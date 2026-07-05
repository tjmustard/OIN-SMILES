"""MetalloGen 3D-generation backend adapter.

Bridges OIN-SMILES (``ParsedOIN``) to the vendored MetalloGen engine in
``oinsmiles.generator3d`` (dummy-metal + RDKit ``CoordMap`` embed + constrained
MMFF/UFF cleanup). Selected via ``OIN3DGenerator(engine="metallogen")``; the
legacy Molassembler/stitch backend remains the default.

Slot mapping: OIN encodes an explicit per-fragment coordination vector, while
MetalloGen assigns m-SMILES fragments to fixed ``globalvars`` coordinate slots by
atom-map number. Matching each OIN vector to its nearest MetalloGen slot vector is
what preserves geometric isomerism (e.g. keeps cisplatin cis, not trans).
"""

import contextlib
import logging
import re
import sys

import numpy as np
from rdkit import Chem

from ..generator3d import generate_3d_structures, get_xyz_string, globalvars
from .molassembler_adapter import GeneratedStructure
from .oin_parser import OINParser, ParsedOIN

logger = logging.getLogger(__name__)

OIN_TO_METALLOGEN_GEO = {
    "LIN": "2_linear",
    "TPL": "3_trigonal_planar",
    "SQP": "4_square_planar",
    "SPL": "4_square_planar",
    "TET": "4_tetrahedral",
    "TPY": "4_trigonal_pyramidal",
    "SPY": "5_square_pyramidal",
    "TBP": "5_trigonal_bipyramidal",
    "OCT": "6_octahedral",
    "PBP": "7_pentagonal_bipyramidal",
}


def convert_parsed_to_msmiles(parsed: ParsedOIN) -> str:
    """Convert a ``ParsedOIN`` into MetalloGen's ``metal|lig1|...|geo`` m-SMILES.

    Each ligand binding atom is tagged with the atom-map number of the nearest
    MetalloGen coordinate slot (isomerism-preserving nearest-vector match).
    """
    geo = OIN_TO_METALLOGEN_GEO.get(parsed.geo_code, "")
    if not geo:
        raise ValueError(
            f"Geometry code '{parsed.geo_code}' not supported by MetalloGen mapping."
        )

    metallogen_vectors = globalvars.known_geometries_vector_dict[geo]
    num_slots = len(metallogen_vectors)
    ligand_parts = [None] * num_slots

    # Metal fragment: strip OIN annotations (e.g. ``[Pt_SPL]`` -> ``[Pt]``).
    metal_frag = parsed.fragments[parsed.metal_fragment_idx]
    metal_frag = re.sub(r"_[A-Z0-9]+", "", metal_frag)
    metal_frag = re.sub(r"@SP[0-9]+", "", metal_frag)

    for i, frag_smiles in enumerate(parsed.fragments):
        if i == parsed.metal_fragment_idx:
            continue

        mol = Chem.MolFromSmiles(frag_smiles, sanitize=False)
        if mol is None:
            raise ValueError(f"Failed to parse fragment {i}: {frag_smiles}")
            
        mol.UpdatePropertyCache(strict=False)

        # Fix kekulization for neutral radicals (like Cp)
        try:
            Chem.Kekulize(mol)
        except Exception:
            # If it fails to kekulize, it's likely an aromatic ring missing a charge (like Cp)
            # Try adding a -1 charge to one atom in each 5-membered aromatic ring
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    mol.GetAtomWithIdx(ring[0]).SetFormalCharge(-1)
            try:
                Chem.Kekulize(mol)
            except Exception:
                pass # If it still fails, let it be

        frag_vectors = [v for v in parsed.vectors if v.fragment_idx == i]

        for v in frag_vectors:
            target_vec = np.array(v.vector)
            dists = np.linalg.norm(metallogen_vectors - target_vec, axis=1)
            mg_slot_idx = int(np.argmin(dists))
            
            # Heuristic to strip implicit Hs from C#N or C#O
            atom = mol.GetAtomWithIdx(v.atom_in_fragment_idx)
            if atom.GetSymbol() == 'C':
                if any(b.GetBondType() == Chem.BondType.TRIPLE for b in atom.GetBonds()):
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(0)
                elif atom.GetIsAromatic():
                    # Count how many coordinating atoms are in the same ring
                    ring_info = mol.GetRingInfo()
                    coordinating_indices = [vec.atom_in_fragment_idx for vec in frag_vectors]
                    is_haptic = False
                    for ring in ring_info.AtomRings():
                        if atom.GetIdx() in ring:
                            coord_in_ring = sum(1 for idx in ring if idx in coordinating_indices)
                            if coord_in_ring > 2: # Cp has 5, benzene has 6. If >2 it's definitely haptic
                                is_haptic = True
                                break
                    if not is_haptic:
                        atom.SetNoImplicit(True)
                        atom.SetNumExplicitHs(0)
            
            # MetalloGen map numbers are 1-based (slot index + 1).
            atom.SetAtomMapNum(mg_slot_idx + 1)

        mapped_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        # Place the fragment at its first binding slot (monodentate = its only slot;
        # multidentate carries all its map numbers within this one fragment string).
        first_slot = int(
            np.argmin(
                np.linalg.norm(metallogen_vectors - np.array(frag_vectors[0].vector), axis=1)
            )
        )
        ligand_parts[first_slot] = mapped_smiles

    msmiles_parts = [metal_frag] + [p for p in ligand_parts if p is not None]
    return "|".join(msmiles_parts) + f"|{geo}"


def convert_oin_to_msmiles(oin_string: str) -> str:
    """Parse an OIN string and convert it to MetalloGen m-SMILES."""
    return convert_parsed_to_msmiles(OINParser().parse(oin_string))


def _oin_fragment_templates(parsed: ParsedOIN) -> list:
    """Heavy-atom RDKit templates (correct bond orders + stereo) per OIN ligand."""
    templates = []
    for k, frag in enumerate(parsed.fragments):
        if k == parsed.metal_fragment_idx:
            continue
        t = Chem.MolFromSmiles(frag, sanitize=False)
        if t is None:
            continue
        try:
            Chem.SanitizeMol(t)
        except Exception:
            t.UpdatePropertyCache(strict=False)
        try:
            t = Chem.RemoveHs(t)
        except Exception:
            pass
        templates.append(t)
    return templates


def build_contract_mol(parsed: ParsedOIN, mg_mol) -> "Chem.Mol | None":
    """Build a contract-compliant RDKit mol from a MetalloGen result.

    Connectivity + coordinates come from MetalloGen (``adj_matrix``, atom_list);
    bond orders + aromaticity come from the OIN ligand-fragment SMILES (via
    ``AssignBondOrdersFromTemplate``); stereo is perceived from the 3D geometry.
    Metal is at its native index with DATIVE metal->donor bonds. Returns None on
    any failure (caller falls back to coordinate-only output).
    """
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    from ..utils.xyz2mol import TRANSITION_METALS_NUM

    try:
        syms = [a.get_element() for a in mg_mol.atom_list]
        coords = [a.get_coordinate() for a in mg_mol.atom_list]
        adj = np.array(mg_mol.adj_matrix)
        n = len(syms)
        if adj.shape != (n, n) or n == 0:
            return None

        rw = Chem.RWMol()
        for s in syms:
            rw.AddAtom(Chem.Atom(s))
        for i in range(n):
            for j in range(i + 1, n):
                if adj[i, j] > 0:
                    rw.AddBond(i, j, Chem.BondType.SINGLE)
        conf = Chem.Conformer(n)
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        rw.AddConformer(conf, assignId=True)

        metal_idx = next(
            (i for i in range(n) if rw.GetAtomWithIdx(i).GetAtomicNum() in TRANSITION_METALS_NUM),
            None,
        )
        if metal_idx is None:
            return None
        donors = [b.GetOtherAtomIdx(metal_idx) for b in rw.GetAtomWithIdx(metal_idx).GetBonds()]

        frag_rw = Chem.RWMol(rw)
        for d in donors:
            frag_rw.RemoveBond(metal_idx, d)
        mapping = []
        frag_mols = Chem.GetMolFrags(
            frag_rw, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping
        )

        templates = _oin_fragment_templates(parsed)
        used = [False] * len(templates)

        for fi, fm in enumerate(frag_mols):
            orig = mapping[fi]
            if metal_idx in orig:
                continue
            heavy_local = [
                k for k in range(fm.GetNumAtoms()) if fm.GetAtomWithIdx(k).GetAtomicNum() != 1
            ]
            q = Chem.RWMol()
            q2g = []
            loc = {}
            for k in heavy_local:
                loc[k] = q.AddAtom(Chem.Atom(fm.GetAtomWithIdx(k).GetAtomicNum()))
                q2g.append(orig[k])
            for b in fm.GetBonds():
                a, bb = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                if a in loc and bb in loc:
                    q.AddBond(loc[a], loc[bb], Chem.BondType.SINGLE)
            qmol = q.GetMol()
            for ti, t in enumerate(templates):
                if used[ti] or t.GetNumAtoms() != qmol.GetNumAtoms():
                    continue
                try:
                    fixed = AllChem.AssignBondOrdersFromTemplate(t, qmol)
                except Exception:
                    continue
                for b in fixed.GetBonds():
                    rb = rw.GetBondBetweenAtoms(q2g[b.GetBeginAtomIdx()], q2g[b.GetEndAtomIdx()])
                    rb.SetBondType(b.GetBondType())
                    rb.SetIsAromatic(b.GetIsAromatic())
                for a in range(fixed.GetNumAtoms()):
                    if fixed.GetAtomWithIdx(a).GetIsAromatic():
                        rw.GetAtomWithIdx(q2g[a]).SetIsAromatic(True)
                used[ti] = True
                break

        for d in donors:
            rw.GetBondBetweenAtoms(metal_idx, d).SetBondType(Chem.BondType.DATIVE)

        mol = rw.GetMol()
        # No full sanitize: the OIN encoder allows non-standard valences (C#O,
        # charge-less Cp). Perceive rings + 3D stereo leniently.
        for step in (
            lambda: mol.UpdatePropertyCache(strict=False),
            lambda: Chem.GetSymmSSSR(mol),
            lambda: Chem.AssignStereochemistryFrom3D(mol),
        ):
            try:
                step()
            except Exception:
                pass
        return mol
    except Exception:
        return None


class MetalloGenAdapter:
    """Generation backend mirroring ``MolassemblerAdapter.generate(parsed)``."""

    def __init__(
        self,
        timeout: int = 60,
        dg_strategy: str = "single",
        ensemble_size: int = 1,
        optimizer: str | None = None,
    ) -> None:
        self.timeout = timeout
        self.dg_strategy = dg_strategy
        self.ensemble_size = ensemble_size
        # optimizer=None -> FF-relaxed geometry only (default; always available).
        # optimizer="xtb" -> refine the FF pool with GFN2-xTB and energy-rank
        # (requires xtb-python + ase; degrades gracefully to FF if unavailable).
        self.optimizer = optimizer

    def generate(self, parsed: ParsedOIN) -> GeneratedStructure:
        msmiles = convert_parsed_to_msmiles(parsed)
        logger.debug("OIN %r -> m-SMILES %r", parsed.original_oin, msmiles)

        # The MetalloGen engine prints progress/geometry to stdout; redirect it to
        # stderr so the oin2xyz CLI's stdout (the XYZ block) stays clean.
        with contextlib.redirect_stdout(sys.stderr):
            mols = generate_3d_structures(
                msmiles,
                num_conformers=self.ensemble_size,
                optimizer=self.optimizer,
            )
        if not mols:
            raise ValueError(
                f"MetalloGen failed to generate any conformers for m-SMILES {msmiles!r}"
            )

        xyz_str = get_xyz_string(mols[0])
        # Contract mol: MetalloGen connectivity+coords, OIN bond orders + 3D stereo.
        # None on failure -> callers fall back to coordinate re-perception.
        mol = build_contract_mol(parsed, mols[0])
        return GeneratedStructure(xyz=xyz_str, mol=mol)


class OIN3DGeneratorMetallogen:
    """Standalone parse+generate wrapper.

    Retained for direct use; prefer ``OIN3DGenerator(engine="metallogen")`` for
    the integrated seam.
    """

    def __init__(
        self,
        timeout: int = 60,
        ensemble_size: int = 1,
        dg_strategy: str = "single",
        optimizer: str | None = None,
    ) -> None:
        self.parser = OINParser()
        self.adapter = MetalloGenAdapter(
            timeout=timeout,
            dg_strategy=dg_strategy,
            ensemble_size=ensemble_size,
            optimizer=optimizer,
        )

    def generate(self, oin_string: str) -> GeneratedStructure:
        return self.adapter.generate(self.parser.parse(oin_string))
