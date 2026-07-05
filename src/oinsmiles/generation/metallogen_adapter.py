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
import os
import re
import sys

import numpy as np
from rdkit import Chem

from ..generator3d import generate_3d_structures, get_xyz_string, globalvars
from .molassembler_adapter import GeneratedStructure
from .oin_parser import OINParser, ParsedOIN

logger = logging.getLogger(__name__)

# Force-field convergence presets for the constrained MMFF/UFF cleanup. Each maps
# to TMCOptimizer kwargs: ff_max_iters / ff_force_tol / ff_energy_tol feed
# RDKit's ForceField.Minimize(); d_converge is the scan-level geometry tolerance
# (Angstrom); num_relaxation is the number of pull-in scan cycles. Tighter presets
# relax the ligand backbones closer to the FF minimum at the cost of runtime.
# Select via OIN3DGenerator(engine="metallogen", ff_preset="tight"), the
# OIN_FF_PRESET env var, or verify_roundtrip.py --ff-preset.
FF_PRESETS = {
    "loose": {"ff_max_iters": 100, "ff_force_tol": 1e-3, "ff_energy_tol": 1e-4, "d_converge": 0.10},
    "default": {"ff_max_iters": 200, "ff_force_tol": 1e-4, "ff_energy_tol": 1e-6, "d_converge": 0.05},
    "tight": {
        "ff_max_iters": 2000,
        "ff_force_tol": 1e-5,
        "ff_energy_tol": 1e-7,
        "d_converge": 0.02,
        "num_relaxation": 8,
    },
    "very_tight": {
        "ff_max_iters": 10000,
        "ff_force_tol": 1e-6,
        "ff_energy_tol": 1e-8,
        "d_converge": 0.01,
        "num_relaxation": 12,
    },
}


def _resolve_ff_params(ff_preset=None, ff_params=None):
    """Merge a named FF preset (arg or OIN_FF_PRESET env) with explicit overrides."""
    preset = ff_preset or os.environ.get("OIN_FF_PRESET")
    resolved = {}
    if preset:
        if preset not in FF_PRESETS:
            raise ValueError(
                f"Unknown ff_preset {preset!r}; choose from {sorted(FF_PRESETS)}"
            )
        resolved.update(FF_PRESETS[preset])
    if ff_params:
        resolved.update(ff_params)  # explicit kwargs win over the preset
    return resolved or None


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
            elif atom.GetSymbol() in ("O", "S"):
                # Anionic / oxo chalcogen donor (enolate, alkoxide, oxo, thiolate):
                # the covalent metal bond replaces what would be an implicit H.
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
            elif atom.GetSymbol() == "N":
                # Amide / imide N (>= 2 heavy neighbours) is an anionic X-type donor;
                # strip its H. Dative amines (NH3, NHR2) have < 2 heavy neighbours -> keep.
                heavy = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
                if heavy >= 2:
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
        # Label CIP so build_contract_mol can carry encoded sp3 stereo (the
        # template's @/@@ is the ground truth vs the stochastic embed handedness).
        try:
            Chem.AssignStereochemistry(t, cleanIt=True, force=True)
        except Exception:
            pass
        templates.append(t)
    return templates


def _flatten_template(t):
    """Connectivity-only copy of a template (all single bonds, no aromatic/charge).

    Used only for substructure matching, so bond orders can be transferred even
    from templates that never sanitize (e.g. C#O, O valence 3)."""
    ft = Chem.RWMol(t)
    for b in ft.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in ft.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
    m = ft.GetMol()
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return m


def build_contract_mol(parsed: ParsedOIN, mg_mol) -> "Chem.Mol | None":
    """Build a contract-compliant RDKit mol from a MetalloGen result.

    Connectivity + coordinates come from MetalloGen (``adj_matrix``, atom_list);
    bond orders + aromaticity + formal charges are transferred per fragment from
    the OIN ligand-fragment SMILES via a connectivity substructure match (robust
    to non-sanitizable templates such as C#O); stereo is perceived from the 3D
    geometry. Metal is at its native index with DATIVE metal->donor bonds. Returns
    None on any failure (caller falls back to coordinate-only output).
    """
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
        flats = [_flatten_template(t) for t in templates]
        used = [False] * len(templates)
        # global contract-atom idx -> encoded CIP code for sp3 carbon stereocentres
        carbon_stereo_targets: dict[int, str] = {}

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
            try:
                qmol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
            for ti, t in enumerate(templates):
                if used[ti] or t.GetNumAtoms() != qmol.GetNumAtoms():
                    continue
                # Connectivity-only match, then copy real bond orders / aromaticity /
                # charge from the template (works even when the template can't sanitize).
                match = qmol.GetSubstructMatch(flats[ti])
                if not match or len(match) != t.GetNumAtoms():
                    continue
                for b in t.GetBonds():
                    rb = rw.GetBondBetweenAtoms(
                        q2g[match[b.GetBeginAtomIdx()]], q2g[match[b.GetEndAtomIdx()]]
                    )
                    rb.SetBondType(b.GetBondType())
                    rb.SetIsAromatic(b.GetIsAromatic())
                for ai in range(t.GetNumAtoms()):
                    rwa = rw.GetAtomWithIdx(q2g[match[ai]])
                    ta = t.GetAtomWithIdx(ai)
                    if ta.GetIsAromatic():
                        rwa.SetIsAromatic(True)
                    rwa.SetFormalCharge(ta.GetFormalCharge())
                    # Record encoded CIP for sp3 carbon stereocentres so we can
                    # override the (stochastic) embed handedness below. P/N Zone-A
                    # donors are left to ChiralityRecoveryUtility.
                    if ta.GetAtomicNum() == 6 and ta.HasProp("_CIPCode"):
                        carbon_stereo_targets[q2g[match[ai]]] = ta.GetProp("_CIPCode")
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

        # Carry ENCODED sp3-carbon stereo. The embed picks a random handedness
        # at backbone stereocentres, so 3D-perceived stereo can be the enantiomer
        # of what the OIN fragment SMILES encodes. Where the geometry-derived CIP
        # disagrees with the template CIP, flip the tag -- the perceive-then-flip
        # pattern ChiralityRecoveryUtility already uses for Zone-A P (a bare tag
        # flip does not survive get_oin_string's fragment rebuild; a CIP-validated
        # one does). Bounded fixed point: flipping one centre can change another's
        # CIP priority ranking.
        if carbon_stereo_targets:
            _CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
            _CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            for _ in range(3):
                try:
                    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                except Exception:
                    break
                changed = False
                for gidx, want in carbon_stereo_targets.items():
                    a = mol.GetAtomWithIdx(gidx)
                    cur = a.GetPropsAsDict().get("_CIPCode")
                    tag = a.GetChiralTag()
                    if cur and cur != want and tag in (_CW, _CCW):
                        a.SetChiralTag(_CCW if tag == _CW else _CW)
                        changed = True
                if not changed:
                    break
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
        ff_preset: str | None = None,
        ff_params: dict | None = None,
    ) -> None:
        self.timeout = timeout
        self.dg_strategy = dg_strategy
        self.ensemble_size = ensemble_size
        # optimizer=None -> FF-relaxed geometry only (default; always available).
        # optimizer="xtb" -> refine the FF pool with GFN2-xTB and energy-rank
        # (requires xtb-python + ase; degrades gracefully to FF if unavailable).
        self.optimizer = optimizer
        # FF convergence knobs (named preset + optional explicit overrides).
        self.ff_params = _resolve_ff_params(ff_preset, ff_params)

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
                ff_params=self.ff_params,
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
        ff_preset: str | None = None,
        ff_params: dict | None = None,
    ) -> None:
        self.parser = OINParser()
        self.adapter = MetalloGenAdapter(
            timeout=timeout,
            dg_strategy=dg_strategy,
            ensemble_size=ensemble_size,
            optimizer=optimizer,
            ff_preset=ff_preset,
            ff_params=ff_params,
        )

    def generate(self, oin_string: str) -> GeneratedStructure:
        return self.adapter.generate(self.parser.parse(oin_string))
