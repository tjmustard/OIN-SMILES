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
        # Phase B reconstructs a contract-compliant bonded mol (metal@0, DATIVE,
        # single conformer). Until then the mol channel is None.
        return GeneratedStructure(xyz=xyz_str, mol=None)


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
