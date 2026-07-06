import itertools
import logging
from collections import defaultdict

import numpy as np
from scipy.spatial.transform import Rotation

from ..oin.winding import signed_circulation

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
except ImportError:
    pass  # RDKit required for Sanitizer, but Aligner is pure numpy


# ==========================================
# PART 1: SMILES SANITIZER (The V2.4 Fix)
# ==========================================
class OINSanitizer:
    """Generate stable, drift-free canonical SMILES for ligand fragments."""

    @staticmethod
    def generate_robust_smiles(ligand_mol, binding_indices_in_ligand):
        """Generate a canonical SMILES with all binding atoms in explicit brackets.

        Forcing explicit brackets prevents 'drift' where c1 becomes [cH]1 or
        vice versa between runs.

        Args:
            ligand_mol: RDKit molecule object of the ligand fragment.
            binding_indices_in_ligand: Atom indices (int) in the ligand that
                bond to the metal.
        """
        # Create a modifiable copy
        rw_mol = Chem.RWMol(ligand_mol)

        # 0. Restore aromatic bond type on ring bonds between aromatic-flagged
        # atoms that are still typed SINGLE. The OIN→XYZ generator de-aromatizes
        # rings for ETKDG (aromatic → SINGLE) and never restores them, so a
        # re-encode of the generator's bonded mol would otherwise serialize an
        # aromatic Cp as [cH]-[cH]-... . Guard on IsInRing() so genuine biaryl
        # single bonds (e.g. Ir(ppy)3 phenyl-pyridine, BINAP binaphthyl) are
        # left untouched. This is a no-op on fragments perceived from real 3D
        # geometry (already AROMATIC-typed), so pass-1 output is unchanged.
        # A staged SanitizeMol was rejected here: a charge-less Cp anion raises
        # KekulizeException, and SANITIZE_ALL ^ KEKULIZE + SetAromaticity leaves
        # the bond types SINGLE (see spec/worklog/TASK-42 probe).
        for bond in rw_mol.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.SINGLE
                and bond.IsInRing()
                and bond.GetBeginAtom().GetIsAromatic()
                and bond.GetEndAtom().GetIsAromatic()
            ):
                bond.SetBondType(Chem.BondType.AROMATIC)
                bond.SetIsAromatic(True)

        # 1. Force Explicit H attributes on Zone A atoms
        for idx in binding_indices_in_ligand:
            atom = rw_mol.GetAtomWithIdx(idx)

            # Get current H count (explicit + implicit)
            # We want to lock the current state so it doesn't drift.
            total_h = atom.GetTotalNumHs()

            # Force this to be Explicit.
            # This forces RDKit to write brackets like [cH] or [CH3]
            # instead of c or C.
            atom.SetNumExplicitHs(total_h)
            atom.SetNoImplicit(True)

            # Update property cache to check valence state
            try:
                atom.GetOwningMol().UpdatePropertyCache(strict=False)
            except Exception:
                pass

            # Check for valence deficit (to handle radicals like [Cl] vs Cl->HCl)
            # Only do this if not already set (though usually 0)
            if atom.GetNumRadicalElectrons() == 0:
                pt = Chem.GetPeriodicTable()
                default_val = pt.GetValenceList(atom.GetAtomicNum())[0]
                # Use GetTotalValence() which includes our set ExplicitHs after UpdatePropertyCache
                current_val = atom.GetTotalValence()

                if current_val < default_val:
                    deficit = default_val - current_val
                    atom.SetNumRadicalElectrons(deficit)

        # 2. Generate Canonical SMILES
        # isomericSmiles=True ensures we keep stereochem info if present
        kmol = rw_mol.GetMol()
        smiles = Chem.MolToSmiles(kmol, isomericSmiles=True, canonical=True)
        return smiles, kmol


# ==========================================
# PART 2: DISCRETE ALIGNER (V2.4 Logic)
# ==========================================


def normalize_template(arr):
    """Return the template with each row normalized to a unit vector."""
    return arr / np.linalg.norm(arr, axis=1)[:, None]


# Define Templates with Explicit Indices (Order matters!)


def normalize_cols(arr):
    """Return the array with each row normalized to a unit vector."""
    return arr / np.linalg.norm(arr, axis=1)[:, None]


TEMPLATE_SPECS = {
    "LIN": {0: {"pos": [0, 0, 1], "ref": [1, 0, 0]}, 1: {"pos": [0, 0, -1], "ref": [1, 0, 0]}},
    "TPL": {
        0: {"pos": [0, 1, 0], "ref": [0, 0, 1]},
        1: {"pos": [0.8660254, -0.5, 0], "ref": [0, 0, 1]},
        2: {"pos": [-0.8660254, -0.5, 0], "ref": [0, 0, 1]},
    },
    "SPL": {
        0: {"pos": [1, 0, 0], "ref": [0, 0, 1]},
        1: {"pos": [0, 1, 0], "ref": [0, 0, 1]},
        2: {"pos": [-1, 0, 0], "ref": [0, 0, 1]},
        3: {"pos": [0, -1, 0], "ref": [0, 0, 1]},
    },
    "TET": {
        0: {"pos": [1, 1, 1], "ref": [-1, 1, 0]},
        1: {"pos": [1, -1, -1], "ref": [0, 1, -1]},
        2: {"pos": [-1, 1, -1], "ref": [1, 1, 0]},
        3: {"pos": [-1, -1, 1], "ref": [1, 0, 1]},
    },
    "TPY": {
        0: {"pos": [0, 0, 1], "ref": [1, 0, 0]},
        1: {"pos": [0, 1, 0], "ref": [0, 0, 1]},
        2: {"pos": [0.8660254, -0.5, 0], "ref": [0, 0, 1]},
        3: {"pos": [-0.8660254, -0.5, 0], "ref": [0, 0, 1]},
    },
    "TBP": {
        0: {"pos": [0, 0, 1], "ref": [1, 0, 0]},
        1: {"pos": [0, 0, -1], "ref": [1, 0, 0]},
        2: {"pos": [0, 1, 0], "ref": [0, 0, 1]},
        3: {"pos": [0.8660254, -0.5, 0], "ref": [0, 0, 1]},
        4: {"pos": [-0.8660254, -0.5, 0], "ref": [0, 0, 1]},
    },
    "SPY": {
        # Real square-pyramidal complexes lift the metal toward the apical ligand,
        # so the four basal donors sit BELOW the metal's equatorial plane
        # (apex-basal angle ~105 deg, trans-basal ~150 deg), not coplanar (90/180).
        # An idealized 90 deg template loses the angular RMSD to TBP for puckered
        # cases like the vanadyl VOacac2 (apex-basal 104-107, trans-basal 146-152).
        # Basal pos = [sin(105)*cos(phi), sin(105)*sin(phi), cos(105)] keeps the
        # +-x/+-y directions (so slot numbering / apex-at-{0} is unchanged) while
        # matching real geometry. Genuine TBP still wins its own cases via a true
        # ~180 deg axial pair, which a puckered SPY (max trans ~150) cannot match.
        0: {"pos": [0, 0, 1], "ref": [1, 0, 0]},
        1: {"pos": [0.9659258, 0, -0.2588190], "ref": [0, 0, 1]},
        2: {"pos": [-0.9659258, 0, -0.2588190], "ref": [0, 0, 1]},
        3: {"pos": [0, 0.9659258, -0.2588190], "ref": [0, 0, 1]},
        4: {"pos": [0, -0.9659258, -0.2588190], "ref": [0, 0, 1]},
    },
    "OCT": {
        0: {"pos": [0, 0, 1], "ref": [1, 0, 0]},
        1: {"pos": [0, 0, -1], "ref": [1, 0, 0]},
        2: {"pos": [1, 0, 0], "ref": [0, 0, 1]},
        3: {"pos": [-1, 0, 0], "ref": [0, 0, 1]},
        4: {"pos": [0, 1, 0], "ref": [0, 0, 1]},
        5: {"pos": [0, -1, 0], "ref": [0, 0, 1]},
    },
    "PBP": {
        0: {"pos": [0, 0, 1], "ref": [1, 0, 0]},
        1: {"pos": [0, 0, -1], "ref": [1, 0, 0]},
        2: {"pos": [1, 0, 0], "ref": [0, 0, 1]},
        3: {"pos": [0.309017, 0.951057, 0], "ref": [0, 0, 1]},
        4: {"pos": [-0.809017, 0.587785, 0], "ref": [0, 0, 1]},
        5: {"pos": [-0.809017, -0.587785, 0], "ref": [0, 0, 1]},
        6: {"pos": [0.309017, -0.951057, 0], "ref": [0, 0, 1]},
    },
}

# List of Symmetric Ligands that should have a fixed Heading Atom (The first binding atom in SMILES)
# This overrides the geometric alignment which can be unstable for symmetric ligands.
SYMMETRIC_LIGANDS = {
    "C=C",
    "[CH2]=[CH2]",
    "c1cccc1",
    "C1=C-C=C-[CH-]1",
    "[cH]1[cH][cH][cH][cH]1",
    "c1ccccc1",
    "[cH]1[cH][cH][cH][cH][cH]1",
    "c1cccccc1",
    "[cH]1[cH][cH][cH][cH][cH][cH]1",
}

# Generate Legacy TEMPLATES (Pos only) for Aligner compatibility
TEMPLATES = {}
for geo, specs in TEMPLATE_SPECS.items():
    sorted_slots = sorted(specs.keys())
    pos_vecs = np.array([specs[idx]["pos"] for idx in sorted_slots], dtype=float)
    if geo == "TET":
        pos_vecs = normalize_cols(pos_vecs)
        for i, vec in enumerate(pos_vecs):
            TEMPLATE_SPECS[geo][sorted_slots[i]]["pos"] = vec.tolist()
    TEMPLATES[geo] = pos_vecs


class OINDiscreteAligner:
    """Assign ligand binding atoms to discrete coordination-geometry slots."""

    def __init__(self, metal_idx, ligands):
        """Initialize the aligner with the metal index and ligand data.

        Args:
            metal_idx: Index of the metal center atom.
            ligands: List of dicts, each of which MUST contain 'smiles'
                generated by OINSanitizer and a 'binding_atoms' list of
                tuples/lists: [(global_idx, mass, coords, local_idx)].
        """
        self.metal_idx = metal_idx
        self.ligands = ligands

    def generate_canonical_vectors(self):
        """Reduce hapticity and assign binding atoms to canonical slot vectors."""
        # 1. Haptic Reduction
        virtual_atoms = self._reduce_hapticity()

        # 2. Competitive Geometry Detection
        n_eff = len(virtual_atoms)
        # Handle cases where N_eff is small or unexpected, defaulting to LIN
        # if needed or raising error
        # Assuming n_eff corresponds to one of the templates for now

        # Determine candidates based on N_eff
        # This mapping logic might need to be more robust for N=1 etc, but following PRD guidance
        if n_eff < 2:
            # Fallback or specific handling for N=1? usually LIN with one empty slot or just LIN
            # For now let's assume valid N >= 2
            tmpl_name = "LIN"
            tmpl_vectors = TEMPLATES["LIN"]
            mapping = [None, None]
            R_mat = None
        else:
            best_res = self._find_best_geometry_match(n_eff, virtual_atoms)
            if best_res:
                tmpl_name, tmpl_vectors, mapping, R_mat = best_res
            else:
                # Fallback if no match found (e.g. geometry too distorted)
                return "g:NON|w:NON"

        # 3. Canonicalize (Maximization + Homogeneous Sort + Heading Syntax)
        canonical_str = self._permute_and_serialize(
            mapping, tmpl_vectors, geometry_name=tmpl_name, alignment_rotation=R_mat
        )

        return f"g:{tmpl_name}|w:{canonical_str}"

    def _reduce_hapticity(self):
        virtual_atoms = []
        for i, lig in enumerate(self.ligands):
            # Skip Metal Center (It has no binding atoms to itself)
            if i == self.metal_idx:
                continue

            # Skip fragments with no binding atoms (non-coordinating solvents?)
            if not lig.get("binding_atoms"):
                continue

            # Sort Key uses the SANITIZED smiles
            # ligand 'binding_atoms' structure: [(global_idx, mass, coords_array, local_idx)]
            # We need the mass of the first binding atom? Or heaviest?
            # The PRD says "Atomic Mass of the Binding Atom".
            # Let's assume 'binding_atoms' is already sorted or we pick the
            # first one which is usually primary.
            # Using the first one from the supplied list.

            first_binding_atom_mass = lig["binding_atoms"][0][1]

            binding_coords = np.array([ba[2] for ba in lig["binding_atoms"]])
            # We need to track which atom is which. The 'binding_atoms' list has entries.
            # We need LOCAL indices for the output string (Rank.LocalIdx:Slot).
            # Let's store the full atom info.

            # Retrieve Metal Center coordinates
            metal_frag = self.ligands[self.metal_idx]
            metal_origin = metal_frag.get("metal_coords")
            if metal_origin is None:
                # Fallback to origin if not found (should not happen with updated xyz2mol)
                metal_origin = np.array([0.0, 0.0, 0.0])

            zone_a_info = lig["binding_atoms"]  # list of (idx, mass, coords, local_idx)
            n_b = len(binding_coords)

            groups = []
            visited = set()
            for j in range(n_b):
                if j in visited:
                    continue
                stack = [j]
                component = []
                while stack:
                    curr = stack.pop()
                    if curr in visited:
                        continue
                    visited.add(curr)
                    component.append(curr)
                    for k in range(n_b):
                        if k in visited:
                            continue
                        if np.linalg.norm(binding_coords[curr] - binding_coords[k]) < 1.6:
                            stack.append(k)
                groups.append(component)

            for grp in groups:
                # Sort group by local_idx to ensure alignment between coords and indices
                grp.sort(key=lambda k: zone_a_info[k][3])

                grp_coords = binding_coords[grp]
                if len(grp_coords) == 0:
                    continue
                centroid = np.mean(grp_coords, axis=0)

                # Representative Atom: The one with the lowest local index (now first)
                rep_idx = zone_a_info[grp[0]][3]

                # Get ALL local indices for this group (for w-tag expansion)
                constituent_indices = sorted([zone_a_info[k][3] for k in grp])

                virtual_atoms.append(
                    {
                        "rank": i,  # Ligand Rank
                        "local_idx": rep_idx,
                        # List of all atoms in this haptic group
                        "constituent_indices": constituent_indices,
                        "coords": centroid - metal_origin,
                        "group_coords": grp_coords
                        - metal_origin,  # Store centered group coords for heading calc
                        # Identity key for Homogeneous Sorting: (Mass, SMILES)
                        "chem_id": (first_binding_atom_mass, lig["smiles"]),
                    }
                )

        return virtual_atoms

    def _find_best_geometry_match(self, n, virtual_atoms):
        candidates = []
        # Exhaustive Candidate List based on Coordination Number (N)
        if n == 2:
            candidates = ["LIN"]
        elif n == 3:
            candidates = ["TPL"]
        elif n == 4:
            # Check ALL 4-coordinate geometries
            candidates = ["SPL", "TET", "TPY"]
            logger.debug(f"  N=4: Checking all candidates: {candidates}")
        elif n == 5:
            # Check ALL 5-coordinate geometries
            candidates = ["TBP", "SPY"]
            logger.debug(f"  N=5: Checking all candidates: {candidates}")
        elif n == 6:
            candidates = ["OCT"]
        elif n == 7:
            candidates = ["PBP"]
        else:
            # Fallback or robust handling
            if n > 7:
                candidates = ["OCT"]  # Best effort
            else:
                candidates = ["LIN"]

        min_rmsd = float("inf")
        best_result = None

        logger.debug(f"Geometry Selection (N={n})")

        for name in candidates:
            vectors = TEMPLATES.get(name)
            if vectors is None:
                continue

            # If N > Slots, we can't map one-to-one without leaving some out or overlapping.
            # OIN assumes N <= Slots usually. If N < Slots, we have empty slots.
            if n > len(vectors):
                continue

            mapping, rmsd, R_mat = self._map_to_template(virtual_atoms, vectors)
            logger.debug(f"  Candidate: {name}, RMSD: {rmsd:.4f}")

            if mapping is not None and rmsd < min_rmsd:
                min_rmsd = rmsd
                best_result = (name, vectors, mapping, R_mat)

        if best_result:
            logger.debug(f"  Selected: {best_result[0]} (RMSD {min_rmsd:.4f})")

        return best_result

    def _map_to_template(self, virtual_atoms, template_vectors):
        n_atoms = len(virtual_atoms)
        n_slots = len(template_vectors)

        if n_atoms == 0:
            return [None] * n_slots, 0.0, Rotation.from_matrix(np.eye(3))

        input_vecs = np.array([a["coords"] for a in virtual_atoms])
        input_norms = input_vecs / (np.linalg.norm(input_vecs, axis=1)[:, None] + 1e-9)

        best_rmsd = float("inf")
        best_mapping = None
        best_R = Rotation.from_matrix(np.eye(3))

        import itertools

        perm_iterator = itertools.permutations(range(n_slots), n_atoms)

        for slot_indices in perm_iterator:
            target_vecs = template_vectors[list(slot_indices)]

            try:
                R, rmsd = Rotation.align_vectors(target_vecs, input_norms)
            except Exception:
                continue

            if rmsd < best_rmsd:
                best_rmsd = rmsd
                current_mapping = [None] * n_slots
                for atom_idx, slot_idx in enumerate(slot_indices):
                    current_mapping[slot_idx] = virtual_atoms[atom_idx]

                best_mapping = current_mapping
                best_R = R

        return best_mapping, best_rmsd, best_R

    @staticmethod
    def _fragment_mol_for_canonicalization(smiles):
        """Parse a fragment SMILES into an index-aligned mol for canonicalization.

        Atom indices line up 1:1 with `local_idx`/`constituent_indices`.
        `Chem.MolFromSmiles(smiles, sanitize=False)` is the exact same
        construction `xyz2mol.py` uses to derive `local_idx`
        (`smiles_mol = Chem.MolFromSmiles(sanitized_smiles, sanitize=False)`
        at xyz2mol.py:952), so no substructure-match re-mapping is needed --
        the indices are already the same numbering by construction (RT-2).
        Sanitizes everything except kekulization: Zone-A atoms carry forced
        explicit-H counts (`OINSanitizer`) that can defeat RDKit's kekulizer,
        but ring/aromaticity perception (needed for a meaningful canonical
        rank) does not require kekulization. Returns None on any failure
        (RT-3 fail-safe -- callers must fall back to existing behavior).
        """
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                return None
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
            return mol
        except Exception:
            return None

    @classmethod
    def _canonical_heading_atom(cls, smiles, constituent_indices):
        """Return the constituent atom with the lowest canonical rank (RC2).

        Uses `Chem.CanonicalRankAtoms(breakTies=True)`; returns None if it
        cannot be computed/mapped (RT-2/RT-3 fail-safe).
        """
        mol = cls._fragment_mol_for_canonicalization(smiles)
        if mol is None or mol.GetNumAtoms() <= max(constituent_indices, default=-1):
            return None
        try:
            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
        except Exception:
            return None
        return min(constituent_indices, key=lambda idx: ranks[idx])

    @classmethod
    def _canonical_ring_signature(cls, smiles):
        """Return a heading-independent canonical SMILES for the fragment (RC1).

        Returns None if it cannot be computed (RT-2/RT-3 fail-safe).
        """
        mol = cls._fragment_mol_for_canonicalization(smiles)
        if mol is None:
            return None
        try:
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            return None

    def _permute_and_serialize(
        self, slot_assignment, tmpl_vectors, geometry_name=None, alignment_rotation=None
    ):
        # 1. Symmetry Permutations
        symmetries = self._brute_force_symmetries(tmpl_vectors)
        best_sequence = None
        best_final_map = None  # Stores {Rank: SlotIndex}

        # If slot_assignment is None or empty, return empty
        if not slot_assignment:
            return ""

        # 2. Find Canonical View (Maximization)
        for perm in symmetries:
            current_view_map = []

            for old_slot_idx, atom in enumerate(slot_assignment):
                if atom is None:
                    continue

                # Permutation maps OLD slot to NEW slot
                new_slot_idx = perm[old_slot_idx]
                ideal_vec = tmpl_vectors[new_slot_idx]

                # Sanitize Vector
                v_clean = []
                for val in ideal_vec:
                    if abs(val) < 1e-9:
                        val = 0.0
                    v_clean.append(val)
                ideal_vec = tuple(v_clean)

                current_view_map.append(
                    {
                        "rank": atom["rank"],
                        "local_idx": atom["local_idx"],
                        "constituent_indices": atom.get("constituent_indices", [atom["local_idx"]]),
                        "group_coords": atom.get("group_coords"),  # Propagate!
                        "chem_id": atom["chem_id"],
                        "vec": ideal_vec,
                        "slot": new_slot_idx,
                    }
                )

            if not current_view_map:
                continue

            # 3. Homogeneous Sorting
            grouped = defaultdict(list)
            for item in current_view_map:
                grouped[item["chem_id"]].append(item)

            final_sorted_view = []
            for chem_id, items in grouped.items():
                # Group items by Rank (Fragment) to preserve connectivity
                frag_groups = defaultdict(list)
                for it in items:
                    frag_groups[it["rank"]].append(it)

                target_ranks = sorted(list(frag_groups.keys()))

                available_sets = []
                for rank in target_ranks:
                    f_items = frag_groups[rank]
                    f_items.sort(key=lambda x: x["local_idx"])
                    vec_set = tuple([x["vec"] for x in f_items])
                    available_sets.append({"vec_set": vec_set, "items": f_items})

                # Sort the available sets descending (Maximization)
                available_sets.sort(key=lambda x: x["vec_set"], reverse=True)

                # Assign Best Sets to Lowest Ranks
                for i, target_rank in enumerate(target_ranks):
                    assigned_set = available_sets[i]
                    for atom_data in assigned_set["items"]:
                        final_sorted_view.append(
                            {
                                "rank": target_rank,
                                "local_idx": atom_data["local_idx"],
                                "constituent_indices": atom_data["constituent_indices"],
                                "group_coords": atom_data["group_coords"],
                                "vec": atom_data["vec"],
                                "slot": atom_data["slot"],
                                "chem_id": chem_id,
                            }
                        )

            # Sort by Rank/LocalIdx for comparison
            final_sorted_view.sort(key=lambda x: (x["rank"], x["local_idx"]))
            current_sequence = [x["vec"] for x in final_sorted_view]

            # Maximization Check
            if best_sequence is None or current_sequence > best_sequence:
                best_sequence = current_sequence
                best_final_map = final_sorted_view

        if not best_final_map:
            return "error"

        # 3b. RC1 -- scoped eta-only rank swap (content-canonical fragment
        # order). Only same-mass, single-virtual-atom haptic (eta) fragments
        # -- e.g. two differently-substituted Cp rings -- are eligible;
        # every non-eta fragment and the metal (rank 0) keep their exact
        # rank (US-001.3, RT-1/Option A). This permutes WHICH content occupies
        # the SAME set of rank slots the eligible fragments already occupy;
        # it never introduces, removes, or touches any other rank.
        by_rank_for_swap = defaultdict(list)
        for x in best_final_map:
            by_rank_for_swap[x["rank"]].append(x)

        eta_fragment_ranks_by_mass = defaultdict(list)
        for rank, items in by_rank_for_swap.items():
            if len(items) == 1 and len(items[0].get("constituent_indices", [])) > 1:
                eta_fragment_ranks_by_mass[items[0]["chem_id"][0]].append(rank)

        for _mass, ranks in eta_fragment_ranks_by_mass.items():
            if len(ranks) < 2:
                continue  # only one eta fragment at this mass -- nothing to swap

            candidates = []
            fail_safe = False
            for rank in ranks:
                item = by_rank_for_swap[rank][0]
                signature = self._canonical_ring_signature(item["chem_id"][1])
                if signature is None:
                    fail_safe = True
                    break
                candidates.append((item, signature))

            if fail_safe:
                # RT-3: never a silent partial reorder -- leave this whole
                # same-mass group at its existing (arrival-order) ranks.
                continue

            def _winding_tiebreak(entry):
                # Content-identical-ring tiebreak ONLY (D-RC1). Winding
                # character is start-invariant (RT-4), so any constituent
                # atom works as the provisional star -- the final heading
                # atom need not be known yet.
                item, _sig = entry
                slot_def = TEMPLATE_SPECS.get(geometry_name, {}).get(item["slot"])
                if slot_def is None or item.get("group_coords") is None:
                    return ""
                c_indices = sorted(item["constituent_indices"])
                return self._determine_winding(
                    grp_coords=item["group_coords"],
                    star_idx=c_indices[0],
                    constituent_indices=c_indices,
                    slot_z=np.array(slot_def["pos"]),
                    slot_x_ref=np.array(slot_def["ref"]),
                    alignment_rotation=alignment_rotation,
                )

            candidates.sort(
                key=lambda entry: (
                    entry[1],  # canonical_ring_smiles (heading-independent)
                    _winding_tiebreak(entry),
                    # Final deterministic tiebreak. Not the fixture's original
                    # XYZ atom index (unstable across a generate/re-encode
                    # round trip) -- the ring's own lowest constituent
                    # local_idx, which is derived from canonical SMILES atom
                    # order and is therefore itself embedding-independent.
                    min(entry[0]["constituent_indices"]),
                )
            )

            for new_rank, (item, _sig) in zip(sorted(ranks), candidates):
                item["rank"] = new_rank

        best_final_map.sort(key=lambda x: (x["rank"], x["local_idx"]))

        # 4. Heading Atom Selection (V3.4)
        heading_local_indices = set()

        if geometry_name and alignment_rotation is not None and geometry_name in TEMPLATE_SPECS:
            # Group by rank to handle full haptic fragments
            by_rank = defaultdict(list)
            for x in best_final_map:
                by_rank[x["rank"]].append(x)

            template_spec = TEMPLATE_SPECS[geometry_name]

            # 4a. Content-canonical heading (RC2): for substituted/asymmetric
            # eta rings (NOT in SYMMETRIC_LIGANDS), choose the ring atom with
            # the lowest Chem.CanonicalRankAtoms rank -- a topological
            # property, invariant to 3D embedding orientation (US-002).
            # Falls back to the geometric best_idx loop below if the
            # canonical rank can't be computed/mapped (RT-2/RT-3 fail-safe).
            content_canonical_ranks_handled = set()
            for rank, items in by_rank.items():
                first_item = items[0]
                smiles = first_item["chem_id"][1]
                grp_coords = first_item.get("group_coords")
                constituent_indices = first_item.get("constituent_indices", [])

                if smiles in SYMMETRIC_LIGANDS:
                    continue
                if grp_coords is None or len(grp_coords) < 2:
                    continue

                canonical_idx = self._canonical_heading_atom(smiles, constituent_indices)
                if canonical_idx is None:
                    continue

                heading_local_indices.add((rank, canonical_idx))
                content_canonical_ranks_handled.add(rank)

            for rank, items in by_rank.items():
                if rank in content_canonical_ranks_handled:
                    continue

                first_item = items[0]
                slot_idx = first_item["slot"]

                # Skip if slot not in specs or no ref vector
                if slot_idx not in template_spec or "ref" not in template_spec[slot_idx]:
                    continue

                ref_vec = np.array(template_spec[slot_idx]["ref"])

                # Get Group Coords
                grp_coords = first_item.get("group_coords")
                if grp_coords is None:
                    continue

                # Calculate Centroid
                centroid = np.mean(grp_coords, axis=0)

                if len(grp_coords) < 2:
                    continue

                best_dot = -float("inf")
                best_idx = -1

                ordered_indices = first_item["constituent_indices"]

                for k, coord in enumerate(grp_coords):
                    # Vector from Centroid -> Atom (Molecular Frame)
                    v_mol = coord - centroid

                    # Transform to Template Frame using Alignment Rotation R
                    v_tmpl = alignment_rotation.apply(v_mol)

                    # Normalize
                    norm = np.linalg.norm(v_tmpl)
                    if norm > 1e-6:
                        v_tmpl_n = v_tmpl / norm
                        dot = np.dot(v_tmpl_n, ref_vec)

                        if dot > best_dot:
                            best_dot = dot
                            best_idx = ordered_indices[k]

                if best_idx != -1:
                    heading_local_indices.add((rank, best_idx))

            # 4b. Symmetric Ligand Override
            # For known symmetric ligands, we force the heading atom to be the
            # first one (Index 0 in local SMILES)
            # This ensures deterministic output (e.g. [CH2]{^}=[CH2] instead of arbitrary).
            for rank, items in by_rank.items():
                first_item = items[0]
                smiles = first_item["chem_id"][1]

                if smiles in SYMMETRIC_LIGANDS:
                    # Remove any existing heading assignment for this ligand
                    # (from geometric step above). We want to replace it, not
                    # add to it (though usually only one heading per ligand group)
                    # But 'heading_local_indices' is a set of (rank, idx).
                    # We should clear entries for this rank first to be safe?.
                    # Actually, the geometric block above loop over ranks.
                    # It's cleaner to check SYMMETRIC_LIGANDS *inside* the loop
                    # above, but separating logic is also fine.
                    # Let's just Enforce it here.

                    # Find min index among constituent_indices
                    ordered_indices = sorted(first_item["constituent_indices"])
                    forced_idx = ordered_indices[0]

                    # Remove any other heading for this rank just in case
                    # (e.g. if geometric picked another)
                    to_remove = [idx for r, idx in heading_local_indices if r == rank]
                    for idx in to_remove:
                        heading_local_indices.remove((rank, idx))

                    heading_local_indices.add((rank, forced_idx))
                    logger.debug(
                        f"Forced Heading Atom for Symmetric Ligand {smiles}: Index {forced_idx}"
                    )

        # 5. Serialize Index-Based Format: w:Rank.Idx:Slot
        parts = []
        for x in best_final_map:
            slot = x["slot"]
            rank = x["rank"]

            indices = x.get("constituent_indices", [x["local_idx"]])

            for idx in indices:
                tag = f"{rank}.{idx}:{slot}"

                # Append Heading Marker if this is the chosen atom
                if (rank, idx) in heading_local_indices:
                    # Calculate Winding Direction using V3.6 Algorithm
                    direction_char = ">"  # Default

                    if geometry_name in TEMPLATE_SPECS and slot in TEMPLATE_SPECS[geometry_name]:
                        slot_def = TEMPLATE_SPECS[geometry_name][slot]

                        # Need constituent_indices (ordered) to find next atom
                        c_indices = sorted(x.get("constituent_indices", [idx]))

                        direction_char = self._determine_winding(
                            grp_coords=x.get("group_coords"),
                            star_idx=idx,  # This is the local_idx of the star
                            constituent_indices=c_indices,
                            slot_z=np.array(slot_def["pos"]),
                            slot_x_ref=np.array(slot_def["ref"]),
                            alignment_rotation=alignment_rotation,
                        )

                    tag += direction_char

                parts.append(tag)

        return ";".join(parts)

    def _determine_winding(
        self, grp_coords, star_idx, constituent_indices, slot_z, slot_x_ref, alignment_rotation=None
    ):
        """Determine the OIN V3.6 winding direction of a haptic group.

        Determines if the SMILES winding (Star -> Next) is Clockwise (>) or
        Counter-Clockwise (<) relative to the template Slot Z vector.

        Builds the coords/star index in SMILES-fragment order and delegates
        the sign computation to `signed_circulation` (`oin/winding.py`), the
        single source of truth shared with the generation-side haptic-face
        correction (Stereo Phase 3). The `n < 3` degenerate default now lives
        in that helper, not here.
        """
        # 1. Identify Star Index in the List. `constituent_indices` and
        # `grp_coords` are already in SMILES/fragment order (ascending
        # local_idx -- see `_reduce_hapticity`), matching
        # `signed_circulation`'s ordering contract.
        try:
            star_local_idx = constituent_indices.index(star_idx)
        except ValueError:
            return ">"

        # 2. `slot_z` is the TEMPLATE-FRAME vector from the metal (template
        # origin) outward to this slot -- already metal->centroid outward by
        # construction in TEMPLATE_SPECS. Assert/normalize that convention
        # here so a malformed template can't silently flip the helper's sign.
        slot_z_norm = np.linalg.norm(slot_z)
        if slot_z_norm < 1e-9:
            raise ValueError("slot_z must be a nonzero, outward-facing (metal->centroid) vector")
        axis_template = slot_z / slot_z_norm

        # 3. `grp_coords` is in the MOLECULAR frame (centered on the metal --
        # see `_reduce_hapticity`: 'group_coords': grp_coords - metal_origin).
        # Rather than rotate every coordinate into the template frame, rotate
        # the axis into the molecular frame: rotations preserve dot products,
        # so dot(R.a, b) == dot(a, R^-1.b).
        axis_mol = (
            alignment_rotation.inv().apply(axis_template)
            if alignment_rotation is not None
            else axis_template
        )

        return signed_circulation(grp_coords, star_local_idx, axis_mol)

    def _brute_force_symmetries(self, vectors):
        n = len(vectors)
        valid = set()
        steps = [0, 90, 120, 180, 240, 270]
        for rx, ry, rz in itertools.product(steps, repeat=3):
            R = Rotation.from_euler("xyz", [rx, ry, rz], degrees=True)
            rot = R.apply(vectors)
            perm = [-1] * n
            matches = 0

            # Check if this rotation maps the template to itself
            # Each vector in 'rot' must match a vector in 'vectors'

            for i in range(n):
                dists = np.linalg.norm(vectors - rot[i], axis=1)
                best = np.argmin(dists)
                if dists[best] < 0.1:
                    perm[i] = best
                    matches += 1
            if matches == n:
                valid.add(tuple(perm))
        return sorted(list(valid))


def classify_coordination_geometry(donor_vectors):
    """Return the best-matching OIN geometry code for a set of donor vectors.

    ``donor_vectors`` is an iterable of length-3 array-likes, each the
    metal-centered vector ``donor_pos - metal_pos`` of a single coordinating atom
    (no hapticity reduction -- pass one vector per discrete donor). Returns the
    OIN geometry code (e.g. ``"SPL"``, ``"TET"``, ``"TPY"``, ``"OCT"``) chosen by
    the same discrete-template matcher the XYZ->OIN encoder uses, or ``None`` if
    no template matched.

    This wraps ``OINDiscreteAligner._find_best_geometry_match``, which is pure
    w.r.t. instance state (it reads only its arguments and the module-level
    ``TEMPLATES``), so a throwaway aligner with empty ligands is sufficient. Note
    the candidate set is chosen purely by coordination number (len of the input),
    and for ``n > 7`` the matcher falls back to ``"OCT"`` as a best effort -- so
    callers that need an eta/haptic guard must gate on the expected coordination
    number themselves rather than relying on a ``None`` return.
    """
    virtual_atoms = [{"coords": np.asarray(v, dtype=float)} for v in donor_vectors]
    aligner = OINDiscreteAligner(0, [])  # ligands unused by the matcher
    result = aligner._find_best_geometry_match(len(virtual_atoms), virtual_atoms)
    return result[0] if result else None
