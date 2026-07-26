import itertools
import logging
from collections import defaultdict
from functools import lru_cache

import numpy as np
from scipy.spatial.transform import Rotation

from ..oin.canonical_slots import derive_rotation_group
from ..oin.hydrogen import h_faithful_smiles
from ..oin.levers import lever_enabled
from ..oin.winding import signed_circulation

logger = logging.getLogger(__name__)

#: Canonicalize the winding character across interchangeable, automorphic eta groups
#: (`OINDiscreteAligner._canonical_eta_winding`). Default OFF: it is the one Lane 3
#: change that is NOT byte-identical, because it rewrites the emitted winding for the
#: meso arrangement of a symmetric ansa-metallocene -- two spellings of one achiral
#: compound that the encoder currently picks between by geometry.
#: Read at import, so a test must patch this attribute rather than the environment.
#: The DEFAULT now lives in oin.levers (promoted to ON in v0.4.5) rather than being
#: spelled here -- the old inline `"0") != "0"` form is exactly the drift lever_enabled
#: exists to prevent, since a sibling site used bare truthiness where "0" ENABLED a lever.
CANONICAL_ETA_WINDING = lever_enabled("OIN_CANONICAL_ETA_WINDING")

try:
    from rdkit import Chem
except ImportError:
    pass  # RDKit required for Sanitizer, but Aligner is pure numpy


# ==========================================
# PART 0: ETA-GROUP ORIENTATION SYMMETRY
# ==========================================
# A haptic group's winding character ({n>} / {n<}) says which way its SMILES
# atom order circulates when you look at the ring from the metal -- i.e. which
# FACE of the ring the metal sees. That is only a property of the *structure*
# when the ring cannot be turned over onto itself.
#
# Turning an eta ring over is a 180 deg rotation about an in-plane axis through
# its centroid: a PROPER rotation, which leaves the metal and every other ligand
# where they are and maps the ring's atoms onto the ring's own positions in
# reverse cyclic order. So:
#
#   winding is meaningless (``orientation-free``) for an eta group
#   <=> some automorphism of its ligand fragment reverses the group's cyclic order.
#
# Cp*, benzene, mesitylene and a BPh4- phenyl all satisfy this, so their
# geometric winding is an artifact of whichever face happened to point at the
# metal in this particular embedding, and it flips at random between the input
# structure and a regenerated one. An ansa-bis(indenyl) ring does not satisfy it,
# and there its winding is exactly what tells rac from meso.
#
# Note this is strictly weaker than "every ring atom is in one symmetry class":
# mesitylene's ring has two classes and an arm-substituted Cp* has four, yet both
# are orientation-free. Testing symmetry classes would silently leave those broken.
#
# The test avoids enumerating the automorphism group. A canonical SMILES of a
# vertex-labeled graph is a complete invariant of that labeled graph, so labeling
# the eta cycle 1..n and canonicalizing answers "is the reversed labeling the same
# labeled graph?" directly. Comparing canonical SMILES computed inside one
# interpreter run also keeps this stable across RDKit versions, where a specific
# CanonicalRankAtoms ordering would not be.


def _orientation_symmetry_graph(mol):
    """Rebuild `mol` as the bare graph the reversal test needs.

    Built from scratch rather than edited in place: an inherited atom carries
    aromaticity, chiral tags, radicals and implicit-H flags that leak into both
    `CanonicalRankAtoms` and the SMILES writer, and the writer will silently
    prefer its own implicit-valence guess over a `SetNumExplicitHs` we stored.
    Two graph-identical phenyls then serialize differently and the whole test
    quietly fails open. Only four things survive the rebuild:

    * **Atomic number** and **formal charge** -- real chemical identity.
    * **Bonds, all SINGLE** -- so a Kekule structure's alternating bonds cannot
      fake an asymmetry that the delocalized ring does not have.
    * **Total-H count, carried in the isotope field** as ``H + 1``. H count must
      be part of the graph (else cyclopentadiene's sp3 CH2 is conflated with a
      Cp- ring carbon), and the isotope is the one per-atom field that both the
      canonical ranker and the SMILES writer are obliged to honor.

    Aromatic flags, chiral tags and radicals are dropped: `generate_robust_smiles`
    adds radicals to binding atoms only, which would make a coordinated ring look
    unlike its own uncoordinated twin.

    Valence is tolerated throughout (`UpdatePropertyCache(strict=False)` +
    `FastFindRings`, never `SanitizeMol`): a BPh4- boron carries four bonds and no
    charge in a fragment SMILES and trips RDKit's valence check -- the same failure
    that makes `_fragment_mol_for_canonicalization` return None for borates.

    Returns None on any failure -- callers must fall back to geometric winding.
    """
    try:
        src = Chem.Mol(mol)
        src.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(src)

        out = Chem.RWMol()
        for atom in src.GetAtoms():
            fresh = Chem.Atom(atom.GetAtomicNum())
            fresh.SetFormalCharge(atom.GetFormalCharge())
            fresh.SetNoImplicit(True)
            fresh.SetNumExplicitHs(0)
            fresh.SetIsotope(atom.GetTotalNumHs() + 1)
            out.AddAtom(fresh)
        for bond in src.GetBonds():
            out.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)

        graph = out.GetMol()
        graph.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(graph)
        return graph
    except Exception:
        return None


def _eta_traversal_order(mol, constituent_indices):
    """Order a haptic group's atoms along its own connectivity.

    Returns ``("ring", cycle)`` for a closed cycle (Cp, arene), ``("path", chain)``
    for an open chain (eta3-allyl, eta4-diene), or ``(None, None)`` when the atoms
    are not a single cycle or chain -- e.g. a pincer's three mutually non-bonded
    donors, or one fragment carrying two separate eta rings. Those must never be
    treated as one orientable group.
    """
    atoms = sorted(constituent_indices)
    if len(atoms) < 2:
        return None, None
    members = set(atoms)
    try:
        adj = {
            a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in members]
            for a in atoms
        }
    except Exception:
        return None, None

    degrees = {a: len(neigh) for a, neigh in adj.items()}

    if all(deg == 2 for deg in degrees.values()):
        start = atoms[0]
        cycle = [start]
        prev, cur = None, start
        while True:
            nxt = [x for x in adj[cur] if x != prev]
            if not nxt:
                return None, None
            step = nxt[0]
            if step == start:
                break
            cycle.append(step)
            prev, cur = cur, step
            if len(cycle) > len(atoms):
                return None, None
        return ("ring", cycle) if len(cycle) == len(atoms) else (None, None)

    ends = [a for a, deg in degrees.items() if deg == 1]
    if len(ends) == 2 and all(deg <= 2 for deg in degrees.values()):
        path = [ends[0]]
        prev, cur = None, ends[0]
        while True:
            nxt = [x for x in adj[cur] if x != prev]
            if not nxt:
                break
            path.append(nxt[0])
            prev, cur = cur, nxt[0]
            if len(path) > len(atoms):
                return None, None
        return ("path", path) if len(path) == len(atoms) else (None, None)

    return None, None


def _labeled_canonical_smiles(flat, order):
    """Canonical SMILES of `flat` with `order` labeled 1..n via atom map numbers."""
    try:
        labeled = Chem.RWMol(flat)
        for atom in labeled.GetAtoms():
            atom.SetAtomMapNum(0)
        for position, idx in enumerate(order):
            labeled.GetAtomWithIdx(idx).SetAtomMapNum(position + 1)
        return Chem.MolToSmiles(labeled.GetMol(), canonical=True)
    except Exception:
        return None


def _labeled_forms(flat, order, kind):
    """All canonical labelings of `order` that describe the same traversal sense.

    For a ring every rotation is the same circulation, so the invariant is the SET
    over the n rotations; for a path only the chain itself. Returns None if any
    labeling fails to canonicalize.
    """
    if kind == "ring":
        n = len(order)
        forms = {_labeled_canonical_smiles(flat, order[k:] + order[:k]) for k in range(n)}
    else:
        forms = {_labeled_canonical_smiles(flat, order)}
    return None if None in forms else forms


def _orientation_free_from_mol(mol, constituent_indices):
    """True if reversing this haptic group's cyclic order is a fragment automorphism.

    True  -> winding is notation, not structure; emit a fixed character.
    False -> winding is load-bearing (rac/meso, coordinated face); keep the geometry.
    None  -> undecidable; caller keeps today's geometric behavior.
    """
    if len(constituent_indices) < 3:
        # eta2 and below have no circulation; signed_circulation already returns '>'.
        return True

    flat = _orientation_symmetry_graph(mol)
    if flat is None:
        return None

    kind, order = _eta_traversal_order(mol, constituent_indices)
    if order is None:
        return None

    forward = _labeled_forms(flat, order, kind)
    reverse = _labeled_forms(flat, order[::-1], kind)
    if forward is None or reverse is None:
        return None
    return bool(forward & reverse)


@lru_cache(maxsize=1024)
def _winding_is_orientation_free(smiles, constituent_key):
    """Memoized `_orientation_free_from_mol` keyed on a fragment SMILES.

    `constituent_key` must be a tuple (hashable). Atom indices are the fragment's
    own SMILES atom order, matching `local_idx` / `constituent_indices`.
    """
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    except Exception:
        return None
    if mol is None or mol.GetNumAtoms() <= max(constituent_key, default=-1):
        return None
    return _orientation_free_from_mol(mol, list(constituent_key))


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

        # 1b. Close the shell on NON-binding chalcogen donors that carry a
        # valence deficit (a croconate/oxo ring O, a nitrito -O). The binding
        # loop above leaves every non-binding atom with default implicit-H
        # behaviour. When such an atom has a genuine deficit -- one bond where
        # oxygen wants two, and no H, charge or radical -- MolToSmiles emits it
        # BARE (`c(O)`, `ON=O`), because SetNoImplicit alone does not force a
        # bracket. The MetalloGen adapter's MolFromSmiles then re-adds the
        # implicit H, so the regenerated 3D structure gains a phantom hydrogen
        # the input never had -> round-trip atom-count mismatch (COLWIK
        # croconate 55->58, ACOXEX nitrito 75->77).
        #
        # Fix by charging the atom by its deficit: a valence-1 O becomes [O-]
        # (phenolate/alkoxide) -- 0 H, closed-shell, unambiguous in the string
        # AND embeddable. A neutral 0-H oxygen is a radical, which the adapter
        # drops back to bare O (so it stays phantom) and which UFF cannot type;
        # a formal charge survives the adapter and embeds. Restricted to O/S:
        # an aqua/hydroxo/carbonyl/ether O sits at full valence (no deficit)
        # and a real O-H shows a bonded H, so neither is ever touched. Binding
        # atoms are skipped -- their {slot} marker restores the metal bond on
        # parse. N is intentionally excluded: a bare non-binding N overlaps the
        # nitride/ammine notation ambiguity owned by oin/inline.py.
        binding_set = set(binding_indices_in_ligand)
        pt = Chem.GetPeriodicTable()
        try:
            rw_mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        for atom in rw_mol.GetAtoms():
            if atom.GetIdx() in binding_set:
                continue
            if atom.GetAtomicNum() not in (8, 16):  # O, S only
                continue
            if atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() != 0:
                continue
            if atom.GetTotalNumHs() != 0:
                continue
            if any(nb.GetAtomicNum() == 1 for nb in atom.GetNeighbors()):
                continue  # a real O-H / S-H keeps its hydrogen
            default_val = pt.GetValenceList(atom.GetAtomicNum())[0]
            deficit = default_val - atom.GetTotalValence()
            if deficit > 0:
                atom.SetNoImplicit(True)
                atom.SetFormalCharge(-deficit)

        # 2. Generate Canonical SMILES
        # isomericSmiles=True ensures we keep stereochem info if present
        #
        # h_faithful_smiles, not MolToSmiles: everything above this line works to get
        # each atom's hydrogen count right, and the writer can still throw it away. A
        # 0-H atom whose valence sits between two allowed values -- a 3-valent thiophene
        # sulfur, between sulfur's 2 and 4 -- serializes BARE, and a bare symbol
        # re-reads as "fill to the next allowed valence with H", so the count step 1
        # just froze comes back one too high. Step 1b above is a narrower, per-motif
        # version of the same repair. Default-OFF lever; see oin/hydrogen.py.
        kmol = rw_mol.GetMol()
        smiles = h_faithful_smiles(kmol, isomericSmiles=True, canonical=True)
        return smiles, kmol

    @staticmethod
    def canonical_donor_representative(ligand_mol, binder_idx):
        """Return the canonical representative of a binder's symmetric donor set.

        A monodentate carboxylate binds through one of two O's that are
        equivalent by resonance; in any single Kekule structure they are
        distinguishable (=O vs -O), so which O the 3D bond perception picked --
        and thus which carries the {slot} marker -- drifts between the two
        round-trip directions (ABAZIO: ``O{n}C(=O)`` vs ``OC(=O{n})``).
        Canonicalize which atom carries the slot:

        1. Find the binder's donor SET on a bond-order-*flattened* copy (every
           bond -> SINGLE, aromatic flags cleared): the two carboxyl O's become
           genuine graph automorphs and share a canonical-rank class. Nitro
           ({O,O}) and sulfonate/phosphonate ({O,O,O}) sets fall out for free;
           an ester's two inequivalent O's stay in separate classes (left
           alone -- no over-collapse).
        2. Pick the member with the lowest bond-order-*aware*
           ``CanonicalRankAtoms(breakTies=True)`` rank. Both round-trip mols
           carry the same fragment graph (only the binder *flag* differs, and
           the flag is not part of the graph), so this returns the same graph
           position in both -- a different physical atom, but the canonical one.

        Returns ``binder_idx`` unchanged on any failure, or when the binder has
        no symmetric partner (fail-safe: the emitted string is then
        byte-identical to today). Operates in ``ligand_mol``'s own atom-index
        space; call it before ``generate_robust_smiles`` so the force-bracket,
        the radical fill, and the {slot} marker all land on the same rep.
        """
        try:
            n_atoms = ligand_mol.GetNumAtoms()
            if not (0 <= binder_idx < n_atoms):
                return binder_idx

            # Bond-order-AWARE view (real C=O), for the final tie-break rank.
            # Mask KEKULIZE like _fragment_mol_for_canonicalization: forced
            # explicit-H binders can defeat the kekulizer, but aromaticity
            # perception (needed for a meaningful rank) does not require it.
            work = Chem.RWMol(ligand_mol)
            Chem.SanitizeMol(work, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)

            # Bond-order-AGNOSTIC view, for donor-set membership: flatten every
            # bond to SINGLE and drop aromatic flags so resonance-equivalent
            # donors collapse into one symmetry class.
            flat = Chem.RWMol(ligand_mol)
            for bond in flat.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
            for atom in flat.GetAtoms():
                atom.SetIsAromatic(False)
            Chem.SanitizeMol(flat, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)

            classes = list(Chem.CanonicalRankAtoms(flat, breakTies=False))
            donor_set = [i for i in range(n_atoms) if classes[i] == classes[binder_idx]]
            if len(donor_set) <= 1:
                return binder_idx

            ranks = list(Chem.CanonicalRankAtoms(work, breakTies=True))
            return min(donor_set, key=lambda idx: ranks[idx])
        except Exception:
            return binder_idx

    @staticmethod
    def canonical_eta_set_representative(ligand_mol, binder_indices):
        """Canonicalize WHICH of several equivalent rings carries an eta slot.

        The set-valued sibling of `canonical_donor_representative`. Tetraphenyl-
        borate binds a metal through one of its four interchangeable phenyls; the
        3D bond perception marks whichever phenyl physically faces the metal, so
        an input structure and its regenerated twin mark *different* rings and the
        OIN strings differ although the structure does not (SOJMIQ). Pick the
        canonical member of the ring's automorphic family instead.

        Returns ``{original_local_idx: canonical_local_idx}``, or ``{}`` when there
        is nothing to canonicalize or anything at all goes wrong (fail-safe: the
        emitted string is then byte-identical to today). Call it before
        `generate_robust_smiles`, so the forced brackets, the radical fill and the
        {slot} markers all land on the same ring -- otherwise brackets would mark
        the physical ring while {slot} marks the canonical one.

        Four guards, all required:

        1. the binders form a single closed ring. A fragment carrying two eta
           rings (ansa-metallocene) or a ring plus a pendant donor (a
           Cp-tethered phosphine) is not one orientable group and is left alone.
        2. the ring is orientation-free. Relabeling a ring whose winding is
           load-bearing could hide a genuine coordinated-face difference.
        3. at least one *other* ring is an automorphic image of it -- otherwise
           the marker is already on the only ring it could be on.
        4. no candidate overlaps the ring itself (fused systems).
        """
        try:
            binders = sorted(binder_indices)
            if len(binders) < 3 or len(set(binders)) != len(binders):
                return {}
            if max(binders) >= ligand_mol.GetNumAtoms():
                return {}

            kind, order = _eta_traversal_order(ligand_mol, binders)
            if kind != "ring":
                return {}  # guard 1
            if _orientation_free_from_mol(ligand_mol, binders) is not True:
                return {}  # guard 2

            flat = _orientation_symmetry_graph(ligand_mol)
            if flat is None:
                return {}

            def dihedral_signature(cycle):
                """(min canonical SMILES, the labeling that achieves it)."""
                n = len(cycle)
                best = None
                for base in (list(cycle), list(cycle)[::-1]):
                    for k in range(n):
                        labeling = base[k:] + base[:k]
                        canon = _labeled_canonical_smiles(flat, labeling)
                        if canon is None:
                            return None, None
                        if best is None or (canon, labeling) < best:
                            best = (canon, labeling)
                return best

            actual_sig, actual_labeling = dihedral_signature(order)
            if actual_sig is None:
                return {}

            actual_set = set(order)
            candidates = [(order, actual_labeling)]
            for ring in Chem.GetSymmSSSR(flat):
                ring = list(ring)
                if len(ring) != len(order) or set(ring) == actual_set:
                    continue
                if set(ring) & actual_set:
                    continue  # guard 4: fused, not an independent image
                ring_kind, ring_order = _eta_traversal_order(flat, ring)
                if ring_kind != "ring":
                    continue
                sig, labeling = dihedral_signature(ring_order)
                if sig == actual_sig:
                    candidates.append((ring_order, labeling))

            if len(candidates) < 2:
                return {}  # guard 3

            ranks = list(Chem.CanonicalRankAtoms(flat, breakTies=True))
            chosen_order, chosen_labeling = min(
                candidates, key=lambda entry: sorted(ranks[i] for i in entry[0])
            )
            if set(chosen_order) == actual_set:
                return {}  # already canonical -- emit exactly what we emit today

            # Both labelings realize the same canonical labeled graph, so mapping
            # them position-wise carries each marked atom onto its automorphic image.
            return dict(zip(actual_labeling, chosen_labeling))
        except Exception:
            return {}


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


# INTENTIONALLY SEPARATE from the generation-side generation/oin_parser.py::TEMPLATES.
# This is the *encoder*-side geometry table (positions + reference vectors used for
# alignment); the parser's TEMPLATES is positions-only and uses different slot
# conventions on some geometries (notably SPY), plus TET here is column-normalized
# below. Do NOT merge the two tables -- the encode and generate pipelines are
# deliberately independent. See the matching note in generation/oin_parser.py.
TEMPLATE_SPECS: dict[str, dict[int, dict]] = {
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
    # CN 8 -- square antiprismatic. `pos` mirrors MetalloGen's
    # `8_squre_antiprismatic` (globalvars.known_geometries_vector_dict), the same
    # /sqrt(3) vectors the generator places donors on, so encoder slot i and
    # generator slot i point the same direction (round-trip slot correspondence).
    # Two staggered squares: slots 0-3 top (+z), slots 4-7 bottom (-z), rotated
    # 45 deg. `ref` is vestigial for winding (only haptic groups use it, via the
    # actual ring centroid) but the schema requires a non-parallel vector.
    "SQA": {
        0: {"pos": [-0.5773503, 0.5773503, 0.5773503], "ref": [0, 0, 1]},
        1: {"pos": [0.5773503, 0.5773503, 0.5773503], "ref": [0, 0, 1]},
        2: {"pos": [-0.5773503, -0.5773503, 0.5773503], "ref": [0, 0, 1]},
        3: {"pos": [0.5773503, -0.5773503, 0.5773503], "ref": [0, 0, 1]},
        4: {"pos": [-0.8141210, 0.0, -0.5773503], "ref": [0, 0, 1]},
        5: {"pos": [0.0, -0.8141210, -0.5773503], "ref": [0, 0, 1]},
        6: {"pos": [0.8141210, 0.0, -0.5773503], "ref": [0, 0, 1]},
        7: {"pos": [0.0, 0.8141210, -0.5773503], "ref": [0, 0, 1]},
    },
    # CN 9 -- tricapped trigonal prismatic. `pos` is MetalloGen's
    # `9_tricapped_trigonal_prismatic` (globalvars.known_geometries_vector_dict)
    # unit-normalized, in the same slot order, so encoder slot i and generator
    # slot i point the same direction (round-trip slot correspondence, like SQA).
    # Slots 0-5 are the trigonal-prism vertices (two triangles on +z/-z), slots
    # 6-8 the three equatorial caps. `ref` is vestigial for winding (CN9 donors
    # are monodentate here; only haptic groups use ref, via the actual ring
    # centroid) but the schema requires a vector non-parallel to `pos` -- [0,0,1]
    # is non-parallel to every slot (max |z-component| 0.655).
    "TCT": {
        0: {"pos": [0.7555736, 0.0000000, 0.6550638], "ref": [0, 0, 1]},
        1: {"pos": [-0.3780710, 0.6546229, 0.6546229], "ref": [0, 0, 1]},
        2: {"pos": [-0.3780710, -0.6546229, 0.6546229], "ref": [0, 0, 1]},
        3: {"pos": [0.7555736, 0.0000000, -0.6550638], "ref": [0, 0, 1]},
        4: {"pos": [-0.3780710, 0.6546229, -0.6546229], "ref": [0, 0, 1]},
        5: {"pos": [-0.3780710, -0.6546229, -0.6546229], "ref": [0, 0, 1]},
        6: {"pos": [-1.0000000, 0.0000000, 0.0000000], "ref": [0, 0, 1]},
        7: {"pos": [0.5821873, -0.8130547, 0.0000000], "ref": [0, 0, 1]},
        8: {"pos": [0.5821873, 0.8130547, 0.0000000], "ref": [0, 0, 1]},
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


# Slack for the geometry-matcher's batched-numpy prefilter (see
# OINDiscreteAligner._candidate_permutations). The batched Kabsch rssd differs
# from scipy's per-permutation result by <1e-7 even on the worst-conditioned
# (ideal, symmetric) spheres, so 1e-6 keeps the candidate set a strict superset
# of the true scipy argmin while staying tight (~symmetry-order, not n!).
_MATCH_CANDIDATE_MARGIN = 1e-6


# --- Electronic (ligand-field) geometry prior --------------------------------------
# d-block group numbers, for d-electron count = group - oxidation state.
GROUP_NUMBER = {
    "Sc": 3,
    "Y": 3,
    "La": 3,
    "Ti": 4,
    "Zr": 4,
    "Hf": 4,
    "V": 5,
    "Nb": 5,
    "Ta": 5,
    "Cr": 6,
    "Mo": 6,
    "W": 6,
    "Mn": 7,
    "Tc": 7,
    "Re": 7,
    "Fe": 8,
    "Ru": 8,
    "Os": 8,
    "Co": 9,
    "Rh": 9,
    "Ir": 9,
    "Ni": 10,
    "Pd": 10,
    "Pt": 10,
    "Cu": 11,
    "Ag": 11,
    "Au": 11,
    "Zn": 12,
    "Cd": 12,
    "Hg": 12,
}

# Ambiguity band: the electronic prior fires only when the two competing low-CN
# templates BOTH fit poorly AND nearly-equally (the sphere is genuinely between
# them). Clean spheres (e.g. cisplatin: SPL 0.14 / TET 1.22) sit far outside this
# band, so a decisive geometric assignment is never overridden. Calibrated on
# cisplatin/BEPCAC (see tests/unit/test_geometry_electronic_prior.py).
_GEO_PRIOR_MIN_RMSD = 0.35  # the better-fitting template is still poor (distorted)
_GEO_PRIOR_TIE_TOL = 0.25  # the two templates fit within this RMSD of each other

# (coordination number, d-electron count) -> LFSE-preferred template, for cells with
# a strong, spin-independent electronic preference. d8 (Ni2+/Pd2+/Pt2+...) is square
# planar; d0/d10 (no LFSE preference) default to tetrahedral; d8 CN5 favours square
# pyramidal. Cells not listed leave the geometric winner untouched.
_ELECTRONIC_GEOMETRY_PRIOR = {
    (4, 8): "SPL",
    (4, 0): "TET",
    (4, 10): "TET",
    (5, 8): "SPY",
}


def metal_d_electron_count(mol):
    """d-electron count of the complex's metal, or None if not derivable.

    d = group_number - oxidation_state, with the oxidation state read from the metal
    atom's formal charge (``get_tmc_mol`` assigns the metal its formal oxidation
    state, e.g. Pt(II)/Ni(II) -> +2). Returns None for non-d-block centres or
    out-of-range counts, so the electronic prior simply does not fire.
    """
    for atom in mol.GetAtoms():
        group = GROUP_NUMBER.get(atom.GetSymbol())
        if group is None:
            continue
        dn = group - atom.GetFormalCharge()
        return dn if 0 <= dn <= 10 else None
    return None


def _electronic_geometry_override(n, current, rmsds_by_name, metal_dn):
    """LFSE-preferred template to override an *ambiguous* geometric match, else None.

    Overrides only when (a) a d-count is known, (b) the geometric winner and the
    electronic preference both fit poorly and nearly-equally (the ambiguity band),
    and (c) the (CN, d-count) cell has a strong preference. Otherwise trusts geometry.
    """
    if metal_dn is None:
        return None
    pref = _ELECTRONIC_GEOMETRY_PRIOR.get((n, metal_dn))
    if pref is None or pref == current:
        return None
    a = rmsds_by_name.get(current)
    b = rmsds_by_name.get(pref)
    if a is None or b is None:
        return None
    if min(a, b) > _GEO_PRIOR_MIN_RMSD and abs(a - b) < _GEO_PRIOR_TIE_TOL:
        return pref
    return None


class OINDiscreteAligner:
    """Assign ligand binding atoms to discrete coordination-geometry slots."""

    def __init__(self, metal_idx, ligands, metal_dn=None):
        """Initialize the aligner with the metal index and ligand data.

        Args:
            metal_idx: Index of the metal center atom.
            ligands: List of dicts, each of which MUST contain 'smiles'
                generated by OINSanitizer and a 'binding_atoms' list of
                tuples/lists: [(global_idx, mass, coords, local_idx)].
            metal_dn: d-electron count of the metal (from ``metal_d_electron_count``),
                or None. When set, an ambiguous low-CN geometric match is broken by the
                electronic (ligand-field) prior so the discrete geometry label is stable
                across distorted conformers. None disables the prior.
        """
        self.metal_idx = metal_idx
        self.ligands = ligands
        self.metal_dn = metal_dn

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
            best_res, rmsds_by_name = self._match_geometry_candidates(n_eff, virtual_atoms)
            if best_res:
                tmpl_name, tmpl_vectors, mapping, R_mat = best_res
                # Electronic tie-break: if the geometric match is ambiguous (distorted
                # sphere between two low-CN templates), let the metal's d-count decide,
                # so force-field conformer noise cannot flip the discrete label.
                override = _electronic_geometry_override(
                    n_eff, tmpl_name, rmsds_by_name, self.metal_dn
                )
                if override and override in TEMPLATES:
                    new_mapping, _rmsd, new_R = self._map_to_template(
                        virtual_atoms, TEMPLATES[override]
                    )
                    if new_mapping is not None:
                        tmpl_name, tmpl_vectors = override, TEMPLATES[override]
                        mapping, R_mat = new_mapping, new_R
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
        """Tightest-fitting template for ``n`` donors: ``(name, vectors, mapping, R)`` or None.

        Thin wrapper over ``_match_geometry_candidates`` -- returns only the best
        result, exactly as before the single-pass refactor.
        """
        best_result, _ = self._match_geometry_candidates(n, virtual_atoms)
        return best_result

    def _match_geometry_candidates(self, n, virtual_atoms):
        """Match ``virtual_atoms`` against every CN-``n`` candidate template in ONE pass.

        Returns ``(best_result, rmsds_by_name)`` where ``best_result`` is
        ``(name, vectors, mapping, R_mat)`` for the tightest-fitting template (or
        ``None``) -- byte-identical to what the pre-refactor ``_find_best_geometry_match``
        produced -- and ``rmsds_by_name`` maps each *evaluated* candidate template
        name to the RMSD its ``_map_to_template`` returned. Candidates skipped
        because the template has fewer slots than there are donors are absent from
        the map, mirroring ``coordination_geometry_fit``'s ``inf`` guard. This lets a
        caller read a target template's fit straight from the map without a second
        ``_map_to_template`` pass (the redundancy ``classify_and_fit`` removes).
        """
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
        elif n == 8:
            candidates = ["SQA"]
        elif n == 9:
            candidates = ["TCT"]
        else:
            # Fallback or robust handling
            if n > 9:
                candidates = ["OCT"]  # Best effort
            else:
                candidates = ["LIN"]

        min_rmsd = float("inf")
        best_result = None
        rmsds_by_name = {}

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
            rmsds_by_name[name] = rmsd
            logger.debug(f"  Candidate: {name}, RMSD: {rmsd:.4f}")

            if mapping is not None and rmsd < min_rmsd:
                min_rmsd = rmsd
                best_result = (name, vectors, mapping, R_mat)

        if best_result:
            logger.debug(f"  Selected: {best_result[0]} (RMSD {min_rmsd:.4f})")

        return best_result, rmsds_by_name

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

        # Exhaustively assign donors to slots, but score only the permutations that
        # can plausibly win. The per-permutation scipy Kabsch (Rotation.align_vectors)
        # dominates cost ~45x over the SVD math itself, so a vectorised batched-numpy
        # Kabsch first ranks every permutation and nominates the ones whose approximate
        # rssd is within _MATCH_CANDIDATE_MARGIN of the best (_candidate_permutations).
        # scipy then scores ONLY those, in lexicographic order, keeping the first strict
        # minimum -- exactly as the full sweep did. Because the margin dwarfs the
        # numpy-vs-scipy rssd gap, the true scipy argmin (and every permutation tied
        # with it) is provably among the candidates, so best_mapping / best_rmsd /
        # best_R are byte-identical to scoring all n! permutations.
        perms = list(itertools.permutations(range(n_slots), n_atoms))
        candidates = self._candidate_permutations(perms, template_vectors, input_norms)

        for cand_idx in candidates:
            slot_indices = perms[cand_idx]
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
    def _candidate_permutations(perms, template_vectors, input_norms):
        """Permutations worth scoring with scipy, via a batched-numpy Kabsch prefilter.

        Returns the indices into ``perms`` (ascending, i.e. lexicographic order) of
        every permutation whose approximate rssd is within ``_MATCH_CANDIDATE_MARGIN``
        of the best. The approximate rssd comes from a single vectorised Kabsch over
        all permutations at once; it differs from scipy's per-permutation result by
        <1e-7 even on the worst-conditioned (ideal, symmetric) coordination spheres,
        so a 1e-6 margin is a strict superset of the true scipy argmin and of any
        permutation tied with it. The caller re-scores the returned candidates with
        scipy for a byte-identical result -- this only prunes permutations that
        provably cannot win.

        ``template_vectors`` and ``input_norms`` are unit (or near-unit) vectors, so
        ``||a - R b||^2 = 2 - 2 a . R b`` and the minimised residual reduces to the
        singular values of the correlation matrix (the standard Kabsch identity).
        """
        idx = np.asarray(perms)  # (P, n_atoms)
        targets = template_vectors[idx]  # (P, n_atoms, 3)
        corr = np.einsum("pij,ik->pjk", targets, input_norms)
        u, s, vt = np.linalg.svd(corr)
        d = np.sign(np.linalg.det(u @ vt))
        alignment = s[:, 0] + s[:, 1] + d * s[:, 2]
        a2 = (targets**2).sum(axis=(1, 2))
        b2 = (input_norms**2).sum()
        approx = np.sqrt(np.maximum(a2 + b2 - 2.0 * alignment, 0.0))
        threshold = approx.min() + _MATCH_CANDIDATE_MARGIN
        return np.nonzero(approx <= threshold)[0]

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

    @staticmethod
    def _item_orientation_free(item):
        """True if this haptic group's winding is notation rather than structure.

        See the module header. False for anything that is not a haptic group, and
        for any group where the test is undecidable -- both keep today's geometric
        winding.
        """
        constituents = item.get("constituent_indices") or []
        if len(constituents) < 2:
            return False  # monodentate donor: no circulation to speak of
        return bool(_winding_is_orientation_free(item["chem_id"][1], tuple(sorted(constituents))))

    @classmethod
    def _topological_heading_atom(cls, smiles, constituent_indices):
        """Embedding-independent heading atom for an orientation-free eta group.

        Prefers the strict RC2 rank. Falls back to the valence-tolerant symmetry
        graph when strict sanitization fails -- a BPh4- borate is the case in
        practice, and today it silently drops through to the *geometric* heading,
        which is exactly what makes the marker wander between round-trip directions.
        Last resort is the lowest constituent index, which is still derived from
        canonical SMILES atom order and so is embedding-independent.
        """
        canonical_idx = cls._canonical_heading_atom(smiles, constituent_indices)
        if canonical_idx is not None:
            return canonical_idx

        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            flat = _orientation_symmetry_graph(mol) if mol is not None else None
            if flat is not None and flat.GetNumAtoms() > max(constituent_indices, default=-1):
                ranks = list(Chem.CanonicalRankAtoms(flat, breakTies=True))
                return min(constituent_indices, key=lambda idx: ranks[idx])
        except Exception:
            pass

        return min(constituent_indices, default=None)

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
                if self._item_orientation_free(item):
                    # A free ring's geometric winding depends on which face the
                    # embedding happened to present, so it is worthless as a
                    # *canonical* sort key. Collapse it and let the lowest
                    # constituent index below decide.
                    return ">"
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
            # Group by rank, but assign a heading to EACH haptic group (slot)
            # within a rank -- a single fragment can carry more than one eta
            # ring. A silane-bridged bis-indenyl ansa-metallocene is one
            # connected fragment (one rank) that occupies two eta slots; heading
            # (hence winding) must be assigned per haptic SLOT, not once per
            # rank, or the second ring loses its winding marker and rac/meso
            # diastereomers become indistinguishable in the OIN string. The
            # len(grp_coords) >= 2 guard below keeps this scoped to haptic
            # groups -- monodentate/polydentate single-atom donors are skipped
            # exactly as before, so non-eta ligands are unaffected.
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
            content_canonical_slots_handled = set()
            for rank, items in by_rank.items():
                for item in items:
                    smiles = item["chem_id"][1]
                    grp_coords = item.get("group_coords")
                    constituent_indices = item.get("constituent_indices", [])

                    if smiles in SYMMETRIC_LIGANDS:
                        continue
                    if grp_coords is None or len(grp_coords) < 2:
                        continue

                    # The heading of an eta group must never come from geometry, whether
                    # or not the ring is orientation-free. For a free ring the
                    # constituent indices may have been remapped onto a canonical
                    # automorphic ring (see
                    # OINSanitizer.canonical_eta_set_representative); for any ring the
                    # geometric pick tracks the embedding, so it lands on a different
                    # ring atom in an input structure than in its regenerated twin and
                    # the star wanders (GIPDEQ: a boron ylide whose fragment will not
                    # kekulize, so the strict rank is unavailable and the old code fell
                    # straight through to geometry). `_topological_heading_atom` has its
                    # own three tiers -- strict canonical rank, the valence-tolerant
                    # symmetry graph, then the lowest constituent index -- and every one
                    # of them is a function of the graph alone.
                    canonical_idx = self._topological_heading_atom(smiles, constituent_indices)
                    if canonical_idx is None:
                        continue

                    heading_local_indices.add((rank, canonical_idx))
                    content_canonical_slots_handled.add(item["slot"])

            for rank, items in by_rank.items():
                for item in items:
                    if item["slot"] in content_canonical_slots_handled:
                        continue

                    slot_idx = item["slot"]

                    # Skip if slot not in specs or no ref vector
                    if slot_idx not in template_spec or "ref" not in template_spec[slot_idx]:
                        continue

                    ref_vec = np.array(template_spec[slot_idx]["ref"])

                    # Get Group Coords
                    grp_coords = item.get("group_coords")
                    if grp_coords is None:
                        continue

                    # Calculate Centroid
                    centroid = np.mean(grp_coords, axis=0)

                    if len(grp_coords) < 2:
                        continue

                    best_dot = -float("inf")
                    best_idx = -1

                    ordered_indices = item["constituent_indices"]

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
            # For known symmetric ligands, force the heading atom to the lowest
            # local SMILES index (deterministic output, e.g. [CH2]{^}=[CH2]).
            # Applied per haptic group so an ansa-bis-Cp still marks both rings.
            for rank, items in by_rank.items():
                for item in items:
                    smiles = item["chem_id"][1]

                    if smiles not in SYMMETRIC_LIGANDS:
                        continue

                    grp_coords = item.get("group_coords")
                    if grp_coords is None or len(grp_coords) < 2:
                        continue

                    constituent = item["constituent_indices"]

                    # Find min index among this group's constituent_indices
                    ordered_indices = sorted(constituent)
                    forced_idx = ordered_indices[0]

                    # Remove any other heading picked for THIS group (e.g. by the
                    # geometric step above); match on the group's own constituent
                    # set so a sibling eta ring in the same rank keeps its heading.
                    to_remove = [
                        (r, idx)
                        for (r, idx) in heading_local_indices
                        if r == rank and idx in constituent
                    ]
                    for entry in to_remove:
                        heading_local_indices.discard(entry)

                    heading_local_indices.add((rank, forced_idx))
                    logger.debug(
                        f"Forced Heading Atom for Symmetric Ligand {smiles}: Index {forced_idx}"
                    )

        # 4c. Canonical winding across interchangeable automorphic eta groups (Lane 3).
        # Gated OFF by default: it rewrites the emitted winding character for the meso
        # arrangement of a symmetric ansa-metallocene, so it is not byte-identical.
        winding_override = {}
        if CANONICAL_ETA_WINDING and geometry_name in TEMPLATE_SPECS:
            winding_override = self._canonical_eta_winding(
                best_final_map,
                heading_local_indices,
                geometry_name,
                tmpl_vectors,
                alignment_rotation,
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

                    # An orientation-free ring (Cp*, an arene, a BPh4- phenyl) can
                    # be turned over onto itself by a proper rotation, so its
                    # geometric winding records only which face this particular
                    # embedding presented to the metal -- it flips at random between
                    # an input structure and a regenerated one. Emit the degenerate
                    # '>' instead, and the OIN string becomes a function of the
                    # structure alone. Rings that cannot be turned over (an
                    # ansa-bis(indenyl)'s rac vs meso) keep the geometric sign.
                    if (
                        geometry_name in TEMPLATE_SPECS
                        and slot in TEMPLATE_SPECS[geometry_name]
                        and not self._item_orientation_free(x)
                    ):
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
                        direction_char = winding_override.get((rank, slot, idx), direction_char)

                    tag += direction_char

                parts.append(tag)

        return ";".join(parts)

    @classmethod
    def _eta_automorphism_class(cls, smiles, constituent_indices):
        """Identity shared by eta groups related by a ligand-graph automorphism.

        The canonical SMILES of the fragment with every constituent atom carrying the
        SAME atom-map number: two eta groups get the same id exactly when some
        automorphism of the fragment graph carries one constituent set onto the other.
        Two eta rings on *different* fragments compare equal when the fragment strings
        are identical and the rings sit at corresponding positions, which is what makes
        an unbridged bis(indenyl) (two separate, identical fragments) come out as one
        class alongside a silane-bridged one (two rings inside a single fragment).

        Sanitization is deliberately skipped -- an unkekulizable borate or boron ylide
        must still get an id -- so this is a *string* identity, deterministic for a
        given fragment SMILES, which is all an equality test needs. Returns None when
        it cannot be computed, and a None id never joins any class.
        """
        if not smiles or not constituent_indices:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None or mol.GetNumAtoms() <= max(constituent_indices):
                return None
            rw = Chem.RWMol(mol)
            for atom in rw.GetAtoms():
                atom.SetAtomMapNum(0)
            for idx in constituent_indices:
                rw.GetAtomWithIdx(int(idx)).SetAtomMapNum(1)
            return Chem.MolToSmiles(rw.GetMol(), canonical=True)
        except Exception:
            return None

    def _canonical_eta_winding(
        self,
        best_final_map,
        heading_local_indices,
        geometry_name,
        tmpl_vectors,
        alignment_rotation,
    ):
        """Canonical winding assignment across interchangeable, automorphic eta groups.

        Returns ``{(rank, slot, heading_idx): winding_char}`` -- overrides only, empty
        when nothing needs canonicalizing.

        WHY THIS IS NEEDED
        ------------------
        Two OIN spellings denote the SAME compound exactly when they are related by a
        proper rotation of the coordination polyhedron composed with a graph automorphism
        of the ligand assembly. Both operations carry each ring's winding character along
        **unchanged**: a proper rotation preserves handedness, and an automorphism
        mapping one eta ring onto another preserves the cyclic sense (an automorphism
        that *reverses* it is precisely the orientation-free case, which is already
        forced to a fixed ``'>'`` and is excluded here). So the pair
        ``(ring identity, winding char)`` travels as a unit and the only freedom left is
        WHICH interchangeable slot each unit is written at. Today geometry decides that,
        which is why a symmetric ansa-metallocene's meso diastereomer encodes as
        ``{0<}...{1>}`` from one structure and ``{0>}...{1<}`` from its regenerated twin.

        WHY IT CANNOT DESTROY STEREOCHEMISTRY
        -------------------------------------
        No reflection is ever applied -- only proper rotations and graph automorphisms,
        the two operations that preserve molecular identity by definition. Concretely,
        on an ansa-metallocene: the **rac** diastereomer carries the same character on
        both rings, so sorting is a no-op and its mirror image (the other character on
        both rings) stays a different string; only **meso**, whose two mirror-related
        spellings really are one achiral compound, collapses. A Cp/fluorenyl pair is not
        automorphic, so each sits in a singleton orbit and nothing moves -- correctly,
        because for an unsymmetrical bridge the two spellings are genuine enantiomers.
        This is the guard the v0.4.4 axial wave lacked when it sorted a token by sign and
        silently made it reflection-invariant.
        """
        template_spec = TEMPLATE_SPECS[geometry_name]

        # 1. The eta groups whose winding is load-bearing, with their geometric sign.
        eta = {}
        for x in best_final_map:
            slot = x["slot"]
            cons = sorted(x.get("constituent_indices") or [])
            if len(cons) < 2 or slot not in template_spec:
                continue
            if self._item_orientation_free(x):
                continue  # winding is notation, already a fixed character
            heading = next((idx for idx in cons if (x["rank"], idx) in heading_local_indices), None)
            if heading is None or x.get("group_coords") is None:
                continue
            slot_def = template_spec[slot]
            eta[slot] = {
                "rank": x["rank"],
                "heading": heading,
                "cons": cons,
                "smiles": x["chem_id"][1],
                "class": self._eta_automorphism_class(x["chem_id"][1], cons),
                "winding": self._determine_winding(
                    grp_coords=x["group_coords"],
                    star_idx=heading,
                    constituent_indices=cons,
                    slot_z=np.array(slot_def["pos"]),
                    slot_x_ref=np.array(slot_def["ref"]),
                    alignment_rotation=alignment_rotation,
                ),
            }
        if len(eta) < 2:
            return {}

        # 2. Colour every occupied slot. Winding is deliberately NOT part of the colour:
        # we are asking which slots may be *relabelled*, and the answer must not depend
        # on the very characters we are about to reassign.
        colour = {}
        for x in best_final_map:
            slot = x["slot"]
            cons = sorted(x.get("constituent_indices") or [])
            colour[slot] = (x["chem_id"], self._eta_automorphism_class(x["chem_id"][1], cons))

        # 3. The proper rotations of this polyhedron that preserve every slot colour.
        # `_brute_force_symmetries` is built from `Rotation.from_euler`, so every element
        # is proper -- no reflection can enter here.
        allowed = [
            perm
            for perm in self._brute_force_symmetries(tmpl_vectors)
            if all(colour.get(perm[s]) == c for s, c in colour.items())
        ]
        if not allowed:
            return {}

        # 4. Orbits of eta slots under that group.
        orbit_of = {}
        for slot in eta:
            reachable = frozenset(perm[slot] for perm in allowed if perm[slot] in eta) | frozenset(
                [slot]
            )
            orbit_of[slot] = reachable
        orbits = {frozenset().union(*(orbit_of[s] for s in orbit_of[slot])) for slot in orbit_of}

        override = {}
        for orbit in orbits:
            members = sorted(orbit)
            # Scoped to a 2-orbit: that is every case the corpus presents (a metallocene's
            # two eta slots), and for a larger orbit the rotation group may realize only
            # SOME rearrangements, so the reachability argument below would need the full
            # induced-group computation. Fail safe to today's behaviour instead.
            if len(members) != 2:
                continue
            classes = {eta[s]["class"] for s in members}
            if len(classes) != 1 or None in classes:
                continue  # not automorphic -- these spellings are NOT interchangeable

            s0, s1 = members
            a, b = eta[s0], eta[s1]
            eps = self._eta_swap_sense(a, b)
            if eps is None:
                continue

            # Is this arrangement ACHIRAL -- i.e. does the mirror spelling denote the same
            # compound? Reflecting the structure flips BOTH characters. Applying the
            # slot-swapping proper rotation together with the ligand automorphism maps
            # (w0, w1) to (eps*w1, eps*w0), so that composite equals the mirror spelling
            # (-w0, -w1) exactly when:
            #     eps == +1  and  w0 != w1      (two separate but identical fragments:
            #                                    each ring is canonically ordered the same
            #                                    way, so the swap preserves cyclic sense)
            #     eps == -1  and  w0 == w1      (two rings inside ONE bridged fragment:
            #                                    the canonical SMILES traverses out along
            #                                    one ring and back along the other, so the
            #                                    swap reverses cyclic sense)
            # When that holds the mirror is reachable by proper operations alone, the
            # compound is achiral, and its two spellings must collapse to one. When it
            # does NOT hold the two spellings are genuine ENANTIOMERS and folding them
            # would destroy exactly the stereochemistry the winding marker exists to carry
            # -- an ansa-metallocene's rac/meso pair. Both branches were checked against
            # the independent geometric oracle in `tools/eta_core_chirality.py`
            # (agreement 4/4); the naive "sort the two characters" rule this replaces got
            # every bridged case backwards.
            same = a["winding"] == b["winding"]
            achiral = (eps == 1 and not same) or (eps == -1 and same)
            if not achiral:
                continue

            flip = {">": "<", "<": ">"}
            current = (a["winding"], b["winding"])
            mirrored = (flip[a["winding"]], flip[b["winding"]])
            if mirrored < current:  # '<' precedes '>': deterministic canonical choice
                for slot, char in zip(members, mirrored):
                    override[(eta[slot]["rank"], slot, eta[slot]["heading"])] = char
        return override

    @classmethod
    def _eta_swap_sense(cls, a, b):
        """+1/-1: does swapping these two eta groups preserve or reverse cyclic sense?

        The winding character is read off each ring's OWN ascending-index (canonical
        SMILES) order, so an automorphism carrying one ring onto the other does not
        simply hand the character over -- it hands it over *possibly reversed*, and which
        it is decides whether the achiral arrangement is the same-sign or the
        opposite-sign one. Getting this backwards silently makes the encoder
        reflection-invariant, which is the v0.4.4 axial wave's failure exactly.

        Two eta groups on *different* fragments with the same fragment SMILES are two
        copies of one ligand, canonically ordered identically, so the swap is
        sense-preserving: +1. Two groups inside ONE fragment need the real automorphism.
        Returns None when the groups are not automorphic at all (a Cp/fluorenyl ansa
        pair), which means their spellings are not interchangeable and nothing may move.
        """
        if a["smiles"] != b["smiles"] or not a["smiles"]:
            return None
        if a["rank"] != b["rank"]:
            return 1  # two separate, identical fragments

        A, B = list(a["cons"]), list(b["cons"])
        if len(A) != len(B):
            return None
        try:
            mol = Chem.MolFromSmiles(a["smiles"], sanitize=False)
            if mol is None or mol.GetNumAtoms() <= max(A + B):
                return None
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
            n = len(B)
            rev = B[::-1]
            for match in mol.GetSubstructMatches(
                mol, uniquify=False, useChirality=False, maxMatches=20000
            ):
                if len(match) != mol.GetNumAtoms():
                    continue
                image = [match[i] for i in A]
                if set(image) != set(B):
                    continue
                for k in range(n):
                    if image == B[k:] + B[:k]:
                        return 1
                    if image == rev[k:] + rev[:k]:
                        return -1
        except Exception:
            return None
        return None

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

        # 2. Winding axis = this ring's ACTUAL metal->centroid direction, in the
        # molecular frame. `grp_coords` is metal-centered (see `_reduce_hapticity`:
        # 'group_coords': grp_coords - metal_origin), so the centroid of its atoms
        # IS the metal->centroid (outward) vector. Using the real per-ring axis --
        # rather than the idealized template slot position `slot_z` -- keeps the
        # winding sign robust when the coordination sphere is distorted: an ansa-
        # metallocene's bridged bite angle squeezes its two eta rings well inside
        # the ideal tetrahedral slot directions, so a single global
        # alignment_rotation of the template axis lands on the wrong side of the
        # (second) ring plane and flips its sign. The actual centroid axis can't.
        # It also matches the generation-side convention (haptic-face winding is
        # measured against the actual ring centroid).
        grp = np.asarray(grp_coords, dtype=float)
        axis_mol = grp.mean(axis=0)

        if np.linalg.norm(axis_mol) < 1e-9:
            # Degenerate (ring centroid coincident with metal): fall back to the
            # template slot direction so the sign is still deterministic.
            slot_z_norm = np.linalg.norm(slot_z)
            if slot_z_norm < 1e-9:
                raise ValueError(
                    "slot_z must be a nonzero, outward-facing (metal->centroid) vector"
                )
            axis_template = slot_z / slot_z_norm
            axis_mol = (
                alignment_rotation.inv().apply(axis_template)
                if alignment_rotation is not None
                else axis_template
            )

        return signed_circulation(grp_coords, star_local_idx, axis_mol)

    def _brute_force_symmetries(self, vectors):
        """Proper-rotation vertex permutations of this coordination template.

        v0.4.5 Lane 2: delegates to ``oin.canonical_slots.derive_rotation_group``, now the
        single source of truth for the polyhedron rotation groups -- the encoder here, the
        comparison key in ``compare.py`` and the canonical-slot post-pass all read one
        table and one group derivation (open debt TD-005).

        This used to brute-force Euler triples from the fixed grid
        ``[0, 90, 120, 180, 240, 270]``. That grid **cannot express a 72-degree five-fold
        rotation**, so on PBP it found 2 of the 10 proper rotations and the encoder could
        not canonicalize a pentagonal-bipyramidal equatorial labeling at all. Measured
        against the derived group it agreed on the other 10 geometries and never invented
        a non-rotation, so unifying is a no-op everywhere except PBP, where it is a fix.
        ``tests/unit/test_canonical_slots.py`` pins both halves of that statement.

        The name is kept because it is what the call site and the tests already say;
        nothing is brute-forced any more.
        """
        key = tuple(tuple(np.round(np.asarray(v, dtype=float), 6)) for v in vectors)
        return _cached_rotation_group(key)


@lru_cache(maxsize=64)
def _cached_rotation_group(vectors_key):
    """``derive_rotation_group`` memoized on a rounded vertex tuple.

    ``_permute_and_serialize`` asks for the group on every molecule, always for one of the
    eleven fixed templates, so the derivation runs once per geometry per process.
    """
    return derive_rotation_group([list(v) for v in vectors_key])


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
    with ``n == 8`` mapping to ``"SQA"`` (square antiprismatic) and ``n == 9`` to
    ``"TCT"`` (tricapped trigonal prismatic); for ``n > 9`` the matcher falls back
    to ``"OCT"`` as a best effort (which then fails the slot-count check, yielding
    ``None``) -- so
    callers that need an eta/haptic guard must gate on the expected coordination
    number themselves rather than relying on a ``None`` return.
    """
    virtual_atoms = [{"coords": np.asarray(v, dtype=float)} for v in donor_vectors]
    aligner = OINDiscreteAligner(0, [])  # ligands unused by the matcher
    result = aligner._find_best_geometry_match(len(virtual_atoms), virtual_atoms)
    return result[0] if result else None


def coordination_geometry_fit(donor_vectors, geo_code):
    """RMSD of the best assignment of ``donor_vectors`` to a geometry's ideal template.

    Uses the template for ``geo_code`` (from ``TEMPLATES``); returns
    ``float('inf')`` if the code is unknown or the template has fewer slots than
    there are donors. Lower is a tighter fit to the ideal geometry. Classification
    alone only tells
    you which template a coordination sphere is *closest* to -- a heavily puckered
    square-plane can still be labelled ``"SPL"`` because it beats ``"TET"``/``"TPY"``
    -- so this score lets callers pick the *cleanest* conformer of a target
    geometry, not merely the first one that classifies as it. Uses the same
    rotation-invariant Kabsch/permutation matcher as the encoder, so scores are
    comparable across conformers.
    """
    template = TEMPLATES.get(geo_code)
    if template is None or len(donor_vectors) > len(template):
        return float("inf")
    virtual_atoms = [{"coords": np.asarray(v, dtype=float)} for v in donor_vectors]
    aligner = OINDiscreteAligner(0, [])
    mapping, rmsd, _ = aligner._map_to_template(virtual_atoms, template)
    if mapping is None:
        return float("inf")
    return rmsd


def classify_and_fit(donor_vectors, geo_code):
    """Best-match OIN code and the fit RMSD of ``donor_vectors`` to ``geo_code``, in ONE pass.

    Single-pass equivalent of ``classify_coordination_geometry(donor_vectors)``
    followed by ``coordination_geometry_fit(donor_vectors, geo_code)``. A single
    ``_match_geometry_candidates`` call yields both the classified label and, from
    the same candidate loop, the RMSD of ``donor_vectors`` against ``geo_code``'s
    template -- so when the two are consumed together (the geometry-aware conformer
    selector), the identical permutation/Kabsch search runs once instead of twice.

    Returns ``(label, fit)``:

    * ``label`` -- the tightest-fitting template's OIN code (``None`` if none matched),
      exactly ``classify_coordination_geometry(donor_vectors)``.
    * ``fit`` -- the RMSD to ``geo_code``'s template. When ``geo_code`` is among the
      coordination-number's candidates it comes free from the match; otherwise it
      falls back to ``coordination_geometry_fit`` (``float('inf')`` for an unknown
      code or a template with fewer slots than there are donors). Byte-identical to
      the old two-call sequence -- it is the same computation.
    """
    virtual_atoms = [{"coords": np.asarray(v, dtype=float)} for v in donor_vectors]
    aligner = OINDiscreteAligner(0, [])
    best_result, rmsds_by_name = aligner._match_geometry_candidates(
        len(virtual_atoms), virtual_atoms
    )
    label = best_result[0] if best_result else None
    if geo_code in rmsds_by_name:
        return label, rmsds_by_name[geo_code]
    return label, coordination_geometry_fit(donor_vectors, geo_code)
