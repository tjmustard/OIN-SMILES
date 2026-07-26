"""MetalloGen 3D-generation backend adapter.

Bridges OIN-SMILES (``ParsedOIN``) to the vendored MetalloGen engine in
``oinsmiles.generator3d`` (dummy-metal + RDKit ``CoordMap`` embed + constrained
MMFF/UFF cleanup). This is the 3D-generation backend used by ``OIN3DGenerator``.

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

from ..generator3d import clash, generate_3d_structures, get_xyz_string, globalvars, om
from ..oin.axial import mol_axial_token, parse_axial_token
from ..oin.compare import canonical_roundtrip_key
from . import _telemetry
from .oin_parser import OINParser, ParsedOIN
from .structure import GeneratedStructure

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
    "default": {
        "ff_max_iters": 200,
        "ff_force_tol": 1e-4,
        "ff_energy_tol": 1e-6,
        "d_converge": 0.05,
    },
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
            raise ValueError(f"Unknown ff_preset {preset!r}; choose from {sorted(FF_PRESETS)}")
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
    # Value MUST byte-match the vendored globalvars key, which is MISSPELLED
    # "squre" (see generator3d/globalvars.known_geometries_vector_dict). Do not
    # "fix" it -- the lookup KeyErrors otherwise. Sibling key "8_sqaure_prismatic"
    # is likewise misspelled upstream.
    "SQA": "8_squre_antiprismatic",
    # CN 9 -- tricapped trigonal prismatic (e.g. Y/Ln with 3 bidentate + 1
    # tridentate donors). Value byte-matches the vendored globalvars key.
    "TCT": "9_tricapped_trigonal_prismatic",
}


class UncoordinatedFragmentError(ValueError):
    """An OIN fragment carries no binding slot.

    Outer-sphere counterions and uncoordinated solvent (a free water, a borate
    anion) are emitted as fragments by the encoder but have no bond to the
    metal. MetalloGen's ``metal|lig1|...|geo`` m-SMILES has no way to express
    such a fragment, so the structure cannot be generated.

    Subclasses ``ValueError`` so existing callers that catch ``ValueError``
    keep working.
    """


def _prepare_ligand_fragments(parsed: ParsedOIN):
    """Shared per-fragment preparation for the m-SMILES and OIN-direct paths (SL2).

    Returns ``(metal_frag, ligand_specs, geo)`` where ``ligand_specs`` is a list of
    ``(mapped_smiles, winding)`` ordered exactly as the m-SMILES join (ascending
    first coordination slot). Each ``mapped_smiles`` is the canonical per-fragment
    string the m-SMILES path emits: kekulized (Cp), bare-donor hydrogen reconciled,
    each binding atom tagged with the atom-map number of the nearest MetalloGen
    coordinate slot (isomerism-preserving nearest-vector match).

    ``winding`` is ``None`` unless the fragment is an eta ring whose heading atom
    carries a ``>``/``<`` marker, in which case it is
    ``(output_order, star_frag_idx, char)``: ``output_order`` is RDKit's
    ``_smilesAtomOutputOrder`` from the canonicalizing ``MolToSmiles``
    (``output_order[canonical_pos] = original_fragment_index``), ``star_frag_idx``
    is the heading atom's original fragment index, and ``char`` is the requested
    winding. This lets ``om.get_om_from_parsed`` recover the encoder's ring order
    and star atom after canonicalization.

    Raises:
        UncoordinatedFragmentError: a ligand fragment has no binding slot.
    """
    geo = OIN_TO_METALLOGEN_GEO.get(parsed.geo_code, "")
    if not geo:
        raise ValueError(f"Geometry code '{parsed.geo_code}' not supported by MetalloGen mapping.")

    metallogen_vectors = globalvars.known_geometries_vector_dict[geo]
    num_slots = len(metallogen_vectors)
    specs_by_slot: list = [None] * num_slots

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
                pass  # If it still fails, let it be

        frag_vectors = [v for v in parsed.vectors if v.fragment_idx == i]
        if not frag_vectors:
            raise UncoordinatedFragmentError(
                f"Fragment {i} ({frag_smiles!r}) has no binding slot; uncoordinated "
                "fragments are not representable in MetalloGen m-SMILES."
            )

        for v in frag_vectors:
            target_vec = np.array(v.vector)
            dists = np.linalg.norm(metallogen_vectors - target_vec, axis=1)
            mg_slot_idx = int(np.argmin(dists))

            # Reconcile the binding atom's hydrogen count with its donor role.
            #
            # An EXPLICIT H count from the OIN (a bracket atom -- [NH], [OH2],
            # [CH2] -- has NoImplicit set) is authoritative: the encoder already
            # decided this is a neutral dative L-type donor that KEEPS its H
            # (secondary amine, aqua, sigma-alkyl/benzyl). Only a BARE binding
            # atom carries a phantom implicit H that the metal bond replaces, so
            # only bare atoms are reinterpreted below.
            #
            # The bracket/bare split is decided by ``replace_map`` in
            # ``oin/inline.py``: it de-brackets a binding atom only when the
            # bracket content is a bare organic-subset symbol (C, N, O, n, ...).
            # An H-bearing donor therefore always serializes bracketed -- the
            # sole exception is ammine NH3, force-de-bracketed to ``N{n}``.
            atom = mol.GetAtomWithIdx(v.atom_in_fragment_idx)
            if not atom.GetNoImplicit():
                sym = atom.GetSymbol()
                heavy = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
                strip = False
                if sym == "C":
                    if any(b.GetBondType() == Chem.BondType.TRIPLE for b in atom.GetBonds()):
                        strip = True  # C#O / C#N carbon has no H
                    elif atom.GetIsAromatic():
                        # sigma-aryl carbanion strips its H; a haptic ring carbon
                        # (Cp/arene, >2 coordinating atoms in the ring) keeps it.
                        ring_info = mol.GetRingInfo()
                        coordinating_indices = [vec.atom_in_fragment_idx for vec in frag_vectors]
                        is_haptic = any(
                            atom.GetIdx() in ring
                            and sum(1 for idx in ring if idx in coordinating_indices) > 2
                            for ring in ring_info.AtomRings()
                        )
                        strip = not is_haptic
                        # Lock a bare 0-H haptic carbon (ipso / ring-fusion /
                        # substituted) so the anionic-aromatic re-parse in
                        # get_ligand_from_smiles cannot re-derive a phantom
                        # implicit H on it. The bare `c{n}` these fragments use
                        # for an eta ipso/fusion carbon has no explicit H, so
                        # AddHs(explicitOnly=False) protonates it to `[cH]`
                        # (ARONEA +4, BOXJUU +6). Freezing it at its
                        # correctly-perceived 0-H count reproduces the input.
                        # A genuine C-H haptic carbon keeps H via its explicit
                        # `[cH]` bracket (NoImplicit already set -> outer guard
                        # skips it), so Cp/arene passers stay byte-identical.
                        if is_haptic and atom.GetTotalNumHs() == 0:
                            atom.SetNumExplicitHs(0)
                            atom.SetNoImplicit(True)
                    elif heavy >= 2:
                        # Non-aromatic sigma-carbon donor already bonded to >=2 heavy
                        # atoms: an NHC carbene (C between two ring N) or a carbanion.
                        # It cannot carry implicit H AND a metal bond without exceeding
                        # valence, so it is a 0-H donor -- drop the phantom CH2 hydrogens.
                        strip = True
                elif sym in ("O", "S"):
                    # Bare chalcogen donor = anionic alkoxide / thiolate / oxo -> 0 H.
                    # (A dative aqua/hydroxo/alcohol keeps its H via the explicit branch.)
                    strip = True
                elif sym == "N":
                    # A bare N with any heavy neighbour is a 0-H anionic X-type
                    # donor: amido, anilide, silylamide, azide, phosphinimide
                    # (P=N). Every N donor that keeps a hydrogen serializes as
                    # [NH2]/[NH] and took the explicit branch above; the only
                    # bare H-bearing N is ammine NH3, which has no heavy
                    # neighbour. Hence `heavy >= 1` -- exact, not a heuristic.
                    #
                    # A bare heavy==0 N is genuinely ambiguous: nitride [N] and
                    # ammine [NH3] both serialize to `N{n}`. Resolving it needs
                    # an OIN format change (oin/inline.py), so it is left alone
                    # here and the ammine reading wins.
                    strip = heavy >= 1
                if strip:
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(0)

            # MetalloGen map numbers are 1-based (slot index + 1).
            atom.SetAtomMapNum(mg_slot_idx + 1)

        mapped_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        # RDKit records the canonical output ordering as a side effect of
        # MolToSmiles: output_order[canonical_pos] = original fragment atom index.
        # Needed to map the heading (star) atom and ring order back through
        # canonicalization for deterministic winding construction.
        output_order = list(
            mol.GetPropsAsDict(includePrivate=True, includeComputed=True)["_smilesAtomOutputOrder"]
        )
        # Place the fragment at its first binding slot (monodentate = its only slot;
        # multidentate carries all its map numbers within this one fragment string).
        first_slot = int(
            np.argmin(np.linalg.norm(metallogen_vectors - np.array(frag_vectors[0].vector), axis=1))
        )

        # Eta-ring winding target: the heading atom is the sole ring atom carrying a
        # >/< marker (OINVector.winding); its absence means no winding to construct.
        winding = None
        heading = next((v for v in frag_vectors if v.winding is not None), None)
        if heading is not None:
            winding = (output_order, heading.atom_in_fragment_idx, heading.winding)

        specs_by_slot[first_slot] = (mapped_smiles, winding)

    ligand_specs = [s for s in specs_by_slot if s is not None]
    return metal_frag, ligand_specs, geo


def convert_parsed_to_msmiles(parsed: ParsedOIN) -> str:
    """Convert a ``ParsedOIN`` into MetalloGen's ``metal|lig1|...|geo`` m-SMILES.

    Each ligand binding atom is tagged with the atom-map number of the nearest
    MetalloGen coordinate slot (isomerism-preserving nearest-vector match). Thin
    wrapper over ``_prepare_ligand_fragments`` -- byte-identical to the pre-SL2
    output (winding metadata is discarded here; the OIN-direct path consumes it).

    Raises:
        UncoordinatedFragmentError: a ligand fragment has no binding slot.
    """
    metal_frag, ligand_specs, geo = _prepare_ligand_fragments(parsed)
    msmiles_parts = [metal_frag] + [mapped_smiles for mapped_smiles, _ in ligand_specs]
    return "|".join(msmiles_parts) + f"|{geo}"


def convert_oin_to_msmiles(oin_string: str) -> str:
    """Parse an OIN string and convert it to MetalloGen m-SMILES."""
    return convert_parsed_to_msmiles(OINParser().parse(oin_string))


# Atom property carrying a binding atom's OIN coordination slot (-1 = does not bind the
# metal) through the bond-order-transfer substructure match, and the mol property carrying
# a template's OIN fragment index.
_SLOT_PROP = "_oinSlot"
_FRAG_IDX_PROP = "_oinFragIdx"
_NON_DONOR = -1

# Bond orders that a wrong automorphism can mislocalize, and how strongly a candidate
# match is penalized for placing them on a long bond (a triple is shorter than a double).
_ORDER_WEIGHT = {Chem.BondType.DOUBLE: 1.0, Chem.BondType.TRIPLE: 2.0}

# Cap on candidate matches enumerated when disambiguating bond-order transfer. Across the
# dataset's ligand templates the automorphism count has median 4 and p90 72, and the eta
# ligands this exists for top out around 48 -- so the cap only ever bites on ligands whose
# automorphisms come from freely-permuting substituents (a tBu/aryl-substituted porphyrin
# reaches 3e4+). Truncating there is safe: every enumerated candidate already satisfies the
# donor constraint, which is more than the old unconstrained single match did, and RDKit's
# enumeration order is deterministic so the pick stays reproducible. Keep it small -- the
# porphyrins are 85-atom graphs and enumeration, not scoring, is what costs (an 85-atom
# template runs ~4.5 ms per call here against ~800 ms at MATCH_MAX = 2000).
MATCH_MAX = 512

# How much worse (in Angstrom, summed over a template's multiple bonds) the legacy
# unconstrained match may score before it is treated as wrong and replaced. Misplacing one
# C=C moves it from ~1.34 A onto a ~1.50 A single bond, so a genuine error costs >= 0.15
# and COD's two double bonds cost >= 0.30. Anything under this is symmetry-equivalent maps
# separated by force-field noise, where re-picking would churn the re-encoded OIN string.
SCORE_TOL = 0.05

# Above this many binding atoms in ONE fragment, skip the exact donor->slot assignment
# (its cost is 2**n) and leave the fragment to the unanchored fallback. A fragment's donor
# count is a coordination number, so this is never reached in practice.
_ASSIGN_MAX_DONORS = 12


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
            # sanitize=False: RemoveHs' internal sanitize KEKULIZES rings it cannot
            # leave aromatic and strips their aromatic flags. For a neutral-radical
            # Cp (e.g. TiCat2, whose SanitizeMol above already failed) that corrupts
            # an aromatic c1cccc1 into an aliphatic C1=CC=CC1, so the re-encoded OIN
            # string mismatches the aromatic input. Skipping the sanitize preserves
            # the aromatic perception; it still removes the explicit H atoms, so the
            # heavy-atom count the contract-mol match relies on is unchanged. No-op
            # for fragments that sanitized cleanly above.
            t = Chem.RemoveHs(t, sanitize=False)
        except Exception:
            pass
        # Label CIP so build_contract_mol can carry encoded sp3 stereo (the
        # template's @/@@ is the ground truth vs the stochastic embed handedness).
        try:
            Chem.AssignStereochemistry(t, cleanIt=True, force=True)
        except Exception:
            pass
        # Unparseable fragments are skipped above, so a template's position in the
        # returned list is NOT its OIN fragment index. build_contract_mol needs the
        # real index to look up the fragment's coordination vectors.
        t.SetIntProp(_FRAG_IDX_PROP, k)
        templates.append(t)
    return templates


def _template_donor_slots(t, parsed: ParsedOIN):
    """``(atom_idx -> slot colour, slot -> unit vector)`` for a template's binding atoms.

    The colour is the OIN coordination slot, which is what separates a chelate's donors
    from one another: a porphyrin's four N sit in four distinct slots, so colouring by
    slot forbids the macrocycle rotations that a bare donor/non-donor colour allows.
    Atoms of one haptic group (eta2-alkene, eta3-allyl, Cp) share a slot, hence a slot.

    Degrades to a single colour 0 -- i.e. plain donor/non-donor -- when the OIN slots or
    their vectors are unusable, so a fragment is never left unconstrained.
    """
    if not t.HasProp(_FRAG_IDX_PROP):
        return {}, {}
    fk = t.GetIntProp(_FRAG_IDX_PROP)
    vecs = [v for v in parsed.vectors if v.fragment_idx == fk]
    if not vecs:
        return {}, {}
    binary: tuple[dict, dict] = ({v.atom_in_fragment_idx: 0 for v in vecs}, {})
    if any(v.slot < 0 for v in vecs):
        return binary
    slots = {v.atom_in_fragment_idx: v.slot for v in vecs}
    unit = {}
    for v in vecs:
        if v.slot not in unit:
            a = np.asarray(v.vector, dtype=float)
            n = float(np.linalg.norm(a))
            if n > 0:
                unit[v.slot] = a / n
    # Every slot needs a direction, else a generated donor cannot be assigned to one.
    return (slots, unit) if set(unit) == set(slots.values()) else binary


def _generated_donor_slots(donors_local, tslots, slot_unit, carr, metal_idx):
    """Assign each generated binding atom the OIN slot it occupies, as a whole.

    Not a per-atom nearest-vector lookup: a haptic group straddles its slot vector over a
    wide arc, so an eta3-allyl terminus can point nearer a *neighbouring* slot than its own
    (FIKXIJ, whose allyl and phosphine chelate one metal). Solving it globally -- assign the
    fragment's donors to its slots, with each slot taking exactly as many donors as the
    template gives it, maximizing total alignment -- keeps every donor in the right group.

    Exact via a bitmask DP; the donor count per fragment is a coordination number, so the
    2**n stays tiny. Returns None if it cannot be assigned.
    """
    n = len(donors_local)
    if not slot_unit:
        return {qi: 0 for qi, _ in donors_local}  # degraded: donor/non-donor only
    slots = sorted(tslots.values())
    if n == 0 or len(slots) != n or n > _ASSIGN_MAX_DONORS:
        return None

    cost = []
    for _, gidx in donors_local:
        u = carr[gidx] - carr[metal_idx]
        norm = float(np.linalg.norm(u))
        if norm == 0.0:
            return None
        u = u / norm
        cost.append([-float(np.dot(u, slot_unit[s])) for s in slots])

    size = 1 << n
    inf = float("inf")
    dp = [inf] * size
    pick = [-1] * size
    dp[0] = 0.0
    for mask in range(size):
        if dp[mask] == inf:
            continue
        i = bin(mask).count("1")  # donors 0..i-1 are placed
        if i == n:
            continue
        for j in range(n):
            if mask >> j & 1:
                continue
            nxt = mask | (1 << j)
            c = dp[mask] + cost[i][j]
            if c < dp[nxt]:
                dp[nxt] = c
                pick[nxt] = j
    out = {}
    mask = size - 1
    for i in range(n - 1, -1, -1):
        j = pick[mask]
        if j < 0:
            return None
        out[donors_local[i][0]] = slots[j]
        mask ^= 1 << j
    return out


def _flatten_template(t, donor_slots=None):
    """Connectivity-only copy of a template (all single bonds, no aromatic/charge).

    Used only for substructure matching, so bond orders can be transferred even
    from templates that never sanitize (e.g. C#O, O valence 3).

    Radical electrons / explicit-H valence must also be cleared: an OIN ligand
    atom that binds the (now-stripped) metal is under-valent, so its template
    atom carries a radical (e.g. the three metal-bound carbons of an eta3-allyl,
    ``[CH2][CH]=[CH]...``). RDKit's substructure matcher treats radical-electron
    count as a match constraint, so a radical-bearing query never matches the
    generated fragment (whose atoms are H-saturated with 0 radicals) -- the
    match silently returns empty and no bond orders/aromaticity transfer, which
    dearomatizes the ligand in the round trip. Normalizing to a plain
    connectivity graph (0 radicals, implicit H) fixes the match without ever
    loosening it (heavy-atom count + connectivity + element still gate it).

    ``donor_slots`` stamps ``_oinSlot`` on EVERY atom (the binding atom's OIN slot, or
    -1), so the match can be restricted to maps that send a template atom binding slot
    *s* onto a generated atom binding that same slot. RDKit's ``atomProperties`` treats
    an ABSENT property as a non-match, so the stamp must cover every atom of query and
    target alike, not just the donors.
    """
    ft = Chem.RWMol(t)
    for b in ft.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in ft.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNumRadicalElectrons(0)
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
        if donor_slots is not None:
            a.SetIntProp(_SLOT_PROP, donor_slots.get(a.GetIdx(), _NON_DONOR))
    m = ft.GetMol()
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return m


def _localizable_bonds(t) -> list:
    """``(begin, end, weight)`` for each non-aromatic double/triple bond of a template.

    Only such bonds can be *mislocalized* by picking the wrong automorphism: an
    all-single or purely aromatic template transfers identically under every map.
    """
    out = []
    for b in t.GetBonds():
        weight = _ORDER_WEIGHT.get(b.GetBondType())
        if weight is not None and not b.GetIsAromatic():
            out.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), weight))
    return out


def _transfer_score(bonds, match, q2g, dmat) -> float:
    """Total generated bond length lying under the template's multiple bonds.

    Every candidate match is an isomorphism onto the same generated fragment, so a
    given template bond's candidate images all join the same pair of elements --
    the sums are directly comparable across matches. Minimizing puts double/triple
    bonds on the SHORT edges, which is where MetalloGen actually embedded them
    (it built the geometry from an m-SMILES that still had the true bond orders).
    Aromatic bonds are excluded: they are ~1.39 A whichever way they map.
    """
    return sum(w * dmat[q2g[match[i]], q2g[match[j]]] for i, j, w in bonds)


def _slot_valid(flat, qmol, match) -> bool:
    """True if ``match`` sends every template atom onto one binding the same OIN slot."""
    return all(
        flat.GetAtomWithIdx(i).GetIntProp(_SLOT_PROP)
        == qmol.GetAtomWithIdx(match[i]).GetIntProp(_SLOT_PROP)
        for i in range(flat.GetNumAtoms())
    )


def _select_match(qmol, flat, t, q2g, dmat):
    """Donor-anchored, geometry-scored substructure match of ``flat`` into ``qmol``.

    ``qmol`` (the generated fragment) and ``flat`` are both heavy-atom, all-single
    connectivity graphs with equal atom counts, so a plain ``GetSubstructMatch`` is an
    automorphism search that returns an ARBITRARY map. For a symmetric eta ligand most
    of those maps move the bond orders: 1,5-COD's flattened 8-ring has 16 automorphisms
    and only 4 keep the C=C on the metal-bound carbons, which is why generated COD
    re-encodes with its double bonds two positions around the ring.

    Two constraints identify a good map. The ``_oinSlot`` colour forbids sending a
    template atom that binds slot *s* onto a generated atom that does not. That kills the
    COD ring rotations, the ring swap that moved PIJCAO's alkyne onto a para-ethyl, and
    the macrocycle rotations of a porphyrin's four distinct N slots. It cannot separate
    the two ends of a symmetric eta3-allyl -- all three carbons share one slot, and both
    maps put the C=C on opposite sides -- so candidates are ranked by ``_transfer_score``
    against the embedded 3D geometry, where the true C=C is ~0.1 A shorter.

    The legacy unconstrained match is PREFERRED when it is slot-valid and scores within
    ``SCORE_TOL`` of the best candidate. Bond orders are not the only thing the caller
    transfers through this map -- formal charges and ``_CIPCode`` stereo ride along -- so
    among equally-good maps the choice is arbitrary but observable downstream. Repairing
    only demonstrably-wrong maps keeps this a strict repair rather than a re-pick.
    """
    params = Chem.SubstructMatchParameters()
    params.uniquify = False
    params.atomProperties = [_SLOT_PROP]
    # An all-single/aromatic template transfers identically under every automorphism, so
    # stop at the first slot-valid map. This is also what keeps the common case cheap:
    # most big automorphism groups belong to tBu/phenyl-rich phosphines with no
    # localizable bond, and they never reach the enumeration below.
    bonds = _localizable_bonds(t)
    params.maxMatches = MATCH_MAX if bonds else 1
    candidates = qmol.GetSubstructMatches(flat, params)
    if not candidates:
        return ()
    if len(candidates) == 1 or not bonds:
        return candidates[0]

    def score(m):
        return _transfer_score(bonds, m, q2g, dmat)

    best = min(candidates, key=score)
    legacy = qmol.GetSubstructMatch(flat)
    if (
        legacy
        and len(legacy) == flat.GetNumAtoms()
        and _slot_valid(flat, qmol, legacy)
        and score(legacy) <= score(best) + SCORE_TOL
    ):
        return legacy
    return best


def _template_lp_label(t, ai: int) -> "str | None":
    """Return the lone-pair CIP ('R'/'S') for a Zone-A P donor template atom.

    Deliberately uses ``rdCIPLabeler`` -- NOT the legacy ``Chem.AssignStereochemistry``
    that ``_oin_fragment_templates`` stamps as ``_CIPCode`` for the backbone paths.
    The two labelers disagree for a 3-coordinate P (ACUWUT: legacy 'R' vs
    rdCIPLabeler 'S'), and ``ChiralityRecoveryUtility.recover``'s Zone-A lone-pair
    branch recomputes with ``rdCIPLabeler`` on the metal-free fragment. The template
    is likewise metal-absent and 3-coordinate -- the same configuration recover()
    sees -- so an rdCIPLabeler label taken here round-trips against recover()'s own
    recomputation. Returns None if the label can't be assigned.
    """
    # Same aromatic-preserving fresh re-parse as _template_sp3_label: a P donor
    # bonded to an aromatic ring (GUXPIA) has a representation-sensitive rdCIPLabeler
    # label, and recover()'s lone-pair branch reads it on the same re-parse.
    return _template_sp3_label(t, ai)


def _template_sp3_label(t, ai: int) -> "str | None":
    """Aromatic-preserving rdCIPLabeler label for a specified sp3 C/Si/S template atom.

    rdCIPLabeler gives OPPOSITE R/S for a stereocentre bonded to an aromatic haptic
    ring (Cp, indenyl, fluorenyl) depending on whether that ring is left aromatic vs
    charged/kekulized -- and the metal-free fragment ``recover()`` re-orients against
    emits the ring aromatic. So the label is taken with kekulization SKIPPED
    (``SANITIZE_ALL ^ SANITIZE_KEKULIZE``) to match the convention recover() reads.

    The label is computed on a FRESH re-parse of the template's SMILES, not the
    ``_oin_fragment_templates`` object directly: that object's aromatic state is
    corrupted by its ``RemoveHs(sanitize=False)`` for a FUSED haptic ring (indenyl,
    fluorenyl), which flips rdCIPLabeler for the adjacent centre (KAGXUM) even under
    the aromatic-preserving sanitize. An atom-map probe survives the re-parse so the
    target atom is found regardless of SMILES atom re-ordering. The label is read
    through the SAME shared reparse helper (``_reparse_cip_label_once``, fill-first)
    that ``ChiralityRecoveryUtility.recover`` uses, so the stamp and the fragment
    comparison agree on the donor-normalised convention (a metal-adjacent carbene/
    alkene/oxo donor's open valence is H-filled identically on both sides -- ORIHUU/
    XILZID/JEKQAS). Returns None on failure.
    """
    from ..core.chirality import _reparse_cip_label_once

    _PROBE = 99
    try:
        tagged = Chem.Mol(t)
        tagged.GetAtomWithIdx(ai).SetAtomMapNum(_PROBE)
        smiles = Chem.MolToSmiles(tagged)
    except Exception:
        return None
    for fill_deficit in (True, False):
        label = _reparse_cip_label_once(smiles, _PROBE, fill_deficit)
        if label is not None:
            return label
    return None


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

    from ..core.chirality import _LP_CIP_PROP, _SP3_CIP_PROP
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
        mapping: list = []
        frag_mols = Chem.GetMolFrags(
            frag_rw, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping
        )

        donor_set = set(donors)
        # Precomputed once: _select_match scores up to MATCH_MAX candidates per fragment
        # and would otherwise recompute the same bond lengths for every one of them.
        carr = np.asarray(coords, dtype=float)
        dmat = np.linalg.norm(carr[:, None, :] - carr[None, :, :], axis=-1)
        templates = _oin_fragment_templates(parsed)
        template_slots = [_template_donor_slots(t, parsed) for t in templates]
        flats = [_flatten_template(t, s) for t, (s, _) in zip(templates, template_slots)]
        used = [False] * len(templates)
        # global contract-atom idx -> encoded CIP code for sp3 stereocentres that
        # recover() leaves untouched (carbon, silicon, sulfur); oriented by the
        # perceive-then-flip loop below. Backbone phosphorus is handled separately
        # via an _OIN_CIPCode stamp (see the template loop).
        sp3_stereo_targets: dict[int, str] = {}
        # global contract-atom idx -> rdCIPLabeler lone-pair CIP for a Zone-A P
        # *donor* (metal-bonded, stereogenic lone pair). Stamped as _OIN_CIPCode_LP
        # and tag-seeded AFTER 3D perception below so recover()'s lone-pair
        # verify-and-flip keeps it -- the same end state a forward-pass CIPAssigner
        # mol reaches. Kept separate from sp3_stereo_targets because the label is
        # the rdCIPLabeler convention (recover() recomputes with rdCIPLabeler),
        # not the legacy _CIPCode the backbone paths carry.
        zone_a_lp_targets: dict[int, str] = {}
        # global contract-atom idx of every non-N/P sp3 centre the PARSED OIN
        # specified with an @ (a specified template chiral tag). A superset of
        # sp3_stereo_targets (which additionally requires a legacy _CIPCode).
        # Used below to clear geometry-invented tags on centres the OIN left
        # UNspecified, so AssignStereochemistryFrom3D cannot fabricate an @ the
        # forward encoder never emitted.
        oin_specified_sp3: set[int] = set()

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
            donors_local = [(qi, g) for qi, g in enumerate(q2g) if g in donor_set]
            for qi in range(qmol.GetNumAtoms()):
                qmol.GetAtomWithIdx(qi).SetIntProp(_SLOT_PROP, _NON_DONOR)

            # Pass 1 anchors the match on each binding atom's OIN slot and ranks the
            # survivors against the 3D geometry. Pass 2 is the pre-existing unanchored
            # match, reached only when no template survives pass 1 (e.g. MetalloGen bonded
            # a haptic ring only partially, so the donor counts disagree) -- so this is
            # never worse than the arbitrary-automorphism behaviour it replaces.
            chosen = None
            for anchored in (True, False):
                for ti, t in enumerate(templates):
                    if used[ti] or t.GetNumAtoms() != qmol.GetNumAtoms():
                        continue
                    if anchored:
                        tslots, slot_unit = template_slots[ti]
                        # Equal heavy-atom counts are not enough to pair two ligands: a
                        # sigma-ethyl would otherwise consume an eta2-ethylene template.
                        if len(tslots) != len(donors_local):
                            continue
                        gslots = _generated_donor_slots(
                            donors_local, tslots, slot_unit, carr, metal_idx
                        )
                        if gslots is None:
                            continue
                        for qi, slot in gslots.items():
                            qmol.GetAtomWithIdx(qi).SetIntProp(_SLOT_PROP, slot)
                    # Connectivity-only match, then copy real bond orders / aromaticity /
                    # charge from the template (works even when the template can't sanitize).
                    match = (
                        _select_match(qmol, flats[ti], t, q2g, dmat)
                        if anchored
                        else qmol.GetSubstructMatch(flats[ti])
                    )
                    if not match or len(match) != t.GetNumAtoms():
                        continue
                    chosen = (ti, t, match)
                    break
                if chosen:
                    break
            if chosen:
                ti, t, match = chosen
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
                    # Carry encoded sp3 stereo from the template. Backbone carbon
                    # (and Si/S, which recover() leaves untouched) is oriented by
                    # the perceive-then-flip loop below. A backbone phosphorus is an
                    # exception: ChiralityRecoveryUtility.recover() would clear its
                    # tag as a stray (its 4-neighbour, no-_OIN_CIPCode else-branch),
                    # so instead stamp the encoded CIP as _OIN_CIPCode -- recover()
                    # then keeps and orients it, exactly as for a forward-pass mol
                    # that went through CIPAssigner.
                    gidx = q2g[match[ai]]
                    anum = ta.GetAtomicNum()
                    if anum in (6, 14, 16) and (
                        ta.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    ):
                        oin_specified_sp3.add(gidx)
                        # Stamp the metal-free template's rdCIPLabeler label so
                        # recover() re-orients this centre on the metal-free
                        # fragment (the metal-present flip loop below mis-orients
                        # a metal-/eta-adjacent centre whose CIP label flips when
                        # the metal is stripped). Aromatic-preserving so a centre
                        # bonded to a haptic Cp/indenyl ring is labelled in the
                        # same convention recover() reads.
                        sp3_label = _template_sp3_label(t, ai)
                        if sp3_label is not None:
                            rwa.SetProp(_SP3_CIP_PROP, sp3_label)
                    elif (
                        anum == 7
                        and ta.GetTotalDegree() == 4
                        and ta.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    ):
                        # Genuine quaternary ammonium N+: a real tetrahedral
                        # stereocentre. Route it through the same metal-free
                        # _SP3_CIP_PROP re-orientation recover() applies to C/Si/S
                        # (POYJIX). Without a stamp the generated N carries no CIP
                        # prop, and recover()'s 4-neighbour no-_OIN_CIPCode fallback
                        # clears it ([N@@+] -> [N+]). Degree-4 gate keeps a trivalent
                        # amine N (RDKit-cleared inversion) unstamped and deferred.
                        n_label = _template_sp3_label(t, ai)
                        if n_label is not None:
                            rwa.SetProp(_SP3_CIP_PROP, n_label)
                    if anum == 15 and gidx in donors:
                        # Zone-A P donor: stereogenic lone pair. recover()'s
                        # lone-pair branch needs _OIN_CIPCode_LP (rdCIPLabeler
                        # convention) + a seeded tag, both applied after 3D
                        # perception below. Gate on the parsed [P@] chiral tag --
                        # reliable even when legacy AssignStereochemistry declines
                        # to CIP-label a 3-coordinate P -- and take the label from
                        # rdCIPLabeler, NOT the legacy _CIPCode (they disagree).
                        if ta.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            label = _template_lp_label(t, ai)
                            if label is not None:
                                zone_a_lp_targets[gidx] = label
                    elif ta.HasProp("_CIPCode"):
                        if anum in (6, 14, 16):
                            sp3_stereo_targets[gidx] = ta.GetProp("_CIPCode")
                        elif anum == 15 and gidx not in donors:
                            rwa.SetProp("_OIN_CIPCode", ta.GetProp("_CIPCode"))
                used[ti] = True

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

        # Carry ENCODED sp3 stereo (carbon, silicon, sulfur). The embed picks a
        # random handedness at backbone stereocentres, so 3D-perceived stereo can
        # be the enantiomer of what the OIN fragment SMILES encodes. Where the
        # geometry-derived CIP disagrees with the template CIP, flip the tag -- the
        # perceive-then-flip pattern ChiralityRecoveryUtility already uses for
        # Zone-A P (a bare tag flip does not survive get_oin_string's fragment
        # rebuild; a CIP-validated one does). Bounded fixed point: flipping one
        # centre can change another's CIP priority ranking. (Backbone P is oriented
        # downstream by recover() via its _OIN_CIPCode stamp, not here.)
        if sp3_stereo_targets:
            _CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
            _CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            for _ in range(3):
                try:
                    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                except Exception:
                    break
                changed = False
                for gidx, want in sp3_stereo_targets.items():
                    a = mol.GetAtomWithIdx(gidx)
                    cur = a.GetPropsAsDict().get("_CIPCode")
                    tag = a.GetChiralTag()
                    if cur and cur != want and tag in (_CW, _CCW):
                        a.SetChiralTag(_CCW if tag == _CW else _CW)
                        changed = True
                if not changed:
                    break

        # Honor the OIN as the source of truth for which sp3 centres are
        # stereo-specified. AssignStereochemistryFrom3D above stamps a tag on
        # every tetrahedral-looking centre in the embed geometry, including ones
        # the OIN left UNspecified -- a non-stereogenic carbon frozen into a
        # chiral-looking conformation, or a centre the fully-sanitised forward
        # encoder (get_tmc_mol) drops where this lenient contract mol keeps it.
        # That invents an @ the input never emitted (KAPCEM: 0 in -> 4 out), so
        # the round trip fails on a centre that was never specified. Clear the
        # geometry-derived tag on any non-N/P sp3 centre the parsed template did
        # not specify. N/P are oriented by recover(); Zone-A P donors are seeded
        # just below; specified centres (kept, and oriented by the flip loop) are
        # a superset of sp3_stereo_targets so this never undoes the carry above.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in (7, 15):
                continue
            if atom.GetChiralTag() not in (
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            ):
                continue
            if atom.GetIdx() not in oin_specified_sp3:
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # Zone-A P donor lone-pair stereo. Make the generated donor arrive at
        # get_oin_string in the same state a forward-pass CIPAssigner mol would:
        # the rdCIPLabeler lone-pair label stamped as _OIN_CIPCode_LP plus a
        # specified chiral tag. recover()'s lone-pair branch then verifies and
        # flips to match the stored label on the metal-free fragment. Must run
        # AFTER 3D perception and the sp3 flip loop above, both of which clear or
        # overwrite chiral tags (the dative metal->P bond makes
        # AssignStereochemistryFrom3D return CHI_UNSPECIFIED here anyway). The
        # seeded handedness is arbitrary -- recover() recomputes and reorients --
        # it only needs to be non-UNSPECIFIED for the flip to have a target.
        for gidx, label in zone_a_lp_targets.items():
            p = mol.GetAtomWithIdx(gidx)
            p.SetProp(_LP_CIP_PROP, label)
            p.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        # Drop E/Z from double bonds a metal-containing ring holds rigid, exactly as
        # the XYZ->OIN convert path does before its own get_oin_string call
        # (translator.XYZToSMILES.convert). AssignStereochemistryFrom3D above stamps a
        # directional marker on every localized ring double bond -- including a
        # porphyrin's meso C=C/C=N bridges, which chelate the metal and are physically
        # ring-locked. The forward encode strips those (they have no free E/Z), but the
        # fast re-encode here feeds this mol straight to get_oin_string, so without the
        # same clear the generated OIN carries slashes the input never had and the
        # round trip fails on a bond that was never stereogenic.
        from ..core.translator import _clear_chelate_locked_bond_stereo

        _clear_chelate_locked_bond_stereo(mol)
        return mol
    except Exception:
        return None


# Number of energy-ranked conformers to consider for geometry-aware selection.
# Matches generator3d's default FF pool size, so with an optimizer set (which
# already MACE-optimizes the whole pool) this adds no extra optimizer cost -- it
# only changes *which* of the already-optimized conformers is returned.
DEFAULT_SELECT_POOL = 5

# Wider conformer pool when the OIN encodes eta-ring winding. A ring's
# coordinated face (hence an ansa-metallocene's rac/meso diastereomer and its
# enantiomer) is set stochastically by each embed, so more candidates are needed
# to reliably sample the requested winding for winding-aware selection to pick.
# Re-encoding is cheap (contract-mol fast path), so the extra candidates cost
# mainly embedding, not perception.
ETA_SELECT_POOL = 16


def _norm_geo_code(code):
    """Fold equivalent square-planar codes together (``SQP`` == ``SPL``)."""
    return "SPL" if code == "SQP" else code


def _expected_coordination_number(geo_code):
    """Donor count implied by an OIN geo code.

    Taken from the MetalloGen name prefix (e.g. ``SPL`` -> ``4_square_planar``
    -> 4), or None if unknown.
    """
    name = OIN_TO_METALLOGEN_GEO.get(geo_code)
    if not name:
        return None
    try:
        return int(name.split("_", 1)[0])
    except (ValueError, IndexError):
        return None


_HAPTIC_GROUP_CUTOFF = 1.6  # A -- same threshold as oin_aligner._reduce_hapticity


def _reduce_haptic_positions(donor_positions, expected_n):
    """Cluster raw donor positions into haptic groups (one centroid per group).

    Returns one centroid per group only if the number of groups equals
    ``expected_n``. Mirrors the XYZ->OIN encoder's hapticity reduction
    (``oin_aligner.OINDiscreteAligner._reduce_hapticity``): binding atoms within
    ``_HAPTIC_GROUP_CUTOFF`` of each other (transitively, so a whole Cp ring is one
    group) collapse to a single coordination point at their centroid. This lets an
    eta complex whose raw donor count (e.g. TiCat's 12 ring C + 2 methyl) would gate
    it out of geometry-aware selection reduce to its true coordination number (4
    here: 2 Cp centroids + 2 methyl) so it can be classified/ranked like a discrete
    geometry. Returns ``None`` when the group count doesn't match ``expected_n`` --
    keeping selection strictly non-regressive for genuinely non-matching spheres.
    """
    n = len(donor_positions)
    visited = set()
    groups = []
    for j in range(n):
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
            for k in range(n):
                if k in visited:
                    continue
                if (
                    np.linalg.norm(donor_positions[curr] - donor_positions[k])
                    < _HAPTIC_GROUP_CUTOFF
                ):
                    stack.append(k)
        groups.append(component)
    if len(groups) != expected_n:
        return None
    return [np.mean([donor_positions[k] for k in g], axis=0) for g in groups]


def _coordination_vectors(contract_mol, expected_n=None):
    """Metal-centered donor vectors (``donor_pos - metal_pos``) for a contract mol.

    Returns None when no metal is found. When ``expected_n`` is given and the
    metal's raw donor count doesn't match it, first attempt hapticity reduction
    (``_reduce_haptic_positions``) -- an eta ligand bonds the metal to many ring
    atoms that collapse to a single coordination point, so a bent metallocene's 14
    raw donors reduce to 4 centroid donors. If reduction still can't reach
    ``expected_n`` the sphere has no discrete template and None is returned.
    """
    from ..utils.xyz2mol import TRANSITION_METALS_NUM

    metal_idx = next(
        (
            i
            for i in range(contract_mol.GetNumAtoms())
            if contract_mol.GetAtomWithIdx(i).GetAtomicNum() in TRANSITION_METALS_NUM
        ),
        None,
    )
    if metal_idx is None:
        return None
    donors = [
        b.GetOtherAtomIdx(metal_idx) for b in contract_mol.GetAtomWithIdx(metal_idx).GetBonds()
    ]
    conf = contract_mol.GetConformer()
    m = conf.GetAtomPosition(metal_idx)
    metal_pos = np.array([m.x, m.y, m.z])
    donor_positions = [
        np.array(
            [
                conf.GetAtomPosition(d).x,
                conf.GetAtomPosition(d).y,
                conf.GetAtomPosition(d).z,
            ]
        )
        for d in donors
    ]
    if expected_n is not None and len(donor_positions) != expected_n:
        # More raw donors than sites => try collapsing eta groups to centroids.
        if len(donor_positions) > expected_n:
            reduced = _reduce_haptic_positions(donor_positions, expected_n)
            if reduced is not None:
                return [tuple(c - metal_pos) for c in reduced]
        return None
    return [tuple(dp - metal_pos) for dp in donor_positions]


def _perceive_geo_code(contract_mol, expected_n=None):
    """Best-matching OIN geometry code for a contract mol's coordination sphere.

    Returns None when perception is not possible (see ``_coordination_vectors``
    for the guards).
    """
    from ..utils.oin_aligner import classify_coordination_geometry

    vecs = _coordination_vectors(contract_mol, expected_n)
    if vecs is None:
        return None
    return classify_coordination_geometry(vecs)


# Winding heading marker on an eta-ring atom: {n>} (clockwise) or {n<} (ccw).
_ETA_WINDING_RE = re.compile(r"\{\d+([<>])\}")


def _eta_winding_multiset(oin_string):
    """Sorted list of eta-ring winding characters in an OIN string.

    An ansa-metallocene's rac/meso (and a ring's coordinated face / enantiomer)
    is captured by the multiset of per-ring winding markers -- e.g. a rac
    bis-indenyl is ``['>', '>']`` while its meso diastereomer is ``['<', '>']``.
    Using the multiset (not per-slot) makes the comparison robust to which
    identical ring is assigned to which slot on re-encoding.
    """
    if not oin_string:
        return []
    return sorted(_ETA_WINDING_RE.findall(oin_string))


def _reencode_oin_fast(contract_mol):
    """Fast XYZ->OIN re-encode of an already-perceived contract mol.

    Skips the expensive bond perception in ``get_tmc_mol`` by feeding the
    contract mol (metal + dative bonds + bond orders + 3D conformer, from
    ``build_contract_mol``) straight into ``get_oin_string`` -- the same
    aligner/serializer the full encoder uses downstream of perception. Verified
    to yield the same eta winding as the full ``XYZToSMILES().convert`` path,
    which makes a wide winding-selection pool affordable. Returns None on
    failure (caller falls back to the full XYZ re-encode).
    """
    from ..utils.xyz2mol import get_oin_string

    if contract_mol is None:
        return None
    try:
        conf = contract_mol.GetConformer()
        coords = np.array(
            [
                [
                    conf.GetAtomPosition(a).x,
                    conf.GetAtomPosition(a).y,
                    conf.GetAtomPosition(a).z,
                ]
                for a in range(contract_mol.GetNumAtoms())
            ]
        )
        return get_oin_string(contract_mol, coords)
    except Exception:
        logger.debug("fast winding re-encode failed for a conformer", exc_info=True)
        return None


def _reencode_oin(mol):
    """Re-encode a generated conformer's 3D structure back to an OIN string.

    Full-fidelity fallback for ``_reencode_oin_fast``: uses the same XYZ->OIN
    path the round-trip verification uses (write XYZ -> ``XYZToSMILES().convert``),
    so the winding it reports is exactly what the round trip will compare
    against. Returns None on any failure.
    """
    import os
    import tempfile

    from ..core.translator import XYZToSMILES

    tmp_path = None
    try:
        xyz = get_xyz_string(mol)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as tmp_file:
            tmp_file.write(xyz)
            tmp_path = tmp_file.name
        return XYZToSMILES().convert(tmp_path)
    except Exception:
        logger.debug("winding re-encode failed for a conformer", exc_info=True)
        return None
    finally:
        if tmp_path is not None:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


def _reencode_key_matches(parsed, m, target_key, cmol=None, require_no_stretch=False, cache=None):
    """True when conformer ``m`` INDEPENDENTLY re-encodes to ``target_key`` (SL1 accept stamp).

    The honest round-trip acceptance test for the generate-until-key-exact early-exit
    (``OIN_EARLY_EXIT`` / ``ff_params["early_exit"]``). ``target_key`` is
    ``canonical_roundtrip_key(input_oin)`` (SL0's fac/mer key). Two-stage:

    1. **Cheap pre-filter** -- ``_reencode_oin_fast(cmol)`` re-serializes the *generated
       geometry* through the generator's own contract-mol connectivity. Its key reflects the
       generated donor arrangement (fac/mer, winding), so a MISMATCH here is a reliable
       "geometry is wrong" signal -> reject without paying the full path. A fast MATCH is NOT
       trusted for acceptance: it shares the generator's connectivity and is blind to
       bond-order / connectivity regressions (circular), so we confirm below.
    2. **Independent confirm** -- ``_reencode_oin(m)`` writes XYZ and re-perceives everything
       via ``XYZToSMILES().convert``. This is the accept decision.

    When ``require_no_stretch`` (i.e. ``clash.STRETCHED_BOND_ENABLED``) is set, a conformer that
    round-trips but carries any stretched bond is rejected. Returns ``False`` on any failure.

    ``cache``: optional dict, keyed on ``id(m)``, memoizing the expensive step-2 re-encode
    (measured 48-57s per call on an eta/haptic conformer -- a full ``XYZToSMILES().convert()``
    round trip). This function is called on the SAME mol object from two sites that can both
    run within one ``MetalloGenAdapter.generate()`` call: the pool-fill loop's ``accept_fn``
    (``generator3d/__init__.py``) and ``_select_by_geometry_impl``'s own early-exit re-scan
    below. Every mol reaching the second site already had this exact predicate evaluated once
    at the first (a mol only skips that first test when it comes from the untested
    ``stereo_rejects`` fallback, or when ``accept_fn`` itself was never constructed -- both
    leave no cache entry, so the lookup simply misses and this recomputes exactly as before).
    ``m``'s geometry is immutable between the two sites, so a cache hit reproduces the prior
    result byte-for-byte; this changes only which of two call sites pays the cost, never the
    verdict. Default ``None`` -> no cache is consulted or written -> identical to pristine for
    any caller that does not pass one.
    """
    try:
        if cmol is None:
            cmol = build_contract_mol(parsed, m)
        fast = _reencode_oin_fast(cmol)
        if fast is not None and canonical_roundtrip_key(fast) != target_key:
            return False
        if cache is not None and id(m) in cache:
            full = cache[id(m)]
        else:
            full = _reencode_oin(m)
            if cache is not None:
                cache[id(m)] = full
        if full is None or canonical_roundtrip_key(full) != target_key:
            return False
        if require_no_stretch and clash.mol_stretched_bond_count(m) > 0:
            return False
    except Exception:
        logger.debug("accept-first re-encode/key comparison failed for a conformer", exc_info=True)
        return False
    return True


def _axial_narrow(candidates, to_mol, target_axial):
    """Keep candidates whose axial token equals ``target_axial``.

    Returns ``(kept, n_blind)`` where ``n_blind`` counts candidates whose token could not be
    PERCEIVED at all (``mol_axial_token`` -> ``None``) as opposed to perceived-and-different.

    That distinction is the point. A genuine miss (every conformer perceived, none matching)
    and a perception failure (nothing perceived, so nothing can match) both produce an empty
    result and both fall through to the unfiltered pool -- silently returning a plausible
    structure of the WRONG atropisomer. This is exactly how the axial pass was defeated
    during development: raw MetalloGen conformers are unsanitized, ``CanonicalRankAtoms``
    raised, every token came back ``None``, and the filter compared ``None`` against ``'-'``.
    Counting the blind ones lets the caller say which failure it was.
    """
    kept = []
    n_blind = 0
    for c in candidates:
        mol = to_mol(c)
        token = mol_axial_token(mol) if mol is not None else None
        if token is None:
            n_blind += 1
        elif token == target_axial:
            kept.append(c)
    return kept, n_blind


def _axial_report_miss(n_candidates, n_blind, target_axial, *, key):
    """Record an axial narrowing miss, loudly when it was really a PERCEPTION failure."""
    if n_candidates and n_blind == n_candidates:
        # Nothing was perceivable -- the filter never had a chance. This is a defect in the
        # axial path itself, not a legitimately absent conformer, so warn rather than
        # whisper: the run will otherwise return the wrong atropisomer and still look fine.
        logger.warning(
            "axial-aware selection: axial perception FAILED on all %d candidates "
            "(requested %r); falling through -- the generated atropisomer may be wrong",
            n_candidates,
            target_axial,
        )
        _telemetry.record("adapter.axial_perception_failed", **{key: n_candidates})
        return
    logger.debug(
        "axial-aware selection: no candidate matched token %r (%d candidates, %d unperceived)",
        target_axial,
        n_candidates,
        n_blind,
    )
    _telemetry.record("adapter.axial_miss", n_blind=n_blind, **{key: n_candidates})


def _select_by_geometry_impl(
    parsed, mols, honor_winding=True, early_exit=False, reencode_cache=None
):
    """Choose the conformer that best realizes the requested coordination geometry.

    Falls back to the lowest-energy conformer. Three levels of preference over
    plain energy ranking:
      1. Only conformers whose coordination sphere *classifies* as the target OIN
         code are eligible (a distorted geometry that best matches a different
         template is rejected).
      2. When the OIN encodes eta-ring winding (an ansa-metallocene's rac/meso, a
         ring's coordinated face/enantiomer), prefer -- among the geometry-
         eligible conformers, best geometry first -- the one whose re-encoded
         winding multiset matches the request. The embed produces the ring face
         stochastically, so without this the diastereomer/enantiomer is left to
         chance; matching it here is what makes the eta round trip reproducible.
      3. Otherwise pick the tightest fit to the ideal target template -- not
         merely the lowest-energy one. Classification is only a nearest-template
         label, so an energetically-competitive but heavily puckered square-plane
         can still read as ``SPL`` yet sit far from the input geometry; ranking by
         template fit selects the cleanest realization.
    Energy (pool order) breaks ties.

    ``mols`` is assumed sorted best (lowest energy) first. Returns
    ``(chosen_mol, chosen_contract_mol)``. Strictly non-regressive: selection is
    skipped -- returning ``mols[0]`` -- when the target geometry has no discrete
    template, when no pooled conformer classifies as the target, or on any
    perception failure. The winding pass likewise falls back to the best-geometry
    (then lowest-energy) conformer when no winding is requested or none matches.

    ``reencode_cache``: forwarded to :func:`_reencode_key_matches` (see its docstring) so the
    early-exit re-scan below reuses whatever the pool-fill loop's ``accept_fn`` already
    computed for a mol, instead of paying the expensive re-encode a second time. ``None``
    (default) -> no caching, byte-identical to pristine.
    """
    from ..utils.oin_aligner import classify_and_fit

    # The atropisomer token requested by the OIN, if any (Y2 P2). None -- the case for every
    # OIN encoded without OIN_EMIT_AXIAL -- disables every axial-aware branch below, so this
    # function stays byte-identical to pristine for them.
    target_axial = parse_axial_token(getattr(parsed, "original_oin", None))

    # SL1 accept-first (opt-in, OIN_EARLY_EXIT / ff_params["early_exit"]): accept the first
    # conformer that INDEPENDENTLY re-encodes to the requested OIN's fac/mer key, short-
    # circuiting AHEAD of the geometry-classification gate below. The round-trip key is the
    # ground truth, so a key-exact conformer is correct even if its coordination sphere
    # classifies as a neighbouring template. Strictly non-regressive: on no match we fall
    # through to the geometry / winding / lowest-energy logic unchanged. Off-flag this whole
    # block is skipped, so selection is byte-identical to pristine.
    if early_exit:
        try:
            target_key = canonical_roundtrip_key(getattr(parsed, "original_oin", "") or "")
        except Exception:
            target_key = None
        if target_key is not None:
            require_no_stretch = clash.STRETCHED_BOND_ENABLED
            for m in mols:
                cmol = build_contract_mol(parsed, m)
                if cmol is None:
                    continue
                if _reencode_key_matches(
                    parsed,
                    m,
                    target_key,
                    cmol=cmol,
                    require_no_stretch=require_no_stretch,
                    cache=reencode_cache,
                ) and (target_axial is None or mol_axial_token(cmol) == target_axial):
                    logger.debug("accept-first: conformer re-encodes to the requested fac/mer key")
                    _telemetry.record("adapter.early_exit_hit")
                    return m, cmol
            _telemetry.record("adapter.early_exit_miss", n_mols=len(mols))

    target = _norm_geo_code(parsed.geo_code)
    expected_n = _expected_coordination_number(parsed.geo_code)

    scored = []  # (clash_vdw, fit_rmsd, energy_rank, mol, contract_mol)
    if target and expected_n is not None:
        for rank, m in enumerate(mols):
            try:
                cmol = build_contract_mol(parsed, m)
                if cmol is None:
                    continue
                vecs = _coordination_vectors(cmol, expected_n)
                if vecs is None:
                    continue
                label, fit = classify_and_fit(vecs, target)
                if _norm_geo_code(label) != target:
                    continue
                # Whole-complex vdW clash as the primary sort key -- gated OFF by default
                # (clash.VDW_ACCEPTANCE_ENABLED). The elimination study found coordination-
                # sphere template fit uncorrelated with real distortion (rho=0.22), so when
                # enabled we pick the least-clashing realization first (fit/energy break
                # ties). Off by default because even among conformers that classify as the
                # target, the least-clashing one has donors splayed to the edge of the gate
                # and re-perceives as detached -> round-trip regression (see clash.py). When
                # disabled, clash=0 for all -> (0, fit, rank) == the pre-A3 (fit, rank) order.
                clash_vdw = clash.mol_clash_count(m) if clash.VDW_ACCEPTANCE_ENABLED else 0
                scored.append((clash_vdw, fit, rank, m, cmol))
            except Exception:
                logger.debug("geometry perception failed for a conformer", exc_info=True)
                continue
        scored.sort(key=lambda t: (t[0], t[1], t[2]))

    # Axial-aware narrowing (Y2 P2): when the requested OIN carries an atropisomer token
    # (` |ax:-|`), keep only conformers whose own axial token reproduces it. The embed sets
    # a biaryl torsion stochastically, so without this the atropisomer is left to chance --
    # the same problem, and the same "selection beats construction" remedy, as the eta
    # winding pass below. Self-gating: an OIN encoded without OIN_EMIT_AXIAL carries no
    # token, so `target_axial` is None and this is a no-op. Strictly non-regressive: if no
    # conformer matches, the pool is left exactly as it was and we fall through unchanged.
    # The token is canonical (graph-derived), which is what makes it comparable between the
    # requested OIN and a freshly embedded conformer with different atom numbering.
    if target_axial:
        if scored:
            # t[4] is the CONTRACT mol: raw pool conformers are unsanitized, so axial
            # perception must run on the same prepared mol the re-encode uses.
            keep, blind = _axial_narrow(scored, lambda t: t[4], target_axial)
            if keep:
                logger.debug(
                    "axial-aware selection: %d/%d conformers match token %r",
                    len(keep),
                    len(scored),
                    target_axial,
                )
                scored = keep
            else:
                _axial_report_miss(len(scored), blind, target_axial, key="n_scored")
        else:
            keep_mols, blind = _axial_narrow(
                mols, lambda m: build_contract_mol(parsed, m), target_axial
            )
            if keep_mols:
                logger.debug(
                    "axial-aware selection: %d/%d pool conformers match token %r",
                    len(keep_mols),
                    len(mols),
                    target_axial,
                )
                mols = keep_mols
            else:
                _axial_report_miss(len(mols), blind, target_axial, key="n_mols")

    # Winding-aware pick: prefer the conformer whose re-encoded eta-ring winding
    # matches the requested OIN. Search geometry-eligible conformers first (best
    # geometry first); if geometry perception found none, fall back to searching
    # the whole energy-ranked pool so winding can still be honored.
    # The OIN-direct path constructs eta winding deterministically at placement
    # time, so the (expensive) re-encode-and-filter pass is unnecessary there.
    target_windings = _eta_winding_multiset(getattr(parsed, "original_oin", None))
    if honor_winding and target_windings:
        if scored:
            candidates = [(m, cmol) for (_clash, _fit, _rank, m, cmol) in scored]
        else:
            candidates = [(m, None) for m in mols]
        for m, cmol in candidates:
            if cmol is None:
                cmol = build_contract_mol(parsed, m)
            # Fast contract-mol re-encode; fall back to the full XYZ path.
            oin = _reencode_oin_fast(cmol) or _reencode_oin(m)
            if oin is None:
                continue
            if _eta_winding_multiset(oin) == target_windings:
                logger.debug("winding-aware selection: matched eta winding %s", target_windings)
                return m, cmol if cmol is not None else build_contract_mol(parsed, m)
        logger.debug(
            "winding-aware selection: no conformer matched eta winding %s; "
            "falling back to geometry/energy",
            target_windings,
        )
        _telemetry.record(
            "adapter.winding_fallthrough",
            n_targets=len(target_windings),
            n_candidates=len(candidates),
        )

    if scored:
        clash_vdw, fit, rank, m, cmol = scored[0]
        logger.debug(
            "geometry-aware selection: chose energy-rank %d as %s "
            "(vdW clashes %d, template fit %.4f, target %s) from %d matching conformer(s)",
            rank,
            target,
            clash_vdw,
            fit,
            target,
            len(scored),
        )
        return m, cmol

    # Fallback: lowest-energy conformer -- the pre-selection default behavior.
    _telemetry.record("adapter.geometry_select_fallthrough", target=target, n_mols=len(mols))
    return mols[0], build_contract_mol(parsed, mols[0])


def _verify_axial_honored(parsed, chosen_mol, cmol):
    """Warn when the requested atropisomer is NOT the one we are about to return.

    Every axial branch above is deliberately non-regressive: on any miss it falls through to
    the unfiltered pool. That safety has a cost -- the run then yields a structurally fine
    conformer of the WRONG atropisomer, whose round-trip key still matches (the key folds the
    axial token), so nothing downstream complains. Verifying once at the single exit turns
    that silent wrong answer into an observable event.
    """
    target = parse_axial_token(getattr(parsed, "original_oin", None))
    if not target:
        return
    probe = cmol if cmol is not None else build_contract_mol(parsed, chosen_mol)
    got = mol_axial_token(probe) if probe is not None else None
    if got == target:
        return
    logger.warning(
        "axial NOT honored: requested %r but returning %r -- the round-trip key folds the "
        "axial token, so this will not surface as a round-trip failure",
        target,
        got,
    )
    _telemetry.record("adapter.axial_not_honored", requested=target, got=got)


def _select_by_geometry(parsed, mols, honor_winding=True, early_exit=False, reencode_cache=None):
    """Conformer choice (see :func:`_select_by_geometry_impl`) plus an axial self-check.

    Kept as the module's entry point so every caller -- including the tests and
    ``tools/benchmark_generation.py``, which monkeypatch this name -- gets the verification.

    ``reencode_cache``: forwarded as-is; see :func:`_select_by_geometry_impl` and
    :func:`_reencode_key_matches`. ``None`` (default) -> byte-identical to pristine.
    """
    chosen, cmol = _select_by_geometry_impl(
        parsed,
        mols,
        honor_winding=honor_winding,
        early_exit=early_exit,
        reencode_cache=reencode_cache,
    )
    try:
        _verify_axial_honored(parsed, chosen, cmol)
    except Exception:  # a diagnostic must never break generation
        logger.debug("axial verification raised", exc_info=True)
    return chosen, cmol


class MetalloGenAdapter:
    """MetalloGen 3D-generation backend. Exposes ``generate(parsed)``."""

    def __init__(
        self,
        timeout: int = 300,
        dg_strategy: str = "single",
        ensemble_size: int = 1,
        optimizer: str | None = None,
        ff_preset: str | None = None,
        ff_params: dict | None = None,
        seed: int = 42,
    ) -> None:
        """Configure the backend (timeout, DG strategy, optimizer, FF knobs, seed)."""
        self.timeout = timeout
        self.dg_strategy = dg_strategy
        self.ensemble_size = ensemble_size
        # Base ETKDG seed threaded to generate_3d_structures (which offsets it
        # per embed attempt). A fixed seed keeps generation reproducible.
        self.seed = seed

        # Treat "FF" or "none" (case-insensitive) as None
        if optimizer is not None and optimizer.lower() in ("ff", "none"):
            optimizer = None

        # optimizer=None -> FF-relaxed geometry only (default; always available).
        # optimizer="xtb" -> refine the FF pool with standard g-xTB and energy-rank
        # (requires g-xTB binary + ase; degrades gracefully to FF if unavailable).
        self.optimizer = optimizer
        # FF convergence knobs (named preset + optional explicit overrides).
        self.ff_params = _resolve_ff_params(ff_preset, ff_params)

    def _oin_direct_enabled(self) -> bool:
        """Whether the OIN-direct assembly + winding path is active (SL2, default off).

        ``ff_params["oin_direct"]`` wins; otherwise the ``OIN_DIRECT_ASSEMBLY`` env
        var (``"1"`` to enable). Unset -> the byte-identical m-SMILES path.
        """
        if self.ff_params and self.ff_params.get("oin_direct"):
            return True
        return os.environ.get("OIN_DIRECT_ASSEMBLY", "0") == "1"

    def _direct_dg_enabled(self) -> bool:
        """Whether OIN-direct assembly feeds the DG embed -- the DEFAULT path (v0.4.4).

        Builds the ``MetalComplex`` straight from ``ParsedOIN`` (preserving winding +
        metal chirality, no lossy m-SMILES bridge) and runs it through the SAME distance-
        geometry embed + winding-search + early-exit/selection body as before -- so OIN
        format changes reach 3D generation without a translation layer. An A/B on a
        38-molecule stratified sample matched the m-SMILES path byte-for-byte per molecule
        (0 regressions, 0 gains), so it ships as the default; the m-SMILES bridge remains
        only as a fallback (see ``generate``). Opt out with ``OIN_DIRECT_DG=0`` or
        ``ff_params["direct_dg"]=False``. Distinct from ``oin_direct`` (SL2), which rigidly
        *constructs* haptic winding (kept opt-in, regresses eta).
        """
        if self.ff_params and "direct_dg" in self.ff_params:
            return bool(self.ff_params["direct_dg"])
        return os.environ.get("OIN_DIRECT_DG", "1") != "0"

    def generate(self, parsed: ParsedOIN) -> GeneratedStructure:
        """Generate a 3D structure for a parsed OIN via the MetalloGen engine."""
        if self._oin_direct_enabled():
            direct = self._maybe_generate_oin_direct(parsed)
            if direct is not None:
                return direct
            # A haptic face whose winding could not be constructed (a rare
            # mixed-denticity eta arm, or a non-inline OIN) fell through -- honor it
            # via the searched m-SMILES path below rather than ship an arbitrary face.

        # DEFAULT path (v0.4.4): build the complex straight from the OIN (no lossy m-SMILES
        # bridge) and run it through the SAME DG-embed + winding-search + early-exit body.
        # Keeping the metal + all ligand info (winding on the ligand, metal chirality) in one
        # representation means OIN-format changes reach 3D generation without a translation
        # layer. The m-SMILES bridge is retained ONLY as a fallback: if direct assembly raises
        # on some edge case, drop to the (winding-lossy) m-SMILES path rather than hard-fail.
        # Opt out of the direct path entirely with OIN_DIRECT_DG=0.
        prebuilt_complex = None
        msmiles = None
        if self._direct_dg_enabled():
            try:
                metal_frag, ligand_specs, geo = _prepare_ligand_fragments(parsed)
                prebuilt_complex = om.get_om_from_parsed(metal_frag, ligand_specs, geo)
                logger.debug("OIN %r -> direct MetalComplex (%s)", parsed.original_oin, geo)
            except Exception:
                logger.debug("OIN-direct assembly failed; falling back to m-SMILES", exc_info=True)
                prebuilt_complex = None
        if prebuilt_complex is None:
            msmiles = convert_parsed_to_msmiles(parsed)
            logger.debug("OIN %r -> m-SMILES %r", parsed.original_oin, msmiles)

        # Build the full energy-ranked conformer pool so geometry-aware selection
        # has candidates to choose among. The pool width is driven by
        # ``uff_pool_size`` (the UFF pre-pool); ``num_conformers`` asks the callee
        # to return the whole ranked list (not just the top-1). With an optimizer
        # set the pool is MACE-optimized regardless, so this adds no optimizer
        # cost over the previous fixed pool of 5.
        #
        # When the OIN encodes eta-ring winding, widen the pool (and the UFF
        # pre-pool that feeds it) so the requested ring face / diastereomer is
        # actually sampled -- winding-aware selection can only pick a winding that
        # exists in the pool.
        needs_winding = bool(_eta_winding_multiset(getattr(parsed, "original_oin", None)))
        base_pool = ETA_SELECT_POOL if needs_winding else DEFAULT_SELECT_POOL
        pool_n = max(self.ensemble_size, base_pool)

        # The MetalloGen engine prints progress/geometry to stdout; redirect it to
        # stderr so the oin2xyz CLI's stdout (the XYZ block) stays clean.
        # Extract deduplication params from ff_params if provided, else use defaults
        uff_pool_size = self.ff_params.get("uff_pool_size", 10) if self.ff_params else 10
        if needs_winding:
            # Diastereomer diversity comes from the UFF pre-pool; make sure it is
            # at least as wide as the selection pool.
            uff_pool_size = max(uff_pool_size, 2 * pool_n)
        rmsd_threshold = self.ff_params.get("rmsd_threshold", 0.5) if self.ff_params else 0.5
        energy_threshold = self.ff_params.get("energy_threshold", 2.0) if self.ff_params else 2.0

        clean_ff_params = {
            k: v
            for k, v in (self.ff_params or {}).items()
            if k
            not in [
                "uff_pool_size",
                "rmsd_threshold",
                "energy_threshold",
                "oin_direct",
                "early_exit",
            ]
        }

        # SL1 generate-until-key-exact early-exit. When enabled, hand the engine an
        # ``accept_fn`` that returns True as soon as an embedded conformer INDEPENDENTLY
        # re-encodes to the requested OIN's fac/mer key (SL0's ``canonical_roundtrip_key``),
        # so the attempt loop stops building the pool -- a pool of 1 when conformer-1
        # round-trips. The same flag turns on the adapter-side accept-first pick below.
        #
        # PROMOTED to default-ON in v0.4.4 (promote A/B --
        # tmCAT-tmPHOTO_xyz_dataset/results-v0.4.4-promote-ab/VALIDATION.md): on the
        # stratified worst-cohort sample it lifted byte-exact 44.7->60.5%, key-match
        # 55.3->73.7% with ZERO regressions, and ran ~5x faster (it picks the conformer
        # reproducing the requested fac/mer isomer instead of leaving it to geometry-
        # classification luck, and stops as soon as one appears). Non-regressive by
        # construction: if no pooled conformer matches the key it falls through to the prior
        # selection. Opt OUT via OIN_EARLY_EXIT=0 or ff_params["early_exit"]=False.
        if self.ff_params is not None and "early_exit" in self.ff_params:
            early_exit = bool(self.ff_params["early_exit"])
        else:
            early_exit = os.environ.get("OIN_EARLY_EXIT", "1") != "0"
        accept_fn = None
        # Per-generation memo for the expensive step-2 re-encode inside
        # _reencode_key_matches (a full XYZToSMILES().convert() round trip -- measured
        # 48-57s/call on an eta/haptic conformer). accept_fn below and
        # _select_by_geometry's own early-exit re-scan test the SAME mol objects against
        # the SAME predicate: every mol _select_by_geometry sees either already ran
        # through accept_fn during the pool-fill loop (a miss there is why the loop kept
        # going) or arrived via the never-tested stereo_rejects fallback -- the cache
        # handles both correctly (a genuine miss just recomputes once, same as before).
        # Fresh dict per generate() call -- no cross-molecule staleness, same pattern as
        # generator3d's alt_cache/PuLP memos. Harmless (empty, unused) when early_exit is
        # off, since accept_fn is never built and nothing ever writes to it.
        _reencode_cache: dict = {}
        if early_exit:
            try:
                _target_key = canonical_roundtrip_key(getattr(parsed, "original_oin", "") or "")
            except Exception:
                _target_key = None
            if _target_key is not None:
                _require_no_stretch = clash.STRETCHED_BOND_ENABLED
                # The round-trip key deliberately FOLDS the axial token (so the batch
                # harness is unaffected by OIN_EMIT_AXIAL). That makes the key alone an
                # unsound acceptance test once an axial token is requested: it would accept
                # the first key-matching conformer of EITHER atropisomer and stop the pool,
                # so the axial-aware pick in _select_by_geometry would never see an
                # alternative. Require both. None when the OIN carries no token, so this is
                # a no-op for every OIN encoded without the flag.
                _target_axial = parse_axial_token(getattr(parsed, "original_oin", None))

                def accept_fn(
                    mg_mol,
                    _pk=_target_key,
                    _rns=_require_no_stretch,
                    _ax=_target_axial,
                    _cache=_reencode_cache,
                ):
                    if _ax is None:
                        # byte-identical to the pre-axial predicate for every OIN that
                        # carries no token -- i.e. everything, absent OIN_EMIT_AXIAL.
                        return _reencode_key_matches(
                            parsed, mg_mol, _pk, require_no_stretch=_rns, cache=_cache
                        )
                    # Perceive on the CONTRACT mol (raw pool conformers are unsanitized, so
                    # axial perception throws on them); build it once and share it with the
                    # key check rather than letting each build its own.
                    cmol = build_contract_mol(parsed, mg_mol)
                    if cmol is None:
                        return False
                    return (
                        _reencode_key_matches(
                            parsed, mg_mol, _pk, cmol=cmol, require_no_stretch=_rns, cache=_cache
                        )
                        and mol_axial_token(cmol) == _ax
                    )

        with contextlib.redirect_stdout(sys.stderr):
            mols = generate_3d_structures(
                msmiles,
                num_conformers=pool_n,
                optimizer=self.optimizer,
                ff_params=clean_ff_params,
                uff_pool_size=uff_pool_size,
                rmsd_threshold=rmsd_threshold,
                energy_threshold=energy_threshold,
                timeout=self.timeout,
                # Bound the FF-only attempt loop by wall-clock. self.timeout is the
                # same per-molecule budget the harness already plumbs (300 s full /
                # 30 s quick); reuse it so a pathological embed (ZIHGEE ~1696 s) fails
                # fast instead of running all 250 attempts. Kept distinct from the ASE
                # `timeout` semantics, which cap a single optimizer call.
                embed_time_budget=self.timeout,
                seed=self.seed,
                accept_fn=accept_fn,
                metal_complex=prebuilt_complex,
            )
        if not mols:
            raise ValueError(
                f"MetalloGen failed to generate any conformers for m-SMILES {msmiles!r}"
            )

        # Prefer the conformer whose re-perceived coordination geometry matches the
        # requested OIN code (fixes floppy-donor cases where a distorted geometry is
        # energetically competitive), falling back to the lowest-energy conformer. With
        # early-exit on, the accept-first pass inside also short-circuits to a key-exact
        # conformer (belt-and-suspenders with the engine's accept_fn: the returned pool may be
        # size 1 already, but a full-pool fallback still gets the honest accept-first pick).
        chosen_mol, mol = _select_by_geometry(
            parsed, mols, early_exit=early_exit, reencode_cache=_reencode_cache
        )

        xyz_str = get_xyz_string(chosen_mol)
        # Contract mol: MetalloGen connectivity+coords, OIN bond orders + 3D stereo.
        # None on failure -> callers fall back to coordinate re-perception.
        return GeneratedStructure(xyz=xyz_str, mol=mol)

    def _maybe_generate_oin_direct(self, parsed: ParsedOIN):
        """OIN-direct assembly with deterministic eta-ring winding (SL2).

        Builds the ``MetalComplex`` straight from ``ParsedOIN`` (bypassing the
        winding-lossy m-SMILES bridge) with the winding target attached to each
        haptic ligand. Haptic complexes are placed rigidly (option 3) so the ring
        face -- and thus the winding -- is *constructed* by ``_place_haptic`` rather
        than sampled, which retires the wide winding-search pool: the ``ETA_SELECT_POOL``
        widening and its UFF doubling are dropped, and the re-encode winding filter in
        ``_select_by_geometry`` is skipped. Non-haptic complexes take the same DG
        embed as the m-SMILES path (equivalent geometry; the only difference is the
        assembly plumbing).

        Returns ``None`` (so the caller uses the searched m-SMILES path) when a haptic
        face has no constructed winding -- a mixed-denticity eta arm that takes the
        chelate branch, or a non-inline OIN -- so such rings are never shipped with an
        arbitrary face and no fallback search.
        """
        metal_frag, ligand_specs, geo = _prepare_ligand_fragments(parsed)
        metal_complex = om.get_om_from_parsed(metal_frag, ligand_specs, geo)

        haptic_ligs = [
            lig for lig in metal_complex.ligands if any(len(bi[0]) > 1 for bi in lig.binding_infos)
        ]
        if haptic_ligs and not all(lig.winding is not None for lig in haptic_ligs):
            return None
        is_haptic = bool(haptic_ligs)
        logger.debug("OIN %r -> direct MetalComplex (%s)", parsed.original_oin, geo)

        # Winding is constructed, not searched -> the default (narrow) pool suffices;
        # no ETA_SELECT_POOL widening, no UFF pre-pool doubling.
        pool_n = max(self.ensemble_size, DEFAULT_SELECT_POOL)
        uff_pool_size = self.ff_params.get("uff_pool_size", 10) if self.ff_params else 10
        rmsd_threshold = self.ff_params.get("rmsd_threshold", 0.5) if self.ff_params else 0.5
        energy_threshold = self.ff_params.get("energy_threshold", 2.0) if self.ff_params else 2.0

        clean_ff_params = {
            k: v
            for k, v in (self.ff_params or {}).items()
            if k
            not in [
                "uff_pool_size",
                "rmsd_threshold",
                "energy_threshold",
                "oin_direct",
                "early_exit",
            ]
        }
        if is_haptic:
            # Route haptic complexes through the rigid placer (option 3), where
            # _place_haptic constructs the winding from ligand.winding.
            clean_ff_params["kabsch_only"] = True

        with contextlib.redirect_stdout(sys.stderr):
            mols = generate_3d_structures(
                None,
                num_conformers=pool_n,
                optimizer=self.optimizer,
                ff_params=clean_ff_params,
                uff_pool_size=uff_pool_size,
                rmsd_threshold=rmsd_threshold,
                energy_threshold=energy_threshold,
                timeout=self.timeout,
                embed_time_budget=self.timeout,
                seed=self.seed,
                metal_complex=metal_complex,
            )
        if not mols:
            raise ValueError(
                f"MetalloGen (OIN-direct) failed to generate any conformers for OIN "
                f"{parsed.original_oin!r}"
            )

        # Winding is already constructed -> skip the re-encode winding filter.
        chosen_mol, mol = _select_by_geometry(parsed, mols, honor_winding=False)
        xyz_str = get_xyz_string(chosen_mol)
        return GeneratedStructure(xyz=xyz_str, mol=mol)


class OIN3DGeneratorMetallogen:
    """Standalone parse+generate wrapper.

    Retained for direct use; prefer ``OIN3DGenerator(engine="metallogen")`` for
    the integrated seam.
    """

    def __init__(
        self,
        timeout: int = 300,
        ensemble_size: int = 1,
        dg_strategy: str = "single",
        optimizer: str | None = None,
        ff_preset: str | None = None,
        ff_params: dict | None = None,
    ) -> None:
        """Build the OIN parser and the underlying ``MetalloGenAdapter``."""
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
        """Parse an OIN string and generate its 3D structure."""
        return self.adapter.generate(self.parser.parse(oin_string))
