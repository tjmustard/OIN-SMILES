"""Module for generating rdkit molobj/smiles/molecular graph from free atoms.

Main implementation by Jan H. Jensen, based on the paper

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
    to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
    DOI: 10.1002/bkcs.10334

Modified by Maria Harris Rasmussen 2024
"""

import contextlib
import copy
import itertools

try:
    from rdkit.Chem import rdEHTTools  # requires RDKit 2019.9.1 or later
except ImportError:
    rdEHTTools = None

import logging
import os
import sys
from collections import OrderedDict, defaultdict

import networkx as nx
import numpy as np
from rdkit import Chem

# Single-source lever registry. This was a guarded import with a local fallback while the
# registry did not yet exist on every branch; now that it ships, the fallback is a TRAP --
# it defaulted every lever to OFF, which would silently revert the six levers promoted to
# default-ON in v0.4.5 if the import ever broke. A missing registry should be a loud
# ImportError, not six quiet behaviour changes.
#
# The registry's purpose is to close a sense-inversion trap that is live in older code:
# os.environ.get("X") is truthy for the string "0", so the obvious way to opt out of a
# bare-truthiness lever turns it ON.
from ..oin.levers import lever_enabled as _lever_enabled

logger = logging.getLogger(__name__)

global __ATOM_LIST__
__ATOM_LIST__ = [
    "h",
    "he",
    "li",
    "be",
    "b",
    "c",
    "n",
    "o",
    "f",
    "ne",
    "na",
    "mg",
    "al",
    "si",
    "p",
    "s",
    "cl",
    "ar",
    "k",
    "ca",
    "sc",
    "ti",
    "v",
    "cr",
    "mn",
    "fe",
    "co",
    "ni",
    "cu",
    "zn",
    "ga",
    "ge",
    "as",
    "se",
    "br",
    "kr",
    "rb",
    "sr",
    "y",
    "zr",
    "nb",
    "mo",
    "tc",
    "ru",
    "rh",
    "pd",
    "ag",
    "cd",
    "in",
    "sn",
    "sb",
    "te",
    "i",
    "xe",
    "cs",
    "ba",
    "la",
    "ce",
    "pr",
    "nd",
    "pm",
    "sm",
    "eu",
    "gd",
    "tb",
    "dy",
    "ho",
    "er",
    "tm",
    "yb",
    "lu",
    "hf",
    "ta",
    "w",
    "re",
    "os",
    "ir",
    "pt",
    "au",
    "hg",
    "tl",
    "pb",
    "bi",
    "po",
    "at",
    "rn",
    "fr",
    "ra",
    "ac",
    "th",
    "pa",
    "u",
    "np",
    "pu",
]


global atomic_valence
global atomic_valence_electrons

atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3, 4]
atomic_valence[6] = [4, 2]
atomic_valence[7] = [3, 4]
atomic_valence[8] = [2, 1, 3]  # [2,1,3]
atomic_valence[9] = [1]
atomic_valence[13] = [3, 4]
atomic_valence[14] = [4]
atomic_valence[15] = [3, 5]  # [5,4,3]
atomic_valence[16] = [2, 4, 6]  # [6,3,2]
atomic_valence[17] = [1]
atomic_valence[18] = [0]
atomic_valence[32] = [4]
atomic_valence[33] = [5, 3]
atomic_valence[35] = [1]
atomic_valence[34] = [2]
atomic_valence[52] = [2]
atomic_valence[53] = [1]

atomic_valence[21] = [20]
atomic_valence[22] = [20]
atomic_valence[23] = [20]
atomic_valence[24] = [20]
atomic_valence[25] = [20]
atomic_valence[26] = [20]
atomic_valence[27] = [20]
atomic_valence[28] = [20]
atomic_valence[29] = [20]
atomic_valence[30] = [20]

atomic_valence[39] = [20]
atomic_valence[40] = [20]
atomic_valence[41] = [20]
atomic_valence[42] = [20]
atomic_valence[43] = [20]
atomic_valence[44] = [20]
atomic_valence[45] = [20]
atomic_valence[46] = [20]
atomic_valence[47] = [20]
atomic_valence[48] = [20]


atomic_valence[57] = [20]
atomic_valence[72] = [20]
atomic_valence[73] = [20]
atomic_valence[74] = [20]
atomic_valence[75] = [20]
atomic_valence[76] = [20]
atomic_valence[77] = [20]
atomic_valence[78] = [20]
atomic_valence[79] = [20]
atomic_valence[80] = [20]


atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[13] = 3
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[18] = 8
atomic_valence_electrons[32] = 4
atomic_valence_electrons[33] = 5
atomic_valence_electrons[35] = 7
atomic_valence_electrons[34] = 6
atomic_valence_electrons[52] = 6
atomic_valence_electrons[53] = 7


def str_atom(atom):
    """Convert integer atom to string atom."""
    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    return atom


def int_atom(atom):
    """Convert str atom to integer atom."""
    global __ATOM_LIST__
    # print(atom)
    atom = atom.lower()
    return __ATOM_LIST__.index(atom) + 1


def get_UA(maxValence_list, valence_list):
    """Unsaturated atoms and their degree of unsaturation.

    Same loop as before, with the subtraction done once instead of twice and the
    ``append`` bound outside. Deliberately *not* vectorised: at ~60 atoms a numpy
    round trip costs more than the Python loop it replaces (measured -- see
    ``docs/ENCODER_PERF_v0.4.5.md``).
    """
    UA = []
    DU = []
    add_UA = UA.append
    add_DU = DU.append
    for i, (maxValence, valence) in enumerate(zip(maxValence_list, valence_list)):
        d = maxValence - valence
        if d > 0:
            add_UA(i)
            add_DU(d)
    return UA, DU


def get_BO(AC, UA, DU, valences, UA_pairs, use_graph=True):
    """"""
    BO = AC.copy()
    DU_save = []

    while DU_save != DU:
        for i, j in UA_pairs:
            BO[i, j] += 1
            BO[j, i] += 1

        # .tolist() not list(): see valences_not_too_large -- boxed np.int64 makes
        # every downstream subtraction in get_UA go through numpy scalar dispatch.
        BO_valence = BO.sum(axis=1).tolist()
        DU_save = copy.copy(DU)
        UA, DU = get_UA(valences, BO_valence)
        UA_pairs = get_UA_pairs(UA, AC, DU, use_graph=use_graph)[0]
    return BO


def valences_not_too_large(BO, valences):
    """Whether no atom carries more bonds than its assigned valence.

    Same loop as before. The one change is ``.tolist()`` in place of ``list()``:
    ``list()`` over a numpy array yields ``np.int64`` *objects*, and every subsequent
    comparison then goes through numpy's scalar machinery, which is several times
    slower than the equivalent on a Python ``int``. ``.tolist()`` converts once.
    """
    number_of_bonds_list = BO.sum(axis=1).tolist()
    for valence, number_of_bonds in zip(valences, number_of_bonds_list):
        if number_of_bonds > valence:
            return False

    return True


def charge_is_OK(
    BO,
    AC,
    charge,
    DU,
    atomic_valence_electrons,
    atoms,
    valences,
    allow_charged_fragments=True,
    allow_carbenes=True,
):
    """Whether the formal charges implied by ``BO`` sum to the requested ``charge``.

    Same loop, same branch ladder, same order-dependent carbon corrections. This was
    the second-largest cost in the encoder (measured 14.4-16.7 s of a 49-58 s encode on
    QIDKUL_comp_0 / QIDKIZ_comp_0) and essentially all of it was interpreter overhead,
    removed here in three ways:

    * ``list(BO[i, :]).count(1)`` allocated a fresh 59-element list of boxed
      ``np.int64`` per carbon atom, per call. One vectorised ``(BO == 1).sum(axis=1)``
      gives the same counts for every atom at once. (Counting row entries equal to 1 is
      exactly what ``.count(1)`` did.)
    * :func:`get_atomic_charge` was one Python call per atom -- 4.4 million calls in a
      single encode of QIDKUL_comp_0. Inlined here as the identical ``elif`` ladder.
    * ``list(BO.sum(axis=1))`` yields boxed ``np.int64`` objects, so every downstream
      ``1 - v``, ``v == 2`` and ``Q += q`` ran through numpy's scalar machinery.
      ``.tolist()`` converts once and the arithmetic is then plain Python ``int``.

    The corrections are left as a Python loop on purpose: they read the *running* ``Q``
    (``if ns == 3 and Q + 1 < charge``), so they cannot be reordered or vectorised
    without changing the answer.

    The original also built a ``q_list`` that it never read; that dead accumulator is
    dropped.
    """
    # total charge
    Q = 0

    if allow_charged_fragments:
        BO_valences = BO.sum(axis=1).tolist()
        # Only carbons consult it, so pay for it only when there is a carbon.
        n_single_list = (BO == 1).sum(axis=1).tolist() if 6 in atoms else None
        for i, atom in enumerate(atoms):
            v = BO_valences[i]
            # get_atomic_charge, inlined verbatim
            if atom == 1:
                q = 1 - v
            elif atom == 5:
                q = 3 - v
            elif atom == 6 and v == 2:
                q = 0
            elif atom == 13:
                q = 3 - v
            elif atom == 15 and v == 5:
                q = 0
            elif atom == 16 and v == 6:
                q = 0
            elif atom == 16 and v == 4:  # testing for sulphur
                q = 0
            elif atom == 16 and v == 5:
                q = 1
            else:
                q = atomic_valence_electrons[atom] - 8 + v
            Q += q
            if atom == 6:
                number_of_single_bonds_to_C = n_single_list[i]
                if not allow_carbenes and number_of_single_bonds_to_C == 2 and v == 2:
                    logger.debug("found illegal carbene")
                    Q += 1
                if number_of_single_bonds_to_C == 3 and Q + 1 < charge:
                    Q += 2
    return charge == Q


def BO_is_OK(
    BO,
    AC,
    charge,
    DU,
    atomic_valence_electrons,
    atoms,
    valences,
    allow_charged_fragments=True,
    allow_carbenes=True,
):
    """Sanity of bond-orders.

    Args:
        BO -
        AC -
        charge -
        DU -


    optional
        allow_charges_fragments -


    Returns:
        boolean - true of molecule is OK, false if not
    """
    if not valences_not_too_large(BO, valences):
        return False

    # Left as-is on purpose. `BO.sum() - AC.sum()` avoids the NxN temporary and is
    # algebraically identical, but measured 0.72x -- two reductions cost more than one
    # subtract-plus-reduction at this matrix size. See docs/ENCODER_PERF_v0.4.5.md.
    check_sum = (BO - AC).sum() == sum(DU)
    check_charge = charge_is_OK(
        BO,
        AC,
        charge,
        DU,
        atomic_valence_electrons,
        atoms,
        valences,
        allow_charged_fragments,
        allow_carbenes=True,
    )

    if check_charge and check_sum:
        return True

    return False


def get_atomic_charge(atom, atomic_valence_electrons, BO_valence):
    """"""
    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 6 and BO_valence == 2:
        charge = 0
    elif atom == 13:
        charge = 3 - BO_valence
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    elif atom == 16 and BO_valence == 4:  # testing for sulphur
        charge = 0
    elif atom == 16 and BO_valence == 5:
        charge = 1

    else:
        charge = atomic_valence_electrons - 8 + BO_valence

    return charge


def BO2mol(
    mol,
    BO_matrix,
    atoms,
    atomic_valence_electrons,
    mol_charge,
    allow_charged_fragments=True,
    use_atom_maps=True,
):
    """Based on code written by Paolo Toscani.

    From bond order, atoms, valence structure and total charge, generate an
    rdkit molecule.

    Args:
        mol - rdkit molecule
        BO_matrix - bond order matrix of molecule
        atoms - list of integer atomic symbols
        atomic_valence_electrons -
        mol_charge - total charge of molecule

    optional:
        allow_charged_fragments - bool - allow charged fragments

    Returns:
        mol - updated rdkit molecule with bond connectivity
    """
    length_bo = len(BO_matrix)
    length_atoms = len(atoms)
    BO_valences = list(BO_matrix.sum(axis=1))

    if length_bo != length_atoms:
        raise RuntimeError(
            "sizes of adjMat ({0:d}) and Atoms {1:d} differ".format(length_bo, length_atoms)
        )

    rwMol = Chem.RWMol(mol)

    bondTypeDict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }

    for i in range(length_bo):
        for j in range(i + 1, length_bo):
            bo = int(round(BO_matrix[i, j]))
            if bo == 0:
                continue
            bt = bondTypeDict.get(bo, Chem.BondType.SINGLE)
            rwMol.RemoveBond(i, j)  # added this for TMC procedure
            rwMol.AddBond(i, j, bt)

    mol = rwMol.GetMol()

    if allow_charged_fragments:
        mol = set_atomic_charges(
            mol,
            atoms,
            atomic_valence_electrons,
            BO_valences,
            BO_matrix,
            mol_charge,
            use_atom_maps=use_atom_maps,
        )
    else:
        mol = set_atomic_radicals(
            mol,
            atoms,
            atomic_valence_electrons,
            BO_valences,
            use_atom_maps=use_atom_maps,
        )

    Chem.SanitizeMol(mol)

    return mol


def set_atomic_charges(
    mol,
    atoms,
    atomic_valence_electrons,
    BO_valences,
    BO_matrix,
    mol_charge,
    use_atom_maps=True,
):
    """"""
    q = 0
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        if use_atom_maps:
            a.SetAtomMapNum(i + 1)
        charge = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])
        q += charge
        if atom == 6:
            number_of_single_bonds_to_C = list(BO_matrix[i, :]).count(1)
            if BO_valences[i] == 2:
                # q += 1
                a.SetNumRadicalElectrons(2)
                charge = 0
            if number_of_single_bonds_to_C == 3 and q + 1 < mol_charge:
                q += 2
                charge = 1

        if abs(charge) > 0:
            a.SetFormalCharge(int(charge))

    # mol = clean_charges(mol)

    return mol


def set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences, use_atom_maps=True):
    """The number of radical electrons = absolute atomic charge."""
    atomic_valence[8] = [2, 1]
    atomic_valence[7] = [3, 2]
    atomic_valence[6] = [4, 2]

    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        if use_atom_maps:
            a.SetAtomMapNum(i + 1)
        charge = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])

        if abs(charge) > 0:
            a.SetNumRadicalElectrons(abs(int(charge)))

    return mol


# --- memo for AC2BO's candidate-generation loop -------------------------------------
#
# ``xyz2mol.py::_select_lig_mol`` runs a charge/carbene ladder that calls
# ``AC2mol`` -> ``AC2BO`` up to five times on the *same* adjacency matrix, and
# ``_rescue_unusable_perception`` sweeps up to eight more charges over one AC of its own.
# Only ``BO_is_OK`` / ``charge_is_OK`` read ``charge``: the candidate-generation half of
# ``AC2BO``'s loop (``get_UA`` -> ``get_UA_pairs`` -> ``get_BO``) therefore recomputes
# bit-identical results on every arm. Measured on QIDKUL_comp_0: 60 001 matching-running
# ``get_UA_pairs`` calls for 26 668 distinct results, and 180 005 ``get_bonds`` calls for
# 4 007 distinct results.
#
# A small LRU over adjacency matrices rather than a single slot, because a round trip
# re-perceives *the same ligand connectivity* many times -- the input structure, then
# every generated conformer the SL1 accept-first check re-encodes -- and each of those is
# a fresh AC array with identical contents.
#
# Keyed on the exact bytes (plus shape and dtype) of the matrix, and only *read* when the
# caller's ``AC`` **is** an array the cache holds a live reference to. That identity test
# is what makes it safe: no ``id()`` recycling can alias a stale entry, and an unrelated
# caller simply misses and recomputes exactly as before.
_AC2BO_SLOTS = OrderedDict()

# How many distinct adjacency matrices to keep. Small: the working set of one round trip
# is a handful of ligand fragments.
_AC2BO_MEMO_SLOTS = 6

# Total cached entries across all slots. Above this the cache stops growing and calls just
# recompute -- slower, never wrong. A blow-up guard: one capped AC2BO search produces
# ~29 000 entries, so this leaves room for the working set without unbounded growth.
_AC2BO_MEMO_MAX = 200_000


def _ac2bo_memo_for(AC):
    """The memo dicts for ``AC``, or ``(None, None)`` if no slot holds this array."""
    slot = _AC2BO_SLOTS.get(id(AC))
    # `is` and not `==`: the id() lookup is only a fast index into slots whose arrays we
    # keep alive, and this check is what rejects a recycled id.
    if slot is not None and slot["ac"] is AC:
        return slot["bonds"], slot["uap"]
    return None, None


def _ac2bo_memo_entries():
    """Total cached entries across all slots."""
    return sum(len(s["bonds"]) + len(s["uap"]) for s in _AC2BO_SLOTS.values())


def _ac2bo_memo_anchor(AC):
    """Register ``AC`` with the cache, reusing entries when its contents are already known.

    Called once per ``AC2BO`` entry -- five or so times per encode -- so hashing the matrix
    here costs nothing measurable. Shape and dtype join the byte content in the tag because
    ``tobytes`` alone would not distinguish differently-shaped matrices.
    """
    key = id(AC)
    slot = _AC2BO_SLOTS.get(key)
    if slot is not None and slot["ac"] is AC:
        _AC2BO_SLOTS.move_to_end(key)
        return

    tag = (AC.shape, AC.dtype.str, AC.tobytes())
    # Snapshot the items: this loop deletes from the mapping it walks. Nothing in the
    # library calls AC2BO from more than one thread today (the only ThreadPoolExecutor is
    # around xtb optimisation, which does not perceive), but a snapshot makes the walk
    # immune to a concurrent anchor resizing the dict, and at <= 6 slots it is free.
    for other_key, other in list(_AC2BO_SLOTS.items()):
        if other["tag"] == tag:
            # Same connectivity, different array object (a re-perceived conformer, or the
            # charge sweep rebuilding its AC): adopt the entries and re-anchor the identity
            # test onto the object this caller will pass down.
            del _AC2BO_SLOTS[other_key]
            other["ac"] = AC
            _AC2BO_SLOTS[key] = other
            return

    if _ac2bo_memo_entries() >= _AC2BO_MEMO_MAX:
        _AC2BO_SLOTS.clear()
    _AC2BO_SLOTS[key] = {"ac": AC, "tag": tag, "bonds": {}, "uap": {}}
    while len(_AC2BO_SLOTS) > _AC2BO_MEMO_SLOTS:
        _AC2BO_SLOTS.popitem(last=False)


def _ac2bo_memo_clear():
    """Drop every memo slot (frees its memory). Used by tests."""
    _AC2BO_SLOTS.clear()


def get_bonds(UA, AC):
    """Bonds of ``AC`` whose both ends are unsaturated.

    A function of ``(AC, UA)`` alone, so it is memoised against the current AC slot.
    Callers mutate the list they get back (``get_UA_pairs`` appends the virtual-node
    edges to it), so a hit is rehydrated into a fresh list and the snapshot stored in
    the memo is an immutable tuple taken before any caller can touch it.
    """
    memo, _ = _ac2bo_memo_for(AC)
    key = None
    if memo is not None:
        key = tuple(UA)
        hit = memo.get(key)
        if hit is not None:
            return list(hit)

    bonds = []

    for k, i in enumerate(UA):
        for j in UA[k + 1 :]:
            if AC[i, j] == 1:
                bonds.append(tuple(sorted([i, j])))

    if key is not None and _ac2bo_memo_entries() < _AC2BO_MEMO_MAX:
        memo[key] = tuple(bonds)
    return bonds


# Env lever, default OFF (unset == nx.max_weight_matching, the historical matcher, so the
# default path is byte-identical). The graph built in ``get_UA_pairs`` carries **no weight
# attributes**, so ``nx.max_weight_matching`` is solving maximum *cardinality* matching with
# a general (Blossom) weighted matcher. Alternatives can return a *different* matching of the
# same size, which is a different Kekule structure and therefore a different perceived BO --
# hence the lever. Measured in docs/VALENCE_SEARCH_v0.4.5.md.
_MATCHER_ENV = "OIN_VALENCE_MATCHER"


def _maximum_matching(G):
    """Maximum matching of an unweighted graph, matcher selected by ``OIN_VALENCE_MATCHER``.

    ``nx`` (default) is the historical ``nx.max_weight_matching``. ``maxcard`` asks the same
    Blossom implementation for the max-cardinality variant. ``greedy`` is
    ``nx.maximal_matching``, which is *maximal* but not *maximum* -- it is included only so
    the cost/fidelity trade-off can be measured, and it can return a smaller matching.
    """
    which = os.environ.get(_MATCHER_ENV) or "nx"
    if which == "nx":
        return nx.max_weight_matching(G)
    if which == "maxcard":
        return nx.max_weight_matching(G, maxcardinality=True)
    if which == "greedy":
        return nx.maximal_matching(G)
    logger.warning("%s=%r is not a known matcher; using nx", _MATCHER_ENV, which)
    return nx.max_weight_matching(G)


def get_UA_pairs(UA, AC, DU, use_graph=True):
    """Maximum matching over the unsaturated-atom bond graph.

    Memoised against the current AC slot on ``(UA, du > 1 pattern)``. That key is
    sufficient -- not merely convenient -- because ``DU`` is read **only** through the
    ``du > 1`` predicate below, which decides how many virtual matching nodes get
    allocated and to which atoms. Everything else the result depends on comes from
    ``UA`` and ``AC``. So two calls agreeing on that key are handed the identical edge
    list in the identical insertion order, and ``nx.max_weight_matching`` is
    deterministic on identical input; the memo can only make the same answer arrive
    sooner.

    Only the ``use_graph`` result is cached. The ``len(bonds) == 0`` early-out is
    already cheap (its only real work, ``get_bonds``, is memoised separately) and the
    ``use_graph=False`` combinatorial branch is not reached by this codebase.
    """
    _, memo = _ac2bo_memo_for(AC)
    key = None
    if memo is not None and use_graph:
        key = (tuple(UA), tuple(du > 1 for du in DU))
        hit = memo.get(key)
        if hit is not None:
            # Fresh list per call: the caller is free to mutate what it gets back.
            return [list(hit)]

    N_UA = 10000
    matching_ids = dict()
    matching_ids2 = dict()
    for i, du in zip(UA, DU):
        if du > 1:
            matching_ids[i] = N_UA
            matching_ids2[N_UA] = i
            N_UA += 1

    bonds = get_bonds(UA, AC)
    for i, j in bonds:
        if i in matching_ids:
            bonds.append(tuple(sorted([matching_ids[i], j])))

        elif j in matching_ids:
            bonds.append(tuple(sorted([i, matching_ids[j]])))

    if len(bonds) == 0:
        return [()]

    if use_graph:
        G = nx.Graph()
        G.add_edges_from(bonds)
        AC2BO_STATS["matching_calls"] += 1
        UA_pairs = [list(_maximum_matching(G))]
        UA_pair = UA_pairs[0]

        remove_pairs = []
        add_pairs = []
        for i, j in UA_pair:
            if i in matching_ids2 and j in matching_ids2:
                remove_pairs.append(tuple([i, j]))
                add_pairs.append(tuple([matching_ids2[i], matching_ids2[j]]))
                # UA_pair.remove(tuple([i,j]))
                # UA_pair.append(tuple([matching_ids2[i], matching_ids2[j]]))
            elif i in matching_ids2:
                # UA_pair.remove(tuple([i,j]))
                remove_pairs.append(tuple([i, j]))
                add_pairs.append(tuple([matching_ids2[i], j]))
                # UA_pair.append(tuple([matching_ids2[i],j]))
            elif j in matching_ids2:
                remove_pairs.append(tuple([i, j]))
                add_pairs.append(tuple([i, matching_ids2[j]]))

                # UA_pair.remove(tuple([i,j]))
                # UA_pair.append(tuple([i,matching_ids2[j]]))
        for p1, p2 in zip(remove_pairs, add_pairs):
            UA_pair.remove(p1)
            UA_pair.append(p2)
        if key is not None and _ac2bo_memo_entries() < _AC2BO_MEMO_MAX:
            memo[key] = tuple(UA_pair)
        return [UA_pair]

    max_atoms_in_combo = 0
    UA_pairs = [()]
    for combo in list(itertools.combinations(bonds, int(len(UA) / 2))):
        flat_list = [item for sublist in combo for item in sublist]
        atoms_in_combo = len(set(flat_list))
        if atoms_in_combo > max_atoms_in_combo:
            max_atoms_in_combo = atoms_in_combo
            UA_pairs = [combo]

        elif atoms_in_combo == max_atoms_in_combo:
            UA_pairs.append(combo)

    return UA_pairs


# Above this many total per-atom valence combinations, AC2BO skips the O/N/C/P/S
# valence-ordering heuristic (which would materialise the full Cartesian product and
# hang on a large conjugated ligand) and iterates a bounded lazy product instead.
# Chosen well above any ligand that currently perceives quickly and far below the
# exponential blow-ups of the encode-fail timeout cohort, so sub-cap ligands (every
# currently-encodable structure) are byte-identical.
_VALENCE_COMBO_CAP = 500_000

# In the over-cap fallback the main loop iterates the *unsorted* lazy product and
# early-returns on the first valid assignment; this bounds how many candidates it will
# try before giving up (returning best_BO, which downstream perception then judges).
# Far smaller than the sort cap because each iteration runs the full BO/charge check --
# enough to catch a chemically-sensible assignment near the front of the product order,
# small enough that a genuinely unperceivable large ligand fails fast instead of hanging.
_VALENCE_FALLBACK_TRIES = 20_000

# Env lever, default OFF (unset == the historical 20 000, so the default path is
# byte-identical). Set OIN_VALENCE_FALLBACK_TRIES=<int> to change how many candidate
# valence assignments the over-cap fallback grinds before giving up and returning
# ``best_BO``. This is NOT byte-identical for over-cap ligands: ``best_BO`` is updated
# with ``BO.sum() >= best_BO.sum()`` (note ``>=``), so the *last* candidate attaining the
# maximum sum wins and a smaller budget can select a different one. See
# docs/VALENCE_SEARCH_v0.4.5.md for the measured stability table.
_FALLBACK_TRIES_ENV = "OIN_VALENCE_FALLBACK_TRIES"


def _fallback_tries():
    """How many over-cap candidates to try. Env-overridable; default is unchanged."""
    raw = os.environ.get(_FALLBACK_TRIES_ENV)
    if not raw:
        return _VALENCE_FALLBACK_TRIES
    try:
        value = int(raw)
    except ValueError:
        logger.warning(
            "%s=%r is not an integer; using %d", _FALLBACK_TRIES_ENV, raw, _VALENCE_FALLBACK_TRIES
        )
        return _VALENCE_FALLBACK_TRIES
    if value < 1:
        logger.warning(
            "%s=%d is < 1; using %d", _FALLBACK_TRIES_ENV, value, _VALENCE_FALLBACK_TRIES
        )
        return _VALENCE_FALLBACK_TRIES
    return value


# Deterministic counters for the valence search. Wall clock is unusable on this host
# (the release sweep keeps load above 12), so every claim in
# docs/VALENCE_SEARCH_v0.4.5.md is made in these instead.
AC2BO_STATS = {
    "ac2bo_calls": 0,  # AC2BO invocations
    "over_cap_calls": 0,  # ... of which took the bounded lazy-product branch
    "over_cap_ordered_calls": 0,  # ... of which enumerated in heuristic order (lever ON)
    "over_cap_filtered_calls": 0,  # ... of which enumerated charge-feasible-only (lever ON)
    "over_cap_infeasible": 0,  # ... of which had NO feasible candidate at all (charge is wrong)
    "over_cap_filter_unsupported": 0,  # ... of which the filter declined (metal in the fragment)
    "candidates": 0,  # candidate valence assignments examined, all branches
    "over_cap_candidates": 0,  # ... of which on the over-cap branch
    "found_valid": 0,  # AC2BO calls that early-returned a valid Lewis structure
    "over_cap_found_valid": 0,  # ... of which were over-cap
    "over_cap_exhausted": 0,  # over-cap calls that fell through and returned best_BO
    "over_cap_best_bo_improved": 0,  # best_BO reassignments on the over-cap branch
    "matching_calls": 0,  # nx.max_weight_matching invocations
}


def reset_ac2bo_stats():
    """Zero the deterministic counters. Call before a measured encode."""
    for key in AC2BO_STATS:
        AC2BO_STATS[key] = 0


def possible_valences(AC_valence, atoms, allow_carbenes=True):
    """The per-atom candidate valence lists whose Cartesian product ``AC2BO`` searches.

    Extracted verbatim from ``AC2BO`` so the size of that product can be computed
    without running the search -- ``tools/valsearch_scan.py`` uses this to find the
    over-cap ligand population for free. ``AC2BO`` calls it, so the two cannot drift.
    """
    valences_list_of_lists = []
    for i, (atomicNum, valence) in enumerate(zip(atoms, AC_valence)):
        # valence can't be smaller than number of neighbourgs
        possible_valence = [x for x in atomic_valence[atomicNum] if x >= valence]
        if atomicNum == 6 and valence == 1:
            possible_valence.remove(2)
        if atomicNum == 6 and not allow_carbenes and valence == 2:
            possible_valence.remove(2)
        if atomicNum == 6 and valence == 2:
            possible_valence.append(3)
        if atomicNum == 16 and valence == 1:
            possible_valence = [1, 2]

        if not possible_valence:
            logger.debug(
                "%s %s %s %s %s %s %s",
                "Valence of atom",
                i,
                "is",
                valence,
                "which bigger than allowed max",
                max(atomic_valence[atomicNum]),
                ". Stopping",
            )
            sys.exit()
        valences_list_of_lists.append(possible_valence)
    return valences_list_of_lists


def valence_combo_size(valences_list_of_lists, cap=_VALENCE_COMBO_CAP):
    """Product of the per-atom valence-list lengths, short-circuited above ``cap``.

    Mirrors ``AC2BO``'s own early break, so the returned number is exactly what the
    encoder compares against ``_VALENCE_COMBO_CAP``: the true product at or below the
    cap, and merely "greater than the cap" above it.
    """
    combo_size = 1
    for _vll in valences_list_of_lists:
        combo_size *= len(_vll)
        if combo_size > cap:
            break
    return combo_size


# The element grouping _ordered_valences sorts by, in its own priority order:
# O slowest-varying, then N, then C, then P, with S fastest. iter_ordered_valences
# reproduces this exactly, so the constant must stay in sync with the group extraction
# in _ordered_valences below (test_valence_order.py asserts the two orders are equal).
_HEURISTIC_ELEMENTS = (8, 7, 6, 15, 16)  # O, N, C, P, S

# Env lever, default OFF. When enabled, the over-cap branch of AC2BO enumerates candidate
# valence assignments in the SAME order _ordered_valences produces -- lazily, so it does not
# materialise the exponential product -- instead of the raw itertools.product order. Read
# only inside `if over_cap:`, so sub-cap ligands (99.8% of the corpus) are byte-identical by
# construction, exactly as OIN_VALENCE_FALLBACK_TRIES is. See docs/VALENCE_ORDER_v0.4.5.md.
_ORDERED_FALLBACK_ENV = "OIN_VALENCE_ORDERED_FALLBACK"


def iter_ordered_valences(valences_list_of_lists, atoms):
    """Yield candidates in exactly ``_ordered_valences``' order, lazily.

    ``_ordered_valences`` materialises the whole Cartesian product twice (once to score,
    once to sort), which is why ``AC2BO`` cannot use it above ``_VALENCE_COMBO_CAP``. But
    its sort key is ``(order_idx, full_tuple)`` where ``order_idx`` is the position in
    ``itertools.product(O_sums, N_sums, C_sums, P_sums, S_sums)`` -- i.e. a plain nested
    enumeration with O outermost and S innermost. Ties in ``order_idx`` can only differ at
    atoms that are none of O/N/C/P/S (the key pins every O/N/C/P/S valence), and
    ``sorted()`` breaks them lexicographically on the valence tuple in atom order, which
    is ascending per position.

    So six nested ``itertools.product``s over small per-element groups reproduce the
    order **exactly** in O(1) memory. ``tests/unit/test_valence_order.py`` asserts the two
    orders are element-for-element equal, so this is an equality claim, not an
    approximation.
    """
    groups = [[i for i, num in enumerate(atoms) if num == el] for el in _HEURISTIC_ELEMENTS]
    grouped = {i for group in groups for i in group}
    other = [i for i in range(len(atoms)) if i not in grouped]

    group_lists = [[valences_list_of_lists[i] for i in group] for group in groups]
    # The sorted() tie-break is ascending per position; the grouped elements keep their
    # own list order (which is atomic_valence's preference order, not ascending).
    other_lists = [sorted(valences_list_of_lists[i]) for i in other]

    out = [0] * len(atoms)

    def scatter(idxs, values):
        for i, value in zip(idxs, values):
            out[i] = value

    for o_vals in itertools.product(*group_lists[0]):
        scatter(groups[0], o_vals)
        for n_vals in itertools.product(*group_lists[1]):
            scatter(groups[1], n_vals)
            for c_vals in itertools.product(*group_lists[2]):
                scatter(groups[2], c_vals)
                for p_vals in itertools.product(*group_lists[3]):
                    scatter(groups[3], p_vals)
                    for s_vals in itertools.product(*group_lists[4]):
                        scatter(groups[4], s_vals)
                        for other_vals in itertools.product(*other_lists):
                            scatter(other, other_vals)
                            yield tuple(out)


_CHARGE_FILTER_ENV = "OIN_VALENCE_CHARGE_FILTER"


def charge_filter_supported(atoms):
    """Can the charge filter reason about this fragment at all?

    ``get_atomic_charge`` is called as ``get_atomic_charge(z, atomic_valence_electrons[z], v)``,
    and the 30 transition metals have an ``atomic_valence`` entry (``[20]``) but **no**
    ``atomic_valence_electrons`` entry -- so that lookup is a ``KeyError`` for them. The
    default path hits the same wall later (via ``charge_is_OK``), but the filter would hit it
    *before examining any candidate*, which would turn "crashes at candidate 1" into "crashes
    at candidate 0" and, in the exotic case where no candidate ever reaches ``charge_is_OK``,
    a working default into a crash. So an unsupported fragment declines the filter and takes
    the historical enumeration instead. Strict dominance again: the lever must never make
    anything worse.
    """
    return all(z in atomic_valence_electrons for z in atoms)


def _charge_feasible_suffix_counts(order_atoms, order_choices):
    """``counts[k]`` = the ``(dq, dpar)`` states reachable using variables ``k..n-1``.

    ``dq`` is the total ``get_atomic_charge`` contribution and ``dpar`` the parity of the
    valence sum. Both are additive over atoms, which is the whole reason this works.
    """
    n = len(order_atoms)
    counts = [None] * (n + 1)
    counts[n] = {(0, 0)}
    for k in range(n - 1, -1, -1):
        z = order_atoms[k]
        table = set()
        for val in order_choices[k]:
            dq = get_atomic_charge(z, atomic_valence_electrons[z], val)
            dpar = val & 1
            for q, par in counts[k + 1]:
                table.add((q + dq, par ^ dpar))
        counts[k] = table
    return counts


def iter_charge_feasible_valences(valences_list_of_lists, atoms, charge, ac_valence):
    """Yield only the candidates that can possibly satisfy ``AC2BO``'s own predicate.

    **This is a subsequence of ``itertools.product(*valences_list_of_lists)``, in the same
    relative order** -- it skips candidates, it never reorders them. That matters: every
    candidate a valid perception could use is still yielded, at its original relative
    position, so the *first valid candidate found is the same one an unbounded raw search
    would have found*. What changes is that the budget is no longer spent on candidates
    that provably cannot be valid.

    Why the skip is sound. A valid return from ``AC2BO`` forces ``BO.sum(axis=1) ==
    valences`` exactly: either ``UA`` is empty and ``BO = AC``, or ``valences_not_too_large``
    bounds every atom above by ``valences`` while ``(BO - AC).sum() == sum(DU)`` fixes the
    total, so every per-atom slack is zero. ``charge_is_OK`` therefore evaluates
    ``get_atomic_charge`` on the candidate's own valences, giving
    ``Q0 = sum_i get_atomic_charge(z_i, v_i)``; its only correction is ``Q += 2`` per
    trivalent single-bonded carbon and only while running below the target, so
    ``Q_final = Q0 + 2k`` with ``k >= 0``. Hence:

    * **C1** ``Q0 <= charge`` and ``charge - Q0`` even;
    * **C2** ``sum(valences) - sum(ac_valence)`` even, since every added bond raises
      ``BO.sum()`` by exactly 2.

    Both are necessary and both are additive, so a suffix DP prunes whole subtrees. The
    conditions are **not sufficient** -- ``get_BO``'s matching may still fail to saturate --
    so survivors are still handed to the real predicate by the caller.

    ``tests/unit/test_valence_order.py`` brute-forces small ligands and asserts that every
    candidate the real predicate accepts passes this filter, so the derivation above is
    pinned by measurement rather than by argument.
    """
    n = len(atoms)
    choices = [list(lst) for lst in valences_list_of_lists]
    counts = _charge_feasible_suffix_counts(atoms, choices)
    want_par = int(sum(ac_valence)) & 1
    targets = {
        (q, par)
        for q, par in counts[0]
        if par == want_par and q <= charge and (charge - q) % 2 == 0
    }
    if not targets:
        return  # provably no valid Lewis structure at this charge -- yield nothing

    out = [0] * n

    def walk(k, acc_q, acc_par):
        if k == n:
            yield tuple(out)
            return
        z = atoms[k]
        for val in choices[k]:
            q2 = acc_q + get_atomic_charge(z, atomic_valence_electrons[z], val)
            par2 = acc_par ^ (val & 1)
            reachable = counts[k + 1]
            if not any((tq - q2, tp ^ par2) in reachable for tq, tp in targets):
                continue
            out[k] = val
            yield from walk(k + 1, q2, par2)

    yield from walk(0, 0, 0)


def _ordered_valences(valences_list_of_lists, atoms):
    """Sort candidate valence assignments by the O/N/C/P/S grouping heuristic.

    Extracted verbatim from ``AC2BO``'s inner loop so the caller can bypass it when
    the Cartesian product of per-atom valences is too large to materialise. The
    heuristic only reorders which valid assignment the main loop finds first, so
    skipping it changes nothing a sub-cap (currently-encodable) ligand relies on.

    ``iter_ordered_valences`` produces the identical order without materialising
    anything; this function is kept as the sub-cap path (so that path stays
    byte-identical by construction) and as the reference the equality test compares to.
    """
    valences_list = itertools.product(*valences_list_of_lists)

    O_valences = [
        v_list for v_list, atomicNum in zip(valences_list_of_lists, atoms) if atomicNum == 8
    ]
    N_valences = [
        v_list for v_list, atomicNum in zip(valences_list_of_lists, atoms) if atomicNum == 7
    ]
    C_valences = [
        v_list for v_list, atomicNum in zip(valences_list_of_lists, atoms) if atomicNum == 6
    ]
    P_valences = [
        v_list for v_list, atomicNum in zip(valences_list_of_lists, atoms) if atomicNum == 15
    ]
    S_valences = [
        v_list for v_list, atomicNum in zip(valences_list_of_lists, atoms) if atomicNum == 16
    ]

    O_sums = []
    for v_list in itertools.product(*O_valences):
        O_sums.append(v_list)

    N_sums = []
    for v_list in itertools.product(*N_valences):
        N_sums.append(v_list)

    C_sums = []
    for v_list in itertools.product(*C_valences):
        C_sums.append(v_list)

    P_sums = []
    for v_list in itertools.product(*P_valences):
        P_sums.append(v_list)

    S_sums = []
    for v_list in itertools.product(*S_valences):
        S_sums.append(v_list)

    order_dict = dict()
    for i, v_list in enumerate(itertools.product(*[O_sums, N_sums, C_sums, P_sums, S_sums])):
        order_dict[v_list] = i

    valence_order_list = []
    for valence_list in valences_list:
        C_sum = []
        N_sum = []
        O_sum = []
        P_sum = []
        S_sum = []
        for v, atomicNum in zip(valence_list, atoms):
            if atomicNum == 6:
                C_sum.append(v)
            if atomicNum == 7:
                N_sum.append(v)
            if atomicNum == 8:
                O_sum.append(v)
            if atomicNum == 15:
                P_sum.append(v)
            if atomicNum == 16:
                S_sum.append(v)

        order_idx = order_dict[
            (tuple(O_sum), tuple(N_sum), tuple(C_sum), tuple(P_sum), tuple(S_sum))
        ]
        valence_order_list.append(order_idx)

    return [
        y
        for x, y in sorted(
            zip(valence_order_list, list(itertools.product(*valences_list_of_lists)))
        )
    ]


#: Process-local kill switch for ``OIN_CANONICAL_PERCEPTION``. ``get_tmc_mol`` sets it to
#: retry an encode in input order after the canonical perception produced a molecule the
#: encoder cannot use. Module-global rather than a parameter because the perception is
#: reached through several layers of the vendored xyz2mol call chain; the encoder already
#: forks for isolation where concurrency matters.
_SUPPRESS_CANONICAL_PERCEPTION = False


@contextlib.contextmanager
def suppress_canonical_perception():
    """Force input-order perception inside this block, whatever the env lever says."""
    global _SUPPRESS_CANONICAL_PERCEPTION
    previous = _SUPPRESS_CANONICAL_PERCEPTION
    _SUPPRESS_CANONICAL_PERCEPTION = True
    try:
        yield
    finally:
        _SUPPRESS_CANONICAL_PERCEPTION = previous


def _canonical_atom_permutation(AC, atoms):
    """Canonical atom ordering of the connectivity graph, or ``None``.

    Returns ``perm`` with ``perm[new_position] = old_index``: a total order on the atoms
    that depends only on the *graph* (elements + adjacency), never on the order the atoms
    happened to appear in the input file.

    The order is taken from ``_smilesAtomOutputOrder`` after writing the all-single-bond
    graph with ``MolToSmiles`` -- **not** from ``CanonicalRankAtoms(breakTies=True)``.
    That distinction is load-bearing and was measured: over 20 random renumberings of
    ``CC(N)=NC``, ``CanonicalRankAtoms(breakTies=True)`` returned a different ranking 18
    times (it settles ties between symmetry-equivalent atoms on the input index), while
    ``MolToSmiles`` returned one single string every time. The canonical *string* is the
    invariant RDKit actually guarantees, so the write order is what to build on. Two atoms
    swapped by a graph automorphism may still land in swapped positions, which is harmless:
    the perceived bond orders then come out as the automorphic image, and the automorphic
    image serializes to the same canonical SMILES.

    Bond orders are deliberately not used: at this point in perception they do not exist
    yet -- deciding them is what the caller is about to do. ``strict=False`` plus
    ``FastFindRings`` keeps an over-valent transition metal (six single bonds to a d-block
    centre) from failing the property cache.
    """
    try:
        n = len(atoms)
        if n == 0 or AC.shape[0] != n:
            return None
        rw = Chem.RWMol()
        for z in atoms:
            rw.AddAtom(Chem.Atom(int(z)))
        for i in range(n):
            for j in range(i + 1, n):
                if AC[i, j]:
                    rw.AddBond(i, j, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        Chem.MolToSmiles(mol)
        raw = mol.GetProp("_smilesAtomOutputOrder")
        perm = [int(x) for x in raw.strip("[]").rstrip(",").split(",") if x != ""]
        if sorted(perm) != list(range(n)):
            return None  # not a total order -- refuse rather than guess
        return perm
    except Exception:
        return None


def _valence_search_is_truncated(AC, atoms, allow_carbenes=True):
    """Whether ``_AC2BO_core`` will cap its valence walk, making its answer order-sensitive.

    Mirrors the per-atom ``possible_valence`` construction and the ``_VALENCE_COMBO_CAP``
    test at the top of ``_AC2BO_core``, and reads nothing but ``(AC, atoms)`` -- so the
    answer is the same however the atoms were numbered, which is what lets ``AC2BO`` use it
    as a renumbering-invariant switch. Errs toward ``True`` (do not canonicalize) on
    anything unexpected.
    """
    try:
        combo = 1
        for atomicNum, valence in zip(atoms, list(AC.sum(axis=1))):
            possible = [x for x in atomic_valence[atomicNum] if x >= valence]
            if atomicNum == 6 and valence == 1 and 2 in possible:
                possible.remove(2)
            if atomicNum == 6 and not allow_carbenes and valence == 2 and 2 in possible:
                possible.remove(2)
            if atomicNum == 6 and valence == 2:
                possible.append(3)
            if atomicNum == 16 and valence == 1:
                possible = [1, 2]
            if not possible:
                return True
            combo *= len(possible)
            if combo > _VALENCE_COMBO_CAP:
                return True
        return False
    except Exception:
        return True


def AC2BO(
    AC,
    atoms,
    charge,
    allow_charged_fragments=True,
    use_graph=True,
    allow_carbenes=True,
):
    """Bond orders from atomic connectivity, optionally made renumbering-invariant.

    ``_AC2BO_core`` returns, in its own words, "an arbitrary resonance form": it walks
    candidate per-atom valence assignments and returns the **first** one that validates,
    and both the walk order (``_ordered_valences`` -> ``itertools.product`` over per-atom
    lists in *input index* order) and the Kekule double-bond placement inside it
    (``get_UA_pairs`` -> ``nx.max_weight_matching``, whose result depends on the order
    edges were inserted) are functions of the input atom numbering. Permuting the atoms in
    the XYZ file therefore yields a genuinely different bond-order assignment for the same
    molecule -- an amidinate flips ``N=C(N-)`` to ``N-C(=N)``, a dioxime flips to its
    nitroso-enamine form -- and the emitted OIN string moves with it. Measured on the
    rotation/renumbering probe, this is the dominant cause of ``rdkit_canonical`` drift,
    and in about a fifth of drifting molecules it moves the comparison KEY too, not just
    the string. Re-serializing cannot repair it: ``MolToSmiles`` is faithful to whichever
    resonance form it is handed.

    ``OIN_CANONICAL_PERCEPTION`` (default OFF -> byte-identical) closes it by conjugation:
    relabel the atoms into a canonical order that depends only on the graph, do the whole
    perception there, and map the result back. That makes **every** index-order dependence
    inside the core -- the valence walk, the matching, the ``best_BO`` tie-break -- a
    function of the canonical labelling instead of the input numbering, in one place,
    rather than hardening each of them separately.

    **Canonicality never outranks perception quality.** The core's valence walk is
    *truncated* above ``_VALENCE_COMBO_CAP``: it tries only the first
    ``_VALENCE_FALLBACK_TRIES`` members of an exponential product and returns ``best_BO``
    if none validates. Which members those are depends on the atom order, so relabelling
    can move a valid assignment out of the searched window. Those molecules keep the
    un-permuted path, decided by ``_valence_search_is_truncated``, which reads only
    ``(AC, atoms)`` and so is itself renumbering-invariant.

    Even below the cap the walk can find a *different but equally valid* Lewis structure --
    on ``tests/fixtures/AGUFEN.xyz`` (a PPN counter-cation) the canonical order drops the
    total bond order 128 -> 126 and yields a pentavalent carbon that
    ``kekulize_safe_sanitize`` rejects, turning a working encode into an ``OINEncodeError``.
    ``utils.xyz2mol.get_tmc_mol`` catches that case: if the encode raises, it retries the
    whole perception with ``suppress_canonical_perception()``. A right answer that drifts
    beats a reproducible wrong one.

    Both fallbacks are decided from the CANONICAL result alone, never by comparing it with
    the un-permuted one. Comparing would silently re-import the order-dependence this
    closes, since the un-permuted result is a function of the input numbering -- that was
    tried, and it broke the ``NAXDOI`` invariance guard.

    Any failure falls through to the un-permuted path, so behaviour is unchanged whenever
    the canonical order cannot be computed.
    """
    plain = lambda: _AC2BO_core(  # noqa: E731
        AC,
        atoms,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph,
        allow_carbenes=allow_carbenes,
    )
    if _SUPPRESS_CANONICAL_PERCEPTION or not _lever_enabled("OIN_CANONICAL_PERCEPTION"):
        return plain()

    perm = _canonical_atom_permutation(AC, atoms)
    if perm is None or _valence_search_is_truncated(AC, atoms, allow_carbenes):
        return plain()

    idx = np.asarray(perm)
    BO_c, atomic_valence_electrons_out = _AC2BO_core(
        AC[np.ix_(idx, idx)],
        [atoms[i] for i in perm],
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph,
        allow_carbenes=allow_carbenes,
    )
    # atomic_valence_electrons is keyed by atomic NUMBER, not by atom index, so it needs
    # no un-permuting; only the BO matrix does.
    BO = np.zeros_like(BO_c)
    BO[np.ix_(idx, idx)] = BO_c
    return BO, atomic_valence_electrons_out


def _AC2BO_core(
    AC, atoms, charge, allow_charged_fragments=True, use_graph=True, allow_carbenes=True
):
    """Implemenation of algorithm shown in Figure 2.

    UA: unsaturated atoms

    DU: degree of unsaturation (u matrix in Figure)

    best_BO: Bcurr in Figure
    """
    global atomic_valence
    global atomic_valence_electrons

    # Anchor the single-slot memo (see _AC2BO_MEMO) on this AC. Entries survive a
    # repeat call on the same adjacency matrix -- which is exactly what
    # _select_lig_mol's charge/carbene ladder does -- and are dropped otherwise.
    _ac2bo_memo_anchor(AC)

    AC2BO_STATS["ac2bo_calls"] += 1

    # make a list of valences, e.g. for CO: [[4],[2,1]]
    # MERGE NOTE (encspeed x valsearch): both lanes changed this block. encspeed added the
    # memo anchor above; valsearch extracted the per-atom valence construction into
    # possible_valences() and added the call counter. Both are kept -- the memo and the
    # counter are orthogonal to the extraction, and the extracted function carries the same
    # logic the inline block had.
    AC_valence = list(AC.sum(axis=1))
    valences_list_of_lists = possible_valences(AC_valence, atoms, allow_carbenes=allow_carbenes)

    # convert [[4],[2,1]] to [[4,2],[4,1]]
    best_BO = AC.copy()

    # The O/N/C/P/S valence-ordering heuristic (in _ordered_valences) materialises the
    # full Cartesian product of per-atom valences to try chemically-sensible assignments
    # first. For a large conjugated ligand that product is exponential and AC2BO hangs
    # building it (the encode-fail timeout cohort: BENVOG et al. never return). Only sort
    # when the product is small enough to materialise cheaply; above the cap iterate the
    # lazy (unsorted) product, capped at _VALENCE_FALLBACK_TRIES candidates, so the main
    # loop early-returns on the first valid assignment (or fails fast) instead of hanging.
    # Every currently-encodable ligand is well under the cap and takes the byte-identical
    # sorted path.
    combo_size = valence_combo_size(valences_list_of_lists)

    over_cap = combo_size > _VALENCE_COMBO_CAP
    if over_cap:
        AC2BO_STATS["over_cap_calls"] += 1
        # Default: the raw product order. Two default-OFF levers change what the bounded
        # prefix contains -- OIN_VALENCE_ORDERED_FALLBACK reorders it into the sub-cap
        # heuristic's order (measured WORSE; see docs/VALENCE_ORDER_v0.4.5.md), and
        # OIN_VALENCE_CHARGE_FILTER keeps the order but drops candidates that provably
        # cannot be valid. The filter wins if both are set, since it subsumes the question.
        want_filter = _lever_enabled(_CHARGE_FILTER_ENV)
        if want_filter and not charge_filter_supported(atoms):
            AC2BO_STATS["over_cap_filter_unsupported"] += 1
            candidate_source = itertools.product(*valences_list_of_lists)
        elif want_filter:
            AC2BO_STATS["over_cap_filtered_calls"] += 1
            feasible = iter_charge_feasible_valences(
                valences_list_of_lists, atoms, charge, AC_valence
            )
            first = next(feasible, None)
            if first is None:
                # Provably nothing here can be valid, so there is no valid structure to
                # find -- but best_BO is still the value downstream judges, and it is built
                # from candidates the filter would drop. Fall back to the historical
                # enumeration so this case stays byte-identical to the default: the lever's
                # blast radius is then exactly "a guess becomes a real Lewis structure",
                # and never "a guess becomes a different guess". (The DP has already proved
                # the grind cannot succeed, so short-circuiting it is a separate, safe
                # optimisation -- see docs/VALENCE_ORDER_v0.4.5.md.)
                AC2BO_STATS["over_cap_infeasible"] += 1
                candidate_source = itertools.product(*valences_list_of_lists)
            else:
                candidate_source = itertools.chain([first], feasible)
        elif _lever_enabled(_ORDERED_FALLBACK_ENV):
            AC2BO_STATS["over_cap_ordered_calls"] += 1
            candidate_source = iter_ordered_valences(valences_list_of_lists, atoms)
        else:
            candidate_source = itertools.product(*valences_list_of_lists)
        sorted_valences_list = itertools.islice(candidate_source, _fallback_tries())
    else:
        # Same order, same first valid candidate, same best_BO -- generated instead of
        # materialised. `_ordered_valences` builds the full Cartesian product TWICE (once
        # to score, once inside `sorted(zip(...))`) plus a dict over the five-group
        # product, and the loop below measurably consumes only a prefix of it: on
        # `NOCGAN_comp_0` one call consumed **1** candidate of a materialised 20 736
        # (1.956 s of that call's 1.972 s), and on `UDITAD_comp_0` 9 of 11 664 (94.2%).
        # `iter_ordered_valences` yields the identical sequence in O(1) memory --
        # `tests/unit/test_valence_order.py` asserts element-for-element equality on 400
        # random configs AND on corpus-shaped ones, so this is dead-work removal against
        # an equality, not a heuristic swap. Byte-identity therefore holds whether the
        # loop early-returns (same first valid candidate) or exhausts (same candidates in
        # the same order, so the same `>=` tie-break survives into `best_BO`).
        sorted_valences_list = iter_ordered_valences(valences_list_of_lists, atoms)

    for valences in sorted_valences_list:  # valences_list:
        AC2BO_STATS["candidates"] += 1
        if over_cap:
            AC2BO_STATS["over_cap_candidates"] += 1
        UA, DU_from_AC = get_UA(valences, AC_valence)
        check_len = len(UA) == 0
        if check_len:
            check_bo = BO_is_OK(
                AC,
                AC,
                charge,
                DU_from_AC,
                atomic_valence_electrons,
                atoms,
                valences,
                allow_charged_fragments=allow_charged_fragments,
                allow_carbenes=allow_carbenes,
            )
        else:
            check_bo = None

        if check_len and check_bo:
            AC2BO_STATS["found_valid"] += 1
            if over_cap:
                AC2BO_STATS["over_cap_found_valid"] += 1
            return AC, atomic_valence_electrons

        UA_pairs_list = get_UA_pairs(UA, AC, DU_from_AC, use_graph=use_graph)
        for UA_pairs in UA_pairs_list:
            BO = get_BO(AC, UA, DU_from_AC, valences, UA_pairs, use_graph=use_graph)
            status = BO_is_OK(
                BO,
                AC,
                charge,
                DU_from_AC,
                atomic_valence_electrons,
                atoms,
                valences,
                allow_charged_fragments=allow_charged_fragments,
                allow_carbenes=allow_carbenes,
            )
            if status:
                AC2BO_STATS["found_valid"] += 1
                if over_cap:
                    AC2BO_STATS["over_cap_found_valid"] += 1
                return BO, atomic_valence_electrons
            # `charge_is_OK` was computed eagerly above this branch, then consumed only
            # here -- behind two cheaper predicates that already short-circuit, and even
            # when `status` had already returned. Evaluating it in place is a pure
            # dead-work removal: the value used is the same value, just not computed when
            # it cannot be read. (It is a pure predicate; its only side effect is a
            # DEBUG log line reachable only on the allow_carbenes=False arm, so the
            # emitted OIN is untouched.)
            elif (
                BO.sum() >= best_BO.sum()
                and valences_not_too_large(BO, valences)
                and charge_is_OK(
                    BO,
                    AC,
                    charge,
                    DU_from_AC,
                    atomic_valence_electrons,
                    atoms,
                    valences,
                    allow_charged_fragments=allow_charged_fragments,
                    allow_carbenes=allow_carbenes,
                )
            ):
                best_BO = BO.copy()
                if over_cap:
                    AC2BO_STATS["over_cap_best_bo_improved"] += 1

    if over_cap:
        AC2BO_STATS["over_cap_exhausted"] += 1
    return best_BO, atomic_valence_electrons


def AC2mol(
    mol,
    AC,
    atoms,
    charge,
    allow_charged_fragments=True,
    use_graph=True,
    use_atom_maps=True,
    allow_carbenes=True,
):
    """"""

    # convert AC matrix to bond order (BO) matrix
    BO, atomic_valence_electrons = AC2BO(
        AC,
        atoms,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph,
        allow_carbenes=allow_carbenes,
    )
    # add BO connectivity and charge info to mol object
    mol = BO2mol(
        mol,
        BO,
        atoms,
        atomic_valence_electrons,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_atom_maps=use_atom_maps,
    )

    # print(Chem.GetFormalCharge(mol), charge)
    # If charge is not correct don't return mol
    if Chem.GetFormalCharge(mol) != charge:
        return None

    # BO2mol returns an arbitrary resonance form. Let's make the rest

    # mols = rdchem.ResonanceMolSupplier(mol)
    # mols = [mol for mol in mols]
    # print(mols)

    return mol


def get_proto_mol(atoms):
    """"""
    mol = Chem.MolFromSmarts("[#" + str(atoms[0]) + "]")
    rwMol = Chem.RWMol(mol)
    for i in range(1, len(atoms)):
        a = Chem.Atom(atoms[i])
        rwMol.AddAtom(a)

    mol = rwMol.GetMol()

    return mol


def read_xyz_file(filename, look_for_charge=True):
    """"""
    atomic_symbols = []
    xyz_coordinates = []
    charge = 0

    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                int(line)
            elif line_number == 1:
                if "charge=" in line:
                    charge = int(line.split("=")[1])
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    atoms = [int_atom(atom) for atom in atomic_symbols]

    return atoms, charge, xyz_coordinates


def xyz2AC(atoms, xyz, charge, use_huckel=False, use_obabel=False):
    """Atoms and coordinates to atom connectivity (AC)

    Args:
        atoms - int atom types
        xyz - coordinates
        charge - molecule charge

    optional:
        use_huckel - Use Huckel method for atom connecitivty
        use_obabel - Use Opne Babel method for atom connectivity

    Returns:
        ac - atom connectivity matrix
        mol - rdkit molecule
    """
    if use_huckel:
        return xyz2AC_huckel(atoms, xyz, charge)
    elif use_obabel:
        return xyz2AC_obabel(atoms, xyz)
    else:
        return xyz2AC_vdW(atoms, xyz)


def xyz2AC_vdW(atoms, xyz):
    # Get mol template
    mol = get_proto_mol(atoms)

    # Set coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)

    AC = get_AC(mol)

    return AC, mol


def get_AC(mol, covalent_factor=1.3):
    """Generate adjacent matrix from atoms and coordinates.

    AC is a (num_atoms, num_atoms) matrix with 1 being covalent bond and 0 is not


    covalent_factor - 1.3 is an arbitrary factor

    Args:
        mol - rdkit molobj with 3D conformer

    optional
        covalent_factor - increase covalent bond length threshold with facto

    Returns:
        AC - adjacent matrix
    """
    # Calculate distance matrix
    dMat = Chem.Get3DDistanceMatrix(mol)

    pt = Chem.GetPeriodicTable()
    num_atoms = mol.GetNumAtoms()
    AC = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum()) * covalent_factor
        for j in range(i + 1, num_atoms):
            a_j = mol.GetAtomWithIdx(j)
            Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum()) * covalent_factor
            if dMat[i, j] <= Rcov_i + Rcov_j:
                AC[i, j] = 1
                AC[j, i] = 1

    return AC


def xyz2AC_huckel(atomicNumList, xyz, charge, tolerance=0.2):
    """Args.

        atomicNumList - atom type list
        xyz - coordinates
        charge - molecule charge

    optional
        tolerance - Huckel bond cutoff

    Returns:
        ac - atom connectivity
        mol - rdkit molecule
    """
    # print(charge)
    mol = get_proto_mol(atomicNumList)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)

    num_atoms = len(atomicNumList)
    AC = np.zeros((num_atoms, num_atoms)).astype(int)

    mol_huckel = Chem.Mol(mol)
    mol_huckel.GetAtomWithIdx(0).SetFormalCharge(charge)  # mol charge arbitrarily added to 1st atom

    passed, result = rdEHTTools.RunMol(mol_huckel)
    opop = result.GetReducedOverlapPopulationMatrix()
    tri = np.zeros((num_atoms, num_atoms))
    tri[np.tril(np.ones((num_atoms, num_atoms), dtype=bool))] = (
        opop  # lower triangular to square matrix
    )
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            pair_pop = abs(tri[j, i])
            if pair_pop >= tolerance:  # arbitry cutoff for bond. May need adjustment
                AC[i, j] = 1
                AC[j, i] = 1

    dMat = Chem.Get3DDistanceMatrix(mol)
    pt = Chem.GetPeriodicTable()

    # filter adjacency matrix if max valence is exceeded
    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        N_con = np.sum(AC[i, :])
        # print(a_i.GetAtomicNum(), N_con)
        while N_con > max(atomic_valence[a_i.GetAtomicNum()]):
            # print("removing longest bond")
            AC = remove_weakest_bond(mol, i, AC, dMat, pt)
            N_con = np.sum(AC[i, :])

    return AC, mol


def boron_cage_vertices(atoms, AC):
    """Indices of boron atoms sitting on a deltahedral (3c-2e) cage vertex.

    The motif is a **B-B-B triangle**: a boron with two boron neighbours that are
    themselves bonded. Every closo/nido borane and carborane vertex is on at least
    one such triangular face, and nothing else in this corpus is:

    * ``BPh4-`` / ``BF4-`` borates have **no B-B bond at all**;
    * a diborane or diboryl ``B-B`` bond has no third boron to close a triangle;
    * an ordinary boronic ester / boroxine ring alternates B-O-B, so again no
      B-B-B triangle.

    So this recognises the actual bonding motif rather than "molecule contains
    boron", which is what keeps the pruning exemption in ``xyz2AC_obabel`` from
    reaching any of the ~6,600 non-cage molecules.

    Args:
        atoms: list of atomic numbers, indexed like ``AC``.
        AC: adjacency matrix (symmetric, 0/1).

    Returns:
        set[int]: indices of cage-vertex borons (empty for every non-cage input).
    """
    borons = [i for i, z in enumerate(atoms) if z == 5]
    if len(borons) < 3:
        return set()
    bset = set(borons)
    b_nbrs = {i: [j for j in bset if j != i and AC[i][j]] for i in borons}
    cage = set()
    for i in borons:
        nb = b_nbrs[i]
        for a_i in range(len(nb)):
            for b_i in range(a_i + 1, len(nb)):
                if AC[nb[a_i]][nb[b_i]]:
                    cage.add(i)
                    cage.add(nb[a_i])
                    cage.add(nb[b_i])
                    break
            if i in cage:
                break
    return cage


def remove_weakest_bond(mol, atom_idx, AC, dMat, pt):
    extra_bond_lengths = []
    bond_atoms = np.nonzero(AC[atom_idx, :])[0]
    # print(bond_atoms)
    a_i = mol.GetAtomWithIdx(atom_idx)
    # print(a_i.GetAtomicNum())
    rcovi = pt.GetRcovalent(a_i.GetAtomicNum())
    for j in bond_atoms:
        # print(j)
        a_j = mol.GetAtomWithIdx(int(j))
        # print(a_j.GetAtomicNum())
        rcovj = pt.GetRcovalent(a_j.GetAtomicNum())
        extra_bond_length = dMat[atom_idx, j] - rcovj - rcovi
        extra_bond_lengths.append(extra_bond_length)

    longest_bond_index = bond_atoms[np.argmax(extra_bond_lengths)]
    AC[atom_idx, longest_bond_index] = 0
    AC[longest_bond_index, atom_idx] = 0

    return AC


def xyz2AC_obabel(atoms, xyz, tolerance=0.45):
    """Generate adjacent matrix from atoms and coordinates in a way similar to
    open babels.

    AC is a (num_atoms, num_atoms) matrix with 1 being covalent bond and 0 is not


    tolerance - 0.45Å is from the open babel paper

    Args:
        mol - rdkit molobj with 3D conformer

    optional
        tolerance - atoms connected if distance is shorter than sum of atomic
        radii + tolerance. If too many bonds to an atom; break longest bond

    Returns:
        AC - adjacency matrix
    """
    global atomic_valence
    # atomic_valence[8] = [2,1]
    # atomic_valence[7] = [3,2]
    atomic_valence[6] = [4, 2]

    # Get mol template
    mol = get_proto_mol(atoms)

    # Set coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)
    # Calculate distance matrix
    dMat = Chem.Get3DDistanceMatrix(mol)

    pt = Chem.GetPeriodicTable()
    num_atoms = mol.GetNumAtoms()
    AC = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum())
        for j in range(i + 1, num_atoms):
            a_j = mol.GetAtomWithIdx(j)
            Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum())
            if dMat[i, j] <= Rcov_i + Rcov_j + tolerance:
                AC[i, j] = 1
                AC[j, i] = 1

    # Filter the adjacency matrix where max valence is exceeded.
    #
    # TWO INDEPENDENT DEFECTS LIVE IN THIS ONE LOOP, and v0.4.5 fixed both. They compose:
    # one changes WHICH atoms are capped, the other changes the ORDER they are capped in.
    #
    # (a) OIN_BORON_CAGE -- WHICH. The loop deletes an atom's longest bonds while its
    #     connectivity exceeds max(atomic_valence[Z]). For boron that cap is 4, while a
    #     closo/nido deltahedral vertex has 5-6 neighbours -- so on a carborane it amputates
    #     7-19 B-B cage edges, shattering an intact, correctly-perceived cage into sub-cages
    #     plus loose [H]B. Measured on all 34 `boron_cluster` encode_fail molecules: the
    #     distance criterion above recovers textbook-exact topologies (o-carborane 10 B / 21
    #     edges, closo-B12 12/30, dicarbollide 9/18) and this loop then destroys them. The
    #     downstream "no 2c-2e Lewis structure for a 3c-2e cage" error is a CONSEQUENCE of
    #     this pruning, not an independent model limit. Exemption is scoped to element B in a
    #     B-B-B triangle motif and computed from the PRE-pruning AC, so it cannot be
    #     triggered by pruning itself; borates, diboranes and boroxines are untouched.
    #
    # (b) OIN_STABLE_METAL_AC -- ORDER. By default the loop runs in input atom order, so
    #     perception depends on how the XYZ happened to be numbered: capping atom i removes a
    #     bond, lowering some atom j's count, so whether j still needs capping depends on
    #     whether i came first. The distance pass above is order-free, making this loop the
    #     ONLY order-dependent step in AC perception. Measured on DUDREA_comp_0 (a Y
    #     borohydride): the bridging hydride is bonded to both B and Y, exceeding H's valence
    #     of 1 -- cap Y first and Y-H survives (degree 5, SPY), cap that H first and Y-H drops
    #     (degree 4, TET). Renumbering flipped the AC in 3 of 8 trials, flipping the emitted
    #     geometry tag. The replacement order is heaviest-element-first (so a metal claims its
    #     bridging hydrides before H's valence rule discards them), then a per-atom
    #     fingerprint of rotation- and permutation-invariant scalars only.
    #
    # See docs/BORON_CAGE_v0.4.5.md and docs/RENUMBERING_INSTABILITY_v0.4.5.md.
    # Both default OFF; with neither set this loop is byte-identical to pre-v0.4.5.
    exempt = set()
    if _lever_enabled("OIN_BORON_CAGE"):
        atomic_nums = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in range(num_atoms)]
        exempt = boron_cage_vertices(atomic_nums, AC)

    if _lever_enabled("OIN_STABLE_METAL_AC"):

        def _cap_key(i):
            nbrs = np.nonzero(AC[i, :])[0]
            dists = tuple(sorted(round(float(dMat[i, j]), 4) for j in nbrs))
            return (-mol.GetAtomWithIdx(int(i)).GetAtomicNum(), -len(nbrs), dists)

        cap_order = sorted(range(num_atoms), key=_cap_key)
    else:
        cap_order = range(num_atoms)

    for i in cap_order:
        if i in exempt:
            continue
        a_i = mol.GetAtomWithIdx(int(i))
        N_con = np.sum(AC[i, :])
        while N_con > max(atomic_valence[a_i.GetAtomicNum()]):
            # print("removing longest bond")
            AC = remove_weakest_bond(mol, int(i), AC, dMat, pt)
            N_con = np.sum(AC[i, :])

    # print(Chem.MolToSmiles(mol))

    return AC, mol


def chiral_stereo_check(mol):
    """Find and embed chiral information into the model based on the
    coordinates.

    Args:
        mol - rdkit molecule, with embeded conformer
    """
    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol, -1)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(mol, -1)

    return


def xyz2mol(
    atoms,
    coordinates,
    charge=0,
    allow_charged_fragments=True,
    use_graph=True,
    use_huckel=False,
    use_obabel=False,
    embed_chiral=True,
    use_atom_maps=True,
):
    """Generate a rdkit molobj from atoms, coordinates and a total_charge.

    Args:
        atoms - list of atom types (int)
        coordinates - 3xN Cartesian coordinates
        charge - total charge of the system (default: 0)

    optional:
        allow_charged_fragments - alternatively radicals are made
        use_graph - use graph (networkx)
        use_huckel - Use Huckel method for atom connectivity prediction
        embed_chiral - embed chiral information to the molecule

    Returns:
        mols - list of rdkit molobjects
    """
    # Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge,
    # and mol object with no connectivity information
    AC, mol = xyz2AC(atoms, coordinates, charge, use_huckel=use_huckel, use_obabel=use_obabel)
    # Convert AC to bond order matrix and add connectivity and charge info to
    # mol object
    new_mol = AC2mol(
        mol,
        AC,
        atoms,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph,
        use_atom_maps=use_atom_maps,
    )

    # Check for stereocenters and chiral centers
    if embed_chiral:
        chiral_stereo_check(new_mol)

    return new_mol


def canonicalize_smiles(structure_smiles):
    """Remove all structural info an atom mapping information."""
    mol = Chem.MolFromSmiles(structure_smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    return canonical_smiles


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(usage="%(prog)s [options] molecule.xyz")
    parser.add_argument("structure", metavar="structure", type=str)
    parser.add_argument("-s", "--sdf", action="store_true", help="Dump sdf file")
    parser.add_argument("--ignore-chiral", action="store_true", help="Ignore chiral centers")
    parser.add_argument(
        "--no-charged-fragments", action="store_true", help="Allow radicals to be made"
    )
    parser.add_argument(
        "--no-graph",
        action="store_true",
        help="Run xyz2mol without networkx dependencies",
    )

    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    parser.add_argument(
        "--use-huckel",
        action="store_true",
        help="Use Huckel method for atom connectivity",
    )
    parser.add_argument(
        "--use-obabel",
        action="store_true",
        help="Use Open Babel way of obtaining atom connectivity; recommended for radicals",
    )
    parser.add_argument(
        "-o",
        "--output-format",
        action="store",
        type=str,
        help="Output format [smiles,sdf] (default=sdf)",
    )
    parser.add_argument(
        "-c",
        "--charge",
        action="store",
        metavar="int",
        type=int,
        help="Total charge of the system",
    )
    parser.add_argument(
        "--use-atom-maps",
        action="store_true",
        help="Label atoms with map numbers according to their order in the .xyz file",
    )

    args = parser.parse_args()

    # read xyz file
    filename = args.structure

    # allow for charged fragments, alternatively radicals are made
    charged_fragments = not args.no_charged_fragments

    # quick is faster for large systems but requires networkx
    # if you don't want to install networkx set quick=False and
    # uncomment 'import networkx as nx' at the top of the file
    quick = not args.no_graph

    # chiral comment
    embed_chiral = not args.ignore_chiral

    # read atoms and coordinates. Try to find the charge
    atoms, charge, xyz_coordinates = read_xyz_file(filename)

    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    use_huckel = args.use_huckel

    use_obabel = args.use_obabel

    # if explicit charge from args, set it
    if args.charge is not None:
        charge = int(args.charge)

    use_atom_maps = args.use_atom_maps
    if not charged_fragments:
        atomic_valence[8] = [2, 1]
        atomic_valence[7] = [3, 2]
        atomic_valence[6] = [4, 2]

    # Get the molobjs
    mols = xyz2mol(
        atoms,
        xyz_coordinates,
        charge=charge,
        use_graph=quick,
        allow_charged_fragments=charged_fragments,
        embed_chiral=embed_chiral,
        use_huckel=use_huckel,
        use_obabel=use_obabel,
        use_atom_maps=use_atom_maps,
    )

    # Print output
    for mol in [mols]:
        if args.output_format == "sdf":
            txt = Chem.MolToMolBlock(mol)
            logger.debug(txt)

        else:
            # Canonical hack
            isomeric_smiles = not args.ignore_chiral
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles)
            # m = Chem.MolFromSmiles(smiles, sanitize=False)
            # smiles = Chem.MolToSmiles(m, isomericSmiles=isomeric_smiles)

            smiles = canonicalize_smiles(smiles)
