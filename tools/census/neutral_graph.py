"""The neutral ruler: compare two XYZ structures with NO perception and NO encoder.

WHY THIS EXISTS
---------------
Every accuracy verdict in this project routes through the encoder: ``XYZ -> E -> OIN -> G ->
XYZ' -> E -> OIN'``. A round-trip failure therefore cannot say whether ``G`` built the wrong
structure or ``E`` re-read a right one differently -- ``G`` always emits a new atom order and a
new conformer, so every round trip is also an invariance test of ``E``. This module judges
``XYZ'`` against ``XYZ`` directly. It imports NOTHING from ``oinsmiles`` on purpose: a ruler
that shares code with the thing it measures inherits its blind spots, and an instrument with no
``oinsmiles`` import cannot silently run against the wrong source tree.

WHAT IT COMPARES
----------------
1. **Graph.** Element-labelled heavy-atom adjacency from distances alone (covalent-radius
   cutoffs; no bond orders, no charges, no aromaticity). Hydrogens are counted per heavy atom
   rather than kept as nodes, because the generator is known to change H counts.
2. **Metal spheres.** Donor directions about each metal, superposed by proper and by improper
   rotation over every donor matching that preserves symmetry class and chelate linkage. This
   one test covers cis/trans, fac/mer and Delta/Lambda, with no slot numbering anywhere.
3. **Tetrahedral-type centres.** Signed volumes where all neighbours are symmetry-distinct.
   A mirror image negates every one of them.

Everything is compared through symmetry CLASSES, never through raw atom indices, so the verdict
does not depend on which of several graph isomorphisms the matcher happens to return.

KNOWN BLIND SPOTS (stated, not hidden)
--------------------------------------
* E/Z about double bonds -- the ruler has no bond orders, so it cannot tell a rotatable single
  bond from a locked double one.
* Atropisomers, helicenes and the planar chirality of a substituted eta ring (it is reduced to
  its centroid) -- no local centre to sign.
* Ring cis/trans where both ring paths are symmetry-equivalent -- the centre has tied neighbours
  and is skipped.
* Two symmetry classes that each carry mixed signs AND sit at equal graph distance from each
  other -- ``_relative`` separates the rest; ``n_mixed_classes`` is recorded per molecule.
A blind spot can only produce a false ``SAME``, never a false ``MIRROR``/``DIFFERENT``.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
from rdkit import Chem

_PT = Chem.GetPeriodicTable()

# fmt: off
#: Duplicated from ``oinsmiles.core.constants`` ON PURPOSE -- see the module docstring.
METALS = frozenset([
    21, 22, 23, 24, 25, 26, 27, 57, 28, 29, 30, 39, 40, 41,
    42, 43, 44, 45, 46, 47, 48, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
])
# fmt: on

#: non-metal pair bonded when d < COV_FACTOR * (r_i + r_j)  (the classic xyz2mol value)
COV_FACTOR = 1.30
#: metal pair bonded when d < r_M + r_X + METAL_SLACK  (``oin/coordination.py`` convention)
METAL_SLACK = 0.45
#: a heavy pair within this fraction of its cutoff is counted as "marginal" (reported only)
MARGIN = 0.06
#: a metal contact X is SECOND-SHELL (dropped) when it is covalently bonded to another contact Y
#: of the same metal that sits this much closer, as a fraction of its own cutoff. Measured need:
#: the generator builds M-L bonds ~0.35 A short, which drags chelate-backbone atoms (the C of a
#: carboxylate, the C next to an N donor) inside the metal cutoff -- 2194 "gained" C contacts on
#: passing molecules before this filter. Eta rings survive it: their carbons sit at near-equal
#: distance, so none is "much closer" than its ring neighbour.
SECOND_SHELL = 0.10

#: |normalised signed volume| below this is "flat": the sign carries no information
FLAT = 0.20
#: a centre must be at least this pyramidal in the INPUT to count as a stereo element
PYRAMIDAL = 0.35
#: donor-metal-donor angle bins (degrees)
TRANS_MIN, CIS_MAX = 150.0, 130.0


def parse_xyz(path):
    """``(symbols, coords)``; ``([], empty)`` on anything malformed."""
    try:
        lines = Path(path).read_text().splitlines()
        n = int(lines[0].strip())
        syms, xyz = [], []
        for ln in lines[2 : 2 + n]:
            p = ln.split()
            syms.append(p[0].capitalize())
            xyz.append([float(v) for v in p[1:4]])
        if len(syms) != n:
            return [], np.zeros((0, 3))
        return syms, np.asarray(xyz, dtype=float)
    except (OSError, IndexError, ValueError):
        return [], np.zeros((0, 3))


class NeutralGraph:
    """Heavy-atom graph of one structure.

    ``metal`` picks how generous the METAL contact definition is; ``cov`` scales the covalent
    cutoff. One fixed cutoff is the wrong model for a metal sphere: a semi-coordinated donor or
    a chelate-backbone atom sits near it and lands on either side in the crystal and in the
    generated geometry independently. So each structure is read three ways --

    * ``min``: firm contacts only (ratio < 0.90, second-shell atoms dropped out to 2 bonds)
    * ``std``: ratio < 1.00, second-shell atoms dropped out to 1 bond
    * ``max``: everything within ratio < 1.10, nothing dropped

    -- and ``compare`` accepts a match between ANY pair of readings, labelling anything other
    than std-vs-std as marginal rather than pretending the call was clean.
    """

    _METAL_CUT = {"min": 0.90, "std": 1.00, "max": 1.10}

    def __init__(self, symbols, coords, metal: str = "std", cov: float = 1.0):
        self.symbols, self.coords = list(symbols), np.asarray(coords, dtype=float)
        z = np.array([_PT.GetAtomicNumber(s) for s in self.symbols])
        r = np.array([_PT.GetRcovalent(int(a)) for a in z])
        is_m = np.isin(z, list(METALS))
        d = np.linalg.norm(self.coords[:, None, :] - self.coords[None, :, :], axis=2)
        either_m = is_m[:, None] | is_m[None, :]
        cut = np.where(
            either_m, r[:, None] + r[None, :] + METAL_SLACK, COV_FACTOR * (r[:, None] + r[None, :])
        )
        ratio = d / cut
        np.fill_diagonal(ratio, np.inf)

        self.heavy = [i for i in range(len(z)) if z[i] > 1]
        self.z, self.is_metal = z, is_m
        hpos = {a: k for k, a in enumerate(self.heavy)}
        heavy_mask = z > 1
        mcut = self._METAL_CUT[metal]
        iu = np.triu(np.ones_like(ratio, dtype=bool), 1) & heavy_mask[:, None] & heavy_mask[None, :]
        self.n_marginal = int((iu & (np.abs(ratio - 1.0) < MARGIN)).sum())
        cov_adj: dict = {}
        for a, b in np.argwhere(iu & ~either_m & (ratio < cov)):
            cov_adj.setdefault(int(a), set()).add(int(b))
            cov_adj.setdefault(int(b), set()).add(int(a))
        self.edges = {(hpos[a], hpos[b]) for a in cov_adj for b in cov_adj[a] if a < b}
        self.n_second_shell = 0
        for a, b in np.argwhere(iu & either_m & (ratio < mcut)):
            a, b = int(a), int(b)
            m, x = (a, b) if is_m[a] else (b, a)
            if metal != "max" and not is_m[x]:
                near = set(cov_adj.get(x, ()))
                if metal == "min":
                    near |= {w for y in list(near) for w in cov_adj.get(y, ())} - {x}
                if any(
                    ratio[m, y] < mcut and ratio[m, y] < ratio[m, x] - SECOND_SHELL for y in near
                ):
                    self.n_second_shell += 1
                    continue
            self.edges.add((hpos[min(a, b)], hpos[max(a, b)]))
        # each H goes to its relatively-nearest heavy atom; H with no partner is "free"
        self.n_h = [0] * len(self.heavy)
        self.n_free_h = 0
        hv = np.array(self.heavy, dtype=int)
        for i in np.where(z == 1)[0]:
            if len(hv) == 0:
                self.n_free_h += 1
                continue
            k = int(np.argmin(ratio[i, hv]))
            if ratio[i, hv[k]] < 1.0:
                self.n_h[k] += 1
            else:
                self.n_free_h += 1

        self.nbrs = [[] for _ in self.heavy]
        for a, b in self.edges:
            self.nbrs[a].append(b)
            self.nbrs[b].append(a)
        self.xyz = self.coords[hv] if len(hv) else np.zeros((0, 3))
        self.hz = [int(z[a]) for a in self.heavy]
        self.metals = [k for k, a in enumerate(self.heavy) if is_m[a]]

    def formula(self, with_h: bool = False) -> str:
        c = Counter(_PT.GetElementSymbol(a) for a in self.hz)
        if with_h:
            c["H"] = sum(self.n_h) + self.n_free_h
        return "".join(f"{k}{v}" for k, v in sorted(c.items()))

    def rdmol(self) -> Chem.Mol:
        """All-single-bond, unsanitized mol: a graph container, not a chemical claim."""
        rw = Chem.RWMol()
        for a in self.hz:
            at = Chem.Atom(a)
            at.SetNoImplicit(True)
            rw.AddAtom(at)
        for a, b in sorted(self.edges):
            rw.AddBond(a, b, Chem.BondType.SINGLE)
        m = rw.GetMol()
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        return m

    def ligand_components(self):
        """Connected components of the graph with the metals removed (= ligands)."""
        seen, comps = set(self.metals), []
        for s in range(len(self.heavy)):
            if s in seen:
                continue
            stack, comp = [s], []
            seen.add(s)
            while stack:
                u = stack.pop()
                comp.append(u)
                for v in self.nbrs[u]:
                    if v not in seen:
                        seen.add(v)
                        stack.append(v)
            comps.append(sorted(comp))
        return comps


def isomorphism(g_a: NeutralGraph, g_b: NeutralGraph, m_a=None, m_b=None):
    """Atom map ``a -> b`` (tuple) if the two heavy graphs are isomorphic, else ``None``."""
    if len(g_a.heavy) != len(g_b.heavy) or len(g_a.edges) != len(g_b.edges):
        return None
    sig = lambda g: sorted((g.hz[k], len(g.nbrs[k])) for k in range(len(g.heavy)))  # noqa: E731
    if sig(g_a) != sig(g_b):
        return None
    m_a = m_a or g_a.rdmol()
    m_b = m_b or g_b.rdmol()
    match = m_b.GetSubstructMatch(m_a, useChirality=False)
    # same atom count + same bond count + subgraph match  =>  isomorphism
    return tuple(match) if len(match) == len(g_a.heavy) else None


def _ligand_multisets_equal(g_a: NeutralGraph, g_b: NeutralGraph) -> bool:
    """Do the two structures have the same multiset of ligand skeletons (metals removed)?"""
    ca, cb = g_a.ligand_components(), g_b.ligand_components()
    if len(ca) != len(cb):
        return False

    def sub(g, comp):
        s = NeutralGraph.__new__(NeutralGraph)
        pos = {a: k for k, a in enumerate(comp)}
        s.heavy, s.hz = comp, [g.hz[a] for a in comp]
        s.edges = {(pos[a], pos[b]) for a, b in g.edges if a in pos and b in pos}
        s.nbrs = [[pos[v] for v in g.nbrs[a] if v in pos] for a in comp]
        return s

    left = [sub(g_b, c) for c in cb]
    for c in ca:
        sa = sub(g_a, c)
        hit = next((k for k, sb in enumerate(left) if isomorphism(sa, sb) is not None), None)
        if hit is None:
            return False
        left.pop(hit)
    return True


def _vol(c, p1, p2, p3) -> float:
    """Normalised signed volume of three neighbours about ``c``: 0 flat, ~0.77 tetrahedral."""
    v = [p - c for p in (p1, p2, p3)]
    n = [np.linalg.norm(x) for x in v]
    if min(n) < 1e-6:
        return 0.0
    return float(np.dot(v[0], np.cross(v[1], v[2])) / (n[0] * n[1] * n[2]))


def _sgn(x: float) -> int:
    return 0 if abs(x) < FLAT else (1 if x > 0 else -1)


def stereo_signature(g: NeutralGraph, cls, elements=None):
    """Signed volumes at tetrahedral-type centres, keyed by symmetry class (label-free).

    ``cls[k]`` is the class of heavy atom ``k``. For the generated structure the classes are
    the INPUT's classes carried over the isomorphism, so both sides speak one vocabulary.
    ``elements`` = the INPUT's list of stereo elements; when given, exactly those are evaluated
    (so a centre that went flat in the generated geometry reads 0 instead of vanishing).
    Returns ``(instances, elements)`` with ``instances[key] = [(sign, atom), ...]``.
    """
    find = elements is None
    elements = [] if find else elements
    if find:
        for c in range(len(g.heavy)):
            if g.is_metal[g.heavy[c]]:
                continue
            nb = sorted(g.nbrs[c], key=lambda a: cls[a])
            ncls = [cls[a] for a in nb]
            if len(set(ncls)) != len(ncls):
                continue
            if len(nb) == 4 or (len(nb) == 3 and _tet3_ok(g, c, nb)):
                key = ("tet", cls[c], tuple(ncls[:3]))
                if (
                    key not in elements
                    and abs(_vol(g.xyz[c], *(g.xyz[a] for a in nb[:3]))) >= PYRAMIDAL
                ):
                    elements.append(key)
    by_cls: dict = {}
    for c in range(len(g.heavy)):
        by_cls.setdefault(cls[c], []).append(c)
    inst: dict = {}
    for key in elements:
        _, ccls, ncls = key
        out = []
        for c in by_cls.get(ccls, []):
            pick = [next((a for a in g.nbrs[c] if cls[a] == k), None) for k in ncls]
            out.append((0 if None in pick else _sgn(_vol(g.xyz[c], *(g.xyz[a] for a in pick))), c))
        inst[key] = out or [(0, -1)]
    return inst, elements


def spheres(g: NeutralGraph, cls):
    """Per metal: ``(metal_class, [(donor_class, ligand_id, unit_vector), ...])``.

    Donors bonded to each other (an eta ring, a side-on alkene) collapse to ONE pseudo-donor at
    their centroid. Vectors are normalised: the generator's M-L bonds are ~0.35 A short, and a
    sphere comparison must not read a uniform compression as a different arrangement.
    """
    lig_of = {a: li for li, comp in enumerate(g.ligand_components()) for a in comp}
    out = []
    for m in g.metals:
        donors = [a for a in g.nbrs[m] if not g.is_metal[g.heavy[a]]]
        dset, seen, pd = set(donors), set(), []
        for d in donors:
            if d in seen:
                continue
            stack, grp = [d], []
            seen.add(d)
            while stack:
                u = stack.pop()
                grp.append(u)
                for v in g.nbrs[u]:
                    if v in dset and v not in seen:
                        seen.add(v)
                        stack.append(v)
            v = g.xyz[grp].mean(axis=0) - g.xyz[m]
            n = np.linalg.norm(v)
            pd.append(
                (
                    tuple(sorted(cls[a] for a in grp)),
                    lig_of.get(grp[0], -1),
                    v / n if n > 1e-6 else v,
                )
            )
        out.append((cls[m], pd))
    return out


#: cap on donor matchings tried per metal (a homoleptic CN-8 sphere is 8! = 40320)
MAX_MATCHINGS = 60000
#: unit-vector RMSD: below = the spheres superpose; a cis<->trans swap of two of six is ~0.8
SPHERE_FIT = 0.40
#: proper and improper fits must differ by at least this to call a handedness
SPHERE_GAP = 0.12


def _rot_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    """RMSD after the best PROPER rotation about the origin (the metal) taking ``q`` onto ``p``."""
    u, _, vt = np.linalg.svd(p.T @ q)
    dmat = np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt)) or 1.0])
    return float(np.sqrt(np.mean(np.sum((p - q @ (u @ dmat @ vt).T) ** 2, axis=1))))


def sphere_fit(a, b):
    """Best (proper, improper) RMSD over donor matchings preserving class AND chelate linkage.

    Whole-molecule superposition confuses conformation with configuration; the first
    coordination sphere does not flex that way, so rigid superposition is sound HERE. The
    matching is enumerated rather than guessed, so the answer cannot depend on which graph
    isomorphism the matcher returned. Returns ``(proper, improper, n_matchings, capped)``.
    """
    if len(a) != len(b) or sorted(x[0] for x in a) != sorted(x[0] for x in b):
        return np.inf, np.inf, 0, False
    if not a:
        return 0.0, 0.0, 1, False
    pa = np.array([x[2] for x in a])
    pb = np.array([x[2] for x in b])
    mirror = pb * np.array([1.0, 1.0, -1.0])
    best = [np.inf, np.inf]
    count = [0]
    order = sorted(
        range(len(a)), key=lambda i: sum(x[0] == a[i][0] for x in a)
    )  # rarest class first

    def rec(k, used, lmap, rmap, assign):
        if count[0] >= MAX_MATCHINGS:
            return
        if k == len(order):
            count[0] += 1
            idx = [assign[i] for i in range(len(a))]
            best[0] = min(best[0], _rot_rmsd(pa, pb[idx]))
            best[1] = min(best[1], _rot_rmsd(pa, mirror[idx]))
            return
        i = order[k]
        for j in range(len(b)):
            if j in used or b[j][0] != a[i][0]:
                continue
            la, lb = a[i][1], b[j][1]
            if lmap.get(la, lb) != lb or rmap.get(lb, la) != la:
                continue
            new_l, new_r = la not in lmap, lb not in rmap
            lmap[la], rmap[lb] = lb, la
            assign[i] = j
            rec(k + 1, used | {j}, lmap, rmap, assign)
            if new_l:
                del lmap[la]
            if new_r:
                del rmap[lb]

    rec(0, frozenset(), {}, {}, {})
    return best[0], best[1], count[0], count[0] >= MAX_MATCHINGS


def sphere_verdict(sa, sb):
    """``(state, detail)`` for one metal: same / inverted / achiral / arr_diff / undecided."""
    p, q, n, capped = sphere_fit(sa, sb)
    c, _, _, _ = sphere_fit(sa, [(x[0], x[1], x[2] * np.array([1.0, 1.0, -1.0])) for x in sa])
    detail = {
        "proper": round(p, 3),
        "improper": round(q, 3),
        "self_mirror": round(c, 3),
        "matchings": n,
    }
    if capped or n == 0:
        return "undecided", detail
    if min(p, q) > SPHERE_FIT:
        return "arr_diff", detail
    if c < SPHERE_FIT - SPHERE_GAP or abs(p - q) < SPHERE_GAP:
        return (
            "achiral",
            detail,
        )  # the sphere superposes on its own mirror: no handedness to get wrong
    return ("same" if p < q else "inverted"), detail


def _absolute(inst) -> dict:
    """Per-element sorted sign multiset. A reflection negates every entry."""
    return {k: tuple(sorted(x for x, _ in v)) for k, v in inst.items()}


def _relative(g: NeutralGraph, inst) -> Counter:
    """RELATIVE configuration between stereo centres: (keys, sign product, graph distance).

    The per-class multisets cannot tell ``(A+ B- | C+ D-)`` from ``(A+ B- | C- D+)`` when A~B
    and C~D are symmetry-equivalent -- two genuinely different diastereomers. The sign PRODUCT
    of a pair, tagged with how far apart the two centres sit in the graph, can. It is invariant
    under reflection (both signs flip), so it separates diastereomers without confusing them
    with enantiomers. Residual blind spot: pairs at equal graph distance.
    """
    flat = [(k, x, a) for k, v in inst.items() for x, a in v if x != 0 and a >= 0]
    dist = {}
    for a in {a for _, _, a in flat}:
        d, frontier = {a: 0}, [a]
        while frontier:
            nxt = []
            for u in frontier:
                for w in g.nbrs[u]:
                    if w not in d:
                        d[w] = d[u] + 1
                        nxt.append(w)
            frontier = nxt
        dist[a] = d
    rel: Counter = Counter()
    for (k1, s1, a1), (k2, s2, a2) in combinations(flat, 2):
        rel[(tuple(sorted((repr(k1), repr(k2)))), s1 * s2, dist[a1].get(a2, -1))] += 1
    return rel


def _tet3_ok(g: NeutralGraph, c: int, nb) -> bool:
    """A 3-heavy-neighbour centre is a stereo candidate only if its 4th site is fixed.

    One H (C-H, metal-bound N-H) or a configurationally stable lone pair (P, As, S, Se). A free
    amine N inverts at room temperature and a bare 3-coordinate C is sp2 -- both are skipped.
    """
    z = g.hz[c]
    if g.n_h[c] == 1:
        return True
    if g.n_h[c] == 0 and z in (15, 33, 16, 34):
        return True
    return False


_LADDER = [(mi, mg, cg) for cg in (1.0, 0.88) for mi, mg in (
    ("std", "std"), ("std", "min"), ("min", "std"), ("min", "min"), ("std", "max"), ("max", "std"), ("max", "max"),
    ("min", "max"), ("max", "min"),
)]  # fmt: skip


def compare(inp_path, gen_path=None, *, gen=None):
    """Full ruler verdict for one molecule. ``gen=(symbols, coords)`` overrides ``gen_path``."""
    s_in, x_in = parse_xyz(inp_path)
    s_g, x_g = gen if gen is not None else parse_xyz(gen_path)
    if not s_in or not s_g:
        return {"graph": "UNREADABLE", "stereo": "NA"}
    g_std = NeutralGraph(s_in, x_in)
    out = {
        "n_heavy_in": len(g_std.heavy),
        "n_atoms_in": len(s_in),
        "n_atoms_gen": len(s_g),
        "n_metals": len(g_std.metals),
        "marginal_in": g_std.n_marginal,
    }
    # Readings are tried from the cleanest outward; the FIRST match wins and is recorded, so a
    # marginal verdict always says exactly which concession bought it. ``cov=0.88`` is the clash
    # rescue: a non-bonded heavy-atom pair squeezed inside bonding range (C...C < 1.98 A).
    ins, gens = {"std": g_std}, {}
    g_in = m_in = g_gen = amap = None
    graph, combo = "NON_ISO", None
    for mi, mg, cg in _LADDER:
        if mi not in ins:
            ins[mi] = NeutralGraph(s_in, x_in, metal=mi)
        if (mg, cg) not in gens:
            gens[(mg, cg)] = NeutralGraph(s_g, x_g, metal=mg, cov=cg)
        if not hasattr(ins[mi], "_mol"):
            ins[mi]._mol = ins[mi].rdmol()
        amap = isomorphism(ins[mi], gens[(mg, cg)], m_a=ins[mi]._mol)
        if amap is not None:
            g_in, m_in, g_gen, combo = ins[mi], ins[mi]._mol, gens[(mg, cg)], f"{mi}/{mg}/{cg}"
            graph = (
                "ISO"
                if (mi, mg, cg) == ("std", "std", 1.0)
                else ("ISO_CLASH" if cg != 1.0 else "ISO_MARGINAL")
            )
            break
    if amap is None:
        g_in, m_in, g_gen = g_std, None, gens[("std", 1.0)]
    out["combo"] = combo
    out["marginal_gen"] = g_gen.n_marginal
    out["formula_in"], out["formula_gen"] = g_in.formula(), g_gen.formula()
    out["h_in"], out["h_gen"] = sum(g_in.n_h) + g_in.n_free_h, sum(g_gen.n_h) + g_gen.n_free_h
    if amap is None:
        if out["formula_in"] != out["formula_gen"]:
            graph = "FORMULA_DIFF"
        elif _ligand_multisets_equal(g_in, g_gen):
            graph = "SPHERE_DIFF"  # same ligands, attached to the metal differently
        else:
            graph = "LIGAND_DIFF"  # a covalent bond inside a ligand was made or broken
        out["graph"], out["stereo"] = graph, "NA"
        out["cn_in"] = sorted(len(g_in.nbrs[m]) for m in g_in.metals)
        out["cn_gen"] = sorted(len(g_gen.nbrs[m]) for m in g_gen.metals)
        return out
    out["graph"] = graph

    cls_in = list(
        Chem.CanonicalRankAtoms(
            m_in, breakTies=False, includeChirality=False, includeIsotopes=False
        )
    )
    cls_gen = [0] * len(cls_in)
    for a, b in enumerate(amap):
        cls_gen[b] = cls_in[a]
    # H decoration, per class (an index-wise compare would depend on WHICH isomorphism we got)
    h_in = Counter((cls_in[k], g_in.n_h[k]) for k in range(len(cls_in)))
    h_gen = Counter((cls_gen[k], g_gen.n_h[k]) for k in range(len(cls_gen)))
    out["h_decoration"] = "SAME" if h_in == h_gen else "DIFF"

    sg_in, elements = stereo_signature(g_in, cls_in)
    sg_gen, _ = stereo_signature(g_gen, cls_gen, elements=elements)
    out["n_tet"] = len(elements)

    # --- metal spheres: pair each input metal with the best same-class generated metal --------
    sp_in, sp_gen = spheres(g_in, cls_in), spheres(g_gen, cls_gen)
    states, details, left = [], [], list(range(len(sp_gen)))
    for mc, sa in sp_in:
        cand = [(sphere_verdict(sa, sp_gen[j][1]), j) for j in left if sp_gen[j][0] == mc]
        if not cand:
            states.append("undecided")
            continue
        (st, det), j = min(cand, key=lambda t: min(t[0][1]["proper"], t[0][1]["improper"]))
        left.remove(j)
        states.append(st)
        details.append(det)
    out["sphere_states"], out["sphere_detail"] = states, details
    out["sphere"] = (
        "ARRANGEMENT_DIFF"
        if "arr_diff" in states
        else ("UNDECIDED" if "undecided" in states else "SAME")
    )

    # --- tetrahedral-type centres -----------------------------------------------------------
    ab_in, ab_gen = _absolute(sg_in), _absolute(sg_gen)
    neg = {k: tuple(sorted(-x for x in v)) for k, v in ab_in.items()}
    flat_keys = {k for k, v in ab_gen.items() if v.count(0) > ab_in[k].count(0)}
    out["n_flattened"] = len(flat_keys)
    out["n_mixed_classes"] = sum(1 for v in ab_in.values() if 1 in v and -1 in v)
    for k in ab_in:
        if k in flat_keys:
            continue
        if ab_in[k] == neg[k]:
            states.append("achiral")  # balanced class (meso-like): its mirror is itself
        else:
            states.append(
                "same" if ab_gen[k] == ab_in[k] else ("inverted" if ab_gen[k] == neg[k] else "diff")
            )
    rel_same = bool(flat_keys) or _relative(g_in, sg_in) == _relative(g_gen, sg_gen)

    chiral = [x for x in states if x in ("same", "inverted")]
    out["input_chiral_by_ruler"] = bool(chiral) or "diff" in states
    part = "_PARTIAL" if flat_keys else ""
    if "arr_diff" in states or "diff" in states or not rel_same:
        stereo = "DIFFERENT"
    elif not chiral:
        stereo = (
            "UNDECIDED" if ("undecided" in states or (flat_keys and not elements)) else "NO_STEREO"
        )
    elif all(x == "same" for x in chiral):
        stereo = "SAME" + part
    elif all(x == "inverted" for x in chiral):
        stereo = "MIRROR" + part
    else:
        stereo = "DIFFERENT"  # some elements inverted, others kept: a diastereomer
    if stereo in ("SAME", "MIRROR") and "undecided" in states:
        stereo += "_PARTIAL"
    out["stereo"] = stereo
    return out
