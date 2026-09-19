"""I3 of the census: is the STRING sufficient? Parse-back, collision scan, perception flags.

    python tools/census/string_sufficiency.py parseback --control delete    # must be 0% ISO
    python tools/census/string_sufficiency.py parseback --control reattach  # must be 0% ISO
    python tools/census/string_sufficiency.py parseback                     # smiles_1 vs input
    python tools/census/string_sufficiency.py parseback --side gen          # smiles_2 vs generated
    python tools/census/string_sufficiency.py collide                       # natural twins
    python tools/census/string_sufficiency.py pflags --charge-probe         # perception flags
    python tools/census/string_sufficiency.py summary                       # crosstabs by bucket

THREE QUESTIONS, NO EMBEDDING
-----------------------------
(a) **Parse-back.** Read ``smiles_1`` back into a graph -- at the PARSER level (the inline
    handler's fragments and slot markers, nothing else) and at the ADAPTER level (the ligand
    SMILES the generator actually hands to MetalloGen, after its kekulisation and donor-H
    reconciliation) -- and compare each to the input's neutral graph from ``neutral_graph.py``
    (heavy-atom isomorphism; H counts per symmetry class reported beside it). No 3D on the
    string side, so this isolates the serializer/parser from the embedding. It is the only
    instrument that can say anything about the 266 ``hard_fail`` that produced no structure.
(b) **Collision scan.** Over the 5,000 inputs, cluster by heavy-graph isomorphism (the I1
    ruler's graph), then within a cluster settle the stereo relation with the ruler itself
    (input vs input). Same string + ruler MIRROR/DIFFERENT = a collision (E non-injective in the
    wild). Different string + ruler SAME = non-canonical across independent crystal structures.
    The cat/photo twins are byte-identical files (1,033/1,033, checked) and E is deterministic
    (C2), so the twin scan reduces to ``_comp_n`` siblings and refcode-level duplicates.
(c) **Perception flags.** ``XYZToSMILES.convert`` hardcodes ``charge = 0`` and never reads the
    ``Charge:`` field of the comment line; perception then sets the metal to
    ``0 - sum(ligand charges)``. So for every stated-charged input the perceived oxidation
    state is off by exactly the stated charge. ``pflags`` records the stated charge, the
    perceived metal charge, whether it is in ``ALLOWED_OXIDATION_STATES`` (defined in
    ``perception_tmc.py`` and used nowhere), the ligand charge sum and soft valence anomalies.
    ``--charge-probe`` re-encodes each stated-charged input with the charge HONOURED (the
    perception entry point is patched in the worker) and asks whether the string changes.

WHAT A BROKEN VERSION WOULD PRINT
---------------------------------
A parse-back that ignored the slot markers would still read ISO on every metal-free ligand
graph, so the ``reattach`` control moves ONE metal edge to a neighbouring atom and must read
non-ISO; ``delete`` drops a whole ligand and must read non-ISO. Both are run and read before
the real numbers. Both FAILED on the full cohort the first time, and both failures were the
graph's stated blind spots: a heavy-atom graph cannot see a hydride ``[H]{n}`` leave (79 read
ISO, all 79 read H-DIFF), and moving the metal to an automorphic neighbour in a bond-order-free
graph (unsubstituted phenyl, N2, tetrazolate, triazolyl N3/N4) is not a change (11 read ISO).
The controls now exclude exactly those cases; one residual (``XEWSID``) matches only at the
``max`` metal reading. For the charge probe the control is inside every row: the probe first encodes
at charge 0 through the patched path and that string must equal ``smiles_1`` byte for byte;
the charged encode must then move the perceived metal charge by exactly the stated charge, or
the patch did not fire and the row is void. Every table prints its denominator.

Paths default to the MAIN checkout (absolute): the dataset and the sweep results are gitignored
and do not exist in a worktree. ``oinsmiles.__file__`` is printed before any worker forks.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import multiprocessing as mp
import os
import re
import signal
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling module only
import neutral_graph as ng  # noqa: E402

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
DEF_COHORT = MAIN / "cohort-v0.4.5-5k"
DEF_SWEEP = MAIN / "results-v0.4.14-sweep"
DEF_OUT = MAIN / "results-census"

ISO_SET = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
CHARGE_RE = re.compile(r"Charge:\s*(-?\d+)?")
PASS_STEREO = ("SAME", "NO_STEREO", "SAME_PARTIAL")


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _comment_charge(path) -> "int | None":
    try:
        line = Path(path).read_text().splitlines()[1]
    except (OSError, IndexError):
        return None
    m = CHARGE_RE.search(line)
    if not m or m.group(1) is None:
        return None
    return int(m.group(1))


def _load_rows(sweep):
    rows = json.loads((sweep / "bucket_report_honest.json").read_text())
    return {r["molecule"]: r for r in rows}


def _check_cohort(cohort):
    mols = sorted(p.stem for p in cohort.glob("*.xyz"))
    dangling = [m for m in mols if not (cohort / f"{m}.xyz").exists()]
    if dangling:
        sys.exit(f"ABORT: {len(dangling)} dangling cohort symlinks -- restore the dataset first")
    return mols


# ----------------------------------------------------------------------------------------------
# a graph built from SMILES fragments, shaped like the ruler's NeutralGraph (heavy atoms only)
# ----------------------------------------------------------------------------------------------


def graph_from_ligands(metal_z: int, lig_mols, donor_edges) -> ng.NeutralGraph:
    """A ruler-shaped heavy-atom graph: one metal + ligand mols + metal->donor edges.

    ``donor_edges`` is a list of ``(ligand_index, atom_index)``. Hydrogen atoms that appear as
    nodes (``[H]`` hydride fragments, explicit ``[H]`` in a ligand) are folded into the H count
    of the atom they bond to -- the metal for a hydride -- exactly as the ruler folds them.
    """
    g = ng.NeutralGraph.__new__(ng.NeutralGraph)
    hz, n_h, edges = [metal_z], [0], set()
    free_h = 0
    donor_set = set(donor_edges)
    maps = []
    for li, mol in enumerate(lig_mols):
        idx = {}
        for a in mol.GetAtoms():
            if a.GetAtomicNum() > 1:
                idx[a.GetIdx()] = len(hz)
                hz.append(a.GetAtomicNum())
                n_h.append(a.GetTotalNumHs())
        for a in mol.GetAtoms():
            if a.GetAtomicNum() != 1:
                continue
            heavy_nb = [n for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
            if heavy_nb:
                n_h[idx[heavy_nb[0].GetIdx()]] += 1
            elif (li, a.GetIdx()) in donor_set:
                n_h[0] += 1  # hydride on the metal
            else:
                free_h += 1
        for b in mol.GetBonds():
            p, q = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if p in idx and q in idx:
                edges.add((min(idx[p], idx[q]), max(idx[p], idx[q])))
        maps.append(idx)
    for li, ai in donor_edges:
        if li < len(maps) and ai in maps[li]:
            edges.add((0, maps[li][ai]))
    n = len(hz)
    g.heavy, g.hz, g.edges, g.n_h, g.n_free_h = list(range(n)), hz, edges, n_h, free_h
    g.z = np.array(hz)
    g.is_metal = np.array([k == 0 for k in range(n)])
    g.metals = [0]
    g.nbrs = [[] for _ in range(n)]
    for a, b in edges:
        g.nbrs[a].append(b)
        g.nbrs[b].append(a)
    g.n_marginal = g.n_second_shell = 0
    g.xyz = np.zeros((n, 3))
    g.symbols, g.coords = [], np.zeros((0, 3))
    return g


def compare_to_input(s_in, x_in, gp: ng.NeutralGraph) -> dict:
    """Heavy-graph verdict of a SMILES-built graph against an XYZ input (metal-reading ladder)."""
    # Same ladder as the ruler's XYZ side: three metal readings, then the clash rescue
    # (``cov=0.88``: a non-bonded heavy pair squeezed inside bonding range, which generated
    # geometries do and crystal structures do not). First match wins and is labelled.
    ins = {}
    mp_ = gp.rdmol()
    verdict, amap, g_in = "NON_ISO", None, None
    for cov in (1.0, 0.88):
        for reading in ("std", "min", "max"):
            ins[(reading, cov)] = ng.NeutralGraph(s_in, x_in, metal=reading, cov=cov)
            ins[(reading, cov)]._mol = ins[(reading, cov)].rdmol()
            amap = ng.isomorphism(ins[(reading, cov)], gp, m_a=ins[(reading, cov)]._mol, m_b=mp_)
            if amap is not None:
                g_in = ins[(reading, cov)]
                verdict = (
                    "ISO_CLASH" if cov != 1.0 else "ISO" if reading == "std" else "ISO_MARGINAL"
                )
                break
        if amap is not None:
            break
    if amap is None:
        g_in = ins[("std", 1.0)]
    pt = Chem.GetPeriodicTable()
    out = {
        "graph": verdict,
        "n_heavy_in": len(g_in.heavy),
        "n_heavy_str": len(gp.heavy),
        "h_in": sum(g_in.n_h) + g_in.n_free_h,
        "h_str": sum(gp.n_h) + gp.n_free_h,
        "cn_in": sorted(len(g_in.nbrs[m]) for m in g_in.metals),
        "cn_str": [len(gp.nbrs[0])],
        "n_metals_in": len(g_in.metals),
        "marginal_in": g_in.n_marginal,
        "second_shell_in": g_in.n_second_shell,
        "donors_in": sorted(
            pt.GetElementSymbol(g_in.hz[a]) for m in g_in.metals for a in g_in.nbrs[m]
        ),
        "donors_str": sorted(pt.GetElementSymbol(gp.hz[a]) for a in gp.nbrs[0]),
        "n_components_in": len(g_in.ligand_components()),
        "n_components_str": len(gp.ligand_components()),
    }
    if amap is None:
        if g_in.formula() != gp.formula():
            out["graph"] = "FORMULA_DIFF"
        elif ng._ligand_multisets_equal(g_in, gp):
            out["graph"] = "SPHERE_DIFF"
        else:
            out["graph"] = "LIGAND_DIFF"
        return out
    cls = list(
        Chem.CanonicalRankAtoms(
            g_in._mol, breakTies=False, includeChirality=False, includeIsotopes=False
        )
    )
    h_in = Counter((cls[k], g_in.n_h[k]) for k in range(len(cls)))
    h_str = Counter((cls[a], gp.n_h[b]) for a, b in enumerate(amap))
    out["h_decoration"] = "SAME" if h_in == h_str else "DIFF"
    out["h_diff_atoms"] = sum((h_str - h_in).values())
    return out


def _parser_level(oin: str):
    """Fragments + slot assignments straight from the inline handler; nothing else touched."""
    from oinsmiles.oin.inline import OINInlineHandler

    m = OINInlineHandler.METAL_REGEX.search(oin)
    if not m:
        raise ValueError("no [M_GEO] token")
    metal_z = Chem.GetPeriodicTable().GetAtomicNumber(m.group(1))
    clean, geo, assigns = OINInlineHandler.parse_inline_string(oin)
    frags = clean.split(".")[1:]
    ligs = []
    for f in frags:
        mol = Chem.MolFromSmiles(f, sanitize=False)
        if mol is None:
            raise ValueError(f"fragment unparsable: {f[:60]}")
        mol.UpdatePropertyCache(strict=False)
        ligs.append(mol)
    donors = [(a.lig_rank - 1, a.atom_idx) for a in assigns]
    slotted = {a.lig_rank - 1 for a in assigns}
    info = {
        "geo": geo,
        "n_lig": len(ligs),
        "n_unslotted": sum(1 for i in range(len(ligs)) if i not in slotted),
        "n_slots": len({a.slot for a in assigns}),
        "n_donor_atoms": len(set(donors)),
        "str_charge": sum(a.GetFormalCharge() for mol in ligs for a in mol.GetAtoms()),
    }
    # soft valence check on the ligands AS WRITTEN: which fragments RDKit refuses to sanitize
    bad = []
    for f in frags:
        try:
            Chem.SanitizeMol(Chem.MolFromSmiles(f, sanitize=False))
        except Exception as exc:  # AtomValenceException / KekulizeException
            bad.append(type(exc).__name__)
    info["str_valence_fail"] = len(bad)
    info["str_valence_exc"] = sorted(set(bad))
    return metal_z, ligs, donors, info


def _adapter_level(oin: str, metal_z: int):
    """The ligand SMILES the generator hands to MetalloGen, read back as a graph."""
    from oinsmiles.generation.metallogen_adapter import _prepare_ligand_fragments
    from oinsmiles.generation.oin_parser import OINParser

    parsed = OINParser().parse(oin)
    _, specs, _ = _prepare_ligand_fragments(parsed)
    ligs, donors, unsan = [], [], 0
    for li, (sm, _w) in enumerate(specs):
        mol = Chem.MolFromSmiles(sm, sanitize=False)
        if mol is None:
            raise ValueError(f"adapter fragment unparsable: {sm[:60]}")
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            unsan += 1
            mol = Chem.MolFromSmiles(sm, sanitize=False)
            mol.UpdatePropertyCache(strict=False)
        for a in mol.GetAtoms():
            if a.GetAtomMapNum() > 0:
                donors.append((li, a.GetIdx()))
        ligs.append(mol)
    return graph_from_ligands(metal_z, ligs, donors), unsan, specs


def _control_string(oin: str, mode: str) -> "str | None":
    if mode != "delete":
        return oin
    parts = oin.split(".")
    # a hydride ``[H]{n}`` is not a heavy atom: deleting it is invisible to a heavy-atom graph
    # BY DESIGN (79/4,985 read ISO on the first run), so the control deletes a heavy ligand
    slotted = [
        i for i in range(1, len(parts)) if "{" in parts[i] and not parts[i].startswith("[H]")
    ]
    if not slotted:
        return None
    del parts[slotted[-1]]
    return ".".join(parts)


def _reattach(gp: ng.NeutralGraph) -> bool:
    """Move ONE metal edge from its donor to a non-donor heavy neighbour of that donor.

    The neighbour must NOT be symmetry-equivalent to the donor in the metal-free ligand graph
    (the other donors of the ligand marked): the ruler's graph has no bond orders, so in an
    unsubstituted phenyl, N2, a tetrazolate or a 1,2,4-triazolyl (N3/N4) the move would be an
    automorphism and ISO would be the CORRECT reading, not a miss (11/4,972 on the first run).
    """
    donors = set(gp.nbrs[0])
    free = Chem.RWMol()
    for a in gp.hz:
        at = Chem.Atom(int(a))
        at.SetNoImplicit(True)
        free.AddAtom(at)
    for a, b in sorted(gp.edges):
        if a != 0 and b != 0:
            free.AddBond(a, b, Chem.BondType.SINGLE)
    for d in sorted(donors):
        for x in donors - {d}:
            free.GetAtomWithIdx(x).SetIsotope(1)
        m = free.GetMol()
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        cls = list(Chem.CanonicalRankAtoms(m, breakTies=False, includeIsotopes=True))
        for x in donors - {d}:
            free.GetAtomWithIdx(x).SetIsotope(0)
        for w in gp.nbrs[d]:
            if w != 0 and w not in donors and cls[w] != cls[d]:
                gp.edges.discard((0, d))
                gp.edges.add((0, w))
                gp.nbrs[0].remove(d)
                gp.nbrs[d].remove(0)
                gp.nbrs[0].append(w)
                gp.nbrs[w].append(0)
                return True
    return False


def parseback_one(item) -> dict:
    """One molecule: parser-level and adapter-level parse-back of its string vs its XYZ."""
    rec = {"molecule": item["molecule"]}
    oin = item["oin"]
    if not oin:
        rec["parser_graph"] = rec["adapter_graph"] = "NO_STRING"
        return rec
    s_in, x_in = ng.parse_xyz(item["xyz"])
    if not s_in:
        rec["parser_graph"] = rec["adapter_graph"] = "UNREADABLE"
        return rec
    control = item.get("control")
    oin_c = _control_string(oin, control) if control else oin
    if oin_c is None:
        rec["parser_graph"] = rec["adapter_graph"] = "CONTROL_NA"
        return rec
    try:
        metal_z, ligs, donors, info = _parser_level(oin_c)
        rec.update(info)
        gp = graph_from_ligands(metal_z, ligs, donors)
        if control == "reattach" and not _reattach(gp):
            rec["parser_graph"] = rec["adapter_graph"] = "CONTROL_NA"
            return rec
        p = compare_to_input(s_in, x_in, gp)
        rec.update(
            {
                f"parser_{k}" if k in ("graph", "h_decoration", "h_diff_atoms") else k: v
                for k, v in p.items()
            }
        )
    except Exception as exc:
        rec["parser_graph"] = "PARSE_FAIL"
        rec["parser_error"] = f"{type(exc).__name__}: {exc}"[:200]
        rec["adapter_graph"] = "NA"
        return rec
    if control:
        rec["adapter_graph"] = "NA"
        return rec
    try:
        ga, unsan, _ = _adapter_level(oin, metal_z)
        a = compare_to_input(s_in, x_in, ga)
        rec["adapter_graph"] = a["graph"]
        rec["adapter_h_decoration"] = a.get("h_decoration")
        rec["adapter_h_diff_atoms"] = a.get("h_diff_atoms")
        rec["adapter_h_str"] = a["h_str"]
        rec["adapter_cn_str"] = a["cn_str"]
        rec["adapter_unsanitizable"] = unsan
    except Exception as exc:
        rec["adapter_graph"] = "ADAPTER_FAIL"
        rec["adapter_error"] = f"{type(exc).__name__}: {exc}"[:200]
    return rec


# ----------------------------------------------------------------------------------------------
# worker pool (forked from the MAIN thread, polled; SIGALRM at SIG_DFL kills a stuck C++ match)
# ----------------------------------------------------------------------------------------------


def _worker(items, part, fn, init, timeout, timeout_rec):
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    if init is not None:
        init()
    with open(part, "a") as fh:
        for item in items:
            fh.write(f"#START {item['molecule']}\n")
            fh.flush()
            signal.alarm(timeout)
            try:
                rec = fn(item)
            except Exception as exc:  # a crash is a verdict about the instrument; record it
                rec = timeout_rec(item["molecule"], f"{type(exc).__name__}: {exc}"[:200])
            signal.alarm(0)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
        os.fsync(fh.fileno())
    os._exit(0)  # skip interpreter finalizers (BLAS thread pools) on the way out


def _reconcile(part, todo, timeout_rec):
    done, started = set(), None
    if part.exists():
        for ln in part.read_text().splitlines():
            if ln.startswith("#START "):
                started = ln[7:]
            else:
                done.add(json.loads(ln)["molecule"])
    if started is not None and started not in done:
        with open(part, "a") as fh:
            fh.write(json.dumps(timeout_rec(started, "TIMEOUT")) + "\n")
        done.add(started)
    return [it for it in todo if it["molecule"] not in done]


def run_pool(items, fn, out, cpu, timeout, timeout_rec, init=None):
    """Run ``fn`` over ``items`` in ``cpu`` forked workers; JSONL to ``out`` + ``#DONE``."""
    out.parent.mkdir(parents=True, exist_ok=True)
    parts = [out.with_name(f".{out.stem}.part{k}") for k in range(cpu)]
    for p in parts:
        p.unlink(missing_ok=True)
    todo = {k: items[k::cpu] for k in range(cpu)}
    procs: dict = {}
    t0 = time.time()
    while any(todo.values()) or procs:
        for k in range(cpu):
            p = procs.get(k)
            if p is not None and p.is_alive():
                continue
            if p is not None:
                p.join()
                del procs[k]
                todo[k] = _reconcile(parts[k], todo[k], timeout_rec)
            if todo[k]:
                procs[k] = mp.Process(
                    target=_worker, args=(todo[k], parts[k], fn, init, timeout, timeout_rec)
                )
                procs[k].start()
        time.sleep(0.2)
    recs = []
    for p in parts:
        if p.exists():
            recs += [json.loads(ln) for ln in p.read_text().splitlines() if not ln.startswith("#")]
            p.unlink()
    recs.sort(key=lambda r: r["molecule"])
    with open(out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
        fh.write(f"#DONE {len(recs)}\n")
    print(f"wrote {out}  n={len(recs)}  {time.time() - t0:.0f}s")
    return recs


def _read_jsonl(path):
    return [json.loads(ln) for ln in path.read_text().splitlines() if not ln.startswith("#")]


def _print_oinsmiles():
    import oinsmiles

    print(f"oinsmiles: {oinsmiles.__file__}")


# ----------------------------------------------------------------------------------------------
# (a) parse-back
# ----------------------------------------------------------------------------------------------


def cmd_parseback(args):
    """Parse ``smiles_1`` (or ``smiles_2``) back to a graph and compare to the XYZ."""
    _print_oinsmiles()
    rows = _load_rows(args.sweep)
    mols = _check_cohort(args.cohort)
    items = []
    for m in mols:
        r = rows.get(m)
        if r is None:
            continue
        if args.side == "gen":
            xyz = args.sweep / "structures" / f"{m}_generated.xyz"
            if not xyz.exists():
                continue
            oin = r.get("smiles_2")
        else:
            xyz, oin = args.cohort / f"{m}.xyz", r.get("smiles_1")
        items.append({"molecule": m, "oin": oin, "xyz": str(xyz), "control": args.control})
    if args.n:
        items = list(np.random.default_rng(args.seed).choice(items, args.n, replace=False))
    name = "parseback" + ("_gen" if args.side == "gen" else "")
    if args.control:
        name += f"_control-{args.control}"
    print(f"cohort={len(mols)}  items={len(items)}  side={args.side}  control={args.control}")

    def trec(mol, err):
        return {
            "molecule": mol,
            "parser_graph": "TIMEOUT" if err == "TIMEOUT" else "ERROR",
            "adapter_graph": "NA",
            "parser_error": err,
        }

    recs = run_pool(items, parseback_one, args.out / f"{name}.jsonl", args.cpu, args.timeout, trec)
    summarize_parseback(recs, args.control)


def summarize_parseback(recs, control):
    """Print denominators and verdict counts; judge a control."""
    n = len(recs)
    pg = Counter(r.get("parser_graph") for r in recs)
    ag = Counter(r.get("adapter_graph") for r in recs)
    print(f"\nDENOMINATOR n={n}")
    print("parser_graph :", dict(pg.most_common()))
    print("adapter_graph:", dict(ag.most_common()))
    print("parser_h     :", dict(Counter(r.get("parser_h_decoration") for r in recs).most_common()))
    print(
        "adapter_h    :", dict(Counter(r.get("adapter_h_decoration") for r in recs).most_common())
    )
    if not control:
        return
    usable = [
        r for r in recs if r.get("parser_graph") not in ("CONTROL_NA", "NO_STRING", "UNREADABLE")
    ]
    bad = [r for r in usable if r.get("parser_graph") in ISO_SET]
    print(
        f"CONTROL {control}: usable={len(usable)}  FAIL(ISO)={len(bad)}  ->  "
        f"{'PASS' if usable and not bad else 'FAIL'}"
    )
    for r in bad[:12]:
        print("   ", r["molecule"], r.get("parser_graph"), r.get("cn_in"), r.get("cn_str"))


# ----------------------------------------------------------------------------------------------
# (b) collision scan
# ----------------------------------------------------------------------------------------------


def cmd_collide(args):
    """Cluster the inputs by heavy-graph isomorphism; judge strings within each cluster."""
    rows = _load_rows(args.sweep)
    mols = _check_cohort(args.cohort)
    t0 = time.time()
    graphs, by_formula = {}, defaultdict(list)
    for m in mols:
        s, x = ng.parse_xyz(args.cohort / f"{m}.xyz")
        if not s:
            continue
        g = ng.NeutralGraph(s, x)
        g._mol = g.rdmol()
        graphs[m] = g
        by_formula[g.formula(with_h=True)].append(m)
    print(
        f"inputs read: {len(graphs)}/{len(mols)}  formula groups: {len(by_formula)}  "
        f"({time.time() - t0:.0f}s)"
    )
    # union-find over isomorphic pairs inside each formula group
    parent = {m: m for m in graphs}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    n_pairs_tested = 0
    for f, members in by_formula.items():
        if len(members) < 2:
            continue
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                if find(a) == find(b):
                    continue
                n_pairs_tested += 1
                if ng.isomorphism(graphs[a], graphs[b], m_a=graphs[a]._mol, m_b=graphs[b]._mol):
                    parent[find(a)] = find(b)
    clusters = defaultdict(list)
    for m in graphs:
        clusters[find(m)].append(m)
    multi = {k: sorted(v) for k, v in clusters.items() if len(v) > 1}
    print(
        f"isomorphism pairs tested: {n_pairs_tested}  clusters(size>1): {len(multi)}  "
        f"molecules in them: {sum(len(v) for v in multi.values())}  ({time.time() - t0:.0f}s)"
    )

    def refcode(m):
        return m.split("_comp_")[0]

    pairs, cluster_out = [], []
    for k, members in multi.items():
        strings = Counter(rows[m].get("smiles_1") for m in members)
        first = members[0]
        rel = {}
        # star: every member vs the first (stereo relation by the ruler, input vs input)
        for m in members[1:]:
            v = ng.compare(args.cohort / f"{first}.xyz", args.cohort / f"{m}.xyz")
            rel[(first, m)] = v
        # plus every pair whose strings DIFFER (the non-canonical question needs the pair itself)
        for i in range(1, len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                if rows[a].get("smiles_1") != rows[b].get("smiles_1"):
                    rel[(a, b)] = ng.compare(args.cohort / f"{a}.xyz", args.cohort / f"{b}.xyz")
        for (a, b), v in rel.items():
            sa, sb = rows[a].get("smiles_1"), rows[b].get("smiles_1")
            st = v.get("stereo", "NA")
            if sa is None or sb is None:
                cls = "NO_STRING"
            elif sa == sb:
                cls = (
                    "CONSISTENT"
                    if st in PASS_STEREO
                    else "COLLISION_ENANTIOMER"
                    if st.startswith("MIRROR")
                    else "COLLISION_DIASTEREOMER"
                    if st == "DIFFERENT"
                    else "UNDECIDED_SAME_STRING"
                )
            else:
                cls = (
                    "NON_CANONICAL"
                    if st in PASS_STEREO
                    else "DISTINGUISHED"
                    if st.startswith("MIRROR") or st == "DIFFERENT"
                    else "UNDECIDED_DIFF_STRING"
                )
            pairs.append(
                {
                    "a": a,
                    "b": b,
                    "cluster": k,
                    "class": cls,
                    "stereo": st,
                    "graph": v.get("graph"),
                    "sphere": v.get("sphere"),
                    "same_string": sa == sb,
                    "sibling": refcode(a) == refcode(b),
                    "bucket_a": rows[a]["bucket"],
                    "bucket_b": rows[b]["bucket"],
                    "h_decoration": v.get("h_decoration"),
                }
            )
        cluster_out.append(
            {
                "cluster": k,
                "size": len(members),
                "members": members,
                "n_strings": len([s for s in strings if s]),
                "formula": graphs[first].formula(True),
                "n_heavy": len(graphs[first].heavy),
            }
        )
    # same string across DIFFERENT clusters: the string cannot tell two graphs apart
    by_string = defaultdict(set)
    for m in graphs:
        s = rows[m].get("smiles_1")
        if s:
            by_string[s].add(find(m))
    graph_collisions = []
    for s, cl in by_string.items():
        if len(cl) > 1:
            mem = [m for m in graphs if rows[m].get("smiles_1") == s]
            graph_collisions.append(
                {
                    "string": s,
                    "n_clusters": len(cl),
                    "members": sorted(mem),
                    "formulas": sorted({graphs[m].formula(True) for m in mem}),
                    "n_heavy": sorted({len(graphs[m].heavy) for m in mem}),
                }
            )
    summary = {
        "n_inputs": len(graphs),
        "n_formula_groups": len(by_formula),
        "n_iso_pairs_tested": n_pairs_tested,
        "n_clusters_multi": len(multi),
        "n_molecules_in_multi": sum(len(v) for v in multi.values()),
        "n_sibling_refcodes": sum(1 for v in Counter(refcode(m) for m in graphs).values() if v > 1),
        "n_pairs_judged": len(pairs),
        "pair_classes": dict(Counter(p["class"] for p in pairs).most_common()),
        "pair_classes_siblings": dict(
            Counter(p["class"] for p in pairs if p["sibling"]).most_common()
        ),
        "clusters_with_multiple_strings": sum(1 for c in cluster_out if c["n_strings"] > 1),
        "n_strings_shared_across_clusters": len(graph_collisions),
        "n_molecules_in_graph_collisions": sum(len(c["members"]) for c in graph_collisions),
    }
    out = args.out / "collisions.json"
    out.write_text(
        json.dumps(
            {
                "summary": summary,
                "clusters": cluster_out,
                "pairs": pairs,
                "graph_collisions": graph_collisions,
            },
            indent=1,
        )
    )
    print(json.dumps(summary, indent=1))
    print(f"wrote {out}  ({time.time() - t0:.0f}s)")


# ----------------------------------------------------------------------------------------------
# (c) perception flags + charge probe
# ----------------------------------------------------------------------------------------------

_PROBE = {"charge": None, "mol": None}


def _pflags_init():
    import oinsmiles.utils.perception_tmc as pt

    orig = pt.get_tmc_mol

    def patched(xyz_file, overall_charge, with_stereo=False):
        q = _PROBE["charge"] if _PROBE["charge"] is not None else overall_charge
        tmc_mol, xyz = orig(xyz_file, q, with_stereo=with_stereo)
        _PROBE["mol"] = tmc_mol
        return tmc_mol, xyz

    pt.get_tmc_mol = patched


def _mol_flags(mol) -> dict:
    from oinsmiles.utils.perception_tmc import TRANSITION_METALS_NUM

    metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    ligs = [a for a in mol.GetAtoms() if a.GetAtomicNum() not in TRANSITION_METALS_NUM]
    return {
        "metal": [a.GetSymbol() for a in metals],
        "os_perceived": [a.GetFormalCharge() for a in metals],
        "lig_charge_sum": sum(a.GetFormalCharge() for a in ligs),
        "total_charge": Chem.GetFormalCharge(mol),
        "n_charged_atoms": sum(1 for a in ligs if a.GetFormalCharge() != 0),
        "n_abs_charge_ge2": sum(1 for a in ligs if abs(a.GetFormalCharge()) >= 2),
        "n_charged_carbon": sum(1 for a in ligs if a.GetAtomicNum() == 6 and a.GetFormalCharge()),
        "n_radical_atoms": sum(1 for a in ligs if a.GetNumRadicalElectrons() > 0),
        "n_atoms": mol.GetNumAtoms(),
    }


def pflags_one(item) -> dict:
    """Perceive one input the way the encoder does (charge 0); optionally probe the true charge."""
    import oinsmiles.utils.perception_tmc as pt
    from oinsmiles.utils.perception_tmc import ALLOWED_OXIDATION_STATES

    rec = {"molecule": item["molecule"], "charge_stated": _comment_charge(item["xyz"])}
    q = rec["charge_stated"]
    _PROBE["charge"] = None
    try:
        with _silence_fds():
            mol, _ = pt.get_tmc_mol(Path(item["xyz"]), 0)
        rec.update(_mol_flags(mol))
        rec["p_status"] = "OK"
    except Exception as exc:
        rec["p_status"] = f"{type(exc).__name__}: {exc}"[:200]
        return rec
    allowed = [ALLOWED_OXIDATION_STATES.get(s) for s in rec["metal"]]
    rec["os_in_allowed"] = all(
        a is not None and o in a for a, o in zip(allowed, rec["os_perceived"])
    )
    if q is not None:
        # perception put the whole stated charge on the metal (tm_ox = 0 - lig charges), so
        # honouring it would move the metal by exactly q -- if the ligand assignment holds
        rec["os_if_honoured"] = [o + q for o in rec["os_perceived"]]
        rec["os_if_honoured_in_allowed"] = all(
            a is not None and o in a for a, o in zip(allowed, rec["os_if_honoured"])
        )
    if not item.get("probe") or q in (None, 0):
        return rec
    from oinsmiles import XYZToSMILES

    try:
        _PROBE["charge"] = 0
        with _silence_fds():
            oin0 = XYZToSMILES().convert(item["xyz"])
        rec["probe_control_equal"] = oin0 == item["oin"]
        _PROBE["charge"] = q
        with _silence_fds():
            oinq = XYZToSMILES().convert(item["xyz"])
        fq = _mol_flags(_PROBE["mol"])
        rec["probe_os_perceived"] = fq["os_perceived"]
        rec["probe_lig_charge_sum"] = fq["lig_charge_sum"]
        rec["probe_total_charge"] = fq["total_charge"]
        rec["probe_fired"] = fq["total_charge"] == q
        rec["probe_equal"] = oinq == item["oin"]
        if not rec["probe_equal"]:
            rec["probe_oin"] = oinq
    except Exception as exc:
        rec["probe_error"] = f"{type(exc).__name__}: {exc}"[:200]
    finally:
        _PROBE["charge"] = None
    return rec


def cmd_pflags(args):
    """Perception plausibility flags for every input; ``--charge-probe`` re-encodes charged ones."""
    _print_oinsmiles()
    rows = _load_rows(args.sweep)
    mols = _check_cohort(args.cohort)
    items = [
        {
            "molecule": m,
            "xyz": str(args.cohort / f"{m}.xyz"),
            "oin": rows.get(m, {}).get("smiles_1"),
            "probe": args.charge_probe,
        }
        for m in mols
    ]
    if args.n:
        items = list(np.random.default_rng(args.seed).choice(items, args.n, replace=False))
    stated = Counter(_comment_charge(it["xyz"]) for it in items)
    print(f"cohort={len(mols)}  items={len(items)}  stated charge: {dict(stated.most_common())}")

    def trec(mol, err):
        return {"molecule": mol, "p_status": err}

    recs = run_pool(
        items,
        pflags_one,
        args.out / "pflags.jsonl",
        args.cpu,
        args.timeout,
        trec,
        init=_pflags_init,
    )
    summarize_pflags(recs)


def summarize_pflags(recs):
    """Denominators for the perception flags and the charge probe."""
    n = len(recs)
    ok = [r for r in recs if r.get("p_status") == "OK"]
    print(
        f"\nDENOMINATOR n={n}  perceived OK={len(ok)}  "
        f"status: {dict(Counter(r.get('p_status', '')[:30] for r in recs).most_common(6))}"
    )
    print("stated charge   :", dict(Counter(r.get("charge_stated") for r in recs).most_common()))
    print("os_in_allowed   :", dict(Counter(r.get("os_in_allowed") for r in ok).most_common()))
    charged = [r for r in ok if r.get("charge_stated") not in (None, 0)]
    print(
        f"stated-charged  : {len(charged)}  os_if_honoured_in_allowed: "
        f"{dict(Counter(r.get('os_if_honoured_in_allowed') for r in charged).most_common())}"
    )
    probed = [r for r in charged if "probe_equal" in r or "probe_error" in r]
    if probed:
        fired = sum(1 for r in probed if r.get("probe_fired"))
        ctrl = sum(1 for r in probed if r.get("probe_control_equal"))
        print(
            f"charge probe    : probed={len(probed)}  control(charge0==smiles_1)={ctrl}  "
            f"fired(total==stated)={fired}  string_changed="
            f"{sum(1 for r in probed if r.get('probe_equal') is False)}  "
            f"errors={sum(1 for r in probed if 'probe_error' in r)}"
        )


# ----------------------------------------------------------------------------------------------
# summary: join with the bucket, the C1 ruler verdict and the C2 mirror column
# ----------------------------------------------------------------------------------------------


def _xtab(recs, key, rows, label):
    tab = defaultdict(Counter)
    for r in recs:
        b = rows.get(r["molecule"], {}).get("bucket", "?")
        tab[b][str(r.get(key))] += 1
    print(f"\n{label} x bucket")
    for b in sorted(tab, key=lambda x: -sum(tab[x].values())):
        tot = sum(tab[b].values())
        print(f"  {b:17s} n={tot:5d}  ", dict(tab[b].most_common()))
    return {b: dict(c) for b, c in tab.items()}


def cmd_summary(args):
    """Cross-tabulate every I3 output against the honest bucket; write the summary JSON."""
    rows = _load_rows(args.sweep)
    out = {"n_rows": len(rows)}
    pb = args.out / "parseback.jsonl"
    if pb.exists():
        recs = _read_jsonl(pb)
        out["parseback_n"] = len(recs)
        out["parser_graph_x_bucket"] = _xtab(recs, "parser_graph", rows, "parser_graph")
        out["adapter_graph_x_bucket"] = _xtab(recs, "adapter_graph", rows, "adapter_graph")
        out["parser_h_x_bucket"] = _xtab(recs, "parser_h_decoration", rows, "parser H decoration")
        out["adapter_h_x_bucket"] = _xtab(
            recs, "adapter_h_decoration", rows, "adapter H decoration"
        )
        out["unslotted_x_bucket"] = _xtab(
            [dict(r, unslotted=bool(r.get("n_unslotted"))) for r in recs],
            "unslotted",
            rows,
            "has an unslotted fragment",
        )
        out["str_valence_x_bucket"] = _xtab(
            [dict(r, vf=bool(r.get("str_valence_fail"))) for r in recs],
            "vf",
            rows,
            "ligand as written fails RDKit sanitize",
        )
        errs = Counter(
            (r.get("adapter_error") or "")[:70]
            for r in recs
            if r.get("adapter_graph") == "ADAPTER_FAIL"
        )
        out["adapter_errors"] = dict(errs.most_common(12))
        print("\nadapter errors:", json.dumps(out["adapter_errors"], indent=1))
        # cross-check with C1: byte_exact passes whose GENERATED graph was not ISO to the input
        gv = args.out / "g_verdict.jsonl"
        if gv.exists():
            g1 = {r["molecule"]: r for r in _read_jsonl(gv)}
            xc = Counter()
            for r in recs:
                v = g1.get(r["molecule"])
                if v is None or rows[r["molecule"]]["bucket"] != "byte_exact":
                    continue
                if v.get("graph") in ("ISO", "ISO_MARGINAL", "ISO_CLASH"):
                    continue
                xc[(v.get("graph"), r.get("parser_graph"), r.get("adapter_graph"))] += 1
            out["byte_exact_C1_nonISO_x_parseback"] = {
                " | ".join(map(str, k)): v for k, v in xc.items()
            }
            print("\nbyte_exact passes with a C1 non-ISO generated graph, by parse-back verdict:")
            for k, v in xc.most_common():
                print(f"  C1={k[0]:13s} parser={k[1]:13s} adapter={k[2]:13s} {v}")
    pg = args.out / "parseback_gen.jsonl"
    if pg.exists():
        recs = _read_jsonl(pg)
        out["parseback_gen_n"] = len(recs)
        out["gen_parser_graph_x_bucket"] = _xtab(
            recs, "parser_graph", rows, "GEN side parser_graph"
        )
        out["gen_adapter_h_x_bucket"] = _xtab(
            recs, "adapter_h_decoration", rows, "GEN side adapter H"
        )
    pf = args.out / "pflags.jsonl"
    if pf.exists():
        recs = _read_jsonl(pf)
        out["pflags_n"] = len(recs)
        out["charge_stated_x_bucket"] = _xtab(recs, "charge_stated", rows, "stated charge")
        out["os_in_allowed_x_bucket"] = _xtab(
            recs, "os_in_allowed", rows, "perceived OS in allowed set"
        )
        charged = [r for r in recs if r.get("charge_stated") not in (None, 0)]
        out["os_if_honoured_x_bucket"] = _xtab(
            charged,
            "os_if_honoured_in_allowed",
            rows,
            "stated-charged: OS if honoured in allowed set",
        )
        probed = [r for r in charged if "probe_equal" in r]
        if probed:
            out["probe_equal_x_bucket"] = _xtab(
                probed, "probe_equal", rows, "charge probe: string unchanged"
            )
            out["probe_fired"] = sum(1 for r in probed if r.get("probe_fired"))
            out["probe_control_equal"] = sum(1 for r in probed if r.get("probe_control_equal"))
            out["probe_n"] = len(probed)
        out["radical_x_bucket"] = _xtab(
            [dict(r, rad=bool(r.get("n_radical_atoms"))) for r in recs],
            "rad",
            rows,
            "ligand radical atoms present",
        )
        out["charged_carbon_x_bucket"] = _xtab(
            [dict(r, cc=bool(r.get("n_charged_carbon"))) for r in recs],
            "cc",
            rows,
            "charged carbon present",
        )
        # C2's mirror column beside the charge flag: is perception fragility charge-correlated?
        es = args.out / "e_selfconsistency.jsonl"
        if es.exists():
            e2 = {r["molecule"]: r for r in _read_jsonl(es)}
            frag = Counter()
            for r in recs:
                e = e2.get(r["molecule"])
                if not e or not e.get("base"):
                    continue
                rel = e.get("rel", {})
                fragile = any(
                    rel.get(k) not in (None, "byte_exact")
                    for k in ("renum0", "renum1", "renum2", "noise0", "noise1", "noise2")
                )
                frag[(r.get("charge_stated") not in (None, 0), fragile)] += 1
            out["stated_charged_x_C2_fragile"] = {
                f"charged={k[0]} fragile={k[1]}": v for k, v in frag.items()
            }
            print("\nstated-charged x C2 perception-fragile:", out["stated_charged_x_C2_fragile"])
    co = args.out / "collisions.json"
    if co.exists():
        out["collisions"] = json.loads(co.read_text())["summary"]
    (args.out / "string_sufficiency_summary.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {args.out / 'string_sufficiency_summary.json'}")


def main():
    """CLI."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["parseback", "collide", "pflags", "summary"])
    ap.add_argument("--cohort", type=Path, default=DEF_COHORT)
    ap.add_argument("--sweep", type=Path, default=DEF_SWEEP)
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument("--control", choices=["delete", "reattach"], default=None)
    ap.add_argument("--side", choices=["input", "gen"], default="input")
    ap.add_argument("--charge-probe", action="store_true")
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=120, help="seconds per molecule (SIGKILL)")
    ap.add_argument("--n", type=int, default=0, help="random subsample (0 = all)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    for p in (args.cohort, args.sweep):
        if not p.exists():
            sys.exit(f"missing: {p}")
    {
        "parseback": cmd_parseback,
        "collide": cmd_collide,
        "pflags": cmd_pflags,
        "summary": cmd_summary,
    }[args.cmd](args)


if __name__ == "__main__":
    main()
