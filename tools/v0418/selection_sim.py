"""Selection lane, step 0: can an INPUT-FREE predicate choose the right eta structure? Reads only.

The L2 A/B left three generated structures per eta molecule -- arms ``off`` (upstream target),
``on`` (the promoted pair) and ``on3`` (+ OIN_ETA_TARGET_UNSCALED) -- each with a known verdict.
A per-molecule oracle over them is +222 VERIFIED with no losses, against +150 for the best single
configuration. The generator never sees the input xyz, so a real selector can only compare a
candidate against THE STRING. This simulates exactly that:

    declared   per element, how many atoms smiles_1 tags as binding ({n}, {n>}, {n<})
    measured   per element, how many atoms of the generated structure sit inside the encoder's
               contact cutoff of the metal (covalent radii + 0.45 A) -- oin.coordination
    deficit    sum over elements of max(0, declared - measured)      a donor that did not arrive
    surplus    sum over elements of max(0, measured - declared)      an atom that should not bind

and asks: if each molecule took the arm the predicate prefers (ties -> the shipped default), how
many VERIFIED molecules would that be -- and how many would it LOSE? Unlike an offline re-score
this CAN express a loss: every arm's verdict is known per molecule.

WHAT IT CANNOT SAY. One pool holding both policies is not three separate runs: the fill order, the
dedup and the early exit differ. This BOUNDS the lane and ranks predicates; it is not the A/B.

WHAT A BROKEN VERSION WOULD PRINT. A tokenizer that misses the binding tags declares 0 everywhere,
every structure scores deficit 0 and the predicate never leaves the default -- so the script prints
how many molecules the predicate DISCRIMINATES on, and the declared total must equal the input's own
contact count for most verified passes (printed as the control line).

    <main>/.venv/bin/python tools/v0418/selection_sim.py [--out selection_sim.json]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "src"))  # beats PYTHONPATH, which is the point

from oinsmiles.oin.coordination import metal_contacts, metal_indices, parse_xyz  # noqa: E402

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
AB = DATA / "results-v0.4.18-eta-detached"
COHORT = DATA / "cohort-v0.4.18-eta"
TABLE = DATA / "results-v0.4.17-reattribution" / "attribution_table.tsv"
ARMS = ("on", "on3", "off")  # tie-break order: the shipped default first
ISO = ("ISO", "ISO_MARGINAL", "ISO_CLASH")
STRING_FAULTS = ("P_E1_COVERAGE", "DATA_MULTI", "P_DETACHED", "E1_GRAPH", "E1_HCOUNT")
BAD_STEREO = ("MIRROR", "MIRROR_PARTIAL", "DIFFERENT")

#: an atom token immediately followed by a binding tag: [cH]{0>}  N{1}  c{2}  [CH3]{0}  Cl{3}
_TAGGED = re.compile(r"(\[[^\]]+\]|Cl|Br|[BCNOPSFI]|[bcnops])\{\d+[<>]?\}")
_ELEM = re.compile(r"^\[?\d*([A-Z][a-z]?|[a-z]{1,2})")


def declared(oin: str) -> Counter:
    out = Counter()
    for tok in _TAGGED.findall(oin):
        m = _ELEM.match(tok)
        sym = m.group(1) if m else ""
        sym = sym.capitalize() if sym.islower() else sym  # aromatic c/n/se -> C/N/Se
        if sym in ("Ch", "Nh", "Oh", "Bh", "Ph", "Sh"):  # [cH] read as two lower-case letters
            sym = sym[0]
        out[sym] += 1
    return out


def measured(xyz_text: str):
    sym, xyz = parse_xyz(xyz_text)
    metals = metal_indices(sym) if sym else []
    if len(metals) != 1:
        return None
    counts, _ = metal_contacts(sym, xyz, metals[0])
    return Counter(counts)


def score(dec: Counter, mea: Counter):
    els = set(dec) | set(mea)
    return (
        sum(max(0, dec[e] - mea[e]) for e in els),
        sum(max(0, mea[e] - dec[e]) for e in els),
    )


def _jsonl(path):
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not rows[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    return {r["molecule"]: r for r in map(json.loads, rows[:-1])}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    T = {r["molecule"]: r for r in csv.DictReader(open(TABLE), delimiter="\t")}
    B, G = {}, {}
    for a in ARMS:
        B[a] = {
            r["molecule"].removesuffix(".xyz"): r
            for r in json.loads((AB / f"ab_{a}" / "bucket_report_honest.json").read_text())
        }
        G[a] = _jsonl(AB / f"ab_{a}" / "g_verdict.jsonl")
    mols = sorted(p.stem for p in COHORT.glob("*.xyz"))
    if len(mols) != 1146:
        sys.exit(f"ABORT: cohort has {len(mols)} molecules")

    def verified(a, m):
        if B[a][m]["bucket"] != "byte_exact" or T[m]["fault"] in STRING_FAULTS:
            return False
        g = G[a].get(m)
        return bool(g) and g.get("graph") in ISO and g.get("stereo") not in BAD_STEREO

    ver = {a: {m: verified(a, m) for m in mols} for a in ARMS}
    print(f"DENOMINATOR {len(mols)} eta-bound molecules, three arms each")
    print("  verified per arm:", {a: sum(ver[a].values()) for a in ARMS})
    oracle = sum(any(ver[a][m] for a in ARMS) for m in mols)
    print(f"  oracle (best arm per molecule): {oracle}")

    # ---- the predicate on every structure ---------------------------------------------------
    S, unscored, ctrl_ok, ctrl_n = {}, 0, 0, 0
    for m in mols:
        oin = B["on"][m].get("smiles_1") or ""
        dec = declared(oin)
        if not dec:
            unscored += 1
            continue
        mi = measured((COHORT / f"{m}.xyz").read_text())
        if mi is None:
            unscored += 1
            continue
        if ver["off"][m]:  # control: on a verified pass the string must describe the INPUT
            ctrl_n += 1
            ctrl_ok += score(dec, mi) == (0, 0)
        S[m] = {}
        for a in ARMS:
            p = AB / f"ab_{a}" / "structures" / f"{m}_generated.xyz"
            mea = measured(p.read_text()) if p.exists() else None
            S[m][a] = score(dec, mea) if mea is not None else None
    print(
        f"  scored {len(S)} / {len(mols)}  (unscored {unscored}: no tags parsed or not single-metal)"
    )
    print(
        f"  CONTROL declared == the INPUT's own contacts on verified passes: {ctrl_ok}/{ctrl_n} "
        f"({100 * ctrl_ok / max(ctrl_n, 1):.1f}%)  -- the rest are BOUNDARY contacts / tag parsing"
    )

    # ---- how well does (deficit, surplus) == (0, 0) track the verdict? -------------------------
    print("\n1. DOES THE PREDICATE TRACK THE VERDICT?   per structure, all arms pooled")
    tab = Counter()
    for m, by in S.items():
        for a, sc in by.items():
            if sc is not None:
                tab[(sc == (0, 0), ver[a][m])] += 1
    for clean in (True, False):
        v, nv = tab[(clean, True)], tab[(clean, False)]
        print(
            f"   {'clean (0,0)  ' if clean else 'not clean    '} verified {v:5d}   not verified {nv:5d}"
            f"   P(verified) = {100 * v / max(v + nv, 1):.1f}%"
        )

    # ---- simulated selectors ---------------------------------------------------------------
    def pick(m, key):
        by = S.get(m)
        if not by:
            return "on"
        cands = [a for a in ARMS if by[a] is not None]
        return min(cands, key=lambda a: (key(by[a]), ARMS.index(a))) if cands else "on"

    print("\n2. SIMULATED SELECTORS   (ties -> on, the shipped default)     base: on =", end=" ")
    base = sum(ver["on"].values())
    print(base)
    rep = {"base_on": base, "oracle": oracle}
    for name, key, arms in (
        ("deficit, then surplus   over on/on3/off", lambda s: s, ARMS),
        ("deficit only            over on/on3/off", lambda s: (s[0],), ARMS),
        ("deficit, then surplus   over on/on3    ", lambda s: s, ("on", "on3")),
        ("deficit only            over on/on3    ", lambda s: (s[0],), ("on", "on3")),
    ):

        def choose(m, key=key, arms=arms):
            by = S.get(m)
            if not by:
                return "on"
            cands = [a for a in arms if by[a] is not None]
            return min(cands, key=lambda a: (key(by[a]), arms.index(a))) if cands else "on"

        ch = {m: choose(m) for m in mols}
        got = sum(ver[ch[m]][m] for m in mols)
        gains = sum(ver[ch[m]][m] and not ver["on"][m] for m in mols)
        losses = sum(ver["on"][m] and not ver[ch[m]][m] for m in mols)
        moved = Counter(ch[m] for m in mols if ch[m] != "on")
        print(
            f"   {name}: {got}  ({got - base:+d} vs on: +{gains} / -{losses})   left the default on "
            f"{sum(moved.values())} molecules {dict(moved)}"
        )
        rep[name.strip()] = {"verified": got, "gains": gains, "losses": losses}

    # ---- where is the oracle's headroom, in predicate terms? --------------------------------
    print("\n3. THE ORACLE'S HEADROOM   molecules NOT verified in `on` but verified in another arm")
    head = [m for m in mols if not ver["on"][m] and any(ver[a][m] for a in ("on3", "off"))]
    seen = Counter()
    for m in head:
        by = S.get(m)
        if not by or by["on"] is None:
            seen["unscored / no `on` structure"] += 1
            continue
        better = [a for a in ("on3", "off") if ver[a][m] and by[a] is not None]
        if any(by[a] < by["on"] for a in better):
            seen["predicate PREFERS the verified arm"] += 1
        elif any(by[a] == by["on"] for a in better):
            seen["predicate TIES (cannot tell)"] += 1
        else:
            seen["predicate prefers the WRONG arm"] += 1
    print(f"   {len(head)} molecules:", dict(seen))
    rep["headroom"] = dict(seen)
    if args.out:
        args.out.write_text(json.dumps(rep, indent=1, sort_keys=True))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
