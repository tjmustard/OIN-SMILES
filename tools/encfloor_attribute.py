#!/usr/bin/env python
"""Per-``AC2BO``-call cost attribution for the encode floor (v0.4.7, lane L3-encfloor).

``docs/agentic-notes/v0.4.5/ENCODER_PERF_v0.4.5.md`` attributes 99.8 % of a slow encode to ``AC2BO`` and stops
there. That attribution predates the default-ON ``OIN_CANONICAL_PERCEPTION`` wrapper, so
two functions that now run on **every** ``AC2BO`` call -- ``_canonical_atom_permutation``
and ``_valence_search_is_truncated`` -- have never appeared in any measurement, and neither
has ``_ordered_valences``, the sub-cap enumeration that 99.8 % of the corpus takes.

This tool opens ``AC2BO`` up. For each call it records, in one pass:

  * ``combo_size`` and whether the call was over cap;
  * inclusive wall of ``_AC2BO_core``, of the ``_ordered_valences`` materialisation inside
    it, and of the two wrapper-level helpers;
  * how many candidate valence assignments the loop actually consumed, and how many
    ``nx.max_weight_matching`` calls it made (deltas of ``AC2BO_STATS``, which are exact
    and load-independent);
  * the **redundancy** question for LEAD 2: total calls vs distinct arguments for
    ``_canonical_atom_permutation`` / ``_valence_search_is_truncated`` / ``possible_valences``,
    keyed on the arguments they actually read.

Wall clock on this host is contended, so every *ratio* below is computed within a single
process and single run; absolute seconds are ADVISORY. The call counts and candidate counts
are not -- they are exact.

    $V tools/encfloor_attribute.py --dataset $DS --molecules QIDKIZ_comp_0,XIRMER_comp_0
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    if os.path.exists(molecule):
        return molecule
    for sub in ("cat", "photo", "regression_inputs", "."):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    p = os.path.join("tests/fixtures", f"{molecule}.xyz")
    if os.path.exists(p):
        return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}")


def _ac_key(AC, atoms, *extra):
    h = hashlib.sha256()
    h.update(AC.tobytes())
    h.update(repr(list(atoms)).encode())
    h.update(repr(extra).encode())
    return h.hexdigest()[:16]


class Recorder:
    """Instruments perception_core in place; ``restore()`` puts it back."""

    def __init__(self, loc):
        self.loc = loc
        self.restores = []
        self.calls = []  # one dict per _AC2BO_core call
        self.helper = {}  # name -> {"calls": n, "wall": s, "keys": set}
        self._pending_ordered = None

    def _slot(self, name):
        return self.helper.setdefault(name, {"calls": 0, "wall": 0.0, "keys": set()})

    def _patch(self, name, fn):
        orig = getattr(self.loc, name)
        self.restores.append((name, orig))
        setattr(self.loc, name, fn(orig))

    def install(self):
        loc = self.loc

        def wrap_core(orig):
            def core(AC, atoms, charge, **kw):
                self._pending_ordered = None
                c0 = dict(loc.AC2BO_STATS)
                t0 = time.perf_counter()
                try:
                    return orig(AC, atoms, charge, **kw)
                finally:
                    dt = time.perf_counter() - t0
                    c1 = loc.AC2BO_STATS
                    ordered = self._pending_ordered
                    self.calls.append(
                        {
                            "n_atoms": len(atoms),
                            "charge": charge,
                            "allow_carbenes": kw.get("allow_carbenes", True),
                            "ac_key": _ac_key(AC, atoms),
                            "wall_core": dt,
                            "candidates": c1["candidates"] - c0["candidates"],
                            "matching": c1["matching_calls"] - c0["matching_calls"],
                            "over_cap": (c1["over_cap_calls"] - c0["over_cap_calls"]) > 0,
                            "found_valid": (c1["found_valid"] - c0["found_valid"]) > 0,
                            "ordered_wall": None if ordered is None else ordered[0],
                            "ordered_len": None if ordered is None else ordered[1],
                        }
                    )

            return core

        def wrap_ordered(orig):
            def ordered(vll, atoms):
                t0 = time.perf_counter()
                out = orig(vll, atoms)
                dt = time.perf_counter() - t0
                self._pending_ordered = (dt, len(out))
                s = self._slot("_ordered_valences")
                s["calls"] += 1
                s["wall"] += dt
                return out

            return ordered

        def wrap_keyed(name, keyfn):
            def maker(orig):
                def fn(*a, **kw):
                    t0 = time.perf_counter()
                    try:
                        return orig(*a, **kw)
                    finally:
                        s = self._slot(name)
                        s["calls"] += 1
                        s["wall"] += time.perf_counter() - t0
                        try:
                            s["keys"].add(keyfn(*a, **kw))
                        except Exception:  # noqa: BLE001 - key is diagnostic only
                            pass

                return fn

            return maker

        self._patch("_AC2BO_core", wrap_core)
        self._patch("_ordered_valences", wrap_ordered)
        self._patch(
            "_canonical_atom_permutation",
            wrap_keyed("_canonical_atom_permutation", lambda AC, atoms: _ac_key(AC, atoms)),
        )
        self._patch(
            "_valence_search_is_truncated",
            wrap_keyed(
                "_valence_search_is_truncated",
                lambda AC, atoms, allow_carbenes=True: _ac_key(AC, atoms, allow_carbenes),
            ),
        )
        self._patch(
            "possible_valences",
            wrap_keyed(
                "possible_valences",
                lambda AC_valence, atoms, allow_carbenes=True: _ac_key(
                    __import__("numpy").asarray(list(AC_valence)), atoms, allow_carbenes
                ),
            ),
        )
        self._patch(
            "valence_combo_size",
            wrap_keyed("valence_combo_size", lambda vll: repr(vll)[:200] + str(len(vll))),
        )

    def install_resonance(self, x2m):
        """Instrument the SL5 forked resonance path in ``utils.perception_tmc``.

        A sub-cap ligand fragment can still take tens of minutes to encode with ``AC2BO``
        contributing ~nothing (``HACYEQ_comp_0``: 54.89 s encode, ``_AC2BO_core`` 0.13 s =
        0.2 %). ``lig_checks`` forks a CPU-budgeted child per call
        (``_RESONANCE_CPU_BUDGET_S`` = 120 CPU seconds, 900 s wall backstop), and the
        charge ladder plus ``_rescue_unusable_perception`` can call it many times over --
        so the floor for that cohort is ``k x 120`` CPU seconds, not valence combinatorics.
        The ``status`` histogram is the load-independent evidence: ``timeout`` means the
        child burned the whole budget.
        """
        self.resonance = {"calls": 0, "wall": 0.0, "status": {}, "lig_checks": 0}
        self.res_restores = []

        orig_iso = getattr(x2m, "_resonance_candidates_isolated", None)
        if orig_iso is not None:

            def iso(*a, **kw):
                t0 = time.perf_counter()
                out = orig_iso(*a, **kw)
                self.resonance["calls"] += 1
                self.resonance["wall"] += time.perf_counter() - t0
                st = out[0] if isinstance(out, tuple) else "?"
                self.resonance["status"][st] = self.resonance["status"].get(st, 0) + 1
                return out

            x2m._resonance_candidates_isolated = iso
            self.res_restores.append((x2m, "_resonance_candidates_isolated", orig_iso))

        orig_lc = getattr(x2m, "lig_checks", None)
        if orig_lc is not None:

            def lc(*a, **kw):
                self.resonance["lig_checks"] += 1
                return orig_lc(*a, **kw)

            x2m.lig_checks = lc
            self.res_restores.append((x2m, "lig_checks", orig_lc))

    def restore(self):
        for name, orig in reversed(self.restores):
            setattr(self.loc, name, orig)
        self.restores = []
        for obj, attr, orig in reversed(getattr(self, "res_restores", [])):
            setattr(obj, attr, orig)
        self.res_restores = []


def run_one(path, name):
    from oinsmiles import XYZToSMILES
    from oinsmiles.utils import perception_core as loc
    from oinsmiles.utils import perception_tmc as x2m

    rec = Recorder(loc)
    loc.reset_ac2bo_stats()
    clear = getattr(loc, "_ac2bo_memo_clear", None)
    if clear is not None:
        clear()
    rec.install()
    rec.install_resonance(x2m)
    t0 = time.perf_counter()
    err = None
    oin = None
    try:
        oin = XYZToSMILES().convert(path)
    except Exception as exc:  # noqa: BLE001 - a failing encode is still an attribution
        err = f"{type(exc).__name__}: {exc}"
    finally:
        wall = time.perf_counter() - t0
        rec.restore()

    helper = {
        k: {"calls": v["calls"], "wall": v["wall"], "distinct": len(v["keys"])}
        for k, v in rec.helper.items()
    }
    return {
        "molecule": name,
        "path": path,
        "wall_total": wall,
        "error": err,
        "oin_sha": None if oin is None else hashlib.sha256(oin.encode()).hexdigest()[:16],
        "ac2bo_calls": rec.calls,
        "helpers": helper,
        "ac2bo_stats": dict(loc.AC2BO_STATS),
        "resonance": rec.resonance,
    }


def report(res):
    print(f"\n=== {res['molecule']}  ({res['path']}) ===")
    # The A/B payload, printed per molecule rather than only written to --json-out at the
    # end: a run killed partway (these molecules can take 20 minutes each) must still leave
    # a comparable arm behind. Counters accompany the sha because they localise a
    # divergence -- a different sha with identical `candidates` is a downstream difference,
    # not a valence-search one.
    st = res["ac2bo_stats"]
    print(
        f"AB oin_sha={res['oin_sha']} candidates={st['candidates']} "
        f"matching={st['matching_calls']} found_valid={st['found_valid']} "
        f"ac2bo_calls={st['ac2bo_calls']} over_cap={st['over_cap_calls']}"
    )
    print(f"wall_total (ADVISORY, host contended) = {res['wall_total']:.2f}s  err={res['error']}")
    calls = res["ac2bo_calls"]
    core_sum = sum(c["wall_core"] for c in calls)
    ord_sum = sum(c["ordered_wall"] or 0.0 for c in calls)
    print(
        f"_AC2BO_core calls={len(calls)}  sum(wall_core)={core_sum:.2f}s "
        f"({100 * core_sum / res['wall_total'] if res['wall_total'] else 0:.1f}% of encode)"
    )
    print(
        f"{'#':>3} {'atoms':>6} {'chg':>4} {'carb':>5} {'cap':>4} {'cands':>8} {'match':>8} "
        f"{'valid':>5} {'ord_len':>9} {'ord_s':>8} {'core_s':>8} {'ord%':>6}"
    )
    for i, c in enumerate(calls, 1):
        ol = "-" if c["ordered_len"] is None else str(c["ordered_len"])
        ow = 0.0 if c["ordered_wall"] is None else c["ordered_wall"]
        pct = 100 * ow / c["wall_core"] if c["wall_core"] else 0.0
        print(
            f"{i:>3} {c['n_atoms']:>6} {c['charge']:>4} {str(c['allow_carbenes']):>5} "
            f"{'OVER' if c['over_cap'] else 'sub':>4} {c['candidates']:>8} {c['matching']:>8} "
            f"{str(c['found_valid']):>5} {ol:>9} {ow:>8.3f} {c['wall_core']:>8.3f} {pct:>5.1f}%"
        )
    print(f"\n{'helper':<34} {'calls':>7} {'distinct':>9} {'wall_s':>9} {'% encode':>9}")
    for k in sorted(res["helpers"], key=lambda k: -res["helpers"][k]["wall"]):
        v = res["helpers"][k]
        pc = 100 * v["wall"] / res["wall_total"] if res["wall_total"] else 0.0
        print(f"{k:<34} {v['calls']:>7} {v['distinct']:>9} {v['wall']:>9.3f} {pc:>8.1f}%")
    print(
        f"\n_ordered_valences total = {ord_sum:.3f}s "
        f"({100 * ord_sum / res['wall_total'] if res['wall_total'] else 0:.2f}% of encode, "
        f"{100 * ord_sum / core_sum if core_sum else 0:.2f}% of _AC2BO_core)"
    )
    r = res.get("resonance") or {}
    if r:
        pc = 100 * r["wall"] / res["wall_total"] if res["wall_total"] else 0.0
        print(
            f"\nFORKED RESONANCE: lig_checks={r['lig_checks']} forks={r['calls']} "
            f"status={r['status']} wall={r['wall']:.2f}s ({pc:.1f}% of encode)"
        )
    print("AC2BO_STATS:", json.dumps(res["ac2bo_stats"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tmCAT-tmPHOTO_xyz_dataset")
    ap.add_argument("--molecules", required=True, help="comma-separated, or @file")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    raw = (
        open(args.molecules[1:]).read().split()
        if args.molecules.startswith("@")
        else args.molecules.split(",")
    )
    names = [x.strip() for x in raw if x.strip()]
    print(f"loadavg_at_start={os.getloadavg()}", flush=True)

    out = []
    for name in names:
        try:
            path = find_xyz(args.dataset, name)
        except FileNotFoundError as exc:
            print(f"SKIP {name}: {exc}", flush=True)
            continue
        res = run_one(path, name)
        out.append(res)
        report(res)
        print(f"#DONE {len(out)}", flush=True)

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"\nwrote {args.json_out}")
    print(f"#DONE_ALL {len(out)}/{len(names)}", flush=True)


if __name__ == "__main__":
    main()
