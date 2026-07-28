"""How often does the cheap acceptance prefilter veto a conformer the strict test would take?
(v0.4.13 Lane 1)

THE DEFECT
==========
``_reencode_key_matches`` step 1 re-serializes the generated geometry through the GENERATOR'S OWN
contract-mol connectivity and rejects on a key mismatch, justified in-code as *"a MISMATCH here is
a reliable 'geometry is wrong' signal"*.

``AROHIA_comp_0`` falsifies that: the cheap test matches **0 of 48** conformers while the strict
independent test matches **16 of 48**. Because the cheap ``return False`` fires first, those 16
are unreachable **in both arms of every A/B ever run on this molecule**. The error runs in the
**pessimistic** direction -- it makes the project look worse than it is -- which is why it
survived five releases without anyone chasing it.

``OIN_PREFILTER_ADVISORY`` makes step 1 advisory. This tool measures what that is worth, and what
it costs.

WHAT A BROKEN VERSION OF THIS WOULD PRINT
=========================================
**Zero overrides** -- which is also what a working run prints when the population is genuinely
empty. The two are indistinguishable from the output alone, so this tool refuses to let them be:

* Every arm reports ``overridden`` **beside** ``confirmed`` and ``cheap_pass``. An arm where all
  three are zero did not measure a small effect; it did not run the predicate at all.
* ``--molecule AROHIA_comp_0`` is the two-point fixture, where the answer is known independently
  (cheap 0/48, strict 16/48). A zero there means the lever is not wired, and nothing else this
  tool prints can be trusted.

THE COST IS A DELIVERABLE, NOT A FOOTNOTE
=========================================
*Do not fix a prefilter by removing it.* The prefilter exists to make acceptance cheap; every
cheap-veto this lever overrides now pays a full ``XYZToSMILES().convert()`` round trip, measured
at 48-57 s per call on an eta/haptic conformer. A correct-but-slow prefilter moves the cost into
v0.4.12's territory, so both arms are timed and the delta is reported next to the recovery count.

⚠ **Do not run this on a loaded machine.** v0.4.12 measured ``UQUXAG_comp_0`` at 17.93 s under
load and 11.06 s clean -- a 62% inflation, comparable to the whole effect being measured -- and
threw the loaded run away rather than banking it.

Usage
-----
    V=<main>/.venv/bin/python; export PYTHONPATH=$PWD/src
    $V tools/prefilter_prevalence.py --xyz <main>/tmCAT-tmPHOTO_xyz_dataset/cat/AROHIA_comp_0.xyz
    $V tools/prefilter_prevalence.py --cohort <dir> --out prevalence.json
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

SITES = (
    "adapter.prefilter_veto_overridden",
    "adapter.prefilter_veto_confirmed",
    "adapter.prefilter_cheap_pass",
    "adapter.prefilter_veto_advisory",
)


def run_one(xyz_path, advisory, timeout):
    """Generate once and return (counts, elapsed_s, passed, error)."""
    os.environ["OIN_PREFILTER_ADVISORY"] = "1" if advisory else "0"
    os.environ["OIN_TELEMETRY"] = "1"

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation import _telemetry
    from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen
    from oinsmiles.oin.compare import canonical_roundtrip_key

    t0 = time.monotonic()
    counts, passed, err = {}, None, None
    try:
        with _telemetry.collecting():
            oin_in = XYZToSMILES().convert(xyz_path)
            gen = OIN3DGeneratorMetallogen(
                optimizer=None, ensemble_size=1, timeout=timeout, ff_params=None
            )
            res = gen.generate(oin_in)
            counts = dict(_telemetry.counts())
        # ⚠ SCORED HONESTLY, and the obvious spelling is the WRONG one.
        #
        # The tempting line here is
        #     canonical_roundtrip_key(get_oin_string(res.mol, coords))
        # which is what the older A/B tools in this directory do -- and it is CIRCULAR: it scores
        # the round trip with `res.mol`, the GENERATOR'S OWN bond graph, i.e. exactly the artifact
        # that would have to be wrong for the test to fail. v0.4.8 measured that at 61/633 = 9.6%
        # false positives and replaced it with OIN_INDEP_SCORE.
        #
        # It matters specifically HERE, more than in a normal A/B: this lane's whole subject is a
        # disagreement between the cheap (generator-connectivity) predicate and the strict
        # (independently re-perceived) one. Scoring the outcome with the cheap predicate would
        # judge the lever by the very test it exists to override, and would report "recovered
        # nothing" by construction. So write the XYZ and re-perceive it, the same single call
        # tools/honest_rescore.py makes.
        if getattr(res, "xyz", None):
            with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as fh:
                fh.write(res.xyz)
                tmp = fh.name
            try:
                indep = XYZToSMILES().convert(tmp)
                passed = canonical_roundtrip_key(oin_in) == canonical_roundtrip_key(indep)
            finally:
                os.unlink(tmp)
    except Exception as exc:  # noqa: BLE001 -- a failed molecule is data, not a crash
        err = f"{type(exc).__name__}: {exc}"
    return counts, round(time.monotonic() - t0, 2), passed, err


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xyz", action="append", default=[], help="one input .xyz, repeatable")
    ap.add_argument("--cohort", help="directory of *.xyz")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--out", help="write per-molecule JSON here")
    args = ap.parse_args()

    paths = list(args.xyz)
    if args.cohort:
        paths += sorted(glob.glob(os.path.join(args.cohort, "*.xyz")))
    if args.limit:
        paths = paths[: args.limit]
    if not paths:
        sys.exit("FATAL: no input. Pass --xyz or --cohort.")

    rows = []
    tot = {"off": dict.fromkeys(SITES, 0), "on": dict.fromkeys(SITES, 0)}
    t_off = t_on = 0.0
    n_recovered = 0

    for i, p in enumerate(paths, 1):
        mol = os.path.splitext(os.path.basename(p))[0]
        c_off, e_off, pass_off, err_off = run_one(p, advisory=False, timeout=args.timeout)
        c_on, e_on, pass_on, err_on = run_one(p, advisory=True, timeout=args.timeout)
        for s in SITES:
            tot["off"][s] += c_off.get(s, 0)
            tot["on"][s] += c_on.get(s, 0)
        t_off += e_off
        t_on += e_on
        recovered = bool(pass_on) and not bool(pass_off)
        n_recovered += recovered
        rows.append(
            {
                "molecule": mol,
                "off": {"counts": c_off, "elapsed_s": e_off, "passed": pass_off, "error": err_off},
                "on": {"counts": c_on, "elapsed_s": e_on, "passed": pass_on, "error": err_on},
                "recovered": recovered,
            }
        )
        ov = c_on.get("adapter.prefilter_veto_overridden", 0)
        cf = c_on.get("adapter.prefilter_veto_confirmed", 0)
        print(
            f"[{i}/{len(paths)}] {mol:22s} overridden={ov:4d} confirmed={cf:4d} "
            f"pass {pass_off}->{pass_on}  {e_off:7.2f}s -> {e_on:7.2f}s",
            flush=True,
        )

    ov = tot["on"]["adapter.prefilter_veto_overridden"]
    cf = tot["on"]["adapter.prefilter_veto_confirmed"]
    cp = tot["on"]["adapter.prefilter_cheap_pass"]
    adv = tot["on"]["adapter.prefilter_veto_advisory"]

    print(f"\n## Lever ON, over {len(paths)} molecules — every counter, with its denominator\n")
    print(f"  cheap veto, strict ACCEPTS (overridden) : {ov:6d}   <- the AROHIA shape")
    print(f"  cheap veto, strict rejects (confirmed)  : {cf:6d}")
    print(f"  cheap veto total (advisory fired)       : {adv:6d}   (= overridden + confirmed)")
    print(f"  cheap PASSED, veto never consulted      : {cp:6d}")
    print(f"  molecules recovered (fail -> pass)      : {n_recovered:6d} / {len(paths)}")
    print(f"\n  wall clock: OFF {t_off:.1f}s   ON {t_on:.1f}s   delta {t_on - t_off:+.1f}s")
    print(
        "\n  ⚠ THESE DENOMINATORS ARE NOT POOL SIZES. The first override ACCEPTS, which stops the\n"
        "    pool filling, so the lever-ON arm evaluates FEWER conformers than the lever-OFF arm\n"
        "    and far fewer than the pool would hold. `AROHIA_comp_0`'s documented 0/48-vs-16/48 was\n"
        "    measured with the pool FORCED FULL (tools/probe_accept_gap.py); a count of 3 here is\n"
        "    the same defect observed until it stopped mattering, not a smaller one. Read\n"
        "    `overridden` as 'conformers recovered before early exit', never as a rate over 48.\n"
        "    For the same reason a NEGATIVE latency delta is expected and is not a speed claim."
    )

    print()
    if adv == 0 and cp == 0:
        print("  🔴 EVERY COUNTER IS ZERO. The predicate did not run -- this measured NOTHING.")
        print("     Check the lever is wired and telemetry is enabled before reading anything.")
        verdict = "INSTRUMENT_DEAD"
    elif adv == 0:
        print(f"  ⚠ The cheap prefilter never vetoed anything across {len(paths)} molecules")
        print(f"    ({cp} cheap passes). There is no population here to recover -- a real null.")
        verdict = "NO_POPULATION"
    elif ov == 0:
        print(f"  ✅ MEASURED NULL: the prefilter vetoed {adv} conformers and the strict test")
        print("     agreed with it every time. The veto is not costing conformers on this")
        print("     sample, and the lane should say so rather than ship a fix.")
        verdict = "PREFILTER_VINDICATED"
    else:
        print(f"  🔴 CONFIRMED: {ov} conformers were vetoed by the cheap prefilter that the")
        print(f"     strict independent test ACCEPTS ({100 * ov / adv:.1f}% of its {adv} vetoes).")
        print(f"     Those are unreachable in the shipped default. Cost: {t_on - t_off:+.1f}s.")
        verdict = "DEFECT_CONFIRMED"

    if args.out:
        json.dump(
            {
                "verdict": verdict,
                "n_molecules": len(paths),
                "totals": tot,
                "n_recovered": n_recovered,
                "elapsed_off_s": round(t_off, 2),
                "elapsed_on_s": round(t_on, 2),
                "rows": rows,
            },
            open(args.out, "w"),
            indent=2,
        )
        print(f"\nwrote {args.out}")
    print(f"\n#DONE {len(paths)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
