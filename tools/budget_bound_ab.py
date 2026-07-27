"""A/B the enforced generation bound: asked for B seconds, how many were spent?

WHY NOT MEASURE THIS THROUGH THE SWEEP HARNESS
==============================================
``tools/test_dataset_roundtrip.py`` reports ``metrics.elapsed_s``, which is a **sum over
up to three separately SIGKILLed attempts** (PASS 1 ``UFF_1``, PASS 2 ``tier1``,
``tier5``). A bound is a property of ONE generation call, so measuring it through a field
that adds three of them together -- under an outer kill that would mask an overrun anyway
-- cannot answer the question. That conflation is what produced this release's refuted
premise; see ``docs/agentic-notes/v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md``.

So this calls ``OIN3DGenerator(timeout=B).generate(oin)`` directly, once per molecule per
arm, and records the wall-clock actually spent. ``eps = spent - B`` is then exactly the
quantity the release has to name.

ONE SUBPROCESS PER (MOLECULE, ARM)
==================================
``OIN_ENFORCE_BUDGET`` is read once per generation, but several sibling ``OIN_*`` levers
and module-level constants are frozen at IMPORT time, and the generator carries
process-lifetime caches (the PuLP topology memo, embed pool state). Flipping a lever
inside one interpreter and re-running would measure the caches as much as the lever. A
fresh interpreter per measurement is the only isolation that does not depend on
enumerating every such cache correctly -- the same discipline ``tools/gate_v047.sh`` ARM 2
already follows.

THE OFF ARM NEEDS ITS OWN KILL
==============================
With the lever off there is by construction no bound, so the OFF arm is run under an
outer ``SIGKILL`` (``--off-arm-kill``, default 6x the budget). A row killed that way is
recorded as ``killed`` with the kill time -- **a lower bound on its true spend, not its
spend.** Reporting it as if it were the real number would understate the very defect
being measured, so the summary counts and labels it separately.

Emits newline-delimited JSON, one object per (molecule, arm), each written with a real
append+flush so a kill mid-run cannot discard output, followed by ``#DONE <n>``.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/budget_bound_ab.py \
        --cohort-dir <DS>/cohort-v049-strata --budget 30 --limit 60 \
        --out results-v0.4.9-bound/ab.jsonl
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

_CHILD = r"""
import json, os, sys, time
sys.path.insert(0, os.environ["OIN_SRC"])
from oinsmiles import XYZToSMILES
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen as G

xyz, budget = sys.argv[1], float(sys.argv[2])
rec = {"molecule": os.path.basename(xyz)[:-4], "budget_s": budget}
t0 = time.monotonic()
try:
    oin = XYZToSMILES().convert(xyz)
except Exception as e:
    rec.update(outcome="encode_fail", spent_s=round(time.monotonic() - t0, 3),
               error=f"{type(e).__name__}: {e}"[:200])
    print(json.dumps(rec), flush=True); raise SystemExit(0)
rec["encode_s"] = round(time.monotonic() - t0, 3)

t1 = time.monotonic()
try:
    res = G(optimizer=None, ensemble_size=1, timeout=budget, ff_params=None).generate(oin)
    rec.update(outcome="ok" if getattr(res, "xyz", None) else "empty",
               xyz_len=len(getattr(res, "xyz", "") or ""))
except Exception as e:
    rec.update(outcome="raised", error_type=type(e).__name__,
               error=f"{type(e).__name__}: {e}"[:200])
rec["gen_s"] = round(time.monotonic() - t1, 3)
rec["spent_s"] = round(time.monotonic() - t0, 3)
print(json.dumps(rec), flush=True)
"""


def run_one(py, src, xyz, budget, enforce, kill_after):
    env = dict(os.environ, OIN_SRC=src)
    # lever_enabled() semantics: "0" DISABLES. Never set it to "0" expecting truthiness.
    env["OIN_ENFORCE_BUDGET"] = "1" if enforce else "0"
    t0 = time.monotonic()
    try:
        out = subprocess.run(
            [py, "-c", _CHILD, xyz, str(budget)],
            capture_output=True,
            text=True,
            timeout=kill_after,
            env=env,
        )
        for line in out.stdout.splitlines():
            line = line.strip()
            if line.startswith("{"):
                return json.loads(line)
        return {
            "molecule": os.path.basename(xyz)[:-4],
            "budget_s": budget,
            "outcome": "no_output",
            "spent_s": round(time.monotonic() - t0, 3),
            "stderr_tail": (out.stderr or "")[-300:],
        }
    except subprocess.TimeoutExpired:
        # NOT the molecule's spend -- a lower bound on it. Labelled so the summary
        # cannot quietly average it in with real measurements.
        return {
            "molecule": os.path.basename(xyz)[:-4],
            "budget_s": budget,
            "outcome": "killed",
            "spent_s": round(time.monotonic() - t0, 3),
            "spent_is_lower_bound": True,
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--budget", type=float, default=30.0, help="OIN3DGenerator(timeout=)")
    ap.add_argument("--limit", type=int, default=0, help="first N molecules (0 = all)")
    ap.add_argument("--names-file", default=None, help="restrict to these basenames")
    ap.add_argument(
        "--off-arm-kill",
        type=float,
        default=None,
        help="outer SIGKILL for the UNBOUNDED arm (default 6x budget). The ON arm gets "
        "the same kill so both arms are measured under identical outer conditions -- an "
        "arm measured without a kill is not comparable to one measured with it.",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    kill_after = args.off_arm_kill or (6.0 * args.budget)
    src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
    py = sys.executable

    names = sorted(f[:-4] for f in os.listdir(args.cohort_dir) if f.endswith(".xyz"))
    if args.names_file:
        keep = {
            ln.strip().removesuffix(".xyz")
            for ln in open(args.names_file)
            if ln.strip() and not ln.startswith("#")
        }
        names = [n for n in names if n in keep]
    if args.limit:
        names = names[: args.limit]
    if not names:
        sys.exit("error: 0 molecules selected -- refusing an empty A/B")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    print(
        f"# {len(names)} molecules x 2 arms, budget={args.budget}s, outer kill={kill_after}s",
        file=sys.stderr,
    )
    n = 0
    with open(args.out, "w") as log:
        for i, name in enumerate(names, 1):
            xyz = os.path.join(args.cohort_dir, name + ".xyz")
            for enforce in (False, True):
                rec = run_one(py, src, xyz, args.budget, enforce, kill_after)
                rec["arm"] = "ON" if enforce else "OFF"
                log.write(json.dumps(rec) + "\n")
                log.flush()
                os.fsync(log.fileno())
                n += 1
            print(
                f"[{i}/{len(names)}] {name}: "
                f"OFF/ON recorded ({rec.get('outcome')}, {rec.get('spent_s')}s)",
                file=sys.stderr,
            )
        log.write(f"#DONE {n}\n")
    print(f"#DONE {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
