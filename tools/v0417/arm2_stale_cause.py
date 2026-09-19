"""Which lever left an ARM 2 golden row stale? One lever off at a time, fresh process each.

``arm2_field2_audit.py`` says a row's field 2 was already wrong before v0.4.17 (class ``STALE``,
``STALE_AND_MOVED`` or ``HEALED``). It cannot say WHEN. The stored sweeps bracket it -- reproduced by
the v0.4.8 sweeps, by none from v0.4.14 on -- and several levers were promoted in between.

This settles it by intervention. The base is the encoder as it was before this lane
(``OIN_EXACT_DONOR_FOLD=0``); then each default-ON lever is set to "0" ALONE, plus the one pair that
was promoted together and is documented as coupled (donor fold + parity veto). A lever is the cause
of a row iff switching it off brings the golden's hash back.

WHAT A BROKEN VERSION WOULD PRINT: "no lever restores it" for every row -- which is also what a
worker importing the wrong ``oinsmiles`` prints. So the base arm must reproduce the audit's ``off``
hash for every row, and the run aborts if it does not.

    <main>/.venv/bin/python tools/v0417/arm2_stale_cause.py --out <jsonl>
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve()
SRC = HERE.parents[2] / "src"
AUDIT_TOOL = HERE.with_name("arm2_field2_audit.py")
DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
OUT = DATA / "results-v0.4.17-exactfold" / "arm2_refreeze"
COHORTS = {"v047": DATA / "cohort-v047-slow100", "v049": DATA / "cohort-v049-strata"}
BASE = {"OIN_EXACT_DONOR_FOLD": "0"}
PAIR = ("OIN_CANONICAL_DONOR_FOLD", "OIN_FOLD_PARITY_VETO")


def _levers():
    sys.path.insert(0, str(SRC))
    from oinsmiles.oin.levers import default_on

    return sorted(set(default_on()) - set(BASE))


def _encode(xyz: Path, extra: dict):
    env = {k: v for k, v in os.environ.items() if not k.startswith("OIN_")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    env.update(BASE)
    env.update(extra)  # written as "0", never deleted: unset means the shipped default
    p = subprocess.run(
        [sys.executable, str(AUDIT_TOOL), "--one", str(xyz)],
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith("{")]
    return json.loads(lines[-1]) if lines else {"error": p.stderr.strip()[-200:]}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()

    rows = {}
    for tag in COHORTS:
        lines = (OUT / f"field2_audit_{tag}.jsonl").read_text().splitlines()
        if not lines[-1].startswith("#DONE"):
            sys.exit(f"ABORT: the {tag} audit has no #DONE trailer")
        for r in map(json.loads, lines[:-1]):
            if r["class"] in ("STALE", "STALE_AND_MOVED", "HEALED"):
                rows.setdefault(r["molecule"], (tag, r))  # JUXPID is in both; one input file
    if not rows:
        sys.exit("ABORT: no stale rows in the audits -- refusing an empty denominator")

    levers = _levers()
    arms = (
        [("base", {})]
        + [(lv, {lv: "0"}) for lv in levers]
        + [("+".join(PAIR), dict.fromkeys(PAIR, "0"))]
    )
    jobs = [(m, name, extra) for m in sorted(rows) for name, extra in arms]

    def run(job):
        m, name, extra = job
        return m, name, _encode(COHORTS[rows[m][0]] / f"{m}.xyz", extra)

    with ThreadPoolExecutor(args.jobs) as pool:
        got = list(pool.map(run, jobs))

    foreign = {r.get("src") for _m, _n, r in got} - {None, str(SRC)}
    if foreign:
        sys.exit(f"ABORT: a worker imported oinsmiles from {foreign}, not {SRC}")

    by = {}
    for m, name, r in got:
        by.setdefault(m, {})[name] = r.get("sha") or r.get("error")
    cause = Counter()
    with open(args.out, "w") as fh:
        for m in sorted(by):
            audit = rows[m][1]
            if by[m]["base"] != audit["off"]["sha"]:
                sys.exit(f"ABORT: {m}: the base arm does not reproduce the audit's `off` hash")
            restores = [n for n, h in by[m].items() if h == audit["golden_sha1"]]
            moves = [n for n, h in by[m].items() if n != "base" and h != by[m]["base"]]
            cause[" | ".join(restores) or "NOTHING RESTORES IT"] += 1
            fh.write(
                json.dumps(
                    {
                        "molecule": m,
                        "audit_class": audit["class"],
                        "golden_sha1": audit["golden_sha1"],
                        "restored_by": restores,
                        "moved_by": moves,
                        "sha_by_arm": by[m],
                    }
                )
                + "\n"
            )
            print(
                f"  {m:16s} {audit['class']:16s} restored by: {', '.join(restores) or '-- nothing --'}"
            )
        fh.write(f"#DONE {len(by)}\n")
    print(f"DENOMINATOR {len(by)} stale rows x {len(arms)} arms   src={SRC}")
    for k, c in cause.most_common():
        print(f"  {c:3d}  {k}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
