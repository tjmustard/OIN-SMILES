"""Did a GENERATOR-side promotion fire inside the ARM 2 gate run -- and did the gate see it? Reads only.

``/refreeze-goldens`` for v0.4.18 re-froze 0 of 425 rows. A zero is a finding only once the lever is
shown to FIRE, and ARM 2's goldens cannot show it: they hold sha256(smiles_1) and sha256(smiles_2),
and the runner computes smiles_2 through the generator's OWN bond graph. But every FRESH row also
carries ``xyz_sha`` (an observation column), and the v0.4.17 re-freeze left a full run of the same
428 rows on pre-lever code. Generation is seeded, so:

    structure sha differs between the two full runs   -> the levers changed what was BUILT
    field 3 differs between the two full runs         -> the gate SAW it

    <main>/.venv/bin/python tools/v0418/arm2_lever_fired.py [--now <dir>] [--before <dir>]
        # defaults: the L2 promotion's full run against the v0.4.17 one; the release re-freeze
        # passes --now results-v0.4.18-release-arm2-refreeze --before results-v0.4.18-arm2-refreeze
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DATA = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
NOW = DATA / "results-v0.4.18-arm2-refreeze" / "on_rows.tsv"
BEFORE = DATA / "results-v0.4.17-exactfold" / "arm2_refreeze" / "on_rows.tsv"


def _rows(path: Path):
    lines = path.read_text().splitlines()
    if not lines[-1].startswith("#DONE"):
        sys.exit(f"ABORT: {path} has no #DONE trailer")
    out = {}
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split("\t")  # golden, name, sha_in, sha_out, len_in, len_out, eta, xyz_sha, status
        out[(f[0], f[1])] = f
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--now", type=Path, default=NOW.parent)
    ap.add_argument("--before", type=Path, default=BEFORE.parent)
    args = ap.parse_args()
    now_p, before_p = args.now / "on_rows.tsv", args.before / "on_rows.tsv"
    now, before = _rows(now_p), _rows(before_p)
    if set(now) != set(before) or not now:
        sys.exit(f"ABORT: the two runs do not hold the same rows ({len(now)} vs {len(before)})")
    ok = [
        k for k in now if now[k][8] == "OK" and before[k][8] == "OK" and now[k][7] and before[k][7]
    ]
    print(f"DENOMINATOR {len(now)} rows in both full runs; built OK in both: {len(ok)}")
    print(f"  pre-lever run : {before_p.read_text().splitlines()[0].split('--')[-1].strip()}")
    print(f"  promoted run  : {now_p.read_text().splitlines()[0].split('--')[-1].strip()}")
    for label, ks in (
        ("eta    ", [k for k in ok if now[k][6] == "eta"]),
        ("non-eta", [k for k in ok if now[k][6] != "eta"]),
    ):
        built = sum(now[k][7] != before[k][7] for k in ks)
        f3 = sum(now[k][3] != before[k][3] for k in ks)
        f2 = sum(now[k][2] != before[k][2] for k in ks)
        print(
            f"  {label} n={len(ks):3d}   structure differs {built:3d}   field 3 (gated) differs {f3}"
            f"   field 2 differs {f2}"
        )
    status = sum(now[k][8].split(":")[0] != before[k][8].split(":")[0] for k in now)
    print(f"  rows whose status class changed between the runs: {status}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
