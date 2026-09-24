"""What the ARM 2 honest observation column says, over a full gate run. Reads only.

Every fresh ``OK`` row now ends with ``honest`` = sha256 of an independent ``XYZToSMILES().convert``
of the generated coordinates -- the harness's own predicate -- beside the gated field 3, which is
the same coordinates serialised through the generator's OWN bond graph. Nothing compares the new
column yet; this reads it, over ``on_rows.tsv`` (tools/v0417/arm2_refreeze.py rows):

    gated pass      sha_out == sha_in         (what ARM 2 has always scored)
    honest pass     honest  == sha_in         (what the harness scores)

and prints the 2x2, by eta, plus the rows the gate PASSES that the harness would FAIL -- the class
of defect ARM 2 could not see at the L2 promotion. Rows where the encode was not attempted (near
the hard timeout) or errored are counted, never dropped.

    <main>/.venv/bin/python tools/v0418/arm2_honest_column.py --run <results-…-arm2-refreeze>
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--run", type=Path, required=True, help="a re-freeze run dir holding on_rows.tsv"
    )
    args = ap.parse_args()
    lines = (args.run / "on_rows.tsv").read_text().splitlines()
    if not lines[-1].startswith("#DONE"):
        sys.exit("ABORT: on_rows.tsv has no #DONE trailer")
    rows = [ln.split("\t") for ln in lines if ln.strip() and not ln.startswith("#")]
    # golden, name, sha_in, sha_out, len_in, len_out, eta, xyz_sha, status, honest
    ok = [r for r in rows if len(r) > 8 and r[8] == "OK"]
    print(f"DENOMINATOR {len(rows)} rows; OK {len(ok)}")
    kinds = Counter(
        "column absent"
        if len(r) < 10 or not r[9]
        else "not attempted"
        if r[9].startswith("HONEST_NOT_ATTEMPTED")
        else "error"
        if r[9].startswith("HONEST_ERROR")
        else "measured"
        for r in ok
    )
    print("  honest column:", dict(kinds))
    meas = [r for r in ok if len(r) > 9 and r[9] and not r[9].startswith("HONEST_")]
    for label, sel in (
        ("all", meas),
        ("eta", [r for r in meas if r[6] == "eta"]),
        ("non-eta", [r for r in meas if r[6] != "eta"]),
    ):
        t = Counter((r[3] == r[2], r[9] == r[2]) for r in sel)
        print(
            f"  {label:7s} n={len(sel):3d}   gated pass & honest pass {t[(True, True)]:3d}   "
            f"gated PASS & honest FAIL {t[(True, False)]:3d}   gated fail & honest pass {t[(False, True)]:3d}   "
            f"both fail {t[(False, False)]:3d}"
        )
    blind = sorted(r[1] for r in meas if r[3] == r[2] and r[9] != r[2])
    print(f"  rows the gate passes that the harness would fail: {len(blind)}")
    print("   ", " ".join(m.replace("_comp_0", "") for m in blind))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
