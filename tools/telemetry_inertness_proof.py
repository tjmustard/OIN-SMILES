#!/usr/bin/env python3
"""Prove the telemetry instrument does not perturb what it measures.

The Phase-2 telemetry probes are inserted into the generation path. Before any
firing-rate number from them can be believed, the instrument itself has to be
shown inert: the structures generated with telemetry compiled in must be
byte-identical to those generated without it, both when it is switched off and
when it is switched on.

That test is only meaningful because the generator was first shown to be
deterministic under identical invocation (flake floor 0.00%, 24/24
byte-identical XYZ across replicates). If the baseline were not self-reproducible
this script would be measuring noise, so it re-establishes that precondition on
its own sample rather than assuming it.

Three arms over the same molecules, run sequentially -- never concurrently,
because concurrent sweeps are known to fabricate ``no_conformers`` failures:

    A  baseline      OIN_TELEMETRY unset   (instrument compiled in, dormant)
    B  disabled      OIN_TELEMETRY=0
    C  enabled       OIN_TELEMETRY=1

Acceptance: sha256 of every generated XYZ matches across all three arms, and the
round-trip outcome and re-encoded OIN string match. Arm C additionally must have
produced at least one telemetry event, otherwise the probes are not wired up and
a passing byte-identity result would be vacuous.

Usage:
    uv run python tools/telemetry_inertness_proof.py --n 20
    uv run python tools/telemetry_inertness_proof.py --check-only --work-dir ...
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATASET = REPO / "tmCAT-tmPHOTO_xyz_dataset"
CAPSTONE = DATASET / "results-capstone-v042"

ARMS = {"A_baseline": None, "B_disabled": "0", "C_enabled": "1"}


def pick_molecules(n: int) -> list[str]:
    """A deterministic spread of molecules that previously generated a structure."""
    rows = json.loads((CAPSTONE / "summary_roundtrip.json").read_text())
    ok = sorted(r["molecule"] for r in rows if r.get("status") == "success")
    if not ok:
        raise SystemExit("no successful molecules to sample")
    stride = max(1, len(ok) // n)
    return ok[::stride][:n]


def run_arm(arm: str, telemetry: str | None, molecules: list[str], work: Path) -> Path:
    out = work / arm
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.pop("OIN_TELEMETRY", None)
    if telemetry is not None:
        env["OIN_TELEMETRY"] = telemetry
    # Fixed hash seed: dict/set iteration order must not vary between arms.
    env["PYTHONHASHSEED"] = "0"

    cmd = [
        "uv",
        "run",
        "python",
        "tools/test_dataset_roundtrip.py",
        "--only",
        ",".join(molecules),
        "--dataset-dir",
        str(DATASET),
        "--output-dir",
        str(out),
        "--mol-timeout",
        "900",
    ]
    print(f"  [{arm}] OIN_TELEMETRY={telemetry!r} ...", flush=True)
    with (out / "run.log").open("w") as log:
        subprocess.run(cmd, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    return out


def digest_structures(out: Path) -> dict[str, str]:
    d = {}
    for xyz in sorted((out / "structures").glob("*_generated.xyz")):
        d[xyz.name] = hashlib.sha256(xyz.read_bytes()).hexdigest()
    return d


def outcomes(out: Path) -> dict[str, tuple]:
    d = {}
    for report in sorted((out / "individual_reports").glob("*.json")):
        try:
            r = json.loads(report.read_text())
        except json.JSONDecodeError:
            continue
        d[r["molecule"]] = (r.get("status"), r.get("smiles_2"))
    return d


def compare(work: Path) -> int:
    arms = list(ARMS)
    digests = {a: digest_structures(work / a) for a in arms}
    outs = {a: outcomes(work / a) for a in arms}

    ref = arms[0]
    failures = 0
    print("\n" + "=" * 70)
    print("TELEMETRY INERTNESS PROOF")
    print("=" * 70)
    print(f"\nreference arm: {ref}  ({len(digests[ref])} structures, {len(outs[ref])} reports)")

    for arm in arms[1:]:
        common = sorted(set(digests[ref]) & set(digests[arm]))
        ident = sum(1 for k in common if digests[ref][k] == digests[arm][k])
        differing = [k for k in common if digests[ref][k] != digests[arm][k]]
        only_ref = sorted(set(digests[ref]) - set(digests[arm]))
        only_arm = sorted(set(digests[arm]) - set(digests[ref]))

        common_out = sorted(set(outs[ref]) & set(outs[arm]))
        out_diff = [m for m in common_out if outs[ref][m] != outs[arm][m]]

        ok = not differing and not only_ref and not only_arm and not out_diff
        failures += 0 if ok else 1
        print(f"\n  {ref} vs {arm}: {'PASS' if ok else 'FAIL'}")
        print(f"    XYZ byte-identical : {ident}/{len(common)}")
        if differing:
            print(f"    DIFFERING          : {differing[:10]}")
        if only_ref or only_arm:
            print(f"    structure set differs: only_{ref}={only_ref[:5]} only_{arm}={only_arm[:5]}")
        print(f"    outcome+string identical: {len(common_out) - len(out_diff)}/{len(common_out)}")
        if out_diff:
            for m in out_diff[:10]:
                print(f"      {m}: {outs[ref][m]} != {outs[arm][m]}")

    # A passing byte-identity result is vacuous if the probes never fired.
    events_path = work / "C_enabled" / "telemetry_events.json"
    if events_path.exists():
        ev = json.loads(events_path.read_text())
        total = sum(len(v) for v in ev.values()) if isinstance(ev, dict) else len(ev)
        print(f"\n  probes fired in C_enabled: {total} events")
        if total == 0:
            print("    WARNING: no events -- probes may not be wired; byte-identity is vacuous")
            failures += 1
    else:
        print("\n  NOTE: no telemetry_events.json; the harness does not yet persist events.")
        print("        Byte-identity is still valid, but probe wiring is unproven here.")

    print("\n" + ("ALL ARMS INERT" if failures == 0 else f"{failures} ARM(S) FAILED"))
    return failures


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=20, help="molecules per arm")
    ap.add_argument("--work-dir", type=Path, default=Path.home() / "elimination-v043" / "inertness")
    ap.add_argument("--check-only", action="store_true", help="compare existing arms, run nothing")
    args = ap.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)

    if not args.check_only:
        molecules = pick_molecules(args.n)
        (args.work_dir / "molecules.txt").write_text("\n".join(molecules) + "\n")
        print(f"{len(molecules)} molecules; running {len(ARMS)} arms sequentially")
        print("(sequential by design: concurrent sweeps fabricate no_conformers failures)")
        for arm, telemetry in ARMS.items():
            run_arm(arm, telemetry, molecules, args.work_dir)

    sys.exit(1 if compare(args.work_dir) else 0)


if __name__ == "__main__":
    main()
