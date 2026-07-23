#!/usr/bin/env python
"""SL5 encoder-robustness sub-triage.

Runs each ``encode_fail`` molecule from the v0.4.4 BASELINE through the XYZ->OIN
encoder in an *isolated subprocess with an OS-level timeout*, because the timeout
failures hang inside a C-level RDKit call (``ResonanceMolSupplier`` enumeration in
``lig_checks``) that a Python ``signal.alarm`` cannot preempt.

Classifies each molecule into a (stage, reason) bucket and freezes a markdown +
JSON table that is the SL5 worklist.

Usage:
    # driver (runs all 48, writes tools/sl5_triage.{md,json}):
    PYTHONPATH=src python tools/sl5_triage.py --dataset-dir <abs path to tmCAT-tmPHOTO_xyz_dataset>
    # single-molecule worker (internal; prints one JSON line):
    PYTHONPATH=src python tools/sl5_triage.py --worker <MOL> --dataset-dir <...>
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback

# The 48-molecule encode_fail cohort frozen in spec/handoffs/v0.4.4/BASELINE.md.
COHORT = [
    "ASISAX", "AVOFIB", "BEKLUA", "BEKMIP", "BENVOG", "CAKBEW", "CAKBOG", "COZCEZ",
    "FAQYUU", "GANYEZ", "GOHWOQ", "HAXJAS", "HAXJOG", "HICLAG", "HOCVAY", "HOHKUL",
    "HUCNAU", "ICEZIC", "JABGAX", "JAFMIP", "JAFTAO", "JAFTES", "KAXVOX", "KAXWAK",
    "KEMTED", "KESWUB", "LEZWAO", "MAFSIY", "MODZUA", "NAKLET", "OZAREO", "PAQBOZ",
    "PAQCAM", "PAYTUH", "RANCIU", "RANMUR", "RAWJEG", "RIWKAK", "RIWKEO", "RONPES",
    "RONQET", "RONQOD", "RULBUV", "ULODUU", "WEFZAL", "XUKRIF", "YIBZIV", "YIVLAQ",
]

PER_MOL_TIMEOUT_S = 90


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def boron_count(xyz_path: str) -> int:
    n = 0
    with open(xyz_path) as fh:
        for i, line in enumerate(fh):
            if i < 2:
                continue
            tok = line.split()
            if tok and tok[0] == "B":
                n += 1
    return n


def run_worker(mol: str, dataset_dir: str) -> dict:
    """Convert one molecule; return a structured result dict (no timeout here)."""
    import warnings

    warnings.filterwarnings("ignore")
    from oinsmiles import XYZToSMILES

    path = find_xyz(mol, dataset_dir)
    if path is None:
        return {"mol": mol, "status": "NO_FILE"}
    try:
        oin = XYZToSMILES().convert(path)
        return {"mol": mol, "status": "OK", "oin": oin}
    except Exception as e:  # noqa: BLE001 - triage captures everything
        tb = traceback.extract_tb(sys.exc_info()[2])
        last = None
        for fr in tb:
            if "oinsmiles" in fr.filename:
                last = fr
        loc = f"{os.path.basename(last.filename)}:{last.lineno}:{last.name}" if last else "?"
        return {
            "mol": mol,
            "status": "FAIL",
            "etype": type(e).__name__,
            "loc": loc,
            "msg": str(e).replace("\n", " "),
        }


def classify(rec: dict, nB: int) -> str:
    if rec["status"] == "OK":
        return "encodes_now"
    if rec["status"] == "TIMEOUT":
        return "resonance_timeout"
    if rec["status"] == "NO_FILE":
        return "no_file"
    msg = rec.get("msg", "")
    loc = rec.get("loc", "")
    if nB >= 3:
        return "boron_cluster"
    if "get_lig_mol failed" in msg:
        return "perception_charge_gap"
    if "aligner" in loc or "oin_aligner" in loc or "inline" in loc:
        return "aligner_serializer"
    return "other"


def drive(dataset_dir: str, out_stem: str) -> None:
    results = []
    for mol in COHORT:
        path = find_xyz(mol, dataset_dir)
        nB = boron_count(path) if path else 0
        cmd = [
            sys.executable, os.path.abspath(__file__),
            "--worker", mol, "--dataset-dir", dataset_dir,
        ]
        env = dict(os.environ)
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=PER_MOL_TIMEOUT_S, env=env,
            )
            line = next(
                (ln for ln in proc.stdout.splitlines() if ln.startswith("{")), None,
            )
            rec = json.loads(line) if line else {
                "mol": mol, "status": "FAIL", "etype": "NoOutput",
                "loc": "?", "msg": (proc.stderr[-200:] or "no json line"),
            }
        except subprocess.TimeoutExpired:
            rec = {"mol": mol, "status": "TIMEOUT"}
        rec["nB"] = nB
        rec["bucket"] = classify(rec, nB)
        results.append(rec)
        print(f"{mol:14s} nB={nB:<3} {rec['bucket']:22s} {rec['status']}", flush=True)

    with open(f"{out_stem}.json", "w") as fh:
        json.dump(results, fh, indent=2)

    from collections import Counter

    tally = Counter(r["bucket"] for r in results)
    lines = ["# SL5 sub-triage — encode_fail cohort (48)\n"]
    lines.append("| bucket | count |")
    lines.append("|---|---:|")
    for b, n in tally.most_common():
        lines.append(f"| `{b}` | {n} |")
    lines.append(f"| **total** | **{len(results)}** |\n")
    lines.append("| molecule | nB | bucket | status | detail |")
    lines.append("|---|---:|---|---|---|")
    for r in sorted(results, key=lambda r: (r["bucket"], r["mol"])):
        detail = r.get("oin", "") if r["status"] == "OK" else r.get("msg", r["status"])
        detail = str(detail)[:80].replace("|", r"\|")
        lines.append(
            f"| {r['mol']} | {r['nB']} | `{r['bucket']}` | {r['status']} | {detail} |"
        )
    with open(f"{out_stem}.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n=== TALLY ===", flush=True)
    for b, n in tally.most_common():
        print(f"{n:3d}  {b}", flush=True)
    print(f"\nWrote {out_stem}.md and {out_stem}.json", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", metavar="MOL", default=None)
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--out-stem", default="tools/sl5_triage")
    args = ap.parse_args()
    if args.worker:
        print(json.dumps(run_worker(args.worker, args.dataset_dir)), flush=True)
    else:
        drive(args.dataset_dir, args.out_stem)


if __name__ == "__main__":
    main()
