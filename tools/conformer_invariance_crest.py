#!/usr/bin/env python3
"""Verify conformer invariance of OIN-SMILES using CREST-generated ensembles.

For each structure in the curated test set, run CREST to generate a conformer
ensemble, encode every conformer to OIN-SMILES, and check that they all collapse
to one canonical string. Unlike the hermetic rotation test (which only rotates a
fixed geometry), this exercises *genuinely different geometries* of the same
molecule.

CREST is an OPTIONAL external binary and NOT a dependency of this project. If the
``crest`` binary is not on PATH this script prints a notice and exits 0 (it never
fails a run just because CREST is absent), mirroring the ``shutil.which`` gate used
for the optional xtb optimizer.

Match rule: the primary gate is a **byte-identical** raw OIN string. As diagnostics
we also compute the tolerant ``winding_canonical_key`` and an independent
covalent-radius connectivity signature, so a divergence caused by CREST altering
metal-ligand bonding (expected chemistry) is reported separately from a genuine
canonicalization failure (same connectivity, different string).

Usage::

    # install CREST first (see README "Optional: CREST conformer cross-check"), then:
    python tools/conformer_invariance_crest.py --only CisPlatin,TransPlatin
    python tools/conformer_invariance_crest.py            # whole manifest
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "tests" / "fixtures" / "conformer_set" / "manifest.json"
DEFAULT_OUT = REPO / "conformer_invariance_results"

# Covalent radii (Angstrom, Cordero 2008 subset) for an encoder-independent
# connectivity fingerprint. Unknown elements fall back to 1.5.
COVALENT_R = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
    "Se": 1.20,
    "As": 1.19,
    "Sc": 1.70,
    "Ti": 1.60,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.61,
    "Fe": 1.52,
    "Co": 1.50,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Y": 1.90,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "Hf": 1.75,
    "Ta": 1.70,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
}


def read_frames(path: Path):
    """Parse a (possibly multi-frame) XYZ file into [(comment, elements, coords)]."""
    lines = path.read_text().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].strip())
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        elements, coords = [], []
        for ln in lines[i + 2 : i + 2 + n]:
            parts = ln.split()
            elements.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        frames.append((comment, elements, np.asarray(coords, dtype=float)))
        i += 2 + n
    return frames


def write_xyz(path: Path, comment, elements, coords):
    with open(path, "w") as f:
        f.write(f"{len(elements)}\n{comment}\n")
        for el, c in zip(elements, coords):
            f.write(f"{el} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n")


def connectivity(elements, coords, tol=1.3):
    """Encoder-independent bond set from covalent radii (frozenset of index pairs)."""
    n = len(elements)
    radii = np.array([COVALENT_R.get(e, 1.5) for e in elements])
    edges = set()
    for i in range(n):
        d = np.linalg.norm(coords[i + 1 :] - coords[i], axis=1)
        thr = tol * (radii[i] + radii[i + 1 :])
        for off in np.where(d <= thr)[0]:
            edges.add((i, i + 1 + int(off)))
    return frozenset(edges)


def _as_text(x):
    if x is None:
        return ""
    return x.decode(errors="replace") if isinstance(x, bytes) else x


def preopt(src: Path, pdir: Path, charge: int, mult: int, args, log_path=None) -> tuple:
    """Pre-optimize the input geometry with xtb; return (xtbopt.xyz | None, note).

    Falls back to the original geometry (returns None) if xtb fails/times out. Atom
    order is preserved by xtb, which the OIN canonicalization relies on.
    """
    pdir.mkdir(parents=True, exist_ok=True)
    inp = pdir / "preopt_input.xyz"
    shutil.copy(src, inp)
    method = {"gxtb": ["--gxtb"], "gfn2": ["--gfn", "2"], "gfnff": ["--gfnff"]}[args.preopt]
    cmd = [
        args.xtb_bin,
        "preopt_input.xyz",
        "--opt",
        "--chrg",
        str(charge),
        "--uhf",
        str(max(0, mult - 1)),
        *method,
    ]
    try:
        proc = subprocess.run(
            cmd, cwd=pdir, capture_output=True, text=True, timeout=args.preopt_timeout
        )
        stdout, stderr, note = (
            _as_text(proc.stdout),
            _as_text(proc.stderr),
            f"{args.preopt} rc={proc.returncode}",
        )
    except subprocess.TimeoutExpired as exc:
        stdout, stderr, note = (
            _as_text(exc.stdout),
            _as_text(exc.stderr),
            f"{args.preopt} timeout {args.preopt_timeout}s",
        )
    except FileNotFoundError:
        return None, f"{args.preopt} xtb-bin not found: {args.xtb_bin}"
    if log_path is not None:
        log_path.write_text(
            f"# cmd: {' '.join(cmd)}\n# result: {note}\n\n===== STDOUT =====\n{stdout}\n"
            f"\n===== STDERR =====\n{stderr}\n"
        )
    opt = pdir / "xtbopt.xyz"
    if opt.exists():
        return opt, note
    return None, note + " (no xtbopt.xyz; used original input)"


def run_crest(src: Path, workdir: Path, charge: int, mult: int, args, log_path=None) -> tuple:
    """Run CREST in workdir; return (conformers_path | None, note).

    Always writes CREST's full stdout/stderr to log_path (if given) — including on
    timeout — so failures are never silently discarded even without --keep-tmp.
    """
    inp = workdir / "input.xyz"
    shutil.copy(src, inp)
    cmd = ["crest", "input.xyz", "--chrg", str(charge), "--uhf", str(max(0, mult - 1))]
    if args.method:
        cmd.append(f"--{args.method}")
    if args.threads:
        cmd += ["-T", str(args.threads)]
    if args.crest_args:
        cmd += args.crest_args.split()
    try:
        proc = subprocess.run(
            cmd, cwd=workdir, capture_output=True, text=True, timeout=args.crest_timeout
        )
        stdout, stderr, note = _as_text(proc.stdout), _as_text(proc.stderr), f"rc={proc.returncode}"
    except subprocess.TimeoutExpired as exc:
        # Partial output captured up to the kill is on the exception.
        stdout, stderr, note = (
            _as_text(exc.stdout),
            _as_text(exc.stderr),
            f"timeout after {args.crest_timeout}s",
        )
    if log_path is not None:
        log_path.write_text(
            f"# cmd: {' '.join(cmd)}\n# result: {note}\n\n===== STDOUT =====\n{stdout}\n"
            f"\n===== STDERR =====\n{stderr}\n"
        )
    out = workdir / "crest_conformers.xyz"
    if not out.exists():
        found = list(workdir.rglob("crest_conformers.xyz"))
        out = found[0] if found else None
    if out is None:
        tail = (stderr or stdout or "").strip().splitlines()[-3:]
        return None, f"{note}, no ensemble ({' | '.join(tail)})"
    return out, note


def reopt_geom(elements, coords, charge, mult, args):
    """Re-optimize one conformer with xtb (--reopt); return (elements, coords, ok).

    Falls back to the input geometry (ok=False) on failure. Force fields distort metal
    coordination spheres (e.g. GFN-FF flattens/puckers square-planar d8), which desyncs
    the discrete OIN geometry label across conformers; an xtb (g-xTB) re-opt restores the
    physical geometry so the label is conformer-stable.
    """
    method = {"gxtb": ["--gxtb"], "gfn2": ["--gfn", "2"], "gfnff": ["--gfnff"]}[args.reopt]
    with tempfile.TemporaryDirectory() as td:
        write_xyz(Path(td) / "c.xyz", "conformer", elements, coords)
        cmd = [
            args.xtb_bin,
            "c.xyz",
            "--opt",
            "--chrg",
            str(charge),
            "--uhf",
            str(max(0, mult - 1)),
            *method,
        ]
        try:
            subprocess.run(cmd, cwd=td, capture_output=True, text=True, timeout=args.reopt_timeout)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return elements, coords, False
        opt = Path(td) / "xtbopt.xyz"
        if opt.exists():
            frames = read_frames(opt)
            if frames and len(frames[0][1]) == len(elements):
                return frames[0][1], frames[0][2], True
    return elements, coords, False


def analyze(structure, converter, conformers_path, args):
    """Encode each conformer; classify invariance vs chemistry vs canonicalization bug.

    With --reopt, each conformer is re-optimized with xtb (default g-xTB) before encoding.
    """
    charge, mult = structure.get("charge", 0), structure.get("mult", 1)
    frames = read_frames(conformers_path)[: args.max_confs]
    per = []
    n_reopt = 0
    for comment, elements, coords in frames:
        if args.reopt != "none":
            elements, coords, ok = reopt_geom(elements, coords, charge, mult, args)
            n_reopt += int(ok)
        with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as tf:
            tp = Path(tf.name)
        try:
            write_xyz(tp, comment, elements, coords)
            oin = converter.convert(str(tp))
        except Exception as exc:  # noqa: BLE001 - record, don't abort the whole run
            oin = f"<encode-error: {type(exc).__name__}: {exc}>"
        finally:
            tp.unlink(missing_ok=True)
        per.append(
            {
                "oin": oin,
                "key": str(winding_canonical_key(normalize_oin_for_comparison(oin)))
                if oin.startswith("[")
                else None,
                "conn": connectivity(elements, coords),
            }
        )

    raw = Counter(p["oin"] for p in per)
    keys = {p["key"] for p in per}
    # Group by connectivity; a canonicalization failure = one connectivity, >1 OIN.
    by_conn = defaultdict(set)
    for p in per:
        by_conn[p["conn"]].add(p["oin"])
    same_conn_diff_oin = any(len(oins) > 1 for oins in by_conn.values())

    if len(per) < 2:
        # CREST produced <2 conformers (rigid/small molecule): nothing to compare, so
        # "invariant" would be vacuous. Report honestly instead.
        verdict = "single-conformer"
    elif len(raw) == 1:
        verdict = "invariant"
    elif len(keys) == 1:
        verdict = "notation-drift"
    elif same_conn_diff_oin:
        # Same covalent connectivity, different OIN. Either a canonicalization bug
        # or a coordination/metal-stereo isomer the covalent signature can't see.
        # Flag for human review rather than auto-declaring a bug.
        verdict = "review-divergent"
    else:
        verdict = "connectivity-varied"

    return {
        "molecule": structure["molecule"],
        "metal": structure["metal"],
        "n_heavy": structure["n_heavy"],
        "n_conformers": len(per),
        "n_reopt": n_reopt,
        "n_distinct_oin": len(raw),
        "n_distinct_keys": len(keys),
        "n_distinct_connectivity": len(by_conn),
        "verdict": verdict,
        "distinct_oin": [{"oin": o, "count": c} for o, c in raw.most_common()]
        if len(raw) > 1
        else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--only", help="comma-separated molecule names to restrict to")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--max-confs", type=int, default=10, help="conformers encoded per structure")
    ap.add_argument(
        "--method", default="gfnff", help="CREST method flag (e.g. gfnff, gfn2); '' to omit"
    )
    ap.add_argument("--threads", type=int, default=None, help="CREST -T value")
    ap.add_argument("--crest-args", default="", help="extra args appended to the crest command")
    ap.add_argument(
        "--crest-timeout", type=int, default=900, help="per-structure wall-clock seconds"
    )
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument(
        "--preopt",
        choices=["none", "gxtb", "gfn2", "gfnff"],
        default="none",
        help="pre-optimize each input geometry with xtb before CREST (default none).",
    )
    ap.add_argument(
        "--xtb-bin",
        default=os.environ.get(
            "OIN_XTB_BIN", "/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/xtb"
        ),
        help="xtb binary for --preopt/--reopt (default: $OIN_XTB_BIN, else the g-xTB build "
        "in the main .venv). Must be a full path since the CREST env's plain xtb shadows it "
        "on PATH; shutil.which('xtb') would pick that GFN2 build over the g-xTB one.",
    )
    ap.add_argument(
        "--preopt-timeout", type=int, default=1200, help="per-structure pre-opt wall-clock seconds"
    )
    ap.add_argument(
        "--reopt",
        choices=["none", "gxtb", "gfn2", "gfnff"],
        default="none",
        help="re-optimize EACH CREST conformer with xtb before encoding (default none). "
        "Removes force-field metal-geometry distortions that desync the OIN geometry label.",
    )
    ap.add_argument(
        "--reopt-timeout", type=int, default=600, help="per-conformer re-opt wall-clock seconds"
    )
    args = ap.parse_args()

    if shutil.which("crest") is None:
        print(
            "CREST is not installed (no 'crest' binary on PATH). This is an optional\n"
            "external tool, not a dependency. See the README section\n"
            "'Optional: CREST conformer cross-check' to install it, then re-run.\n"
            "Nothing to do; exiting cleanly."
        )
        return 0

    if not args.manifest.exists():
        print(
            f"ERROR: manifest not found: {args.manifest}\n"
            "Generate the set with tools/select_conformer_test_set.py.",
            file=sys.stderr,
        )
        return 2

    structures = json.loads(args.manifest.read_text())["structures"]
    if args.only:
        want = {n.strip() for n in args.only.split(",")}
        structures = [s for s in structures if s["molecule"] in want]
    if args.limit:
        structures = structures[: args.limit]
    if not structures:
        print("No structures selected.", file=sys.stderr)
        return 2

    converter = XYZToSMILES()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for idx, s in enumerate(structures, 1):
        src = REPO / s["path"]
        print(
            f"[{idx}/{len(structures)}] {s['molecule']} ({s['metal']}, {s['n_heavy']} heavy) ...",
            flush=True,
        )
        base_tmp = args.out_dir / f"crest_{s['molecule']}" if args.keep_tmp else None
        ctx = (base_tmp, None) if args.keep_tmp else (None, tempfile.TemporaryDirectory())
        workdir = base_tmp or Path(ctx[1].name)
        if base_tmp:
            base_tmp.mkdir(parents=True, exist_ok=True)
        log_path = args.out_dir / f"{s['molecule']}.crest.log"
        charge, mult = s.get("charge", 0), s.get("mult", 1)
        crest_input, preopt_note = src, "no-preopt"
        if args.preopt != "none":
            opt, preopt_note = preopt(
                src,
                workdir / "preopt",
                charge,
                mult,
                args,
                log_path=args.out_dir / f"{s['molecule']}.preopt.log",
            )
            if opt is not None:
                crest_input = opt
            print(f"    preopt: {preopt_note}", flush=True)
        try:
            conformers, note = run_crest(
                crest_input, workdir, charge, mult, args, log_path=log_path
            )
            if conformers is None:
                res = {
                    "molecule": s["molecule"],
                    "metal": s["metal"],
                    "n_heavy": s["n_heavy"],
                    "verdict": "crest-failed",
                    "note": note,
                    "preopt": preopt_note,
                }
            else:
                res = analyze(s, converter, conformers, args)
                res["note"] = note
                res["preopt"] = preopt_note
        finally:
            if ctx[1] is not None:
                ctx[1].cleanup()
        results.append(res)
        extra = (
            f" [{res['n_distinct_oin']} distinct OIN / {res['n_conformers']} confs]"
            if "n_distinct_oin" in res
            else ""
        )
        print(f"    -> {res['verdict']}{extra}", flush=True)

    write_report(args.out_dir, results, args)
    tally = Counter(r["verdict"] for r in results)
    print("\n=== verdict tally ===")
    for v, c in tally.most_common():
        print(f"  {v}: {c}")
    # Exit non-zero when any structure needs review (same connectivity, different OIN).
    return 1 if tally.get("review-divergent") else 0


def write_report(out_dir: Path, results, args):
    (out_dir / "conformer_invariance_report.json").write_text(
        json.dumps(
            {
                "settings": {
                    "method": args.method,
                    "max_confs": args.max_confs,
                    "crest_timeout": args.crest_timeout,
                    "crest_args": args.crest_args,
                    "preopt": args.preopt,
                    "reopt": args.reopt,
                },
                "results": results,
            },
            indent=2,
        )
        + "\n"
    )
    tally = Counter(r["verdict"] for r in results)
    lines = [
        "# CREST conformer-invariance report",
        "",
        f"method `{args.method or 'default'}` · preopt `{args.preopt}` · "
        f"crest-args `{args.crest_args or 'none'}` · up to {args.max_confs} conformers/structure "
        f"· {args.crest_timeout}s timeout",
        "",
        "Verdicts: **invariant** (>=2 conformers, all → one OIN), **single-conformer** "
        "(CREST produced <2 conformers — rigid molecule, invariance not exercised), "
        "**notation-drift** (same canonical key, cosmetic string diff), "
        "**connectivity-varied** (CREST changed bonding — expected chemistry, not an "
        "encoder bug), **review-divergent** (same covalent connectivity but different OIN — "
        "either a canonicalization bug or a coordination/metal-stereo isomer the covalent "
        "signature cannot distinguish; review manually), **crest-failed** (no ensemble "
        "produced).",
        "",
        "## Tally",
        "",
        "| verdict | count |",
        "|---|---|",
    ]
    for v, c in tally.most_common():
        lines.append(f"| {v} | {c} |")
    lines += [
        "",
        "## Per structure",
        "",
        "| molecule | metal | heavy | confs | distinct OIN | distinct keys | verdict |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['molecule']} | {r['metal']} | {r['n_heavy']} | {r.get('n_conformers', '-')} | "
            f"{r.get('n_distinct_oin', '-')} | {r.get('n_distinct_keys', '-')} | {r['verdict']} |"
        )
    flagged = [
        r
        for r in results
        if r["verdict"] in ("review-divergent", "notation-drift") and r.get("distinct_oin")
    ]
    if flagged:
        lines += ["", "## Divergent strings", ""]
        for r in flagged:
            lines.append(f"### {r['molecule']} ({r['verdict']})")
            for d in r["distinct_oin"]:
                lines.append(f"- ×{d['count']}: `{d['oin']}`")
            lines.append("")
    (out_dir / "conformer_invariance_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
