"""Cohort A/B for the P2 axial fix -- turns the 2/2 fixture result into a RATE.

``axial_roundtrip_ab`` measured one fixture pair. That is directional, not a rate, and it
never exercised a multi-axis token (the corpus has structures with 2 and 4 hindered
stereogenic axes, e.g. ``-++-``). This runs the same experiment over the emitting cohort
found by ``axial_population --mirror-check``.

For each structure, and for each of its two enantiomers (deposited and z-mirrored):

1. encode with ``OIN_EMIT_AXIAL=1`` -> an OIN carrying an axial token
2. generate 3D from that OIN
3. read the generated structure's own canonical axial token
4. ``match`` iff it equals the requested one

Two arms. ``--no-select`` strips the token just before generation, disabling the adapter's
axial-aware pass, and measures how often the right atropisomer appears by chance. The delta
is what the fix buys; the selection arm's rate is what gates the ON-by-default decision.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.axial_cohort_ab [--no-select] [--limit N]
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.oin.axial import AXIAL_TOKEN_RE, mol_axial_token, parse_axial_token  # noqa: E402

OUT_DIR = REPO / "results-injectivity-y2"
COHORT_JSON = OUT_DIR / "axial_population.json"
SEED = 42


@contextlib.contextmanager
def _silence():
    with open(os.devnull, "w") as devnull:
        old = os.dup(2)
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            os.dup2(old, 2)
            os.close(old)


def load_cohort(limit: int = 0) -> list[dict]:
    if not COHORT_JSON.exists():
        raise SystemExit(
            f"cohort not found: {COHORT_JSON}\n"
            "run: python -m tools.injectivity.axial_population --n 1500 --mirror-check"
        )
    rows = json.loads(COHORT_JSON.read_text()).get("emitting_examples", [])
    rows = [r for r in rows if r.get("path") and Path(r["path"]).exists()]
    rows.sort(key=lambda r: r["name"])
    return rows[:limit] if limit else rows


def _mirror(path: Path, dst: Path) -> Path:
    lines = path.read_text().splitlines()
    n = int(lines[0])
    out = [lines[0], lines[1]]
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        out.append(f"{p[0]}  {p[1]}  {p[2]}  {-float(p[3]):.6f}")
    dst.write_text("\n".join(out) + "\n")
    return dst


def _generated_token(oin: str, timeout: int):
    from oinsmiles.generation.engine import OIN3DGenerator

    try:
        with _silence():
            res = OIN3DGenerator(optimizer="ff", seed=SEED, timeout=timeout).generate(oin)
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"
    mol = getattr(res, "mol", None)
    if mol is None:
        return None, "generator returned no mol"
    return mol_axial_token(mol), None


def run(cohort: list[dict], no_select: bool, timeout: int) -> dict:
    from oinsmiles import XYZToSMILES

    rows = []
    os.environ["OIN_EMIT_AXIAL"] = "1"
    try:
        for i, entry in enumerate(cohort, 1):
            src = Path(entry["path"])
            with tempfile.TemporaryDirectory() as d:
                twins = {"base": src, "mirror": _mirror(src, Path(d) / "mirror.xyz")}
                for label, path in twins.items():
                    try:
                        with _silence():
                            oin = XYZToSMILES().convert(str(path))
                    except Exception as e:
                        rows.append(
                            {
                                "name": entry["name"],
                                "twin": label,
                                "error": f"encode: {type(e).__name__}: {e}",
                                "match": None,
                            }
                        )
                        continue
                    requested = parse_axial_token(oin)
                    gen_in = AXIAL_TOKEN_RE.sub("", oin).strip() if no_select else oin
                    got, err = _generated_token(gen_in, timeout)
                    rows.append(
                        {
                            "name": entry["name"],
                            "twin": label,
                            "n_axes": entry.get("n_axes"),
                            "requested": requested,
                            "generated": got,
                            "match": (got == requested) if err is None else None,
                            "error": err,
                        }
                    )
            done = [r for r in rows if r["match"] is not None]
            hit = sum(1 for r in done if r["match"])
            print(
                f"  [{i}/{len(cohort)}] {entry['name']}: running {hit}/{len(done)}",
                flush=True,
            )
    finally:
        os.environ.pop("OIN_EMIT_AXIAL", None)

    done = [r for r in rows if r["match"] is not None]
    multi = [r for r in done if (r.get("n_axes") or 1) > 1]
    return {
        "arm": "no-select (baseline)" if no_select else "axial-aware selection",
        "seed": SEED,
        "timeout_s": timeout,
        "n_structures": len(cohort),
        "n_evaluated": len(done),
        "n_match": sum(1 for r in done if r["match"]),
        "n_errors": sum(1 for r in rows if r["match"] is None),
        "multi_axis_evaluated": len(multi),
        "multi_axis_match": sum(1 for r in multi if r["match"]),
        "rows": rows,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-select", action="store_true", help="baseline arm: strip the token")
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--limit", type=int, default=0, help="cap the cohort (0 = all)")
    args = ap.parse_args(argv)

    cohort = load_cohort(args.limit)
    arm = "no-select (baseline)" if args.no_select else "axial-aware selection"
    print(f"== arm: {arm} == cohort: {len(cohort)} structures ({2 * len(cohort)} enantiomers)")
    res = run(cohort, args.no_select, args.timeout)
    OUT_DIR.mkdir(exist_ok=True)
    name = "axial_cohort_baseline.json" if args.no_select else "axial_cohort_selection.json"
    (OUT_DIR / name).write_text(json.dumps(res, indent=2) + "\n")
    pct = 100.0 * res["n_match"] / res["n_evaluated"] if res["n_evaluated"] else 0.0
    print(
        f"\n{arm}: {res['n_match']}/{res['n_evaluated']} matched ({pct:.1f}%), "
        f"{res['n_errors']} unevaluable; "
        f"multi-axis {res['multi_axis_match']}/{res['multi_axis_evaluated']}"
    )
    print(f"-> {OUT_DIR / name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
