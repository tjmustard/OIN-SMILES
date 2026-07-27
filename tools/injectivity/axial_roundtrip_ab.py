"""A/B the P2 axial fix end-to-end: does an axial-token OIN regenerate the right atropisomer?

This is the measurement that decides whether ``OIN_EMIT_AXIAL`` can be gated ON by default.
Emitting a stereo descriptor the generator cannot reproduce would convert a silent
round-trip FALSE POSITIVE into a generator-caused FALSE NEGATIVE, so the emit is only safe
once generation honours the token.

For each fixture and each enantiomer (the structure and its z-mirror):

1. encode with ``OIN_EMIT_AXIAL=1`` -> an OIN carrying ``|ax:+|`` or ``|ax:-|``
2. generate 3D from that OIN
3. read the generated structure's own canonical axial token
4. ``match`` iff it equals the requested one

``--no-select`` disables the adapter's axial-aware pass (by stripping the token from the
OIN just before generation), giving the baseline: how often the right atropisomer comes out
by chance. The delta between the two arms is the value of the fix.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.axial_roundtrip_ab [--no-select]
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
SEED = 42
FIXTURES = ["PdCl2-R-BINAP.xyz"]


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
    """Generate 3D from ``oin`` and return (generated_axial_token, error)."""
    from rdkit import Chem

    from oinsmiles.generation.engine import OIN3DGenerator

    try:
        with _silence():
            gen = OIN3DGenerator(optimizer="ff", seed=SEED, timeout=timeout)
            res = gen.generate(oin)
    except Exception as e:  # generation failure is a datum, not a crash
        return None, f"{type(e).__name__}: {e}"
    mol = getattr(res, "mol", None)
    if mol is None:
        xyz = getattr(res, "xyz", None)
        if not xyz:
            return None, "no mol and no xyz returned"
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "gen.xyz"
            p.write_text(xyz)
            from oinsmiles.utils.perception_tmc import get_tmc_mol

            try:
                with _silence():
                    mol, _ = get_tmc_mol(p, 0, with_stereo=False)
                Chem.SanitizeMol(mol)
            except Exception as e:
                return None, f"reperception failed: {type(e).__name__}: {e}"
    return mol_axial_token(mol), None


def run(no_select: bool, timeout: int) -> dict:
    from oinsmiles import XYZToSMILES

    rows = []
    os.environ["OIN_EMIT_AXIAL"] = "1"
    try:
        for fx in FIXTURES:
            src = REPO / "tests" / "fixtures" / fx
            with tempfile.TemporaryDirectory() as d:
                twins = {"base": src, "mirror": _mirror(src, Path(d) / "mirror.xyz")}
                for label, path in twins.items():
                    with _silence():
                        oin = XYZToSMILES().convert(str(path))
                    requested = parse_axial_token(oin)
                    gen_input = AXIAL_TOKEN_RE.sub("", oin).strip() if no_select else oin
                    got, err = _generated_token(gen_input, timeout)
                    rows.append(
                        {
                            "fixture": fx,
                            "twin": label,
                            "requested": requested,
                            "generated": got,
                            "match": (got == requested) if err is None else None,
                            "error": err,
                        }
                    )
                    print(
                        f"  {fx} [{label}]: requested={requested!r} generated={got!r} "
                        f"match={rows[-1]['match']}" + (f" err={err}" if err else ""),
                        flush=True,
                    )
    finally:
        os.environ.pop("OIN_EMIT_AXIAL", None)

    done = [r for r in rows if r["match"] is not None]
    return {
        "arm": "no-select (baseline)" if no_select else "axial-aware selection",
        "seed": SEED,
        "timeout_s": timeout,
        "rows": rows,
        "n_evaluated": len(done),
        "n_match": sum(1 for r in done if r["match"]),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-select", action="store_true", help="baseline arm: strip the token")
    ap.add_argument("--timeout", type=int, default=300)
    args = ap.parse_args(argv)

    print(f"== arm: {'no-select (baseline)' if args.no_select else 'axial-aware selection'} ==")
    res = run(args.no_select, args.timeout)
    OUT_DIR.mkdir(exist_ok=True)
    name = "axial_ab_baseline.json" if args.no_select else "axial_ab_selection.json"
    (OUT_DIR / name).write_text(json.dumps(res, indent=2) + "\n")
    print(f"\n{res['n_match']}/{res['n_evaluated']} matched  -> {OUT_DIR / name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
