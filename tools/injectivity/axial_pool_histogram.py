"""Where does a multi-axis atropisomer request die? Pool-generation, or selection?

The cohort A/B says multi-axis structures round-trip 0/2 in BOTH arms, and the generated
structure returns an EMPTY axial token. Empty admits three very different causes, with
opposite fixes:

* **pool-generation** -- no embedded conformer ever holds two hindered axes at once, so
  selection has nothing to pick (remedy: widen the pool / constrain the embed);
* **selection** -- some pooled conformer *does* carry the requested token but a later stage
  returns a different one (remedy: fix the filter / its ordering);
* **relaxation** -- a pooled conformer carried it and the FF then flattened the torsion out
  of the hindered window (remedy: guard the relaxation).

This tool separates them by dumping, for every conformer in the pool the adapter actually
sees, the full per-axis detail behind its token: the signed dihedral, whether it passed the
``hindered`` window (20-160 deg) and whether it is ``stereogenic``. A structure whose axes
are all detected but not ``hindered`` has *relaxed flat*; one whose axes are missing
entirely has a graph/perception problem. Aggregated as a histogram over the pool.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.axial_pool_histogram \
          tests/fixtures/YESKOZ.xyz [--timeout 300] [--pool N]
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.oin.axial import detect_axial_axes, mol_axial_token, parse_axial_token  # noqa: E402

OUT_DIR = REPO / "results-injectivity-y2"
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


def _axis_detail(mol) -> list[dict]:
    """Per-axis diagnostics for one conformer, or ``[]`` when nothing is detectable."""
    try:
        axes = detect_axial_axes(mol)
    except Exception as e:
        return [{"error": f"{type(e).__name__}: {e}"}]
    return [
        {
            "a1": ax.a1,
            "a2": ax.a2,
            "dihedral": ax.dihedral_deg,
            "hindered": ax.hindered,
            "stereogenic": ax.stereogenic,
            "sign": ax.sign,
        }
        for ax in axes
    ]


def probe(xyz: Path, timeout: int, pool: int | None) -> dict:
    """Encode ``xyz`` with a token, generate, and report the pool the adapter saw."""
    from oinsmiles import XYZToSMILES
    from oinsmiles.generation import metallogen_adapter as ma
    from oinsmiles.generation.engine import OIN3DGenerator

    os.environ["OIN_EMIT_AXIAL"] = "1"
    try:
        with _silence():
            oin = XYZToSMILES().convert(str(xyz))
    finally:
        os.environ.pop("OIN_EMIT_AXIAL", None)
    requested = parse_axial_token(oin)

    captured: dict = {}
    real_select = ma._select_by_geometry

    def spy(parsed, mols, honor_winding=True, early_exit=False):
        # Snapshot the pool BEFORE any narrowing so we see what the embed produced,
        # not what selection kept.
        rows = []
        for i, m in enumerate(mols):
            cmol = ma.build_contract_mol(parsed, m)
            rows.append(
                {
                    "rank": i,
                    "contract_mol": cmol is not None,
                    "token": mol_axial_token(cmol) if cmol is not None else None,
                    "axes": _axis_detail(cmol) if cmol is not None else [],
                }
            )
        captured.setdefault("pool", []).extend(rows)
        captured["n_pool"] = len(mols)
        return real_select(parsed, mols, honor_winding=honor_winding, early_exit=early_exit)

    ma._select_by_geometry = spy
    ff_params = {"uff_pool_size": 2 * pool} if pool else None
    try:
        with _silence():
            res = OIN3DGenerator(
                optimizer="ff",
                seed=SEED,
                timeout=timeout,
                ensemble_size=pool or 1,
                ff_params=ff_params,
            ).generate(oin)
        final_token = mol_axial_token(res.mol) if getattr(res, "mol", None) is not None else None
        final_axes = _axis_detail(res.mol) if getattr(res, "mol", None) is not None else []
        error = None
    except Exception as e:
        final_token, final_axes, error = None, [], f"{type(e).__name__}: {e}"
    finally:
        ma._select_by_geometry = real_select

    pool_rows = captured.get("pool", [])
    hist = Counter(r["token"] if r["token"] is not None else "<unperceived>" for r in pool_rows)
    # Was the requested token EVER in the pool? That single bit decides the remedy.
    ever = any(r["token"] == requested for r in pool_rows)
    return {
        "structure": xyz.stem,
        "oin": oin,
        "requested": requested,
        "n_pool_seen": len(pool_rows),
        "pool_token_histogram": dict(hist),
        "requested_in_pool": ever,
        "final_token": final_token,
        "final_axes": final_axes,
        "error": error,
        "pool": pool_rows,
        "verdict": (
            "SELECTION/RELAXATION -- the pool contained the requested token"
            if ever and final_token != requested
            else "OK -- requested token returned"
            if final_token == requested
            else "POOL-GENERATION -- no embedded conformer ever held the requested token"
        ),
    }


def render(r: dict) -> str:
    lines = [
        f"### {r['structure']}  requested `{r['requested']}`",
        f"- pool conformers seen by the adapter: **{r['n_pool_seen']}**",
        f"- requested token present anywhere in the pool: **{r['requested_in_pool']}**",
        f"- final returned token: `{r['final_token']}`",
        f"- verdict: **{r['verdict']}**",
        "",
        "| pool token | count |",
        "|---|---:|",
    ]
    for tok, n in sorted(r["pool_token_histogram"].items(), key=lambda kv: -kv[1]):
        lines.append(f"| `{tok}` | {n} |")
    lines.append("")
    if r["pool"]:
        lines += [
            "Per-axis detail (why a token is empty: `hindered=False` means the torsion "
            "relaxed out of the 20-160 deg window):",
            "",
            "| rank | token | axes (dihedral / hindered / stereogenic) |",
            "|---:|---|---|",
        ]
        for row in r["pool"]:
            axes = (
                ", ".join(
                    f"{a.get('dihedral')}/{'H' if a.get('hindered') else 'h'}"
                    f"/{'S' if a.get('stereogenic') else 's'}"
                    for a in row["axes"]
                )
                or "(none detected)"
            )
            lines.append(f"| {row['rank']} | `{row['token']}` | {axes} |")
        lines.append("")
    if r["error"]:
        lines += [f"- generation error: `{r['error']}`", ""]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("xyz", nargs="+", type=Path)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--pool", type=int, default=None, help="force an ensemble/pool width")
    ap.add_argument("--out", default="axial_pool_histogram")
    args = ap.parse_args(argv)

    results = []
    md = ["# Axial pool token histogram", ""]
    for x in args.xyz:
        if not x.exists():
            print(f"missing: {x}", file=sys.stderr)
            return 2
        print(f"== probing {x.stem} ==", flush=True)
        r = probe(x, args.timeout, args.pool)
        results.append(r)
        block = render(r)
        md.append(block)
        print(block, flush=True)

    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / f"{args.out}.json").write_text(json.dumps(results, indent=2) + "\n")
    (OUT_DIR / f"{args.out}.md").write_text("\n".join(md))
    print(f"-> {OUT_DIR / (args.out + '.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
