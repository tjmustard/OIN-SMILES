"""Population-at-risk for the P2 axial blind spot (Y2) -- a SOUND, conformation-independent count.

The Wave-1 rigid-mirror chirality scan is conformation-inflated and can only give an upper
bound (``report.py`` labels it as such). Motif counting does not have that problem: whether a
structure *contains* a hindered, stereogenic biaryl axis is a property of the graph plus the
torsion, not of which conformer was deposited. This is the same reasoning that made the P3
metal-bound-amine figure (~2.85 % of corpus) trustworthy.

Reports, over a seeded sample of the dataset:

* ``any_biaryl``    -- has an inter-ring aromatic single bond at all
* ``hindered``      -- that bond is twisted and ortho-walled (atropisomer *candidate*)
* ``emitting``      -- hindered AND stereogenic -> ``OIN_EMIT_AXIAL`` would emit a token.
  **This is the blast radius**: the set of structures whose OIN string changes when the flag
  is turned on, and the set currently suffering the P2 encoder blind spot.
* ``hindered_not_stereogenic`` -- the over-sensitivity that the stereogenicity gate prevents:
  sterically twisted but achiral (a ring end with symmetry-equivalent ortho neighbours).

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.axial_population --n 1500 [--dataset DIR]
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import random
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.oin.axial import detect_axial_axes  # noqa: E402

DEFAULT_DATASET = Path("/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset")
OUT_DIR = REPO / "results-injectivity-y2"
SAMPLE_SEED = 42


@contextlib.contextmanager
def _silence():
    """Mute xyz2mol's C-level 'suspicious distance' chatter (crystal C-H lengths)."""
    with open(os.devnull, "w") as devnull:
        old = os.dup(2)
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            os.dup2(old, 2)
            os.close(old)


def _iter_dataset(dataset: Path) -> list[Path]:
    files: list[Path] = []
    for sub in ("cat", "photo"):
        d = dataset / sub
        if d.is_dir():
            files.extend(sorted(d.glob("*.xyz")))
    if not files:
        files = sorted(dataset.rglob("*.xyz"))
    return files


def _mirror_flip_check(path: Path) -> dict:
    """Does this structure's axial token FLIP for its mirror image?

    The sign convention is only trustworthy if it tracks *handedness*. If the reference
    ortho neighbour instead correlated with something structural, signs could be
    systematically skewed and the token would be quietly wrong at scale -- the emitted
    signs do run ~3:1 toward ``+`` across the corpus, which has a benign chemical reading
    ((R)-BINAP-like ligands dominate catalysis sets) and a malign one. Mirroring each
    structure separates them: a sound convention flips EVERY time, whatever the skew.
    """
    from oinsmiles.oin.axial import axial_token
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    lines = path.read_text().splitlines()
    n_at = int(lines[0])
    body = []
    for ln in lines[2 : 2 + n_at]:
        p = ln.split()
        body.append(f"{p[0]}  {p[1]}  {p[2]}  {-float(p[3]):.6f}")
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        mp = Path(d) / "mirror.xyz"
        mp.write_text(f"{lines[0]}\n{lines[1]}\n" + "\n".join(body) + "\n")
        with _silence():
            base_mol, _ = get_tmc_mol(path, 0, with_stereo=False)
            mir_mol, _ = get_tmc_mol(mp, 0, with_stereo=False)
    base_tok = axial_token(base_mol)
    mir_tok = axial_token(mir_mol)
    expected = base_tok.translate(str.maketrans("+-", "-+"))
    return {
        "base": base_tok,
        "mirror": mir_tok,
        "expected_mirror": expected,
        "flips": bool(base_tok) and mir_tok == expected,
    }


def scan(dataset: Path, n: int, mirror_check: bool = False) -> dict:
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    files = _iter_dataset(dataset)
    rng = random.Random(SAMPLE_SEED)
    sample = files if n <= 0 or n >= len(files) else rng.sample(files, n)
    sample = sorted(sample)

    counts = Counter()
    emitting_examples: list[dict] = []
    oversensitive_examples: list[dict] = []
    failures = 0

    for i, path in enumerate(sample, 1):
        try:
            with _silence():
                mol, _ = get_tmc_mol(path, 0, with_stereo=False)
            axes = detect_axial_axes(mol)
        except Exception:
            failures += 1
            continue
        counts["scanned"] += 1
        if not axes:
            continue
        counts["any_biaryl"] += 1
        hindered = [a for a in axes if a.hindered]
        if not hindered:
            continue
        counts["hindered"] += 1
        emitting = [a for a in hindered if a.stereogenic]
        if emitting:
            counts["emitting"] += 1
            row = {
                "name": path.stem,
                "path": str(path),
                "n_axes": len(emitting),
                "dihedrals": [a.dihedral_deg for a in emitting],
                "token": "".join("+" if a.sign > 0 else "-" for a in emitting),
            }
            if mirror_check:
                try:
                    row["mirror"] = _mirror_flip_check(path)
                except Exception as e:
                    row["mirror"] = {"error": f"{type(e).__name__}: {e}", "flips": None}
                if row["mirror"].get("flips") is True:
                    counts["mirror_flips"] += 1
                elif row["mirror"].get("flips") is False:
                    counts["mirror_NO_flip"] += 1
                else:
                    counts["mirror_error"] += 1
            emitting_examples.append(row)
        else:
            counts["hindered_not_stereogenic"] += 1
            if len(oversensitive_examples) < 25:
                oversensitive_examples.append(
                    {"name": path.stem, "dihedrals": [a.dihedral_deg for a in hindered]}
                )
        if i % 100 == 0:
            print(f"  ... {i}/{len(sample)}  emitting={counts['emitting']}", flush=True)

    scanned = max(counts["scanned"], 1)
    return {
        "dataset": str(dataset),
        "sample_seed": SAMPLE_SEED,
        "requested": n,
        "counts": dict(counts),
        "load_failures": failures,
        "fractions": {
            k: round(counts[k] / scanned, 4)
            for k in ("any_biaryl", "hindered", "emitting", "hindered_not_stereogenic")
        },
        "emitting_examples": emitting_examples,
        "oversensitive_examples": oversensitive_examples,
    }


def render(res: dict) -> str:
    c, f = res["counts"], res["fractions"]
    scanned = c.get("scanned", 0)
    lines = [
        "# P2 axial population-at-risk (Y2)",
        "",
        "Conformation-independent motif count -- a SOUND rate, unlike the Wave-1 rigid-mirror",
        "chirality scan (an upper bound only). Seeded sample; re-running reproduces it.",
        "",
        f"- dataset: `{res['dataset']}`",
        f"- sampled / scanned: {res['requested']} / **{scanned}**"
        f" (load failures: {res['load_failures']})",
        "",
        "| population | count | fraction of scanned |",
        "|---|---:|---:|",
        f"| has any inter-ring aromatic single bond | {c.get('any_biaryl', 0)} | {f['any_biaryl']:.2%} |",
        f"| hindered (twisted + ortho-walled) | {c.get('hindered', 0)} | {f['hindered']:.2%} |",
        f"| **emitting (hindered + stereogenic)** | **{c.get('emitting', 0)}** | **{f['emitting']:.2%}** |",
        f"| hindered but NOT stereogenic (gate suppresses) | {c.get('hindered_not_stereogenic', 0)} "
        f"| {f['hindered_not_stereogenic']:.2%} |",
        "",
        "**Blast radius** = the *emitting* row: the structures whose OIN string changes when",
        "`OIN_EMIT_AXIAL` is on, and equally the structures currently exposed to the P2 blind",
        "spot (their enantiomers encode byte-identically today).",
        "",
        "The last row is the over-sensitivity the stereogenicity gate prevents: sterically",
        "twisted biaryls that are nonetheless achiral (a ring end with symmetry-equivalent",
        "ortho neighbours). Without the gate the encoder would claim chirality that isn't there.",
        "",
    ]
    if "mirror_flips" in c or "mirror_NO_flip" in c:
        flips = c.get("mirror_flips", 0)
        noflip = c.get("mirror_NO_flip", 0)
        err = c.get("mirror_error", 0)
        total = flips + noflip + err
        verdict = "SOUND" if noflip == 0 and flips else "**BROKEN**"
        lines += [
            "## Sign-convention audit (does the token flip for the mirror?)",
            "",
            "The emitted signs skew toward `+`. That has a benign reading (catalysis sets are",
            "rich in single-enantiomer ligands) and a malign one (the reference-neighbour rule",
            "tracks structure rather than handedness). Mirroring every emitting structure",
            "separates them: a sound convention flips **every** time, whatever the skew.",
            "",
            f"- flips correctly: **{flips}/{total}**",
            f"- does NOT flip: **{noflip}**",
            f"- errored: {err}",
            "",
            f"**Verdict: the sign convention is {verdict}.**",
            "",
        ]
    if res["emitting_examples"]:
        lines += [
            "## Emitting examples",
            "",
            "| structure | axes | dihedrals | token |",
            "|---|---:|---|---|",
        ]
        for e in res["emitting_examples"]:
            lines.append(
                f"| {e['name']} | {e['n_axes']} | {', '.join(f'{d:+.1f}' for d in e['dihedrals'])} "
                f"| `{e['token']}` |"
            )
        lines.append("")
    if res["oversensitive_examples"]:
        lines += [
            "## Suppressed by the stereogenicity gate (achiral)",
            "",
            "| structure | dihedrals |",
            "|---|---|",
        ]
        for e in res["oversensitive_examples"]:
            lines.append(f"| {e['name']} | {', '.join(f'{d:+.1f}' for d in e['dihedrals'])} |")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=1500, help="sample size (<=0 for the whole set)")
    ap.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    ap.add_argument(
        "--mirror-check",
        action="store_true",
        help="verify each emitting structure's token FLIPS for its mirror (sign-convention audit)",
    )
    args = ap.parse_args(argv)

    if not args.dataset.exists():
        print(f"dataset not found: {args.dataset}", file=sys.stderr)
        return 2
    res = scan(args.dataset, args.n, mirror_check=args.mirror_check)
    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / "axial_population.json").write_text(json.dumps(res, indent=2) + "\n")
    md = render(res)
    (OUT_DIR / "axial_population.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
