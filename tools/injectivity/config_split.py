"""Split the unknown-unknown hunt's ambiguous residue into conformational vs configurational.

The Wave-3 hunt (:mod:`tools.injectivity.uu_hunt`) left an ambiguous residue: structures the
encoder collapses with its mirror, where **InChI agrees with the collapse** and no modelled
stereo element (CIP centre, P1 metal, P2 axial, P3 bound amine) explains the chirality the
*rigid* oracle reported. Two very different things live there:

* **conformational chirality** -- the mirror is a different conformer of the same isomer, so
  collapsing it is correct (that is the conformer invariance the encoder is supposed to have);
* **jointly-blind axes** -- InChI is itself blind to metal Delta/Lambda and atropisomerism, so a
  genuinely new axis that InChI also misses would hide in exactly this bucket.

:mod:`tools.injectivity.torsion_oracle` is the instrument that separates them. This module is
the driver: it reads the residue out of ``results-injectivity-y2/uu_hunt.json`` and applies the
oracle to every member, then reports the split with each verdict's supporting numbers.

**What a ``configurational`` row means, and what it does not.** It means: genuinely chiral,
the encoder emits one string for both twins, InChI cannot see the difference either, no
modelled stereo element accounts for it, and the mirror is not reachable by rotating bonds.
That is a *candidate* new blind spot and the acceptance bar is hand inspection of every one --
not a finding on its own. ``inconclusive`` rows are honest failures of the search, not weak
positives.

Run:
  PYTHONPATH=$PWD/src python -m tools.injectivity.config_split
  PYTHONPATH=$PWD/src python -m tools.injectivity.config_split --limit 5 --restarts 4
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem

from .config_oracle import bound_amine_centers, metal_stereo_descriptors
from .oracle import load_mol
from .torsion_oracle import N_RESTARTS, configurational_verdict, rotatable_torsions

REPO = Path(__file__).resolve().parents[2]
UU_JSON = REPO / "results-injectivity-y2" / "uu_hunt.json"
OUT_DIR = REPO / "results-injectivity-y3"


def load_residue(path: Path) -> list[dict]:
    """The FULL ambiguous residue. Falls back to the capped sample on an older json."""
    h = json.loads(path.read_text())
    rows = h.get("residue") or h.get("residue_examples") or []
    if not h.get("residue") and rows:
        print(
            f"  ! {path.name} predates the full-residue field; only {len(rows)} of "
            f"{h.get('n_residue_ambiguous')} residue rows are available. "
            "Re-run uu_hunt to split all of them.",
        )
    return rows


def _diagnostics(path: Path) -> dict:
    """Extra context a hand inspection needs, gathered once per structure."""
    try:
        mol, _coords = load_mol(path)
    except Exception as e:
        return {"diag_error": repr(e)[:150]}
    try:
        skel = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))
    except Exception:
        skel = ""
    root = max(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetAtomicNum())
    return {
        "skeleton_smiles": skel,
        "metal_descriptors": [
            {"shape": m.shape, "permutation": m.permutation} for m in metal_stereo_descriptors(mol)
        ],
        "bound_amines": [{"sign": a.sign, "volume": a.volume} for a in bound_amine_centers(mol)],
        "n_rings": mol.GetRingInfo().NumRings(),
        "torsion_bonds": [
            f"{mol.GetAtomWithIdx(a).GetSymbol()}{a}-{mol.GetAtomWithIdx(b).GetSymbol()}{b}"
            for a, b, _m, _d in rotatable_torsions(mol, root)
        ],
    }


def split(rows: list[dict], *, restarts: int, limit: int = 0) -> dict:
    out, counts = [], Counter()
    todo = rows[:limit] if limit else rows
    for i, r in enumerate(todo, 1):
        path = Path(r["path"])
        try:
            v = configurational_verdict(path, restarts=restarts, name=r["name"])
            row = v.to_dict()
        except Exception as e:
            row = {"name": r["name"], "verdict": "error", "note": repr(e)[:200]}
        row["path"] = str(path)
        row["uu_rigid_mirror_rmsd"] = r.get("mirror_rmsd")
        row["uu_rotatable_bonds"] = r.get("rotatable_bonds")
        row["oin"] = r.get("oin", "")
        if row["verdict"] == "configurational":
            row.update(_diagnostics(path))
        counts[row["verdict"]] += 1
        out.append(row)
        print(
            f"  [{i}/{len(todo)}] {row['name']:22s} {row['verdict']:16s} "
            f"d_mirror={row.get('d_mirror')} d_control={row.get('d_control')}",
            flush=True,
        )
    # the ones that must be inspected by hand, hardest evidence first
    survivors = sorted(
        (r for r in out if r["verdict"] == "configurational"),
        key=lambda r: -float(r.get("d_mirror") or 0.0),
    )
    return {
        "n_residue": len(todo),
        "restarts": restarts,
        "counts": dict(counts),
        "survivors_needing_hand_inspection": survivors,
        "rows": out,
    }


def render(s: dict) -> str:
    n = max(s["n_residue"], 1)
    c = s["counts"]

    def line(k, label, meaning):
        return f"| {label} | {c.get(k, 0)} | {c.get(k, 0) / n:.1%} | {meaning} |"

    lines = [
        "# Splitting the Wave-3 ambiguous residue (Lane 7 / Task A)",
        "",
        "The unknown-unknown hunt left an ambiguous residue: mirror twins the encoder collapses,",
        "where InChI agrees with the collapse and no modelled stereo element explains the",
        "chirality the *rigid* oracle reported. The torsion-aware oracle",
        "(`tools/injectivity/torsion_oracle.py`) asks the question the rigid one cannot: **is the",
        "mirror reachable from the structure by rotating bonds?** If it is, the two are conformers",
        "of one isomer and collapsing them is correct. If it is not, the difference is",
        "configurational.",
        "",
        f"- residue members probed: **{s['n_residue']}**",
        f"- search budget: {s['restarts']} random restarts plus the dihedral-negating seed,",
        "  4 sweeps, 30 deg coarse grid; match threshold 0.5 A",
        "- every verdict is paired with a positive control on the same molecule and budget: a",
        "  randomly-torsioned copy of the structure, reachable **by construction**. A failure to",
        "  reach the mirror only counts when the control succeeded.",
        "",
        "| verdict | count | share | meaning |",
        "|---|---:|---:|---|",
        line(
            "conformational", "conformational", "mirror is a conformer -- **collapse is CORRECT**"
        ),
        line(
            "rigid_achiral",
            "rigid achiral",
            "mirror superimposes with no torsion change -- correct invariance; the rigid "
            "oracle over-reported",
        ),
        line(
            "configurational",
            "**configurational**",
            "**candidate new blind spot -- must be inspected by hand**",
        ),
        line("inconclusive", "inconclusive", "the search failed its own control -- no evidence"),
        line("error", "error", "structure could not be probed"),
        "",
    ]
    surv = s["survivors_needing_hand_inspection"]
    if surv:
        lines += [
            "## Survivors requiring hand inspection",
            "",
            "| structure | d_mirror | d_control | rigid RMSD | torsions | worst residual atom |",
            "|---|---:|---:|---:|---:|---|",
        ]
        for r in surv:
            worst = r.get("worst_atoms") or [{}]
            w = worst[0]
            lines.append(
                f"| {r['name']} | {r.get('d_mirror')} | {r.get('d_control')} "
                f"| {r.get('rigid_mirror_rmsd')} | {r.get('n_torsions')} "
                f"| {w.get('element', '-')} {w.get('deviation', '-')} A |"
            )
        lines += [
            "",
            "`d_mirror` is the best mirror RMSD over the whole torsion orbit; `d_control` is the",
            "same optimiser recovering a reachable target. A survivor with a large `d_mirror` and",
            "a tiny `d_control` is the strong form of the claim.",
            "",
        ]
    else:
        lines += [
            "## No survivors",
            "",
            "Every member of the residue is either a conformer of its own mirror or was shown to",
            "be achiral once the automorphism set was complete. On this sample the residue holds",
            "**no unnamed configurational axis** -- a negative result for the sample size, not a",
            "proof of injectivity.",
            "",
        ]
    lines += [
        "## Reproduce",
        "",
        "```",
        "PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 300   # writes the residue",
        "PYTHONPATH=$PWD/src python -m tools.injectivity.config_split",
        "```",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uu-json", type=Path, default=UU_JSON)
    ap.add_argument("--restarts", type=int, default=N_RESTARTS)
    ap.add_argument("--limit", type=int, default=0, help="probe only the first N (smoke test)")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)

    if not args.uu_json.exists():
        raise SystemExit(
            f"{args.uu_json} not found -- run: python -m tools.injectivity.uu_hunt --n 300"
        )
    rows = load_residue(args.uu_json)
    if not rows:
        raise SystemExit("no residue rows in the uu_hunt json")
    s = split(rows, restarts=args.restarts, limit=args.limit)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "config_split.json").write_text(json.dumps(s, indent=2) + "\n")
    md = render(s)
    (args.out / "config_split.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
