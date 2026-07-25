"""Wave 3 -- the unknown-unknown hunt: find chirality axes nobody has named.

Waves 1-2 confirmed three *named* blind spots (metal Delta/Lambda, axial, metal-bound amine).
This looks for the ones no one has thought of, using the same generator-free mirror-twin
method but driven by a triage that isolates the interesting residue.

For each structure we mirror it, encode both twins, and ask three independent questions:

* is it chiral at all?           -- geometric oracle (Wave 1)
* does the encoder separate them? -- ``raw_equal`` from ``convert()``
* if not, is the cause KNOWN?     -- RDKit CIP stereocentres, or one of the P1/P2/P3 axes
  from the configurational oracle (Wave 2)

A collapse whose cause is known is a re-detection of an existing finding. A collapse with
**no known cause** is the signal: a chiral structure the encoder cannot separate, and whose
chirality is not explained by any stereo element we currently model. That is an unnamed axis.

**Known weakness, stated up front.** The geometric oracle is a rigid-superposition test, so
on flexible structures it over-reports chirality: a floppy achiral molecule's mirror is a
non-superimposable *conformer* and reads as chiral. Wave 1 hit this and refused to publish a
rate. The same caution applies here -- the candidate list is an UPPER BOUND and a triage
queue, not a set of confirmed findings. To keep the noise down we report a rigidity proxy
(rotatable-bond count) with every candidate and sort the least flexible first, because those
are the ones where the oracle is trustworthy.

Run:  PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 400
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import random
import sys
import tempfile
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
for _p in (str(REPO), str(REPO / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

from oinsmiles.oin.axial import detect_axial_axes  # noqa: E402
from tools.injectivity.config_oracle import (  # noqa: E402
    bound_amine_centers,
    metal_stereo_descriptors,
)
from tools.injectivity.oracle import geometric_chirality  # noqa: E402

DEFAULT_DATASET = Path("/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset")
OUT_DIR = REPO / "results-injectivity-y2"
SAMPLE_SEED = 42

#: at or below this many rotatable bonds the rigid mirror test is trustworthy enough that a
#: "chiral" verdict is unlikely to be pure conformational noise.
RIGID_ROTB_MAX = 4

VERDICTS = (
    "achiral_ok",  # not chiral; encoder collapses -- correct invariance
    "distinguished",  # chiral and the encoder separates the twins -- correct
    "known_cip",  # collapse, but RDKit sees stereocentres (a modelled stereo type)
    "known_axis",  # collapse explained by a named P1/P2/P3 axis
    "CONFIRMED_BLIND",  # collapse while InChI SEPARATES the twins -- rigorous, no confound
    "conformational_or_joint_blind",  # collapse, InChI agrees -- ambiguous residue
)


def _inchi(mol):
    """Configuration-based, conformation-INDEPENDENT identity. None when unavailable."""
    from rdkit.Chem import inchi

    try:
        probe = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(probe)
        return inchi.MolToInchi(probe) or None
    except Exception:
        return None


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


def _iter_dataset(dataset: Path) -> list[Path]:
    files: list[Path] = []
    for sub in ("cat", "photo"):
        d = dataset / sub
        if d.is_dir():
            files.extend(sorted(d.glob("*.xyz")))
    return files or sorted(dataset.rglob("*.xyz"))


def _mirror_file(path: Path, dst: Path) -> Path:
    lines = path.read_text().splitlines()
    n = int(lines[0])
    out = [lines[0], lines[1]]
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        out.append(f"{p[0]}  {p[1]}  {p[2]}  {-float(p[3]):.6f}")
    dst.write_text("\n".join(out) + "\n")
    return dst


def _known_causes(mol: Chem.Mol) -> list[str]:
    """Named stereo elements that would explain a chirality the encoder lost."""
    causes = []
    try:
        probe = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(probe)
        if any(a.HasProp("_CIPCode") for a in probe.GetAtoms()):
            causes.append("cip_stereocentre")
    except Exception:
        pass
    try:
        if metal_stereo_descriptors(mol):
            causes.append("P1_metal")
    except Exception:
        pass
    try:
        if any(a.emits for a in detect_axial_axes(mol)):
            causes.append("P2_axial")
    except Exception:
        pass
    try:
        if bound_amine_centers(mol):
            causes.append("P3_bound_amine")
    except Exception:
        pass
    return causes


def probe(path: Path) -> dict | None:
    from oinsmiles import XYZToSMILES
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    with tempfile.TemporaryDirectory() as d:
        mirror = _mirror_file(path, Path(d) / "mirror.xyz")
        try:
            with _silence():
                mol, _ = get_tmc_mol(path, 0, with_stereo=False)
                Chem.SanitizeMol(mol)
                mirror_mol, _ = get_tmc_mol(mirror, 0, with_stereo=False)
                Chem.SanitizeMol(mirror_mol)
                oin_a = XYZToSMILES().convert(str(path))
                oin_b = XYZToSMILES().convert(str(mirror))
        except Exception:
            return None

    try:
        rmsd, _n_auto, chiral = geometric_chirality(mol, mol.GetConformer().GetPositions())
    except Exception:
        return None

    raw_equal = oin_a == oin_b
    causes = _known_causes(mol) if (chiral and raw_equal) else []

    # InChI is configuration-based and conformation-INDEPENDENT. If it separates the twins
    # while the OIN does not, the encoder has lost something standard cheminformatics keeps
    # -- a blind spot with NO conformational confound. That is the rigorous signal.
    inchi_a = inchi_b = None
    inchi_separates = None
    if raw_equal and chiral:
        with _silence():
            inchi_a, inchi_b = _inchi(mol), _inchi(mirror_mol)
        if inchi_a and inchi_b:
            inchi_separates = inchi_a != inchi_b

    if not chiral:
        verdict = "achiral_ok" if raw_equal else "distinguished"
    elif not raw_equal:
        verdict = "distinguished"
    elif inchi_separates:
        verdict = "CONFIRMED_BLIND"
    elif "cip_stereocentre" in causes:
        verdict = "known_cip"
    elif causes:
        verdict = "known_axis"
    else:
        verdict = "conformational_or_joint_blind"

    try:
        rot_b = Descriptors.NumRotatableBonds(mol)
    except Exception:
        rot_b = -1

    return {
        "name": path.stem,
        "path": str(path),
        "verdict": verdict,
        "mirror_rmsd": round(float(rmsd), 3),
        "raw_equal": raw_equal,
        "known_causes": causes,
        "inchi_separates": inchi_separates,
        "rotatable_bonds": rot_b,
        "rigid": 0 <= rot_b <= RIGID_ROTB_MAX,
        "oin": oin_a,
    }


def hunt(dataset: Path, n: int) -> dict:
    files = _iter_dataset(dataset)
    rng = random.Random(SAMPLE_SEED)
    sample = sorted(files if n <= 0 or n >= len(files) else rng.sample(files, n))

    rows, counts, failures = [], Counter(), 0
    for i, path in enumerate(sample, 1):
        r = probe(path)
        if r is None:
            failures += 1
            continue
        rows.append(r)
        counts[r["verdict"]] += 1
        if i % 50 == 0:
            print(
                f"  ... {i}/{len(sample)}  confirmed_blind={counts['CONFIRMED_BLIND']}",
                flush=True,
            )

    candidates = [r for r in rows if r["verdict"] == "CONFIRMED_BLIND"]
    residue = [r for r in rows if r["verdict"] == "conformational_or_joint_blind"]
    # least flexible first: the oracle is trustworthy exactly where the molecule is rigid
    candidates.sort(key=lambda r: (not r["rigid"], r["rotatable_bonds"], -r["mirror_rmsd"]))
    rigid_candidates = [r for r in candidates if r["rigid"]]

    return {
        "dataset": str(dataset),
        "sample_seed": SAMPLE_SEED,
        "requested": n,
        "n_probed": len(rows),
        "load_failures": failures,
        "verdicts": {v: counts.get(v, 0) for v in VERDICTS},
        "n_candidates": len(candidates),
        "n_residue_ambiguous": len(residue),
        "residue_examples": residue[:25],
        "n_candidates_rigid": len(rigid_candidates),
        "rigidity_threshold_rotb": RIGID_ROTB_MAX,
        "candidates": candidates[:60],
    }


def render(h: dict) -> str:
    v = h["verdicts"]
    n = max(h["n_probed"], 1)
    lines = [
        "# Unknown-unknown hunt (Y3) -- chirality axes nobody has named",
        "",
        "Generator-free mirror twins over a seeded dataset sample, triaged so that only the",
        "unexplained residue is reported. A collapse the encoder makes for a KNOWN reason (a",
        "CIP stereocentre, or one of the P1/P2/P3 axes) is a re-detection. The signal is a",
        "collapse that **InChI** -- configuration-based and conformation-independent -- does",
        "NOT make: there the encoder has lost something standard cheminformatics keeps, with",
        "no conformational confound to explain it away.",
        "",
        f"- dataset: `{h['dataset']}`  (seed {h['sample_seed']})",
        f"- probed: **{h['n_probed']}** (load/encode failures: {h['load_failures']})",
        "",
        "| verdict | count | share | meaning |",
        "|---|---:|---:|---|",
        f"| achiral, encoder collapses | {v['achiral_ok']} | {v['achiral_ok'] / n:.1%} |"
        " correct invariance |",
        f"| chiral, encoder separates | {v['distinguished']} | {v['distinguished'] / n:.1%} |"
        " correct injectivity |",
        f"| collapse, CIP stereocentre present | {v['known_cip']} | {v['known_cip'] / n:.1%} |"
        " modelled stereo type |",
        f"| collapse, named P1/P2/P3 axis | {v['known_axis']} | {v['known_axis'] / n:.1%} |"
        " known blind spot |",
        f"| **collapse, InChI SEPARATES them** | **{v['CONFIRMED_BLIND']}** |"
        f" {v['CONFIRMED_BLIND'] / n:.1%} | **confirmed blind spot** |",
        f"| collapse, InChI agrees | {v['conformational_or_joint_blind']} |"
        f" {v['conformational_or_joint_blind'] / n:.1%} | ambiguous residue |",
        "",
        "## The two residues, and why they differ in strength",
        "",
        "**Confirmed blind spots** are collapses where InChI -- which is configuration-based",
        "and conformation-INDEPENDENT -- separates the two twins while the OIN does not. There",
        "is no conformational confound in that comparison, so it is a rigorous finding: the",
        "encoder has lost something standard cheminformatics keeps.",
        "",
        "**The ambiguous residue** is collapses where InChI agrees. Two very different things",
        "live there and the geometric oracle cannot separate them:",
        "",
        "* *conformational chirality* -- a twisted but freely-rotating fragment (an",
        "  unsubstituted biphenyl, a floppy arm). The mirror is a different CONFORMER of the",
        "  same isomer, so collapsing it is CORRECT: it is conformer invariance, not a blind",
        "  spot. Manual triage of the top-ranked rigid candidates found exactly this -- e.g.",
        "  `EDOQIZ` is an unsubstituted biphenyl on linear Au, whose axis the stereogenicity",
        "  gate already (correctly) refuses to call stereogenic.",
        "* *jointly-blind axes* -- P1/P2/P3 live here too, because InChI is itself blind to",
        "  metal Delta/Lambda and to atropisomerism. So a genuinely new axis that InChI also",
        "  misses would hide in this bucket.",
        "",
        "Separating those two requires a torsion-aware configurational test rather than a rigid",
        "one; until then this residue is a triage queue, not a finding.",
        "",
    ]
    if h["candidates"]:
        lines += [
            "## Triage queue (rigid first)",
            "",
            "| structure | rot. bonds | rigid? | mirror RMSD | OIN |",
            "|---|---:|---|---:|---|",
        ]
        for c in h["candidates"]:
            oin = c["oin"] if len(c["oin"]) <= 70 else c["oin"][:67] + "..."
            lines.append(
                f"| {c['name']} | {c['rotatable_bonds']} | {'yes' if c['rigid'] else 'no'} "
                f"| {c['mirror_rmsd']} | `{oin}` |"
            )
        lines.append("")
    else:
        lines += [
            "No unexplained collapses in this sample: every chiral structure the encoder",
            "collapsed was explained by a stereo element we already model. That is a negative",
            "result for the sample size, not proof that no unnamed axis exists.",
            "",
        ]
    lines += [
        "## Reproduce",
        "",
        "```",
        "PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 400",
        "```",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=400, help="sample size (<=0 for the whole set)")
    ap.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    args = ap.parse_args(argv)

    if not args.dataset.exists():
        raise SystemExit(f"dataset not found: {args.dataset}")
    h = hunt(args.dataset, args.n)
    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / "uu_hunt.json").write_text(json.dumps(h, indent=2) + "\n")
    md = render(h)
    (OUT_DIR / "uu_hunt.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
