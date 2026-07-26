"""Task B -- is the round-trip key blind to a donor swap? Answered from real geometry.

**The unsound probe this replaces, and why it was unsound.** The Y1 target map listed
"symmetric-donor swap" as a candidate key blind spot. The obvious test is to hand-write two OIN
strings with the donors in swapped slots and compare their keys; the key folds them, which
looks like a blind spot. It is not evidence of anything. That probe never established that its
two strings denote *different isomers* -- and for a square-planar complex with two identical
ancillary ligands, swapping slots 0<->1 is a reflection of a **planar, hence achiral** complex,
so the two strings name the same molecule and folding them is correct. (The octahedral variant
was worse: it put a short bidentate across *trans* slots, which is geometrically impossible.)

The rule the audit adopted for encoder claims applies to key claims too: **drive the comparison
from real geometry, with something independent certifying that the two are distinct isomers.**
This module does that two ways, so the verdict does not rest on one construction:

1. **Corpus pairs.** Two real crystal structures with the *same constitution* (identical
   canonical SMILES of the perceived complex) but a *different donor arrangement*. The
   arrangement is read off the geometry as the multiset of **trans donor pairs** -- which
   donors sit opposite which. For a fixed constitution that multiset is a configurational
   invariant: two complexes whose trans-pair multisets differ cannot be the same isomer, no
   matter how either one is oriented or how its ligands rotate. No OIN string is involved in
   that certification.

2. **Operator pairs.** :func:`tools.injectivity.twin_operators.swap_donor` applied to a real
   structure, exchanging two donors' coordination sites by a rigid motion. Distinctness is
   certified by the torsion-orbit oracle, and geometry by the vdW clash gate.

Run:
  PYTHONPATH=$PWD/src python -m tools.injectivity.positional_isomers --scan
  PYTHONPATH=$PWD/src python -m tools.injectivity.positional_isomers
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem

from oinsmiles.core.constants import TRANSITION_METALS_NUM

from .oracle import load_mol
from .twin_collision import _key, _silence_fds
from .twin_operators import enumerate_donor_swaps, probe_operator, swap_donor

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DATASET = Path("/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset")
OUT_DIR = REPO / "results-injectivity-y3"

#: donor-M-donor angle above which two donors count as trans.
TRANS_ANGLE_DEG = 150.0

#: A trans-pair multiset read off a real crystal geometry is only meaningful if it does not
#: depend on where the cutoff was put. Every signature is recomputed at each of these angles and
#: a pair only counts when the two members differ at ALL of them. Without this the first "pair"
#: the scan produced was BOCYEA_comp_0 vs BOCYEA_comp_1 -- two crystallographically independent
#: copies of the SAME compound, separated only by one donor-M-donor angle sitting either side of
#: 150 deg. It would have been reported as a folded positional isomer.
TRANS_ANGLE_LADDER = (140.0, 150.0, 160.0)

#: structures larger than this are skipped: perception cost, and big flexible ligands make the
#: "same constitution" test noisy without adding cases.
MAX_ATOMS = 70


def _metal(mol: Chem.Mol) -> int | None:
    metals = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM]
    return metals[0] if len(metals) == 1 else None


def arrangement_signature(mol: Chem.Mol, coords, cutoff: float = TRANS_ANGLE_DEG) -> tuple:
    """The multiset of trans donor pairs, as sorted element pairs.

    A configurational invariant of the coordination sphere: it is unchanged by any rotation of
    the whole complex and by any rotation of any bond, and it differs between cis/trans and
    between fac/mer. Two structures of the same constitution whose signatures differ are
    therefore different isomers -- certified without reference to the OIN.
    """
    m = _metal(mol)
    if m is None:
        return ()
    donors = [n.GetIdx() for n in mol.GetAtomWithIdx(m).GetNeighbors()]
    out = []
    for i, a in enumerate(donors):
        for b in donors[i + 1 :]:
            va, vb = coords[a] - coords[m], coords[b] - coords[m]
            cos = float(va @ vb) / (float((va @ va) ** 0.5) * float((vb @ vb) ** 0.5))
            if math.degrees(math.acos(max(-1.0, min(1.0, cos)))) > cutoff:
                out.append(
                    tuple(
                        sorted(
                            (
                                mol.GetAtomWithIdx(a).GetSymbol(),
                                mol.GetAtomWithIdx(b).GetSymbol(),
                            )
                        )
                    )
                )
    return tuple(sorted(out))


def signature_ladder(mol: Chem.Mol, coords) -> list[tuple]:
    """The signature at every cutoff in :data:`TRANS_ANGLE_LADDER` (see that constant)."""
    return [arrangement_signature(mol, coords, c) for c in TRANS_ANGLE_LADDER]


def _robustly_different(la: list[tuple], lb: list[tuple]) -> bool:
    """The two arrangements differ, and differ for a reason other than the cutoff.

    Three conditions, each removing one way a crystal geometry can fake an isomer difference:

    * they differ at **every** cutoff, so no single borderline angle is doing the work;
    * the number of trans pairs matches at every cutoff -- a differing *count* means one donor
      pair straddled the boundary or the two have different coordination numbers, neither of
      which is a positional isomerism;
    * the count is stable across the ladder within each structure, so neither signature is
      itself sitting on a boundary.
    """
    if len({len(s) for s in la}) != 1 or len({len(s) for s in lb}) != 1:
        return False
    if len(la[0]) != len(lb[0]):
        return False
    return all(a != b for a, b in zip(la, lb, strict=True))


def constitution_key(mol: Chem.Mol) -> str:
    """Canonical SMILES of the perceived complex, stereo stripped -- 'which compound is this'."""
    probe = Chem.Mol(mol)
    Chem.RemoveStereochemistry(probe)
    return Chem.MolToSmiles(Chem.RemoveHs(probe))


def _prefilter(path: Path) -> tuple | None:
    """Cheap, rdkit-free triage: (metal, n_atoms, formula) for single-metal structures."""
    lines = path.read_text().splitlines()
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        return None
    if n > MAX_ATOMS:
        return None
    syms = [ln.split()[0] for ln in lines[2 : 2 + n] if ln.split()]
    if len(syms) != n:
        return None
    pt = Chem.GetPeriodicTable()
    metals = [s for s in syms if pt.GetAtomicNumber(s) in TRANSITION_METALS_NUM]
    if len(metals) != 1:
        return None
    return (metals[0], "".join(f"{k}{v}" for k, v in sorted(Counter(syms).items())))


def scan(dataset: Path, limit: int = 0, checkpoint: Path | None = None) -> dict:
    """Find real positional-isomer pairs: same constitution, different trans-pair signature."""
    files = sorted(p for sub in ("cat", "photo") for p in (dataset / sub).glob("*.xyz"))
    if limit:
        files = files[:limit]

    # pass 1 -- rdkit-free: only formulas that occur more than once can hold a pair
    buckets: dict[tuple, list[Path]] = defaultdict(list)
    for p in files:
        try:
            k = _prefilter(p)
        except Exception:
            continue
        if k:
            buckets[k].append(p)
    shortlist = [p for group in buckets.values() if len(group) > 1 for p in group]
    print(
        f"  pass 1: {len(files)} files -> {len(buckets)} formula buckets -> "
        f"{len(shortlist)} structures in a bucket with >1 member",
        flush=True,
    )

    # pass 2 -- perceive only the shortlist, group by constitution.
    #
    # Checkpointed to a JSONL sidecar and resumable, because pass 2 is hours of perception on a
    # loaded machine and the whole result used to be written only at the very end: two runs were
    # lost at 2200/3116 and 150/300 to a harness task timeout, each discarding every perception
    # it had done. One append per structure costs nothing and makes a restart free.
    groups: dict[str, list[dict]] = defaultdict(list)
    perceive_fail = 0
    done: dict[str, dict] = {}
    if checkpoint and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:  # a torn final line from a hard kill
                continue
            done[rec["path"]] = rec
        print(f"  pass 2: resuming, {len(done)} structures already perceived", flush=True)
    fh = checkpoint.open("a") if checkpoint else None
    try:
        for i, p in enumerate(shortlist, 1):
            if i % 200 == 0:
                print(f"  pass 2: {i}/{len(shortlist)}", flush=True)
            rec = done.get(str(p))
            if rec is None:
                try:
                    with _silence_fds():
                        mol, coords = load_mol(p)
                        Chem.SanitizeMol(mol)
                        ck = constitution_key(mol)
                        ladder = signature_ladder(mol, coords)
                except Exception:
                    rec = {"path": str(p), "failed": True}
                else:
                    rec = {
                        "path": str(p),
                        "name": p.stem,
                        "refcode": p.stem.split("_")[0],
                        "constitution": ck,
                        "signature": [list(s) for s in ladder[1]],
                        "ladder": [[list(s) for s in sig] for sig in ladder],
                    }
                if fh:
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
            if rec.get("failed"):
                perceive_fail += 1
                continue
            if not rec["signature"]:
                continue
            groups[rec["constitution"]].append(
                {k: rec[k] for k in ("path", "name", "refcode", "signature", "ladder")}
            )
    finally:
        if fh:
            fh.close()

    pairs, rejected = [], 0
    for ck, members in groups.items():
        for i, ma in enumerate(members):
            for mb in members[i + 1 :]:
                # two crystallographically independent copies of one compound are not isomers
                if ma["refcode"] == mb["refcode"]:
                    continue
                la = [tuple(tuple(x) for x in s) for s in ma["ladder"]]
                lb = [tuple(tuple(x) for x in s) for s in mb["ladder"]]
                if la[1] == lb[1]:
                    continue
                if not _robustly_different(la, lb):
                    rejected += 1
                    continue
                pairs.append(
                    {
                        "constitution": ck,
                        "a": ma,
                        "b": mb,
                        "signature_a": ma["signature"],
                        "signature_b": mb["signature"],
                    }
                )
    return {
        "n_files": len(files),
        "n_shortlist": len(shortlist),
        "n_perceive_failures": perceive_fail,
        "n_constitution_groups": len(groups),
        "n_rejected_cutoff_artifacts": rejected,
        "n_pairs": len(pairs),
        "pairs": pairs,
    }


def compare_pair(pair: dict) -> dict:
    """Encode both members from geometry; compare raw strings and round-trip keys."""
    from oinsmiles import XYZToSMILES

    with _silence_fds():
        oin_a = XYZToSMILES().convert(pair["a"]["path"])
        oin_b = XYZToSMILES().convert(pair["b"]["path"])
    raw_equal = oin_a == oin_b
    key_equal = _key(oin_a) == _key(oin_b)
    return {
        "a": pair["a"]["name"],
        "b": pair["b"]["name"],
        "constitution": pair["constitution"],
        "signature_a": pair["signature_a"],
        "signature_b": pair["signature_b"],
        # the independent certification: same constitution, different trans-pair multiset
        "oracle_distinct": pair["signature_a"] != pair["signature_b"],
        "oin_a": oin_a,
        "oin_b": oin_b,
        "raw_equal": raw_equal,
        "key_equal": key_equal,
        "verdict": "KEY FOLDS a distinct positional isomer" if key_equal else "key separates them",
    }


def operator_pairs(paths: list[Path], *, restarts: int = 8) -> list[dict]:
    """Second line of evidence: swap_donor applied to real geometry."""
    out = []
    for p in paths:
        try:
            with _silence_fds():
                mol, coords = load_mol(p)
        except Exception as e:
            out.append({"name": p.stem, "error": repr(e)[:150]})
            continue
        for a, b in enumerate_donor_swaps(mol):
            o = probe_operator(p, swap_donor(mol, coords, a, b), restarts=restarts)
            out.append(o.to_dict())
            print(
                f"  {p.stem:24s} {o.detail:18s} geom_ok={o.geometry_ok} "
                f"distinct={o.oracle_distinct} raw_eq={o.raw_equal} key_eq={o.key_equal} "
                f"-> {o.verdict}",
                flush=True,
            )
    return out


def render(scan_res: dict | None, compared: list[dict], ops: list[dict]) -> str:
    folds = [c for c in compared if c["oracle_distinct"] and c["key_equal"]]
    op_scored = [o for o in ops if o.get("geometry_ok") and o.get("oracle_distinct")]
    op_folds = [o for o in op_scored if o.get("key_equal")]
    if compared or op_scored:
        verdict = "CONFIRMED" if (folds or op_folds) else "REFUTED"
    else:
        verdict = "UNDETERMINED"
    lines = [
        "# Donor swap and the round-trip key (Lane 7 / Task B)",
        "",
        f'## Verdict: **{verdict}** -- "the key folds a donor swap"',
        "",
        "Both lines of evidence are driven from real 3D geometry, and in both the claim that the",
        "two structures are *different isomers* is certified without reference to any OIN string.",
        "",
        "### Line 1 -- real positional-isomer pairs from the corpus",
        "",
        "Same constitution (identical canonical SMILES of the perceived complex), different",
        "**trans-donor multiset**. For a fixed constitution that multiset is a configurational",
        "invariant -- unchanged by any rotation of the complex or of any bond -- so two members",
        "whose multisets differ cannot be the same isomer.",
        "",
    ]
    if scan_res:
        lines += [
            f"- corpus scanned: {scan_res['n_files']} structures",
            f"- shortlisted (formula shared with another structure): {scan_res['n_shortlist']}",
            f"- constitution groups formed: {scan_res['n_constitution_groups']}",
            f"- **positional-isomer pairs found: {scan_res['n_pairs']}**",
            "",
        ]
    if compared:
        lines += [
            "| A | B | trans(A) | trans(B) | distinct? | raw equal | key equal |",
            "|---|---|---|---|---|---|---|",
        ]
        for c in compared:
            ta = " ".join("-".join(s) for s in c["signature_a"]) or "(none)"
            tb = " ".join("-".join(s) for s in c["signature_b"]) or "(none)"
            lines.append(
                f"| {c['a']} | {c['b']} | {ta} | {tb} | {c['oracle_distinct']} "
                f"| {c['raw_equal']} | {c['key_equal']} |"
            )
        lines.append("")
    else:
        lines += [
            "No pair survived. That is an honest negative for this corpus, not evidence that the",
            "key is sound -- see line 2.",
            "",
        ]
    lines += [
        "### Line 2 -- swap_donor on real geometry",
        "",
        "`swap_donor` exchanges two donors' coordination sites by a rigid motion of each ligand,",
        "so every ligand's internal geometry is untouched and only *which site each donor",
        "occupies* changes. Symmetry-equivalent donor pairs are skipped -- exchanging those is",
        "the identity on the isomer, and that is precisely the mistake the hand-written probe",
        "made. Distinctness comes from the torsion-orbit oracle, geometry from the vdW clash gate.",
        "",
        "| structure | swap | geometry ok | distinct isomer | raw equal | key equal | verdict |",
        "|---|---|---|---|---|---|---|",
    ]
    for o in ops:
        if "error" in o:
            continue
        lines.append(
            f"| {o['name']} | {o['detail']} | {o['geometry_ok']} | {o['oracle_distinct']} "
            f"| {o['raw_equal']} | {o['key_equal']} | {o['verdict']} |"
        )
    lines += [
        "",
        f"Distinct-isomer swaps scored: **{len(op_scored)}**; of those, folded by the key: "
        f"**{len(op_folds)}**.",
        "",
        "A swap of two donors that sit *trans* to each other is a 180 deg rotation of the whole",
        "complex, so it is the same isomer and the oracle correctly reports `distinct=False`.",
        "Those rows are the built-in negative control: an instrument that called them distinct",
        "would be manufacturing blind spots.",
        "",
        "## Reproduce",
        "",
        "```",
        "PYTHONPATH=$PWD/src python -m tools.injectivity.positional_isomers --scan",
        "```",
        "",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    ap.add_argument("--scan", action="store_true", help="run the corpus scan (slow)")
    ap.add_argument("--limit", type=int, default=0, help="only the first N corpus files")
    ap.add_argument("--max-pairs", type=int, default=40)
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)

    scan_res = None
    compared: list[dict] = []
    cache = args.out / "positional_isomer_scan.json"
    if args.scan or not cache.exists():
        args.out.mkdir(parents=True, exist_ok=True)
        scan_res = scan(args.dataset, args.limit, checkpoint=args.out / "pass2_checkpoint.jsonl")
        args.out.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(scan_res, indent=2) + "\n")
    else:
        scan_res = json.loads(cache.read_text())
    print(f"  pairs found: {scan_res['n_pairs']}", flush=True)
    for pair in scan_res["pairs"][: args.max_pairs]:
        try:
            compared.append(compare_pair(pair))
        except Exception as e:
            print(f"  ! {pair['a']['name']}/{pair['b']['name']}: {e!r}"[:160])
        if compared:
            c = compared[-1]
            print(f"  {c['a']} vs {c['b']}: raw_eq={c['raw_equal']} key_eq={c['key_equal']}")

    fx = REPO / "tests" / "fixtures"
    ops = operator_pairs(
        [
            p
            for p in (
                fx / "PtMeNH3ClBr-Cis.xyz",
                fx / "CisPlatin.xyz",
                fx / "JEGKOW.xyz",
                fx / "TransPlatin.xyz",
                fx / "FeH2(CO)4.xyz",
            )
            if p.exists()
        ]
    )

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "positional_isomers.json").write_text(
        json.dumps({"scan": scan_res, "pairs": compared, "operator_pairs": ops}, indent=2) + "\n"
    )
    md = render(scan_res, compared, ops)
    (args.out / "positional_isomers.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
