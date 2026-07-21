#!/usr/bin/env python3
"""E2c: static invariant sweeps over the real OIN corpus.

Parse-only and single-core. Runs the real OIN parser and the real
``convert_parsed_to_msmiles`` over every stored ``.oin`` string, but never
generates 3D, so it is cheap enough to run on a loaded machine.

Where the synthetic battery (``tools/synthetic_oin_battery.py``) *detects*
handoff defects with unambiguous attribution, this *sizes* them on real
chemistry and cross-tabs them against the round-trip outcome.

Checks:

  H2b  slot drop     -- slot markers present in the OIN string vs coordination
                        vectors surviving into ParsedOIN. The parser drops a slot
                        silently when no geometry template resolves it.
  H3b  slot bijection -- on real molecules, does the nearest-vector match ever
                        put two ligands in the same ``ligand_parts`` index (which
                        would overwrite one of them)?
  H3   heavy formula -- heavy-atom composition of the OIN fragments vs the emitted
                        m-SMILES. Heavy atoms are unambiguous: unlike hydrogen,
                        no phantom-implicit convention applies, so a mismatch is a
                        real loss rather than a convention difference.
  H3x  handoff raise -- how often the handoff refuses outright, and with what.
  H4f  actinide      -- complexes whose metal is rewritten to a lanthanide analogue.
  H4e  multiplicity  -- complexes whose metal admits a high-spin state, all of which
                        are forced to singlet/doublet by ``om.py``.

Usage:
    uv run python tools/elimination_static_sweeps.py
    uv run python tools/elimination_static_sweeps.py --out results/e2c.json --limit 500
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    OIN_TO_METALLOGEN_GEO,
    convert_parsed_to_msmiles,
)
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.generator3d import globalvars  # noqa: E402

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
DEFAULT_RESULTS = REPO / "tmCAT-tmPHOTO_xyz_dataset" / "results-capstone-v042"

SLOT_RE = re.compile(r"\{(\d+)[><]?\}")
METAL_TAG_RE = re.compile(r"\[([A-Z][a-z]?)(?:_([A-Z]{2,3}))?")

ACTINIDES = {
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
}
# First-row TMs plus the f-block: the metals that routinely carry S > 1/2 and so
# are mis-specified by a multiplicity forced to singlet/doublet.
HIGH_SPIN_CAPABLE = {
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Mo",
    "Ru",
    "Gd",
    "Eu",
    "Sm",
    "Nd",
    "Dy",
    "Tb",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Ce",
    "Pr",
}


TETRAHEDRAL_RE = re.compile(r"@{1,2}(?![A-Z])")
NAMED_STEREO_RE = re.compile(r"@(?:SP|TB|OH|AL)\d*")
DIRECTIONAL_RE = re.compile(r"[/\\]")


def stereo_profile(smiles_fragments: list[str]) -> dict[str, int]:
    """Count stereo markers as they appear in SMILES text.

    This is a *screen*, not a verdict. A textual count over-reports: RDKit's
    canonical output legitimately re-expresses the same stereochemistry with a
    different ``@``/``@@`` distribution when atom ordering changes, so a
    non-zero delta here does not imply information was lost. Every textual
    mismatch is adjudicated by :func:`perceived_stereocentres`, which compares
    what the graph actually encodes.
    """
    joined = " ".join(smiles_fragments)
    return {
        "tetrahedral": len(TETRAHEDRAL_RE.findall(joined)),
        "named": len(NAMED_STEREO_RE.findall(joined)),
        "directional": len(DIRECTIONAL_RE.findall(joined)),
    }


def perceived_stereocentres(smiles_fragments: list[str]) -> int:
    """Number of assigned stereocentres RDKit perceives across the fragments."""
    total = 0
    for frag in smiles_fragments:
        mol = Chem.MolFromSmiles(frag, sanitize=False)
        if mol is None:
            continue
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol, catchErrors=True)
            total += len(
                Chem.FindMolChiralCenters(
                    mol, includeUnassigned=False, useLegacyImplementation=False
                )
            )
        except Exception:
            continue
    return total


def heavy_formula(smiles_fragments: list[str]) -> Counter:
    """Heavy-atom composition of a list of SMILES fragments."""
    out: Counter = Counter()
    for frag in smiles_fragments:
        mol = Chem.MolFromSmiles(frag, sanitize=False)
        if mol is None:
            out["__unparseable__"] += 1
            continue
        for atom in mol.GetAtoms():
            out[atom.GetSymbol()] += 1
    return out


def analyse(molecule: str, oin: str, status: str) -> dict[str, Any]:
    row: dict[str, Any] = {"molecule": molecule, "status": status, "oin_len": len(oin)}

    m = METAL_TAG_RE.match(oin)
    row["metal"] = m.group(1) if m else None
    row["geo"] = m.group(2) if m and m.group(2) else None
    row["actinide"] = row["metal"] in ACTINIDES
    row["high_spin_capable"] = row["metal"] in HIGH_SPIN_CAPABLE

    declared_slots = {int(s) for s in SLOT_RE.findall(oin)}
    row["n_slots_declared"] = len(declared_slots)

    try:
        parsed = OINParser().parse(oin)
    except Exception as exc:
        row.update({"stage": "parse", "error": f"{type(exc).__name__}: {exc}"})
        return row

    resolved_slots = {v.slot for v in parsed.vectors}
    row["n_slots_resolved"] = len(resolved_slots)
    row["slots_dropped"] = sorted(declared_slots - resolved_slots)
    row["n_slots_dropped"] = len(row["slots_dropped"])
    row["geo_code"] = parsed.geo_code
    row["geo_supported"] = parsed.geo_code in OIN_TO_METALLOGEN_GEO
    row["n_fragments"] = len(parsed.fragments)

    # H3b on real data: would two ligands land in the same ligand_parts index?
    if parsed.geo_code in OIN_TO_METALLOGEN_GEO:
        mg = globalvars.known_geometries_vector_dict[OIN_TO_METALLOGEN_GEO[parsed.geo_code]]
        firsts = []
        for frag_idx in range(len(parsed.fragments)):
            if frag_idx == parsed.metal_fragment_idx:
                continue
            fv = [v for v in parsed.vectors if v.fragment_idx == frag_idx]
            if not fv:
                continue
            firsts.append(int(np.argmin(np.linalg.norm(mg - np.array(fv[0].vector), axis=1))))
        dupes = {k: c for k, c in Counter(firsts).items() if c > 1}
        row["ligand_slot_collision"] = bool(dupes)
        row["n_ligands_lost_to_collision"] = sum(c - 1 for c in dupes.values())

    try:
        msmiles = convert_parsed_to_msmiles(parsed)
    except Exception as exc:
        row.update({"stage": "handoff", "error": f"{type(exc).__name__}: {exc}"})
        return row

    row["stage"] = "ok"
    parts = msmiles.split("|")
    ligand_parts = parts[1:-1]
    row["n_ligands_msmiles"] = len(ligand_parts)
    row["n_ligands_oin"] = len(parsed.fragments) - 1

    oin_heavy = heavy_formula(
        [f for i, f in enumerate(parsed.fragments) if i != parsed.metal_fragment_idx]
    )
    mg_heavy = heavy_formula(ligand_parts)
    row["heavy_formula_match"] = oin_heavy == mg_heavy
    if oin_heavy != mg_heavy:
        diff = {k: mg_heavy.get(k, 0) - oin_heavy.get(k, 0) for k in set(oin_heavy) | set(mg_heavy)}
        row["heavy_formula_diff"] = {k: v for k, v in diff.items() if v}

    # H3a: stereo annotation carried through the handoff. Ligand fragments only --
    # the metal fragment is compared separately because it is the one the
    # @SP-stripping regex in convert_parsed_to_msmiles actually touches.
    oin_ligands = [f for i, f in enumerate(parsed.fragments) if i != parsed.metal_fragment_idx]
    oin_stereo = stereo_profile(oin_ligands)
    mg_stereo = stereo_profile(ligand_parts)
    row["stereo_oin"] = oin_stereo
    row["stereo_msmiles"] = mg_stereo
    row["stereo_match"] = oin_stereo == mg_stereo
    if not row["stereo_match"]:
        row["stereo_diff"] = {
            k: mg_stereo[k] - oin_stereo[k] for k in oin_stereo if mg_stereo[k] != oin_stereo[k]
        }
        # Adjudicate: a textual delta is only a real loss if the graph lost a
        # stereocentre too. Perception is the expensive path, so it runs only
        # on the molecules the cheap screen flagged.
        oin_centres = perceived_stereocentres(oin_ligands)
        mg_centres = perceived_stereocentres(ligand_parts)
        row["stereo_centres_oin"] = oin_centres
        row["stereo_centres_msmiles"] = mg_centres
        row["stereo_real_loss"] = mg_centres < oin_centres

    metal_frag = parsed.fragments[parsed.metal_fragment_idx]
    row["metal_named_stereo"] = len(NAMED_STEREO_RE.findall(metal_frag))
    row["metal_part_msmiles"] = parts[0]
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--limit", type=int, help="only the first N molecules (smoke test)")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    summary_path = args.results_dir / "summary_roundtrip.json"
    status_by_mol = {}
    if summary_path.exists():
        for r in json.loads(summary_path.read_text()):
            status_by_mol[r["molecule"]] = r.get("status", "?")

    oin_files = sorted((args.results_dir / "structures").glob("*.oin"))
    if args.limit:
        oin_files = oin_files[: args.limit]
    if not oin_files:
        raise SystemExit(f"no .oin files under {args.results_dir / 'structures'}")

    rows = []
    for path in oin_files:
        mol = path.stem
        try:
            oin = path.read_text().strip()
        except OSError:
            continue
        rows.append(analyse(mol, oin, status_by_mol.get(mol, "?")))

    n = len(rows)
    ok = [r for r in rows if r.get("stage") == "ok"]

    def pct(k: int) -> str:
        return f"{k:6d}  ({100.0 * k / n:5.2f}%)"

    print("=" * 78)
    print(f"E2c STATIC INVARIANT SWEEPS -- {n} OIN strings from {args.results_dir.name}")
    print("=" * 78)

    print("\nSTAGE REACHED")
    for stage, c in Counter(r.get("stage", "parse_error") for r in rows).most_common():
        print(f"  {stage:<12} {pct(c)}")

    print("\nH2b  SLOT DROP (parser silently discards an unresolvable slot)")
    dropped = [r for r in rows if r.get("n_slots_dropped")]
    print(f"  molecules with >=1 dropped slot: {pct(len(dropped))}")
    if dropped:
        print(f"  total slots dropped: {sum(r['n_slots_dropped'] for r in dropped)}")
        by_geo = Counter(r.get("geo_code") for r in dropped)
        print(f"  by geometry: {dict(by_geo.most_common(10))}")
        for r in dropped[:5]:
            print(f"    {r['molecule']:<22} geo={r.get('geo_code')} dropped={r['slots_dropped']}")

    print("\nH3b  LIGAND SLOT COLLISION (one ligand overwrites another)")
    coll = [r for r in rows if r.get("ligand_slot_collision")]
    print(f"  molecules with a collision: {pct(len(coll))}")
    if coll:
        print(f"  total ligands lost: {sum(r['n_ligands_lost_to_collision'] for r in coll)}")
        for r in coll[:10]:
            print(f"    {r['molecule']:<22} geo={r.get('geo_code')}")

    print("\nH3   HEAVY-ATOM FORMULA PRESERVED THROUGH THE HANDOFF")
    bad = [r for r in ok if not r.get("heavy_formula_match")]
    print(f"  of {len(ok)} that reached m-SMILES, mismatches: {len(bad)}")
    if bad:
        kinds: Counter = Counter()
        for r in bad:
            kinds[json.dumps(r.get("heavy_formula_diff", {}), sort_keys=True)] += 1
        for k, c in kinds.most_common(10):
            print(f"    {c:5d}  {k}")
        for r in bad[:5]:
            print(f"    e.g. {r['molecule']} {r.get('heavy_formula_diff')}")

    print("\n     LIGAND COUNT PRESERVED")
    lost = [r for r in ok if r.get("n_ligands_msmiles") != r.get("n_ligands_oin")]
    print(f"  ligand-count mismatches: {len(lost)} of {len(ok)}")
    for r in lost[:10]:
        print(f"    {r['molecule']:<22} oin={r['n_ligands_oin']} msmiles={r['n_ligands_msmiles']}")

    print("\nH3a  STEREO ANNOTATION CARRIED THROUGH THE HANDOFF")
    with_stereo = [
        r for r in ok if r.get("stereo_oin") and any(v for v in r["stereo_oin"].values())
    ]
    stereo_bad = [r for r in ok if r.get("stereo_match") is False]
    print(f"  molecules carrying any stereo annotation: {pct(len(with_stereo))}")
    print(f"  textual-screen mismatches: {len(stereo_bad)} of {len(ok)}")
    if stereo_bad:
        kinds = Counter(json.dumps(r.get("stereo_diff", {}), sort_keys=True) for r in stereo_bad)
        for k, c in kinds.most_common(8):
            print(f"    {c:5d}  {k}")
        real = [r for r in stereo_bad if r.get("stereo_real_loss")]
        print(f"    ADJUDICATED as real stereocentre loss: {len(real)}")
        for r in real[:10]:
            print(
                f"      {r['molecule']:<22} centres "
                f"{r.get('stereo_centres_oin')} -> {r.get('stereo_centres_msmiles')}"
            )
        if not real:
            print("      (all textual deltas are canonicalisation, not loss)")
        st = Counter(r["status"] for r in stereo_bad)
        fr = 100.0 * (len(stereo_bad) - st.get("success", 0)) / len(stereo_bad)
        print(f"    fail_rate among textual mismatches: {fr:5.1f}%  {dict(st)}")

    metal_stereo = [r for r in ok if r.get("metal_named_stereo")]
    print(f"  metal fragments carrying @SP/@TB/@OH (the stripped tag): {pct(len(metal_stereo))}")
    if metal_stereo:
        for r in metal_stereo[:5]:
            print(f"      {r['molecule']:<22} -> m-SMILES metal part {r['metal_part_msmiles']!r}")

    print("\nH3x  HANDOFF REFUSALS")
    refused = [r for r in rows if r.get("stage") == "handoff"]
    print(f"  molecules the handoff rejected: {pct(len(refused))}")
    kinds = Counter((r.get("error") or "").split(":")[0] for r in refused)
    for k, c in kinds.most_common(10):
        print(f"    {c:5d}  {k}")

    print("\nH4f  ACTINIDE (silently rewritten to a lanthanide analogue)")
    act = [r for r in rows if r.get("actinide")]
    print(f"  complexes with an actinide centre: {pct(len(act))}")
    if act:
        print(f"  metals: {dict(Counter(r['metal'] for r in act))}")

    print("\nH4e  MULTIPLICITY FORCED TO SINGLET/DOUBLET")
    hs = [r for r in rows if r.get("high_spin_capable")]
    print(f"  complexes on a high-spin-capable metal: {pct(len(hs))}")
    print(f"  metals: {dict(Counter(r['metal'] for r in hs).most_common(10))}")

    print("\nCROSS-TAB: defect vs round-trip outcome")
    for label, subset in (
        ("slot dropped", dropped),
        ("slot collision", coll),
        ("heavy formula mismatch", bad),
        ("handoff refused", refused),
    ):
        if not subset:
            print(f"  {label:<26} n=0")
            continue
        st = Counter(r["status"] for r in subset)
        fail_rate = 100.0 * (len(subset) - st.get("success", 0)) / len(subset)
        print(f"  {label:<26} n={len(subset):5d}  fail_rate={fail_rate:5.1f}%  {dict(st)}")
    base = Counter(r["status"] for r in rows)
    print(
        f"  {'BASELINE (all molecules)':<26} n={n:5d}  "
        f"fail_rate={100.0 * (n - base.get('success', 0)) / n:5.1f}%"
    )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({"rows": rows}, indent=2, default=str))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
