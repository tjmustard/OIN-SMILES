#!/usr/bin/env python3
"""C4: synthetic OIN battery -- tests the OIN -> m-SMILES handoff in isolation.

Read-only with respect to the library; no 3D is generated.

Every case here is a hand-constructed OIN string whose correct composition,
charge, donor set and coordination number are known *by construction*. That
removes every other source of error -- no crystal perception, no template
matching, no embedding, no force field -- so a discrepancy can only come from
``convert_parsed_to_msmiles``.

This is deliberately the *detector* for the handoff hypothesis; the corpus-wide
audit is the *sizer*. If the handoff corrupts information, it shows up here with
unambiguous attribution and no statistics required.

Two independent checks:

  1. SLOT BIJECTION (pure geometry, no molecules). For every supported geometry,
     map each OIN template slot to its nearest MetalloGen slot and assert the
     assignment is injective. A collision means two ligands resolve to the same
     ``ligand_parts`` index and one is silently overwritten.

  2. COMPOSITION FIDELITY. Build an OIN string from known parts, run the real
     handoff, and compare the m-SMILES against the by-construction truth:
     heavy-atom formula, hydrogen count, total charge, ligand count, donor
     elements, and the set of occupied slots.

Usage:
    uv run python tools/synthetic_oin_battery.py
    uv run python tools/synthetic_oin_battery.py --out results/c4_battery.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    OIN_TO_METALLOGEN_GEO,
    convert_oin_to_msmiles,
)
from oinsmiles.generation.oin_parser import TEMPLATES  # noqa: E402
from oinsmiles.generator3d import globalvars  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Coordination number implied by each geometry code.
GEO_CN = {code: len(vecs) for code, vecs in TEMPLATES.items()}


# --------------------------------------------------------------------------
# Ligand catalogue -- OIN fragment text plus its known-by-construction content.
#
# ``heavy`` is the heavy-atom formula and ``h`` the hydrogen count the OIN
# *intends*: a bare donor atom's implicit H is phantom (the metal bond replaces
# it), so a bare ``O``/``S``/``N`` donor contributes 0 H even though RDKit would
# read the free fragment as water/H2S/ammonia. That reinterpretation is exactly
# what the handoff must get right, so it is declared here by hand.
#
# Formal charge is NOT declared: in the OIN/m-SMILES convention an X-type donor
# is written neutral and the complex charge lives on the metal (``[Cl]``, not
# ``[Cl-]``). Charge is explicit and unambiguous in SMILES, so the expected value
# is read back off the OIN fragment itself rather than asserted from chemical
# intuition -- asserting "chloride is -1" tests the convention, not the code.
# --------------------------------------------------------------------------
LIGANDS: dict[str, dict[str, Any]] = {
    "chloride": {"frag": "[Cl]{{{s}}}", "heavy": {"Cl": 1}, "h": 0, "donor": "Cl"},
    "ammine": {"frag": "N{{{s}}}", "heavy": {"N": 1}, "h": 3, "donor": "N"},
    "aqua": {"frag": "[OH2]{{{s}}}", "heavy": {"O": 1}, "h": 2, "donor": "O"},
    "carbonyl": {"frag": "C{{{s}}}#O", "heavy": {"C": 1, "O": 1}, "h": 0, "donor": "C"},
    "cyanide": {"frag": "[C-]{{{s}}}#N", "heavy": {"C": 1, "N": 1}, "h": 0, "donor": "C"},
    "hydroxo": {"frag": "[OH]{{{s}}}", "heavy": {"O": 1}, "h": 1, "donor": "O"},
    "methyl": {"frag": "[CH3]{{{s}}}", "heavy": {"C": 1}, "h": 3, "donor": "C"},
    # --- adversarial: bare donors whose H count the adapter must reinterpret ---
    "oxo_bare_O": {"frag": "O{{{s}}}", "heavy": {"O": 1}, "h": 0, "donor": "O"},
    "thiolate_bare_S": {"frag": "S{{{s}}}", "heavy": {"S": 1}, "h": 0, "donor": "S"},
    "nitride_bare_N": {"frag": "N{{{s}}}", "heavy": {"N": 1}, "h": 0, "donor": "N"},
    "amido_bare_N": {"frag": "CN{{{s}}}", "heavy": {"C": 1, "N": 1}, "h": 3, "donor": "N"},
}


def oin_fragment_charge(key: str) -> int:
    """Formal charge as literally written in the OIN fragment (markers stripped)."""
    frag = LIGANDS[key]["frag"].format(s=0).replace("{0}", "")
    mol = Chem.MolFromSmiles(frag, sanitize=False)
    if mol is None:
        raise ValueError(f"catalogue fragment for {key!r} does not parse: {frag}")
    mol.UpdatePropertyCache(strict=False)
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


METALS = ["Pt", "Pd", "Fe", "Ru", "Ir", "Zr", "Y", "Ti"]


def build_oin(metal: str, geo: str, ligand_keys: list[str]) -> str:
    """Assemble an inline OIN string placing each ligand at consecutive slots."""
    frags = [f"[{metal}_{geo}]"]
    for slot, key in enumerate(ligand_keys):
        frags.append(LIGANDS[key]["frag"].format(s=slot))
    return ".".join(frags)


def expected_content(ligand_keys: list[str]) -> dict[str, Any]:
    heavy: Counter = Counter()
    h = 0
    charge = 0
    donors: Counter = Counter()
    for key in ligand_keys:
        lig = LIGANDS[key]
        heavy.update(lig["heavy"])
        h += lig["h"]
        charge += oin_fragment_charge(key)
        donors[lig["donor"]] += 1
    return {
        "heavy": dict(heavy),
        "h": h,
        "ligand_charge": charge,
        "donors": dict(donors),
        "n_ligands": len(ligand_keys),
    }


def msmiles_content(msmiles: str) -> dict[str, Any]:
    """Extract composition from the emitted m-SMILES ``metal|lig...|geo``."""
    parts = msmiles.split("|")
    metal_part, ligand_parts = parts[0], parts[1:-1]

    heavy: Counter = Counter()
    h = 0
    charge = 0
    donors: Counter = Counter()
    slots: list[int] = []

    for lig in ligand_parts:
        mol = Chem.MolFromSmiles(lig, sanitize=False)
        if mol is None:
            return {"parse_error": lig}
        mol.UpdatePropertyCache(strict=False)
        for atom in mol.GetAtoms():
            heavy[atom.GetSymbol()] += 1
            h += atom.GetTotalNumHs()
            charge += atom.GetFormalCharge()
            if atom.GetAtomMapNum() > 0:
                slots.append(atom.GetAtomMapNum())
                donors[atom.GetSymbol()] += 1

    return {
        "heavy": dict(heavy),
        "h": h,
        "ligand_charge": charge,
        "donors": dict(donors),
        "n_ligands": len(ligand_parts),
        "slots": sorted(slots),
        "metal_part": metal_part,
    }


# --------------------------------------------------------------------------
# Check 1: slot bijection
# --------------------------------------------------------------------------
def check_slot_bijection() -> list[dict[str, Any]]:
    """H3b: does the nearest-vector match ever collapse two OIN slots into one?"""
    results = []
    for oin_code, oin_vecs in sorted(TEMPLATES.items()):
        geo = OIN_TO_METALLOGEN_GEO.get(oin_code)
        if not geo:
            results.append({"geo": oin_code, "status": "unmapped", "injective": None})
            continue
        mg = globalvars.known_geometries_vector_dict[geo]
        mapping = [int(np.argmin(np.linalg.norm(mg - v, axis=1))) for v in oin_vecs]
        counts = Counter(mapping)
        collisions = {k: c for k, c in counts.items() if c > 1}
        results.append(
            {
                "geo": oin_code,
                "metallogen_geo": geo,
                "n_oin_slots": len(oin_vecs),
                "n_mg_slots": len(mg),
                "mapping": mapping,
                "injective": not collisions,
                "collisions": collisions,
                "unused_mg_slots": sorted(set(range(len(mg))) - set(mapping)),
            }
        )
    return results


# --------------------------------------------------------------------------
# Check 2: composition fidelity
# --------------------------------------------------------------------------
def run_case(name: str, metal: str, geo: str, ligand_keys: list[str]) -> dict[str, Any]:
    oin = build_oin(metal, geo, ligand_keys)
    exp = expected_content(ligand_keys)
    row: dict[str, Any] = {
        "case": name,
        "metal": metal,
        "geo": geo,
        "ligands": ligand_keys,
        "oin": oin,
        "expected": exp,
    }
    try:
        msmiles = convert_oin_to_msmiles(oin)
    except Exception as exc:
        row.update({"status": "RAISED", "error": f"{type(exc).__name__}: {exc}"})
        return row

    got = msmiles_content(msmiles)
    row["msmiles"] = msmiles
    row["got"] = got
    if "parse_error" in got:
        row.update({"status": "MSMILES_UNPARSEABLE"})
        return row

    violations = []
    if got["heavy"] != exp["heavy"]:
        violations.append("heavy_formula")
    if got["h"] != exp["h"]:
        violations.append("hydrogen_count")
    if got["ligand_charge"] != exp["ligand_charge"]:
        violations.append("ligand_charge")
    if got["n_ligands"] != exp["n_ligands"]:
        violations.append("ligand_dropped")
    if got["donors"] != exp["donors"]:
        violations.append("donor_set")
    if len(set(got["slots"])) != len(got["slots"]):
        violations.append("slot_collision")

    row["violations"] = violations
    row["status"] = "OK" if not violations else "VIOLATION"
    return row


def build_cases() -> list[tuple[str, str, str, list[str]]]:
    """Systematic sweep plus targeted adversarial probes."""
    cases: list[tuple[str, str, str, list[str]]] = []

    # Systematic: every geometry x every simple monodentate ligand, homoleptic.
    simple = ["chloride", "ammine", "aqua", "carbonyl", "cyanide", "hydroxo", "methyl"]
    for geo, cn in sorted(GEO_CN.items()):
        if geo not in OIN_TO_METALLOGEN_GEO:
            continue
        for lig in simple:
            cases.append((f"homoleptic/{geo}/{lig}", "Pt" if cn == 4 else "Fe", geo, [lig] * cn))

    # Mixed-ligand: exercises multiple slots with different donor elements.
    for geo, cn in sorted(GEO_CN.items()):
        if geo not in OIN_TO_METALLOGEN_GEO or cn < 4:
            continue
        mixed = [simple[i % len(simple)] for i in range(cn)]
        cases.append((f"mixed/{geo}", "Ru", geo, mixed))

    # Metal sweep at a fixed simple geometry.
    for metal in METALS:
        cases.append((f"metal/{metal}", metal, "OCT", ["chloride"] * 6))

    # Adversarial: the bare-donor reinterpretation branches.
    for lig in ["oxo_bare_O", "thiolate_bare_S", "nitride_bare_N", "amido_bare_N"]:
        cases.append((f"adversarial/{lig}/SPL", "Ti", "SPL", [lig] + ["chloride"] * 3))
        cases.append((f"adversarial/{lig}/OCT", "Zr", "OCT", [lig] + ["chloride"] * 5))

    # The documented ammine-vs-nitride ambiguity, stated both ways.
    cases.append(("ambiguity/ammine_reading", "Pt", "SPL", ["ammine"] * 4))
    cases.append(("ambiguity/nitride_reading", "Ti", "SPL", ["nitride_bare_N"] + ["chloride"] * 3))

    return cases


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, help="write full results as JSON")
    ap.add_argument("--verbose", action="store_true", help="print every case, not just failures")
    args = ap.parse_args()

    bijection = check_slot_bijection()
    non_injective = [b for b in bijection if b.get("injective") is False]

    rows = [run_case(*c) for c in build_cases()]
    by_status = Counter(r["status"] for r in rows)
    violation_kinds: Counter = Counter()
    for r in rows:
        for v in r.get("violations", []):
            violation_kinds[v] += 1

    print("=" * 78)
    print("C4 SYNTHETIC OIN BATTERY -- handoff tested in isolation")
    print("=" * 78)

    print(f"\n[1] SLOT BIJECTION over {len(bijection)} geometries")
    for b in bijection:
        if b.get("injective") is None:
            print(f"    {b['geo']:<5} UNMAPPED (no MetalloGen geometry)")
        else:
            flag = "OK " if b["injective"] else "COLLISION"
            extra = f"  collisions={b['collisions']}" if not b["injective"] else ""
            print(
                f"    {b['geo']:<5} {flag} oin_slots={b['n_oin_slots']} "
                f"mg_slots={b['n_mg_slots']}{extra}"
            )
    print(f"    -> non-injective geometries: {len(non_injective)}")

    print(f"\n[2] COMPOSITION FIDELITY over {len(rows)} synthetic complexes")
    for status, n in sorted(by_status.items()):
        print(f"    {status:<20} {n}")
    if violation_kinds:
        print("\n    violation kinds:")
        for k, n in violation_kinds.most_common():
            print(f"      {k:<20} {n}")

    bad = [r for r in rows if r["status"] != "OK"]
    if bad:
        print(f"\n    first {min(15, len(bad))} non-OK cases:")
        for r in bad[:15]:
            detail = r.get("error") or ",".join(r.get("violations", []))
            print(f"      {r['case']:<38} {r['status']:<14} {detail}")
            if args.verbose and r.get("msmiles"):
                print(f"        oin: {r['oin']}")
                print(f"        m  : {r['msmiles']}")
                print(f"        exp: {r['expected']}")
                print(f"        got: {r['got']}")

    if args.verbose:
        print("\n    all OK cases:")
        for r in rows:
            if r["status"] == "OK":
                print(f"      {r['case']:<38} {r['msmiles']}")

    result = {
        "slot_bijection": bijection,
        "non_injective_count": len(non_injective),
        "cases": rows,
        "status_counts": dict(by_status),
        "violation_kinds": dict(violation_kinds),
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2, default=str))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
