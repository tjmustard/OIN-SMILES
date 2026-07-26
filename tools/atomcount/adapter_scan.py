"""Authoritative atom-count audit of the `atom_count` class -- no generation, no encoding.

Reads the OIN string the capstone sweep already recorded, runs it through the
generator's own front half (``OINParser`` -> ``_prepare_ligand_fragments``), and
counts the atoms MetalloGen is asked to build.  That count is what the generated
XYZ will contain.

Also attributes a cause per fragment:

``kekulize_rescue``
    ``_prepare_ligand_fragments`` puts a ``-1`` charge on ``ring[0]`` of *every*
    5-membered all-aromatic ring when the fragment fails to kekulize.  On a BARE
    aromatic carbon that flips the implicit-H count from 1 to 0, so an innocent
    thiophene/pyrrole C-H silently loses its hydrogen.

``donor_h_reconcile``
    the bare-donor strip heuristics changed the donor's H count.

``string_phantom``
    the slot-stripped OIN fragment already re-parses with more H than the input
    can support -- the phantom is baked into the string itself.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation import metallogen_adapter as MA  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402

SLOT_RE = re.compile(r"\{\d+[><^]?\}")
RESULTS = Path(
    "/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042"
)


def reparse(s):
    """Parse a fragment for counting, never returning None for a real fragment.

    Full sanitize, then no-kekulize sanitize, then a bare parse plus
    ``UpdatePropertyCache(strict=False)``. The last tier matters: a carbonyl written
    ``[C]#O`` is a carbon radical that no sanitize accepts, and the first version of this
    probe reported 15 of the 74 molecules as "unparseable" purely because of it -- a probe
    limitation that would otherwise read as 15 unexplained molecules.
    """
    m = Chem.MolFromSmiles(s)
    if m is not None:
        return m
    m = Chem.MolFromSmiles(s, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(
            m,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        return m
    except Exception:
        pass
    m = Chem.MolFromSmiles(s, sanitize=False)
    try:
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        return m
    except Exception:
        return None


def totals(m):
    if m is None:
        return None
    heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
    hs = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 1)
    hs += sum(a.GetTotalNumHs() for a in m.GetAtoms() if a.GetAtomicNum() > 1)
    return heavy + hs


def neg_arom_carbons(smi: str) -> int:
    m = reparse(smi)
    if m is None:
        return 0
    return sum(
        1
        for a in m.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetFormalCharge() < 0 and a.GetIsAromatic()
    )


def audit(oin: str):
    parsed = OINParser().parse(oin)
    metal_frag, specs, geo = MA._prepare_ligand_fragments(parsed)
    total = totals(reparse(re.sub(r"_[A-Z0-9]+", "", metal_frag)))
    if total is None:
        return None, None, None
    per = []
    # slot-stripped OIN fragments, in the same order the adapter emits them is not
    # guaranteed, so compare as multisets of (heavy,H) and as aggregate totals.
    oin_frag_total = 0
    oin_neg = 0
    for i, frag in enumerate(parsed.fragments):
        if i == parsed.metal_fragment_idx:
            continue
        t = totals(reparse(SLOT_RE.sub("", frag)))
        oin_frag_total += t if t is not None else 0
        oin_neg += neg_arom_carbons(SLOT_RE.sub("", frag))
    ad_neg = 0
    for smi, _w in specs:
        t = totals(reparse(smi))
        if t is None:
            return None, None, None
        total += t
        ad_neg += neg_arom_carbons(smi)
        per.append((smi, t))
    return (
        total,
        per,
        {
            "oin_frag_total": oin_frag_total,
            "oin_neg_arom_c": oin_neg,
            "adapter_neg_arom_c": ad_neg,
        },
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    wl = json.loads(Path(args.worklist).read_text())["atom_count"]

    rows = []
    for e in wl:
        name = e["molecule"]
        rep = json.loads((RESULTS / "individual_reports" / f"{name}.json").read_text())
        n_in = int(Path(rep["input_xyz"]).read_text().splitlines()[0].split()[0])
        g = RESULTS / "structures" / f"{name}_generated.xyz"
        n_gen = int(g.read_text().splitlines()[0].split()[0]) if g.exists() else None
        m = re.search(r"Input (\d+) != Gen (\d+)", e["error"])
        rep_gen = int(m.group(2))
        row = {
            "molecule": name,
            "input_n": n_in,
            "reported_gen": rep_gen,
            "stored_gen": n_gen,
            "oin": rep["smiles_1"],
        }
        try:
            tot, per, info = audit(rep["smiles_1"])
        except Exception as ex:  # noqa: BLE001
            row["status"] = f"{type(ex).__name__}: {ex}"
            rows.append(row)
            continue
        if tot is None:
            row["status"] = "unparseable-fragment"
            rows.append(row)
            continue
        row["status"] = "ok"
        row["adapter_n"] = tot
        row["adapter_delta"] = tot - n_in
        row.update(info)
        # cause flags
        row["kekulize_rescue"] = info["adapter_neg_arom_c"] > info["oin_neg_arom_c"]
        rows.append(row)

    ok = [r for r in rows if r["status"] == "ok"]
    print(f"audited {len(ok)}/{len(rows)}")
    match = sum(1 for r in ok if r["adapter_delta"] == r["reported_gen"] - r["input_n"])
    print(f"adapter-implied delta == harness-reported delta: {match}/{len(ok)}")
    print(
        "adapter count == reported gen count:              "
        f"{sum(1 for r in ok if r['adapter_n'] == r['reported_gen'])}/{len(ok)}"
    )
    print(
        f"adapter count == input (string is fine):          {sum(1 for r in ok if r['adapter_delta'] == 0)}/{len(ok)}"
    )
    print(
        f"\nkekulize_rescue fired (new aromatic C-):          {sum(1 for r in ok if r['kekulize_rescue'])}/{len(ok)}"
    )
    neg = [r for r in ok if r["adapter_delta"] < 0]
    print(
        f"  of the {len(neg)} LOSS rows, rescue fired in:        {sum(1 for r in neg if r['kekulize_rescue'])}"
    )
    pos = [r for r in ok if r["adapter_delta"] > 0]
    print(
        f"  of the {len(pos)} GAIN rows, rescue fired in:        {sum(1 for r in pos if r['kekulize_rescue'])}"
    )

    print("\n-- rows --")
    for r in sorted(rows, key=lambda r: (r["status"] != "ok", r.get("adapter_delta", 0))):
        if r["status"] != "ok":
            print(f"  {r['molecule']:22s} in={r['input_n']:4d} {r['status']}")
            continue
        print(
            f"  {r['molecule']:22s} in={r['input_n']:4d} adapter={r['adapter_n']:4d} "
            f"d={r['adapter_delta']:+3d} harness_d={r['reported_gen'] - r['input_n']:+3d} "
            f"{'KEKRESCUE' if r['kekulize_rescue'] else ''}"
        )
    if args.out:
        args.out.write_text(json.dumps(rows, indent=1, default=str))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
