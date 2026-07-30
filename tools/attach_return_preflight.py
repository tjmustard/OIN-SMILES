#!/usr/bin/env python3
"""v0.4.15 Lane 1 pre-flight: would an attachment-preferring RETURN path displace a passing conformer?

Lane 1 makes ``_select_by_geometry_impl``'s fallback return prefer a conformer whose claimed
coordination sites still hold an atom. That is a pure win only if the conformers it *demotes*
are not already round-tripping. The baseline sweep says 52 of the 3858 ``byte_exact`` molecules
carry a ``DETACHED`` verdict from ``oin.coordination.coordination_report`` -- so on those, an
attachment preference could reorder the pool away from the conformer that passes.

🔴 THE TWO PREDICATES ARE NOT THE SAME TEST, and that is the whole reason this probe exists.

* ``coordination_report`` (what ``attach_class_audit.py`` reports) compares the INPUT structure's
  donor set against the GENERATED structure's donor set. It answers "did the donor set change".
* ``attach_check.ligands_attached`` (what the guard uses) asks whether every coordination SITE
  THE OIN CLAIMS still holds at least one atom, measured on the generated geometry alone. It
  answers "did a site go empty".

A molecule can trip the first and pass the second: swapping which hydride the metal binds changes
the donor set without emptying a site. So "52 DETACHED byte_exact" is an upper bound on Lane 1's
exposure, not the exposure. This tool measures the guard's own verdict.

The claim reference is the **OIN's distinct slot count** -- coordination SITES, since an eta ring
writes one slot across every ring atom. In production the guard takes its claim from the contract
mol's metal bonds; the OIN is the same claim expressed in the artifact this run actually stored.
Bonds are never read as evidence of attachment: a detached ligand keeps its bond.

Usage:
    python tools/attach_return_preflight.py --sweep <results dir> --out-json out.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from oinsmiles.generation.attach_check import (  # noqa: E402
    encoder_donor_set,
    group_sites,
)

# Same slot alphabet as ``oin.compare._SLOT_RE`` -- ``^`` is a winding marker too, and a pattern
# that omits it silently under-counts the claim on haptic molecules, which are exactly the
# population Lane 1 targets.
_SLOT_RE = re.compile(r"\{(\d+)[><^]?\}")


def oin_claimed_sites(oin: str | None) -> int:
    """Distinct coordination SLOTS the OIN states, i.e. sites, not donor atoms."""
    if not oin:
        return 0
    return len(set(int(x) for x in _SLOT_RE.findall(oin)))


def read_xyz(path: str):
    """``(atomic_numbers, coords)`` in the file's own order -- which is ``res.xyz`` order.

    Order is load-bearing: it must match the generator mol's atom order. Same reader as
    ``tools/attach_probe.py``.
    """
    from rdkit.Chem import GetPeriodicTable

    pt = GetPeriodicTable()
    lines = open(path).read().splitlines()
    n = int(lines[0].strip())
    z, xyz = [], []
    for i in range(n):
        f = lines[2 + i].split()
        z.append(pt.GetAtomicNumber(f[0]))
        xyz.append([float(v) for v in f[1:4]])
    return z, np.array(xyz)


def score_one(mol: str, sweep: str):
    """The guard's verdict for one molecule's stored generated geometry."""
    xyz_path = os.path.join(sweep, "structures", f"{mol}_generated.xyz")
    oin_path = os.path.join(sweep, "structures", f"{mol}.oin")
    if not os.path.exists(xyz_path) or not os.path.exists(oin_path):
        return {"molecule": mol, "verdict": "NO_STRUCTURE"}
    try:
        znums, coords = read_xyz(xyz_path)
        oin = open(oin_path).read().strip()
        claimed = oin_claimed_sites(oin)
        metal_idx, actual = encoder_donor_set(znums, coords)
        if metal_idx is None:
            return {"molecule": mol, "verdict": "NO_METAL"}
        if claimed == 0:
            # Nothing claimed -> the guard abstains rather than rejects (see ligands_attached).
            return {"molecule": mol, "verdict": "ABSTAIN_NO_CLAIM"}
        actual_sites = len(group_sites(sorted(actual), coords)) if actual else 0
        # A claimed site that holds nothing is what the guard rejects on. With site COUNTS as
        # the observable, `actual < claimed` proves at least one claimed site is empty.
        verdict = "SITE_LOST" if actual_sites < claimed else "SITES_HELD"
        return {
            "molecule": mol,
            "verdict": verdict,
            "claimed_sites": claimed,
            "actual_sites": actual_sites,
            "actual_donors": len(actual),
        }
    except Exception as exc:  # a probe must say so, never abstain silently
        return {"molecule": mol, "verdict": "ERROR", "error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True, help="a results-* dir with structures/ + audit json")
    ap.add_argument("--out-json", required=True)
    ap.add_argument(
        "--control-n",
        type=int,
        default=250,
        help="size of the ATTACHED byte_exact control draw (0 disables)",
    )
    args = ap.parse_args()

    audit = json.load(open(os.path.join(args.sweep, "attach_class_audit.json")))
    table = audit["table"]

    pops = {
        # Lane 1's false-positive exposure: passing today, DETACHED by coordination_report.
        "byte_exact_DETACHED": table["byte_exact"]["DETACHED"],
        # Lane 1's target.
        "structural_DETACHED": table["structural"]["DETACHED"],
    }
    if args.control_n:
        # Control: byte_exact molecules the audit calls INTACT. If the guard says SITE_LOST here
        # it is firing on molecules with no attachment problem at all, and the probe -- or the
        # claim reference -- is wrong. A pre-flight without this arm cannot tell a real exposure
        # from a broken claim count.
        intact = table["byte_exact"]["INTACT"]
        step = max(1, len(intact) // args.control_n)
        pops["byte_exact_INTACT_control"] = intact[::step][: args.control_n]

    out = {"sweep": os.path.abspath(args.sweep), "populations": {}}
    for name, mols in pops.items():
        rows = [score_one(m, args.sweep) for m in mols]
        tally: dict[str, int] = {}
        for r in rows:
            tally[r["verdict"]] = tally.get(r["verdict"], 0) + 1
        out["populations"][name] = {"n": len(rows), "tally": tally, "rows": rows}
        print(f"\n=== {name} (n={len(rows)}) ===", flush=True)
        for k in sorted(tally, key=lambda k: -tally[k]):
            print(f"  {k:18s} {tally[k]:5d}  ({100 * tally[k] / len(rows):.1f}%)", flush=True)

    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {args.out_json}")

    be = out["populations"]["byte_exact_DETACHED"]
    st = out["populations"]["structural_DETACHED"]
    print(
        "\nLane 1 reading: exposure = SITE_LOST among byte_exact_DETACHED "
        f"({be['tally'].get('SITE_LOST', 0)}/{be['n']}); "
        f"target = SITE_LOST among structural_DETACHED "
        f"({st['tally'].get('SITE_LOST', 0)}/{st['n']})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
