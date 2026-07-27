"""Encode-side byte-identity A/B: sha256 of the OIN string for a fixed molecule set.

The encoder's output *is* the product, so any speed change must be SHA-identical. Run
this on both revisions and diff the output:

    git show HEAD:src/oinsmiles/utils/perception_core.py > /tmp/base.py
    trap 'cp /tmp/mine.py src/oinsmiles/utils/perception_core.py' EXIT   # NEVER git stash
    ...

Prints one ``name<TAB>sha256<TAB>len<TAB>eta`` line per molecule, sorted, plus a sha256
over the whole manifest so a single string can be compared. ``eta`` flags a haptic OIN
(``{N>}`` / ``{N<}``) so the set is visibly not just cisplatin.

Each molecule is encoded in a **fresh subprocess-free but cache-cleared** state so a
cross-molecule memo cannot mask a difference.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

HAPTIC = re.compile(r"\{\d+[<>]\}")

# 4 goldens + a stratified dataset sample. The dataset entries deliberately include the
# two molecules whose encode this lane is about (QIDKUL eta / QIDKIZ non-eta), several
# haptic complexes, and ordinary fast ones -- a change that only preserved cisplatin
# would be caught here.
DEFAULT_FIXTURES = [
    "CisPlatin",
    "TransPlatin",
    "Ferrocene",
    "POJJOP",
    "Cis-PtCl2(en)",
    "fac-Ir(ppy)3",
    "mer-Ir(ppy)3",
    "PdCl2-R-BINAP",
    "Zeises_salt",
    "YESKOZ",
    "FeCO5",
    "CuCN2",
]

DEFAULT_DATASET = [
    "QIDKUL_comp_0",
    "QIDKIZ_comp_0",
    "ABERUW_comp_0",
    "ABEXOU_comp_0",
    "ABIFUM_comp_0",
    "ACASOO_comp_0",
    "ACAXIO_comp_0",
    "ACINIM_comp_0",
    "ADAMAT_comp_0",
    "ADANEB_comp_0",
    "AFADOC_comp_0",
    "AGUFEN_comp_0",
]


def resolve(name: str, fixtures_dir: str, dataset_dir: str) -> str | None:
    cands = [os.path.join(fixtures_dir, f"{name}.xyz")]
    for sub in ("regression_inputs", "cat", "photo"):
        cands.append(os.path.join(dataset_dir, sub, f"{name}.xyz"))
    for c in cands:
        if os.path.exists(c):
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures-dir", default="tests/fixtures")
    ap.add_argument("--dataset", default="tmCAT-tmPHOTO_xyz_dataset")
    ap.add_argument("--only", default=None, help="comma-separated subset")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    names = (
        [n.strip() for n in args.only.split(",") if n.strip()]
        if args.only
        else DEFAULT_FIXTURES + DEFAULT_DATASET
    )

    from oinsmiles import XYZToSMILES
    from oinsmiles.utils import perception_core as loc

    lines = []
    missing = []
    for name in sorted(names):
        path = resolve(name, args.fixtures_dir, args.dataset)
        if path is None:
            missing.append(name)
            continue
        # Clear the AC2BO memo between molecules: a cross-molecule cache hit must never
        # be what makes two revisions agree.
        clear = getattr(loc, "_ac2bo_memo_clear", None)
        if clear is not None:
            clear()
        t0 = time.perf_counter()
        try:
            oin = XYZToSMILES().convert(path)
            sha = hashlib.sha256(oin.encode()).hexdigest()
            eta = "eta" if HAPTIC.search(oin) else "-"
            lines.append(f"{name}\t{sha}\t{len(oin)}\t{eta}")
        except Exception as e:
            # An error is part of the contract too: it must be the SAME error.
            lines.append(f"{name}\tERROR:{type(e).__name__}:{e}\t-\t-")
        print(f"# {name} {time.perf_counter() - t0:.2f}s", file=sys.stderr, flush=True)

    manifest = "\n".join(lines)
    print(manifest)
    print(f"\n# molecules={len(lines)} missing={missing}")
    print(f"# MANIFEST_SHA256={hashlib.sha256(manifest.encode()).hexdigest()}  {args.label}")


if __name__ == "__main__":
    main()
