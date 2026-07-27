#!/usr/bin/env python
"""Time GENERATION for a boron-cage molecule -- the arm the boron lane never measured.

``docs/BORON_CAGE_v0.4.5.md`` §5's 34/34 table is nine ENCODER-side checks (encodes,
byte-identical repeat encode, fragments re-parse, atom/bond multisets, key computable and
stable). None of them invokes the 3D generator, so that document's "encode and round-trip"
does not cover OIN -> XYZ -> OIN.

``XIQKOY_comp_0`` shows why it matters: with ``OIN_BORON_CAGE=0`` the encoder emits the cage as
a disconnected fragment and generation raises ``UncoordinatedFragmentError`` in 0.01 s (total
0.87 s); with the lever ON the encode is a correct coordinated B10 cage and generation runs past
340 s. This bounds how much of the 34 behaves that way -- see §10 of that document for results.

Prints one JSON object per run so a shell loop can accumulate JSONL.

Usage:
    GEN_CAP=60 python tools/boron_gen_time.py <MOLNAME>     # e.g. AVOFIB
"""

from __future__ import annotations

import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.generation.metallogen_adapter import OIN3DGeneratorMetallogen  # noqa: E402

#: Dataset root, overridable so this is not pinned to one checkout.
DATASET = os.environ.get(
    "OIN_DATASET_DIR",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "tmCAT-tmPHOTO_xyz_dataset",
    ),
)


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: boron_gen_time.py <MOLNAME>", file=sys.stderr)
        return 2
    name = sys.argv[1]
    cap = float(os.environ.get("GEN_CAP", "60"))
    out: dict = {"mol": name, "cap_s": cap}

    hits = sorted(glob.glob(f"{DATASET}/cat/{name}*.xyz")) + sorted(
        glob.glob(f"{DATASET}/photo/{name}*.xyz")
    )
    if not hits:
        out["error"] = f"no xyz found under {DATASET}"
        print(json.dumps(out))
        return 1

    out["xyz"] = hits[0]
    t0 = time.monotonic()
    try:
        oin = XYZToSMILES().convert(hits[0])
        out["encode_s"] = round(time.monotonic() - t0, 2)
        out["oin"] = oin
        t1 = time.monotonic()
        res = OIN3DGeneratorMetallogen(
            optimizer=None, ensemble_size=1, timeout=cap, ff_params=None
        ).generate(oin)
        # NOTE: this routinely EXCEEDS `cap`. embed_time_budget=timeout bounds the embed attempt
        # loop, not the OIN-direct assembly around it, so only the harness's per-molecule SIGKILL
        # subprocess really enforces a budget. Measured overruns: 60s asked, 60.7-137.9s spent.
        out["gen_s"] = round(time.monotonic() - t1, 2)
        out["got_mol"] = getattr(res, "mol", None) is not None
    except Exception as e:
        out["gen_s"] = round(time.monotonic() - t0 - out.get("encode_s", 0.0), 2)
        out["got_mol"] = False
        out["error"] = f"{type(e).__name__}: {str(e)[:110]}"
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
