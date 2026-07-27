#!/usr/bin/env python
"""Profile ONE molecule's encode and 3D-generation, replicating the harness UFF_1 tier.

The eta runtime question ("can a job finish in <30 s?") has been answered twice by
reasoning about which stage `accept_fn` sees, and both answers were imprecise. This
measures instead: it splits wall-clock into encode vs generate, then reports the
cumulative-time leaders inside generate so the cost has a named owner.

Kwargs are copied from `tools/test_dataset_roundtrip.py`'s PASS 1 UFF fast-pass
(`optimizer=None, ensemble_size=1, timeout=300, ff_params=None`) so the numbers are
comparable to a sweep row rather than to a synthetic configuration.

Usage:
    python tools/profile_eta.py <mol.xyz> [--top 30] [--sort cumtime] [--out prof.txt]
"""

from __future__ import annotations

import argparse
import cProfile
import io
import os
import pstats
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.generation.metallogen_adapter import (  # noqa: E402
    OIN3DGeneratorMetallogen as OIN3DGenerator,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--sort", default="cumtime")
    ap.add_argument("--out")
    ap.add_argument(
        "--timeout",
        type=float,
        default=300.0,
        help="generator timeout; match the sweep budget being explained",
    )
    args = ap.parse_args()

    print(f"=== {os.path.basename(args.xyz)} ===", flush=True)

    t0 = time.monotonic()
    oin = XYZToSMILES().convert(args.xyz)
    t_encode = time.monotonic() - t0
    print(f"encode : {t_encode:8.2f}s", flush=True)
    print(f"OIN    : {oin}", flush=True)

    gen = OIN3DGenerator(optimizer=None, ensemble_size=1, timeout=args.timeout, ff_params=None)

    prof = cProfile.Profile()
    t1 = time.monotonic()
    prof.enable()
    try:
        result = gen.generate(oin)
    except Exception as e:  # a failure still has a profile worth reading
        result = None
        print(f"generate RAISED: {type(e).__name__}: {e}", flush=True)
    finally:
        prof.disable()
    t_gen = time.monotonic() - t1

    print(f"generate: {t_gen:8.2f}s   (result={'ok' if result else 'None'})", flush=True)
    print(f"TOTAL   : {t_encode + t_gen:8.2f}s   encode={100 * t_encode / (t_encode + t_gen):.1f}%")

    buf = io.StringIO()
    pstats.Stats(prof, stream=buf).sort_stats(args.sort).print_stats(args.top)
    text = buf.getvalue()
    print(text, flush=True)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(f"{os.path.basename(args.xyz)} encode={t_encode:.2f}s gen={t_gen:.2f}s\n")
            fh.write(f"OIN: {oin}\n\n")
            fh.write(text)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
