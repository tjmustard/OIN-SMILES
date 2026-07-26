"""Does a perception lever change WHICH coordination geometry the encoder reports?

The veto this exists to enforce
===============================
``OIN_STABLE_METAL_AC`` makes the AC valence-capping loop order-independent by capping the
heaviest atom first. That is a **perception** change, and it is asymmetric: capping the metal
before the light atoms bridging to it can only **keep** metal bonds that the old
atom-order-dependent iteration discarded. More metal bonds means a higher coordination
number, which means the geometric template fit can land on a *different polyhedron*.

So the lever can silently rewrite the `[M_XXX]` tag — and the tag is not cosmetic: it selects
the vertex table, hence the rotation group, hence the canonical slot labelling and the
comparison key's whole vertex signature. A lever that improves byte-stability while quietly
reclassifying geometries across the corpus would be a bad trade made invisibly.

Fixing one molecule cleanly (`DUDREA_comp_0`, `[Y_SPY]` → `[Y_TET]` under renumbering) is not
evidence of safety at scale. This measures the scale.

What it reports
===============
Per molecule, one encode with the lever off and one with it on, then:

* the `[M_XXX]` transition matrix (off-tag → on-tag), so a shift is visible by direction
  rather than only as a count;
* coordination number, taken as the number of distinct `{n}` slots, so a tag change can be
  attributed to the metal genuinely gaining a donor rather than to a template tie flipping;
* whether the change is confined to the tag or the whole string moved.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/geometry_tag_shift.py --n 200 \\
        --lever OIN_STABLE_METAL_AC [--shard 1:4] [--out DIR]
"""

import argparse
import contextlib
import json
import os
import random
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles import XYZToSMILES  # noqa: E402

DEFAULT_DATASET = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "tmCAT-tmPHOTO_xyz_dataset")
)
SEED = 42
_METAL_RE = re.compile(r"\[([A-Z][a-z]?)(?:@[A-Z0-9]+)?_([A-Z]{3})\]")
_SLOT_RE = re.compile(r"\{(\d+)[<>^]?\}")


@contextlib.contextmanager
def _silence_fds():
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _geo(oin):
    m = _METAL_RE.search(oin or "")
    return m.group(2) if m else None


def _cn(oin):
    """Coordination number as the count of DISTINCT slot indices (eta rings share one)."""
    return len({int(x) for x in _SLOT_RE.findall(oin or "")})


def main():
    ap = argparse.ArgumentParser(description="Geometry-tag distribution shift under a lever.")
    ap.add_argument("--lever", default="OIN_STABLE_METAL_AC")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--dataset", default=DEFAULT_DATASET)
    ap.add_argument("--subdirs", default="cat,photo")
    ap.add_argument("--shard", help="1-based I:N")
    ap.add_argument("--out")
    args = ap.parse_args()

    files, seen = [], set()
    for sub in args.subdirs.split(","):
        root = os.path.join(os.path.abspath(args.dataset), sub.strip())
        if not os.path.isdir(root):
            continue
        for dirpath, _d, fns in sorted(os.walk(root)):
            for f in sorted(fns):
                if f.endswith(".xyz") and not f.endswith("_generated.xyz") and f not in seen:
                    seen.add(f)
                    files.append(os.path.join(dirpath, f))
    if not files:
        sys.exit(
            f"error: no .xyz under {os.path.abspath(args.dataset)}/{{{args.subdirs}}} — the "
            "dataset is gitignored, so a worktree lacks it; pass --dataset explicitly"
        )
    files = sorted(random.Random(SEED).sample(files, min(args.n, len(files))))
    if args.shard:
        i, n = (int(x) for x in args.shard.split(":"))
        files = files[i - 1 :: n]

    conv = XYZToSMILES()
    trans = Counter()
    cn_delta = Counter()
    rows = []
    n_ok = n_tag_moved = n_str_moved = 0

    for pos, path in enumerate(files, 1):
        name = os.path.basename(path)[: -len(".xyz")]
        out = {}
        for label, val in (("off", None), ("on", "1")):
            prior = os.environ.get(args.lever)
            if val is None:
                os.environ.pop(args.lever, None)
            else:
                os.environ[args.lever] = val
            try:
                with _silence_fds():
                    out[label] = conv.convert(path)
            except Exception as e:  # noqa: BLE001
                out[label] = None
                out[label + "_err"] = f"{type(e).__name__}: {e}"
            finally:
                if prior is None:
                    os.environ.pop(args.lever, None)
                else:
                    os.environ[args.lever] = prior

        if out.get("off") is None or out.get("on") is None:
            rows.append({"molecule": name, "skipped": "encode failed in one arm"})
            continue
        n_ok += 1
        g_off, g_on = _geo(out["off"]), _geo(out["on"])
        c_off, c_on = _cn(out["off"]), _cn(out["on"])
        if g_off != g_on:
            n_tag_moved += 1
            trans[(g_off, g_on)] += 1
        if out["off"] != out["on"]:
            n_str_moved += 1
        if c_off != c_on:
            cn_delta[(c_off, c_on)] += 1
        if g_off != g_on or c_off != c_on:
            rows.append(
                {
                    "molecule": name,
                    "geo_off": g_off,
                    "geo_on": g_on,
                    "cn_off": c_off,
                    "cn_on": c_on,
                    "oin_off": out["off"],
                    "oin_on": out["on"],
                }
            )
        if pos % 25 == 0:
            print(f"  {pos}/{len(files)}", flush=True)

    print(f"\n{'=' * 68}\nGEOMETRY-TAG SHIFT under {args.lever}   ({n_ok} molecules encoded)")
    print("=" * 68)
    print(f"  string changed at all : {n_str_moved}/{n_ok}")
    print(f"  [M_XXX] tag CHANGED   : {n_tag_moved}/{n_ok}   <-- the veto number")
    if trans:
        print("\n  transitions (off -> on):")
        for (a, b), c in trans.most_common():
            print(f"    {a} -> {b}   {c}")
    if cn_delta:
        print("\n  coordination-number changes (distinct slots, off -> on):")
        for (a, b), c in cn_delta.most_common():
            print(f"    {a} -> {b}   {c}")
        print("  A tag change WITH a CN change is the lever keeping a real metal bond.")
        print("  A tag change WITHOUT one is a template tie flipping -- less defensible.")
    else:
        print("\n  no coordination-number changes at all")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        tag = (args.shard or "all").replace(":", "of")
        with open(os.path.join(args.out, f"geometry_tag_shift_{tag}.json"), "w") as f:
            json.dump(
                {
                    "lever": args.lever,
                    "n_encoded": n_ok,
                    "n_string_changed": n_str_moved,
                    "n_tag_changed": n_tag_moved,
                    "transitions": {f"{a}->{b}": c for (a, b), c in trans.items()},
                    "cn_changes": {f"{a}->{b}": c for (a, b), c in cn_delta.items()},
                    "rows": rows,
                },
                f,
                indent=2,
            )
        print(f"\nWrote {args.out}/geometry_tag_shift_{tag}.json")


if __name__ == "__main__":
    main()
