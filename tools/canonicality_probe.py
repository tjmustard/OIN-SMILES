"""Direct canonicality probe: does the encoder emit ONE string for one molecule?

WHY THIS EXISTS RATHER THAN A SWEEP DIFF
========================================
"Canonical" means: the same complex, presented differently, encodes to the same string.
This probe presents one input structure several ways and asserts the OIN string does not
move. Crucially, the transforms below **hold the molecular graph fixed**, so the expected
answer is not "similar" but **byte-identical** -- a known ground truth, which is what makes
this a sharp instrument.

    rotate    apply a random proper rotation to the coordinates
    renumber  permute the order atoms appear in the XYZ file
    both      renumber, then rotate

Any difference is a canonicality defect, and it is classified with the SAME subclass logic
the sweep bucket report uses (``slot_renumber`` / ``winding_star_drift`` /
``fragment_reorder`` / ``rdkit_canonical``), so a number here is directly comparable to a
``key_equal`` sub-split from a real sweep.

WHAT THIS FIXES ABOUT ``reencode_ab.py``
========================================
``reencode_ab.py`` compares the encoding of a stored INPUT structure against the encoding
of a stored GENERATED structure. That pair is two different geometries, which is the right
shape -- but the round-trip harness does not encode the generated structure the way that
tool must. The harness calls ``get_oin_string(gen_result.mol, coords)`` using the
GENERATOR'S OWN bond graph (``tools/test_dataset_roundtrip.py:185``); only on exception
does it fall back to ``convert()``. Re-encoding a stored generated ``.xyz`` from disk has
no access to ``gen_result.mol``, so it must re-perceive connectivity with
``xyz2AC_obabel`` -- and on a slightly distorted generated geometry that can yield a
genuinely different graph. Measured consequence: ``reencode_ab`` reports ``structural`` at
19.6% where the capstone sweep reports 1.0%. That excess is **connectivity-perception
drift, not serialization drift**, and it is out of scope for the canonicality lanes
(README section 8 names it as a separate perception-hardening problem).

So: use THIS probe as the primary instrument for canonicality lanes. Its transforms cannot
change the perceived graph, so nothing contaminates the signal. ``reencode_ab.py`` remains
useful for the ``key_equal`` sub-split, which tracks the sweep well (12.0% vs 12.3%), but
its ``structural`` and ``byte_exact`` numbers are NOT comparable to a sweep's.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --n 300 --trials 3
    PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --only ABAFUF_comp_0 -v
"""

import argparse
import contextlib
import glob
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.compare import (  # noqa: E402
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
)

DEFAULT_DATASET = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "tmCAT-tmPHOTO_xyz_dataset")
)
SEED = 42
_SLOT_RE = re.compile(r"\{(\d+)([<>^]?)\}")


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
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


def _subclass(s1, s2):
    """Same drift taxonomy as tools/roundtrip_bucket_report.py::_key_equal_subclass."""
    n1 = normalize_oin_for_comparison((s1 or "").strip())
    n2 = normalize_oin_for_comparison((s2 or "").strip())
    f1 = [f for f in n1.split(".") if f]
    f2 = [f for f in n2.split(".") if f]
    if sorted(f1) == sorted(f2) and f1 != f2:
        return "fragment_reorder"
    if _SLOT_RE.sub(r"{\2}", n1) == _SLOT_RE.sub(r"{\2}", n2):
        return "slot_renumber"
    strip_wind = re.compile(r"\{(\d+)[<>^]\}")
    if strip_wind.sub(r"{\1}", n1) == strip_wind.sub(r"{\1}", n2):
        return "winding_star_drift"
    return "rdkit_canonical"


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        p = line.split()
        syms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(coords, dtype=float), (lines[1].rstrip("\n") if len(lines) > 1 else "")


def write_xyz(path, syms, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, c in zip(syms, coords):
            f.write(f"{s:<3} {c[0]:>14.8f} {c[1]:>14.8f} {c[2]:>14.8f}\n")


def random_rotation(rng):
    """A uniformly random PROPER rotation (det = +1).

    Proper, not arbitrary orthogonal: an improper operation would mirror the structure,
    which legitimately changes a chiral molecule's encoding. Built by QR then forced to
    det = +1 -- flipping a column, never a reflection of the result.
    """
    a = np.asarray([[rng.gauss(0, 1) for _ in range(3)] for _ in range(3)], dtype=float)
    q, r = np.linalg.qr(a)
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def probe_one(path, conv, trials, rng, verbose=False):
    """Encode one structure under identity / rotate / renumber / both.

    Returns a per-molecule record. ``base`` failing to encode is reported, not raised --
    encoder coverage is a separate axis (the sweep's ``encode_fail``).
    """
    syms, coords, comment = read_xyz(path)
    rec = {
        "molecule": os.path.basename(path)[: -len(".xyz")],
        "natoms": len(syms),
        "base": None,
        "variants": [],
        "stable": None,
        "drift_subclasses": [],
        "key_stable": None,
    }
    try:
        with _silence_fds():
            rec["base"] = conv.convert(path)
    except Exception as e:  # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
        return rec
    if rec["base"] is None:
        rec["error"] = "encoder returned None"
        return rec

    tmpdir = tempfile.mkdtemp()
    try:
        for t in range(trials):
            for mode in ("rotate", "renumber", "both"):
                order = list(range(len(syms)))
                c = coords
                if mode in ("renumber", "both"):
                    rng.shuffle(order)
                    c = coords[order]
                s = [syms[i] for i in order]
                if mode in ("rotate", "both"):
                    c = c @ random_rotation(rng).T
                p = os.path.join(tmpdir, f"{mode}_{t}.xyz")
                write_xyz(p, s, c, comment)
                try:
                    with _silence_fds():
                        got = conv.convert(p)
                except Exception as e:  # noqa: BLE001
                    got = None
                    err = f"{type(e).__name__}: {e}"
                else:
                    err = None
                same = got == rec["base"]
                entry = {"mode": mode, "trial": t, "equal": same}
                if not same:
                    entry["subclass"] = _subclass(rec["base"], got) if got else "encode_fail"
                    entry["got"] = got
                    entry["error"] = err
                    try:
                        entry["key_equal"] = canonical_roundtrip_key(
                            rec["base"]
                        ) == canonical_roundtrip_key(got)
                    except Exception:  # noqa: BLE001
                        entry["key_equal"] = None
                rec["variants"].append(entry)
    finally:
        for fn in glob.glob(os.path.join(tmpdir, "*")):
            os.unlink(fn)
        os.rmdir(tmpdir)

    bad = [v for v in rec["variants"] if not v["equal"]]
    rec["stable"] = not bad
    rec["drift_subclasses"] = sorted({v["subclass"] for v in bad})
    rec["key_stable"] = all(v.get("key_equal") is not False for v in bad)
    if verbose and bad:
        print(
            f"  {rec['molecule']}: {len(bad)}/{len(rec['variants'])} drifted {rec['drift_subclasses']}"
        )
    return rec


def main():
    ap = argparse.ArgumentParser(description="Rotation/renumbering canonicality probe.")
    ap.add_argument("--n", type=int, default=300, help="Molecules to probe (0 = all)")
    ap.add_argument("--trials", type=int, default=2, help="Random trials per mode")
    ap.add_argument("--dataset", default=DEFAULT_DATASET)
    ap.add_argument("--subdirs", default="cat,photo")
    ap.add_argument("--only", help="Comma-separated molecule ids (bypasses sampling)")
    ap.add_argument("--out", help="Directory for the json/md report")
    ap.add_argument("--shard", help="1-based shard spec I:N")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    files = []
    for sub in args.subdirs.split(","):
        root = os.path.join(os.path.abspath(args.dataset), sub.strip())
        if os.path.isdir(root):
            for dirpath, _d, fns in sorted(os.walk(root)):
                files += [
                    os.path.join(dirpath, f)
                    for f in sorted(fns)
                    if f.endswith(".xyz") and not f.endswith("_generated.xyz")
                ]
    # Dedup by basename (1,033 names exist in both cat/ and photo/), first wins.
    seen, uniq = set(), []
    for p in files:
        b = os.path.basename(p)
        if b not in seen:
            seen.add(b)
            uniq.append(p)
    files = uniq

    if args.only:
        want = {m.strip().removesuffix(".xyz") for m in args.only.split(",") if m.strip()}
        files = [p for p in files if os.path.basename(p)[: -len(".xyz")] in want]
    elif args.n and args.n > 0:
        files = sorted(random.Random(SEED).sample(files, min(args.n, len(files))))
    if args.shard:
        i, n = (int(x) for x in args.shard.split(":"))
        files = files[i - 1 :: n]

    conv = XYZToSMILES()
    rng = random.Random(SEED)
    recs = []
    t0 = time.time()
    for i, path in enumerate(files, 1):
        recs.append(probe_one(path, conv, args.trials, rng, args.verbose))
        if i % 25 == 0:
            print(f"  {i}/{len(files)}  {i / max(time.time() - t0, 1e-9):.1f} mol/s", flush=True)

    encoded = [r for r in recs if r.get("base")]
    stable = [r for r in encoded if r["stable"]]
    drift = [r for r in encoded if not r["stable"]]
    sub = Counter(s for r in drift for s in r["drift_subclasses"])
    by_mode = Counter(v["mode"] for r in drift for v in r["variants"] if not v["equal"])
    key_broken = [r for r in drift if not r["key_stable"]]

    print(f"\n{'=' * 66}\nCANONICALITY PROBE  ({len(recs)} probed, {args.trials} trials/mode)")
    print(f"{'=' * 66}")
    print(f"  encoder failed on input : {len(recs) - len(encoded)}")
    print(
        f"  byte-stable             : {len(stable)}/{len(encoded)}"
        f"  ({100 * len(stable) / max(len(encoded), 1):.2f}%)"
    )
    print(
        f"  DRIFTED                 : {len(drift)}/{len(encoded)}"
        f"  ({100 * len(drift) / max(len(encoded), 1):.2f}%)"
    )
    print("\n  drift by subclass (a molecule can show more than one):")
    for k, v in sub.most_common():
        print(f"    {k:22} {v}")
    print("\n  drift by transform:")
    for k, v in by_mode.most_common():
        print(f"    {k:22} {v}")
    print(f"\n  of the drifted, key ALSO changed (isomer-level, worse): {len(key_broken)}")
    if key_broken:
        print("    " + ", ".join(r["molecule"] for r in key_broken[:15]))

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        try:
            commit = (
                subprocess.check_output(
                    ["git", "rev-parse", "--short", "HEAD"],
                    cwd=os.path.dirname(os.path.abspath(__file__)),
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .strip()
            )
        except Exception:
            commit = "unknown"
        payload = {
            "probed_at": datetime.now().isoformat(timespec="seconds"),
            "commit_id": commit,
            "levers": {k: v for k, v in sorted(os.environ.items()) if k.startswith("OIN_")},
            "trials_per_mode": args.trials,
            "seed": SEED,
            "n_probed": len(recs),
            "n_encoded": len(encoded),
            "n_stable": len(stable),
            "n_drifted": len(drift),
            "drift_by_subclass": dict(sub),
            "drift_by_transform": dict(by_mode),
            "key_broken": [r["molecule"] for r in key_broken],
            "records": recs,
        }
        tag = (args.shard or "all").replace(":", "of")
        with open(os.path.join(args.out, f"canonicality_probe_{tag}.json"), "w") as f:
            json.dump(payload, f, indent=2)
        print(f"\nWrote {args.out}/canonicality_probe_{tag}.json")


if __name__ == "__main__":
    main()
