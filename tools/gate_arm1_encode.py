"""v0.4.7 gate ARM 1: encoder byte-identity over the fixed `tests/fixtures/` set.

For each ``tests/fixtures/*.xyz`` file: encode and emit one line

    name<TAB>sha256(oin)<TAB>len<TAB>eta

sorted by name, then a ``# MANIFEST_SHA256=<...>`` line over the joined manifest text,
then a ``#DONE <n>`` sentinel.

WHY A COUNT, NOT A NAME LIST
=============================
The fixture set is asserted to be *exactly* ``EXPECTED_FIXTURE_COUNT`` files at runtime
(not hardcoded as a name list, unlike ``tools/enc_byte_identity_ab.py``'s curated
subset) -- a fixture added or removed is itself gate-relevant drift, so the count check
fails loudly rather than silently covering a different corpus than the frozen manifest
was built against.

That guard works, and v0.4.10 is the proof: it refused to run for a whole release cycle
after ``ULODUU_comp_0`` was added without a golden row. The cost of the refusal was that
a SECOND drift rode along unseen for just as long -- so when this fires, re-freeze
promptly and diff the pre-existing rows before you do.

MEMO-CLEARING DISCIPLINE (copied from ``tools/enc_byte_identity_ab.py``)
=========================================================================
``_ac2bo_memo_clear()`` runs BETWEEN molecules so a cross-molecule cache hit can
never be what makes two revisions agree.

ERRORS ARE PART OF THE CONTRACT
================================
A fixture that fails to encode still emits a line (``ERROR:<Type>:<msg>`` in the
sha column) and counts toward ``#DONE``. Two revisions must raise the SAME error,
not merely agree when both succeed -- a lane that turns a hard error into a silent
None would otherwise slip through unnoticed.

FLUSH DISCIPLINE
================
Every line is printed with ``flush=True`` immediately after being produced. Python
block-buffers stdout when it is not a tty (e.g. redirected to a file by the caller);
without an explicit flush, a `timeout` kill mid-run discards the buffer and a
downstream ``sort``/count still exits 0 on an empty or truncated file. Combined with
checking the ``#DONE <n>`` sentinel's ``n`` against the expected fixture count
before trusting any line above it, a truncated run cannot be mistaken for agreement.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/gate_arm1_encode.py \
        --fixtures-dir tests/fixtures > /tmp/arm1_out.tsv
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

HAPTIC = re.compile(r"\{\d+[<>]\}")

# v0.4.10: 61 -> 62. ULODUU_comp_0 landed as a fixture in dd51a515 (the boron TET
# correction) and the golden was never extended, so this guard has been hard-refusing
# ever since -- ARM 1 was simply not runnable through the whole of v0.4.9, which is how
# a second, unrelated drift (the v0.4.7 xyz2mol -> perception_tmc rename, visible in
# ASISAX_comp_0's ERROR text) also went unseen. Both are re-frozen in the golden; the
# other 60 rows were verified byte-identical first, so the encoder itself has not moved.
# Keep this count hardcoded rather than derived from the golden's row count: the two are
# asserted INDEPENDENTLY on purpose, and a fixture added without a matching golden row
# is exactly the drift this is here to catch.
EXPECTED_FIXTURE_COUNT = 62


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fixtures-dir", default="tests/fixtures")
    ap.add_argument(
        "--expect-n",
        type=int,
        default=EXPECTED_FIXTURE_COUNT,
        help=f"Fail loudly if the fixtures dir does not have exactly this many .xyz files "
        f"(default {EXPECTED_FIXTURE_COUNT})",
    )
    args = ap.parse_args()

    fixtures_dir = os.path.abspath(args.fixtures_dir)
    paths = sorted(glob.glob(os.path.join(fixtures_dir, "*.xyz")))
    if len(paths) == 0:
        sys.exit(f"error: 0 .xyz files found under {fixtures_dir} -- refusing an empty corpus")
    if len(paths) != args.expect_n:
        sys.exit(
            f"error: {fixtures_dir} has {len(paths)} .xyz files, expected exactly "
            f"{args.expect_n} -- the fixture set has drifted since the golden manifest "
            f"was frozen; update --expect-n deliberately if this is intentional"
        )

    from oinsmiles import XYZToSMILES
    from oinsmiles.utils import perception_core as loc

    clear = getattr(loc, "_ac2bo_memo_clear", None)

    lines = []
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        if clear is not None:
            clear()
        try:
            oin = XYZToSMILES().convert(path)
            sha = hashlib.sha256(oin.encode()).hexdigest()
            eta = "eta" if HAPTIC.search(oin) else "-"
            line = f"{name}\t{sha}\t{len(oin)}\t{eta}"
        except Exception as e:
            # An error is part of the contract too: it must be the SAME error.
            line = f"{name}\tERROR:{type(e).__name__}:{e}\t-\t-"
        lines.append(line)
        print(line, flush=True)

    manifest = "\n".join(lines)
    print(f"# molecules={len(lines)} fixtures_dir={fixtures_dir}", flush=True)
    print(f"# MANIFEST_SHA256={hashlib.sha256(manifest.encode()).hexdigest()}", flush=True)
    print(f"#DONE {len(lines)}", flush=True)


if __name__ == "__main__":
    main()
