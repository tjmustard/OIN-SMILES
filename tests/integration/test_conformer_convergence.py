"""True-conformer convergence of the canonical OIN-SMILES encoding.

A canonical encoding must be conformer-invariant: genuinely *different* 3D geometries of the
*same* isomer (a CREST conformer ensemble -- not merely one geometry rigidly reoriented) must
encode to the *same* canonical key
``winding_canonical_key(normalize_oin_for_comparison(XYZToSMILES().convert(frame)))``.

This is the true-conformer companion to ``test_conformer_invariance.py``. There the gate is
byte-identical raw strings (one geometry, rigidly rotated/translated); here the frames are
distinct minima, so they carry benign rotational slot drift and the gate is the canonical
**key**, not the raw string (that drift is exactly what the key is designed to absorb).

Fixtures are curated, frame-capped CREST ensembles under
``tests/fixtures/conformer_ensembles/`` (see that dir's README; built by
``scratchpad/build_ensemble_fixtures.py`` from the gitignored sweep). Two notable full-tier
fixtures:

* ``BEPCAC_comp_0`` -- a regression guard for the B1 electronic geometry prior (without it this
  Ni d8 complex splits SPL/TET across conformers and this test fails).
* ``CETDAI_comp_0`` (verdict ``notation-drift``) -- its frames span >=2 distinct *raw* OIN
  strings that share *one* canonical key; the test positively asserts the key collapses that
  drift.

Run the fast subset (two small ensembles)::

    uv run python -m unittest tests.integration.test_conformer_convergence

Run the full set::

    OIN_CONVERGENCE_FULL=1 uv run python -m unittest tests.integration.test_conformer_convergence
"""

import contextlib
import json
import os
import tempfile
import unittest
from pathlib import Path

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "tests" / "fixtures" / "conformer_ensembles" / "manifest.json"

FULL = os.environ.get("OIN_CONVERGENCE_FULL") == "1"


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings).

    Only wraps the conversion calls; Python exceptions propagate as objects, so real
    failures are unaffected.
    """
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


def _iter_frame_blocks(path):
    """Yield each frame of a concatenated multi-frame XYZ as a standalone XYZ text block."""
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].strip())
        block = "\n".join(lines[i : i + 2 + n]) + "\n"
        yield block
        i += 2 + n


def _key(oin):
    return winding_canonical_key(normalize_oin_for_comparison(oin))


def _load_ensembles():
    if not MANIFEST.exists():
        return None
    return json.loads(MANIFEST.read_text())["ensembles"]


class ConformerConvergenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ensembles = _load_ensembles()
        if cls.ensembles is None:
            raise unittest.SkipTest(
                f"missing {MANIFEST.relative_to(REPO)}; regenerate with "
                "scratchpad/build_ensemble_fixtures.py"
            )
        cls.converter = XYZToSMILES()

    def _encode_frames(self, ensemble):
        """Return (raws, keys) -- one per frame -- for an ensemble fixture."""
        path = REPO / ensemble["path"]
        self.assertTrue(path.exists(), f"fixture missing: {ensemble['path']}")
        raws, keys = [], []
        with tempfile.TemporaryDirectory() as td:
            for idx, block in enumerate(_iter_frame_blocks(path)):
                fp = os.path.join(td, f"{ensemble['molecule']}_{idx}.xyz")
                Path(fp).write_text(block)
                with _silence_fds():
                    oin = self.converter.convert(fp)
                raws.append(oin)
                keys.append(_key(oin))
        return raws, keys

    def _check(self, ensembles):
        """Every frame of each ensemble must encode to one canonical key."""
        failures = []
        for e in ensembles:
            raws, keys = self._encode_frames(e)
            self.assertGreaterEqual(
                len(keys), 2, f"{e['molecule']}: need >=2 frames to test convergence"
            )
            base = keys[0]
            # NB: _key returns (str, list) -- unhashable -- so compare by equality, not a set.
            for idx, k in enumerate(keys[1:], start=1):
                if k != base:
                    failures.append(
                        f"{e['molecule']} (metal={e['metal']}, verdict={e['verdict']}) "
                        f"frame {idx}: canonical key diverged\n"
                        f"    frame 0: {raws[0]}\n    frame {idx}: {raws[idx]}"
                    )
            # Notation-drift fixtures positively demonstrate the key absorbing raw drift:
            # >=2 distinct RAW strings, still one canonical key.
            if e["verdict"] == "notation-drift":
                distinct_raw = len(set(raws))
                if distinct_raw < 2:
                    failures.append(
                        f"{e['molecule']}: expected >=2 distinct raw OIN encodings to "
                        f"demonstrate notation drift, got {distinct_raw} "
                        f"(all identical: {raws[0]})"
                    )
        if failures:
            self.fail(f"{len(failures)} convergence failure(s):\n" + "\n".join(failures))

    def test_convergence_fast(self):
        """The two small fast-tier ensembles converge to a single canonical key."""
        fast = [e for e in self.ensembles if e["tier"] == "fast"]
        self.assertTrue(fast, "no fast-tier ensembles in manifest")
        self._check(fast)

    @unittest.skipUnless(FULL, "set OIN_CONVERGENCE_FULL=1 to run the full ensemble set")
    def test_convergence_full(self):
        """Every curated ensemble converges (incl. BEPCAC B1-guard and CETDAI drift demo)."""
        self._check(self.ensembles)


if __name__ == "__main__":
    unittest.main()
