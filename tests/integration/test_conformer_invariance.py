"""Conformer / orientation invariance of the canonical OIN-SMILES encoding.

A canonical encoding must be orientation-invariant: the same 3D structure, rigidly
rotated and translated (atom order preserved), must encode to the *same* OIN-SMILES
string. This exercises the PAI orientation canonicalization in
``oinsmiles.utils.perception_tmc._align_to_pai`` across the curated size-stratified set in
``tests/fixtures/conformer_set/`` (built by ``tools/select_conformer_test_set.py``)
plus the CisPlatin / TransPlatin anchors.

Primary gate is **byte-identical** raw strings. As a diagnostic, we also compare
``winding_canonical_key(normalize_oin_for_comparison(...))``: if the raw strings
differ but the canonical keys match, the message flags benign notation drift rather
than a hard canonicalization failure.

This is a rigid-geometry guard; it is NOT a substitute for the true-conformer check
(different geometries of the same molecule) in ``tools/conformer_invariance_crest.py``,
which needs the external CREST binary.

Run the fast subset (2 fixtures + a few small structures)::

    uv run python -m unittest tests.integration.test_conformer_invariance

Run the full ~30-structure set::

    OIN_CONFORMER_FULL=1 uv run python -m unittest tests.integration.test_conformer_invariance
"""

import contextlib
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "tests" / "fixtures" / "conformer_set" / "manifest.json"

FULL = os.environ.get("OIN_CONFORMER_FULL") == "1"
FAST_EXTRA = 3  # smallest dataset structures added to the two fixtures in fast mode


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


def _read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    elements, coords = [], []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return comment, elements, np.asarray(coords, dtype=float)


def _write_xyz(path, comment, elements, coords):
    with open(path, "w") as f:
        f.write(f"{len(elements)}\n{comment}\n")
        for el, c in zip(elements, coords):
            f.write(f"{el} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n")


def _random_rotation(rng):
    """Uniform random proper rotation via QR of a Gaussian matrix."""
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def _key(oin):
    return winding_canonical_key(normalize_oin_for_comparison(oin))


def _load_structures():
    if not MANIFEST.exists():
        return None
    data = json.loads(MANIFEST.read_text())
    return sorted(data["structures"], key=lambda s: (s["n_heavy"], s["molecule"]))


class ConformerInvarianceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structures = _load_structures()
        if cls.structures is None:
            raise unittest.SkipTest(
                f"missing {MANIFEST.relative_to(REPO)}; regenerate with "
                "tools/select_conformer_test_set.py"
            )
        cls.converter = XYZToSMILES()

    def _check(self, structures, n_rotations):
        """Rotate/translate each structure n_rotations times; require identical OIN."""
        failures = []
        with tempfile.TemporaryDirectory() as td:
            for s in structures:
                path = REPO / s["path"]
                self.assertTrue(path.exists(), f"fixture missing: {s['path']}")
                comment, elements, coords = _read_xyz(path)
                # Deterministic per-structure RNG so runs are reproducible
                # (builtin hash() is per-process randomized; use a stable digest).
                digest = hashlib.sha1(s["molecule"].encode()).digest()
                seed = 0xC0FFEE ^ int.from_bytes(digest[:4], "big")
                rng = np.random.default_rng(seed)
                with _silence_fds():
                    base = self.converter.convert(str(path))
                base_key = _key(base)
                for k in range(n_rotations):
                    rot = _random_rotation(rng)
                    trans = rng.uniform(-30.0, 30.0, size=3)
                    moved = coords @ rot.T + trans
                    tp = os.path.join(td, f"{s['molecule']}_{k}.xyz")
                    _write_xyz(tp, comment, elements, moved)
                    with _silence_fds():
                        oin = self.converter.convert(tp)
                    if oin != base:
                        drift = (
                            "KEY-MATCH (notation drift)" if _key(oin) == base_key else "KEY-DIFFERS"
                        )
                        failures.append(
                            f"{s['molecule']} (metal={s['metal']}, heavy={s['n_heavy']}) "
                            f"rotation {k}: {drift}\n    base: {base}\n    got : {oin}"
                        )
        if failures:
            self.fail(
                f"{len(failures)} orientation-invariance divergence(s):\n" + "\n".join(failures)
            )

    def test_rotation_invariance_fast(self):
        """Two Pt fixtures + the smallest few dataset structures encode identically."""
        fixtures = [s for s in self.structures if s["source"] == "fixture"]
        dataset = [s for s in self.structures if s["source"] != "fixture"][:FAST_EXTRA]
        self._check(fixtures + dataset, n_rotations=3)

    @unittest.skipUnless(FULL, "set OIN_CONFORMER_FULL=1 to run the full ~30-structure set")
    def test_rotation_invariance_full(self):
        """Every structure in the curated set encodes identically under rotation."""
        self._check(self.structures, n_rotations=4)


if __name__ == "__main__":
    unittest.main()
