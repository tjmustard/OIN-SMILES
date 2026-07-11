"""Inertness gate for ROADMAP-stereo.md Phase 1 (winding plumbing).

Phase 1 threaded winding direction (`>`/`<`/`^`) through the parse layer
(`oin/inline.py` SlotAssignment) and the generation-side ParsedOIN
(`generation/oin_parser.py` OINVector.winding / winding_by_slot), but must not
change any 3D coordinate, placement, or bonding outcome on the legacy engine --
nothing in `molassembler_adapter.py` reads the new fields. This gate proves that.

The reference values were captured by running `OIN3DGenerator(engine="legacy")`
against `tests/candidate_outputs/*_oin.txt` on the pre-Phase-1 codebase (via
`git stash` of the three modified source files) and confirmed unchanged against the
post-Phase-1 codebase before this test was written.

**Why two gates, not one hash set.** An exact SHA-256 of the emitted XYZ is only a
reliable golden for *rigid* complexes, whose distance-geometry embedding converges
to a single geometry regardless of the runner's BLAS/LAPACK build. A *flexible*
chelate (the en `-CH2CH2-` pucker in cis-PtCl2(en); the bdnn/bdpp backbones) lands
in float-different coordinates across platforms, so its byte hash drifts on the CI
runner even though the molecule built is identical -- a false "regression." Those
fixtures are therefore gated on the platform-stable invariant the plumbing could
actually break: atom count and element multiset. (The metallogen default engine is
stochastic per FF embed and *does* honor winding, so this whole gate is
legacy-specific.)
"""

import hashlib
import os
import unittest
from collections import Counter

from oinsmiles.generation.engine import OIN3DGenerator

_CANDIDATES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../candidate_outputs"))

# Rigid complexes: exact emitted-XYZ hash is a stable, meaningful golden (square-planar
# PtCl2 cores, the ferrocene sandwich, the rigid octahedral tris-chelate Ir complexes).
_RIGID_XYZ_HASHES = {
    "cisplatin_oin.txt": "841c476cf0f20e997e298bd520030233da403e19f7310ad0d1dbbe60ef44e31f",
    "transplatin_oin.txt": "6a9a5ca2907c5be50fcd214d481b3f49a98fd2e49758e96c0c8792468f79de99",
    "ferrocene_oin.txt": "b67ac60a4e6cd9fc0d9bca0380a59c986355443dc4aaf158312dea5d41f62ee2",
    "fac_irppy3_oin.txt": "1e654ec3e991e6255bc950c95837e382d3cc728f46b937e6a1d62f45381c7010",
    "mer_irppy3_oin.txt": "d48014ce1cd1df54e93ebd5af3a09cdf91be812a5d716558560dd56305a0baa1",
}

# Flexible chelates: byte hash is platform-fragile (conformer float-drift), so gate the
# invariant instead -- (atom count, element -> count). Captured on the same reference run.
_FLEXIBLE_FIXTURE_COMPOSITION = {
    "cis_ptcl2en_oin.txt": (15, {"C": 2, "Cl": 2, "H": 8, "N": 2, "Pt": 1}),
    "bdnn_oin.txt": (64, {"C": 29, "Cl": 2, "H": 30, "N": 2, "Pd": 1}),
    "bdpp_oin.txt": (64, {"C": 29, "Cl": 2, "H": 30, "P": 2, "Pd": 1}),
}


def _xyz_composition(xyz):
    """(atom_count, {element: count}) parsed from an XYZ string."""
    lines = xyz.strip().splitlines()
    natoms = int(lines[0].strip())
    counts = Counter(line.split()[0] for line in lines[2 : 2 + natoms])
    return natoms, dict(counts)


class TestWindingPlumbingInertness(unittest.TestCase):
    """Legacy generation output must be unchanged by the winding plumbing."""

    def setUp(self):
        # Pinned to the legacy Molassembler engine: the reference values are its
        # deterministic (rigid) / composition-stable (flexible) output.
        self.generator = OIN3DGenerator(engine="legacy")

    def _generate(self, name):
        with open(os.path.join(_CANDIDATES_DIR, name)) as f:
            return self.generator.generate(f.read().strip()).xyz

    def test_rigid_generation_output_unchanged(self):
        """Rigid fixtures must be byte-identical to the pre-Phase-1 emitted XYZ."""
        mismatches = []
        for name, expected_hash in _RIGID_XYZ_HASHES.items():
            actual_hash = hashlib.sha256(self._generate(name).encode()).hexdigest()
            if actual_hash != expected_hash:
                mismatches.append(f"{name}: expected {expected_hash}, got {actual_hash}")
        self.assertEqual(
            mismatches,
            [],
            "Winding plumbing changed rigid generation output (should be inert):\n"
            + "\n".join(mismatches),
        )

    def test_flexible_generation_composition_stable(self):
        """Flexible chelates: atom count + element multiset must be unchanged.

        The exact coordinates drift across platforms (conformer float-noise), so the
        byte hash is not a portable golden; the composition is what the plumbing could
        actually corrupt (a dropped/duplicated atom, a stray dummy) and is stable.
        """
        mismatches = []
        for name, (exp_n, exp_counts) in _FLEXIBLE_FIXTURE_COMPOSITION.items():
            n, counts = _xyz_composition(self._generate(name))
            if (n, counts) != (exp_n, exp_counts):
                mismatches.append(f"{name}: expected ({exp_n}, {exp_counts}), got ({n}, {counts})")
        self.assertEqual(
            mismatches,
            [],
            "Winding plumbing changed flexible generation composition (should be inert):\n"
            + "\n".join(mismatches),
        )


if __name__ == "__main__":
    unittest.main()
