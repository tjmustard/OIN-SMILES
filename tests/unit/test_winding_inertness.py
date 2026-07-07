"""Inertness gate for ROADMAP-stereo.md Phase 1 (winding plumbing).

Phase 1 threads winding direction (`>`/`<`/`^`) through the parse layer
(`oin/inline.py` SlotAssignment) and the generation-side ParsedOIN
(`generation/oin_parser.py` OINVector.winding / winding_by_slot), but must not
change any 3D coordinate, placement, or bonding outcome — nothing in
`molassembler_adapter.py` reads the new fields yet.

These SHA-256 hashes were captured by running `OIN3DGenerator().generate()`
against `tests/candidate_outputs/*_oin.txt` on the pre-Phase-1 codebase (via
`git stash` of the three modified source files) and confirmed byte-identical
against the post-Phase-1 codebase before this test was written. If this test
ever fails, generation output has drifted — that is a real regression, not a
stale golden.
"""

import hashlib
import os
import unittest

from oinsmiles.generation.engine import OIN3DGenerator

_CANDIDATES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../candidate_outputs"))

# name -> sha256(emitted XYZ) captured on the pre-Phase-1 codebase.
_PRE_PHASE1_XYZ_HASHES = {
    "cisplatin_oin.txt": "841c476cf0f20e997e298bd520030233da403e19f7310ad0d1dbbe60ef44e31f",
    "transplatin_oin.txt": "6a9a5ca2907c5be50fcd214d481b3f49a98fd2e49758e96c0c8792468f79de99",
    "cis_ptcl2en_oin.txt": "5f742fc6735c7ef074a26e1d8bb4b5315af61b5d0537518df013d6906afe865f",
    "ferrocene_oin.txt": "b67ac60a4e6cd9fc0d9bca0380a59c986355443dc4aaf158312dea5d41f62ee2",
    "fac_irppy3_oin.txt": "1e654ec3e991e6255bc950c95837e382d3cc728f46b937e6a1d62f45381c7010",
    "mer_irppy3_oin.txt": "d48014ce1cd1df54e93ebd5af3a09cdf91be812a5d716558560dd56305a0baa1",
    "bdnn_oin.txt": "40becabc091d8a14f3b37aa60996f3bbe0daa15e26b2c935108f42699398a76f",
    "bdpp_oin.txt": "0181f27f472bcdcd5e38651c2540493dce6e1e4ca598fc593a56a8376e4121a5",
}


class TestWindingPlumbingInertness(unittest.TestCase):
    """Emitted XYZ must be byte-identical to pre-Phase-1 output for every fixture."""

    def test_generation_output_unchanged_by_winding_plumbing(self):
        # Pinned to the legacy Molassembler engine: these SHA-256 hashes are of the
        # deterministic legacy output. The metallogen default (now engine="metallogen")
        # is stochastic per FF embed and honors winding (it steers the eta face), so its
        # output is neither hash-stable nor winding-inert -- this gate is legacy-specific.
        generator = OIN3DGenerator(engine="legacy")
        mismatches = []

        for name, expected_hash in _PRE_PHASE1_XYZ_HASHES.items():
            path = os.path.join(_CANDIDATES_DIR, name)
            with open(path) as f:
                oin_string = f.read().strip()

            structure = generator.generate(oin_string)
            actual_hash = hashlib.sha256(structure.xyz.encode()).hexdigest()

            if actual_hash != expected_hash:
                mismatches.append(f"{name}: expected {expected_hash}, got {actual_hash}")

        self.assertEqual(
            mismatches,
            [],
            "Winding plumbing changed generation output (should be inert):\n"
            + "\n".join(mismatches),
        )


if __name__ == "__main__":
    unittest.main()
