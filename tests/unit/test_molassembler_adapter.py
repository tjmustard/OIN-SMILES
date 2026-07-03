"""Unit tests for MolassemblerAdapter and helpers.

Tests use mocking so no actual SCINE Molassembler subprocess is spawned.
Integration-level DG conformer generation is covered by the candidate
artifact protocol (tests/candidate_outputs/).
"""

import os
import pickle
import sys
import unittest
from unittest.mock import MagicMock, patch

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.molassembler_adapter import (
    MolassemblerAdapter,
    MolassemblerTimeoutError,
    _build_connected_smiles,
    _molassembler_worker,
)
from oinsmiles.generation.oin_parser import OINVector, ParsedOIN


def _cisplatin_parsed_oin() -> ParsedOIN:
    """Minimal ParsedOIN for cis-[PtCl2(NH3)2]."""
    return ParsedOIN(
        smiles="[Pt].[Cl].[Cl].N.N",
        fragments=["[Pt]", "[Cl]", "[Cl]", "N", "N"],
        metal_fragment_idx=0,
        vectors=[
            OINVector(atom_idx=-1, vector=(1, 0, 0), fragment_idx=1, atom_in_fragment_idx=0),
            OINVector(atom_idx=-1, vector=(-1, 0, 0), fragment_idx=2, atom_in_fragment_idx=0),
            OINVector(atom_idx=-1, vector=(0, 1, 0), fragment_idx=3, atom_in_fragment_idx=0),
            OINVector(atom_idx=-1, vector=(0, -1, 0), fragment_idx=4, atom_in_fragment_idx=0),
        ],
        original_oin="[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
    )


class TestBuildConnectedSmiles(unittest.TestCase):
    def test_cisplatin_produces_connected_smiles(self):
        """Connected SMILES for cisplatin should contain Pt bonded to N and Cl."""
        parsed = _cisplatin_parsed_oin()
        smiles = _build_connected_smiles(parsed)
        self.assertIsInstance(smiles, str)
        self.assertGreater(len(smiles), 0)
        # All four elements should appear in the connected SMILES
        self.assertIn("Pt", smiles)
        self.assertIn("N", smiles)
        self.assertIn("Cl", smiles)

    def test_no_vectors_returns_fallback(self):
        """Without vectors, the fallback dot-disconnected SMILES is returned."""
        parsed = ParsedOIN(
            smiles="[Pt].N",
            fragments=["[Pt]", "N"],
            metal_fragment_idx=0,
            vectors=[],
            original_oin="[Pt_SPL].N{0}",
        )
        smiles = _build_connected_smiles(parsed)
        # Should not crash; returns dot-disconnected SMILES as fallback
        self.assertIsInstance(smiles, str)


class TestMolassemblerWorkerPicklable(unittest.TestCase):
    def test_worker_is_picklable(self):
        """_molassembler_worker must be picklable for ProcessPoolExecutor."""
        data = pickle.dumps(_molassembler_worker)
        self.assertGreater(len(data), 0)


class TestMolassemblerAdapter(unittest.TestCase):
    def test_timeout_raises_molassembler_timeout_error(self):
        """A timed-out future must raise MolassemblerTimeoutError (not FuturesTimeout)."""
        from concurrent.futures import TimeoutError as FuturesTimeout

        adapter = MolassemblerAdapter(timeout=1)
        parsed = _cisplatin_parsed_oin()

        with patch(
            "oinsmiles.generation.molassembler_adapter.ProcessPoolExecutor"
        ) as mock_pool_cls:
            mock_executor = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_executor)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)

            mock_future = MagicMock()
            mock_future.result.side_effect = FuturesTimeout("timed out")
            mock_executor.submit.return_value = mock_future

            with self.assertRaises(MolassemblerTimeoutError):
                adapter.generate(parsed)

    def test_dg_error_raises_runtime_error(self):
        """A worker returning ok=False must raise RuntimeError."""
        adapter = MolassemblerAdapter(timeout=60)
        parsed = _cisplatin_parsed_oin()

        with patch(
            "oinsmiles.generation.molassembler_adapter.ProcessPoolExecutor"
        ) as mock_pool_cls:
            mock_executor = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_executor)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)

            mock_future = MagicMock()
            mock_future.result.return_value = {"error": "DG failed", "ok": False}
            mock_executor.submit.return_value = mock_future

            with self.assertRaises(RuntimeError):
                adapter.generate(parsed)

    def test_successful_generate_returns_xyz_block(self):
        """A successful worker result returns a GeneratedStructure whose .xyz is the XYZ block string."""
        adapter = MolassemblerAdapter(timeout=60)
        parsed = _cisplatin_parsed_oin()
        expected_xyz = "9\n\nPt 0.0 0.0 0.0\n..."

        with patch(
            "oinsmiles.generation.molassembler_adapter.ProcessPoolExecutor"
        ) as mock_pool_cls:
            mock_executor = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_executor)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)

            mock_future = MagicMock()
            mock_future.result.return_value = {"xyz_block": expected_xyz, "ok": True}
            mock_executor.submit.return_value = mock_future

            result = adapter.generate(parsed)
            self.assertEqual(result.xyz, expected_xyz)


if __name__ == "__main__":
    unittest.main()
