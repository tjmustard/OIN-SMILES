"""Tests for Bug 3 fix: bare exit() in om.py replaced with RuntimeError.

Before the fix, MetalComplex.get_embedding() called exit() when no conformer
could be embedded.  exit() raises SystemExit which cannot be caught by
`except Exception`, so it killed the whole interpreter.

After the fix, RuntimeError is raised instead, which is a normal catchable
exception.
"""

import unittest
from unittest.mock import patch

from oinsmiles.generator3d import embed as embed_module
from oinsmiles.generator3d.om import MetalComplex


class TestGetEmbeddingRaisesNotExits(unittest.TestCase):
    def _make_empty_mc(self):
        return MetalComplex("", center_atom=None, ligands=[], chg=0, multiplicity=1)

    def test_raises_runtime_error_not_system_exit(self):
        """embed failures must raise RuntimeError, not SystemExit."""
        mc = self._make_empty_mc()
        with patch.object(embed_module, "get_embedding", return_value=None):
            with self.assertRaises(RuntimeError):
                mc.get_embedding(num_conformer=1)

    def test_system_exit_not_raised(self):
        """SystemExit must NOT be raised — it would kill the interpreter."""
        mc = self._make_empty_mc()
        with patch.object(embed_module, "get_embedding", return_value=None):
            try:
                mc.get_embedding(num_conformer=1)
            except SystemExit:
                self.fail("get_embedding() raised SystemExit — exit() is still present")
            except RuntimeError:
                pass  # correct behaviour

    def test_exception_is_catchable_by_except_exception(self):
        """The raised exception must be an instance of Exception."""
        mc = self._make_empty_mc()
        caught = False
        with patch.object(embed_module, "get_embedding", return_value=None):
            try:
                mc.get_embedding(num_conformer=1)
            except Exception:
                caught = True
        self.assertTrue(caught, "Exception was not caught by 'except Exception'")

    def test_error_message_is_descriptive(self):
        mc = self._make_empty_mc()
        with patch.object(embed_module, "get_embedding", return_value=None):
            with self.assertRaises(RuntimeError) as ctx:
                mc.get_embedding(num_conformer=1)
        self.assertIn("embedding", str(ctx.exception).lower())


if __name__ == "__main__":
    unittest.main()

