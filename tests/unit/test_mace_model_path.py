"""Guard tests for the MACE weights lookup (R1 audit OSPO-REVIEWER-11).

The MACE optimizer used to call ``load_dotenv()``, which copies every line of a ``.env``
file into ``os.environ``. The ``OIN_*`` levers that change keys are read from
``os.environ`` on each call, so a stray ``.env`` could change keys for the rest of the
process. ``resolve_mace_model_path`` reads only the one weights-path variable instead.
"""

import ast
import importlib.util
import inspect
import os
import tempfile
import unittest
from unittest import mock

from oinsmiles.generator3d import ml_optimizer
from oinsmiles.generator3d.ml_optimizer import resolve_mace_model_path

VAR = "MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"


def _write_env(dirpath, lines):
    path = os.path.join(dirpath, ".env")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


class TestMaceModelPath(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("dotenv"), "python-dotenv is not installed")
    def test_reads_only_its_variable_and_leaves_environ_alone(self):
        with tempfile.TemporaryDirectory() as d, mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(VAR, None)
            os.environ.pop("OIN_DOTENV_PROBE", None)
            os.environ.pop("OIN_XTB_BIN", None)
            env = _write_env(
                d,
                [
                    f'{VAR}="models/mace/w.model"',
                    "OIN_DOTENV_PROBE=0",
                    "OIN_XTB_BIN=/bin/false",
                ],
            )
            before = dict(os.environ)
            self.assertEqual(resolve_mace_model_path(VAR, dotenv_path=env), "models/mace/w.model")
            self.assertEqual(dict(os.environ), before)

    def test_process_environment_wins(self):
        with tempfile.TemporaryDirectory() as d, mock.patch.dict(os.environ, {VAR: "/from/env"}):
            env = _write_env(d, [f"{VAR}=/from/dotenv"])
            self.assertEqual(resolve_mace_model_path(VAR, dotenv_path=env), "/from/env")

    def test_unset_everywhere_is_none(self):
        with tempfile.TemporaryDirectory() as d, mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(VAR, None)
            env = _write_env(d, ["OTHER=1"])
            self.assertIsNone(resolve_mace_model_path(VAR, dotenv_path=env))

    def test_optimizer_never_loads_dotenv_into_environ(self):
        tree = ast.parse(inspect.getsource(ml_optimizer))
        called = {
            getattr(node.func, "id", getattr(node.func, "attr", None))
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
        }
        self.assertNotIn("load_dotenv", called)


if __name__ == "__main__":
    unittest.main()
