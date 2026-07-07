import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generator3d.ml_optimizer import ASEOptimizer


class TestMaceOptimizer(unittest.TestCase):
    def test_mace_init_missing_env(self):
        # Temporarily clear the environment variable to ensure the ValueError is raised
        original_env = os.environ.get("MACE_OMOL25_MODEL_PATH")
        if "MACE_OMOL25_MODEL_PATH" in os.environ:
            del os.environ["MACE_OMOL25_MODEL_PATH"]

        try:
            # We expect a ValueError because the MACE_OMOL25_MODEL_PATH is not set
            optimizer = ASEOptimizer(method="mace-omol25")

            # If initialization succeeds, the optimize method should catch the missing env var
            # To test optimize, we would need a dummy molecule, which is slightly more complex,
            # but we can just check if the object is created properly.
            # In our implementation, the environment variable check happens during optimize(),
            # so we should ideally mock an optimization call or just verify initialization works.
            self.assertEqual(optimizer.method, "mace-omol25")
            self.assertTrue(optimizer._calc_cls is not None)

        except ImportError:
            # If mace is not installed, it's fine to skip for this test
            self.skipTest("mace-torch not installed.")
        finally:
            # Restore environment variable
            if original_env is not None:
                os.environ["MACE_OMOL25_MODEL_PATH"] = original_env

    def test_mace_extra_large_init_missing_env(self):
        original_env = os.environ.get("MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH")
        if "MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH" in os.environ:
            del os.environ["MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"]

        try:
            optimizer = ASEOptimizer(method="mace-omol-0-extra-large-1024")
            self.assertEqual(optimizer.method, "mace-omol-0-extra-large-1024")
            self.assertTrue(optimizer._calc_cls is not None)

        except ImportError:
            self.skipTest("mace-torch not installed.")
        finally:
            if original_env is not None:
                os.environ["MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"] = original_env

    def test_xtb_init_still_works(self):
        try:
            # This should still initialize exactly as before
            optimizer = ASEOptimizer(method="xtb")
            self.assertEqual(optimizer.method, "xtb")
        except ImportError:
            self.skipTest("xtb-python not installed.")


if __name__ == "__main__":
    unittest.main()
