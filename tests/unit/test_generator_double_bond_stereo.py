"""Guard that the MetalloGen 3D generator preserves C=C (cis/trans, E/Z) geometry.

The vendored generator's molecule graph ("ace_mol") is stereo-blind, and the
dummy-metal embed re-perceives bond orders (dropping the alkene C=C to a single
bond), so before this fix the generated double-bond dihedral was random with a
cis bias -- an E-input could embed as Z and vice versa. cis/trans is kinetically
locked, so the generator must reproduce the requested geometry deterministically.

The C=C stereo is threaded from the ligand SMILES (`/C=C/`) through the ace_mol
and enforced on the embed rd_mol (and again on the per-ligand FF mol in the
cleanup), so distance geometry embeds the correct E/Z. These tests feed an E and
a Z alkenyl-amine OIN to the generator, re-encode the generated 3D structure with
stereo perception ON, and assert the perceived C=C stereo matches -- and that E
and Z do not collapse to the same geometry.

Note: `XYZToSMILES.convert` runs the encoder with ``with_stereo=False`` (no 3D
stereo perception), so these tests re-encode via ``get_tmc_mol(with_stereo=True)``
+ ``get_oin_string`` to read the geometry the generator actually produced.
"""

import os
import re
import sys
import tempfile
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol

OIN_E = "[Pt_SPL].C/C=C/CCN{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}"
OIN_Z = r"[Pt_SPL].C/C=C\CCN{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}"


def _perceive_alkene_stereo(xyz):
    """Re-encode a generated XYZ (stereo perception on) and read its C=C stereo.

    Returns ``STEREONONE`` if xyz2mol's re-perception left the double bond's
    stereo undetermined for this particular geometry.
    """
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as tmp_file:
            tmp_file.write(xyz)
            tmp_path = tmp_file.name
        mol, coords = get_tmc_mol(tmp_path, 0, with_stereo=True)
        oin2 = get_oin_string(mol, coords)
    finally:
        if tmp_path is not None:
            os.unlink(tmp_path)
    for frag in oin2.split("."):
        clean = re.sub(r"\{[^}]*\}", "", frag)
        if "C=C" not in clean:
            continue
        m = Chem.MolFromSmiles(clean)
        if m is None:
            continue
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
        for b in m.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and b.GetStereo() != Chem.BondStereo.STEREONONE
            ):
                return b.GetStereo()
    return Chem.BondStereo.STEREONONE


def _generated_alkene_stereo(oin_string):
    """Generate a 3D structure for ``oin_string`` and return its perceived C=C stereo.

    This used to retry up to six times: the MetalloGen embed drew its seed from an
    unseeded ``random.randint``, so generation was sampled rather than
    deterministic and a run could round-trip to an undetermined double bond. That
    made this test the flaky "task #8" case in the worklogs. The embed is now
    seeded (``generate_3d_structures(seed=42)``), so one generation is the whole
    story -- a retry would return the identical structure.
    """
    gen = OIN3DGenerator(engine="metallogen", optimizer="ff")
    return _perceive_alkene_stereo(gen.generate(oin_string).xyz)


class TestGeneratorDoubleBondStereo(unittest.TestCase):
    def test_e_input_generates_e_geometry(self):
        self.assertEqual(_generated_alkene_stereo(OIN_E), Chem.BondStereo.STEREOE)

    def test_z_input_generates_z_geometry(self):
        self.assertEqual(_generated_alkene_stereo(OIN_Z), Chem.BondStereo.STEREOZ)

    def test_e_and_z_do_not_collapse(self):
        """The core guarantee: E and Z inputs must not generate the same geometry."""
        e = _generated_alkene_stereo(OIN_E)
        z = _generated_alkene_stereo(OIN_Z)
        self.assertNotEqual(e, Chem.BondStereo.STEREONONE)
        self.assertNotEqual(z, Chem.BondStereo.STEREONONE)
        self.assertNotEqual(e, z, "E and Z inputs generated the same C=C geometry")

    def test_generation_is_reproducible(self):
        """No retry is needed because the seeded embed returns the same structure."""
        gen = OIN3DGenerator(engine="metallogen", optimizer="ff")
        runs = {gen.generate(OIN_E).xyz for _ in range(3)}
        self.assertEqual(len(runs), 1, "the seeded embed must be reproducible")


class TestDoubleBondStereoDonorFilter(unittest.TestCase):
    """Only geometrically-free double bonds carry an enforced stereo constraint.

    A double bond lying on a chelate ring -- a cycle closed through the metal --
    is already held rigid by the coordination; adding a distance-geometry stereo
    constraint there is unnecessary and over-constrains the metal embed
    (measured: AGULIX conformer yield 9/9 -> 3/9). The ligand builder must drop
    those, while keeping a pendant, freely-rotatable alkene. See
    ``test_chelate_locked_ez.py`` for the predicate itself.
    """

    def test_pendant_alkene_is_enforced(self):
        from oinsmiles.generator3d.ligand import get_ligand_from_smiles

        lig = get_ligand_from_smiles("C/C=C/CC[NH2:1]")
        self.assertTrue(
            lig.molecule.stereo_bonds,
            "a pendant alkene (donor 3 bonds away) should keep its E/Z constraint",
        )

    def test_metal_bound_imine_is_skipped(self):
        from oinsmiles.generator3d.ligand import get_ligand_from_smiles

        # AGULIX ligand: the C=N sits on the S->N chelate path.
        lig = get_ligand_from_smiles(r"CS/C(=N/[N:3]=Cc1ccccc1[O:6])[S:5]")
        self.assertEqual(
            lig.molecule.stereo_bonds,
            [],
            "a C=N on a chelate ring must not be enforced",
        )


if __name__ == "__main__":
    unittest.main()
