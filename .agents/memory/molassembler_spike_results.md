# Molassembler Spike Results

✅ Import `scine_molassembler` successful.
   Version: 3.0.1
   Import path: `scine_molassembler`

✅ `_test_worker` is picklable (40 bytes).
✅ `_cisplatin_conformer` is picklable (48 bytes).

✅ ProcessPoolExecutor minimal test passed: 'Worker output {'test': 123}'

✅ Cisplatin conformer generated successfully via Molassembler.
API usage:
- masm.Molecule()
- mol.add_atom(utils.ElementType.X)
- mol.add_bond(a, b, masm.BondType.Single)
- mol.set_shape_at_atom(pt, masm.shapes.Shape.Square)
- masm.dg.generate_conformation(mol, seed, conf)

XYZ written to tests/candidate_outputs/spike_cisplatin.xyz