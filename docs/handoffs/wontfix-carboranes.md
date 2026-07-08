# Documented limitation: carborane / 3c2e borane clusters (no session this wave)

~14 registry rows (`carborane_unsupported`, e.g. `HIMQIF_comp_0`, `AFOGEK_comp_0`)
fail in the FORWARD encode: `get_lig_mol` cannot build a valid RDKit template
for ligand fragments like `[H]B1[B-]2([H])[B-]3([H])...` — polyhedral boranes
whose 3-center-2-electron bonding has no faithful 2-electron SMILES valence
model.

This is a **notation-design problem, not a bug**: OIN would need an explicit
cluster convention (e.g. treating the cage as an eta-like multi-atom unit or a
pseudo-atom) before these can round-trip. Deferred until a design exists.

Interim expectation: the encoder should fail with a clear
"polyhedral borane ligands unsupported" error rather than a raw RDKit valence
traceback (S3 may route the error message, but the capability itself is out of
scope for this wave). Keep the rows in the registry under `wontfix-docs` so
they are not mistaken for regressions.
