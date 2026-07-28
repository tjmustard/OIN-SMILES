# results-v0.4.8-honest — frozen authoritative honest baseline

Re-scored twin of `results-v0.4.6-sweep` (see `SOURCE`). The generator was NOT re-run: every
`smiles_2_indep` is an independent `XYZToSMILES().convert()` of the *stored* generated XYZ, which
is bit-identical to what `OIN_INDEP_SCORE` computes in the live harness.

| file | what |
|---|---|
| `bucket_report_PASS1_authoritative.{md,json}` | **the frozen record** — scored vs honest side by side, plus the transition matrix (copy of `bucket_report_both.*`) |
| `bucket_report_honest.{md,json}` | honest-only table |
| `honest_rescore.jsonl` | per-molecule transitions, `#DONE 5000` sentinel |
| `encoder_identity.jsonl` | corpus encoder byte-identity gate, `#DONE 4985` |
| `atom_count_provenance.json` | Lane 2's per-atom hydrogen provenance for the 18 |
| `individual_reports/` | the v0.4.6 reports + `smiles_2_indep` / `indep_key_match` / `honest_class` / backfilled `coordination` |
| `structures/` | symlink to the source sweep — geometry is NOT duplicated |

Headline: `byte_exact` 4140 (82.80%) -> **3623 (72.46%)**. Full analysis in
`docs/agentic-notes/v0.4.8/HONEST_BASELINE_v0.4.8.md`.
