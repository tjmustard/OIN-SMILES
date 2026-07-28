# v0.4.13 authoritative table — DERIVED, not swept

Derived from `results-v0.4.8-honest` by `tools/fold_transition_sim.py --arm veto`.
**The generator was NOT re-run** — see `FROZEN.md` for why that is exact here.

## Buckets

| bucket | count | % |
|---|---:|---:|
| byte_exact | 3794 | 75.88% |
| key_equal | 439 | 8.78% |
| facmer_divergent | 16 | 0.32% |
| structural | 417 | 8.34% |
| hard_fail | 319 | 6.38% |
| encode_fail | 15 | 0.30% |
| **total** | **5000** | **100.00%** |

### key_equal sub-split

| subclass | count |
|---|---:|
| rdkit_canonical | 114 |
| slot_renumber | 325 |

## Movement vs the v0.4.8 baseline

| | |
|---|---|
| `byte_exact` before | 3623 (72.46%) |
| `byte_exact` after | **3794 (75.88%)** |
| points | **+3.42** |
| moved in a BAD direction | **0** |
| excluded (drift / unavailable) | **0 / 0** |

- `key_equal/slot_renumber -> byte_exact` — **171**
