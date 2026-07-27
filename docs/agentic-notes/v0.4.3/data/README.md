# Phase 2 telemetry data — v0.4.3 elimination study

Raw inputs backing §8b of `docs/agentic-notes/v0.4.3/FALSIFICATION_v0.4.3_ELIMINATION.md`.

| file | what | produced by |
|---|---|---|
| `telemetry_events.json` | per-molecule silent-degradation events for the 861-molecule stratified run (non-quick, sequential, PASS-1 generator config) | `tools/telemetry_run.py` |
| `telemetry_strata.json` | the four strata: S1 cleanest passes / S2 most-distorted passes / S3 bucket-E gate failures / S4 buckets B+D | built from the distortion MPO + failure buckets |

Regenerate the report tables:

```
uv run python tools/telemetry_analyze.py
```
