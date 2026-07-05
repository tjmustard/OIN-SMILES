# TASK-40: per-complex round-trip filter + failure artifacts (WS-0)

Status: DONE (2026-07-04) — `--only` filter + failure-artifact dump landed;
`--only Ferrocene`=1, `--only TiCat`=4, full run 19/6 with identical failing
set, TiCat3 crash leaves `Ex*_error.txt`. Suites green. See NOTES.md Log.
Depends on: none
Suggested model: Haiku (mechanical, fully specified)
Effort: ~0.25 session · Risk: none (harness-only)

Part of the η-ligand round-trip recovery effort — see
`spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (§3, WS-0). Read that and
`spec/worklog/NOTES.md` first. Baseline (2026-07-04): round-trip 19/25,
6 η-ligand failures (Ferrocene, TiCp2Me2, TiCat1–4).

## Goal

Make `verify_roundtrip.py` iterable on a single complex without the full
25-example run, and leave a forensic artifact when an example crashes. TiCat3/4
currently crash during step-3 re-encode and leave no `step2.oin` and no error
file, so there is nothing to inspect after a run.

## The fix — one file: `tests/integration/verify_roundtrip.py`

Line numbers are as of 2026-07-04; re-confirm the anchor before editing.

### 1. New `--only` argument

In `main()`, immediately after the existing `--limit` argument
(~lines 159-161):

```python
parser.add_argument(
    "--only",
    type=str,
    help="Run only examples whose name contains this substring (case-insensitive).",
)
```

### 2. Filter examples

Immediately after `examples = get_examples()` (~line 176) and **before** the
`if args.limit:` block, so `--only` narrows first and `--limit` can further cap:

```python
if args.only:
    needle = args.only.lower()
    examples = [e for e in examples if needle in e.name.lower()]
    print(f"Filtering to {len(examples)} example(s) matching '{args.only}'.")
```

### 3. Persist a failure artifact

In the outer `except Exception as e:` block at the end of the per-example loop
(~line 433), keep the existing `print` / `traceback.print_exc()` /
`reporter.log_failure(...)`, and add a guarded artifact dump. `i` (loop counter)
and `safe_name` (~line 188) are always defined at this point; `base_name` is
defined inside the `try`, so build the filename from `safe_name`:

```python
except Exception as e:
    print(f"Unified Test FAILED: {e}")
    import traceback

    traceback.print_exc()
    if output_dir:
        try:
            err_path = os.path.join(output_dir, f"Ex{i}_{safe_name}_error.txt")
            with open(err_path, "w") as fh:
                fh.write(f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}")
        except Exception:
            pass
    reporter.log_failure(test_name, f"Exception: {str(e)}")
```

## Acceptance (exact commands + expected)

1. Single-example filter:
   ```
   uv run python tests/integration/verify_roundtrip.py --only Ferrocene 2>&1 | grep -c "Running Round-Trip for Example"
   ```
   Expect: `1` (exactly one example run). `--only TiCat` → `4`.

2. Full run unchanged. Capture before/after and diff the pass/fail set:
   ```
   uv run python tests/integration/verify_roundtrip.py --output-dir /tmp/rt40_after
   python3 -c "import json; d=json.load(open('/tmp/rt40_after/summary_roundtrip.json')); print(d['passed'], d['failed'])"
   ```
   Expect: `19 6` — same example set, same split as the 2026-07-04 baseline in
   `/tmp/baseline_rt/summary_roundtrip.json` (Ferrocene, TiCp2Me2, TiCat1–4 fail).

3. Crash leaves an artifact:
   ```
   uv run python tests/integration/verify_roundtrip.py --only TiCat3 --output-dir /tmp/rt40_err
   ls /tmp/rt40_err/Ex*_TiCat3_error.txt
   ```
   Expect: the error file exists and names the crash
   (`cannot unpack non-iterable NoneType` pre-TASK-41; `get_lig_mol failed`
   post-TASK-41).

4. Suites unaffected (harness-only change):
   ```
   uv run python -m unittest discover tests 2>&1 | tail -3          # 55 OK
   uv run python -m unittest discover tests/unit 2>&1 | tail -3     # OK, skipped=3
   uv run python tests/integration/verify_xyz_to_oin.py 2>&1 | tail -2  # 25 Passed
   ```

## Constraints / DO NOT

- Touch only `tests/integration/verify_roundtrip.py`.
- Do NOT change the example definitions in `verify_xyz_to_oin.py`
  (`get_examples()`), `normalize_oin_for_comparison`, or the RMSD threshold.
- Run `uv run ruff format tests/integration/verify_roundtrip.py` before finishing.

## On completion

Set `Status: DONE`, append a dated Log entry to `spec/worklog/NOTES.md` (the
`--only` arg, the filter placement, the failure-artifact dump, acceptance
results), and note in the handoff status table that WS-0 landed. Do NOT commit
unless asked — leave staged for review (`git add` scoped to this one file).
