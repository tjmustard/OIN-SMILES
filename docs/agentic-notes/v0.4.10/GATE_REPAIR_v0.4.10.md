# v0.4.10 · ARM 1 was dead, and that is how a second drift hid behind it

> **The encoder half of the byte-identity gate refused to run for an entire release cycle.**
> v0.4.9 froze a 328-molecule runtime benchmark, measured a noise floor to 0.28%, and shipped a
> bound — with ARM 1 hard-failing on every invocation. Nobody saw it, because a gate that exits
> non-zero *before comparing anything* looks like a tooling error, not like a missing result.

Found before any v0.4.10 lane ran, by doing the one thing this release's charter asks for first:
run the gate on the pristine tree and see whether anything changed.

## 1. What happened

```
$ bash tools/gate_v047.sh arm1
error: tests/fixtures has 62 .xyz files, expected exactly 61 -- the fixture set has
drifted since the golden manifest was frozen; update --expect-n deliberately if this
is intentional
[gate/arm1] FAIL: gate_arm1_encode.py exited 1
```

`tests/fixtures/ULODUU_comp_0.xyz` landed in **`dd51a515`** (*"CORRECTION(boron): geometry is a
refuted discriminator -- TET joins the safe set"*). `tools/gate_v047_arm1_golden.tsv` was never
extended, so `EXPECTED_FIXTURE_COUNT = 61` refused every run from that commit onward.

**The guard is correct and it worked.** A fixture added or removed *is* gate-relevant drift, and the
alternative — silently gating a different corpus than the manifest was built against — is worse.
The defect is not the refusal. The defect is that the refusal was left standing.

## 2. The cost of leaving it standing

Regenerating at `--expect-n 62` on the pristine post-v0.4.9 tree (`3077282a`) and diffing the
**pre-existing 61 rows** against the frozen golden — *before* accepting anything — gives exactly two
differences:

| row | golden | fresh | cause |
|---|---|---|---|
| `ASISAX_comp_0` | `ERROR:ValueError:`**`xyz2mol`**` failed: get_lig_mol failed for ligand fragment #0 (SMILES: '[H]C1C…')` | `ERROR:ValueError:`**`perception_tmc`**` failed: …` | the **v0.4.7 rename**, `xyz2mol.py` → `perception_tmc.py` |
| `ULODUU_comp_0` | *(absent)* | `3651a131…9701`  192  eta | the v0.4.9 fixture |

**The other 60 rows are byte-identical.** The encoder has not moved.

The `ASISAX` row is the point. The v0.4.7 rename was recorded as **behaviour-neutral**, and it was —
of *behaviour*. It was not neutral for the error **string**, and ARM 1 hashes error strings on
purpose: its own docstring says *"Two revisions must raise the SAME error, not merely agree when
both succeed."* So a real, if benign, gate-relevant change sat undetected for two releases behind a
guard that had already stopped reporting.

> **The transferable finding:** a gate that fails *before* it compares is indistinguishable from a
> gate that is merely inconvenient to run, and it silently stops covering everything else it was
> watching. Both of this session's ARM 1 findings are cheap; the expensive part was that neither was
> visible while the count guard was tripped.

## 3. What shipped

| change | file |
|---|---|
| golden re-frozen at **62 rows** from the pristine `3077282a` tree, both diffs deliberate | `tools/gate_v047_arm1_golden.tsv` |
| `EXPECTED_FIXTURE_COUNT` 61 → 62, with the provenance of both rows in the comment | `tools/gate_arm1_encode.py` |
| `#DONE` sentinel check 61 → 62 | `tools/gate_v047.sh` |
| docstring: *"WHY 61, NOT A NAME LIST"* → *"WHY A COUNT, NOT A NAME LIST"*, plus the standing instruction to **diff the pre-existing rows before re-freezing** | `tools/gate_arm1_encode.py` |

Deliberately **not** changed:

- **The counts stay hardcoded and independently asserted.** Deriving `EXPECTED_FIXTURE_COUNT` from
  the golden's row count would make this class of drift self-heal — which is exactly wrong. A
  fixture added without a matching golden row is the thing the guard exists to catch; making the two
  agree automatically deletes the guard.
- **The error message is not reverted.** `perception_tmc failed:` names the module that actually
  exists. The golden was stale, not the code.

## 4. Reproducing

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src   # rdkit pinned ==2025.9.3

# the refusal, on the pre-repair tree
bash tools/gate_v047.sh arm1                      # exits 1 before comparing anything

# the regeneration + the diff that licenses the re-freeze  (~15 min: several fixtures
# take the forked-resonance path at a 120 CPU-s budget each)
$V tools/gate_arm1_encode.py --fixtures-dir tests/fixtures --expect-n 62 > /tmp/fresh.tsv
diff <(grep -v '^#' tools/gate_v047_arm1_golden.tsv) <(grep -v '^#' /tmp/fresh.tsv)
# expect: exactly the two rows in §2, and nothing else

# post-repair
bash tools/gate_v047.sh arm1                      # PASS, 62 gated
```

⚠ ARM 1 is **not** the "fast" arm its own header claims. Several fixtures route through
`perception_tmc._resonance_candidates_isolated`, a forked child bounded by
`_RESONANCE_CPU_BUDGET_S = 120` **CPU seconds**, so a full run is ~15 minutes rather than the
implied minute or two. That cost is the R3 regime deferred to v0.4.11.
