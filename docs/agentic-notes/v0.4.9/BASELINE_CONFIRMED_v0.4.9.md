# v0.4.9 · Phase 0 — the baseline v0.4.9 is measured against

**Verdict: `byte_exact = 72.46%` stands. No re-baselining. Every v0.4.9 lane may quote it.**

v0.4.8 re-scored the 5000-molecule sweep **offline**, from stored structures, in 334 s instead
of the 55 CPU-h a re-sweep costs. Before building a release on that number, v0.4.9 waited for
the live confirmation run to finish and checked it independently.

The live arm's own verification is already recorded by the lane that ran it — see
`docs/agentic-notes/v0.4.8/HONEST_BASELINE_v0.4.8.md` (commit `ce2265ba`): 275/275
byte-identical `smiles_2_indep` on the same stored conformer, flat bucket distribution on fresh
generation, and a live scored→honest delta of **−11.0 points** against the corpus **−10.34**.
This note records only what that check did not cover.

## 1. Per-molecule agreement, offline vs freshly generated

The v0.4.8 note compares the two arms **bucket by bucket**. Per molecule:

| offline honest → live honest | n |
|---|---:|
| `byte` → `byte` | 214 |
| `other` → `other` | 61 |
| `FAIL` → `FAIL` | 21 |
| `FAIL` → `other` | 2 |
| `byte` → `FAIL` | 1 |
| `other` → `FAIL` | 1 |
| **agreement** | **296 / 300 = 98.7%** |

A flat bucket distribution can hide equal and opposite per-molecule churn; this shows there is
none. **All four disagreements are generation-side and three of the four move toward failure** —
none moves a molecule from failing to passing, which is the only direction that would have
undermined the correction.

**Encoder drift: 0/300.** Every `smiles_1` is byte-identical between the v0.4.6 sweep and live
v0.4.8 code — the corpus-scale encoder-identity gate (4985/4985) reproduced on a fresh draw.

## 2. 🔴 A string-equality count is not a pass count

Counting honest passes as "`honest_class` ends in `->byte`" gives **3631**. The frozen bucket
report says **3623**. The report is right, and the eight-molecule gap is not an error:

```
ALEMOT_comp_0  DOCPAO_comp_0  MIBFEL_comp_0  NEFNER_comp_0
NOYTUS_comp_0  ULOQIX_comp_0  UPABUK_comp_0  XAKCAP_comp_0
```

Every one has `smiles_1 == smiles_2_indep` **byte-for-byte** and `status: "failed"`. The bucket
report applies the harness `status` gate *before* the string comparison. The naive count does not.

**Eight is exactly the population v0.4.8 Lane 2 measured**, and `XAKCAP_comp_0` is one of its
pinned fixtures: these are the molecules the **atom-count gate** catches — they encode to a
string identical to their input while having lost atoms, so **no string comparison can catch
them**. Phase 0 arrived at the same eight from the opposite direction, which is the strongest
available confirmation that the gate is load-bearing.

> **The trap, stated generally:** any instrument that reports `byte_exact` by comparing strings
> alone over-counts by exactly this population. `tools/roundtrip_bucket_report.py` gets it right;
> an ad-hoc `honest_class.endswith("->byte")` does not. This bit an analysis in this very session.

## 3. Runtime, carried into v0.4.9 unchanged

| | median | max | `> 30 s` |
|---|---:|---:|---:|
| live300 (n = 300) | 7.53 s | 448.8 s | 61 (**20.3%**) |
| 5k sweep (n = 5000) | 7.19 s | 759.9 s | 994 (**19.88%**) |

The `> 30 s` rate reproduces to within half a point on an independent draw, so v0.4.9's runtime
baseline is not a single-run artifact. The 448.8 s maximum is the same multi-tier accumulation
documented in `ELAPSED_S_IS_A_SUM_v0.4.9.md` — **not** one attempt overrunning.

## 4. Reproducing §1 and §2

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src   # rdkit pinned ==2025.9.3

$V - <<'PY'
import json, glob, os
from collections import Counter
L = "tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-live300/individual_reports"
H = "tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/individual_reports"
def v(s1, s2):
    if not s1 or not s2: return "FAIL"
    return "byte" if s1.strip() == s2.strip() else "other"
agree, drift = Counter(), 0
for f in glob.glob(os.path.join(L, "*.json")):
    n = os.path.basename(f)[:-5]
    lr, hr = json.load(open(f)), json.load(open(os.path.join(H, n + ".json")))
    agree[(v(hr.get("smiles_1"), hr.get("smiles_2_indep")),
           v(lr.get("smiles_1"), lr.get("smiles_2_indep")))] += 1
    a, b = (lr.get("smiles_1") or "").strip(), (hr.get("smiles_1") or "").strip()
    drift += bool(a and b and a != b)
print(agree.most_common()); print("encoder drift:", drift)
PY

# §2 -- the eight the status gate removes and a string comparison cannot
$V - <<'PY'
import json, glob, os
D = "tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/individual_reports"
odd = [os.path.basename(f)[:-5] for f in glob.glob(os.path.join(D, "*.json"))
       if (lambda r: (r.get("honest_class") or "").endswith("->byte")
                     and r.get("status") != "success")(json.load(open(f)))]
print(len(odd), sorted(odd))
PY
```
