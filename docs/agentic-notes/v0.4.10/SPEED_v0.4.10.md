# v0.4.10 · Cost per attempt — what was measured, what was refuted, what was not done

> **The release changed no answer. That was the point, and it is the part with the strongest
> evidence: ARM 1 62/62 and ARM 2 90/90 byte-identical, for both lanes, on all four gate runs.**
>
> Everything else in this document is smaller than it looks, because **every speed number here is
> bimodal by molecule** — the same change measures −50.2% on one molecule and nothing on the next.

Companion lane documents: `GATE_REPAIR_v0.4.10.md` (L0) · `DEAD_SCAN_v0.4.10.md` (LA) ·
`CIP_MEMO_v0.4.10.md` (LB) · `NEGATIVE_LANES_v0.4.10.md` (LC, LD).

---

## 1. The first finding was that the arbiter was broken

Before any lane ran, the charter's own instruction — run the gate on the pristine tree and see
whether anything changed — produced this:

```
error: tests/fixtures has 62 .xyz files, expected exactly 61
[gate/arm1] FAIL: gate_arm1_encode.py exited 1
```

ARM 1 had been exiting **before comparing anything** since `dd51a515`. **v0.4.9 froze a 328-molecule
benchmark and measured a 0.28% noise floor while the encoder arm of its own gate was non-runnable.**
Behind that refusal sat a second, real drift — the v0.4.7 `xyz2mol` → `perception_tmc` rename, which
was behaviour-neutral but not *string*-neutral, and ARM 1 hashes error strings on purpose. The other
60 rows were byte-identical, so the encoder had not moved. Full account in `GATE_REPAIR_v0.4.10.md`.

> **A gate that fails before it compares is indistinguishable from a gate that is merely
> inconvenient to run, and it silently stops covering everything else it was watching.**

## 2. The charter was re-ranked before it was executed

The two chartered lanes were written before v0.4.9 profiled anything. v0.4.9's close-out named two
targets that are larger, cheaper and better-evidenced than either, and re-sized the original Lane 1
down to third. **Executing the charter as written would have spent the release on its third-best
target.**

| lane | target | outcome |
|---|---|---|
| **L0** | the gate itself | repaired; 2 rows re-frozen, 60 verified identical |
| **LA** | the discarded `.index()` scan | **landed, on by default, no lever** |
| **LB** | `_reparse_cip_label_once` memo | **landed, default OFF, not promoted** |
| **LC** | SVD in `_finalize_positions` (charter Lane 1) | **NEGATIVE — the premise does not reproduce** |
| **LD** | per-attempt redundancy (charter Lane 2) | **already fixed by an earlier release** |

## 3. What got faster, and where it did not

Two changes, four molecules, **and the pattern matters more than any single number**:

| molecule | class | Lane A (`.index()` deleted) | Lane B (CIP memo) |
|---|---|---|---|
| `CAHQEJ_comp_0` | eta, `[Ni_TPL]`, 2 haptic | **−32.9%** | −2.4% |
| `VAFMIA_comp_0` | `[Cu_LIN]`, adamantyl NHC | — | **−86.7%** |
| `FOSNEI_comp_0` | non-eta, boron cage | **+0.3% (nil)** | — |

> ⚠ **Lane A's figures here correct an earlier pair.** A first A/B, taken while four gate processes
> were competing for the box (load **35**), reported **−50.2%** and **+9.6%**. Quiet-box re-measurement
> gives **−32.9%** and **+0.3%** — the gain was over-stated by 17 points and the null was buried in
> 30% within-arm spread. The lane's merge commit quotes the contended figures; these are the ones to
> cite. See §7.

**The two lanes are complementary, and neither touches the third class.** That is the same shape
v0.4.9 reported as "three molecules, three cost regimes", and it is the reason this release states
per-molecule numbers rather than a headline percentage. A corpus figure would require a 5000-molecule
sweep at ~55 CPU-h, which **this release did not run** (see §6).

### The bimodality is attributed, not just observed

Wrapping `Molecule.__eq__` — what `list.index` actually calls — on the pristine tree:

| molecule | `__eq__` calls | cost | % of generate |
|---|---:|---:|---:|
| `CAHQEJ_comp_0` | **99** | **38.52 s** | **38.8%** |
| `FOSNEI_comp_0` | **3** | 0.03 s | **0.0%** |

99 comparisons versus 3. **The `FOSNEI` null is explained rather than excused as noise.** And the
two independent measurements agree on the size of the thing: the attribution puts the scan at
**38.8%** of that run's generation, the quiet A/B at **32.9%** — both *roughly a third of this
molecule*, with the difference inside the 12.7% run-to-run spread arm A shows here.

It also settles the handoff's under-count. v0.4.9 quoted **22%**, which is the `numpy.linalg.eig`
line alone; the whole of `get_c_eig_list` is **38.51 s over 198 calls**, and those 198
eigendecompositions match v0.4.9's count exactly — same work, different attribution boundary. The
difference is the Coulomb matrix built before each decomposition, on both operands.

`VAFMIA` at 10.87 s now clears 30 s. It was the worst of the **eleven** molecules that exceeded
budget in *both* of v0.4.9's arms. **The other ten are untouched, and `max(elapsed_s) < 30 s` is not
delivered.**

## 4. Byte-identity — the primary gate, not a secondary check

| run | arm | verdict |
|---|---|---|
| Lane A | ARM 1 (encoder, 62 fixtures) | **PASS — byte-identical to golden** |
| Lane A | ARM 2 (round trip, 90 molecules, fast band) | **PASS — 90 gated** |
| Lane B | ARM 1 | **PASS — byte-identical to golden** |
| Lane B | ARM 2 | **PASS — 90 gated** |

Plus identical generated-structure fingerprints across every A/B run, reproduced independently on
three separate checkouts.

**Lane B needed both arms and that is not obvious.** `_reparse_cip_label_once` is reached from
`metallogen_adapter._template_sp3_label` (generator) *and* from `chirality.recover()` →
`_reparse_aromatic_cip_label` (**encoder**). A memo transparent for one and not the other would pass
a single-arm gate. The charter frames this release as generator-side work; that framing is what would
have hidden it.

## 5. Predicted vs actual

| | predicted | actual |
|---|---|---|
| `byte_exact` | **FLAT by construction** — "if it moves, the change was not byte-identical; treat that as a stop condition" | **FLAT.** 4/4 gate runs byte-identical. Not re-measured at corpus scale — see §6 |
| median down, p90 down | yes | **not measured at corpus scale.** Per-molecule: −50.2%, −86.7%, −2.4%, nil |
| `> 30 s` count down "by an amount the profiling will predict" | yes | **not re-measured.** One named molecule crosses 30 s |
| `facmer_divergent` does not rise | — | not re-measured (requires the sweep) |
| suite ≥ 877 OK | — | **946 OK** on merged `main` — exactly 930 baseline + 16 new lane tests |

**Two of the four predictions could not be evaluated**, for one reason: they are corpus quantities and
this release did not run a corpus sweep. That is stated rather than papered over — the charter's
"median down; the `> 30 s` count down by an amount the profiling will predict" assumed a sweep that
costs 55 CPU-h.

## 6. What this release did NOT do

- **No 5000-molecule sweep.** `byte_exact` FLAT rests on *construction plus four gate runs*, not on a
  fresh corpus measurement. For a byte-identical release that is a defensible substitute — the gates
  are the proof that no answer changed — but it is not the same claim, and the roadmap's gap table is
  therefore **unchanged** and was not re-derived.
- **`OIN_MEMO_CIP_REPARSE` was not promoted.** The charter permits same-release promotion *"only if
  byte-identity holds on the whole benchmark, including the fast control"*. This release ran the
  **fast band — 90 of 328 molecules**. That is the control, not the whole benchmark.
  **Promotion gate:** the full 328-molecule cohort with the lever on, ~10 CPU-h sharded 6-way
  (`--shard i:6`, **1-based**), against the existing frozen golden. No new baseline run is needed —
  the golden *is* the baseline.
- **No time limit added.** v0.4.9 owns stopping work, shipped no CBC `timeLimit`, and refuted CBC as
  a target (2.1% on `FOSNEI`, 0.35 s on `VAFMIA`). So there was nothing to re-litigate and nothing
  added here.
- **Nothing pushed.** Standing instruction on this repo.

## 7. Method notes worth keeping

- **Byte-identity gates are load-immune; wall-clock is not.** Running four gates beside two A/Bs drove
  a 12-core box to **load 35** and made a −50% and a null statistically indistinguishable. Gates can
  be parallelised freely. Timing runs must be serial, and the box's contention state must be stated
  with every number.
- **`pkill -f` kills the parent, not its multithreaded BLAS children.** Two orphans at 339% CPU each
  survived a kill and silently corrupted every measurement in flight. Verify by PID after any kill.
- **A "pristine" arm pinned to a *path* stops being pristine the moment you merge into it.** The first
  attribution run pointed its baseline at `main` *after* Lane A had been merged there, and duly
  reported zero difference — a null that reads exactly like "no effect". The pristine arm must be a
  worktree pinned to a **commit**.
- **`pgrep -f` in a wait loop matches the waiter's own command line.** A queued job that waited on
  `pgrep -f "run_rest.sh"` never ran, because that string appears in its own `bash -c` body. Wait on a
  marker **file**.
