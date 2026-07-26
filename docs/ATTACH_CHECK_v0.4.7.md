# `OIN_ATTACH_CHECK` — the attachment check for `OIN_ACCEPT_SCORED` (v0.4.7, lane L5-attach)

**Status: falsification PASSED · check BUILT, default-OFF · four gates RE-RUN.**
**Verdict: 2 of the 4 bar criteria met. DO NOT PROMOTE to default-ON (§5.1) — but never run
`OIN_ACCEPT_SCORED` without it (§5.2). The product call is left explicit in §5.3.**

| | bare lever | **lever + check** |
|---|---|---|
| **G3** bytes | 98/98 identical | **97/97 identical** ✅ |
| **G1** structure | clash 77→**105**, severe 5→**14** ❌ | **77→82, severe 5→5** ✅ |
| **G2** independent re-perception | 88/98→62/98, **26 regressions** ❌ | **87/97→78/97, 9 regressions** ⚠ |
| **G4** pass rate | 0 regressions ✅ | 1 regression (timeout-shaped) ⚠ |
| speedup vs default | **3.33×** | **1.43×** (43% retained) ⚠ |

This lane exists because `docs/ACCEPT_SCORED_v0.4.7.md` §5.4 named exactly one route from
`DO NOT PROMOTE` to promotion: *build the attachment check of §6 as the lever's missing safety
condition, then re-run the four gates.* §6.5 left one test outstanding and made it a gate —
**if the predicate cannot separate arm A's accepted conformers from arm B's, §6 is wrong and
the proposal dies.** That test ran first and it is §1 below.

> **Headline, stated before the detail:** the falsification passed, but **not for the
> predicate §6 proposed.** Both candidates in §6.4 — count-based and set-based — wrongly
> reject conformers that round-trip, on 11 of 22 and 3 of 22 respectively. §6.4's stated
> ceiling of "6/8 count-based, 7/8 set-based" is not reachable by those predicates at any
> false-positive rate worth paying. A different reference quantity reaches 7/8 with **zero**
> false positives. The conclusion of §6 survives; its proposed mechanism does not.

---

## 1. THE FALSIFICATION (§6.5) — PASSED, and it corrected the proposal

### 1.1 What was run

`tools/ab_accept_scored.py` gained `--dump-xyz`, which persists each **accepted** conformer's
coordinates plus the generator's claimed metal-donor set. Generation is then paid once and any
number of candidate predicates can be scored offline (`tools/attach_probe.py`) — which is what
made it affordable to test four predicates instead of assuming one.

Cohort: the 22-molecule gap cohort minus `XIQKOY_comp_0` (the DEAD class — produces no
structure in either arm, and 720 s for nothing). **21 molecules, both arms, 40 accepted
conformers.** `--timeout 300 --hard-cap 500 --workers 2`.

**The run reproduced the promote lane exactly**, which is what makes the comparison legitimate:

- the same 8 `indep` regressions, in the same order: `KAQDOV`, `ZITSIE`, `FEXYOZ`, `DAKGON`,
  `RATPEK`, `MEDZUR`, `POVPIA`, `HIDCIH_comp_1`;
- **20 comparable / 20 byte-identical `sha256(smiles_2)` / 0 different**;
- `POVPIA_comp_0` arm A again at `16 clashes / worst 0.4344`, matching both prior runs. That is
  now the **third** independent reproduction of the generator's determinism.

### 1.2 The trap (§6.1), and how the probe stays out of it

`_coordination_vectors` reads donors from `metal.GetBonds()`. **A ligand that has left the
coordination sphere keeps its bond object**, so a check built there certifies exactly the
structures it exists to catch.

In this probe the generator's graph supplies only the **reference** — which atoms to look at.
The **measurement** is `xyz2AC_obabel(Z, coords, tolerance=0.5)` and the metal's row of the
resulting adjacency matrix: coordinates and element numbers, nothing else. A detached ligand
stays in the reference and vanishes from the measurement, and that mismatch is the signal.

An **alignment guard** was added and it is not ceremony: `res.xyz` is written from the raw
MetalloGen mol while `res.mol` is the contract mol, so their shared atom ordering is a
convention (the harness's own `get_oin_string(gen.mol, coords)` relies on it) rather than a
guarantee. If it ever broke, the reference would index into the wrong atoms and the probe
would manufacture a separation out of a numbering bug. **OK on all 40.**

### 1.3 Result — the two predicates §6.4 proposed BOTH fail

Scored against the `indep` oracle. "Separates" = accepts arm A's conformer **and** rejects arm
B's, on the 8 known regressions. "False positives" = conformers that pass independent
re-perception and would nevertheless be rejected — the cost of shipping the check.

| predicate | separates /8 | **false positives** /22 indep-passing |
|---|---|---|
| **P1** `|actual| == |claimed|` — §6.4's *count-based* | 5 | **11** |
| **P2** `claimed ⊆ actual` — §6.4's *set-based* | 5 | **3** |
| **P3** sites(actual) == OIN slot count, raw donor row | 7 | **9** |
| **P4** site coverage, raw donor row | 6 | 0 |
| **P3f** sites == slots, *filtered* donor row | 6 | 0 |
| **P4f** **site coverage, *filtered* donor row** | **7** | **0** |

**Why P1 and P2 fail, and it is one reason:** the number and identity of metal-bonded *atoms*
in the generator's graph is **not conserved** against the encoder's perception, even on
structures that round-trip. `MEDZUR_comp_0`'s contract mol claims 10 metal bonds while the
coordinate path perceives 7 — its Cp ring is already slipped η5→η3 — and it round-trips
`[Ru_TET]` on both sides, because the round-trip key folds that. Anything referenced to the raw
donor-atom count therefore rejects good conformers. §6.4 reasoned about the right *phenomenon*
and picked the wrong *invariant*.

What **is** conserved is the set of coordination **sites**. So the shipped predicate is:

> **every coordination site the generator claims must still retain at least one atom inside
> bonding distance of the metal**

with sites formed by the same 1.6 Å transitive grouping the encoder's hapticity reduction uses,
so a whole Cp ring is one site. Ring slip leaves the site populated; a ring that has drifted out
of the sphere empties it.

### 1.4 The encoder's ring-carbon filter is the mechanism, not a detail

`_get_tmc_mol_impl` drops a coordinating ring carbon when a neighbouring N/O/P/S is also
coordinating — the heteroatom is the real donor. Compare P4 with P4f above: omitting that filter
costs a regression (6/8 → 7/8 when it is applied). The molecule that moves is `DAKGON_comp_0`,
and it is worth stating in full because it explains *why* the filter matters:

- arm B's **filtered** donor set has **6 atoms — the same size as the claim** — but not the same
  atoms;
- the filter drops the NHC carbene carbons **precisely because the adjacent N has come within
  bonding distance**;
- that *is* the C→N donor reassignment §4.5 recorded. Two claimed sites lose every member and
  the conformer is rejected.

So the filter is not a nuisance correction applied to reduce false positives — it is the
donor-reassignment rule, and it is what makes DAKGON detectable at all. §6.4 predicted DAKGON
would need "set not count"; it was right that DAKGON is reachable and wrong about the mechanism.

### 1.5 Two things a narrower test would have gotten wrong

Recorded because this lane's instruction was to hold itself to the sibling lane's standard, and
both were caught only by widening a measurement that had already "passed".

1. **At n=3 the unfiltered site-count predicate looked perfect** — 3/3 separation, 0 false
   positives — and I had already written it up as the winner. At n=21 it wrongly rejects **9 of
   22** round-tripping conformers. The entire difference is the ring-carbon filter, which only
   bites on chelates that the 3-molecule sample did not contain. This is the same correction
   §4.8 had to apply to G1, arrived at independently and from the same cause: a sample too small
   to contain the hard case.
2. **The cheap shortcut is wrong.** Taking the metal's row by the distance criterion alone —
   skipping `xyz2AC_obabel`'s valence-cap pruning loop — runs in **0.17 ms** instead of 24 ms, a
   further 100×. It **disagrees with the real row on 11 of 40 conformers.** The pruning loop can
   drop a metal–donor bond (the code's own `DUDREA_comp_0` bridging-hydride example), so the
   full call is required for the check to measure what the encoder measures. This was measured
   rather than assumed specifically because it was attractive.

A **prediction was recorded before the 21-molecule run completed** (`spec/handoffs/v0.4.7/
PROGRESS-attach.md`): "P3 will separate 7 of 8; POVPIA will not separate." The separation count
was right; the false-positive rate was not anticipated at the magnitude it appeared, and the
predicate had to change as a result.

### 1.6 The residual, stated plainly

`POVPIA_comp_0` is **not** caught. Its metal donor set is intact and every claimed site is
populated; the defect is ligand-internal — a hydrogen detaches and the amine reads as an imine.
No metal-centred predicate can see that, which is exactly what §6.4 said.

**7 of 8 is the ceiling and it is reached. It is not 8 of 8, and this document does not claim it.**

### 1.7 Does it reject real structures? No.

Run on the 21 **crystal inputs** themselves (`tools/attach_probe.py --inputs`): every molecule
detects a full donor set, with metal–donor distances spanning **1.81–2.84 Å**. §6.3 quoted
1.8–2.6 Å on eight molecules; the wider set widens the range slightly (GAVSED to 2.84,
DEJHEF to 2.77) and does not change the verdict.

The stronger version of the same control is the false-positive column of §1.3: **0 of 22**
conformers that pass independent re-perception are rejected by P4f.

### 1.8 Cost

**7–81 ms per conformer** (median **23.6 ms** over 40 evaluations, 47–109 atoms) against the
**48–57 s** strict `_reencode_oin` it stands in for — consistent with §6.3's measurement, and
roughly 1000×. First-call values near 900 ms are import/warm-up and are excluded.

---

## 2. THE IMPLEMENTATION

`src/oinsmiles/generation/attach_check.py`, wired into
`metallogen_adapter._reencode_key_matches` and registered in `oin/levers.py::_HELD_OFF` with its
measured justification, per project convention.

```
if not independent_confirm and fast is not None:          # the OIN_ACCEPT_SCORED branch
    if require_no_stretch and clash.mol_stretched_bond_count(m) > 0:
        return False
    if lever_enabled("OIN_ATTACH_CHECK") and not conformer_ligands_attached(m):
        return False                                      # <-- the check
    return True
```

**Default OFF.** It sits *inside* the lever's own branch, so with `OIN_ACCEPT_SCORED` off it is
unreachable — the default path is byte-identical by construction, and this is verified two ways:
a unit test asserts on the source that the call appears only after the guard, and generating
`ODEWID_comp_0` and `YIYGAP_comp_0` with `OIN_ATTACH_CHECK=1` and `OIN_ACCEPT_SCORED` off
reproduces `f00fdf52371f5cdb` and `816ebda751edf06d` — the promote lane's own recorded shas.

**Structural defence against the §6.1 trap:** `ligands_attached()` takes plain atom indices,
element numbers and coordinates. It cannot be handed connectivity. The only place the module
touches a bond is the convenience wrapper that reads the *reference* donor list off the metal,
and that is one function with the reason written above it.

**Verification against the offline probe:** the production function reproduces the probe's
verdict on **40/40** conformers, so §1's table describes the shipped code and not a prototype.

**17 unit tests** (`tests/unit/test_attach_check.py`). The load-bearing one is
`TestTrapAvoidance::test_detached_ligand_is_caught_even_though_its_bond_survives`: an RDKit mol
whose metal–N bond object is intact while the N sits 6 Å away. A `GetBonds()`-based check passes
that structure; this one must reject it. There is also a ring-slip test (η5→η2 must be
**accepted**) and its mirror (the same rigid ring moved out of the sphere must be **rejected**).

> One fixture was wrong on the first attempt and the fix is documented in place: moving two ring
> carbons out of the ring plane shatters the *ligand* rather than slipping it, so those two
> carbons form a site of their own with nothing bonded in it. The check rejected it and was
> right to. The corrected fixture translates the ring **rigidly**.

### 2.1 Scope limit — an acceptance condition, not a return condition

This guards **acceptance**. On the CHEAP_ONLY class (`docs/eta_accept_gap_cohort.md`, and §3 of
the promote lane's report) *nothing is ever accepted*: the pool fills to completion and
`_select_by_geometry` returns a best-by-geometry conformer that never faced any acceptance test
at all. Those molecules are unaffected by this check, with the lever on or off.

**This is not hypothetical and it is visible in §3's results.** It is the single most important
qualification on everything below.

---
## 2.2 ⚠ The first implementation was a SILENT NO-OP, and a full A/B reported it as a result

Recorded in the body rather than a footnote, because the way it was caught is the transferable
part.

`_reencode_key_matches` was wired to pass **`m`** to the check. `m` is MetalloGen's `Molecule`,
not an `rdkit.Chem.Mol` — it has no `GetAtoms`. Every call raised `AttributeError`, the
abstain-on-error branch swallowed it, and the check never ran.

**A complete 21-molecule arm-C run then reported exactly what a genuine null result looks
like:** identical `sha256(smiles_2)` on all 20 comparable molecules, identical `clash_vdw`,
identical `worst_overlap`, the same 8 `indep` regressions, `worse 0 · better 0 · identical 20`
against arm B. Written up as-is, that is the sentence "the check recovers the structure cost on
zero molecules" — a clean, publishable negative.

What prevented it: **20 of 20 bit-identical is not a plausible null result, it is a plausible
not-wired-up.** A real safety condition that rejects conformers changes *something* somewhere.
A telemetry probe then settled it in one line — `adapter.attach_check_rejected` fired **0**
times while the returned conformer **failed** the check when tested directly.

This is the same defect shape the promote lane had to fix in its own instrument
(`clash.mol_clash_count` returning 0 on `AttributeError`, so `0 over 0` and `0 across 22`
printed identically). I rebuilt it inside the safety check itself, one lane later.

Fixed three ways: pass `cmol`; **record `adapter.attach_check_unevaluable` instead of abstaining
silently**, so "never rejects" and "never runs" can no longer look the same; and two regression
tests, one asserting the source calls `conformer_ligands_attached(cmol)` and not `(m)`.

*(Telemetry is opt-in — `_telemetry.record()` no-ops unless `OIN_TELEMETRY=1`, deliberately, so
the instrument cannot perturb what it measures. Observing an abstention in production therefore
requires setting it. That is a real limitation of this loudness fix and it is stated rather than
glossed.)*

The invalid run is kept as `spec/handoffs/v0.4.7/runs/INVALID_c21_noop.json`. **No number in §3
comes from it.**

---
## 3. THE FOUR GATES, RE-RUN WITH THE CHECK ENABLED

Method: **arm C** = `OIN_ACCEPT_SCORED=1` **+** `OIN_ATTACH_CHECK=1`, run as a `--single-arm`
pass and compared against stored arms A and B. Re-measuring arm A costs more than the rest of
the lane (median 21 s, tail to 400 s); §4.2/§4.4 of the promote lane established generator
determinism at coordinate level, and this lane reproduced it a **third** time, so a stored arm
A is a valid comparator. The premise is checked, not assumed: `sha_in` must agree on every
molecule, and it does.

### 3.1 Gap cohort (n=21, both arms fresh this session) — **G1 and G2 both stop failing**

| gate | question | A → B (bare lever) | **A → C (lever + check)** | verdict |
|---|---|---|---|---|
| **G3** | does `smiles_2` stay byte-identical? | 20/20 identical | **20/20 identical, 0 divergent** | **PASS** |
| **G1** | does structure quality degrade? | clash 16→2, severe 7→0 | **clash 16→0, severe 7→0, worst_overlap min 0.4344→0.75** | **PASS** |
| **G2** | what does dropping independent re-perception cost? | 15/20 → 7/20, **8 regressions** | **15/20 → 13/20, 2 regressions** | **mostly recovered** |
| **G4** | does any passing molecule stop passing? | 19/21, none | **19/21, none** | **PASS** |

Against the bare lever, arm C is a strict improvement on the arm that mattered:
**`indep` 7/20 → 13/20, six fixes, zero regressions.**

#### Recovery, per molecule, on the 8 the bare lever broke

| molecule | A `indep` | B | **C** | A s | B s | **C s** | |
|---|---|---|---|---|---|---|---|
| DAKGON | ✓ | ✗ | **✓** | 17.8 | 3.6 | 20.2 | **RECOVERED** |
| FEXYOZ | ✓ | ✗ | **✓** | 12.7 | 2.9 | 21.7 | **RECOVERED** |
| HIDCIH_comp_1 | ✓ | ✗ | **✓** | 64.2 | 2.0 | 69.0 | **RECOVERED** |
| KAQDOV | ✓ | ✗ | **✓** | 202.6 | 7.1 | 248.9 | **RECOVERED** |
| RATPEK | ✓ | ✗ | **✓** | 53.7 | 9.2 | 69.2 | **RECOVERED** |
| ZITSIE | ✓ | ✗ | **✓** | 188.9 | 3.9 | 210.0 | **RECOVERED** |
| **MEDZUR** | ✓ | ✗ | **✗** | 8.9 | 1.0 | 7.4 | still broken |
| **POVPIA** | ✓ | ✗ | **✗** | 21.2 | 5.4 | 10.5 | still broken |

**6 of 8 recovered.** On all six the check returns arm A's *exact* structure — same
`sha256(smiles_2)`, same `clash_vdw`, same `worst_overlap` to four decimals — at approximately
arm A's runtime.

### 3.2 ⚠ Production recovery is 6/8, not the 7/8 the falsification measured

The offline predicate separates 7 of 8 (§1.3). Production recovers 6. **The two numbers answer
different questions and both belong in the record.**

§1.3 asks *can the predicate tell arm A's conformer from arm B's?* §3.1 asks *does the pipeline
end up returning a good one?* A filter cannot manufacture a conformer that does not exist in
the pool. `MEDZUR_comp_0` is the gap: the check correctly rejected arm B's conformer, the pool
was then filled, and selection returned a **third** structure (`worst_overlap` 0.7611 — neither
arm A's 0.7577 nor arm B's 0.7705) which **passes the attachment check and still fails
independent re-perception.**

### 3.3 Three residual mechanisms, read off arm C's own returned conformers

Scoring the check against what arm C actually returned (`--dump-xyz`, then `attach_probe`)
separates the residual cleanly. Only the first was anticipated.

| molecule | check on C's structure | `indep` | mechanism |
|---|---|---|---|
| POVPIA | passes | ✗ | **Known and predicted (§6.4).** Metal sphere intact; the defect is ligand-internal — a hydrogen detaches and C–N reads as C=N. Unreachable by any metal-centred predicate. |
| MEDZUR | **passes** | ✗ | **New class, not predicted by anyone.** Attachment fully intact and independent re-perception still disagrees. Here the attachment check is simply not the binding constraint. |
| GAVSED | **fails** | ✗ | **The scope limit (§2.1), demonstrated.** Acceptance rejected every conformer, so `_select_by_geometry`'s geometry-ranked fallback returned one regardless — **and that fallback is not attachment-aware.** |

GAVSED is the one to act on: the check guards **acceptance**, and on molecules where it rejects
everything, the structure that ships was chosen by a ranking that never consults it. That is a
named, fixable gap and it is **not** fixed here, because closing it would change arm A's
behaviour too and that is outside this lane's scope.

### 3.4 The cost — and this is where the case weakens

**Runtime is ADVISORY throughout: load averaged 42–65 on 12 cores, with a 5k sweep and two
sibling lanes running.** Totals within one cohort are comparable because both arms met the same
conditions; absolute seconds are not portable.

| gap cohort, n=21 | A default | B bare lever | **C lever + check** |
|---|---|---|---|
| total_s | 1788.9 | **220.6** | 1273.5 |
| median_s | 21.17 | **3.84** | 13.08 |
| `>30s` | 10 | **2** | 7 |
| speedup vs A | — | **8.11×** | **1.40×** |

**The check retains 32.9% of the bare lever's saving on this cohort** (B saves 1568 s against
A; C saves 515 s). The mechanism is not subtle and it is visible per molecule above: when the
check rejects the conformers the score would have accepted, `accept_fn` never fires, the pool
fills to completion, and **pool filling is where arm A's cost lives.** So on precisely the
molecules the lever was fastest on — KAQDOV 7.1 s → 248.9 s, ZITSIE 3.9 s → 210.0 s — the check
hands the runtime straight back.

**This cohort was selected to exhibit the acceptance gap**, so it is the worst case for speedup
retention by construction, in the same way it was the wrong cohort to estimate G1 from (§4.8).
The population number is §3.5 and it is the one to quote.

### 3.5 Guard population (n=100) — **G1 stops failing. G2 does not.**

The gate that inverted the promote lane's verdict (§4.8), re-run. Arm C compared against that
lane's stored arms A and B, same cohort, same `--timeout 150 --hard-cap 240`. `sha_in` control
OK on every molecule.

| gate | A → B (bare lever) | **A → C (lever + check)** | verdict |
|---|---|---|---|
| **G3** | 98/98 identical, 0 divergent | **97/97 identical, 0 divergent** | **PASS** |
| **G1** | clash 77→**105**, severe 5→**14**, worse 23 / better 13 | **clash 77→82, severe 5→5, worst_overlap min unmoved at 0.5101, worse 8 / better 5 / identical 84** | **no longer fails** |
| **G2** | 88/98 → 62/98, **26 regressions**, 0 fixes | **87/97 → 78/97, 9 regressions, 0 fixes** | **still fails** |
| **G4** | 98/100 → 99/100, **0 regressions** | 98/100 → 98/100, **1 regression** (`DURPAH`), 1 fix (`CUCBUZ`) | **marginal** |

Against the bare lever directly, arm C is a strict improvement on every structural axis:
**`indep` 61/98 → 78/98 — 17 fixes, 0 regressions**; **clash 106 → 83**; **severe 14 → 5**;
per-molecule **better 15 / worse 9 / identical 74**.

**Recovery: 17 of the 26 molecules the bare lever broke.** 65%, against 6/8 on the gap cohort.

#### The sharpest cut, in the same form §4.8 used

The lever can only cost something where it changes the returned conformer.

| | bare lever | **lever + check** |
|---|---|---|
| molecules whose conformer changed vs arm A | **35** of 98 | **13** of 98 |
| among those: `clash_vdw` | 13 → **41** (×3.2) | 1 → **6** |
| among those: `clash_severe` | 0 → **9** | **0 → 0** |
| among those: independent re-perception lost | **26 / 35 (74%)** | **9 / 13 (69%)** |

**This is the honest reading and it cuts both ways.** The check does not make the lever's
divergences safer — conditioned on still changing the answer, the damage rate is essentially
unchanged (74% → 69%). What it does is **cut the number of divergences from 35 to 13**. The
population-level G2 rate falls from **29.5% to 11.4%** because the check stops the lever acting
on 22 molecules where it would have acted badly, not because the 13 remaining divergences are
any better than before.

#### G4's single regression is timeout-shaped, and it is real

`DURPAH_comp_0`: arm A passes in 198.0 s, arm B in 112.5 s, **arm C is SIGKILLed at the 240 s
hard cap.** It is not a wrong answer — it is a molecule pushed past the budget by the extra pool
filling. This project has a documented history of reading timeout-shaped deltas as correctness
deltas (v0.4.4: 11 "regressions", all 300 s timeouts, 0 correctness), so the shape is named
explicitly. It is still a molecule that used to pass and no longer does, at a fixed budget.

#### Cost at population scale

| guard population, n=100 | A default | B bare lever | **C lever + check** |
|---|---|---|---|
| total_s | 2480.1 | **745.0** | 1734.7 |
| median_s | 10.51 | **4.70** | 7.23 |
| `>30s` | 20 | **1** | **14** |
| speedup vs A | — | **3.33×** | **1.43×** |

**43.0% of the bare lever's saving survives** — better than the gap cohort's 32.9%, exactly as
expected since the lever is a no-op on most of the population. But `>30s` goes **1 → 14**, which
matters more than the totals for a sweep with a per-molecule budget: that is the mechanism behind
`DURPAH` and it will recur at any cap.

### 3.6 L1's eta-concentrated frozen cohort — **ATTEMPTED, NOT COMPLETED, and the reason is a result**

§6.5 directed this test at L1's frozen slow-100 (`MANIFEST_SHA256=6f61359b…ab794`, `#DONE 100`,
**eta 60/100**) because the attachment claim is specifically about haptic coordination and that
cohort concentrates it.

**It could not be completed in both arms, and the reason is the finding itself.** Arm C on the
slow-100 completed **1 molecule in the time arm B completed 3**, projecting to roughly **25 h at
this machine's load**. Restricted to the first 24 eta molecules of the frozen list, arm C still
completed only 2 of 24 before it had to be stopped to let the population gate finish. **A safety
condition whose cost makes the eta cohort unaffordable to measure is itself evidence about the
speedup question**, and it points the same way as §3.4.

What the two completed pairs show — reported as anecdote, with n=2 stated plainly:

| molecule | B `indep` | **C `indep`** | B s | **C s** | B clash/worst | C clash/worst |
|---|---|---|---|---|---|---|
| ADANEB | ✗ | ✗ | 4.9 | **81.8** (17×) | 0 / 0.8155 | 0 / 0.7519 |
| AGIKUW | ✗ | **✓** | 16.5 | **72.1** (4.4×) | 0 / 0.7563 | **1 / 0.6883** |

One recovery, one not; both several-fold slower; and on the recovered one the vdW numbers got
**worse** (0 → 1 clash, worst overlap 0.7563 → 0.6883), which is a reminder that `indep` and
clash are genuinely different axes (the promote lane's §4.6 named finding).

> **A separate observation from the same run, and it needs its caveat stated louder than the
> number:** the bare lever failed independent re-perception on **11 of 11** completed slow-eta
> molecules. That is far worse than the guard population's rate — but **arm A was never measured
> on this cohort**, and the promote lane's §3 already established that the *default* path also
> ships `indep` failures on hard molecules. So this is **not** attributable to the lever. It is a
> flag that the eta tail may be in poor shape generally, and it is not evidence about
> `OIN_ACCEPT_SCORED` either way.

---

## 4. SCORECARD — the bar, and whether it was met

The brief set the bar explicitly: **"G2 and G1 must stop failing while G3 keeps passing and most
of the speedup survives."**

| criterion | result | met? |
|---|---|---|
| **G1 stops failing** | clash 77→82 (was 77→105), severe 5→5 (was 5→14), worst_overlap min unmoved | **YES** |
| **G2 stops failing** | 9 regressions, 0 fixes (was 26, 0). Population rate 29.5% → **11.4%** | **NO** — reduced by 61%, not eliminated |
| **G3 keeps passing** | 97/97 byte-identical, 0 divergent; 117/117 across both cohorts | **YES** |
| **most of the speedup survives** | **43.0%** retained; 3.33× → **1.43×**; `>30s` 1 → 14 | **NO** |
| *(G4, not in the bar)* | 1 pass regression (`DURPAH`, timeout-shaped), 1 fix | marginal |

**Two of four criteria met.** The check is a large, one-directional improvement on the bare
lever — 17 `indep` fixes and 0 regressions against it, clashes and severe clashes both down —
and it does not reach the bar that was set for promotion.

### 4.1 This is neither a clean win nor a null result, and it should not be rounded to either

- **Not a null result.** The check recovers 17 of 26 population regressions and 6 of 8 on the
  gap cohort, returning arm A's *exact* structure on the recovered ones. G1's population failure
  is genuinely erased. Reporting this as "the check does nothing" would be false — and is
  precisely what the first, broken implementation reported (§2.2).
- **Not a clean win.** G2 still fails, still perfectly one-way (0 fixes against arm A in either
  cohort, across 122 molecules and now three arms). More than half the speedup is gone, `>30s`
  goes 1 → 14, and one currently-passing molecule stops passing at a fixed budget.

### 4.2 What the residual costs, mechanically

The 9 remaining population regressions and the 2 on the gap cohort are not one problem:

1. **Ligand-internal defects** (`POVPIA`) — invisible to any metal-centred predicate. Named in
   advance by §6.4, confirmed here.
2. **Attachment intact, `indep` still disagrees** (`MEDZUR`) — the check is not the binding
   constraint. **This class was not predicted by anyone** and its size is unmeasured.
3. **The unguarded fallback** (`GAVSED`) — when acceptance rejects every conformer, the pool
   fills and `_select_by_geometry` returns a geometry-ranked structure that **never consults the
   check**. Demonstrated: GAVSED's returned conformer *fails* the attachment check and shipped
   anyway.

Class 3 is the only one with an obvious next move, and it is deliberately **not** taken here:
applying the check inside `_select_by_geometry`'s final ranking would change **arm A's** behaviour
too, which is a different change with its own gate, outside this lane's scope.

---

## 5. RECOMMENDATION

### 5.1 On the combined lever: **do not promote to default-ON. Keep both opt-in.**

Two of the four bar criteria are unmet, and the two that fail are the two the promotion case
needed: G2 still regresses independent re-perception on 11.4% of the molecules that survive it
with the lever off, and only 43% of the speedup — the lever's entire reason to exist — survives.
A lever that is 1.43× faster than the default and still loses structure on one molecule in nine
is not obviously worth a default-path behaviour change, and the cost remains **invisible to the
metric that would police it** (`passed` moved by exactly one molecule, in the wrong direction).

`OIN_ATTACH_CHECK` is registered in `_HELD_OFF` with this evidence.

### 5.2 On the pairing, which is a separate and clearer call: **never run `OIN_ACCEPT_SCORED` without it**

This is the one unambiguous finding of the lane. Against the bare lever the check is
**one-directional**: 17 `indep` fixes and **0** regressions, clash 106→83, severe 14→5, byte
identity untouched. Anyone who opts into `OIN_ACCEPT_SCORED` today for throughput work is taking
a structure cost they could largely avoid for 24 ms a conformer and roughly half the speedup.

**Suggested disposition:** make `OIN_ATTACH_CHECK` default-ON *whenever `OIN_ACCEPT_SCORED` is
on* — i.e. the check becomes part of what that lever means — while both stay off by default. That
is a smaller decision than promotion and the evidence for it is not mixed.

### 5.3 The product call, stated so it can be made by someone else

**The choice is not "is the check good" — it is what the round trip is for**, which is the same
fork §5.3 of the promote lane identified, now with the check priced.

- **If the OIN string is the product** (Reading A), the bare lever was already adequate — G3 is
  117/117 byte-identical with or without the check — and the check buys nothing you value while
  costing 57% of the speedup. **Then: ship the bare lever opt-in, skip the check.**
- **If the round trip must reproduce the geometry** (Reading B), the check is what makes the
  lever defensible at all: it erases G1's population failure and recovers two thirds of G2. But
  it does not finish the job, and at 1.43× the default the remaining speedup may not justify a
  second lever. **Then: keep both opt-in and fix the fallback (§4.2 class 3) before revisiting.**

**I am not making this call.** Both readings are priced above with the same measurements.

### 5.4 What would change the recommendation

1. **Close the fallback gap (§4.2 class 3).** `_select_by_geometry` returns structures the check
   rejects. Guarding the *return* rather than only *acceptance* is a bounded change with a
   measurable target, and GAVSED is the fixture.
2. **Size class 2.** `MEDZUR` passes the attachment check and still fails `indep`. Nobody knows
   how many molecules are in that class; until someone does, the ceiling on this approach is
   unknown rather than "7/8".
3. **Recover the speedup.** The cost is entirely pool-filling: rejecting the first conformer
   means paying arm A's embed budget. A cheaper path would be to make the check inform
   *generation* (place ligands that stay attached) rather than filter after the fact.
4. **Re-run on L1's slow-100 in both arms** once (1) or (3) makes arm C affordable there. §3.6
   could not complete, and the eta tail is where both the speedup and the damage live.

---

## 6. Reproduce

```
# falsification (§1): dump both arms' accepted conformers, then score offline
tools/ab_accept_scored.py --cohort spec/handoffs/v0.4.7/cohort_attach21.json \
  --out spec/handoffs/v0.4.7/runs/attach21.json \
  --dump-xyz spec/handoffs/v0.4.7/dump21 --timeout 300 --hard-cap 500 --workers 2
tools/attach_probe.py --dump spec/handoffs/v0.4.7/dump21 --json runs/attach_probe21.json
tools/attach_probe.py --inputs spec/handoffs/v0.4.7/cohort_attach21.json   # crystal-input control

# arm C (§3): the lever PLUS the check, as a single arm
OIN_ATTACH_CHECK=1 tools/ab_accept_scored.py --cohort <cohort>.json \
  --single-arm 1 --label C-lever+check --out runs/c<n>.json

# the scorecard, no re-running
tools/attach_gate_report.py --baseline runs/attach21.json --armc runs/c21.json
tools/attach_gate_report.py --baseline runs/guard100.json --armc runs/cguard100.json
```

All runs under `systemd-run --user -p OOMPolicy=continue -p MemoryMax=6G`, 1–2 workers, on a box
at load 40–65 of 12 cores. **Every second in this document is ADVISORY**; within-cohort ratios
are meaningful because both arms met the same conditions, absolute values are not portable.

Artefacts: `spec/handoffs/v0.4.7/runs/` — `attach21.json`, `c21.json`, `cguard100.json`,
`attach_probe21.json`, `inputs_control.json`, `gates_c21.txt`, `gates_guard100.txt`, and
`INVALID_c21_noop.json` (§2.2, **not** used for any number here).
