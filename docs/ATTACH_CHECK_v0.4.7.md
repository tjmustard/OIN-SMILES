# `OIN_ATTACH_CHECK` — the attachment check for `OIN_ACCEPT_SCORED` (v0.4.7, lane L5-attach)

**Status: falsification PASSED, check BUILT and default-OFF, four gates RE-RUN. Recommendation in §5.**

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
