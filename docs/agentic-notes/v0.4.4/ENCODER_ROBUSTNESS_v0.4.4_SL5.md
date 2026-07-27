# SL5 — Encoder Robustness (v0.4.4)

The v0.4.2 capstone (6,719 molecules) had **48 `encode_fail` molecules** (0.71%): a crystal
XYZ that produced no OIN string at all (`smiles_1 is None`). They can never round-trip, so they
are a hard ceiling on the 99% goal. SL5 sub-triaged the cohort, classified the irreducible part,
and hardened the two perception pathologies that hang the rest — all while keeping every
currently-encodable structure byte-identical (the hard constraint: the encoder feeds the whole
project).

## Sub-triage of the 48 (`tools/sl5_triage.py`)

Each molecule was reproduced through `XYZToSMILES().convert()` in an isolated subprocess with an
OS-level timeout (a Python `signal.alarm` cannot interrupt the C-level hangs these inputs cause).

| bucket | count | root cause |
|---|---:|---|
| `boron_cluster` | **34** | electron-deficient carborane / closo-nido borane cage; RDKit's 2-center-2-electron valence model has no Lewis structure for a 3c-2e cage |
| `resonance_timeout` | **10** | large conjugated ligand hangs in xyz2mol perception |
| aromatic-perception | 3 | quinoid/ylide ring stays invalid after de-aromatizing (`KAXVOX`, `KAXWAK`, `LEZWAO`) |
| `perception_charge_gap` | 1 | `ASISAX`: fused azacage, charge sweep exhausted |

## W1 — the boron-cluster ceiling is classified, not a crash

RDKit cannot perceive a 3c-2e boron cage into a sanitizable Lewis structure, and `get_lig_mol`'s
charge sweep already spans −4..+4, so no charge widening will ever encode these. The encoder now
detects the cage (`_is_electron_deficient_cluster`: ≥3 borons **and** a B–B bond — a `BPh4⁻` borate
has neither) and raises a typed **`OINEncodeError`** naming the cause; `core/translator.py::convert`
re-raises the typed error rather than flattening it, so callers can distinguish a known ceiling from
an unexpected failure. `OINEncodeError` subclasses `ValueError`, so existing handlers are unaffected.
This classifies **34/48** of the cohort (it does not *encode* boron cages — out of scope).

## W3 — the timeout cohort: three hang stages, two bounded

The timeout cohort hangs at **three** independent super-polynomial stages of xyz2mol perception,
each revealed only after bounding the previous (the handoff hypothesis named only the second):

### 1. `AC2BO` valence-order sort (`xyz2mol_local.py`) — bounded, byte-identical
It materialised the full Cartesian product of per-atom valences to sort candidate assignments —
exponential for dozens of multivalent atoms. `_VALENCE_COMBO_CAP` (500k) skips the sort above the
cap and iterates the lazy product bounded to `_VALENCE_FALLBACK_TRIES`; the main loop early-returns
on the first valid assignment. **Provably byte-identical:** a >500k-combo ligand *hangs* on
unbounded `main`, so no currently-encodable molecule reaches the fallback.

### 2. `ResonanceMolSupplier` (`lig_checks`) — forked, CPU-time-bounded (recovers the cohort)
It builds the conjugation-electron groups in a C++ call *before* any enumeration, so `maxStructs`
does not bound it (even `maxStructs=2` hangs), and it explodes on a large conjugated ligand. It
**holds the GIL** inconsistently — a watchdog thread cannot reliably interrupt it — so the
enumeration runs in a **forked child bounded by a CPU-time budget**
(`_resonance_candidates_isolated`, `os.fork` + pipe + `select`, child `RLIMIT_CPU =
_RESONANCE_CPU_BUDGET_S = 120` CPU-seconds):
- **Finishes within the budget →** the child returns the resonance forms (each a property-preserving
  `ToBinary`), reconstructed in the parent. **Byte-identical to the inline path** — verified: e.g.
  `EHADAV` (Co macrocycle, ~21 CPU-s) and `RUTJEW` (72 aromatic atoms) reproduce `main`'s OIN exactly.
- **Burns past the budget →** the kernel kills the child (`SIGXCPU`, enforced even mid-C-call) and
  perception falls back to the single form, **recovering** the otherwise-unencodable molecule
  (e.g. `BENVOG`, `HUCNAU`).

**Why CPU-time, not wall-clock:** a wall-clock timeout would make the outcome depend on machine load
(a starved *completer* could wrongly fall back, changing its OIN). CPU-time is load-independent — a
given ligand's resonance burns ~constant CPU seconds whether the box is idle or saturated — so the
encode stays deterministic. A higher budget only lets *more* completers finish (hangs burn CPU
without bound and are always killed), so the budget is set generously above every observed
completer. `os.fork` (not `multiprocessing`) is deliberate: it works inside the dataset harness's
*daemon* workers (which forbid child processes) and inherits the ligand copy-on-write (no input
pickling). Only large ligands (`_resonance_needs_isolation`: ≥50 heavy or ≥35 aromatic) take this
path; ordinary ligands run inline unchanged. Routing a *completer* here is byte-identical, so the
size gate only trades a fork for a safety net — it never changes an encodable molecule's OIN.

This is what a naive size-based *skip* (tried first, then withdrawn) could not do: that skip changed
the perceived form of passing molecules whose resonance completes quickly (e.g. `EHADAV` regressed —
resonance finds the aromatic form, the skip kept the localized tautomer). The forked bound keeps the
completer's real resonance result and falls back only on a *genuine* hang.

### 3. `get_UA_pairs` → networkx `max_weight_matching` (`xyz2mol_local.py`) — not bounded
O(V³) in the unsaturated-atom graph; dominates on a few large charged conjugated ligands (`FAQYUU`,
`HICLAG`). Some of the cohort hang here rather than in resonance and are not recovered — documented
residual, a candidate for the same fork-timeout treatment if warranted.

## Byte-identity (the hard constraint) — verified

- **Goldens** — `test_regression_stability` passes unchanged.
- **Encoder regression suite** — the 9 XYZ→OIN test modules run **66/66 OK**, including the new
  `test_encoder_robustness` (boron typing + detector, AC2BO cap, forked-resonance recovery).
- **Firing-set A/B** — `tools/sl5_encode_dump.py` / per-molecule `main`-vs-branch on the ligands
  that route through the fork confirms completers are byte-identical (`EHADAV`, `BUJMEZ`, …); the
  only differences are molecules current `main` can no longer encode at all (`NUKXEY`), which the
  branch *recovers*.
- The generation path (`generator3d/`, `generation/`) is untouched.

## Net effect

- **34/48** encode-fails now fail with a **typed, classified** `OINEncodeError` (W1).
- The resonance-hang subset is **recovered byte-identically** by the forked time-bounded resonance;
  a genuine hang now costs a bounded ~90s then degrades instead of hanging forever.
- Remaining: the `get_UA_pairs`-matching hangs, 3 quinoid/ylide de-aromatize cases (already emit a
  typed `OINEncodeError`), and `ASISAX`'s charge gap — documented future work.

Reproduce: `tools/sl5_triage.py` (cohort triage), `tools/sl5_encode_dump.py` / `sl5_reencode_check.py`
(byte-identity gates), `tools/sl5_profile_hang.py` (pinpoint a hang stage).
