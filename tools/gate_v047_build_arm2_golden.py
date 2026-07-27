"""Build the ARM 2 (round-trip) frozen golden manifest for `tools/gate_v047.sh`.

The v0.4.7 slow-100 cohort (see `docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`) is selected from a SINGLE
sweep results dir (`tools/select_slow_byte_exact.py`, top-100 by `metrics.elapsed_s`
within that one dir) -- ``smiles_1``/``smiles_2`` are therefore ALREADY known and
byte-exact there. So the golden values this tool freezes come from that EXISTING
primary source -- no new generation is run here. Re-running the full round trip for
all 100 molecules today would cost the exact wall-clock this cohort was chosen to be
expensive at, which is the quiet-phase sweep's job, not this lane's.

PRIMARY SOURCE, NOT TWO PEERS
=============================
An earlier version of this tool required agreement between TWO sources, treated as
equal peers, and hard-failed on any disagreement or absence in either. That was
wrong: the two source dirs are different, uncorrelated random draws over the
dataset (see ``docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`` -- the "62, not 100" incident) and share only
~100 names total, so requiring both was a near-total sample destruction, not a
robustness check. This tool now takes ONE required primary source
(``--source-a``, the same dir the cohort was selected from -- every cohort name
MUST resolve there and pass the full predicate, or the build fails loudly) and an
OPTIONAL corroboration source (``--source-b``) that is purely informational: where
a name happens to also have a report there, its ``smiles_1`` is cross-checked as an
encoder-determinism sanity check (mismatch => loud warning, never fatal), and
absence is expected and unremarkable.

For each cohort molecule this tool:
  1. Loads its report from the primary source and re-verifies the full selection
     predicate there (status/tier/byte-exact/elapsed) -- a re-check, not a re-trust.
  2. If a corroboration source is given and has a report for this name with a
     non-null ``smiles_1``, compares it to the primary's -- logged, never fatal.
  3. Emits ``name<TAB>sha256(smiles_1)<TAB>sha256(smiles_2)<TAB>len1<TAB>len2<TAB>eta``.

Sorted by name, followed by ``#DONE <n>``. Same TSV shape as ARM 1's manifest so both
arms are diffed the same way by ``gate_v047.sh``.

v0.4.9: ``--predicate any_encoded`` AND WHY A RUNTIME COHORT NEEDS IT
=====================================================================
The default ``--predicate byte_exact_scored`` is the v0.4.7 behaviour, unchanged and still
the default so ``tools/gate_v047_arm2_golden.tsv`` rebuilds bit-identically. It requires
``status=success``, ``tier_passed=UFF_1`` and ``smiles_1 == smiles_2``.

That predicate cannot express the v0.4.9 benchmark, for two independent reasons:

1. **A runtime benchmark must contain the molecules that FAIL.** They are the slow ones --
   the >=300 s band is 291 molecules burning 49.8% of the sweep's entire compute for three
   honest passes. Requiring ``status=success`` deletes exactly the population Goal B is about.
2. **``smiles_1 == smiles_2`` is the SCORED verdict**, which v0.4.8 showed over-states
   ``byte_exact`` by 10.34 points corpus-wide and much more in the slow tail. Of the 100
   molecules the v0.4.7 cohort froze under it, **59 are ``byte->FAIL``** once the round trip
   is re-perceived from coordinates. Selecting on it again would rebuild the same mistake.

``any_encoded`` therefore requires only a report with a non-null ``smiles_1``: the encode
happened, so there is something to freeze.

``NO_STRUCTURE`` -- THE SENTINEL THAT KEEPS THIS GATE FROM BEING A CLOCK
=========================================================================
For a molecule whose source run produced no ``smiles_2``, freezing an empty hash would make
the gate assert "generation fails here". **That is not a property of the code.** Those runs
ended at a wall-clock budget, so the assertion is really about how fast the box was --
``ULODUU`` assembles at a 60 s cap and not at 30 s, which is the whole reason the boron
fast-fail was refuted. A golden built that way turns every future speedup into a MISMATCH.

So column 3 becomes ``NO_STRUCTURE@<budget>s``, carrying the budget it was observed under
(that is what makes the observation interpretable at all), and ``gate_v047.sh`` gates those
rows on ``sha256(smiles_1)`` ONLY, reporting the fresh generation outcome as an observation.
Byte-identity is claimed exactly where byte-identity is a stable property.

``--strata-file`` appends two OBSERVATION columns, ``band`` and ``honest_class``, from
``tools/select_runtime_strata.py``'s manifest. They are never compared. They exist so nobody
reads a green ARM 2 as "the chemistry is right" -- 71 of the 328 rows in the v0.4.9 cohort
are ``byte->FAIL``, and the gate is deliberately blind to that.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/gate_v047_build_arm2_golden.py \
        --names-file cohort_names.txt \
        --source-a <DS>/results-v0.4.5-sweep-partial-2697mols \
        --source-b <DS>/results-v0.4.5-rebaseline \
        --out tools/gate_v047_arm2_golden.tsv

    # v0.4.9 stratified runtime cohort
    PYTHONPATH=src .venv/bin/python tools/gate_v047_build_arm2_golden.py \
        --names-file spec/handoffs/v0.4.9/cohort_v049_names.txt \
        --source-a <DS>/results-v0.4.8-honest \
        --predicate any_encoded --source-budget-s 300 \
        --strata-file spec/handoffs/v0.4.9/cohort_v049_strata.json \
        --allow-unbaselined \
        --out tools/gate_v049_arm2_golden.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys

HAPTIC = re.compile(r"\{\d+[<>]\}")


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def load_report(results_dir: str, name: str) -> dict | None:
    path = os.path.join(results_dir, "individual_reports", f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def read_names(names_file: str) -> list[str]:
    names = []
    with open(names_file) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            names.append(line[: -len(".xyz")] if line.endswith(".xyz") else line)
    return sorted(set(names))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--names-file", required=True, help="one basename (no .xyz) or name.xyz per line"
    )
    ap.add_argument(
        "--source-a",
        required=True,
        help="PRIMARY source -- the dir the cohort was selected from; must pass the full "
        "predicate for every name, or the build fails",
    )
    ap.add_argument(
        "--source-b",
        default=None,
        help="OPTIONAL corroboration source -- informational only, never required, never "
        "blocks the build on absence or mismatch (mismatch is warned loudly instead)",
    )
    ap.add_argument("--out", default=None, help="also write the TSV to this path")
    ap.add_argument(
        "--predicate",
        choices=("byte_exact_scored", "any_encoded"),
        default="byte_exact_scored",
        help="byte_exact_scored (default, v0.4.7 behaviour): status=success, tier UFF_1, "
        "smiles_1 == smiles_2. any_encoded (v0.4.9 runtime cohort): only a non-null "
        "smiles_1 is required, so molecules that FAILED are kept -- they are the slow ones.",
    )
    ap.add_argument(
        "--source-budget-s",
        type=float,
        default=300.0,
        help="the per-attempt wall-clock budget the SOURCE sweep ran under. Written into the "
        "NO_STRUCTURE sentinel, because 'this molecule produced nothing' is only "
        "interpretable alongside the budget it was given (ULODUU assembles at 60 s, not 30).",
    )
    ap.add_argument(
        "--strata-file",
        default=None,
        help="tools/select_runtime_strata.py manifest; appends band + honest_class as "
        "OBSERVATION columns (never compared by the gate)",
    )
    ap.add_argument(
        "--allow-unbaselined",
        action="store_true",
        help="do not fail on cohort names absent from the primary source. They are REPORTED "
        "and EXCLUDED from the golden -- a name with no source row has nothing to freeze, "
        "and inventing a row for it would be the fabrication this flag exists to avoid.",
    )
    args = ap.parse_args()

    names = read_names(args.names_file)
    if not names:
        sys.exit(f"error: 0 names read from {args.names_file}")

    strata = {}
    if args.strata_file:
        with open(args.strata_file) as f:
            strata = {r["name"]: r for r in json.load(f)["molecules"]}

    lines = []
    missing_primary = []
    bad_primary = []
    no_structure = []
    no_structure_det = []
    no_encode = []
    corroborated = 0
    corroboration_mismatches = []
    for name in names:
        ra = load_report(args.source_a, name)
        if ra is None:
            missing_primary.append(name)
            continue
        s1a, s2a = ra.get("smiles_1"), ra.get("smiles_2")
        if args.predicate == "byte_exact_scored":
            if ra.get("status") != "success" or ra.get("tier_passed") != "UFF_1":
                bad_primary.append((name, ra.get("status"), ra.get("tier_passed")))
                continue
            if not s1a or not s2a or s1a.strip() != s2a.strip():
                bad_primary.append((name, "smiles_1/smiles_2 missing or not byte-exact"))
                continue
            oin1, oin2 = s1a.strip(), s2a.strip()
        else:  # any_encoded
            oin1 = s1a.strip() if s1a else None
            oin2 = s2a.strip() if s2a else None
            if oin1 is None:
                # Not an encode *failure* -- an encode TIMEOUT. All 15 such rows in the 5k
                # sweep read "TimeoutException at UFF_1: exceeded 300s while encoding (hard
                # kill)", 1.25 CPU-h of it. That is budget-dependent exactly like
                # NO_STRUCTURE, so it gets a sentinel rather than a frozen hash -- and it is
                # the one population NO generator-side bound can reach, because it never gets
                # as far as the generator. Worth stating loudly in Lane 1's scope.
                no_encode.append(name)

        if args.source_b:
            rb = load_report(args.source_b, name)
            if rb is not None:
                s1b = rb.get("smiles_1")
                if s1b:
                    if s1b.strip() == oin1:
                        corroborated += 1
                    else:
                        corroboration_mismatches.append(name)
                        print(
                            f"warning: {name} smiles_1 DIFFERS in corroboration source "
                            f"{args.source_b} -- possible encoder non-determinism, "
                            f"investigate (NOT blocking the golden build)",
                            file=sys.stderr,
                        )

        if oin1 is None:
            sentinel = f"NO_ENCODE@{args.source_budget_s:g}s"
            row = f"{name}\t{sentinel}\t{sentinel}\t-\t-\t-"
            if strata:
                st = strata.get(name, {})
                row += f"\t{st.get('band') or '-'}\t{st.get('honest_class') or '-'}"
            lines.append(row)
            continue

        eta = "eta" if HAPTIC.search(oin1) else "-"
        if oin2 is None:
            # See the module docstring: an empty hash would assert "generation fails here",
            # which for a run that ended at a wall-clock cap is a fact about the budget.
            #
            # But only 45 of the 50 such rows in the v0.4.9 cohort ended that way. The other
            # 5 raised a real error (GURKUA in 1.2 s on an UncoordinatedFragmentError), and
            # THAT is a code-determined property worth gating strictly -- a lane that
            # accidentally makes one of them generate should be caught, not shrugged at. The
            # source report's own error text separates the two, so use it rather than
            # flattening both into the weaker sentinel.
            if "TimeoutException" in (ra.get("error") or ""):
                sha2 = f"NO_STRUCTURE@{args.source_budget_s:g}s"
                no_structure.append(name)
            else:
                sha2 = "NO_STRUCTURE_DET"
                no_structure_det.append(name)
            len2 = "-"
        else:
            sha2 = sha(oin2)
            len2 = str(len(oin2))
        row = f"{name}\t{sha(oin1)}\t{sha2}\t{len(oin1)}\t{len2}\t{eta}"
        if strata:
            st = strata.get(name, {})
            row += f"\t{st.get('band') or '-'}\t{st.get('honest_class') or '-'}"
        lines.append(row)

    if missing_primary and not args.allow_unbaselined:
        sys.exit(
            f"error: {len(missing_primary)} / {len(names)} names have NO report in the "
            f"primary source {args.source_a} -- every cohort name must resolve there "
            f"(that is where it was selected from): {missing_primary}"
        )
    if missing_primary:
        print(
            f"# UNBASELINED (--allow-unbaselined): {len(missing_primary)} cohort name(s) have "
            f"no report in {args.source_a} and are EXCLUDED from the golden -- they are still "
            f"in the cohort dir and will be measured, they just have nothing frozen yet: "
            f"{sorted(missing_primary)}",
            file=sys.stderr,
        )
    if bad_primary:
        sys.exit(
            f"error: {len(bad_primary)} / {len(names)} names failed the full predicate "
            f"in the primary source (should be impossible -- re-selection would have "
            f"excluded them): {bad_primary}"
        )
    if len(lines) == 0:
        sys.exit("error: 0 golden rows produced -- refusing to write an empty golden manifest")

    print(
        f"# primary={args.source_a} predicate={args.predicate} names={len(names)} "
        f"rows={len(lines)} no_structure={len(no_structure)} "
        f"no_structure_det={len(no_structure_det)} no_encode={len(no_encode)} "
        f"corroborated={corroborated} corroboration_mismatches={len(corroboration_mismatches)}",
        file=sys.stderr,
    )
    if no_structure:
        print(
            f"# {len(no_structure)} row(s) frozen as NO_STRUCTURE@{args.source_budget_s:g}s "
            "(the source run TIMED OUT) -- gated on sha256(smiles_1) only; their generation "
            "outcome is an observation, not a byte-identity claim",
            file=sys.stderr,
        )
    if no_structure_det:
        print(
            f"# {len(no_structure_det)} row(s) frozen as NO_STRUCTURE_DET (the source run "
            "raised, it did not time out) -- gated STRICTLY: producing a structure now is a "
            f"MISMATCH: {sorted(no_structure_det)}",
            file=sys.stderr,
        )
    if no_encode:
        print(
            f"# {len(no_encode)} row(s) frozen as NO_ENCODE@{args.source_budget_s:g}s -- the "
            "encoder itself was hard-killed at the budget, so these rows carry NO "
            f"byte-identity signal at all: {sorted(no_encode)}",
            file=sys.stderr,
        )

    manifest = "\n".join(lines)
    out_text = manifest + f"\n# MANIFEST_SHA256={sha(manifest)}\n#DONE {len(lines)}"
    print(out_text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(out_text + "\n")
        print(f"# wrote {os.path.abspath(args.out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
