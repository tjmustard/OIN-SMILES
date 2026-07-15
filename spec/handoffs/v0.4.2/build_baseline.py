#!/usr/bin/env python3
"""P0: build the v0.4.2 clean-floor BASELINE from existing c7edeeb6-stamped reports.

No new sweep: the live accumulator has already produced a single-commit (c7edeeb6)
sample. We extract (a) the trusted must-not-regress passing set, (b) per-class
goldens confirmed-failing on the exact baseline commit. Zero contention, no pause.
"""
import glob
import json
import os
import re
import sys

ROOT = "/home/tjmustard/Documents/GitHub/OIN-SMILES"
REPORTS = f"{ROOT}/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/individual_reports"
OUTDIR = f"{ROOT}/spec/handoffs/v0.4.2"
BASE = "c7edeeb6"

sys.path.insert(0, f"{ROOT}/tools")
import classify_failures as CF  # noqa: E402

overrides = json.load(open(f"{ROOT}/tools/triage_overrides.json"))
overrides = {k: v for k, v in overrides.items() if not k.startswith("_")}

FIXABLE = [
    "donor_H_atom_count", "H_on_terminal_oxo_imido", "geometry_NON",
    "geometry_or_fragment_change", "winding_flip", "EZ_bond_stereo", "atom_stereo",
    "encode_crash_other", "kekulize_encode_crash", "macrocycle_perception",
    "garbled_aromatic",
]
ARTIFACT = ["timeout", "high_rmsd", "carborane_unsupported", "no_conformers",
            "gen_exception_other", "rmsd_mapping_failed"]
# backlog-named representatives (first golden per class, if present on baseline)
REP = {
    "donor_H_atom_count": "AJIJUY_comp_0", "H_on_terminal_oxo_imido": "BADKAS_comp_0",
    "geometry_NON": "DEKQAN_comp_0", "geometry_or_fragment_change": "AHAQUX_comp_0",
    "winding_flip": "AGOVOK_comp_0", "EZ_bond_stereo": "AHAZOZ_comp_0",
    "atom_stereo": "CUQVUF_comp_0", "encode_crash_other": "CASHOW_comp_0",
    "kekulize_encode_crash": "HANWIE_comp_0", "macrocycle_perception": "FOSNEI_comp_0",
    "garbled_aromatic": "DIXXIS_comp_0",
}

S_STEREO = re.compile(r"\[S@")


def cls_of(rep):
    mol = rep.get("molecule")
    if mol in overrides:
        return overrides[mol]
    return CF.classify(rep)[0]


trusted_pass, untrusted_pass = [], []
fail_by_class = {}          # class -> list of (mol, rep)
s_stereo_rows = []          # [S@SP3] subset of string_mismatch_other on baseline

for f in glob.glob(f"{REPORTS}/*.json"):
    try:
        rep = json.load(open(f))
    except Exception:
        continue
    mol, st, cid = rep.get("molecule"), rep.get("status"), rep.get("commit_id")
    if st == "success":
        (trusted_pass if cid == BASE else untrusted_pass).append(mol)
        continue
    if st != "failed" or cid != BASE:
        continue
    c = cls_of(rep)
    fail_by_class.setdefault(c, []).append((mol, rep))
    if c == "string_mismatch_other":
        s1, s2 = rep.get("smiles_1") or "", rep.get("smiles_2") or ""
        if S_STEREO.search(s2) and not S_STEREO.search(s1):
            s_stereo_rows.append((mol, rep))

trusted_pass = sorted(set(trusted_pass))
untrusted_pass = sorted(set(untrusted_pass) - set(trusted_pass))

# sidecar: trusted passing IDs
with open(f"{OUTDIR}/baseline_pass_c7edeeb6.txt", "w") as fh:
    fh.write("\n".join(trusted_pass) + "\n")

# sidecar: per-class failing members + goldens signatures
fail_json = {}
for c, rows in fail_by_class.items():
    fail_json[c] = sorted(m for m, _ in rows)
with open(f"{OUTDIR}/baseline_fail_c7edeeb6.json", "w") as fh:
    json.dump(fail_json, fh, indent=1, sort_keys=True)


def goldens_for(c, rows, n=8):
    rows = sorted(rows, key=lambda x: x[0])
    rep_mol = REP.get(c)
    ordered = ([r for r in rows if r[0] == rep_mol] +
               [r for r in rows if r[0] != rep_mol])
    return ordered[:n]


def _collapse(s):
    return re.sub(r"\s+", " ", (s or "").strip())


def sig(rep):
    err = _collapse(rep.get("error"))
    if err:
        # drop python traceback tail; keep the message head
        err = err.split("Traceback (most recent call last)")[0].strip()
        return "error: " + err[:160]
    s1, s2 = _collapse(rep.get("smiles_1")), _collapse(rep.get("smiles_2"))
    return f"exp `{s1[:90]}` | got `{s2[:90]}`"


lines = []
w = lines.append
w("# v0.4.2 BASELINE — clean single-commit floor")
w("")
w(f"Baseline commit **`{BASE}`** (= tag `v0.4.1`). Built by P0 **without a new sweep**: the live "
  "accumulator had already produced a single-commit sample on this exact commit, so the floor is "
  "extracted directly (provenance filter + confirmed-on-baseline goldens). It grows as the "
  "accumulator runs — regenerate with `scratchpad/build_baseline.py`.")
w("")
w("## Must-not-regress passing set (trusted provenance)")
w("")
w(f"**{len(trusted_pass)} molecules** stamped `commit_id == {BASE}` and `status == success`. "
  "Full ID list: `spec/handoffs/v0.4.2/baseline_pass_c7edeeb6.txt`. The capstone's blocker gate is "
  "`{{passes on release/v0.4.2}} ⊇ {{this set}}`.")
w("")
w(f"**Untrusted passers ({len(untrusted_pass)}):** stamped on the dirty/older tree "
  "(`5538b722-dirty` etc.), not `c7edeeb6` — a pass there may not hold on the baseline. Treat as "
  "*unverified*; do not add to the must-not-regress set without a baseline re-run.")
w("")
w("## Per-class floor (failures confirmed on `c7edeeb6`)")
w("")
w("Counts are failures **already stamped `c7edeeb6`** (the clean floor sample so far); the wider "
  "mixed-provenance backlog is larger. Goldens are drawn from these baseline-stamped rows, so each "
  "reproduces on the exact baseline commit. Full per-class member lists: "
  "`spec/handoffs/v0.4.2/baseline_fail_c7edeeb6.json`.")
w("")
w("### Fixable classes (wave targets)")
w("")
for c in FIXABLE:
    rows = fail_by_class.get(c, [])
    w(f"#### `{c}` — {len(rows)} on baseline")
    if not rows:
        w("- (none stamped c7edeeb6 yet — pull from the growing backlog / re-run a golden)")
        w("")
        continue
    gs = goldens_for(c, rows)
    w(f"- goldens: {', '.join('`'+m+'`' for m, _ in gs)}")
    m0, r0 = gs[0]
    w(f"- repr `{m0}`: {sig(r0)}")
    w("")

# [S@SP3] subset
w(f"#### `string_mismatch_other` → `[S@SP3]` subset — {len(s_stereo_rows)} on baseline (S6b target)")
if s_stereo_rows:
    gs = sorted(s_stereo_rows, key=lambda x: x[0])[:8]
    w(f"- goldens: {', '.join('`'+m+'`' for m, _ in gs)}")
    w(f"- repr `{gs[0][0]}`: {sig(gs[0][1])}")
else:
    w("- (none isolated on baseline yet; S6b re-triages string_mismatch_other via triage_overrides)")
w("")
w("### Artifact classes — context only (NOT part of the floor; S7/docs own them)")
w("")
for c in ARTIFACT:
    w(f"- `{c}`: {len(fail_by_class.get(c, []))} on baseline")
w("")
other = sorted(set(fail_by_class) - set(FIXABLE) - set(ARTIFACT))
if other:
    w("### Other failed classes on baseline (triage)")
    for c in other:
        w(f"- `{c}`: {len(fail_by_class.get(c, []))}")
    w("")
w("## Method / caveats")
w("")
w("- **No headline percentage.** The floor is a *set of molecule IDs*, per the wave protocol.")
w("- Classes were derived with `tools/classify_failures.py::classify` + `tools/triage_overrides.json` "
  "(same routing the registry uses), over `results-v0.4.0/individual_reports/*.json` filtered to "
  f"`commit_id == {BASE}`.")
w("- Artifact classes (timeout/high_rmsd/carborane/no_conformers) are context only — highest cost, "
  "lowest diagnostic value; S7 triages them at full fidelity, docs records the residual.")

open(f"{OUTDIR}/BASELINE.md", "w").write("\n".join(lines) + "\n")

print("trusted_pass:", len(trusted_pass), "| untrusted_pass:", len(untrusted_pass))
print("baseline failed classes:")
for c in sorted(fail_by_class, key=lambda k: -len(fail_by_class[k])):
    tag = "FIX" if c in FIXABLE else ("ART" if c in ARTIFACT else "?")
    print(f"  [{tag}] {c:32s} {len(fail_by_class[c])}")
print("[S@SP3] subset:", len(s_stereo_rows))
print("wrote BASELINE.md + baseline_pass_c7edeeb6.txt + baseline_fail_c7edeeb6.json")
