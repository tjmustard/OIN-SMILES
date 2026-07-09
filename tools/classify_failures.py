"""Classify round-trip failures into defect classes and emit CASE_REGISTRY.md.

Reads ``individual_reports/*.json`` from a results directory, classifies every
row by error class plus smoking-gun string patterns, maps each class to the
worktree session that owns it (see ``docs/handoffs/``), and writes
``CASE_REGISTRY.md`` + ``case_registry.json`` into the results directory.

Stdlib-only on purpose: it must run anywhere without the chemistry stack.

Usage:
    python tools/classify_failures.py --output-dir tmCAT-tmPHOTO_xyz_dataset/20260707-results
"""

import argparse
import glob
import json
import os
import re
from collections import Counter, defaultdict
from datetime import datetime

SLOT_RE = re.compile(r"\{(\d+)([<>^]?)\}")

# defect class -> owning handoff session (docs/handoffs/S*.md)
SESSION_OF = {
    "donor_H_atom_count": "S1-donor-h",
    "H_on_terminal_oxo_imido": "S1-donor-h",
    "frag_vector_IndexError": "S1-donor-h",
    "eta_diene_localization": "S2-eta-diene",
    "garbled_aromatic": "S3-aromatic-perception",
    "kekulize_encode_crash": "S3-aromatic-perception",
    "xyz2mol_none_crash": "S3-aromatic-perception",
    "encode_crash_other": "S3-aromatic-perception",
    "macrocycle_perception": "S3-aromatic-perception",
    "winding_flip": "S4-eta-winding",
    "eta_slot_placement": "S4-eta-winding",
    "rmsd_sentinel": "S5-metrics",
    "high_rmsd": "S5-metrics",
    "geometry_or_fragment_change": "S5-metrics (triage)",
    "atom_stereo": "S6-stereo",
    "EZ_bond_stereo": "S6-stereo",
    "no_conformers": "S6-stereo (triage)",
    "geometry_NON": "S1-donor-h",
    "carborane_unsupported": "wontfix-docs",
    "timeout": "unassigned-triage",
    "stale_harness_TypeError": "rerun-only (fixed in v0.3.5)",
    "interrupted": "rerun-only",
    "junk_row": "rerun-only",
    "string_mismatch_other": "unassigned-triage",
    "gen_exception_other": "unassigned-triage",
}


def _metal_geo(rep):
    m = re.match(r"\[([A-Za-z]+)_([A-Z0-9]+)\]", rep.get("smiles_1") or "")
    return f"{m.group(1)}_{m.group(2)}" if m else "?"


def _windings(s):
    return sorted(m.group(2) for m in SLOT_RE.finditer(s or "") if m.group(2) in "<>")


def _strip_slots(s):
    return SLOT_RE.sub("", s or "")


def _looks_macrocyclic(s):
    """True if any ligand fragment is a porphyrinoid / large N4 macrocycle.

    String-only heuristic (no RDKit): a big fragment with >=4 ring nitrogens and
    a conjugated imine backbone. Porphyrins/corroles/phthalocyanines re-encode
    with inconsistent macrocycle kekulization, which is an encoder-perception
    issue (S3) even when the net string delta reduces to a few E/Z slashes.
    """
    for frag in _strip_slots(s).split("."):
        n_ring_n = frag.count("n") + frag.count("=N") + frag.count("N=")
        if len(frag) > 55 and (frag.count("N") + frag.count("n")) >= 4 and n_ring_n >= 3:
            return True
    return False


def classify(rep):
    """Return (defect_class, one_line_evidence)."""
    status = rep.get("status")
    err = rep.get("error") or ""
    s1, s2 = rep.get("smiles_1") or "", rep.get("smiles_2") or ""

    if status == "success":
        return "success", ""
    if status not in ("failed", "success"):
        return "interrupted", f"status={status}"
    if not err and not s1:
        return "junk_row", "failed row with no error and no SMILES"

    if "max_attempts" in err and "TMCOptimizer" in err:
        return "stale_harness_TypeError", "harness --quick bug, fixed in v0.3.5"

    if (
        err.startswith(("XYZToSMILES failed", "1D conversion failed", "Input XYZ missing"))
        or "re-encode failed" in err.split("\n")[0]
    ):
        head = err.split("\n")[0]
        if re.search(r"\[H\]B|B-\]|carborane", err):
            return "carborane_unsupported", head[:120]
        if "kekulize" in err.lower():
            return "kekulize_encode_crash", head[:120]
        if "NoneType" in err:
            return "xyz2mol_none_crash", head[:120]
        return "encode_crash_other", head[:120]

    if err.startswith("High RMSD"):
        m = re.search(r"High RMSD at \S+: ([\d.]+)", err)
        val = float(m.group(1)) if m else -1.0
        if val >= 900:
            return "rmsd_sentinel", f"rmsd={val:.0f} (mapping failure, not geometry)"
        return "high_rmsd", f"rmsd={val:.2f}"

    if err.startswith("Atom count mismatch"):
        m = re.search(r"Input (\d+) != Gen (\d+)", err)
        return "donor_H_atom_count", (f"atoms {m.group(1)} -> {m.group(2)}" if m else "")

    if err.startswith("String mismatch"):
        # smoking-gun patterns, most specific first
        base1, base2 = _strip_slots(s1), _strip_slots(s2)
        # exotic main-group cage (borazine, silole): notation limit, not a bug
        if re.search(r"B\d?=N|=\[BH\]|\[si", s1 + s2):
            return "carborane_unsupported", "exotic B/Si cage notation"
        # generated mislocalized a multiple bond onto the ring backbone
        # ([CH2]=[CH2]) or a substituent (=[CH3], [CH2]#) -- bond-order transfer
        if re.search(r"\[CH2?\]=\[CH2?\]|=\[CH3\]|\[CH2\]#", s2) and not re.search(
            r"\[CH2?\]=\[CH2?\]|=\[CH3\]|\[CH2\]#", s1
        ):
            return "eta_diene_localization", "gen mislocalized C=C/C#C (bond-order transfer)"
        # over-valent neutral nitro/nitrite: N bonded to 3 O with no charge, so
        # the fragment is unparseable and the =O position drifts on re-encode
        if re.search(r"N\(O\)=O|N\(=O\)O", s1 + s2):
            return "garbled_aromatic", "over-valent nitro (should be [N+](=O)[O-])"
        if re.search(r"c=|=c[0-9()]", s2) and not re.search(r"c=|=c[0-9()]", s1):
            return "garbled_aromatic", "gen emitted mixed aromatic/double bonds (c=)"
        if re.search(r"=\[[ON]H\d?\]", s2) and not re.search(r"=\[[ON]H\d?\]", s1):
            return "H_on_terminal_oxo_imido", "gen protonated =O / =N donor"
        # porphyrinoid macrocycle: inconsistent kekulization on re-encode. Checked
        # BEFORE atom_stereo/E-Z because the net delta often reduces to a few @/
        # slashes that are symptoms of the localized form, not real stereocenters.
        if _looks_macrocyclic(s1):
            return (
                "macrocycle_perception",
                "porphyrinoid macrocycle re-encoded with inconsistent kekule",
            )
        if base1.replace("@", "") == base2.replace("@", "") and base1 != base2:
            return "atom_stereo", "@ tags differ, skeleton identical"
        stripped1 = base1.replace("/", "").replace("\\", "")
        stripped2 = base2.replace("/", "").replace("\\", "")
        if stripped1 == stripped2 and base1 != base2:
            return "EZ_bond_stereo", "E/Z slashes differ, skeleton identical"
        if _windings(s1) != _windings(s2):
            return "winding_flip", f"winding {_windings(s1)} vs {_windings(s2)}"
        m1 = re.match(r"\[[A-Za-z]+_[A-Z0-9]+\]", s1)
        m2 = re.match(r"\[[A-Za-z]+_[A-Z0-9]+\]", s2)
        if m1 and m2 and m1.group(0) != m2.group(0):
            return (
                "geometry_or_fragment_change",
                f"metal/geo token {m1.group(0)} -> {m2.group(0)} (ligand detached?)",
            )
        if s1.count(".") != s2.count("."):
            return (
                "geometry_or_fragment_change",
                f"fragment count {s1.count('.') + 1} -> {s2.count('.') + 1}",
            )
        return "string_mismatch_other", ""

    if err.startswith("Generation/Verification failed"):
        head = err.split("\n")[0]
        if "IndexError" in head and ("frag_vectors" in err or "list index out of range" in head):
            return "frag_vector_IndexError", "metallogen_adapter convert_parsed_to_msmiles"
        if "failed to generate any conformers" in err:
            return "no_conformers", head[:120]
        if "Geometry code 'NON'" in err:
            return "geometry_NON", "encoder emitted g:NON (no geometry template matched)"
        return "gen_exception_other", head[:120]

    if "TimeoutException" in err:
        return "timeout", ""

    return "string_mismatch_other", err.split("\n")[0][:120]


def _load_overrides():
    """molecule -> hand-triaged class, for cases the stdlib classifier can't route.

    Lives next to this script (tools/triage_overrides.json). Applied to failed
    rows only, so a passing molecule is never forced to a defect class.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "triage_overrides.json")
    try:
        with open(path) as f:
            return {k: v for k, v in json.load(f).items() if not k.startswith("_")}
    except FileNotFoundError:
        return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output-dir", required=True, help="Results dir with individual_reports/")
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    indiv_dir = os.path.join(output_dir, "individual_reports")
    overrides = _load_overrides()

    rows = []
    for fp in sorted(glob.glob(os.path.join(indiv_dir, "*.json"))):
        try:
            with open(fp) as f:
                rep = json.load(f)
        except Exception:
            continue
        cls, evidence = classify(rep)
        if cls != "success" and rep.get("molecule") in overrides:
            cls = overrides[rep["molecule"]]
            evidence = f"{evidence} [manual triage override]".strip()
        rows.append(
            {
                "molecule": rep.get("molecule"),
                "class": cls,
                "session": SESSION_OF.get(cls, "unassigned-triage") if cls != "success" else "",
                "evidence": evidence,
                "metal_geo": _metal_geo(rep),
                "commit_id": rep.get("commit_id", "pre-provenance"),
                "saved_at": rep.get("saved_at")
                or datetime.fromtimestamp(os.path.getmtime(fp)).isoformat(timespec="seconds"),
            }
        )

    n_success = sum(1 for r in rows if r["class"] == "success")
    failures = [r for r in rows if r["class"] != "success"]
    counts = Counter(r["class"] for r in failures)

    by_session = defaultdict(list)
    for r in failures:
        by_session[r["session"]].append(r)

    lines = [
        "# tmCAT/tmPHOTO Round-Trip Case Registry",
        "",
        f"Generated {datetime.now().isoformat(timespec='seconds')} by tools/classify_failures.py",
        f"from {len(rows)} individual reports in `{indiv_dir}`.",
        "",
        f"**{n_success} success / {len(failures)} failure ({100 * n_success / max(len(rows), 1):.1f}% pass).**",
        "",
        "Rows whose `commit_id` is `pre-provenance` predate provenance stamping and may be",
        "stale; re-run before investing debugging effort (`--rerun-failed` or `--only`).",
        "",
        "## Failure classes",
        "",
        "| class | count | owning session |",
        "|---|---|---|",
    ]
    for cls, n in counts.most_common():
        lines.append(f"| {cls} | {n} | {SESSION_OF.get(cls, 'unassigned-triage')} |")

    for session in sorted(by_session):
        rs = by_session[session]
        lines += ["", f"## {session} ({len(rs)} molecules)", ""]
        for r in sorted(rs, key=lambda x: (x["class"], x["molecule"] or "")):
            ev = f" — {r['evidence']}" if r["evidence"] else ""
            lines.append(
                f"- `{r['molecule']}` [{r['metal_geo']}] {r['class']}{ev} "
                f"(code {r['commit_id']}, {r['saved_at']})"
            )

    registry_md = os.path.join(output_dir, "CASE_REGISTRY.md")
    with open(registry_md, "w") as f:
        f.write("\n".join(lines) + "\n")

    registry_json = os.path.join(output_dir, "case_registry.json")
    with open(registry_json, "w") as f:
        json.dump(rows, f, indent=1)

    print(f"{len(rows)} reports: {n_success} success, {len(failures)} failures.")
    for cls, n in counts.most_common():
        print(f"  {n:5d}  {cls:32s} -> {SESSION_OF.get(cls, 'unassigned-triage')}")
    print(f"Wrote {registry_md}")
    print(f"Wrote {registry_json}")


if __name__ == "__main__":
    main()
