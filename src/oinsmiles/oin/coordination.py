"""Coordination-integrity check — does the GENERATED geometry still hold its ligands on?

WHY THIS EXISTS
===============
The round-trip success metric cannot detect a detached ligand. ``tools/test_dataset_roundtrip.py``
scores with ``get_oin_string(gen_result.mol, coords)``, and ``gen_result.mol`` carries **the
generator's own bond graph**. A cyclopentadienyl ring that has drifted off the metal is still
bonded *in that graph*, so the re-encode reproduces the input's coordination and the key matches.

Measured rate on the default path (``docs/METRIC_FALSE_POSITIVES.md``, 633 scored successes of
``results-v0.4.5-rebaseline``): **61/633 = 9.6 % overall, 48/171 = 28.1 % on haptic molecules.**
``FIYHUT_comp_0`` is the clean example — textbook ferrocene Fe–C is ≈2.05 Å and the input matches
to 0.03 Å, while the generated structure puts all ten ring carbons at 2.84–2.96 Å (~0.85 Å too
far, 10 bonded carbons → 0). Both rings are off the iron and it scores as a successful round trip.

WHAT THIS DOES, AND WHY IT IS GEOMETRIC RATHER THAN A RE-ENCODE
===============================================================
The obvious implementation — re-encode the generated XYZ with the full ``XYZToSMILES().convert()``
and diff the strings — costs seconds per molecule and drags in the whole perception stack. This
works on **distances only**: count what is inside covalent contact of the metal in each structure
and compare the multisets. Milliseconds, no perception, and — the point — it consults **neither**
bond graph, so it cannot be fooled the way the metric is.

It deliberately does **not** try to map OIN slot ``{n}`` onto a generated atom index. Slot→atom
correspondence is exactly the fragile part (``xyz2mol.py``'s marker placement has a documented
unreliable ``GetSubstructMatch`` fallback), and it is not needed: an aggregate multiset comparison
is order-independent and already separates the failure from the pass.

CONTRACT
========
``coordination_report`` is a **DIAGNOSTIC**. It returns findings; it does not raise and it does not
gate. Shipping it as a gate immediately would move ~61 molecules from pass to fail in one step,
which is indistinguishable from a regression — the confound this project has been caught by twice
(v0.4.4's 11 "regressions" were all timeouts). Record first, decide later.
"""

from __future__ import annotations

import numpy as np
from rdkit.Chem import GetPeriodicTable

from ..core.constants import TRANSITION_METALS_NUM

_PT = GetPeriodicTable()

#: Bonded-contact slack over the sum of covalent radii, in Angstrom. Matches the criterion
#: ``xyz2AC_obabel`` uses for adjacency (``xyz2mol_local.py``), so "in contact" here means the
#: same thing it means to the encoder. Near this boundary perception is genuinely ambiguous
#: (see ``docs/METRIC_FALSE_POSITIVES.md`` §3.2) — which is why the report also carries the raw
#: distances, so a marginal call can be inspected rather than silently decided.
CONTACT_SLACK = 0.45

#: A contact must move by more than this (Angstrom) before a count change is worth reporting, so
#: thermal-scale jitter at the boundary does not manufacture findings.
MARGINAL_BAND = 0.10


def _z(symbol: str) -> int:
    try:
        return _PT.GetAtomicNumber(symbol)
    except Exception:
        return 0


def parse_xyz(text: str):
    """``(symbols, coords)`` from XYZ content. ``([], empty)`` on anything malformed."""
    lines = text.splitlines()
    if not lines:
        return [], np.zeros((0, 3))
    try:
        n = int(lines[0].strip())
    except ValueError:
        return [], np.zeros((0, 3))
    syms, xyz = [], []
    for i in range(n):
        try:
            parts = lines[2 + i].split()
            syms.append(parts[0])
            xyz.append([float(v) for v in parts[1:4]])
        except (IndexError, ValueError):
            return [], np.zeros((0, 3))
    return syms, np.asarray(xyz, dtype=float)


def metal_indices(symbols) -> list:
    """Indices of transition-metal atoms, in input order.

    The project's invariant is that the metal centre is ``fragments[0]``, but that is a property of
    the *encoded* form; here we are looking at raw coordinates, so the metal is found by element.
    """
    return [i for i, s in enumerate(symbols) if _z(s) in TRANSITION_METALS_NUM]


def metal_contacts(symbols, coords, metal_index: int, slack: float = CONTACT_SLACK):
    """Atoms within covalent contact of ``metal_index``.

    Returns ``(counts, contacts)``: ``counts`` maps element symbol -> how many of that element are
    in contact; ``contacts`` is a list of ``(index, symbol, distance, cutoff)`` sorted by distance,
    retained so a caller can show *why* a verdict was reached rather than only that it was.
    """
    if not symbols or coords.shape[0] != len(symbols):
        return {}, []
    if not (0 <= metal_index < len(symbols)):
        return {}, []
    m_r = _PT.GetRcovalent(symbols[metal_index])
    d = np.linalg.norm(coords - coords[metal_index], axis=1)
    counts: dict = {}
    contacts = []
    for i, s in enumerate(symbols):
        if i == metal_index:
            continue
        cutoff = m_r + _PT.GetRcovalent(s) + slack
        if d[i] < cutoff:
            counts[s] = counts.get(s, 0) + 1
            contacts.append((i, s, float(d[i]), float(cutoff)))
    contacts.sort(key=lambda t: t[2])
    return counts, contacts


def coordination_report(input_xyz_text: str, generated_xyz_text: str, slack: float = CONTACT_SLACK):
    """Compare metal coordination between the input and the generated geometry.

    Returns a JSON-serialisable dict. ``intact`` is True when every metal keeps at least as many
    contacts of each element as it had in the input. ``None`` means *not assessable* (no metal, a
    malformed file, a differing metal count) and must never be read as a pass — the whole point of
    this module is that a silent zero or a silent True is how the existing metric fails.

    ``lost`` is the per-element shortfall. ``marginal`` flags contacts that merely crossed the
    cutoff by less than :data:`MARGINAL_BAND`, so a caller can tell a ring that left the metal from
    a bond length that wandered a few hundredths of an Angstrom over the boundary.
    """
    out: dict = {"intact": None, "reason": None, "metals": []}

    in_syms, in_xyz = parse_xyz(input_xyz_text)
    gen_syms, gen_xyz = parse_xyz(generated_xyz_text)
    if not in_syms or not gen_syms:
        out["reason"] = "unparsable xyz"
        return out

    in_metals, gen_metals = metal_indices(in_syms), metal_indices(gen_syms)
    if not in_metals or not gen_metals:
        out["reason"] = "no transition metal found"
        return out
    if len(in_metals) != len(gen_metals):
        out["reason"] = f"metal count differs: {len(in_metals)} -> {len(gen_metals)}"
        return out

    intact = True
    boundary_only = False
    for k, (mi, mg) in enumerate(zip(in_metals, gen_metals)):
        in_counts, in_contacts = metal_contacts(in_syms, in_xyz, mi, slack)
        gen_counts, gen_contacts = metal_contacts(gen_syms, gen_xyz, mg, slack)

        lost = {}
        for el, n in in_counts.items():
            deficit = n - gen_counts.get(el, 0)
            if deficit > 0:
                lost[el] = deficit

        # Is the loss REAL, or did an atom drift a few hundredths of an Angstrom over the cutoff?
        #
        # This distinction is load-bearing and was added on evidence. Validated against 633
        # molecules with an independent string-level ground truth, a raw loss verdict flagged 45
        # genuine passes -- and 36 of those had a contact within MARGINAL_BAND of the cutoff. The
        # worst were haptic: Ru 9->3, Cr 9->4 on arene rings whose carbons sit right at the
        # boundary, where a per-atom covalent criterion is stricter than the encoder's own ring
        # perception. FIYHUT by contrast has its Fe-C at 2.86 A against a 2.53 A cutoff -- 0.33 A
        # beyond, more than 3x the band -- so the real detachments are untouched by this.
        #
        # A loss made up ENTIRELY of atoms just outside the cutoff is therefore reported as
        # `boundary` rather than as a degradation: inconclusive, not clean, and not a finding.
        beyond = 0
        m_r_gen = _PT.GetRcovalent(gen_syms[mg])
        d_all = np.linalg.norm(gen_xyz - gen_xyz[mg], axis=1)
        for el in lost:
            slack_out = sorted(
                d_all[i] - (m_r_gen + _PT.GetRcovalent(gen_syms[i]) + slack)
                for i in range(len(gen_syms))
                if i != mg
                and gen_syms[i] == el
                and d_all[i] >= (m_r_gen + _PT.GetRcovalent(gen_syms[i]) + slack)
            )
            beyond += sum(1 for v in slack_out[: lost[el]] if v >= MARGINAL_BAND)
        gained = {}
        for el, n in gen_counts.items():
            surplus = n - in_counts.get(el, 0)
            if surplus > 0:
                gained[el] = surplus

        # A contact is "marginal" when the generated distance sits just outside the cutoff -- the
        # regime where the encoder's own perception is ambiguous. Reported separately so a caller
        # never has to guess whether a finding is a detachment or boundary jitter.
        gen_d = np.linalg.norm(gen_xyz - gen_xyz[mg], axis=1)
        gen_cut = np.array(
            [_PT.GetRcovalent(gen_syms[mg]) + _PT.GetRcovalent(s) + slack for s in gen_syms]
        )
        marginal = int(
            sum(
                1
                for i in range(len(gen_syms))
                if i != mg and 0 <= (gen_d[i] - gen_cut[i]) < MARGINAL_BAND
            )
        )

        if lost and beyond > 0:
            intact = False
        elif lost:
            boundary_only = True
        out["metals"].append(
            {
                "metal": in_syms[mi],
                "index_in": mi,
                "index_gen": mg,
                "n_contacts_in": sum(in_counts.values()),
                "n_contacts_gen": sum(gen_counts.values()),
                "counts_in": in_counts,
                "counts_gen": gen_counts,
                "lost": lost,
                "lost_beyond_band": beyond,
                "gained": gained,
                "marginal_gen": marginal,
                "nearest_in": [(s, round(d, 3)) for _i, s, d, _c in in_contacts[:12]],
                "nearest_gen": [(s, round(d, 3)) for _i, s, d, _c in gen_contacts[:12]],
            }
        )
        del k

    out["intact"] = intact
    out["boundary_only"] = bool(boundary_only and intact)
    if not intact:
        out["reason"] = "; ".join(
            f"{m['metal']}: lost {m['lost']} ({m['n_contacts_in']} -> {m['n_contacts_gen']})"
            for m in out["metals"]
            if m["lost"] and m["lost_beyond_band"]
        )
    elif out["boundary_only"]:
        # Contacts were lost, but every one of them sits within MARGINAL_BAND of the cutoff. Not
        # charged as a degradation and not called clean either -- see the note above `beyond`.
        out["reason"] = "; ".join(
            f"{m['metal']}: {sum(m['lost'].values())} contact(s) at the cutoff boundary"
            for m in out["metals"]
            if m["lost"]
        )
    # KNOWN SCOPE LIMIT, measured rather than assumed. This verdict is LOSS-based, so two real
    # failure modes are invisible to it and are left to the caller via `gained` / `counts_*`:
    #   1. over-coordination -- 4 of 61 known false positives GAINED contacts (6->11, 7->12) and
    #      changed geometry tag without losing anything. A gain is not reliably bad, though: a
    #      genuine pass in the same corpus gained 2 (Mo 6->8), so no threshold is asserted here.
    #   2. same-count hapticity rearrangement -- OGARAP goes eta3 -> eta2 with the carbon count
    #      unchanged, which an aggregate per-element multiset cannot see.
    return out


def summarize(report: dict) -> str:
    """One-line human form, for a log or a report field."""
    if report.get("intact") is None:
        return f"coordination: NOT ASSESSED ({report.get('reason')})"
    if report["intact"]:
        n = sum(m["n_contacts_gen"] for m in report["metals"])
        return f"coordination: intact ({n} metal contacts)"
    return f"coordination: DEGRADED -- {report['reason']}"
