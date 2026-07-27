"""Coordinate-only "are the ligands still attached?" test for pool acceptance.

WHY THIS EXISTS
===============
``OIN_ACCEPT_SCORED`` drops the one test in ``_reencode_key_matches`` that does not reuse the
generator's own bond graph. Measured (``docs/ACCEPT_SCORED_v0.4.7.md``), that is a large
speedup and it is byte-safe -- ``sha256(smiles_2)`` is identical on 118/118 molecules -- but
it degrades the structure underneath the string: on a 100-molecule population, +28 vdW
clashes, severe clashes 5->14, and independent re-perception lost on 26 molecules with zero
recoveries. Diffing the failures showed **6 of 8 lose metal-ligand coordination outright**
(the eta ring becomes a free molecule and the metal geometry tag degrades by exactly one
donor: ``[Ru_TET]``->``[Ru_TPL]``, ``[Zr_TET]``->``[Zr_LIN]``, ``[Ti_OCT]``->``[Ti_SPY]``).

This module is that lever's missing safety condition: *accept the first conformer the score
credits **that still has its ligands attached***.

THE TRAP IT AVOIDS -- READ THIS BEFORE CHANGING ANYTHING HERE
=============================================================
``_coordination_vectors`` derives donors from ``metal.GetBonds()``. **A ligand that has
physically left the coordination sphere KEEPS ITS BOND OBJECT.** Any attachment check built
on that path is blind by construction and would certify exactly the structures it exists to
catch. So in this module the generator's graph is allowed to supply only the *reference* --
which atoms to look at, and how they group into coordination sites -- and never the
*measurement*. The measurement is :func:`encoder_donor_set`, which reads coordinates and
element numbers only. That is why :func:`ligands_attached` takes ``claimed_sites`` as plain
index tuples rather than a mol: a caller cannot accidentally pass connectivity in.

THE PREDICATE, AND WHY IT IS THIS ONE
=====================================
Falsified against the ``indep`` oracle on 21 molecules / 40 accepted conformers in both arms
(``tools/attach_probe.py``, ``docs/ATTACH_CHECK_v0.4.7.md`` §1). Four candidates were scored
side by side; the two proposed in §6.4 of the promote lane's report both fail:

===========================================  ===========  ==========================
predicate                                    separates/8  false positives/22
===========================================  ===========  ==========================
count:  ``|actual| == |claimed|``                      5  **11**
subset: ``claimed subseteq actual``                    5  **3**
sites == OIN slot count (raw donor row)                7  **9**
**site coverage on the FILTERED donor set**        **7**  **0**
===========================================  ===========  ==========================

The count- and subset-based predicates fail for one reason: **the number and identity of
metal-bonded atoms in the generator's graph is not conserved against the encoder's perception
even on good structures.** ``MEDZUR_comp_0``'s contract mol claims 10 metal bonds while the
coordinate path perceives 7, and it round-trips (``[Ru_TET]`` both sides) -- its Cp ring is
already slipped eta5->eta3 and the round-trip key forgives that. Anything referenced to the
raw donor-atom count therefore rejects good conformers.

What IS conserved is the set of coordination **sites**. So the shipped predicate is:

    every coordination site the generator claims must still retain at least one atom
    inside bonding distance of the metal

Sites come from the same 1.6 A transitive grouping the encoder's hapticity reduction uses
(``_HAPTIC_GROUP_CUTOFF``), so a whole Cp ring is one site. Ring slip leaves the site
populated; a ring that has drifted out of the sphere empties it.

WHY THE ENCODER'S RING-CARBON FILTER IS LOAD-BEARING, NOT A DETAIL
==================================================================
``_get_tmc_mol_impl`` drops a coordinating ring carbon when a neighbouring N/O/P/S is also
coordinating -- the heteroatom is the real donor. Omitting that filter costs BOTH ways:
without it the predicate wrongly rejects 9 of 22 round-tripping conformers, and it MISSES
``DAKGON_comp_0``. DAKGON is the clean demonstration that the filter is the mechanism rather
than a correction: arm B's filtered donor set has 6 atoms, the same size as the claim, but
not the same atoms -- the filter drops the NHC carbene carbons precisely because the adjacent
N has come within bonding distance, which *is* the C->N donor reassignment. Two claimed sites
lose every member and the conformer is rejected.

WHAT IT DOES NOT CATCH -- SHIP THIS HONESTLY
============================================
``POVPIA_comp_0``: a hydrogen detaches and the amine reads as an imine. The metal donor set is
intact and every site is populated, so this predicate accepts it. That is 1 of 8 known
failures left standing, and it is the case §6.4 already called unreachable by any
metal-centred check. **7 of 8 is the ceiling and it is reached; it is not 8 of 8.**

Scope limit, separate from the residual above: this is an ACCEPTANCE condition, not a RETURN
condition. On the CHEAP_ONLY class (``docs/eta_accept_gap_cohort.md``) nothing is ever
accepted, the pool fills to completion and ``_select_by_geometry`` returns a best-by-geometry
conformer that never faced any acceptance test. Those molecules are unaffected here, with the
lever on or off.

COST
====
7-81 ms per conformer (median 23.6 ms over 40 evaluations, 47-109 atoms), against the 48-57 s
that the dropped strict ``_reencode_oin`` costs on an eta conformer -- roughly 1000x.

A 100x cheaper shortcut was tried and REFUTED: taking the metal's row by the distance
criterion alone, skipping ``xyz2AC_obabel``'s valence-cap pruning loop, runs in 0.17 ms and
disagrees with the real row on 11 of 40 conformers. The pruning loop can drop a metal-donor
bond (the code's own ``DUDREA_comp_0`` bridging-hydride example), so the full call is required
for the check to be measuring what the encoder measures.
"""

from __future__ import annotations

import logging

import numpy as np
from rdkit import Chem

from . import _telemetry

logger = logging.getLogger(__name__)

#: A -- same threshold as ``metallogen_adapter._HAPTIC_GROUP_CUTOFF`` and
#: ``oin_aligner._reduce_hapticity``. Kept in step with them: it defines what counts as ONE
#: coordination site, and the encoder and this check must agree on that or the predicate is
#: testing a different question from the one that fails.
HAPTIC_GROUP_CUTOFF = 1.6

#: Tolerance for ``xyz2AC_obabel``. 0.5 (not obabel's 0.45) is what ``get_basic_mol`` passes,
#: with the comment "Modified tolerance to capture haptic bonds". Using anything else here
#: would make the check disagree with the perception whose change it exists to detect.
AC_TOLERANCE = 0.5

_HETEROATOM_DONORS = frozenset({7, 8, 15, 16})  # N, O, P, S -- _get_tmc_mol_impl's own set


def encoder_donor_set(atomic_nums, coords, tolerance: float = AC_TOLERANCE):
    """The metal's donor set exactly as ``XYZToSMILES().convert()`` derives it.

    Mirrors ``get_basic_mol`` -> ``_get_tmc_mol_impl``: the nonzero metal row of
    ``xyz2AC_obabel``'s adjacency matrix, then the aromatic-ring-carbon filter.

    **Coordinates and element numbers only.** No bond object is read, which is the whole
    point (see the module docstring). Returns ``(metal_idx, donor_index_set)``, or
    ``(None, set())`` when there is no transition metal.
    """
    from ..utils.xyz2mol import TRANSITION_METALS_NUM
    from ..utils.xyz2mol_local import xyz2AC_obabel

    znums = [int(z) for z in atomic_nums]
    metal_idx = next((i for i, z in enumerate(znums) if z in TRANSITION_METALS_NUM), None)
    if metal_idx is None:
        return None, set()

    AC, proto = xyz2AC_obabel(znums, coords, tolerance=tolerance)
    raw = set(int(j) for j in np.nonzero(AC[metal_idx])[0])
    if not raw:
        return metal_idx, raw

    # Rebuild the same mol get_basic_mol builds, purely so RDKit can answer IsInRing().
    rw = Chem.RWMol(proto)
    n = len(AC)
    for i in range(n):
        for j in range(i + 1, n):
            if AC[i, j]:
                rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.GetSymmSSSR(mol)
    except Exception:
        # Ring perception failed: fall back to the unfiltered row. That is the PERMISSIVE
        # direction (more donors -> sites stay populated -> the conformer is accepted), so a
        # failure here can never turn a good conformer into a rejection.
        return metal_idx, raw

    keep = set()
    for idx in raw:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            if any(
                nb.GetAtomicNum() in _HETEROATOM_DONORS and nb.GetIdx() in raw
                for nb in atom.GetNeighbors()
            ):
                continue  # the neighbouring heteroatom is the real donor
        keep.add(idx)
    return metal_idx, keep


def group_sites(indices, coords, cutoff: float = HAPTIC_GROUP_CUTOFF):
    """Transitive single-linkage grouping of ``indices`` at ``cutoff``.

    One group == one coordination site, which is the unit the OIN's slot numbering uses: an
    eta ring writes the same slot number on every ring atom. Same construction as
    ``metallogen_adapter._reduce_haptic_positions``.
    """
    idx = list(indices)
    seen: set[int] = set()
    out: list[list[int]] = []
    for a in range(len(idx)):
        if a in seen:
            continue
        stack, comp = [a], []
        while stack:
            c = stack.pop()
            if c in seen:
                continue
            seen.add(c)
            comp.append(idx[c])
            for k in range(len(idx)):
                if k in seen:
                    continue
                if float(np.linalg.norm(coords[idx[c]] - coords[idx[k]])) < cutoff:
                    stack.append(k)
        out.append(sorted(comp))
    return out


def ligands_attached(claimed_donors, atomic_nums, coords, tolerance: float = AC_TOLERANCE):
    """Does every coordination site the generator claims still have an atom bonded to the metal?

    ``claimed_donors`` is the REFERENCE ONLY -- the atom indices the generator says are metal
    donors. It selects which atoms to look at and (via ``group_sites``) how they partition
    into sites. It is never evidence that anything is attached; a detached ligand keeps its
    bond object and therefore stays in this list, which is exactly how its site is caught
    empty.

    Returns ``(ok, detail)``. ``detail`` names the lost sites for telemetry/debugging.
    """
    claimed = sorted(set(int(i) for i in claimed_donors))
    if not claimed:
        # Nothing claimed -> nothing to verify. Abstain rather than reject: a conformer with
        # no metal bonds at all is a different defect and not this predicate's to judge.
        return True, {"reason": "no claimed donors", "sites_lost": 0}

    _metal_idx, actual = encoder_donor_set(atomic_nums, coords, tolerance)
    if _metal_idx is None:
        return True, {"reason": "no transition metal", "sites_lost": 0}

    sites = group_sites(claimed, coords)
    lost = [s for s in sites if not (set(s) & actual)]
    return not lost, {
        "n_sites_claimed": len(sites),
        "n_donors_claimed": len(claimed),
        "n_donors_actual": len(actual),
        "sites_lost": len(lost),
        "lost_atom_indices": [i for s in lost for i in s],
    }


def conformer_ligands_attached(mol, claimed_donors=None):
    """:func:`ligands_attached` for an RDKit mol that carries a conformer.

    ``claimed_donors`` defaults to the metal's bonded neighbours in ``mol`` -- the ONLY place
    this module touches connectivity, and only as the reference set. Returns ``True`` when the
    check cannot be evaluated, so a perception failure never rejects a conformer that the
    predicate has not actually judged.

    ⚠ **That abstain-on-error branch is a loaded gun and it has already gone off once.** Passed
    the wrong object -- MetalloGen's ``Molecule`` rather than an ``rdkit.Chem.Mol`` -- every
    call raised ``AttributeError``, every call abstained, and the check silently became a no-op
    that a full A/B run reported as "no effect" rather than as "not wired up". This is the same
    shape of defect as ``clash.mol_clash_count`` returning 0 on ``AttributeError``, which the
    promote lane had to fix for the same reason. So the branch now RECORDS
    ``adapter.attach_check_unevaluable`` before abstaining: a check that cannot run must be
    loud, because "never rejects" and "never runs" are otherwise indistinguishable in the
    output.
    """
    from ..utils.xyz2mol import TRANSITION_METALS_NUM

    try:
        conf = mol.GetConformer()
        znums = [a.GetAtomicNum() for a in mol.GetAtoms()]
        coords = np.array(
            [
                [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                for i in range(mol.GetNumAtoms())
            ]
        )
        if claimed_donors is None:
            metal_idx = next(
                (a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM),
                None,
            )
            if metal_idx is None:
                return True
            claimed_donors = [
                b.GetOtherAtomIdx(metal_idx) for b in mol.GetAtomWithIdx(metal_idx).GetBonds()
            ]
        ok, detail = ligands_attached(claimed_donors, znums, coords)
        _telemetry.record(
            "adapter.attach_check_rejected" if not ok else "adapter.attach_check_passed", **detail
        )
        return ok
    except Exception as exc:
        _telemetry.record("adapter.attach_check_unevaluable", error=f"{type(exc).__name__}: {exc}")
        logger.debug("attachment check could not be evaluated for a conformer", exc_info=True)
        return True
