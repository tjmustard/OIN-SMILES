"""Curated metal–ligand bond lengths for the default (MetalloGen) 3D generator.

The default engine places each σ donor at a generic Cordero covalent-radius sum
(``(metal_r + ligand_r) * scale``; ``chem.py`` radii, *Dalton Trans., 2008,
2832-2838*). That sum **systematically overestimates** real dative metal–ligand
distances — most sharply for M–P, M–halide, M–S and early-transition-metal M–O/N —
because a covalent radius does not describe a lone-pair donor bond. This module
supplies hand-curated per-``(metal, ligand)`` distances so the FF-clean scan target
(``clean_geometry.ff_clean``) can pin σ donors at physically realistic separations.

This table originated as a **verbatim copy** of ``_BOND_LENGTHS`` from the legacy
Molassembler backend; that backend was removed in v0.3.7, so this module is now the
authoritative source for these distances (a hand-copied constant is exactly how Sc/Y
once went missing — TD-005).

Unlike legacy ``_bond_length`` (which returns a blanket ``2.10`` for any unlisted
pair), :func:`bond_length` returns ``None`` for an absent pair, so the caller falls
back to the **covalent sum** rather than a flat 2.10 Å that could be worse than it.
"""

from __future__ import annotations

# Typical metal–ligand bond lengths (Angstrom). Originally copied from the legacy
# Molassembler backend's ``_BOND_LENGTHS`` (removed in v0.3.7); now authoritative.
BOND_LENGTHS: dict[str, dict[str, float]] = {
    "Ir": {
        "C": 2.00,
        "N": 2.12,
        "O": 2.05,
        "P": 2.30,
        "Cl": 2.35,
        "Br": 2.50,
        "I": 2.65,
        "S": 2.30,
        "F": 1.95,
    },
    "Pt": {
        "Cl": 2.31,
        "N": 2.10,
        "P": 2.25,
        "C": 2.00,
        "O": 2.05,
        "S": 2.30,
        "Br": 2.40,
        "I": 2.60,
        "F": 2.05,
    },
    "Pd": {
        "Cl": 2.30,
        "N": 2.10,
        "P": 2.30,
        "C": 2.00,
        "O": 2.05,
        "S": 2.30,
        "Br": 2.40,
        "As": 2.40,
        "F": 2.10,
    },
    "Fe": {
        "C": 1.80,
        "N": 2.00,
        "O": 2.00,
        "Cl": 2.20,
        "Br": 2.35,
        "I": 2.50,
        "H": 1.55,
        "P": 2.20,
    },
    "Cr": {"C": 1.85, "N": 2.05, "O": 1.95, "Cl": 2.30},
    "Re": {"C": 1.90, "N": 2.15, "O": 2.00, "Cl": 2.40},
    "Au": {"C": 2.00, "N": 2.10, "Cl": 2.30, "P": 2.25, "S": 2.30},
    "Rh": {"C": 2.00, "N": 2.10, "O": 2.05, "Cl": 2.35, "P": 2.25},
    "Cu": {"C": 1.90, "N": 2.00, "O": 1.95, "Cl": 2.25},
    "Hg": {"I": 2.65, "Cl": 2.40, "Br": 2.50, "N": 2.30},
    "Ti": {"Cl": 2.30, "N": 2.05, "O": 1.80, "C": 2.05},
    "V": {"O": 1.60, "N": 2.00, "Cl": 2.30},
    "Zn": {"N": 2.05, "O": 2.00, "Cl": 2.25, "Br": 2.35, "S": 2.25, "Se": 2.40},
    "Ni": {"N": 2.00, "O": 2.00, "Cl": 2.20, "P": 2.20, "C": 1.90, "S": 2.20},
    "Ag": {"N": 2.20, "O": 2.30, "S": 2.40, "Cl": 2.40},
    "Cd": {"N": 2.30, "O": 2.30, "S": 2.55, "Cl": 2.50},
    "Ru": {"C": 1.90, "N": 2.10, "O": 2.00, "Cl": 2.35, "P": 2.25},
    "La": {"N": 2.60, "O": 2.50, "S": 2.80, "Se": 2.90},
}


# Metals for which the curated table is *applied* by the generator. This is a
# dataset-validated subset of BOND_LENGTHS' 18 metals, not the whole table.
#
# Why a subset. Measured against real input geometries (coordination-sphere mean
# RMSD over the tmCAT-tmPHOTO dataset, ~6 molecules/metal × 3 seeds), the table
# improves the late/post-transition d⁸/d¹⁰ metals — Ni, Pd, Pt (d⁸) and Zn, Cd, Hg,
# Ag (d¹⁰) — whose dative M–L bonds a generic covalent-radius sum systematically
# *over*estimates. It *regresses* several early/mid-transition-metal buckets,
# because some of its entries encode a shorter bond mode than the dative bond
# actually present: Ti–O 1.80 and V–O 1.60 are metal-**oxo** (M=O) distances, and
# the early-TM M–C values are shorter than a real σ M–C. The generator cannot tell
# an oxo from an alkoxide, or a carbene from an alkyl, at this seam. A few metals
# with a strong median but inconsistent per-molecule sign (Rh, Ir) or a flat median
# (Cu, Re, Fe) are also left off — validation could not show they *strictly* help.
#
# Rather than ship a net win that regresses some buckets (the acceptance forbids any
# per-metal regression), the table is gated to the seven metals where validation
# shows it strictly improves fidelity (per-metal median RMSD −0.015 to −0.084 Å;
# every molecule bucket improving). Every other metal keeps the covalent-radius sum
# and is therefore byte-identical to the pre-table generator — the change is
# strictly additive. Expand this set only with the same per-metal RMSD validation,
# never by hand (TD-005).
ENABLED_METALS: frozenset[str] = frozenset({"Ni", "Pd", "Pt", "Zn", "Cd", "Hg", "Ag"})


# σ-donor ``(metal, ligand)`` pairs excluded from the curated table even though the
# metal is in :data:`ENABLED_METALS`. Their curated target sits far enough *below*
# the covalent sum that pinning the donor there (``clean_geometry.ff_clean``'s
# ``maxDispl=0, forceConstant=2000`` constraint) over-strains the UFF clean: the
# ligand-adjacency and "other atoms bind" guards then restore-and-retry far more
# often. Measured on PdCl₂-BINAP (P8, seed 42, deterministic): Pd–P at the 2.30 Å
# table target (vs the 2.46 Å covalent sum) quadrupled ``ff_clean`` failures 3 → 12
# and pushed generation 34 s → 57 s; routing Pd–P alone back to the covalent sum
# restores both (fail=3, 36 s) at an RMSD give-back below P4's ≥0.05 Å acceptance
# floor. Pd–Cl is deliberately *not* listed: with Pd–P exempt it triggers no extra
# retries (fail stays 3), so it keeps its P4 fidelity gain. A σ donor whose pair is
# listed falls back to the covalent sum — byte-identical to the pre-table generator
# for that donor — while every unlisted pair keeps its :data:`ENABLED_METALS`
# behaviour. Expand this set only with the same speed + per-metal RMSD validation
# P8 used, never by hand.
SHORT_PIN_EXEMPT_PAIRS: frozenset[tuple[str, str]] = frozenset({("Pd", "P")})


def sigma_table_applies(metal: str, ligand: str) -> bool:
    """Whether a σ donor of ``(metal, ligand)`` takes the curated bond-length target.

    ``True`` only when ``metal`` is in :data:`ENABLED_METALS` **and** the pair is not
    in :data:`SHORT_PIN_EXEMPT_PAIRS`. A ``False`` result means the caller
    (``clean_geometry._binding_distance``) falls back to the covalent-radius sum,
    exactly as it does for a non-enabled metal — so an exempt pair is byte-identical
    to the pre-table generator for that donor.
    """
    return metal in ENABLED_METALS and (metal, ligand) not in SHORT_PIN_EXEMPT_PAIRS


def bond_length(metal: str, ligand: str) -> float | None:
    """Curated metal–ligand distance (Å), or ``None`` when the pair is absent.

    This is the raw table accessor — it does **not** apply the
    :data:`ENABLED_METALS` gate; the caller (``clean_geometry._binding_distance``)
    does. ``None`` (not a blanket default) is deliberate: an unlisted pair must fall
    back to the covalent-radius sum, which is a better estimate than a flat 2.10 Å.
    """
    return BOND_LENGTHS.get(metal, {}).get(ligand)
