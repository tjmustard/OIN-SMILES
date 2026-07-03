"""Shared constants for the OIN-SMILES core and utils packages.

Single source of truth for the transition-metal atomic-number predicate used
by both the XYZ->OIN encode pipeline (``utils/xyz2mol.py``) and the
chirality module (``core/chirality.py``).

This module exists specifically to avoid a circular import: ``xyz2mol.py``
imports ``ChiralityRecoveryUtility`` from ``core/chirality.py``, so
``core/chirality.py`` cannot import ``TRANSITION_METALS_NUM`` back out of
``utils/xyz2mol.py`` directly. Both modules import from here instead (TD-005:
never duplicate the metal list -- that multiplies the stale-``is_metal``
failure surface).
"""

from __future__ import annotations

# fmt: off
TRANSITION_METALS: list[str] = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "La", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
]

TRANSITION_METALS_NUM: list[int] = [
    21, 22, 23, 24, 25, 26, 27, 57, 28, 29, 30, 39, 40, 41,
    42, 43, 44, 45, 46, 47, 48, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
]
# fmt: on
