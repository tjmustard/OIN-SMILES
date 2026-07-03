"""Shared winding-sign computation for haptic (eta) groups.

Single source of truth for the OIN V3.6 winding-direction sign, used by both
the XYZ->OIN encoder (``utils/oin_aligner.py::_determine_winding``) and the
OIN->XYZ generation-side haptic-face correction (Stereo Phase 3). Extracting
the sign math into one helper prevents the two call sites from silently
diverging on the winding convention.
"""

import numpy as np


def signed_circulation(coords, star_local_idx, axis) -> str:
    """Return the OIN winding character (``'>'`` or ``'<'``) for a haptic group.

    Args:
        coords: Array-like of shape (n, 3) with the group's atom positions,
            in **SMILES/fragment order** (i.e. ``coords[i]`` corresponds to
            the atom at SMILES-order position ``i`` within the group). Coords
            need not be pre-centered -- this function centers them on their
            own centroid internally.
        star_local_idx: Index into ``coords`` (SMILES-order position, not an
            atom identifier) of the "star" atom -- the heading/reference atom
            whose winding relative to its cyclic-next neighbour is measured.
        axis: The reference axis vector, by convention **metal -> centroid,
            outward**. Sign is measured relative to this direction, so callers
            must ensure `axis` actually points outward (see
            `utils/oin_aligner.py::_determine_winding` for the encoder-side
            assertion/normalization of this convention).

    Returns:
        ``'>'`` (clockwise / forward) if ``n < 3`` (degenerate default -- e.g.
        linear/monodentate groups), or if
        ``cross(v_star, v_next) . axis >= 0``. Otherwise ``'<'``
        (counter-clockwise / reverse).
    """
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if n < 3:
        return ">"  # Linear/monodentate groups default to Forward.

    centered = coords - coords.mean(axis=0)

    v_star = centered[star_local_idx]
    next_local_idx = (star_local_idx + 1) % n
    v_next = centered[next_local_idx]

    winding_normal = np.cross(v_star, v_next)
    dot = np.dot(winding_normal, axis)

    return ">" if dot >= 0 else "<"
