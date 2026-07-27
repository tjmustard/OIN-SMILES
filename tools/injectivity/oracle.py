"""Independent isomer-distinctness oracle -- does NOT pass through the OIN encoder.

The injectivity probes need a ground truth for "is this twin a genuinely different
isomer?" that is independent of the thing under test (the OIN string). Two signals:

1. **Geometric enantiomer test** (primary, robust for rigid species). A molecule is
   chiral iff its mirror image cannot be superimposed on it by a *proper* rotation +
   translation, modulo a relabeling of atoms that respects the connectivity graph.
   We mirror the coordinates, enumerate the graph automorphisms, and take the minimum
   proper-rotation RMSD over them. Small ⇒ the mirror is the same structure (achiral);
   large ⇒ the mirror is a distinct enantiomer. This uses only geometry + topology, so
   it certifies distinctness for metal Δ/Λ and biaryl atropisomers alike -- exactly the
   axes where descriptor-based perception (RDKit / InChI) is itself blind.

   Caveat: rigid-superposition conflates conformation with configuration, so it is only
   valid for *rigid* molecules (metal complexes with rigid chelates, biaryls). For a
   flexible achiral molecule the mirror is a non-superimposable conformer and would read
   as "chiral". The Y1 probe fixtures are curated rigid species; the dataset-scale UU
   hunt (Wave 3) needs a configurational comparison instead.

2. **RDKit stereo fingerprint** (secondary, informative). ``FindPotentialStereo`` + CIP
   on the metal-free ligand fragment. Often empty for the hard axes -- which is a
   finding, not a bug: standard cheminformatics does not perceive them either. We do
   NOT toggle the global ``SetUseLegacyStereoPerception`` here; it is process-global and
   would contaminate the encoder's own perception if run in the same process.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from oinsmiles.core.constants import TRANSITION_METALS_NUM
from oinsmiles.utils.perception_tmc import get_tmc_mol

#: min proper-rotation mirror RMSD (Å) above which a rigid structure is called chiral.
#: Validated separation on the Y1 fixtures: achiral controls ~0.05, chiral 3-4.
CHIRALITY_RMSD_THRESHOLD = 0.5

#: cap on the number of graph automorphisms enumerated. Generous: more automorphisms
#: only give MORE chances to find an aligning permutation, so a too-low cap risks
#: calling an achiral-but-symmetric molecule chiral. Logged when hit.
MAX_AUTOMORPHISMS = 4000


@dataclass
class OracleVerdict:
    """Independent verdict on whether ``mirror(structure)`` is a distinct isomer."""

    chiral: bool
    rmsd: float
    n_automorphisms: int
    automorphism_cap_hit: bool
    method: str = "geometric-enantiomer"
    fingerprint: dict = field(default_factory=dict)
    note: str = ""

    @property
    def distinct(self) -> bool:
        """The mirror is a genuinely different isomer (ground truth for the probe)."""
        return self.chiral


def _kabsch_proper_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Min RMSD aligning ``Q`` onto ``P`` by translation + a PROPER rotation (det=+1).

    Standard Kabsch would happily return a reflection when that fits better; that is
    exactly what we must forbid, since a reflection is what distinguishes an enantiomer.
    """
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, _, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Qr = Qc @ R.T
    return float(np.sqrt(((Qr - Pc) ** 2).sum(1).mean()))


def geometric_chirality(
    mol: Chem.Mol,
    coords: np.ndarray,
    *,
    max_autos: int = MAX_AUTOMORPHISMS,
    threshold: float = CHIRALITY_RMSD_THRESHOLD,
) -> tuple[float, int, bool]:
    """Return (min mirror RMSD over proper rotations + automorphisms, n_autos, is_chiral).

    ``coords`` are the atom coordinates in ``mol``'s atom order. The mirror is a z-axis
    reflection (any single-axis reflection is an enantiomer-generating improper isometry).
    """
    coords = np.asarray(coords, dtype=float)
    mirror = coords.copy()
    mirror[:, 2] *= -1.0
    n = mol.GetNumAtoms()
    matches = mol.GetSubstructMatches(mol, uniquify=False, useChirality=False, maxMatches=max_autos)
    best = np.inf
    for m in matches:
        if len(m) != n:
            continue
        perm = np.fromiter(m, dtype=int, count=n)
        best = min(best, _kabsch_proper_rmsd(coords, mirror[perm]))
    if best is np.inf:  # no full-graph automorphism (should not happen: identity matches)
        best = _kabsch_proper_rmsd(coords, mirror)
    return float(best), len(matches), bool(best > threshold)


def metal_free_fragment(mol: Chem.Mol) -> Chem.Mol:
    """A sanitized copy of ``mol`` with every transition-metal atom removed."""
    rw = Chem.RWMol(mol)
    for idx in sorted(
        (a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM),
        reverse=True,
    ):
        rw.RemoveAtom(idx)
    frag = rw.GetMol()
    try:
        Chem.SanitizeMol(frag)
    except Exception:
        try:
            Chem.SanitizeMol(
                frag,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        except Exception:
            pass
    return frag


def stereo_fingerprint(mol: Chem.Mol) -> dict:
    """RDKit potential-stereo + CIP fingerprint on the metal-free ligand fragment.

    Informative secondary signal only. Does NOT mutate global RDKit perception state.
    """
    from rdkit.Chem import FindPotentialStereo, rdCIPLabeler

    frag = metal_free_fragment(mol)
    try:
        si = FindPotentialStereo(frag)
        pot = dict(Counter(str(s.type).split(".")[-1] for s in si))
    except Exception as e:
        pot = {"error": str(e)}
    cips: list[str] = []
    try:
        probe = Chem.Mol(frag)
        Chem.AssignStereochemistryFrom3D(probe)
        rdCIPLabeler.AssignCIPLabels(probe)
        cips = [a.GetPropsAsDict()["_CIPCode"] for a in probe.GetAtoms() if a.HasProp("_CIPCode")]
    except Exception:
        pass
    return {"potential_stereo": pot, "cip_codes": sorted(cips)}


def load_mol(xyz_path: str | Path, charge: int = 0) -> tuple[Chem.Mol, np.ndarray]:
    """Load an XYZ into a bonded RDKit mol + coordinate array (the encoder's own view)."""
    mol, _ = get_tmc_mol(Path(xyz_path), charge, with_stereo=False)
    coords = mol.GetConformer().GetPositions()
    return mol, coords


def is_distinct_enantiomer(
    xyz_path: str | Path,
    *,
    charge: int = 0,
    with_fingerprint: bool = True,
) -> OracleVerdict:
    """High-level oracle: is ``mirror(structure at xyz_path)`` a distinct isomer?

    Independent of the OIN encoder. Uses the geometric enantiomer test; attaches the
    RDKit fingerprint for corroboration.
    """
    mol, coords = load_mol(xyz_path, charge)
    rmsd, n_autos, chiral = geometric_chirality(mol, coords)
    fp = stereo_fingerprint(mol) if with_fingerprint else {}
    cap_hit = n_autos >= MAX_AUTOMORPHISMS
    note = "automorphism cap hit -- chirality verdict may be unreliable" if cap_hit else ""
    return OracleVerdict(
        chiral=chiral,
        rmsd=rmsd,
        n_automorphisms=n_autos,
        automorphism_cap_hit=cap_hit,
        fingerprint=fp,
        note=note,
    )
