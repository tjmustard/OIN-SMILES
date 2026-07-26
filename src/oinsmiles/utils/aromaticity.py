"""Shared handling for aromatic systems that have no valid Kekule structure.

A *stuck ring* is an aromatic ring whose atoms carry an exocyclic double bond --
a quinoid 2-iminopyridine donor, or the ``P=c`` ylide bond that ``AC2mol`` draws
onto a phenyl when it perceives the wrong ligand charge. RDKit cannot kekulize
such a ring, so ``Chem.Kekulize`` and any full ``Chem.SanitizeMol`` raise.

Both directions of the round trip hit this:

* OIN -> XYZ: ``generator3d.process.get_ace_mol_from_rd_mol`` needs Kekule bond
  orders to build its bond-order matrix.
* XYZ -> OIN: ``utils.xyz2mol.get_tmc_mol`` sanitizes the assembled TMC.

The blanket fallback -- treat every aromatic bond as single -- also degrades the
ligand's other, well-behaved rings. Clearing aromaticity on *only* the stuck rings
preserves correct alternating bond orders everywhere else.
"""

from rdkit import Chem

from ..oin.levers import lever_enabled as _lever_enabled


class OINEncodeError(ValueError):
    """The encoder cannot build a valid molecule from this input.

    Subclasses ``ValueError`` so existing ``except ValueError`` handlers (and
    ``core.translator``, which re-raises as ``ValueError``) keep working.
    """


def stuck_ring_atoms(mol):
    """Atom indices of every aromatic ring carrying an exocyclic double bond.

    Requires ring info; the caller should have sanitized at least far enough to
    perceive rings (``SANITIZE_ALL ^ SANITIZE_KEKULIZE`` is enough). Returns an
    empty set when ring info is unavailable rather than raising.
    """
    stuck = set()
    try:
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
    except Exception:
        return stuck

    for ring in rings:
        ring_set = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                continue
            for bond in atom.GetBonds():
                if (
                    bond.GetOtherAtomIdx(idx) not in ring_set
                    and bond.GetBondType() == Chem.BondType.DOUBLE
                ):
                    stuck.update(ring)
    return stuck


def clear_ring_aromaticity(rw_mol, atom_indices):
    """Drop aromatic flags on the given atoms and on the bonds between them.

    An ``AROMATIC``-typed bond inside the cleared set is demoted to ``SINGLE`` so
    nothing downstream trips over a bond type with no Kekule meaning. Mutates
    ``rw_mol`` in place.
    """
    for idx in atom_indices:
        rw_mol.GetAtomWithIdx(idx).SetIsAromatic(False)
    for bond in rw_mol.GetBonds():
        if bond.GetBeginAtomIdx() in atom_indices and bond.GetEndAtomIdx() in atom_indices:
            bond.SetIsAromatic(False)
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetBondType(Chem.BondType.SINGLE)


def dearomatize_stuck_rings(rd_molecule, add_hydrogen):
    """Kekulize a molecule whose aromatic system has no global Kekule structure.

    Clears aromaticity on only the stuck ring(s) and kekulizes the rest normally,
    preserving correct alternating bond orders on the good rings (a blanket
    "every aromatic bond is single" fallback makes MetalloGen embed an innocent
    2,6-dimethylphenyl bridge non-planar, and it then re-encodes as quinoid).
    Returns a best-effort mol -- the original on any failure; the bond-order
    lookup in ``process`` still guards against a stray AROMATIC bond.
    """
    try:
        mol = Chem.AddHs(rd_molecule) if add_hydrogen else Chem.RWMol(rd_molecule).GetMol()
        mol = Chem.RWMol(mol)
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception:
            pass
        stuck = stuck_ring_atoms(mol)
        if not stuck:
            return rd_molecule
        clear_ring_aromaticity(mol, stuck)
        out = mol.GetMol()
        out.UpdatePropertyCache(strict=False)
        Chem.Kekulize(out, clearAromaticFlags=True)
        return out
    except Exception:
        return rd_molecule


def _boron_cage_relaxation_applies(mol):
    """Whether the ``OIN_BORON_CAGE`` valence relaxation should apply to ``mol``.

    Two conditions, both required: the lever is set, and the mol actually contains
    a B-B-B triangle (the deltahedral 3c-2e cage motif). Imported lazily because
    ``xyz2mol_local`` pulls in the heavy perception stack and this module is
    imported by it in turn.
    """
    if not _lever_enabled("OIN_BORON_CAGE"):
        return False
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    if sum(1 for z in atoms if z == 5) < 3:
        return False
    from .xyz2mol_local import boron_cage_vertices

    return bool(boron_cage_vertices(atoms, Chem.rdmolops.GetAdjacencyMatrix(mol)))


def sanitize_allowing_boron_cage(mol):
    """``Chem.SanitizeMol(mol)``, but tolerating a deltahedral boron cage.

    A drop-in replacement for a bare ``Chem.SanitizeMol`` call on the encode path.
    With ``OIN_BORON_CAGE`` unset, or on any mol without a B-B-B triangle, this is
    **exactly** ``Chem.SanitizeMol(mol)`` -- same call, same exceptions, so the
    ~6,600 non-cage molecules are byte-identical. Only a cage mol takes the
    ``SANITIZE_ALL ^ SANITIZE_PROPERTIES`` path, and if that also fails the
    original strict error is raised so the failure is never silently swallowed.
    """
    if not _boron_cage_relaxation_applies(mol):
        Chem.SanitizeMol(mol)
        return mol
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    except Exception:
        Chem.SanitizeMol(mol)
    return mol


def kekulize_safe_sanitize(mol):
    """``Chem.SanitizeMol``, retried with the stuck rings de-aromatized.

    ``AC2mol`` can perceive a ligand charge that forces a double bond from a
    phosphorus onto an aromatic ring carbon (a PPN+ counter-cation, a phosphonium
    ylide). The ring is then unkekulizable and the first full sanitize raises
    ``Can't kekulize mol. Unkekulized atoms: ...`` -- 17 molecules died this way,
    every traceback pointing at whichever sanitize happened to run first.

    Retry once with aromaticity cleared on just the offending rings. If that still
    fails, raise ``OINEncodeError`` naming the atoms, so the caller reports a
    specific limitation instead of a bare kekulize traceback. Relaxing the ring
    cannot rescue a mol that is invalid for another reason -- at the wrong ligand
    charge ``AC2mol`` leaves the ipso carbon pentavalent, and only re-perceiving the
    charge fixes that -- so the retry error is reported verbatim.

    Returns a sanitized mol of the same class as the input (callers pass ``RWMol``
    and go on to call ``GetMol()``). Atom properties, including ``__origIdx``, are
    carried across the copy.

    Under ``OIN_BORON_CAGE`` (default OFF) a mol carrying a deltahedral boron-cage
    vertex sanitizes with ``SANITIZE_PROPERTIES`` skipped, because a 5- or
    6-connected cage boron trips RDKit's valence *rule* while being a perfectly
    well-formed graph. Everything else about the sanitize still runs, and the
    relaxation reaches only a mol that actually contains the cage motif.
    """
    if _boron_cage_relaxation_applies(mol):
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
            return mol
        except Exception:  # noqa: BLE001 - fall through to the normal ladder
            pass

    try:
        Chem.SanitizeMol(mol)
        return mol
    except Chem.KekulizeException as first_error:
        original_error = first_error

    was_rw_mol = isinstance(mol, Chem.RWMol)
    work = Chem.RWMol(mol)
    try:
        # Perceive rings and aromatic flags without demanding a Kekule structure.
        Chem.SanitizeMol(work, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except Exception:
        pass

    stuck = stuck_ring_atoms(work)
    if not stuck:
        # No genuinely quinoid ring to relax, yet the first full sanitize could
        # not kekulize. This is the stale-flag case: ``AC2mol`` left aromatic
        # flags on a fused/charged ring system that DOES have a valid Kekule
        # structure, so kekulization using those flags failed even though a fresh
        # perception succeeds. ``work`` already re-ran ``SetAromaticity`` above
        # (the partial sanitize), so a full sanitize now kekulizes the refreshed
        # flags. Only a mol that still will not sanitize is an honest limitation.
        try:
            Chem.SanitizeMol(work)
            return Chem.RWMol(work) if was_rw_mol else work.GetMol()
        except Exception:
            pass
        raise OINEncodeError(
            f"cannot kekulize molecule and found no quinoid ring to relax: {original_error}"
        ) from original_error

    clear_ring_aromaticity(work, stuck)
    out = work.GetMol()
    out.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(out)
    except Exception as retry_error:
        raise OINEncodeError(
            f"molecule is still invalid after de-aromatizing the quinoid ring(s) at atoms "
            f"{sorted(stuck)}; the perceived bond orders are unusable ({retry_error})"
        ) from retry_error

    return Chem.RWMol(out) if was_rw_mol else out
