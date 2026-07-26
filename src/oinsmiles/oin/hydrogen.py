"""H-faithful SMILES serialization.

RDKit's SMILES writer emits a **bare** organic-subset symbol (``C``, ``n``, ``S``, ...)
whenever it judges brackets unnecessary. A bare symbol is read back as *"fill to the
next allowed valence with hydrogen"*, and ``Atom.SetNoImplicit(True)`` does **not**
force a bracket. So an atom the caller has deliberately set to 0 H serializes bare and
re-reads carrying a hydrogen it never had:

* a thiophene sulfur whose perceived valence is 3 -- strictly between sulfur's allowed
  2 and 4 -- writes ``S`` and re-reads as ``[SH]`` (``CIDDAU_comp_0``, 18 -> 19 atoms);
* a correctly-bracketed ``[C]`` (an NHC carbene donor, a benzylic carbon whose
  hydrogens the input XYZ never had) is silently **de-bracketed** to ``C`` and re-reads
  as ``CH2`` (``INENOF_comp_0``, 58 -> 60).

For a lossless notation this matters directly: the phantom hydrogen is baked into the
OIN string, so the 3D generator faithfully builds a molecule with the wrong number of
atoms and the round trip fails its final gate. It was the largest non-timeout
``hard_fail`` class in the v0.4.5 capstone corpus -- see ``docs/ATOM_COUNT_v0.4.5.md``.

The per-motif alternative has been tried repeatedly and does not converge: the
bare-donor strip heuristics in ``generation/metallogen_adapter.py`` and step 1b of
``utils/oin_aligner.py`` are both compensations for this same asymmetry, each added for
the molecules that happened to be in front of someone at the time (COLWIK, ACOXEX,
ARONEA, BOXJUU). :func:`h_faithful_smiles` instead enforces the property those
heuristics were approximating, and enforces it by *checking* rather than by predicting
which motifs need it.

Gated behind ``OIN_H_FAITHFUL``; with the lever unset the output is byte-identical to
plain ``MolToSmiles``.
"""

from __future__ import annotations

import os

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - rdkit is a hard dependency in practice
    Chem = None  # type: ignore[assignment]

__all__ = ["h_faithful_smiles", "hydrogen_faithfulness_enabled"]


def hydrogen_faithfulness_enabled() -> bool:
    """True when the ``OIN_H_FAITHFUL`` lever is set.

    Read on every call rather than cached at import, so a test can toggle it.
    """
    return bool(os.environ.get("OIN_H_FAITHFUL"))


def _reparse(smiles: str):
    """Parse `smiles` the way every downstream consumer of an OIN fragment does.

    Full sanitize first; on failure retry without kekulization, which is the one step
    a metal-stripped aromatic donor ring reliably fails. Returns ``None`` if the
    fragment cannot be read at all (a borane cluster, ``[C]#O``), in which case the
    caller must leave the string alone -- an unparseable fragment tells us nothing
    about its hydrogen counts.
    """
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
    except Exception:
        return None
    return mol


def _total_h(mol) -> int:
    """Total hydrogen of `mol`, counting both H nodes and folded H counts."""
    total = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            total += 1
        else:
            total += atom.GetTotalNumHs()
    return total


def _output_order(mol):
    try:
        return list(
            mol.GetPropsAsDict(includePrivate=True, includeComputed=True)["_smilesAtomOutputOrder"]
        )
    except Exception:
        return None


def _cached_copy(mol):
    """A property-cached copy of `mol`, or ``None``.

    Callers hand us molecules in wildly different states. ``oin/inline.py`` in
    particular passes the result of ``MolFromSmiles(..., sanitize=False)`` with no
    ``UpdatePropertyCache`` at all, and ``GetTotalNumHs()`` on such an atom does not
    return a wrong answer -- it raises a RDKit pre-condition violation. Work on a copy
    so the caller's molecule is never mutated by our bookkeeping.
    """
    try:
        ref = Chem.Mol(mol)
        ref.UpdatePropertyCache(strict=False)
        return ref
    except Exception:
        return None


def _divergent_atoms(ref, mol, smiles: str) -> list[int]:
    """Indices whose H count changes when `smiles` is read back.

    ``MolToSmiles`` records ``_smilesAtomOutputOrder`` as a side effect:
    ``order[canonical_position] = original atom index``. That is what makes the
    comparison exact rather than a heuristic match.
    """
    order = _output_order(mol)
    if order is None:
        return []
    back = _reparse(smiles)
    if back is None or back.GetNumAtoms() != ref.GetNumAtoms():
        return []
    bad = []
    for pos, orig in enumerate(order):
        want = ref.GetAtomWithIdx(int(orig))
        got = back.GetAtomWithIdx(pos)
        if want.GetAtomicNum() != got.GetAtomicNum():
            return []  # ordering assumption violated; do not touch the string
        if want.GetTotalNumHs() != got.GetTotalNumHs():
            bad.append(int(orig))
    return bad


def h_faithful_smiles(mol, **kwargs) -> str:
    """``Chem.MolToSmiles(mol, **kwargs)``, but the result re-reads with the same H.

    Writes the SMILES, reads it back, and compares hydrogen counts atom by atom. Any
    atom that came back different is given an unpaired electron -- purely to force
    RDKit to bracket it, since SMILES has no radical syntax and so nothing about the
    radical reaches the string -- and the molecule is re-serialized and re-checked.

    Returns the plain ``MolToSmiles`` output unchanged when it was already faithful,
    when the lever is off, when the fragment cannot be re-parsed, or when the repair
    fails to verify. It can therefore only ever change a string that was **measurably
    wrong**, which bounds the blast radius to molecules whose atom count the notation
    was already getting wrong.

    `mol` is not mutated.
    """
    smiles = Chem.MolToSmiles(mol, **kwargs)
    if not hydrogen_faithfulness_enabled():
        return smiles
    try:
        return _repair(mol, smiles, kwargs)
    except Exception:
        # This runs inside the encoder's serialization hot path, and `oin/inline.py`
        # wraps its caller in a bare `except Exception` that falls back to a *different*
        # slot-tagging strategy. An exception escaping here would therefore not surface
        # as an error -- it would silently reroute the encoder. Never raise.
        return smiles


def _repair(mol, smiles: str, kwargs: dict) -> str:
    ref = _cached_copy(mol)
    if ref is None:
        return smiles

    bad = _divergent_atoms(ref, mol, smiles)
    if not bad:
        return smiles

    patched = Chem.RWMol(ref)
    for idx in bad:
        atom = patched.GetAtomWithIdx(idx)
        # Freeze the intended count first: a radical on an atom that is still relying on
        # implicit H would change the very number we are trying to preserve.
        atom.SetNumExplicitHs(ref.GetAtomWithIdx(idx).GetTotalNumHs())
        atom.SetNoImplicit(True)
        if atom.GetNumRadicalElectrons() == 0:
            atom.SetNumRadicalElectrons(1)
    candidate_mol = patched.GetMol()
    candidate = Chem.MolToSmiles(candidate_mol, **kwargs)

    # The atom ORDER must not move. Callers read `_smilesAtomOutputOrder` off the
    # molecule they passed us -- `xyz2mol.get_oin_string` uses it to decide which
    # character position each `{slot}` marker attaches to -- and that property records
    # the FIRST serialization, not the repaired one. An unpaired electron can in
    # principle change the canonical ranking, and if it did, every slot marker would
    # land on the wrong atom: a silently wrong OIN string, which is far worse than the
    # phantom hydrogen we are removing. So decline the repair rather than re-home the
    # property. Coverage lost this way is visible as a still-failing molecule; a
    # mis-slotted string would not be visible at all.
    if _output_order(candidate_mol) != _output_order(mol):
        return smiles

    # Verify rather than assume. The repair is only accepted if the rewritten string
    # actually reads back with the hydrogen budget the caller intended -- otherwise we
    # would be trading one wrong count for another, which is how a "fix" moves a failure
    # instead of closing it.
    back = _reparse(candidate)
    if back is None or _total_h(back) != _total_h(ref):
        return smiles
    return candidate
