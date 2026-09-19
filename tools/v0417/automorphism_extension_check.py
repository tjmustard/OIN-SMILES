"""Does every automorphism found on the PRUNED fragment extend to the FULL fragment?

``canonical_slots._automorphism_query`` strips pendant decoration (t-Bu, CF3, alkyl tails) before
enumerating automorphisms, and pins each surviving atom's full-graph symmetry class as an isotope
so the stripped graph "cannot gain a symmetry". That last clause is an ARGUMENT -- colour
refinement decides isomorphism of pendant trees -- not a proof about RDKit's ranking. The exact
fold's whole safety claim ("every candidate describes the molecule in hand") rests on it, so it is
tested here on every distinct ligand fragment the cohort's strings contain.

For each pruned self-match, the core atoms are pinned to their images with unique isotopes on an
UNPRUNED copy and RDKit is asked for ONE full match (existence only -- nothing is enumerated).
A pruned automorphism with no full extension is a symmetry the ligand does not have: print it.

A broken checker would print 0 failures too, so the NEGATIVE CONTROL pins one core atom to a
wrong-class image and requires the extension to FAIL on every fragment it is tried on.

    PYTHONPATH=$PWD/src <main>/.venv/bin/python tools/v0417/automorphism_extension_check.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

MAIN = Path("/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset")
PER_FRAGMENT = 64  # distinct donor actions checked per fragment; most fragments have <= 8


def _full_skeleton(mol, resonance):
    """The unpruned graph ``_automorphism_query`` flattens, rebuilt by asking it to keep all."""
    from oinsmiles.oin.canonical_slots import _automorphism_query

    return _automorphism_query(mol, resonance, set(range(mol.GetNumAtoms())))


def _extends(full, core_old, images_old, ranks_iso):
    """Is there a full self-match sending ``core_old[k] -> images_old[k]`` for every k?"""
    from rdkit import Chem

    q, t = Chem.RWMol(full), Chem.RWMol(full)
    for a in q.GetAtoms():
        a.SetIsotope(ranks_iso[a.GetIdx()])
    for a in t.GetAtoms():
        a.SetIsotope(ranks_iso[a.GetIdx()])
    for k, (src, dst) in enumerate(zip(core_old, images_old)):
        q.GetAtomWithIdx(src).SetIsotope(100000 + k)
        t.GetAtomWithIdx(dst).SetIsotope(100000 + k)
    return t.GetMol().HasSubstructMatch(q.GetMol(), useChirality=True)


def main():
    from rdkit import RDLogger

    import oinsmiles
    from oinsmiles.oin.canonical_slots import _AUTOMORPHISM_MATCH_CAP, _automorphism_query
    from oinsmiles.oin.compare import _parse_fragment
    from oinsmiles.oin.inline import OINInlineHandler, _count_smiles_atoms_before

    RDLogger.DisableLog("rdApp.*")
    print("oinsmiles from:", oinsmiles.__file__)
    src = MAIN / "results-census/e_selfconsistency.jsonl"
    frags = set()
    for ln in src.read_text().splitlines():
        if ln.startswith("{"):
            base = json.loads(ln)["base"]
            if base:
                frags.update(f for f in base.split(".") if f)
    frags = sorted(f for f in frags if not OINInlineHandler.METAL_REGEX.search(f))
    print(f"distinct non-metal fragments in the cohort's strings: {len(frags)}")

    tally, failures, neg = Counter(), [], Counter()
    for frag in frags:
        slot_atoms: dict[int, set] = {}
        for m in OINInlineHandler.SLOT_REGEX.finditer(frag):
            prefix = OINInlineHandler.SLOT_REGEX.sub("", frag[: m.start()])
            slot_atoms.setdefault(int(m.group(1)), set()).add(
                _count_smiles_atoms_before(prefix, len(prefix))
            )
        if len(slot_atoms) < 2:
            tally["skipped: fewer than 2 slots (never folded)"] += 1
            continue
        mol = _parse_fragment(OINInlineHandler.SLOT_REGEX.sub("", frag))
        donors = {i for atoms in slot_atoms.values() for i in atoms}
        if mol is None or any(i >= mol.GetNumAtoms() for i in donors):
            tally["skipped: unparseable"] += 1
            continue
        for resonance in (True, False):
            try:
                pruned, order = _automorphism_query(mol, resonance, donors)
                full, full_order = _full_skeleton(mol, resonance)
                matches = pruned.GetSubstructMatches(
                    pruned, uniquify=False, useChirality=True, maxMatches=_AUTOMORPHISM_MATCH_CAP
                )
            except Exception as exc:  # noqa: BLE001
                tally[f"skipped: {type(exc).__name__}"] += 1
                continue
            if len(matches) >= _AUTOMORPHISM_MATCH_CAP:
                tally["skipped: match cap"] += 1
                continue
            assert full_order == list(range(mol.GetNumAtoms())), "keep-all must keep all"
            iso = [a.GetIsotope() for a in full.GetAtoms()]
            tally["fragments checked"] += 1
            tally["atoms stripped"] += mol.GetNumAtoms() - len(order)
            pos = {old: k for k, old in enumerate(order)}
            by_action = {}
            for mt in matches:
                by_action.setdefault(tuple(order[mt[pos[d]]] for d in sorted(donors)), mt)
            for mt in list(by_action.values())[:PER_FRAGMENT]:
                tally["pruned automorphisms checked"] += 1
                images = [order[j] for j in mt]
                if not _extends(full, order, images, iso):
                    failures.append((frag, resonance, images))
            # NEGATIVE CONTROL: send the first donor to an atom of a DIFFERENT class.
            d0 = sorted(donors)[0]
            wrong = next((i for i in order if iso[i] != iso[d0]), None)
            if wrong is not None:
                images = list(order)
                images[pos[d0]] = wrong
                neg[_extends(full, order, images, iso)] += 1

    for k, v in tally.most_common():
        print(f"   {v:7d}  {k}")
    print(f"\nPRUNED AUTOMORPHISMS WITH NO FULL EXTENSION: {len(failures)}")
    for frag, resonance, _images in failures[:20]:
        print(f"   resonance={resonance}  {frag}")
    print(f"NEGATIVE CONTROL (wrong-class pin) extended: {neg[True]}  refused: {neg[False]}")
    out = MAIN / "results-v0.4.17-exactfold/automorphism_extension_check.json"
    out.write_text(
        json.dumps(
            {
                "tally": dict(tally),
                "failures": [[f, r] for f, r, _ in failures],
                "negative_control": {"extended": neg[True], "refused": neg[False]},
            },
            indent=1,
        )
    )
    print("wrote", out)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
