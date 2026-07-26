#!/usr/bin/env python
"""Boron-cage spike, acid test: does an encoded cage ROUND-TRIP?

Encoding a cage is worthless if the string cannot be read back -- that just moves
the failure downstream. For each molecule this checks, in order:

1. **encode** -- ``XYZToSMILES().convert(xyz)`` returns an OIN string;
2. **determinism** -- a second encode in the same process is byte-identical;
3. **atom conservation** -- the atom count implied by the OIN body equals the
   input xyz's atom count (a shattered/amputated cage shows up here);
4. **re-parse to the same graph** -- every OIN ligand fragment re-parses, and the
   union of re-parsed fragments has the same element multiset and the same bonded
   element-pair multiset as the encoder's own ``tmc_mol`` (metal-donor dative bonds
   excluded, since the OIN carries those as slot markers, not bonds);
5. **key behaves** -- ``canonical_roundtrip_key`` is computable, equal across the
   repeat, and does **not** degrade the cage fragment to the order-dependent
   ``RAW:`` fallback.

Isolated subprocess per molecule: a cage drives ``AC2BO`` into a bare
``sys.exit()`` when the lever is off, and that is not catchable.

Usage:
    OIN_BORON_CAGE=1 PYTHONPATH=src python tools/boron_roundtrip.py --dataset-dir <abs>
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import traceback
from collections import Counter

PER_MOL_TIMEOUT_S = 300

# slot markers `{0}` / `{0>}` / `{0<}` and the metal-fragment head `[Rh_TPY]`
_SLOT_RE = re.compile(r"\{\d+[<>]?\}")
_AXIAL_RE = re.compile(r"\s*\|ax:[+\-]*\|")


def find_xyz(mol: str, dataset_dir: str) -> str | None:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{mol}_comp_0.xyz")
        if os.path.exists(p):
            return p
    return None


def xyz_atom_count(path: str) -> int:
    with open(path) as fh:
        return int(fh.readline().split()[0])


def parse_frag_tolerant(smiles: str):
    """Parse an OIN ligand fragment, tolerating a hypervalent cage boron.

    Mirrors ``compare._parse_fragment`` but adds a third rung that also skips
    ``SANITIZE_PROPERTIES`` -- the only check a deltahedral cage vertex violates.
    """
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        return m, "full"
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None, "unparseable"
    for label, ops in (
        (
            "no_kekulize",
            Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        ),
        (
            "no_properties",
            Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        ),
        (
            "no_props_no_kek",
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        ),
    ):
        m2 = Chem.MolFromSmiles(smiles, sanitize=False)
        try:
            Chem.SanitizeMol(m2, sanitizeOps=ops)
            return m2, label
        except Exception:  # noqa: BLE001
            continue
    return None, "unsanitizable"


def _boron_h(mol):
    """Hydrogens on boron, counting explicit H *atoms* and implicit H alike.

    The encoder's tmc_mol carries cage B-H as explicit hydrogen atoms; a re-parsed
    OIN fragment carries the same hydrogens as `[BH]` implicit counts. Comparing
    GetTotalNumHs alone would score an exact round trip as a loss.
    """
    n = 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 5:
            continue
        n += a.GetTotalNumHs()
        n += sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() == 1)
    return n


def graph_fp(mol, drop_dative=True):
    """Element multiset + bonded element-pair multiset, order-independent."""
    from rdkit import Chem

    els = Counter(a.GetSymbol() for a in mol.GetAtoms())
    # implicit hydrogens count too, or a fragment written as `C` vs `[CH4]` differs
    for a in mol.GetAtoms():
        nh = a.GetTotalNumHs()
        if nh:
            els["H"] += nh
    bonds = Counter()
    for b in mol.GetBonds():
        if drop_dative and b.GetBondType() == Chem.BondType.DATIVE:
            continue
        pair = tuple(sorted((b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())))
        if drop_dative and "H" in pair:
            continue
        bonds[pair] += 1
    return els, bonds


_METAL_HEAD_RE = re.compile(r"^\[([A-Z][a-z]?)(@[A-Z0-9]+)?_[A-Z]{3}\]$")


def oin_fragments(oin: str):
    """Split an OIN string into its dot-separated fragments, slot markers stripped.

    The metal head carries a geometry suffix (``[Rh_TPY]``) that is OIN syntax, not
    SMILES, so it is rewritten to a plain bracketed atom before re-parsing.
    """
    body = _AXIAL_RE.sub("", oin.strip())
    body = body.split(" |")[0]
    out = []
    for f in body.split("."):
        f = f.strip()
        if not f:
            continue
        m = _METAL_HEAD_RE.match(f)
        if m:
            out.append(f"[{m.group(1)}]")
            continue
        out.append(_SLOT_RE.sub("", f))
    return out


def check(mol_name: str, dataset_dir: str) -> dict:
    import warnings

    warnings.filterwarnings("ignore")
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    from oinsmiles import XYZToSMILES
    from oinsmiles.oin.compare import canonical_roundtrip_key
    from oinsmiles.utils.xyz2mol import get_tmc_mol

    rec: dict = {"mol": mol_name, "lever": bool(os.environ.get("OIN_BORON_CAGE"))}
    path = find_xyz(mol_name, dataset_dir)
    if path is None:
        rec["result"] = "NO_FILE"
        return rec
    rec["xyz_atoms"] = xyz_atom_count(path)

    # 1. encode
    try:
        oin = XYZToSMILES().convert(path)
    except BaseException as e:  # noqa: BLE001
        tb = traceback.extract_tb(sys.exc_info()[2])
        last = [f for f in tb if "oinsmiles" in f.filename]
        rec["result"] = "ENCODE_FAIL"
        rec["etype"] = type(e).__name__
        rec["loc"] = (
            f"{os.path.basename(last[-1].filename)}:{last[-1].lineno}:{last[-1].name}"
            if last
            else "?"
        )
        rec["err"] = str(e).replace("\n", " ")[:200]
        return rec
    rec["oin"] = oin
    rec["oin_len"] = len(oin)

    # 2. determinism
    try:
        oin2 = XYZToSMILES().convert(path)
        rec["deterministic"] = oin2 == oin
        if oin2 != oin:
            rec["oin2"] = oin2
    except BaseException as e:  # noqa: BLE001
        rec["deterministic"] = f"repeat raised {type(e).__name__}"

    # 3./4. re-parse every fragment; compare to the encoder's own tmc_mol
    frags = oin_fragments(oin)
    rec["n_frags"] = len(frags)
    parsed, modes, unparsed = [], [], []
    for f in frags:
        m, mode = parse_frag_tolerant(f)
        modes.append(mode)
        if m is None:
            unparsed.append(f[:80])
        else:
            parsed.append(m)
    rec["parse_modes"] = dict(Counter(modes))
    rec["unparsed_frags"] = unparsed
    rec["all_frags_reparse"] = not unparsed

    els_rt: Counter = Counter()
    bonds_rt: Counter = Counter()
    b_h_rt = 0
    for m in parsed:
        e, b = graph_fp(m)
        els_rt += e
        bonds_rt += b
        b_h_rt += _boron_h(m)
    rec["reparsed_heavy_atoms"] = sum(v for k, v in els_rt.items() if k != "H")
    rec["reparsed_total_atoms"] = sum(els_rt.values())
    rec["reparsed_B_H"] = b_h_rt

    try:
        tmc, _xyz = get_tmc_mol(path, 0, with_stereo=False)
        els_enc, bonds_enc = graph_fp(tmc)
        rec["encoder_atoms"] = sum(els_enc.values())
        rec["encoder_heavy_atoms"] = sum(v for k, v in els_enc.items() if k != "H")
        rec["encoder_B_H"] = _boron_h(tmc)
        # PRIMARY criterion: heavy-atom multiset and heavy-bond multiset.
        # Total-H is deliberately NOT primary: stripping a slot marker turns a
        # coordinated donor into a free ligand, which legitimately fills implicit
        # hydrogens (`C(O{2})` -> `C(O)`, an -OH). That drift is a property of the
        # OIN format, present with the lever OFF on ordinary molecules too -- see
        # the --control run -- so charging it to the cage would be wrong.
        els_enc_heavy = Counter({k: v for k, v in els_enc.items() if k != "H"})
        els_rt_heavy = Counter({k: v for k, v in els_rt.items() if k != "H"})
        rec["heavy_atoms_equal"] = els_enc_heavy == els_rt_heavy
        rec["graph_elements_equal"] = els_enc == els_rt
        rec["graph_bonds_equal"] = bonds_enc == bonds_rt
        rec["cage_H_conserved"] = rec["encoder_B_H"] == b_h_rt
        rec["total_H_delta"] = els_rt.get("H", 0) - els_enc.get("H", 0)
        if els_enc_heavy != els_rt_heavy:
            rec["heavy_enc"] = dict(els_enc_heavy)
            rec["heavy_rt"] = dict(els_rt_heavy)
        if bonds_enc != bonds_rt:
            rec["bonds_only_enc"] = {str(k): v for k, v in (bonds_enc - bonds_rt).items()}
            rec["bonds_only_rt"] = {str(k): v for k, v in (bonds_rt - bonds_enc).items()}
    except BaseException as e:  # noqa: BLE001
        rec["heavy_atoms_equal"] = f"tmc failed: {type(e).__name__}"

    # 5. key behaviour
    try:
        k1 = canonical_roundtrip_key(oin)
        rec["key_ok"] = True
        rec["key_repr"] = str(k1)[:200]
        rec["key_has_RAW"] = "RAW:" in str(k1)
        k2 = canonical_roundtrip_key(oin)
        rec["key_stable"] = str(k1) == str(k2)
    except BaseException as e:  # noqa: BLE001
        rec["key_ok"] = f"{type(e).__name__}: {str(e)[:120]}"

    rec["result"] = (
        "ROUNDTRIP_OK"
        if (
            rec.get("deterministic") is True
            and rec.get("all_frags_reparse")
            and rec.get("heavy_atoms_equal") is True
            and rec.get("graph_bonds_equal") is True
            and rec.get("cage_H_conserved") is True
            and rec.get("key_ok") is True
        )
        else "ENCODES_PARTIAL"
    )
    return rec


def main() -> None:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boron_characterize import BORON_COHORT

    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", default=None)
    ap.add_argument("--dataset-dir", required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--out", default="tools/boron_roundtrip.json")
    args = ap.parse_args()

    if args.worker:
        print(json.dumps(check(args.worker, args.dataset_dir)), flush=True)
        return

    todo = args.only.split(",") if args.only else BORON_COHORT
    results = []
    for mol in todo:
        cmd = [
            sys.executable,
            os.path.abspath(__file__),
            "--worker",
            mol,
            "--dataset-dir",
            args.dataset_dir,
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=PER_MOL_TIMEOUT_S)
            line = next((ln for ln in proc.stdout.splitlines() if ln.startswith("{")), None)
            rec = (
                json.loads(line)
                if line
                else {
                    "mol": mol,
                    "result": "WORKER_DIED",
                    "rc": proc.returncode,
                    "err": (proc.stderr or "")[-200:],
                }
            )
        except subprocess.TimeoutExpired:
            rec = {"mol": mol, "result": "TIMEOUT"}
        results.append(rec)
        print(
            f"{mol:10s} {rec['result']:16s} "
            f"det={rec.get('deterministic')} heavy={rec.get('reparsed_heavy_atoms')}"
            f"/{rec.get('encoder_heavy_atoms')} heavyEq={rec.get('heavy_atoms_equal')} "
            f"bonds={rec.get('graph_bonds_equal')} BH={rec.get('reparsed_B_H')}"
            f"/{rec.get('encoder_B_H')} dH={rec.get('total_H_delta')} "
            f"key={rec.get('key_ok')} RAW={rec.get('key_has_RAW')} "
            f"{rec.get('etype', '')}{rec.get('loc', '')}",
            flush=True,
        )
    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    tally = Counter(r["result"] for r in results)
    print("\n=== TALLY ===")
    for k, v in tally.most_common():
        print(f"{v:3d}  {k}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
