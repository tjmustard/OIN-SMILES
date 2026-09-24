"""v0.4.7 gate ARM 2: one molecule's UFF_1-tier round trip, for `tools/gate_v047.sh`.

Replicates the exact XYZ(1) -> OIN(1) -> XYZ(gen) -> OIN(2) pipeline that
``tools/test_dataset_roundtrip.py::_attempt_generation`` runs for the "UFF_1" tier
(``optimizer=None``, ``ensemble_size=1``, ``ff_params=None`` -- FF-relaxed geometry
only, the tier every molecule in the frozen cohort is known to have passed at,
byte-exact, in two independent v0.4.5 runs).

Prints ONE line to stdout:

    name<TAB>sha256(smiles_1)<TAB>sha256(smiles_2)<TAB>len1<TAB>len2<TAB>eta<TAB>xyz_sha256<TAB>status<TAB>honest

THE GATED RE-ENCODE CANNOT SEE WHAT THE GENERATOR BUILDS (v0.4.18)
===================================================================
``smiles_2`` is ``get_oin_string(result.mol, coords)`` -- the generated coordinates serialised
through the GENERATOR'S OWN bond graph, the scoring path the harness dropped in v0.4.8. Whether an
eta ring sits 2.1 A or 2.9 A from its metal never enters it: the metal-ring bonds are read off
``result.mol``. Measured at the v0.4.18 promotion of the eta construction levers: 138 of 198 eta
rows built a DIFFERENT structure than under v0.4.17, and ``smiles_2`` moved on none of them, while
the harness's honest score moved on +197 / -28 of 5,000. ARM 2 gates the encoder and the
generator's graph + stereo; it passed every release that shipped ferrocene detached.

So the row now ends with ``honest``: sha256 of ``XYZToSMILES().convert(<generated xyz>)`` -- an
INDEPENDENT re-perception of the coordinates, the same predicate as the harness's honest score. It
is an OBSERVATION column (nothing compares it; the gate reads columns 1-3 only) until a full run
has frozen it into the goldens. Because it is one more full encode inside the gate's hard timeout,
it is attempted only when the time already spent plus the input encode's own duration leaves a
margin, and says so otherwise (``HONEST_NOT_ATTEMPTED@<s>``); an encode error is reported as
``HONEST_ERROR:<Type>``. Sentinel rows (no structure / error) carry an empty column.

THE GATE OBJECT IS THE STRING, NOT THE XYZ
===========================================
``sha256(smiles_1)`` and ``sha256(smiles_2)`` are what the gate in ``gate_v047.sh``
diffs against the frozen golden manifest. ``xyz_sha256`` (of the generated XYZ) is
recorded as an OBSERVATION column ONLY -- it is strictly stronger than the notation
contract (it would fail a lane that legitimately picked a different-but-equivalent
conformer) and must never be what fails the gate. ``tools/perf_byte_identity_ab.py``
already covers that stronger, XYZ-level check for lanes that want it.

WHY A SEPARATE PROCESS PER MOLECULE
====================================
Invoked once per molecule by ``gate_v047.sh`` (one subprocess per name), not looped
over a shared interpreter. Several of the encoder/generator's ``OIN_*`` levers
(``oin/levers.py``) and module-level constants (e.g.
``generator3d/clash.py:VDW_ACCEPTANCE_ENABLED``) are frozen at *import* time, and the
generator pipeline carries other process-lifetime caches (PuLP topology memo, embed
pool state). A fresh interpreter per molecule is the only isolation guarantee that
does not depend on enumerating every such cache correctly.

ERRORS ARE PART OF THE CONTRACT
================================
Any exception during encode or generate is caught and reported as
``status=ERROR:<Type>:<msg>`` with empty hash columns -- this can never spuriously
match a golden manifest entry (which is always a real success, by the cohort's own
byte-exact-in-both-runs selection predicate), so a regression that turns a pass into
an error is a loud, visible mismatch rather than a silent skip.

Usage
=====
    PYTHONPATH=src .venv/bin/python tools/gate_arm2_roundtrip_one.py \
        --cohort-dir tmCAT-tmPHOTO_xyz_dataset/cohort-v047-slow100 \
        --molecule MUTYEG_comp_0 --timeout 300
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import sys
import tempfile
import time
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

HAPTIC = re.compile(r"\{\d+[<>]\}")


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--molecule", required=True, help="basename without .xyz")
    ap.add_argument("--timeout", type=float, default=300.0, help="generator timeout budget (s)")
    args = ap.parse_args()

    xyz_path = os.path.join(args.cohort_dir, f"{args.molecule}.xyz")
    if not os.path.exists(xyz_path):
        print(f"{args.molecule}\t\t\t-\t-\t-\t\tERROR:FileNotFoundError:{xyz_path}", flush=True)
        return

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import (
        OIN3DGeneratorMetallogen as OIN3DGenerator,
    )

    t_start = time.monotonic()
    try:
        oin1 = XYZToSMILES().convert(xyz_path)
    except Exception as e:
        print(
            f"{args.molecule}\t\t\t-\t-\t-\t\tERROR:{type(e).__name__}:{e}",
            flush=True,
        )
        return

    encode1_s = time.monotonic() - t_start
    eta = "eta" if HAPTIC.search(oin1) else "-"
    tmp_dir = tempfile.mkdtemp(prefix=f"gate_v047_{args.molecule}_")
    try:
        gen = OIN3DGenerator(optimizer=None, ensemble_size=1, timeout=args.timeout, ff_params=None)
        try:
            result = gen.generate(oin1)
        except Exception as e:
            print(
                f"{args.molecule}\t{sha(oin1)}\t\t{len(oin1)}\t-\t{eta}\t\t"
                f"ERROR:{type(e).__name__}:{e}",
                flush=True,
            )
            return

        if result is None or not getattr(result, "xyz", None):
            print(
                f"{args.molecule}\t{sha(oin1)}\t\t{len(oin1)}\t-\t{eta}\t\tERROR:NoResult:generate returned None/empty",
                flush=True,
            )
            return

        xyz_out = result.xyz
        xyz_sha = sha(xyz_out)
        gen_xyz_path = os.path.join(tmp_dir, "gen.xyz")
        with open(gen_xyz_path, "w") as f:
            f.write(xyz_out)

        mol_gen_bonded = result.mol
        oin2 = None
        if mol_gen_bonded is not None:
            try:
                import numpy as np

                from oinsmiles.utils.perception_tmc import get_oin_string

                with open(gen_xyz_path) as f:
                    xyz_lines = f.readlines()
                natoms = int(xyz_lines[0].strip())
                xyz_coords = np.array(
                    [[float(x) for x in xyz_lines[i].split()[1:4]] for i in range(2, 2 + natoms)]
                )
                oin2 = get_oin_string(mol_gen_bonded, xyz_coords)
            except Exception:
                oin2 = None
        if oin2 is None:
            oin2 = XYZToSMILES().convert(gen_xyz_path)

        # the honest observation column -- see the module docstring. The gate's hard timeout is
        # 1.5 x --timeout and a SIGKILLed row scores as a smiles_1 MISMATCH, so an observation
        # must never be what pushes a row over it.
        spent = time.monotonic() - t_start
        if spent + 1.5 * encode1_s < 1.35 * args.timeout:
            try:
                honest = sha(XYZToSMILES().convert(gen_xyz_path))
            except Exception as e:
                honest = f"HONEST_ERROR:{type(e).__name__}"
        else:
            honest = f"HONEST_NOT_ATTEMPTED@{spent:.0f}s"

        print(
            f"{args.molecule}\t{sha(oin1)}\t{sha(oin2)}\t{len(oin1)}\t{len(oin2)}\t{eta}\t"
            f"{xyz_sha}\tOK\t{honest}",
            flush=True,
        )
    except Exception as e:
        print(
            f"{args.molecule}\t{sha(oin1)}\t\t{len(oin1)}\t-\t{eta}\t\t"
            f"ERROR:{type(e).__name__}:{e}",
            flush=True,
        )
        sys.stderr.write(traceback.format_exc())
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
