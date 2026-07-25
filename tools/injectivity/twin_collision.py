"""Twin-collision harness: perturb one axis, encode both twins, detect a collision.

Given a base structure and a perturbation *operator* that changes exactly one
chemically-meaningful axis, encode both with the OIN encoder and record:

* ``raw_equal``    -- the two raw OIN strings are byte-identical
* ``key_equal``    -- the round-trip equivalence key agrees (what the BATCH harness gates on)
* ``oracle_distinct`` -- the independent oracle says the twin is a different isomer

A collision is ``oracle_distinct and key_equal``: two genuinely different isomers the
round-trip test would call the same. If ``raw_equal`` too, the encoder is *totally* blind
(not even the raw string separates them). No 3D generator is invoked anywhere here, so the
measurement is a pure property of the encoder ``E`` -- immune to MetalloGen / timeout noise.

Run ad-hoc:  ``PYTHONPATH=$PWD/src python -m tools.injectivity.twin_collision <fix.xyz> ...``
"""

from __future__ import annotations

import contextlib
import os
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key

from .oracle import OracleVerdict, is_distinct_enantiomer

# --- classification of a twin pair --------------------------------------------------

#: Verdicts, worst first. The two blind verdicts are round-trip FALSE POSITIVES (the
#: batch gate is the key, and both have key_equal=True); ``encoder_blind`` is strictly
#: worse because not even the raw string separates the isomers.
VERDICT_ENCODER_BLIND = "encoder_blind"  # oracle-distinct, raw_equal  -> total blindness
VERDICT_KEY_BLIND = "key_blind"  # oracle-distinct, raw differs, key_equal -> batch false-positive
VERDICT_DISTINGUISHED = "distinguished"  # oracle-distinct, key differs -> encoder injective here
VERDICT_OVER_SENSITIVE = "over_sensitive"  # NOT oracle-distinct yet raw differs -> false negative
VERDICT_INVARIANT_OK = "invariant_ok"  # NOT oracle-distinct, raw_equal -> correct invariance

_SEVERITY = {
    VERDICT_ENCODER_BLIND: 4,
    VERDICT_KEY_BLIND: 3,
    VERDICT_OVER_SENSITIVE: 2,
    VERDICT_DISTINGUISHED: 0,
    VERDICT_INVARIANT_OK: 0,
}


def classify(raw_equal: bool, key_equal: bool, oracle_distinct: bool) -> str:
    if oracle_distinct:
        if raw_equal:
            return VERDICT_ENCODER_BLIND
        return VERDICT_KEY_BLIND if key_equal else VERDICT_DISTINGUISHED
    return VERDICT_INVARIANT_OK if raw_equal else VERDICT_OVER_SENSITIVE


@dataclass
class ProbeOutcome:
    name: str
    operator: str
    oin_base: str
    oin_twin: str
    raw_equal: bool
    key_equal: bool
    oracle_distinct: bool
    oracle_rmsd: float
    oracle_note: str
    verdict: str
    fingerprint: dict

    @property
    def severity(self) -> int:
        return _SEVERITY.get(self.verdict, 0)

    @property
    def is_collision(self) -> bool:
        """A genuinely different isomer the round-trip key cannot tell from the base."""
        return self.verdict in (VERDICT_ENCODER_BLIND, VERDICT_KEY_BLIND)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["severity"] = self.severity
        d["is_collision"] = self.is_collision
        return d


@dataclass
class TwinProbe:
    """A named probe: a base fixture plus the axis its mirror is meant to flip."""

    name: str
    xyz_path: str
    axis: str  # human label, e.g. "metal Δ/Λ", "axial (atropisomer)"
    operator: str = "mirror_z"
    charge: int = 0


def _key(oin: str):
    return winding_canonical_key(normalize_oin_for_comparison(oin))


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def mirror_z_coords(lines: list[str]) -> list[str]:
    """Return XYZ lines with every z coordinate negated (a mirror-image operator)."""
    out = [lines[0], lines[1]]
    for line in lines[2:]:
        p = line.split()
        if len(p) >= 4:
            out.append(
                f"{p[0]:<3} {float(p[1]):>14.8f} {float(p[2]):>14.8f} {-float(p[3]):>14.8f}\n"
            )
        elif line.strip():
            out.append(line if line.endswith("\n") else line + "\n")
    return out


def _write_mirror(xyz_path: str | Path) -> str:
    with open(xyz_path) as f:
        lines = f.readlines()
    mirrored = mirror_z_coords(lines)
    fd, tmp = tempfile.mkstemp(suffix=".xyz")
    with os.fdopen(fd, "w") as f:
        f.writelines(mirrored)
    return tmp


def probe_mirror(xyz_path: str | Path, *, name: str | None = None, charge: int = 0) -> ProbeOutcome:
    """Encode a structure and its z-mirror; classify the collision, if any."""
    name = name or Path(xyz_path).stem
    conv = XYZToSMILES()
    mirror = _write_mirror(xyz_path)
    try:
        with _silence_fds():
            oin_base = conv.convert(str(xyz_path))
            oin_twin = conv.convert(mirror)
    finally:
        os.unlink(mirror)
    verdict: OracleVerdict = is_distinct_enantiomer(xyz_path, charge=charge)
    raw_equal = oin_base == oin_twin
    key_equal = _key(oin_base) == _key(oin_twin)
    return ProbeOutcome(
        name=name,
        operator="mirror_z",
        oin_base=oin_base,
        oin_twin=oin_twin,
        raw_equal=raw_equal,
        key_equal=key_equal,
        oracle_distinct=verdict.distinct,
        oracle_rmsd=round(verdict.rmsd, 4),
        oracle_note=verdict.note,
        verdict=classify(raw_equal, key_equal, verdict.distinct),
        fingerprint=verdict.fingerprint,
    )


def run_probes(probes: list[TwinProbe]) -> list[ProbeOutcome]:
    outcomes = []
    for p in probes:
        try:
            outcomes.append(probe_mirror(p.xyz_path, name=p.name, charge=p.charge))
        except Exception as e:  # keep the sweep going; record the failure
            outcomes.append(
                ProbeOutcome(
                    name=p.name,
                    operator=p.operator,
                    oin_base="",
                    oin_twin="",
                    raw_equal=False,
                    key_equal=False,
                    oracle_distinct=False,
                    oracle_rmsd=float("nan"),
                    oracle_note=f"probe error: {e}",
                    verdict="error",
                    fingerprint={},
                )
            )
    return outcomes


def _fmt(o: ProbeOutcome) -> str:
    tag = {
        VERDICT_ENCODER_BLIND: "ENCODER-BLIND (total)",
        VERDICT_KEY_BLIND: "KEY-BLIND (batch false-positive)",
        VERDICT_DISTINGUISHED: "distinguished (injective)",
        VERDICT_INVARIANT_OK: "invariant ok (mirror = same isomer)",
        VERDICT_OVER_SENSITIVE: "OVER-SENSITIVE (false negative)",
    }.get(o.verdict, o.verdict)
    return (
        f"\n## {o.name}  [{o.operator}]\n"
        f"  base : {o.oin_base}\n"
        f"  twin : {o.oin_twin}\n"
        f"  oracle: distinct={o.oracle_distinct} (mirror RMSD {o.oracle_rmsd} Å){' -- ' + o.oracle_note if o.oracle_note else ''}\n"
        f"  raw_equal={o.raw_equal}  key_equal={o.key_equal}  ->  {tag}"
    )


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    for path in argv:
        print(_fmt(probe_mirror(path)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
