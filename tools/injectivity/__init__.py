"""Injectivity audit (Y1): tools that probe whether the OIN encoder is *injective*.

The round-trip test can only prove ``E∘G`` is a retraction; it can never prove the
encoder ``E`` maps two distinct isomers to two distinct strings. These tools measure
that missing property directly -- by building a twin that differs from a base
structure along exactly one chemically-meaningful axis and asking whether the
encoder (and the round-trip comparison key) can still tell them apart. None of it
invokes the 3D generator, so it is immune to MetalloGen-inaccuracy / timeout noise.

See ``docs/agentic-notes/injectivity/INJECTIVITY_Y1_OVERVIEW.md`` for the falsification frame and results.
"""

from .oracle import OracleVerdict, geometric_chirality, is_distinct_enantiomer, stereo_fingerprint
from .twin_collision import ProbeOutcome, TwinProbe, mirror_z_coords, probe_mirror

__all__ = [
    "OracleVerdict",
    "geometric_chirality",
    "is_distinct_enantiomer",
    "stereo_fingerprint",
    "ProbeOutcome",
    "TwinProbe",
    "mirror_z_coords",
    "probe_mirror",
]
