"""Silent-degradation telemetry for the v0.4.3 elimination study.

The generation pipeline degrades quietly in many places: it substitutes a
flattened bond-order model when both charge solvers fail, ships the best
*rejected* embedding when none passed validation, returns the wrong
stereoisomer when no conformer matched, falls back to the lowest-energy
conformer when the requested geometry was never produced. Each of these is
logged at DEBUG at most, so a distorted structure is indistinguishable from a
clean one in the sweep output.

This module records which of those paths actually fired, per molecule, so the
study can separate "the fallback fired" from "the fallback caused harm".

Design constraints, because the instrument must not perturb what it measures:

* **Disabled by default.** ``record()`` returns immediately unless the
  ``OIN_TELEMETRY`` environment variable is set to ``1``. Production and the
  A/B control arm therefore execute one env-var-backed boolean test per site.
* **Cannot raise.** The body is wrapped in a bare ``except Exception: pass``.
  A bug in the instrument cannot change the control flow of the code under
  observation.
* **Consumes no randomness.** Nothing here touches ``random``, ``numpy.random``
  or RDKit's RNG. Drawing even one number would shift every subsequent
  stochastic decision and silently invalidate the byte-identity proof.
* **Context-local.** State lives in a :class:`~contextvars.ContextVar`, so the
  harness's ``multiprocessing`` spawn workers each accumulate their own events
  with no cross-talk.

Usage at a degradation site -- one inserted line, passing only locals that have
already been computed::

    from ..generation import _telemetry
    _telemetry.record("embed.best_rejected_returned", score=maximum_value)

Usage by a driver::

    with telemetry.collecting() as events:
        generator.generate(oin)
    # events -> [{"site": "...", "detail": {...}}, ...]
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Any, Iterator

__all__ = ["enabled", "record", "collecting", "snapshot", "reset", "counts"]

_ENV_VAR = "OIN_TELEMETRY"

# None means "not collecting"; a list means events are being accumulated.
_events: ContextVar[list[dict[str, Any]] | None] = ContextVar("_oin_telemetry_events", default=None)


def enabled() -> bool:
    """True when telemetry is switched on for this process."""
    return os.environ.get(_ENV_VAR) == "1"


def record(site: str, **detail: Any) -> None:
    """Note that a silent-degradation site fired.

    A no-op unless telemetry is enabled and a collection is active. Never
    raises, never consumes randomness, and never inspects its arguments beyond
    storing them.

    Args:
        site: Stable dotted identifier, e.g. ``"embed.pulp_and_xyz2mol_both_failed"``.
        **detail: Already-computed locals worth keeping (a score, an exception
            type name, a count). Pass values, not expressions that can raise.
    """
    try:
        if os.environ.get(_ENV_VAR) != "1":
            return
        bucket = _events.get()
        if bucket is None:
            return
        bucket.append({"site": site, "detail": detail} if detail else {"site": site})
    except Exception:
        # The instrument must never alter the behaviour of the code it observes.
        pass


@contextmanager
def collecting() -> Iterator[list[dict[str, Any]]]:
    """Collect events emitted inside the block.

    Yields the (live) event list. Restores any enclosing collection on exit, so
    nesting is safe.
    """
    bucket: list[dict[str, Any]] = []
    token = _events.set(bucket)
    try:
        yield bucket
    finally:
        _events.reset(token)


def snapshot() -> list[dict[str, Any]]:
    """Copy of the events collected so far, or an empty list if not collecting."""
    bucket = _events.get()
    return list(bucket) if bucket is not None else []


def reset() -> None:
    """Discard events collected so far, keeping the collection active."""
    bucket = _events.get()
    if bucket is not None:
        bucket.clear()


def counts() -> dict[str, int]:
    """Per-site firing counts for the active collection."""
    out: dict[str, int] = {}
    for event in snapshot():
        site = event.get("site", "?")
        out[site] = out.get(site, 0) + 1
    return out
