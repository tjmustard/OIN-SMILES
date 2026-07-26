"""Deterministic-counter attribution for the v0.4.5 perf lane.

Runs a single molecule's OIN-DIRECT generation step (the encode is fast; the
generation/embed step is the slow tail this lane owns) with the low-level
engine calls monkeypatched to *count*, not to change behaviour. Counters are
exact regardless of host contention -- unlike wall-clock, which is meaningless
above load ~12 on this box (see spec/handoffs/v0.4.5/RESTART.md SS6).

Counted:
  * ``AllChem.EmbedMolecule`` calls (total, and how many returned -1)
  * ``embed.get_embedding`` calls == outer pool attempts consumed
  * ``embed.get_embeddings_batch`` calls (should be 0 -- num_threads=1 default)
  * CBC/PuLP ``actualSolve`` calls == real solver subprocess-equivalent spawns
  * PuLP topology memo hits/misses (``pulp_cache_stats()``)
  * final pool size (``len(mols)`` returned by ``generate_3d_structures``)
  * ``metallogen_adapter.build_contract_mol`` / ``_reencode_key_matches`` /
    ``_reencode_oin_fast`` / ``_reencode_oin`` calls -- the SL1 accept-first re-encode
    path. This is where the dominant slow-tail cost actually lives (see
    ``docs/PERF_v0.4.5.md``): ``_reencode_oin`` is a full ``XYZToSMILES().convert()``
    re-perception, measured at 48-57s/call on an eta test case, and per-call stderr
    timing lines are printed for each of these so a single run shows exactly which
    calls are expensive and which are cache hits.

Usage:
    PYTHONPATH=src .venv/bin/python tools/perf_attribute.py \
        --dataset /path/to/tmCAT-tmPHOTO_xyz_dataset \
        --molecule QIDKUL_comp_0 --timeout 60
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def find_xyz(dataset_dir: str, molecule: str) -> str:
    for sub in ("cat", "photo"):
        p = os.path.join(dataset_dir, sub, f"{molecule}.xyz")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"{molecule} not found under {dataset_dir}/{{cat,photo}}")


def atom_count(xyz_path: str) -> int:
    with open(xyz_path) as f:
        return int(f.readline().strip())


class Counters:
    def __init__(self):
        self.embed_calls = 0
        self.embed_rc_minus1 = 0
        self.embed_raised = 0
        self.get_embedding_calls = 0
        self.get_embeddings_batch_calls = 0
        self.actual_solve_calls = 0
        self.get_valid_molecule_calls = 0
        self.clean_geometry_calls = 0
        self.build_contract_mol_calls = 0
        self.reencode_key_matches_calls = 0
        self.reencode_oin_fast_calls = 0
        self.reencode_oin_full_calls = 0

    def as_dict(self):
        return dict(self.__dict__)


def instrument(counters: Counters):
    """Monkeypatch the hot call sites; returns a list of (obj, attr, original) to restore."""
    import pulp as pl
    from rdkit.Chem import AllChem

    from oinsmiles.generator3d import chem as chem_mod
    from oinsmiles.generator3d import clean_geometry as cg_mod
    from oinsmiles.generator3d import embed as embed_mod

    restores = []

    orig_embed_molecule = AllChem.EmbedMolecule

    def counted_embed_molecule(*a, **kw):
        counters.embed_calls += 1
        try:
            rc = orig_embed_molecule(*a, **kw)
        except Exception:
            counters.embed_raised += 1
            raise
        if rc == -1:
            counters.embed_rc_minus1 += 1
        return rc

    AllChem.EmbedMolecule = counted_embed_molecule
    restores.append((AllChem, "EmbedMolecule", orig_embed_molecule))
    # embed.py imported `from rdkit.Chem import AllChem` -- same module object,
    # so patching the attribute above is visible there too. No separate patch
    # needed unless a call site imported the *function* by name.

    orig_get_embedding = embed_mod.get_embedding
    _t_start = time.monotonic()

    def counted_get_embedding(*a, **kw):
        counters.get_embedding_calls += 1
        t0 = time.monotonic()
        result = orig_get_embedding(*a, **kw)
        dt = time.monotonic() - t0
        print(
            f"[progress] get_embedding call#{counters.get_embedding_calls} "
            f"took {dt:.2f}s (t+{time.monotonic() - _t_start:.1f}s) "
            f"embed_calls_total={counters.embed_calls} rc-1={counters.embed_rc_minus1}",
            file=sys.stderr,
            flush=True,
        )
        return result

    embed_mod.get_embedding = counted_get_embedding
    restores.append((embed_mod, "get_embedding", orig_get_embedding))

    orig_get_embeddings_batch = embed_mod.get_embeddings_batch

    def counted_get_embeddings_batch(*a, **kw):
        counters.get_embeddings_batch_calls += 1
        return orig_get_embeddings_batch(*a, **kw)

    embed_mod.get_embeddings_batch = counted_get_embeddings_batch
    restores.append((embed_mod, "get_embeddings_batch", orig_get_embeddings_batch))

    # generator3d/__init__.py calls embed.get_embedding via `from . import embed as
    # embed_mod`-style (verify actual import) -- patched above at module level, so
    # any `embed.get_embedding(...)` attribute lookup sees the wrapper.

    orig_solve = pl.LpSolverDefault.actualSolve

    def counted_solve(*a, **kw):
        counters.actual_solve_calls += 1
        return orig_solve(*a, **kw)

    pl.LpSolverDefault.actualSolve = counted_solve
    restores.append((pl.LpSolverDefault, "actualSolve", orig_solve))

    orig_gvm = chem_mod.Molecule.get_valid_molecule

    def counted_gvm(*a, **kw):
        counters.get_valid_molecule_calls += 1
        return orig_gvm(*a, **kw)

    chem_mod.Molecule.get_valid_molecule = counted_gvm
    restores.append((chem_mod.Molecule, "get_valid_molecule", orig_gvm))

    orig_cg = cg_mod.TMCOptimizer.clean_geometry

    def counted_cg(*a, **kw):
        counters.clean_geometry_calls += 1
        t0 = time.monotonic()
        result = orig_cg(*a, **kw)
        print(
            f"[progress] clean_geometry call#{counters.clean_geometry_calls} "
            f"took {time.monotonic() - t0:.2f}s",
            file=sys.stderr,
            flush=True,
        )
        return result

    cg_mod.TMCOptimizer.clean_geometry = counted_cg
    restores.append((cg_mod.TMCOptimizer, "clean_geometry", orig_cg))

    from oinsmiles.generation import metallogen_adapter as adapter_mod

    orig_build_contract = adapter_mod.build_contract_mol

    def counted_build_contract(*a, **kw):
        counters.build_contract_mol_calls += 1
        t0 = time.monotonic()
        result = orig_build_contract(*a, **kw)
        print(
            f"[progress] build_contract_mol call#{counters.build_contract_mol_calls} "
            f"took {time.monotonic() - t0:.2f}s",
            file=sys.stderr,
            flush=True,
        )
        return result

    adapter_mod.build_contract_mol = counted_build_contract
    restores.append((adapter_mod, "build_contract_mol", orig_build_contract))

    orig_key_matches = adapter_mod._reencode_key_matches

    def counted_key_matches(*a, **kw):
        counters.reencode_key_matches_calls += 1
        t0 = time.monotonic()
        result = orig_key_matches(*a, **kw)
        print(
            f"[progress] _reencode_key_matches call#{counters.reencode_key_matches_calls} "
            f"took {time.monotonic() - t0:.2f}s -> {result}",
            file=sys.stderr,
            flush=True,
        )
        return result

    adapter_mod._reencode_key_matches = counted_key_matches
    restores.append((adapter_mod, "_reencode_key_matches", orig_key_matches))

    orig_reencode_fast = adapter_mod._reencode_oin_fast

    def counted_reencode_fast(*a, **kw):
        counters.reencode_oin_fast_calls += 1
        t0 = time.monotonic()
        result = orig_reencode_fast(*a, **kw)
        print(
            f"[progress]   _reencode_oin_fast call#{counters.reencode_oin_fast_calls} "
            f"took {time.monotonic() - t0:.2f}s -> {'None' if result is None else 'oin'}",
            file=sys.stderr,
            flush=True,
        )
        return result

    adapter_mod._reencode_oin_fast = counted_reencode_fast
    restores.append((adapter_mod, "_reencode_oin_fast", orig_reencode_fast))

    orig_reencode_full = adapter_mod._reencode_oin

    def counted_reencode_full(*a, **kw):
        counters.reencode_oin_full_calls += 1
        t0 = time.monotonic()
        result = orig_reencode_full(*a, **kw)
        print(
            f"[progress]   _reencode_oin (FULL XYZToSMILES) call#{counters.reencode_oin_full_calls} "
            f"took {time.monotonic() - t0:.2f}s -> {'None' if result is None else 'oin'}",
            file=sys.stderr,
            flush=True,
        )
        return result

    adapter_mod._reencode_oin = counted_reencode_full
    restores.append((adapter_mod, "_reencode_oin", orig_reencode_full))

    return restores


def restore_all(restores):
    for obj, attr, orig in restores:
        setattr(obj, attr, orig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--molecule", required=True, help="e.g. QIDKUL_comp_0")
    ap.add_argument("--timeout", type=float, default=60.0, help="internal embed_time_budget cap")
    ap.add_argument("--optimizer", default=None)
    ap.add_argument("--ensemble-size", type=int, default=1)
    ap.add_argument(
        "--max-attempts",
        type=int,
        default=None,
        help=(
            "Deterministic cap on the outer pool-fill attempt loop "
            "(ff_params['max_attempts']). Contention-robust: bounds call COUNT, "
            "not wall clock, so it is the preferred triage knob on a busy host."
        ),
    )
    ap.add_argument("--out", default=None, help="optional path to append JSON result")
    args = ap.parse_args()

    from oinsmiles import XYZToSMILES
    from oinsmiles.generation.metallogen_adapter import (
        OIN3DGeneratorMetallogen as OIN3DGenerator,
    )
    from oinsmiles.generator3d.utils.compute_chg_and_bo_pulp import pulp_cache_stats

    xyz_path = find_xyz(args.dataset, args.molecule)
    n_atoms = atom_count(xyz_path)

    xyz_to_smiles = XYZToSMILES()
    t0 = time.monotonic()
    oin_string = xyz_to_smiles.convert(xyz_path)
    encode_s = time.monotonic() - t0
    print(f"[progress] encode_s={encode_s:.2f}", file=sys.stderr, flush=True)

    counters = Counters()
    restores = instrument(counters)
    ff_params = {"max_attempts": args.max_attempts} if args.max_attempts is not None else None
    gen = OIN3DGenerator(
        optimizer=args.optimizer,
        ensemble_size=args.ensemble_size,
        timeout=args.timeout,
        ff_params=ff_params,
    )
    error = None
    pool_size = None
    t0 = time.monotonic()
    try:
        result = gen.generate(oin_string)
        pool_size = 1 if result is not None else 0
    except Exception as e:
        error = f"{type(e).__name__}: {e}"
    gen_s = time.monotonic() - t0
    restore_all(restores)

    cache = pulp_cache_stats()
    row = {
        "molecule": args.molecule,
        "n_atoms": n_atoms,
        "timeout_budget_s": args.timeout,
        "max_attempts": args.max_attempts,
        "optimizer": args.optimizer,
        "ensemble_size": args.ensemble_size,
        "encode_s": round(encode_s, 3),
        "generate_s": round(gen_s, 3),
        "error": error,
        "pool_size_returned": pool_size,
        **counters.as_dict(),
        "pulp_cache_hits": cache["hits"],
        "pulp_cache_misses": cache["misses"],
    }
    print(json.dumps(row, indent=2))
    if args.out:
        with open(args.out, "a") as f:
            f.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
