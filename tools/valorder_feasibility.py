#!/usr/bin/env python
"""Does a valid Lewis structure EXIST for an over-cap ligand -- and where is it in each order?

The 20 000-candidate over-cap search is a prefix of a space with 10^6 to 10^60 members, so
"we tried 20 000 and found none" cannot distinguish *the ligand has no Lewis structure* from
*we looked at the wrong 20 000*. Neither can a bigger budget: you cannot exhaust 10^60.

But ``AC2BO``'s own acceptance predicate has a **necessary condition that is additive over
atoms**, which makes the question exactly decidable:

1. A valid return implies ``BO.sum(axis=1) == valences`` exactly. (Either ``UA`` is empty and
   ``BO = AC``, or ``valences_not_too_large`` bounds every atom above by ``valences`` while
   ``(BO - AC).sum() == sum(DU)`` forces the total, so every per-atom slack is zero.)
2. Therefore ``charge_is_OK`` evaluates ``get_atomic_charge`` on the *candidate* valences, a
   pure per-atom function. Write ``Q0(valences) = sum_i get_atomic_charge(z_i, v_i)``.
3. The only correction ``charge_is_OK`` can apply is ``Q += 2`` per trivalent single-bonded
   carbon, and only while running *below* the target. So ``Q_final = Q0 + 2k, k >= 0``, giving
   **C1: ``Q0 <= charge`` and ``charge - Q0`` even**.
4. Every added bond raises ``BO.sum()`` by 2, so **C2: ``sum(valences) - sum(AC_valence)``
   must be even**.

C1 and C2 are additive, so a suffix DP over ``(Q0, parity)`` counts the survivors in the whole
space exactly -- and because both enumeration orders are *lexicographic in a known variable
order* (raw = atom order; heuristic = O, N, C, P, S, then the rest), the same DP yields the
exact **rank of the first survivor in each order** without iterating to it.

That is the measurement the ordering hypothesis needs:

* **0 survivors** -> no valid Lewis structure exists at this charge/carbene setting. The
  question is closed permanently; no ordering and no budget can help.
* survivors, first one at rank < 20 000 in heuristic order but >> 20 000 in raw order ->
  ``found_valid = 0`` is an ORDERING artefact.
* survivors, both ranks >> 20 000 -> ordering is not the lever either.

C1/C2 are **necessary, not sufficient** (the matching in ``get_BO`` may still fail to
saturate), so the survivors this tool finds are then handed to the real predicate --
``BO_is_OK`` on the actual ``get_BO`` output -- and that verdict is reported separately.

    $V tools/valorder_feasibility.py --cache-dir <dir> --mols HICLAG_comp_0 --verify 40
"""

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import numpy as np  # noqa: E402

from oinsmiles.utils.perception_core import (  # noqa: E402
    _HEURISTIC_ELEMENTS,
    BO_is_OK,
    atomic_valence_electrons,
    get_atomic_charge,
    get_BO,
    get_UA,
    get_UA_pairs,
    possible_valences,
)


def load_inputs(cache):
    data = np.load(cache, allow_pickle=False)
    return {
        "AC": data["AC"],
        "atoms": [int(x) for x in data["atoms"]],
        "charge": int(data["charge"]),
        "allow_carbenes": bool(data["allow_carbenes"]),
    }


def variable_order(atoms, vll, arm):
    """The (index, choices) sequence each enumeration order varies from slowest to fastest.

    Both orders are lexicographic: the FIRST variable is the slowest-varying. For ``raw``
    that is ``itertools.product``'s own convention (atom 0 slowest, atom n-1 fastest). For
    ``ordered`` it is ``iter_ordered_valences``: the O group, then N, C, P, S, then the
    remaining atoms -- and within each group, atom order, with the non-O/N/C/P/S lists
    ascending (which is what ``_ordered_valences``' ``sorted()`` tie-break does).
    """
    if arm == "raw":
        return [(i, list(vll[i])) for i in range(len(atoms))]
    groups = [[i for i, num in enumerate(atoms) if num == el] for el in _HEURISTIC_ELEMENTS]
    grouped = {i for group in groups for i in group}
    order = []
    for group in groups:
        order += [(i, list(vll[i])) for i in group]
    order += [(i, sorted(vll[i])) for i in range(len(atoms)) if i not in grouped]
    return order


def suffix_counts(order, atoms):
    """``counts[k]`` maps ``(dq, dpar)`` contributed by variables ``k..n-1`` to a count."""
    n = len(order)
    counts = [None] * (n + 1)
    counts[n] = {(0, 0): 1}
    for k in range(n - 1, -1, -1):
        idx, choices = order[k]
        z = atoms[idx]
        table = {}
        for val in choices:
            dq = get_atomic_charge(z, atomic_valence_electrons[z], val)
            dpar = val & 1
            for (q, par), c in counts[k + 1].items():
                key = (q + dq, par ^ dpar)
                table[key] = table.get(key, 0) + c
        counts[k] = table
    return counts


def targets(counts0, charge, ac_parity, exact_only):
    """The ``(dq, dpar)`` states that satisfy C1 (or C1 with equality) and C2."""
    want_par = ac_parity & 1
    out = set()
    for q, par in counts0:
        if par != want_par:
            continue
        if q == charge or (not exact_only and q <= charge and (charge - q) % 2 == 0):
            out.add((q, par))
    return out


def count_survivors(counts, target_set):
    return sum(counts[0].get(t, 0) for t in target_set)


def tail_sizes(order):
    """``tail[k]`` = how many candidates share any fixed prefix of length ``k+1``.

    The mixed-radix weight of position ``k``: the enumeration is lexicographic with
    ``order[0]`` slowest, so advancing position ``k`` by one step skips ``tail[k]``
    candidates.
    """
    n = len(order)
    tail = [1] * (n + 1)
    for k in range(n - 1, -1, -1):
        tail[k] = tail[k + 1] * len(order[k][1])
    return tail[1:] + [1]


def first_survivor(order, atoms, counts, target_set):
    """The lexicographically-first survivor in this order.

    Greedy: at each position take the smallest value that still admits a completion
    reaching the target. Feasibility comes from the suffix DP, so no search is needed.
    """
    if count_survivors(counts, target_set) == 0:
        return None
    acc_q, acc_par = 0, 0
    chosen = {}
    for k, (idx, choices) in enumerate(order):
        z = atoms[idx]
        for val in choices:
            q2 = acc_q + get_atomic_charge(z, atomic_valence_electrons[z], val)
            par2 = acc_par ^ (val & 1)
            if any(counts[k + 1].get((tq - q2, tp ^ par2), 0) for tq, tp in target_set):
                chosen[idx] = val
                acc_q, acc_par = q2, par2
                break
        else:  # pragma: no cover - guarded by count_survivors above
            raise AssertionError("no feasible choice despite a non-zero survivor count")
    return tuple(chosen[i] for i in range(len(atoms)))


def rank_of(order, tail, candidate):
    """Exact 0-based position of ``candidate`` in this order's full enumeration.

    This is the number that decides the hypothesis: ``AC2BO`` sees a *prefix* of this
    enumeration, so a candidate at rank >= ``_VALENCE_FALLBACK_TRIES`` is unreachable no
    matter how sensible it is. Mixed-radix, so exact for spaces far beyond 2**64.
    """
    rank = 0
    for k, (idx, choices) in enumerate(order):
        rank += choices.index(candidate[idx]) * tail[k]
    return rank


def verify(candidate, inputs):
    """Run ``AC2BO``'s real acceptance predicate on one candidate."""
    AC, atoms, charge = inputs["AC"], inputs["atoms"], inputs["charge"]
    AC_valence = list(AC.sum(axis=1))
    UA, DU = get_UA(list(candidate), AC_valence)
    if not UA:
        ok = BO_is_OK(
            AC,
            AC,
            charge,
            DU,
            atomic_valence_electrons,
            atoms,
            list(candidate),
            allow_charged_fragments=True,
            allow_carbenes=inputs["allow_carbenes"],
        )
        return {"valid": bool(ok), "path": "UA-empty", "bo_sum": int(AC.sum())}
    for UA_pairs in get_UA_pairs(UA, AC, DU, use_graph=True):
        BO = get_BO(AC, UA, DU, list(candidate), UA_pairs, use_graph=True)
        ok = BO_is_OK(
            BO,
            AC,
            charge,
            DU,
            atomic_valence_electrons,
            atoms,
            list(candidate),
            allow_charged_fragments=True,
            allow_carbenes=inputs["allow_carbenes"],
        )
        if ok:
            return {"valid": True, "path": "matched", "bo_sum": int(BO.sum())}
    return {"valid": False, "path": "matched", "bo_sum": int(BO.sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", required=True, help="npz files from valorder_probe.py")
    ap.add_argument("--mols", required=True)
    ap.add_argument(
        "--verify",
        type=int,
        default=0,
        help="cap on how many survivors to hand to the real predicate, per arm",
    )
    ap.add_argument("--verify-arms", default="", help="e.g. raw,ordered")
    ap.add_argument(
        "--max-rank",
        type=int,
        default=20_000,
        help="stop verifying once a survivor's rank leaves the search's reach "
        "(default = _VALENCE_FALLBACK_TRIES)",
    )
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    names = [x.strip() for x in args.mols.split(",") if x.strip()]
    out = {}
    for name in names:
        cache = Path(args.cache_dir) / f"{name}.npz"
        if not cache.exists():
            print(f"{name}: no cache at {cache} -- run valorder_probe.py first", flush=True)
            continue
        inputs = load_inputs(str(cache))
        atoms, AC, charge = inputs["atoms"], inputs["AC"], inputs["charge"]
        vll = possible_valences(
            list(AC.sum(axis=1)), atoms, allow_carbenes=inputs["allow_carbenes"]
        )
        space = 1
        for lst in vll:
            space *= len(lst)
        ac_parity = int(AC.sum(axis=1).sum()) & 1

        rec = {
            "ligand_atoms": len(atoms),
            "charge": charge,
            "allow_carbenes": inputs["allow_carbenes"],
            "space": space,
            "arms": {},
        }
        # Which atoms actually vary, and by which element. If every multi-choice atom
        # belongs to ONE of _HEURISTIC_ELEMENTS' groups, the heuristic order and the raw
        # product order are the *same sequence*, because single-choice atoms contribute
        # nothing and a group is traversed in atom order. Then the ordering hypothesis is
        # vacuous for this ligand, and saying so is the whole point of printing it.
        varying = {}
        for i, lst in enumerate(vll):
            if len(lst) > 1:
                varying[atoms[i]] = varying.get(atoms[i], 0) + 1
        raw_seq = [(i, tuple(ch)) for i, ch in variable_order(atoms, vll, "raw") if len(ch) > 1]
        ord_seq = [(i, tuple(ch)) for i, ch in variable_order(atoms, vll, "ordered") if len(ch) > 1]
        rec["varying_by_element"] = {str(k): v for k, v in sorted(varying.items())}
        rec["orders_identical"] = raw_seq == ord_seq
        print(
            f"\n=== {name}  {len(atoms)} atoms  charge={charge}  "
            f"carbenes={inputs['allow_carbenes']}  space=10^{len(str(space)) - 1} ({space})",
            flush=True,
        )
        print(
            f"  varying atoms by element: "
            f"{ {k: v for k, v in sorted(varying.items())} }   "
            f"raw and heuristic orders identical: {raw_seq == ord_seq}",
            flush=True,
        )

        for arm in ("raw", "ordered"):
            order = variable_order(atoms, vll, arm)
            tail = tail_sizes(order)
            counts = suffix_counts(order, atoms)
            t_exact = targets(counts[0], charge, ac_parity, exact_only=True)
            t_loose = targets(counts[0], charge, ac_parity, exact_only=False)
            n_exact = count_survivors(counts, t_exact)
            n_loose = count_survivors(counts, t_loose)
            cand = first_survivor(order, atoms, counts, t_exact)
            row = {
                "survivors_exact": str(n_exact),
                "survivors_loose": str(n_loose),
                "first_candidate": list(cand) if cand else None,
            }
            if cand is not None:
                row["first_rank"] = str(rank_of(order, tail, cand))
                row["verify_first"] = verify(cand, inputs)
            rec["arms"][arm] = row
            print(
                f"  {arm:<8s} survivors(Q0==charge) = {n_exact}"
                f"  (10^{len(str(n_exact)) - 1} of 10^{len(str(space)) - 1})"
                f"  survivors(C1) = {n_loose}",
                flush=True,
            )
            if cand is None:
                print("           NO candidate can satisfy the charge condition", flush=True)
            else:
                v = row["verify_first"]
                print(
                    f"           first survivor at rank {row['first_rank']}"
                    f"  real predicate: {'VALID' if v['valid'] else 'rejected'}"
                    f" [{v['path']}, BO.sum={v['bo_sum']}]",
                    flush=True,
                )

        for arm in [a for a in args.verify_arms.split(",") if a.strip()]:
            # Walk the survivors in this order's own sequence and run the REAL predicate
            # until one passes or the rank leaves the search's reach. "First valid rank"
            # is the number that decides whether the shipped budget could ever find it.
            order = variable_order(atoms, vll, arm)
            tail = tail_sizes(order)
            counts = suffix_counts(order, atoms)
            t_exact = targets(counts[0], charge, ac_parity, exact_only=True)
            checked, first_valid, last_rank = 0, None, None
            for cand in iter_survivors(order, atoms, counts, t_exact, limit=args.verify):
                rank = rank_of(order, tail, cand)
                if rank >= args.max_rank:
                    break
                last_rank = rank
                checked += 1
                if verify(cand, inputs)["valid"]:
                    first_valid = rank
                    break
            rec["arms"][arm]["verify_sweep"] = {
                "checked": checked,
                "max_rank": args.max_rank,
                "first_valid_rank": first_valid,
                "last_rank_checked": last_rank,
            }
            print(
                f"  {arm:<8s} verified {checked} survivors with rank < {args.max_rank}: "
                + (
                    f"first VALID at rank {first_valid}"
                    if first_valid is not None
                    else "none valid in reach"
                ),
                flush=True,
            )
        out[name] = rec

    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1, default=str))
        print(f"\nwrote {args.out}", flush=True)


def iter_survivors(order, atoms, counts, target_set, limit):
    """Yield survivors in this order's sequence, skipping infeasible subtrees via the DP."""
    n = len(order)
    yielded = 0
    assign = {}

    def rec(k, acc_q, acc_par):
        nonlocal yielded
        if yielded >= limit:
            return
        if k == n:
            yield tuple(assign[i] for i in range(len(atoms)))
            yielded += 1
            return
        idx, choices = order[k]
        z = atoms[idx]
        for val in choices:
            q2 = acc_q + get_atomic_charge(z, atomic_valence_electrons[z], val)
            par2 = acc_par ^ (val & 1)
            if not any(counts[k + 1].get((tq - q2, tp ^ par2), 0) for tq, tp in target_set):
                continue
            assign[idx] = val
            yield from rec(k + 1, q2, par2)
            if yielded >= limit:
                return

    yield from rec(0, 0, 0)


if __name__ == "__main__":
    main()
