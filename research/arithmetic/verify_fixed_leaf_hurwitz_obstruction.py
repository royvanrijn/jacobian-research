#!/usr/bin/env python3
"""Exact controls for the uniform fixed-leaf Hurwitz strategy obstruction."""
from __future__ import annotations
import argparse
from collections import deque
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
HELPER = HERE / "verify_paired_hurwitz_leaf_extraction.py"
spec = importlib.util.spec_from_file_location("leaf_gates", HELPER)
leaf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(leaf)
m = leaf.m

def owner(t, x):
    moving = [i for i, g in enumerate(t) if g[x] != x]
    return moving[0] if len(moving) == 1 else None

def deletion(t, x):
    at = owner(t, x)
    assert at is not None, "Deletion requires the selected sheet to be a leaf"
    n = len(t[0])
    aux = m.perm(n, [(x, t[at][x])])
    without_aux = m.compose(aux, t[at])
    assert without_aux[x] == x
    assert m.cycle_lengths(without_aux) == [1] * (n - 2) + [2]
    old = [a for a in range(n) if a != x]
    lookup = {a: i for i, a in enumerate(old)}
    small = []
    for i, g in enumerate(t):
        g = without_aux if i == at else g
        assert g[x] == x
        small.append(tuple(lookup[g[a]] for a in old))
    small = tuple(small)
    P = m.product(t)
    Q = m.compose(m.perm(n, [(x, P[x])]), P)
    expected = tuple(lookup[Q[a]] for a in old)
    assert m.product(small) == expected
    assert m.cycle_lengths(expected) == [n - 1]
    return small

def removable_owner(t, x):
    at = owner(t, x)
    if at is None:
        return False
    degree = leaf.degrees(t)
    return all(degree[a] == 1 or degree[b] == 1 for a, b in leaf.edges(t[at]))

def family(r):
    n, x, u = 2 * r + 1, 2 * r, 2 * r - 2
    A = m.perm(n, [(x, u), (0, 1)])
    B = tuple(m.perm(n, [(2 * i - 2, 2 * i), (2 * i - 1, 2 * i + 1)])
              for i in range(1, r))
    return (A,) + B

def all_degree7_tuples():
    D = set(leaf.dt_class(7))
    P = m.perm(7, [tuple(range(7))])
    tuples = []
    for a in D:
        for b in D:
            c = m.compose(m.inverse(m.compose(a, b)), P)
            if c in D:
                tuples.append((a, b, c))
    assert len(tuples) == 392
    return tuples

def complete_equivariance():
    allowed = rejected = chosen_leaves = 0
    for t in all_degree7_tuples():
        for x in range(7):
            if owner(t, x) is None:
                continue
            chosen_leaves += 1
            small = deletion(t, x)
            for i in range(2):
                for sign in (-1, 1):
                    out = m.hurwitz(t, i, sign)
                    if owner(out, x) is None:
                        rejected += 1
                        continue
                    assert deletion(out, x) == m.hurwitz(small, i, sign)
                    allowed += 1
    assert allowed and rejected
    return {"complete_fixed_product_domain": 392,
            "all_chosen_leaf_inputs": chosen_leaves,
            "all_leaf_preserving_moves_checked": allowed,
            "moves_spreading_the_chosen_sheet": rejected}

def family_control(r):
    t = family(r)
    n, x, u = 2 * r + 1, 2 * r, 2 * r - 2
    assert m.cycle_lengths(m.product(t)) == [n]
    degrees = leaf.degrees(t)
    assert [a for a in range(n) if degrees[a] == 1] == [2 * r - 1, x]
    assert all(not all(degrees[a] == 1 or degrees[b] == 1
                       for a, b in leaf.edges(g)) for g in t)
    assert m.power(m.compose(t[0], t[-1]), 2) == m.perm(n, [(x, u - 2, u)])
    lower = deletion(t, x)
    e = m.perm(2 * r, [(0, 1)])
    assert lower == (e,) + tuple(g[:-1] for g in t[1:])
    blocks = [frozenset((2 * j, 2 * j + 1)) for j in range(r)]
    for i, g in enumerate(lower[1:], 1):
        expected = list(blocks)
        expected[i - 1], expected[i] = expected[i], expected[i - 1]
        assert [frozenset(g[a] for a in block) for block in blocks] == expected
    sigma = tuple(range(2 * r))
    flips = [e]
    for i, g in enumerate(lower[1:], 1):
        sigma = m.compose(g, sigma)
        flip = m.conjugate(sigma, e)
        assert flip == m.perm(2 * r, [(2 * i, 2 * i + 1)])
        flips.append(flip)
    assert len(set(flips)) == r
    word = [-i for i in range(1, r - 1)] + [-(r - 1), -(r - 1)]
    assert len(word) == r
    out = t
    owner_counts = []
    for step, code in enumerate(word):
        before = out
        out = m.hurwitz(out, abs(code) - 1, -1)
        owner_counts.append(sum(g[x] != x for g in out))
        if owner(out, x) is not None and owner(before, x) is not None:
            assert deletion(out, x) == m.hurwitz(deletion(before, x), abs(code) - 1, -1)
    assert owner_counts == [1] * (r - 2) + [2, 1]
    assert removable_owner(out, x)
    assert out[-1] == m.perm(n, [(x, u - 2), (2 * r - 3, 2 * r - 1)])
    assert m.product(out) == m.product(t)
    lower_after = deletion(out, x)
    small_group_before = small_group_after = None
    if r <= 4:
        small_group_before = leaf.group_closure(lower)
        assert small_group_before == 2 ** r * factorial(r)
    if r == 3:
        assert leaf.group_closure(t) == factorial(n) // 2
        small_group_after = leaf.group_closure(lower_after)
        assert small_group_after == factorial(2 * r)
    return {"r": r, "degree": n,
            "full_group_scope": "A_N by written primitivity/3-cycle proof",
            "lower_wreath_order_formula": 2 ** r * factorial(r),
            "independent_lower_group_closure_before": small_group_before,
            "independent_lower_group_closure_after": small_group_after,
            "no_initial_removable_DT": True,
            "escape_word_length": len(word),
            "escape_sheet_owner_counts": owner_counts if r <= 6 else "1 repeated r-2 times, then2,1",
            "two_exclusive_leaf_ends_after_escape": True}

def handoff_identity():
    n = 7
    # x=0,u=1,v=2,a=3,b=4,c=5,d=6
    pair = (m.perm(n, [(0, 1), (3, 4)]), m.perm(n, [(1, 2), (5, 6)]))
    middle = m.hurwitz(pair, 0, -1)
    final = m.hurwitz(middle, 0, -1)
    assert final == (m.perm(n, [(1, 2), (3, 4)]),
                     m.perm(n, [(0, 2), (5, 6)]))
    assert owner(pair, 0) is not None and owner(middle, 0) is None and owner(final, 0) is not None
    return {"literal_Q_inverse_squared_identity": True,
            "chosen_sheet_occurs_in_two_intermediate_factors": True}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path,
        default=HERE / "fixed_leaf_hurwitz_obstruction_verification.json")
    args = ap.parse_args()
    equivariance = complete_equivariance()
    controls = [family_control(r) for r in (3, 4, 5, 6, 12, 25, 50)]
    handoff = handoff_identity()
    inputs = ["FIXED_LEAF_HURWITZ_OBSTRUCTION.md",
              "PAIRED_HURWITZ_LEAF_EXTRACTION.md",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py",
              "alternating_hurwitz_source_ledger.json"]
    report = {
        "schema_version": 1, "kind": "fixed_leaf_strategy_obstruction",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "complete_degree7_deletion_equivariance": equivariance,
        "uniform_family_controls": controls,
        "generic_leaf_handoff_identity": handoff,
        "scope": "Written all-degree obstruction to retaining a chosen leaf throughout elementary moves, plus full-alternating witness family and explicit unrestricted escape. No unrestricted orbit distinction or arbitrary-tree leaf-extraction theorem.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()

