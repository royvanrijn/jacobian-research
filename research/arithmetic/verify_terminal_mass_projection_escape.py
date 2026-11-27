#!/usr/bin/env python3
"""Exact controls of fixed-projection escape; the universal proof is written."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import platform

import verify_terminal_mass_hurwitz_descent as v

if not __debug__:
    raise RuntimeError("projection checker requires active assertions")

HERE = Path(__file__).resolve().parent
m, h = v.m, v.h


def trees(t):
    """All labeled trees, via the bijective Pruefer encoding."""
    for code in product(range(t), repeat=t - 2):
        degree = Counter(code)
        edges = []
        for x in code:
            leaf = next(y for y in range(t) if degree[y] == 0)
            edges.append((leaf, x))
            degree[leaf] = -1
            degree[x] -= 1
        leaves = [y for y in range(t) if degree[y] == 0]
        edges.append(tuple(leaves))
        yield tuple(tuple(sorted(e)) for e in edges)


def tm(t, u, w, x=None, y=None):
    return m.perm(3 * t, [(u if x is None else x, w if y is None else y),
                         (2 * t + u, 2 * t + w)])


def projected(g, t):
    edges = h.edges(g)
    middle = [e for e in edges if all(x < 2 * t for x in e)]
    terminal = [e for e in edges if all(2 * t <= x < 3 * t for x in e)]
    assert len(middle) == len(terminal) == 1
    return {x % t for x in middle[0]} == {x - 2 * t for x in terminal[0]}


def conjugated_by(B, C):
    # Independent endpoint substitution, checked against permutation conjugation.
    result = m.perm(len(C), [tuple(B[x] for x in e) for e in h.edges(C)])
    assert result == m.conjugate(B, C)
    return result


def lifted_tree_controls():
    rows = []
    for t in (2, 3, 4):
        tested = split = consistent = 0
        for es in trees(t):
            for choices in product((0, 1), repeat=2 * len(es)):
                factors = [tm(t, u, w, u + t * choices[2 * j],
                              w + t * choices[2 * j + 1])
                           for j, (u, w) in enumerate(es)]
                assert all(projected(C, t) for C in factors)
                lifts = {u: set() for u in range(t)}
                for C in factors:
                    for e in h.edges(C):
                        if e[0] < 2 * t:
                            for x in e:
                                lifts[x % t].add(x)
                tested += 1
                bad = next((u for u, xs in lifts.items() if len(xs) == 2), None)
                failures = [(i, j) for i, C in enumerate(factors)
                            for j, D in enumerate(factors) if i != j
                            and not projected(conjugated_by(C, D), t)]
                assert bool(failures) == (bad is not None)
                if bad is None:
                    consistent += 1
                    assert len(set().union(*lifts.values())) == t
                else:
                    split += 1
        rows.append({"terminal_size": t, "lift_assignments": tested,
                     "consistent": consistent, "split_detected": split})
    return rows


def mm_controls():
    rows = []
    for t in (3, 4, 5):
        tested = fibrewise = touching = 0
        ds = [m.perm(3 * t, h.edges(g)) for g in h.dt_class(2 * t)]
        for es in trees(t):
            factors = [tm(t, u, w) for u, w in es]
            for B in ds:
                tested += 1
                edge_fixed = all(projected(conjugated_by(B, C), t) for C in factors)
                expected = all(B[x] % t == x for x in range(t))
                assert edge_fixed == expected
                if not edge_fixed:
                    continue
                fibrewise += 1
                moved = next((x for x in range(t) if B[x] != x), None)
                if moved is None:
                    continue
                touching += 1
                x = moved
                C = next(C for C in factors if C[x] != x)
                y = C[x]
                Bnew = conjugated_by(C, B)
                assert Bnew[y] == x + t and Bnew[y] % t != y
                witness = next(D for D in factors
                               if not projected(conjugated_by(Bnew, D), t))
                assert projected(witness, t)
        rows.append({"terminal_size": t, "tree_MM_pairs": tested,
                     "edgewise_fibre_preserving": fibrewise,
                     "active_attachments_broken_by_conjugation": touching})
    # The two-vertex tree cannot force alpha to fix individual vertices.
    C = tm(2, 0, 1)
    B = m.perm(6, [(0, 3), (1, 2)])
    assert projected(conjugated_by(B, C), 2)
    assert B[0] % 2 != 0
    return {"all_labeled_tree_controls": rows,
            "two_vertex_edge_exchange_is_a_real_exception": True}


def coordinates(R, root, size):
    answer = {root: 0}
    x = R[root]
    for k in range(1, size):
        assert x not in answer
        answer[x] = k
        x = R[x]
    assert x == root
    return answer


def path_control(t):
    n = 3 * t + 1
    terminal_middle = [m.perm(n, [(i, i + 1), (2 * t + i, 2 * t + i + 1)])
                       for i in range(t - 1)]
    middle_middle = [m.perm(n, [(t - 1 + i, t + i),
                               (t - 1 + i + t // 2, t + i + t // 2)])
                     for i in range(t // 2)]
    rest = tuple(terminal_middle + middle_middle)
    R = m.product(rest)
    mc = coordinates(R, 0, 2 * t)
    uc = coordinates(R, 2 * t, t)
    ui = {k: x for x, k in uc.items()}
    pi = {x: ui[k % t] for x, k in mc.items()}
    c = next(x for x, k in mc.items() if k == t)
    A = m.perm(n, [(3 * t, 0), (c, 2 * t)])
    initial = (A,) + rest
    assert m.cycle_lengths(m.product(initial)) == [n]
    assert v.factor_cut(initial, 0)[0] == v.mu(initial) == t - 1

    def bad_factor(state):
        for j, g in enumerate(state[1:], 1):
            es = h.edges(g)
            ue = [e for e in es if all(x in uc for x in e)]
            if ue:
                me, = [e for e in es if all(x in mc for x in e)]
                if {pi[x] for x in me} != set(ue[0]):
                    return j
        return None

    assert bad_factor(initial) is None
    word = [t]  # Conjugate the final TM by the first MM, keeping A fixed.
    escaped = m.braid(initial, word)
    assert escaped[0] == A and m.product(escaped[1:]) == R
    bad = bad_factor(escaped)
    assert bad is not None
    ue, = [e for e in h.edges(escaped[bad]) if all(x in uc for x in e)]
    me, = [e for e in h.edges(escaped[bad]) if all(x in mc for x in e)]
    missing = min(set(ue) - {pi[x] for x in me})
    phase = uc[missing]
    full_twist = list(range(2, len(initial))) * (len(initial) - 1)
    twisted = m.braid(escaped, full_twist * phase)
    power = m.power(R, phase)
    expected = (A,) + tuple(m.conjugate(m.inverse(power), B) for B in escaped[1:])
    assert twisted == expected
    common = {x for x in range(n) if A[x] != x and twisted[bad][x] != x}
    assert common == {2 * t}
    # Literal bringing and the final pruning move of Lemma 3.
    bringing = v.bring_B_word(0, bad)
    gathered = m.braid(twisted, bringing)
    assert gathered[0] == A and gathered[1] == twisted[bad]
    pruned = m.braid(gathered, [-1])
    assert v.factor_cut(pruned, 1)[0] < t - 1
    assert m.product(pruned) == m.product(initial)
    return {"terminal_size": t, "degree": n,
            "initial_mass": t - 1, "projection_break_word": word,
            "mass_after_projection_break": v.mu(escaped),
            "bad_factor_index_zero_based": bad, "outer_contact_phase": phase,
            "twist_word": full_twist, "twist_repetitions": phase,
            "bringing_word": bringing, "final_prune_word": [-1],
            "final_chosen_factor_mass": v.factor_cut(pruned, 1)[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "terminal_mass_projection_escape_verification.json")
    args = parser.parse_args()
    lifts = lifted_tree_controls()
    middle = mm_controls()
    paths = [path_control(t) for t in (4, 6, 8, 12, 16)]
    inputs = ["TERMINAL_MASS_HURWITZ_DESCENT.md",
              "verify_terminal_mass_hurwitz_descent.py",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt = {
        "schema_version": 1, "kind": "terminal_mass_projection_escape_controls",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
                         for p in inputs},
        "lift_consistency_controls": lifts, "MM_conjugation_controls": middle,
        "literal_saturated_tree_escapes": paths,
        "scope": "Bounded exact algebra controls and five deterministic full-tree "
                 "escape certificates. No new degree census. Lemma 19 and Theorem 20 "
                 "are universal written proofs, not finite extrapolations. Pure-DT "
                 "connectivity additionally uses the retained degree-thirteen base "
                 "and two-leaf absorption. Mixed connectivity remains separate.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "lift_cases": sum(r["lift_assignments"] for r in lifts),
                      "MM_cases": sum(r["tree_MM_pairs"] for r in middle["all_labeled_tree_controls"]),
                      "literal_escapes": len(paths), "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
