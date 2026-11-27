#!/usr/bin/env python3
"""Exact bounded controls and complete small bases for mixed Hurwitz isolation."""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
from itertools import combinations, permutations, product
import json
from math import comb, factorial
from pathlib import Path
import platform

import verify_terminal_mass_hurwitz_descent as v
import verify_terminal_mass_projection_escape as projection

if not __debug__:
    raise RuntimeError("mixed checker requires active assertions")

HERE = Path(__file__).resolve().parent
m, h = v.m, v.h


def support(g):
    return {x for x, y in enumerate(g) if x != y}


def mixed_cut(state):
    n = len(state[0])
    assert n == 2 * len(state) + 1
    indices = [i for i, g in enumerate(state) if m.cycle_lengths(g) == [1] * (n - 3) + [3]]
    assert len(indices) == 1
    at, = indices
    assert all(i == at or m.cycle_lengths(g) == [1] * (n - 4) + [2, 2]
               for i, g in enumerate(state))
    es = [e for g in state for e in h.edges(g)]
    assert len(es) == n - 1 and len(set(map(frozenset, es))) == n - 1
    assert len(v.components(n, es)) == 1
    rest = [e for i, g in enumerate(state) if i != at for e in h.edges(g)]
    cs = v.components(n, rest)
    assert len(cs) == 3 and all(len(c & support(state[at])) == 1 for c in cs)
    return at, cs, n - max(map(len, cs)) - 2


def pair_controls():
    T = m.perm(6, [(0, 1, 2)])
    B = m.perm(6, [(0, 3), (4, 5)])
    plus = m.braid((T, B), [1, 1])
    minus = m.braid((T, B), [-1, -1])
    assert plus == (m.perm(6, [(3, 1, 2)]), m.perm(6, [(0, 2), (4, 5)]))
    assert minus == (m.perm(6, [(0, 3, 2)]), m.perm(6, [(1, 3), (4, 5)]))
    owners_tested = 0
    for old in range(4):
        for leaf in (4, 5):
            raw = list(range(6)); raw[leaf] = old
            values = list(dict.fromkeys(raw))
            owner = [values.index(x) for x in raw]
            def cuts(pair):
                D = next(g for g in pair if 3 not in m.cycle_lengths(g))
                return v.components(5, [(owner[a], owner[b]) for a, b in h.edges(D)])
            initial = cuts((T, B))
            bside = next(c for c in initial if owner[1] in c)
            cside = next(c for c in initial if owner[2] in c)
            assert any(cside < c for c in cuts(plus))
            assert any(bside < c for c in cuts(minus))
            owners_tested += 1
    pair = (m.perm(5, [(0, 1, 2)]), m.perm(5, [(0, 3), (1, 4)]))
    assert m.braid(pair, [1, 1]) == (
        m.perm(5, [(2, 3, 4)]), m.perm(5, [(0, 2), (1, 3)]))
    # The three-component projection obstruction, on a literal pair only.
    B = m.perm(9, [(0, 1), (3, 4)])
    C = m.perm(9, [(0, 2), (6, 8)])
    new = m.conjugate(B, C)
    assert {x % 3 for x in h.edges(new)[0]} != {x % 3 for x in h.edges(new)[1]}
    return {"one_contact_owner_patterns": owners_tested,
            "both_two_move_identities": True, "two_leaf_capture_identity": True,
            "three_equal_projection_break": True}


def star_word(edges, root):
    edges = [tuple(sorted(e)) for e in edges]
    word = []
    def move(j, sign):
        A, B = edges[j:j + 2]
        swap = {a: b for a, b in ((B if sign > 0 else A),)}
        swap.update({b: a for a, b in ((B if sign > 0 else A),)})
        target = A if sign > 0 else B
        changed = tuple(sorted(swap.get(a, a) for a in target))
        edges[j:j + 2] = [B, changed] if sign > 0 else [changed, A]
        word.append(sign * (j + 1))
    def components(vertices, es):
        remaining = set(vertices); answer = []
        while remaining:
            stack = [remaining.pop()]; found = set(stack)
            for x in stack:
                for a, b in es:
                    y = b if a == x else a if b == x else None
                    if y in remaining:
                        remaining.remove(y); stack.append(y); found.add(y)
            answer.append(found)
        return answer
    def rec(start, length, center):
        if length == 0:
            return
        at = next(j for j in range(start, start + length) if center in edges[j])
        while at > start:
            move(at - 1, 1); at -= 1
        e = edges[start]; other = next(x for x in e if x != center)
        vertices = set().union(*(set(q) for q in edges[start:start + length]))
        cs = components(vertices, edges[start + 1:start + length])
        assert len(cs) == 2
        B = next(c for c in cs if other in c); A = vertices - B
        count = len(B) - 1
        for wanted in range(start + 1, start + 1 + count):
            at = next(j for j in range(wanted, start + length) if set(edges[j]) <= B)
            while at > wanted:
                assert set(edges[at]).isdisjoint(edges[at - 1])
                move(at - 1, 1); at -= 1
        rec(start + 1, count, other)
        rec(start + 1 + count, len(A) - 1, center)
        for j in range(start, start + count):
            assert edges[j] == e
            move(j, -1)
        assert all(center in q for q in edges[start:start + length])
    rec(0, len(edges), root)
    return tuple(word), tuple(edges)


def paired_tree(t, edges, orientation=1):
    T = m.perm(2 * t + 1, [(0, 1, t + 1)] if orientation > 0 else [(0, t + 1, 1)])
    return (T,) + tuple(m.perm(2 * t + 1, [(a + 1, b + 1), (a + t + 1, b + t + 1)])
                        for a, b in edges)


def isolation_control(t, edges, orientation=1):
    initial = paired_tree(t, edges, orientation)
    word, star = star_word(edges, 0)
    lifted = [1 + x if x > 0 else x - 1 for x in word]
    normalized = m.braid(initial, lifted)
    assert all(1 in support(D) and t + 1 in support(D) for D in normalized[1:])
    assert normalized[0] == initial[0]
    escaped = m.braid(normalized, [-1, -1, -1, 2, 2])
    at, cs, rho = mixed_cut(escaped)
    assert rho == 0 and sorted(map(len, cs)) == [1, 1, 2 * t - 1]
    leaves = set().union(*(set(c) for c in cs if len(c) == 1))
    assert not leaves & {1, t + 1} and len(leaves) == 2
    smaller = h.delete(escaped, at, leaves)
    assert all(3 not in m.cycle_lengths(g) for g in smaller)
    assert m.product(escaped) == m.product(initial)
    return len(word), [list(e) for e in star]


def star_controls():
    counts = []
    for t in (3, 4, 5):
        cases = steps = 0
        for tree in projection.trees(t):
            for es in permutations(tree):
                length, _ = isolation_control(t, es)
                steps += length; cases += 1
        counts.append({"terminal_size": t, "all_labeled_tree_orders": cases,
                       "total_normalizing_moves": steps})
    larger = []
    for t in (6, 12, 20):
        es = tuple((i, i + 1) for i in range(t - 1))
        for order in (es, es[::-1]):
            for sign in (-1, 1):
                length, _ = isolation_control(t, order, sign)
                larger.append({"terminal_size": t, "orientation": sign,
                               "star_normalizing_word_length": length})
    return {"complete_small_tree_orders": counts, "larger_path_controls": larger}


def unique_largest_controls():
    cases = 0; examples = []
    for s in (2, 3, 4):
        ell = 5; n = 2 * s + ell
        for tree in projection.trees(s):
            for es in permutations(tree):
                pairs = [m.perm(n, [(a, b), (a + s, b + s)]) for a, b in es]
                mm = [m.perm(n, [(2 * s + i, 2 * s + i + 1),
                                 (2 * s + i + 2, 2 * s + i + 3)]) for i in (0, 1)]
                state = (m.perm(n, [(0, s, 2 * s)]), mm[0], *pairs, mm[1])
                initial_rho = mixed_cut(state)[2]
                degree = Counter(x for e in es for x in e)
                leaf = next(x for x in range(1, s) if degree[x] == 1)
                at = next(j for j, D in enumerate(state) if j > 0 and leaf in support(D))
                gather = list(range(-at, -1))
                gathered = m.braid(state, gather)
                D = gathered[1]; root_side = next(x for x in support(D) if x < s and x != leaf)
                R = m.product(gathered[1:]); k = 0; x = 0
                while x != root_side:
                    x = R[x]; k += 1; assert k < s
                twist = list(range(2, len(state))) * (len(state) - 1)
                aligned = m.braid(gathered, twist * k)
                leaf_now = m.power(R, -k)[leaf]
                assert aligned[1] == m.perm(n, [(0, leaf_now), (s, leaf_now + s)])
                escaped = m.braid(aligned, [1, 1])
                _, cs, rho = mixed_cut(escaped)
                assert max(map(len, cs)) == ell + s - 1
                assert rho == initial_rho - (s - 1)
                assert m.product(escaped) == m.product(state)
                cases += 1
                if k and not any(e["smaller_size"] == s for e in examples):
                    examples.append({"smaller_size": s, "largest_size": ell,
                                     "paired_edges": [list(e) for e in es],
                                     "gather_word": gather, "twist_phase": k,
                                     "rho_drop": s - 1})
    return {"all_small_labeled_tree_orders": cases, "nonzero_phase_examples": examples}


def relabel(state, f):
    return tuple(m.conjugate(f, g) for g in state)


def fixed_product(state):
    P = m.product(state); seq = [0]
    for _ in range(len(P) - 1):
        seq.append(P[seq[-1]])
    assert len(set(seq)) == len(P) and P[seq[-1]] == 0
    f = [0] * len(P)
    for i, x in enumerate(seq):
        f[x] = i
    return relabel(state, tuple(f))


def mixed_count(n):
    r = (n - 1) // 2; L = n - 1
    total = Fraction(0)
    for j in range(L + 1):
        x = Fraction(j * (L - j), L * (L - 1))
        total += (-1) ** j * comb(L, j) * (1 - 3 * x) * (1 - 4 * x) ** (r - 1)
    total *= Fraction(2 * comb(n, 3) * (3 * comb(n, 4)) ** (r - 1), factorial(n))
    assert total.denominator == 1
    return int(total)


def pure_small():
    rows = []; representatives = []
    for n in (5, 7):
        r = (n - 1) // 2; D = h.dt_class(n); Ds = set(D)
        P = m.perm(n, [tuple(range(n))]); domain = set()
        for prefix in product(D, repeat=r - 1):
            last = m.compose(m.inverse(m.product(prefix)), P)
            if last in Ds:
                domain.add(prefix + (last,))
        assert len(domain) == (5 if n == 5 else 392)
        remaining = set(domain)
        while remaining:
            seed = min(remaining); orbit = {seed}; queue = deque([seed])
            while queue:
                state = queue.popleft()
                for j in range(r - 1):
                    for sign in (-1, 1):
                        new = m.hurwitz(state, j, sign)
                        assert new in domain
                        if new not in orbit:
                            orbit.add(new); queue.append(new)
            remaining -= orbit
            order = h.group_closure(seed)
            rows.append({"degree": n, "literal_orbit": len(orbit),
                         "absolute_orbit": len(orbit) // n,
                         "group_order": order,
                         "representative": [list(g) for g in seed]})
            if n == 7:
                representatives.append(seed)
    assert sorted((x["literal_orbit"], x["group_order"]) for x in rows if x["degree"] == 7) == [
        (49, 168), (49, 168), (147, 2520), (147, 2520)]
    return rows, representatives


def mixed_bases(pure7):
    rows = []; bridges = []
    for n in (3, 5, 7, 9):
        r = (n - 1) // 2
        if n == 3:
            seed = (m.perm(n, [(0, 1, 2)]),)
        elif n == 5:
            seed = (m.perm(n, [(0, 1, 2)]), m.perm(n, [(0, 3), (1, 4)]))
        else:
            seed = paired_tree(r, tuple((0, i) for i in range(1, r)))
        seed = fixed_product(seed); P = m.product(seed)
        words = {seed: ()}; queue = deque([seed])
        rho_counts = Counter({mixed_cut(seed)[2]: 1})
        while queue:
            state = queue.popleft()
            for j in range(r - 1):
                for sign in (-1, 1):
                    new = m.hurwitz(state, j, sign)
                    if new not in words:
                        assert m.product(new) == P
                        rho_counts[mixed_cut(new)[2]] += 1
                        words[new] = words[state] + (sign * (j + 1),)
                        queue.append(new)
        count = mixed_count(n)
        assert len(words) == r * count
        marked = {s for s in words if 3 in m.cycle_lengths(s[0])}
        assert len(marked) == count
        # Direct closure under the EXACT marked generators, not a larger group.
        seen = {seed}; queue = deque([seed])
        generators = [(j,) for j in range(2, r)] + [(1, 1)] if r > 1 else []
        generators += [tuple(-x for x in reversed(w)) for w in generators]
        while queue:
            state = queue.popleft()
            for word in generators:
                new = m.braid(state, word)
                assert new in marked
                if new not in seen:
                    seen.add(new); queue.append(new)
        assert seen == marked
        if n in (5, 7):
            D = h.dt_class(n); direct = set()
            for rest in product(D, repeat=r - 1):
                T = m.compose(P, m.inverse(m.product(rest)))
                if m.cycle_lengths(T) == [1] * (n - 3) + [3]:
                    direct.add((T,) + rest)
            assert direct == marked
        if n == 5:
            assert {mixed_cut(s)[2] for s in words} == {1}
        assert min(rho_counts) == (1 if n == 5 else 0)
        removable = []
        for at, D in enumerate(seed):
            if 3 in m.cycle_lengths(D):
                continue
            rest_support = set().union(*(support(g) for j, g in enumerate(seed) if j != at))
            if all(set(e) - rest_support for e in h.edges(D)):
                removable.append(at)
        assert bool(removable) == (n >= 5)
        if n == 9:
            for core in pure7:
                for anchor in range(7):
                    for orientation in (-1, 1):
                        extended = tuple(tuple(g) + (7, 8) for g in core) + (
                            m.perm(9, [(anchor, 7, 8)] if orientation > 0 else [(anchor, 8, 7)]),)
                        normalized = fixed_product(extended)
                        word = tuple(-x for x in reversed(words[normalized]))
                        assert m.braid(normalized, word) == seed
                        bridges.append({"pure_group_order": h.group_closure(core),
                                        "pure_representative": [list(g) for g in core],
                                        "anchor": anchor, "orientation": orientation,
                                        "word_to_common_mixed_seed": list(word)})
        digest = hashlib.sha256(b"".join(bytes(x for g in s for x in g)
                                        for s in sorted(words))).hexdigest()
        rows.append({"degree": n, "literal_unmarked_orbit": len(words),
                     "literal_exact_marked_orbit": len(marked),
                     "independent_character_count_marked": count,
                     "sorted_orbit_sha256": digest,
                     "literal_seed": [list(g) for g in seed],
                     "rho_distribution": dict(sorted(rho_counts.items())),
                     "seed_removable_DT_indices": removable,
                     "exact_marked_generators": [list(w) for w in generators],
                     "generated_group_order": h.group_closure(seed) if n <= 7 else None,
                     "group_scope": "exact closure" if n <= 7 else "A9 by written primitivity/Jordan proof"})
    return rows, bridges


def pure_canonical(r):
    n = 2 * r + 1
    return (m.perm(n, [(2 * r, 2), (0, 1)]),) + tuple(
        m.perm(n, [(0, 2 * j), (1, 2 * j + 1)]) for j in range(r - 1, 0, -1))


def anchor_controls():
    tested = 0
    for r in range(2, 9):
        core = pure_canonical(r); n = 2 * r + 1; C = m.product(core)
        expanded = tuple(tuple(g) + (n, n + 1) for g in core)
        f = tuple(C) + (n, n + 1)
        word = list(range(1, r)) * r
        for anchor in range(n):
            for sign in (-1, 1):
                T = m.perm(n + 2, [(anchor, n, n + 1)] if sign > 0 else [(anchor, n + 1, n)])
                state = expanded + (T,)
                moved = relabel(m.braid(state, word), f)
                assert moved[:-1] == expanded
                target = m.perm(n + 2, [(C[anchor], n, n + 1)] if sign > 0 else [(C[anchor], n + 1, n)])
                assert moved[-1] == target and mixed_cut(moved)[2] == 0
                tested += 1
    return {"all_anchors_and_orientations": tested, "pure_factor_counts": list(range(2, 9))}


def clifford_controls():
    def mul(A, B):
        out = {}
        for a, x in A.items():
            for b, y in B.items():
                swaps = sum((b & ((1 << i) - 1)).bit_count() for i in range(4) if a & (1 << i))
                k = a ^ b
                out[k] = out.get(k, 0) + (-1) ** swaps * x * y
        return {k: x for k, x in out.items() if x}
    def half(A):
        return {k: Fraction(x, 2) for k, x in A.items()}
    def vec(a, b):
        return {1 << a: Fraction(1), 1 << b: Fraction(-1)}
    D = half(mul(vec(0, 1), vec(2, 3)))
    H = half(mul(vec(0, 2), vec(1, 3)))
    minus = lambda A: {k: -x for k, x in A.items()}
    assert mul(D, D) == mul(H, H) == {0: -1}
    assert mul(mul(H, D), minus(H)) == minus(D)
    return {"DT_lift_square_minus_one": True,
            "even_centralizer_conjugates_the_two_lift_signs": True}


def negative_controls():
    good = paired_tree(3, ((0, 1), (0, 2)))
    repeated = good[:-1] + (good[1],)
    try:
        mixed_cut(repeated)
    except AssertionError:
        pass
    else:
        raise AssertionError("repeated DT edges were accepted")
    wrong = (m.perm(7, [(0, 1)]),) + good[1:]
    try:
        mixed_cut(wrong)
    except AssertionError:
        pass
    else:
        raise AssertionError("wrong distinguished branch class was accepted")
    # The literal identity alone is insufficient when a purported satellite
    # leaf is also moved by an outside factor. This is an admissible full tree.
    state = tuple(m.perm(9, cs) for cs in (
        ((0, 1, 2),), ((1, 3), (2, 4)), ((1, 5), (2, 6)), ((3, 7), (2, 8))))
    mixed_cut(state)
    changed = m.braid(state, [-1, -1, -1, 2, 2])
    assert mixed_cut(changed)[2] > 0
    return {"repeated_edge_tree_rejected": True, "wrong_branch_class_rejected": True,
            "outside_satellite_contact_prevents_false_isolation": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "mixed_hurwitz_isolation_verification.json")
    args = parser.parse_args()
    local = pair_controls(); stars = star_controls(); unique = unique_largest_controls()
    pure, representatives = pure_small()
    mixed, bridges = mixed_bases(representatives)
    anchors = anchor_controls(); lifting = clifford_controls(); negative = negative_controls()
    inputs = ["MIXED_HURWITZ_ISOLATION_AND_CONNECTIVITY.md",
              "TERMINAL_MASS_HURWITZ_DESCENT.md", "DEGREE_FOURTEEN_HURWITZ_ORBIT.md",
              "UNIFORM_TWO_LEAF_DT_ABSORPTION.md", "UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md",
              "HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md",
              "verify_terminal_mass_hurwitz_descent.py",
              "verify_terminal_mass_projection_escape.py",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt = {"schema_version": 1, "kind": "mixed_hurwitz_isolation_and_small_bases",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "local_identities_and_cut_controls": local, "paired_tree_isolation_controls": stars,
        "unique_largest_phase_alignment_controls": unique,
        "complete_small_pure_orbits": pure, "complete_small_mixed_bases": mixed,
        "all_56_degree7_core_extension_bridges": bridges,
        "one_anchor_restoration_controls": anchors, "lifting_controls": lifting,
        "negative_controls": negative,
        "scope": "Complete fixed-product mixed bases in degrees 3,5,7,9 and pure controls "
                 "in degrees 5,7. All unrestricted isolation/connectivity conclusions use "
                 "the written proof and retained pure-DT bases and absorption. No larger "
                 "census, independent replay, formal verification or external review.",
        "independent_replay": False, "formal_verification": False, "external_review": False}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "mixed_bases": mixed,
                      "pure7_bridges": len(bridges), "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
