#!/usr/bin/env python3
"""Exact controls for arbitrary minimal index-two mixtures and their boundaries."""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from math import comb, factorial
from pathlib import Path
import platform

import verify_mixed_hurwitz_isolation as one

if not __debug__:
    raise RuntimeError("mixture checker requires active assertions")

HERE = Path(__file__).resolve().parent
m, h = one.m, one.h


def kind(g):
    n = len(g); lengths = m.cycle_lengths(g)
    if lengths == [1] * (n - 3) + [3]:
        return "T"
    assert lengths == [1] * (n - 4) + [2, 2]
    return "D"


def cuts(state):
    n = len(state[0]); assert n == 2 * len(state) + 1
    types = [kind(g) for g in state]
    es = [e for g in state for e in h.edges(g)]
    assert len(es) == n - 1 and len(set(map(frozenset, es))) == n - 1
    assert len(one.v.components(n, es)) == 1
    answer = {}
    for i, typ in enumerate(types):
        if typ != "T":
            continue
        rest = [e for j, g in enumerate(state) if i != j for e in h.edges(g)]
        cs = one.v.components(n, rest)
        assert len(cs) == 3 and all(len(c & one.support(state[i])) == 1 for c in cs)
        answer[i] = cs
    return answer


def canonical(k, d):
    assert k >= 1 and d >= 0
    n = 2 * (k + d) + 1
    if d == 0:
        return tuple(m.perm(n, [(0, 2 * i + 1, 2 * i + 2)]) for i in range(k))
    if d == 1:
        core = (m.perm(n, [(0, 1, 2)]), m.perm(n, [(0, 3), (1, 4)]))
        return core + tuple(m.perm(n, [(0, 5 + 2 * i, 6 + 2 * i)]) for i in range(k - 1))
    old = one.pure_canonical(d); q = 2 * d + 1
    core = tuple(tuple(g) + tuple(range(q, n)) for g in old)
    return core + tuple(m.perm(n, [(0, q + 2 * i, q + 2 * i + 1)]) for i in range(k))


def pure_words(r):
    words = []
    for i in range(r - 1):
        for j in range(i + 1, r):
            gather = list(range(j, i + 1, -1))
            word = tuple(gather + [i + 1, i + 1] + [-x for x in reversed(gather)])
            assert m.labels(r, word) == tuple(range(r))
            words += [word, tuple(-x for x in reversed(word))]
    return words


def expected_count(n, k, d):
    L = n - 1; total = Fraction(0)
    for j in range(L + 1):
        x = Fraction(j * (L - j), L * (L - 1))
        total += (-1) ** j * comb(L, j) * (1 - 3 * x) ** k * (1 - 4 * x) ** d
    total *= Fraction((2 * comb(n, 3)) ** k * (3 * comb(n, 4)) ** d, factorial(n))
    assert total.denominator == 1
    assert total == n ** (k + d - 1) * (k + d - 1) ** d
    return int(total)


def comparison_controls(state, cs):
    tested = strict = 0
    for i, parts in cs.items():
        L = max(map(len, parts))
        for j, other in cs.items():
            if i == j:
                continue
            inside, = [V for V in parts if one.support(state[j]) <= V]
            outside = set().union(*(set(V) for V in parts if V != inside))
            big, = [V for V in other if outside <= V]
            root, = one.support(state[i]) & inside
            assert root in big and len(big) >= len(state[0]) - len(inside) + 1
            if len(inside) < L or sum(len(V) == L for V in parts) > 1:
                assert max(map(len, other)) >= L + 2
                strict += 1
            tested += 1
    return tested, strict


def complete_bases():
    rows = []; comparisons = strict_comparisons = 0
    for r in range(1, 5):
        for k in range(1, r + 1):
            d = r - k; n = 2 * r + 1
            seed = one.fixed_product(canonical(k, d)); P = m.product(seed)
            orbit = {seed}; queue = deque([seed]); maximum = 0; isolated = 0
            while queue:
                state = queue.popleft(); cs = cuts(state)
                assert len(cs) == k
                a, b = comparison_controls(state, cs)
                comparisons += a; strict_comparisons += b
                L = max(max(map(len, V)) for V in cs.values())
                maximum = max(maximum, L)
                isolated += any(max(map(len, V)) == n - 2 for V in cs.values())
                for i in range(r - 1):
                    for sign in (-1, 1):
                        new = m.hurwitz(state, i, sign)
                        if new not in orbit:
                            assert m.product(new) == P
                            orbit.add(new); queue.append(new)
            fixed = expected_count(n, k, d)
            assert len(orbit) == comb(r, k) * fixed
            assert bool(isolated) == ((k, d) != (1, 1))
            assert maximum == (2 if (k, d) == (1, 1) else n - 2)
            pure = {seed}; queue = deque([seed]); words = pure_words(r)
            pattern = tuple(map(kind, seed))
            while queue:
                state = queue.popleft()
                for word in words:
                    new = m.braid(state, word)
                    assert tuple(map(kind, new)) == pattern and new in orbit
                    if new not in pure:
                        pure.add(new); queue.append(new)
            assert len(pure) == fixed
            digest = hashlib.sha256(b"".join(bytes(x for g in s for x in g)
                                            for s in sorted(orbit))).hexdigest()
            pure_digest = hashlib.sha256(b"".join(bytes(x for g in s for x in g)
                                                 for s in sorted(pure))).hexdigest()
            rows.append({"degree": n, "three_cycles": k, "DTs": d,
                         "literal_full_orbit": len(orbit), "literal_fully_marked_orbit": len(pure),
                         "class_character_fixed_pattern_total": fixed,
                         "maximum_T_component": maximum, "states_with_isolated_T": isolated,
                         "orbit_sha256": digest, "pure_orbit_sha256": pure_digest,
                         "literal_seed": [list(g) for g in seed],
                         "group_order": h.group_closure(seed) if n <= 7 else None})
            print(json.dumps({"degree": n, "k": k, "m": d, "full": len(orbit),
                              "pure": len(pure)}), flush=True)
    return {"complete_domains": rows, "literal_T_pair_cut_comparisons": comparisons,
            "strict_nonlargest_cut_comparisons": strict_comparisons}


def stabilizer_word(state, i, j, period):
    assert i < j
    gather = list(range(j, i + 1, -1))
    gathered = m.braid(state, gather)
    assert gathered[i] == state[i] and gathered[i + 1] == state[j]
    word = gather + [i + 1] * period + [-x for x in reversed(gather)]
    assert m.braid(state, word) == state
    perm = list(range(len(state))); perm[i], perm[j] = perm[j], perm[i]
    assert m.labels(len(state), word) == tuple(perm)
    # Outside the adjacent selected pair, every permutation stays literal
    # throughout the local loop. Gathering has its separate endpoint scope.
    active = gathered
    for _ in range(period):
        active = m.braid(active, [i + 1])
        assert active[:i] == gathered[:i] and active[i + 2:] == gathered[i + 2:]
    assert active == gathered
    return word


def marking_controls():
    cases = loops = 0; examples = []
    settings = [(k, r - k) for r in range(1, 11) for k in range(1, r + 1)]
    settings += [(12, 8), (8, 12), (20, 0), (1, 20)]
    for k, d in settings:
        state = canonical(k, d); cuts(state)
        Ts = [i for i, g in enumerate(state) if kind(g) == "T"]
        Ds = [i for i, g in enumerate(state) if kind(g) == "D"]
        for j in Ts[1:]:
            word = stabilizer_word(state, Ts[0], j, 5); loops += 1
            if k == 3 and d == 2:
                examples.append({"k": k, "m": d, "pair": [Ts[0], j], "word": word})
        if d >= 2:
            # W_d is the first d factors: A then D_(d-1),...,D_1.
            assert Ds == list(range(d))
            stabilizer_word(state, 0, d - 1, 5); loops += 1
            for j in range(1, d - 1):
                stabilizer_word(state, j, d - 1, 3); loops += 1
        cases += 1
    disjoint = (m.perm(6, [(0, 1, 2)]), m.perm(6, [(3, 4, 5)]))
    assert m.braid(disjoint, [1, 1]) == disjoint
    assert m.braid(disjoint, [1] * 5) != disjoint
    # Test both orientations independently on each petal of a one-contact pair.
    for a in ((0, 1, 2), (0, 2, 1)):
        for b in ((0, 3, 4), (0, 4, 3)):
            pair = (m.perm(5, [a]), m.perm(5, [b]))
            stabilizer_word(pair, 0, 1, 5)
    return {"canonical_mixture_controls": cases, "literal_mark_transposition_loops": loops,
            "separated_T_pair_examples": examples, "four_T_pair_orientations": 4,
            "disjoint_T_pair_does_not_have_odd_local_period": True}


def pure_boundary_controls():
    rows, representatives = one.pure_small(); marked = []
    seed5 = tuple(tuple(g) for g in next(r for r in rows if r["degree"] == 5)["representative"])
    orbit5 = {seed5}; queue = deque([seed5]); words5 = pure_words(2)
    while queue:
        state = queue.popleft()
        for word in words5:
            new = m.braid(state, word)
            if new not in orbit5:
                orbit5.add(new); queue.append(new)
    assert len(orbit5) == 5
    for seed in representatives:
        words = pure_words(3); seen = {seed}; queue = deque([seed])
        while queue:
            state = queue.popleft()
            for word in words:
                new = m.braid(state, word)
                if new not in seen:
                    seen.add(new); queue.append(new)
        order = h.group_closure(seed)
        assert len(seen) == (147 if order == 2520 else 49)
        minimum = min(one.v.mu(s) for s in seen)
        assert minimum == (0 if order == 2520 else 1)
        marked.append({"group_order": order, "literal_pure_braid_orbit": len(seen),
                       "minimum_DT_mass": minimum})
    return {"complete_pure_small_orbits": rows, "literal_pure_braid_degree5_orbit": len(orbit5),
            "fully_marked_degree7_orbits": marked}


def collective_pruning_controls():
    """Test the new pruning step with other Ts, including nontrivial phases."""
    cases = nonzero_phases = 0; examples = []
    for s in (2, 3, 4):
        for tree in one.projection.trees(s):
            for es in permutations(tree):
                degree = Counter(x for edge in es for x in edge)
                leaf = next(x for x in range(1, s) if degree[x] == 1)
                for kt, dt in ((2, 0), (1, 1), (3, 1), (2, 2)):
                    core = canonical(kt, dt); ell = len(core[0]); n = 2 * s + ell
                    assert ell > s
                    large = [m.perm(n, [tuple(x + 2 * s for x in cycle)
                                       for cycle in permutation_cycles(g)]) for g in core]
                    pairs = [m.perm(n, [(a, b), (a + s, b + s)]) for a, b in es]
                    for order in (large, large[::-1]):
                        state = (m.perm(n, [(0, s, 2 * s)]), order[0],
                                 *pairs, *order[1:])
                        assert sorted(map(len, cuts(state)[0])) == [s, s, ell]
                        at = next(i for i, D in enumerate(state) if i and leaf in one.support(D))
                        gather = list(range(-at, -1)); gathered = m.braid(state, gather)
                        assert gathered[0] == state[0]
                        assert all(g[leaf] == leaf and g[leaf + s] == leaf + s
                                   for g in gathered[2:])
                        D = gathered[1]
                        anchor, = [x for x in one.support(D) if x < s and x != leaf]
                        R = m.product(gathered[1:]); phase = 0; x = 0
                        while x != anchor:
                            x = R[x]; phase += 1; assert phase < s
                        twist = list(range(2, len(state))) * (len(state) - 1)
                        aligned = m.braid(gathered, twist * phase)
                        leaf_now = m.power(R, -phase)[leaf]
                        assert aligned[1] == m.perm(n, [(0, leaf_now), (s, leaf_now + s)])
                        assert sorted(map(len, cuts(aligned)[0])) == [s, s, ell]
                        escaped = m.braid(aligned, [1, 1])
                        assert max(map(len, cuts(escaped)[0])) == ell + s - 1
                        assert m.product(escaped) == m.product(state)
                        assert sum(kind(g) == "T" for g in escaped) == kt + 1
                        cases += 1; nonzero_phases += bool(phase)
                        if phase and not any(e["smaller_size"] == s for e in examples):
                            examples.append({"smaller_size": s, "largest_size": ell,
                                             "other_Ts": kt, "large_DTs": dt,
                                             "gather_word": gather, "twist_phase": phase,
                                             "strict_component_gain": s - 1})
    return {"paired_tree_orders_with_other_Ts": cases,
            "nonzero_phase_controls": nonzero_phases, "examples": examples}


def permutation_cycles(g):
    remaining = set(range(len(g))); cycles = []
    while remaining:
        first = min(remaining); cycle = []; x = first
        while x in remaining:
            remaining.remove(x); cycle.append(x); x = g[x]
        assert x == first
        if len(cycle) > 1:
            cycles.append(tuple(cycle))
    return cycles


def clifford_mul(A, B):
    out = {}
    for a, x in A.items():
        for b, y in B.items():
            bits = a; swaps = 0
            while bits:
                low = bits & -bits; swaps += (b & (low - 1)).bit_count(); bits -= low
            z = a ^ b
            out[z] = out.get(z, 0) + (-1) ** swaps * x * y
    return {z: x for z, x in out.items() if x}


def clifford_power(A, exponent):
    answer = {0: Fraction(1)}
    while exponent:
        if exponent & 1:
            answer = clifford_mul(answer, A)
        A = clifford_mul(A, A); exponent >>= 1
    return answer


def T_lift(T):
    (a, b), (bb, c) = h.edges(T); assert b == bb
    u = {1 << a: Fraction(1), 1 << b: Fraction(-1)}
    v = {1 << b: Fraction(1), 1 << c: Fraction(-1)}
    answer = {z: x / 2 for z, x in clifford_mul(u, v).items()}
    assert clifford_power(answer, 3) == {0: 1}
    return answer


def spin_of(state):
    n = len(state[0]); F = {0: Fraction(1)}
    for T in state:
        F = clifford_mul(F, T_lift(T))
    answer = clifford_power(F, n)
    assert answer in ({0: 1}, {0: -1})
    return int(answer[0])


def spin_controls():
    rows = []
    for r in range(1, 5):
        value = spin_of(canonical(r, 0))
        expected = (-1) ** (r * (r - 1) // 2)
        assert value == expected
        rows.append({"three_cycle_factors": r, "degree": 2 * r + 1,
                     "canonical_odd_lift_spin": value})
    assert one.clifford_controls()["even_centralizer_conjugates_the_two_lift_signs"]
    return rows


def higher_index_controls():
    dihedral = []
    for n in (5, 9, 13, 25, 45, 65):
        A = tuple((-x) % n for x in range(n)); B = tuple((-x - 1) % n for x in range(n))
        P = m.perm(n, [tuple(range(n))]); assert m.product((A, B)) == P
        assert m.cycle_lengths(A) == m.cycle_lengths(B) == [1] + [2] * ((n - 1) // 2)
        orbit = set(); state = (A, B)
        while state not in orbit:
            orbit.add(state); state = m.hurwitz(state, 0)
        assert len(orbit) == n and h.group_closure((A, B)) == 2 * n
        divisor = next((d for d in range(2, n) if n % d == 0), None)
        if divisor:
            blocks = {frozenset(range(i, n, divisor)) for i in range(divisor)}
            assert all(frozenset(g[x] for x in V) in blocks for V in blocks for g in (A, B))
        dihedral.append({"degree": n, "factor_index": (n - 1) // 2,
                         "literal_pair_orbit": n, "group_order": 2 * n,
                         "nontrivial_block_residue_modulus": divisor})
    wreath = []
    for b in (3, 5, 7):
        n = 3 * b
        T = m.perm(n, [(0, b, 2 * b)])
        B = m.perm(n, [tuple(range(i * b, (i + 1) * b)) for i in range(3)])
        assert m.product((T, B)) == m.perm(n, [tuple(range(n))])
        assert m.power(B, b) == tuple(range(n))
        blocks = {frozenset((i, i + b, i + 2 * b)) for i in range(b)}
        for i in range(b):
            C = m.conjugate(m.power(B, i), T)
            assert frozenset(one.support(C)) in blocks
        assert all(frozenset(g[x] for x in V) in blocks for V in blocks for g in (T, B))
        order = h.group_closure((T, B)); assert order == b * 3 ** b
        es = [e for g in (T, B) for e in h.edges(g)] if b == 3 else None
        if es is not None:
            assert len(es) == n - 1 and len(one.v.components(n, es)) == 1
        wreath.append({"degree": n, "higher_cycle_length": b,
                       "higher_factor_index": 3 * (b - 1), "group_order": order})
    return {"minimal_dihedral_family": dihedral,
            "minimal_imprimitive_family_with_a_three_cycle": wreath}


def nonminimal_spin_split():
    n = 5; Ts = [m.perm(n, [c]) for v in combinations(range(n), 3)
                 for c in (v, (v[0], v[2], v[1]))]
    Ds = set(Ts); P = m.perm(n, [tuple(range(n))])
    domain = {(A, B, C) for A in Ts for B in Ts
              if (C := m.compose(m.inverse(m.compose(A, B)), P)) in Ds}
    assert len(domain) == 125
    remaining = set(domain); rows = []
    while remaining:
        seed = min(remaining); seen = {seed}; queue = deque([seed])
        while queue:
            state = queue.popleft()
            for j in (0, 1):
                for sign in (-1, 1):
                    new = m.hurwitz(state, j, sign); assert new in domain
                    if new not in seen:
                        seen.add(new); queue.append(new)
        remaining -= seen
        spins = {spin_of(s) for s in seen}; assert len(spins) == 1
        order = h.group_closure(seed); assert order == 60
        digest = hashlib.sha256(b"".join(bytes(x for g in s for x in g)
                                        for s in sorted(seen))).hexdigest()
        rows.append({"literal_orbit": len(seen), "absolute_orbit": len(seen) // 5,
                     "spin": spins.pop(), "group_order": order,
                     "literal_representative": [list(g) for g in seed],
                     "orbit_sha256": digest})
    assert sorted((r["literal_orbit"], r["spin"]) for r in rows) == [(50, -1), (75, 1)]
    try:
        cuts(min(domain))
    except AssertionError:
        pass
    else:
        raise AssertionError("nonminimal full-cycle tuple was admitted to the tree checker")
    return {"degree": 5, "factor_count": 3, "complete_domain": 125,
            "genus_with_cycle_at_infinity": 1, "orbits": rows,
            "rejected_by_minimal_tree_validator": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "index_two_hurwitz_mixtures_verification.json")
    args = parser.parse_args()
    bases = complete_bases(); marks = marking_controls(); pure = pure_boundary_controls()
    pruning = collective_pruning_controls()
    spin = spin_controls(); higher = higher_index_controls(); nonminimal = nonminimal_spin_split()
    inputs = ["INDEX_TWO_HURWITZ_MIXTURES.md", "MIXED_HURWITZ_ISOLATION_AND_CONNECTIVITY.md",
              "verify_mixed_hurwitz_isolation.py", "verify_terminal_mass_hurwitz_descent.py",
              "verify_terminal_mass_projection_escape.py", "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py",
              "sources/muehle-nadeau-alternating-3-cycles-2019.pdf",
              "sources/fried-lifting-invariants-math-0611591v4.pdf"]
    receipt = {"schema_version": 1, "kind": "arbitrary_index_two_hurwitz_mixtures",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "complete_minimal_mixture_controls": bases, "fully_marked_stabilizer_controls": marks,
        "collective_pruning_with_other_Ts": pruning,
        "pure_DT_boundary_controls": pure, "canonical_spin_controls": spin,
        "higher_index_boundary_controls": higher, "nonminimal_lifting_obstruction": nonminimal,
        "scope": "Bounded exact controls, complete minimal mixture domains through N=9, "
                 "and the complete nonminimal N=5 three-cycle triple domain. Universal "
                 "mixture connectivity is a written collective-cut proof and induction "
                 "from the one-T theorem; no unrestricted higher-index or nonminimal "
                 "connectivity is inferred. No independent replay or formal verification.",
        "independent_replay": False, "formal_verification": False, "external_review": False}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "domains": len(bases["complete_domains"]),
                      "T_pair_comparisons": bases["literal_T_pair_cut_comparisons"],
                      "pruning_controls_with_other_Ts": pruning["paired_tree_orders_with_other_Ts"],
                      "marking_loops": marks["literal_mark_transposition_loops"],
                      "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
