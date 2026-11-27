#!/usr/bin/env python3
"""Exact component certificates and block budgets for polynomial passports.

The universal classification remains open. Complete finite coverage here is
proved only when disjoint signed braid closures exhaust an independent exact
class-character total, or a directly enumerated two-factor domain.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, permutations, product
import json
from math import comb, factorial, gcd, lcm
from pathlib import Path
import platform

import verify_index_two_hurwitz_mixtures as v

if not __debug__:
    raise RuntimeError("component checker requires active assertions")

HERE = Path(__file__).resolve().parent
m, h = v.m, v.h


def cycles(g, include_fixed=True):
    remaining = set(range(len(g))); answer = []
    while remaining:
        first = min(remaining); x = first; cycle = []
        while x in remaining:
            remaining.remove(x); cycle.append(x); x = g[x]
        assert x == first
        if include_fixed or len(cycle) > 1:
            answer.append(tuple(cycle))
    return answer


def index(g):
    return len(g) - len(cycles(g))


def validate(state):
    n = len(state[0])
    assert all(len(g) == n and sorted(g) == list(range(n)) and index(g) for g in state)
    assert m.product(state) == m.perm(n, [tuple(range(n))])
    assert sum(map(index, state)) == n - 1
    # H1's incidence tree, independent of a transposition expansion convention.
    nodes = [(i, c) for i, g in enumerate(state) for c in cycles(g, False)]
    edges = [(x, n + j) for j, (_, c) in enumerate(nodes) for x in c]
    assert len(edges) == n + len(nodes) - 1
    assert len(v.one.v.components(n + len(nodes), edges)) == 1


def hook_characters(partition):
    """det(1+x*g | standard representation), by exact polynomial division."""
    n = sum(partition); coefficients = [1]
    for ell in partition:
        new = [0] * (len(coefficients) + ell)
        for j, a in enumerate(coefficients):
            new[j] += a; new[j + ell] -= (-1) ** ell * a
        coefficients = new
    answer = []; previous = 0
    for a in coefficients:
        q = a - previous; answer.append(q); previous = q
    assert answer[-1] == 0 and len(answer) == n + 1
    return answer[:-1]


def character_total(passport):
    n = sum(passport[0]); r = len(passport)
    assert all(sum(p) == n for p in passport)
    assert sum(n - len(p) for p in passport) == n - 1
    characters = [hook_characters(p) for p in passport]
    sizes = []
    for p in passport:
        denominator = 1
        for ell, multiplicity in Counter(p).items():
            denominator *= ell ** multiplicity * factorial(multiplicity)
        sizes.append(factorial(n) // denominator)
    total = Fraction(0)
    for j in range(n):
        term = Fraction((-1) ** j, comb(n - 1, j) ** (r - 1))
        for chi in characters:
            term *= chi[j]
        total += term
    for size in sizes:
        total *= size
    total /= factorial(n)
    assert total.denominator == 1 and total >= 0
    return int(total)


def comparison_total(passport):
    n = sum(passport[0]); answer = Fraction(n ** (len(passport) - 1))
    for p in passport:
        factor = Fraction(factorial(len(p) - 1))
        for count in Counter(p).values():
            factor /= factorial(count)
        answer *= factor
    assert answer.denominator == 1
    return int(answer)


def type_placements(passport):
    return len(set(permutations(tuple(tuple(p) for p in passport))))


def closure(seed, words=None):
    seen = {seed}; queue = deque([seed]); parent = {seed: (None, ())}
    moves = words if words is not None else [
        (sign * (j + 1),) for j in range(len(seed) - 1) for sign in (-1, 1)]
    while queue:
        state = queue.popleft()
        for word in moves:
            new = m.braid(state, word)
            if new not in seen:
                seen.add(new); queue.append(new); parent[new] = (state, tuple(word))
    return seen, parent


def word_to_seed(state, parent):
    word = []
    while parent[state][0] is not None:
        previous, edge = parent[state]
        word.extend(-x for x in reversed(edge)); state = previous
    return word


def rotated(state, shift):
    n = len(state[0])
    return tuple(tuple((g[(x - shift) % n] + shift) % n for x in range(n)) for g in state)


def cyclic_label(state):
    return min(rotated(state, k) for k in range(len(state[0])))


def pair_label(state):
    return min(cyclic_label(state), cyclic_label(m.hurwitz(state, 0)))


def block_divisors(state):
    n = len(state[0])
    return [d for d in range(1, n + 1) if n % d == 0 and
            all(g[(x + d) % n] % d == g[x] % d for g in state for x in range(n))]


def quotient_budget(state, d):
    n = len(state[0]); b = n // d; assert d in block_divisors(state)
    quotient = []; decorations = []; inner = 0
    for g in state:
        q = tuple(g[a] % d for a in range(d)); quotient.append(q); decorated = []
        for O in cycles(q):
            a = min(O); power = m.power(g, len(O))
            ret = tuple((power[a + t * d] - a) // d for t in range(b))
            assert sorted(ret) == list(range(b))
            inner += index(ret); decorated.append((len(O), tuple(m.cycle_lengths(ret))))
        assert index(g) == b * index(q) + sum(
            b - len(return_type) for _, return_type in decorated)
        decorations.append(tuple(sorted(decorated)))
    assert m.product(tuple(quotient)) == m.perm(d, [tuple(range(d))])
    assert sum(map(index, quotient)) == d - 1 and inner == b - 1
    return tuple(quotient), tuple(decorations)


def group_word(state, word):
    value = tuple(range(len(state[0])))
    for letter in word:
        g = state[abs(letter) - 1]
        value = m.compose(value, g if letter > 0 else m.inverse(g))
    return value


def rooted_plane_tree_code(state):
    A, B = state
    def at(edge, color):
        g = (A, B)[color]; children = []; x = g[edge]
        while x != edge:
            children.append(at(x, 1 - color)); x = g[x]
        return color, tuple(children)
    return min((at(e, 0), at(e, 1)) for e in range(len(A)))


def double_threes(n):
    for S in combinations(range(n), 6):
        for rest in combinations(S[1:], 2):
            X = (S[0],) + rest; Y = tuple(x for x in S if x not in X)
            for reverse_x, reverse_y in product((False, True), repeat=2):
                x = (X[0], X[2], X[1]) if reverse_x else X
                y = (Y[0], Y[2], Y[1]) if reverse_y else Y
                yield m.perm(n, [x, y])


def multiply_vector(A, a, b):
    answer = {}
    for mask, coefficient in A.items():
        for at, sign in ((a, 1), (b, -1)):
            target = mask ^ (1 << at)
            value = coefficient * sign * (-1) ** ((mask >> (at + 1)).bit_count())
            answer[target] = answer.get(target, 0) + value
    return {mask: c for mask, c in answer.items() if c}


@lru_cache(None)
def raw_cycle_lift(n):
    """Integer numerator of the adjacent-root lift of P_n, denominator 2^r."""
    r = (n - 1) // 2; assert n % 2
    raw = {0: 1}
    for j in range(n - 1):
        raw = multiply_vector(raw, j, j + 1)
    powered = {0: 1}
    # Multiply by vectors, avoiding a dense Clifford-by-Clifford power.
    for _ in range(n):
        for j in range(n - 1):
            powered = multiply_vector(powered, j, j + 1)
    assert len(powered) == 1 and 0 in powered
    scale = 2 ** (r * n); assert abs(powered[0]) == scale
    sign = powered[0] // scale
    assert sign == (-1) ** (r * (r - 1) // 2)
    return raw, sign


def spin_control(state):
    """Exact standard spin for the retained all-3-cycle-power examples."""
    n = len(state[0]); F = {0: 1}; count = 0
    for g in state:
        for C in cycles(g, False):
            assert len(C) == 3
            F = multiply_vector(F, C[0], C[1])
            F = multiply_vector(F, C[1], C[2]); count += 1
    assert 2 * count == n - 1
    raw, sign = raw_cycle_lift(n)
    ratio = 1 if F == raw else -1
    assert F == {mask: ratio * c for mask, c in raw.items()}
    spin = ratio * sign
    expected = (-1) ** (count + (n * n - 1) // 8)
    assert spin == expected
    return spin


def digest(states):
    return hashlib.sha256(b"".join(bytes(x for g in s for x in g)
                                   for s in sorted(states))).hexdigest()


def pair_controls(inputs):
    rows = []; representatives11 = []; labels_tested = 0
    certificates = {tuple(tuple(g) for g in c["representative"]): c
                    for c in inputs["degree11_pair_certificates"]}
    for n, second_count in ((9, 2), (11, 3)):
        P = m.perm(n, [tuple(range(n))])
        target = [1] * (n - 3 * second_count) + [3] * second_count
        candidates = list(double_threes(n))
        assert len(candidates) == factorial(n) // (3 ** 2 * 2 * factorial(n - 6))
        states = {(A, B) for A in candidates
                  if m.cycle_lengths(B := m.compose(m.inverse(A), P)) == target}
        passport = [m.cycle_lengths(g) for g in min(states)]
        assert len(states) == character_total(passport) == comparison_total(passport)
        domain = states | {m.hurwitz(s, 0) for s in states}
        remaining = set(states); orbits = []; pure_labels = set()
        while remaining:
            seed = min(remaining); validate(seed)
            orbit, _ = closure(seed); assert orbit <= domain
            remaining -= orbit
            normal = pair_label(seed)
            for state in orbit:
                assert pair_label(state) == normal
                assert m.braid(state, [1, 1]) == rotated(state, -1)
                validate(state); labels_tested += 1
            tree_code = rooted_plane_tree_code(seed)
            assert len(block_divisors(seed)) == 2
            pure_labels.update(cyclic_label(s) for s in orbit & states)
            row = {"literal_orbit": len(orbit), "absolute_orbit_size": len({
                       cyclic_label(s) for s in orbit}), "literal_representative": seed,
                   "two_strand_normal_form": normal, "ordered_plane_tree_code": tree_code,
                   "orbit_sha256": digest(orbit), "standard_spin": spin_control(seed)}
            if n == 9:
                row["group_order"] = h.group_closure(seed); assert row["group_order"] == 1512
            else:
                cert = certificates[seed]
                value = group_word(seed, cert["group_three_cycle_word"])
                assert list(value) == cert["word_value"]
                assert m.cycle_lengths(value) == [1] * 8 + [3]
                row.update(group="A11", group_basis="prime-degree primitivity + explicit 3-cycle + Jordan",
                           group_three_cycle_word=cert["group_three_cycle_word"], group_word_value=value)
                representatives11.append(seed)
            orbits.append(row)
        assert len({pair_label(s) for s in states}) == len(orbits)
        assert len({rooted_plane_tree_code(s) for s in states}) == len(pure_labels)
        assert sorted(len(closure(s, v.pure_words(2))[0]) for s in states) == [n] * len(states)
        assert (len(states), len(pure_labels), sorted(r["literal_orbit"] for r in orbits)) == (
            (36, 4, [9, 9, 18]) if n == 9 else (66, 6, [22] * 6))
        rows.append({"degree": n, "ordered_passport": passport,
                     "complete_ordered_domain": len(states), "all_placement_domain": len(domain),
                     "absolute_color_ordered_components": len(pure_labels), "components": orbits})
    return {"complete_two_factor_domains": rows, "normal_form_state_controls": labels_tested}, representatives11


def stabilizer_size(state):
    n = len(state[0])
    return sum(all(g[(x + k) % n] == (g[x] + k) % n
                   for g in state for x in range(n)) for k in range(n))


def imprimitive15():
    n = 15; qa = (1, 2, 0, 3, 4); qb = (0, 1, 3, 4, 2)
    assert h.group_closure((qa, qb)) == 60
    A = tuple(3 * qa[x // 3] + x % 3 for x in range(n))
    B = tuple(3 * qb[x // 3] + x % 3 for x in range(n))
    return v.one.fixed_product((A, B, m.perm(n, [(0, 1, 2)])))


def degree15_controls(inputs):
    primitive = tuple(tuple(g) for g in inputs["degree15_primitive_seed"])
    imprimitive = imprimitive15(); known = set(); rows = []
    passport = [m.cycle_lengths(g) for g in primitive]
    assert passport == [[1] * 6 + [3] * 3, [1] * 6 + [3] * 3, [1] * 12 + [3]]
    ordered_total = character_total(passport)
    assert ordered_total == comparison_total(passport) == 19600
    assert Fraction(factorial(8), factorial(3) * factorial(6)) == Fraction(28, 3)
    for seed, expected, group in ((primitive, 58725, "A15"), (imprimitive, 75, "C3 wr A5")):
        validate(seed); orbit, _ = closure(seed)
        assert not known & orbit and len(orbit) == expected; known |= orbit
        assert {tuple(sorted(tuple(m.cycle_lengths(g)) for g in s)) for s in orbit} == {
            tuple(sorted(tuple(p) for p in passport))}
        stabilizer = 1 if group == "A15" else 3
        for state in orbit:
            assert stabilizer_size(state) == stabilizer
        divisors = block_divisors(seed); budgets = []
        for d in divisors:
            quotient, decorations = quotient_budget(seed, d)
            budgets.append({"block_count": d, "block_size": 15 // d,
                            "quotient": quotient, "return_decorations": decorations})
        assert divisors == ([1, 15] if group == "A15" else [1, 5, 15])
        pure, _ = closure(seed, v.pure_words(3))
        assert pure <= orbit and len(pure) * 3 == len(orbit)
        row = {"group": group, "literal_orbit": len(orbit),
               "literal_pure_braid_fixed_pattern_orbit": len(pure),
               "cyclic_stabilizer_order": stabilizer,
               "absolute_orbit_size": len(orbit) // (15 // stabilizer),
               "literal_representative": seed,
               "standard_spin": spin_control(seed), "block_divisors": divisors,
               "quotient_index_budgets": budgets, "orbit_sha256": digest(orbit),
               "pure_orbit_sha256": digest(pure)}
        if group != "A15":
            assert h.group_closure(seed) == 14580
            quotient, decorations = quotient_budget(seed, 5)
            quotient_orbit, _ = closure(quotient); assert len(quotient_orbit) == 15
            for state in orbit:
                q, dec = quotient_budget(state, 5)
                assert q in quotient_orbit
                assert sorted(dec) == sorted(decorations)
            row.update(group_order=14580, quotient_literal_orbit=15, literal_fibre_size=5)
        rows.append(row)
        print(json.dumps({"degree": 15, "component": group, "full": len(orbit),
                          "pure": len(pure)}), flush=True)
    assert len(known) == ordered_total * type_placements(passport) == 58800
    return {"ordered_passport": passport, "independent_ordered_total": ordered_total,
            "all_placement_total": len(known), "components": rows,
            "every_cyclic_stabilizer_checked": True}


def refinement_controls(representatives):
    refinements = []
    for i, pair in enumerate(representatives):
        A, B = pair; C = cycles(B, False)
        assert len(C) == 3
        for removed in C:
            refined = (A, m.perm(11, [c for c in C if c != removed]), m.perm(11, [removed]))
            validate(refined)
            assert m.product(refined) == m.product(pair)
            refinements.append((i, removed, refined))
    seed = refinements[0][2]; orbit, parents = closure(seed)
    passport = [m.cycle_lengths(g) for g in seed]
    assert len(orbit) == character_total(passport) * type_placements(passport) == 3267
    assert character_total(passport) == comparison_total(passport)
    words = []
    for i, removed, refined in refinements:
        word = word_to_seed(refined, parents)
        assert m.braid(refined, word) == seed
        words.append({"rigid_component_index": i, "split_cycle": removed,
                      "refined_literal_tuple": refined, "braid_word_to_common_seed": word})
    # A second nonrigid all-odd-order case tests the arbitrary hook input.
    n = 13; terms = [[(0, 1, 2), (6, 7, 8)], [(2, 3, 4), (8, 9, 10)],
                     [(4, 5, 6), (10, 11, 12)]]
    seed13 = v.one.fixed_product(tuple(m.perm(n, C) for C in terms)); validate(seed13)
    orbit13, _ = closure(seed13); passport13 = [m.cycle_lengths(g) for g in seed13]
    assert len(orbit13) == character_total(passport13) == comparison_total(passport13) == 10816
    return {"degree11_refined_literal_orbit": len(orbit),
            "degree11_seed": seed, "degree11_orbit_sha256": digest(orbit),
            "all18_rigid_refinement_bridges": words, "degree11_standard_spin": spin_control(seed),
            "degree13_control": {"passport": passport13, "literal_seed": seed13,
                                "literal_orbit": len(orbit13), "orbit_sha256": digest(orbit13),
                                "standard_spin": spin_control(seed13)}}


def boundary_controls():
    rows = []
    for n in (9, 15, 25):
        A = tuple((-x) % n for x in range(n)); B = tuple((-x - 1) % n for x in range(n))
        seed = (A, B); validate(seed); ds = block_divisors(seed)
        for a in ds:
            for b in ds:
                assert gcd(a, b) in ds and lcm(a, b) in ds
        for d in ds:
            quotient_budget(seed, d)
        assert ds == [d for d in range(1, n + 1) if n % d == 0]
        assert block_divisors(m.hurwitz(seed, 0)) == ds
        rows.append({"dihedral_degree": n, "complete_block_divisors": ds})
    bad = (m.perm(5, [(0, 1, 2)]), m.perm(5, [(0, 2, 3)]),
           m.perm(5, [(0, 3, 4)]))
    assert m.product(bad) == m.perm(5, [(0, 1, 2, 3, 4)])
    assert sum(map(index, bad)) == 6
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError("nonminimal genus-one input admitted")
    # A tree's plane order is part of the label, not just its abstract graph.
    return {"block_lattice_controls": rows, "nonminimal_full_cycle_rejected": True,
            "rounding_higher_class_weights_is_not_used": True}


def complete_small_pair_controls():
    """Independent all-cycle-type controls, including even degrees and deck groups."""
    rows = []
    for n in range(3, 7):
        P = m.perm(n, [tuple(range(n))]); grouped = {}; domain = set()
        for A in permutations(range(n)):
            B = m.compose(m.inverse(A), P)
            if not index(A) or not index(B) or index(A) + index(B) != n - 1:
                continue
            state = (A, B); validate(state); domain.add(state)
            key = tuple(tuple(m.cycle_lengths(g)) for g in state)
            grouped.setdefault(key, set()).add(state)
        for passport, states in grouped.items():
            assert len(states) == character_total(passport) == comparison_total(passport)
        labels = {}; symmetric = 0
        for state in domain:
            assert m.braid(state, [1, 1]) == rotated(state, -1)
            label = pair_label(state); labels.setdefault(label, set()).add(state)
            assert pair_label(m.hurwitz(state, 0)) == label
            symmetric += stabilizer_size(state) > 1
        for seed, states in labels.items():
            orbit, _ = closure(seed)
            assert orbit == states
        rows.append({"degree": n, "all_minimal_nonidentity_pairs": len(domain),
                     "ordered_passports": len(grouped), "full_components": len(labels),
                     "pairs_with_nontrivial_cyclic_stabilizer": symmetric})
    assert any(r["pairs_with_nontrivial_cyclic_stabilizer"] for r in rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "polynomial_passport_components_verification.json")
    args = parser.parse_args()
    inputs = json.loads((HERE / "polynomial_passport_components_inputs.json").read_text())
    pairs, representatives = pair_controls(inputs)
    small = complete_small_pair_controls()
    fifteen = degree15_controls(inputs); refinement = refinement_controls(representatives)
    boundaries = boundary_controls()
    paths = ["POLYNOMIAL_PASSPORT_COMPONENTS.md", "polynomial_passport_components_inputs.json",
             "verify_index_two_hurwitz_mixtures.py", "verify_mixed_hurwitz_isolation.py",
             "verify_terminal_mass_hurwitz_descent.py", "verify_terminal_mass_projection_escape.py",
             "verify_paired_hurwitz_leaf_extraction.py", "verify_unmarked_polynomial_hurwitz_reduction.py",
             "HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md",
             "sources/jones-zvonkin-braid-cacti-2002.pdf",
             "sources/fried-lifting-invariants-math-0611591v4.pdf"]
    receipt = {"schema_version": 1, "kind": "polynomial_passport_component_reductions",
               "result": "PASS", "python": platform.python_version(),
               "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in paths},
               "two_factor_classifications": pairs, "degree15_component_classification": fifteen,
               "complete_all_passport_two_factor_controls_through_degree6": small,
               "refinement_controls": refinement, "boundary_controls": boundaries,
               "scope": "Complete finite domains explicitly listed, written all-degree two-factor "
                        "classification and block-budget reductions. The arbitrary-passport full "
                        "component classification remains open; no extrapolation from finite orbits.",
               "independent_replay": False, "formal_verification": False, "external_review": False}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "degree11_rigid_components": 6,
                      "degree15_components": 2, "refinement_bridges": 18,
                      "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
