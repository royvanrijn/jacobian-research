#!/usr/bin/env python3
"""Exact finite controls of written terminal-mass descent lemmas, not global escape."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("terminal-mass checker requires active Python assertions")

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "paired_leaf_helpers", HERE / "verify_paired_hurwitz_leaf_extraction.py")
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
m = h.m


def components(n, edges):
    adjacency = [set() for _ in range(n)]
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    remaining, answer = set(range(n)), []
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        queue = [root]
        for a in queue:
            for b in sorted(adjacency[a] & remaining):
                remaining.remove(b)
                queue.append(b)
        answer.append(frozenset(queue))
    return tuple(answer)


def edge_list(t):
    n = len(t[0])
    assert n == 2 * len(t) + 1
    answer = []
    for g in t:
        assert len(g) == n and sorted(g) == list(range(n))
        assert m.cycle_lengths(g) == [1] * (n - 4) + [2, 2]
        answer.extend(tuple(sorted(e)) for e in h.edges(g))
    assert len(set(answer)) == n - 1
    assert len(components(n, answer)) == 1
    return answer


def factor_cut(t, index):
    es = edge_list(t)
    removed = {tuple(sorted(e)) for e in h.edges(t[index])}
    comp = components(len(t[0]), [e for e in es if e not in removed])
    assert len(comp) == 3
    pair = [set(e) for e in removed]
    middle = next(S for S in comp if all(S & e for e in pair))
    terminal = tuple(S for S in comp if S != middle)
    outer = frozenset(x for e in pair for x in e if x not in middle)
    assert len(outer) == 2
    return len(t[0]) - len(middle) - 2, middle, terminal, outer


def masses(t):
    return tuple(factor_cut(t, i)[0] for i in range(len(t)))


def mu(t):
    return min(masses(t))


def removable(t, index):
    """Independent test by supports, without component sizes."""
    rest = set().union(*(set(i for i, j in enumerate(g) if i != j)
                         for k, g in enumerate(t) if k != index))
    return all(any(a not in rest for a in e) for e in h.edges(t[index]))


def clean_contacts(t):
    mass = mu(t)
    answer = []
    for i, A in enumerate(t):
        val, _, _, outer = factor_cut(t, i)
        if val != mass:
            continue
        for j, B in enumerate(t):
            contact = {a for a in range(len(A)) if A[a] != a and B[a] != a}
            if len(contact) == 1 and contact <= outer:
                answer.append((i, j))
    return answer


def bring_B_word(i, j):
    assert i != j
    return (list(range(j, i + 1, -1)) if i < j
            else [-k for k in range(j + 1, i)])


def prune(t, i, j):
    mass, middle, terminal, outer = factor_cut(t, i)
    assert mass == mu(t) > 0
    A, B = t[i], t[j]
    common = {a for a in range(len(A)) if A[a] != a and B[a] != a}
    assert len(common) == 1 and common <= outer
    a = next(iter(common))
    e = B[a]
    S = next(S for S in terminal if a in S)
    U = next(S for S in terminal if a not in S)
    assert e in S
    initial_sets = frozenset((middle,) + terminal)
    word = bring_B_word(i, j)
    active = t
    current_i = i
    trajectory = [mass]
    for move in word:
        active = m.braid(active, [move])
        assert active[current_i] == A
        cut = factor_cut(active, current_i)
        assert frozenset((cut[1],) + cut[2]) == initial_sets
        assert cut[0] == mass
        trajectory.append(mu(active))
        assert trajectory[-1] <= mass
    final_move = -(i + 1) if i < j else i
    left = i if i < j else i - 1
    assert active[left:left + 2] == ((A, B) if i < j else (B, A))
    # Compute the pruned branch directly in the pre-move full tree.
    es = edge_list(active)
    split = components(len(A), [edge for edge in es
                                if set(edge) <= S and set(edge) != {a, e}])
    Sa = next(V for V in split if a in V)
    Se = S - Sa
    assert e in Se and Se
    active = m.braid(active, [final_move])
    current_i = left + 1 if i < j else left
    assert active[current_i] == A
    cut = factor_cut(active, current_i)
    assert frozenset(cut[2]) == frozenset((Sa, U))
    assert cut[1] == middle | Se
    assert cut[0] == mass - len(Se)
    trajectory.append(mu(active))
    assert trajectory[-1] < mass
    assert m.product(active) == m.product(t)
    full_word = word + [final_move]
    assert len(full_word) == abs(i - j)
    first_drop = next(k for k, val in enumerate(trajectory) if val < mass)
    prefix = full_word[:first_drop]
    assert all(val == mass for val in trajectory[:first_drop])
    assert mu(m.braid(t, prefix)) < mass
    return len(full_word), first_drop, len(Se)


def quotient_mass(pair, index, owners, weights):
    """Compute solely on a weighted five-vertex quotient."""
    es = [tuple(owners[a] for a in e) for g in pair for e in h.edges(g)]
    assert all(a != b for a, b in es)
    assert len({frozenset(e) for e in es}) == 4
    assert len(components(5, es)) == 1
    removed = {frozenset(owners[a] for a in e) for e in h.edges(pair[index])}
    comp = components(5, [e for e in es if frozenset(e) not in removed])
    assert len(comp) == 3
    middle = next(S for S in comp if all(S & e for e in removed))
    return sum(w for i, w in enumerate(weights) if i not in middle) - 2


def rotation_controls(t, at):
    A, B = t[at:at + 2]
    removed = {tuple(sorted(e)) for g in (A, B) for e in h.edges(g)}
    comp = components(len(A), [e for e in edge_list(t) if e not in removed])
    assert len(comp) == 5
    owners = {a: i for i, S in enumerate(comp) for a in S}
    weights = tuple(map(len, comp))
    C = m.compose(A, B)
    ident = tuple(range(len(A)))
    period = next(k for k in range(1, 7) if m.power(C, k) == ident)
    assert period in (2, 3, 4, 5, 6)
    for k in range(period):
        u = m.braid(t, [at + 1] * k)
        expected = (m.compose(A, m.power(C, k)),
                    m.compose(A, m.power(C, k + 1)))
        assert u[at:at + 2] == expected
        assert u[:at] == t[:at] and u[at + 2:] == t[at + 2:]
        for j in (0, 1):
            assert quotient_mass(expected, j, owners, weights) == factor_cut(u, at + j)[0]
    assert m.braid(t, [at + 1] * period) == t
    return period, 2 * period


def capture_controls(t, at):
    A, B = t[at:at + 2]
    common = {a for a in range(len(A)) if A[a] != a and B[a] != a}
    if len(common) != 2:
        return False
    if not all(len(set(e) & common) == 1 for g in (A, B) for e in h.edges(g)):
        return False
    if any(g[a] != a for i, g in enumerate(t) if i not in (at, at + 1) for a in common):
        return False
    u = m.braid(t, [-(at + 1)])
    assert u[at + 1] == A and removable(u, at + 1)
    assert factor_cut(u, at + 1)[0] == 0
    h.delete(u, at + 1, common)
    # Reversed factor order is an equally valid tree input with a different product.
    reversed_t = t[:at] + (B, A) + t[at + 2:]
    v = m.braid(reversed_t, [at + 1])
    assert v[at] == A and removable(v, at)
    h.delete(v, at, common)
    return True


def degree7_controls():
    ds = h.dt_class(7)
    D = set(ds)
    P = m.perm(7, [tuple(range(7))])
    domain = sorted((a, b, c) for a in ds for b in ds
                    if (c := m.compose(m.inverse(m.compose(a, b)), P)) in D)
    assert len(domain) == 392  # Reuses the previously retained finite domain.
    contacts = captures = quotient_checks = 0
    ordering = Counter()
    full_mass = Counter()
    proper_mass = Counter()
    periods = Counter()
    for t in domain:
        vals = masses(t)
        assert all((v == 0) == removable(t, i) for i, v in enumerate(vals))
        full = not h.orthogonal(list(map(h.support, t)))
        (full_mass if full else proper_mass)[min(vals)] += 1
        for i, j in clean_contacts(t) if min(vals) else []:
            prune(t, i, j)
            ordering["B_right" if i < j else "B_left"] += 1
            contacts += 1
        for at in range(2):
            period, checks = rotation_controls(t, at)
            periods[period] += 1
            quotient_checks += checks
            captures += capture_controls(t, at)
    assert full_mass == {0: 210, 1: 42, 2: 42}
    assert proper_mass == {1: 98}
    assert all(ordering[key] for key in ("B_left", "B_right"))
    assert h.group_closure(h.cyclic(3)) == 168
    return {
        "existing_fixed_product_domain": len(domain),
        "outer_contact_descents": contacts,
        "bringing_orders": dict(sorted(ordering.items())),
        "two_contact_captures_checked_in_both_orders": captures,
        "weighted_quotient_mass_checks": quotient_checks,
        "pair_period_distribution": dict(sorted(periods.items())),
        "full_A7_mass_distribution": dict(sorted(full_mass.items())),
        "proper_order168_mass_distribution": dict(sorted(proper_mass.items())),
        "proper_G3_group_order": 168,
    }


def certificate_control():
    data = json.loads((HERE / "terminal_mass_hurwitz_descent_certificate.json").read_text())
    t = tuple(m.perm(data["degree"], es) for es in data["tuple_edges"])
    assert edge_list(t) and masses(t) == tuple(data["factor_masses"])
    assert m.cycle_lengths(m.product(t)) == [15]
    assert (h.support(t[0]) & h.support(t[1])).bit_count() == 1
    three = m.power(m.compose(t[0], t[1]), 2)
    assert three == m.perm(15, [data["three_cycle_from_first_pair_product_squared"]])
    assert not clean_contacts(t)
    assert capture_controls(t, 1)
    trajectory, active = [mu(t)], t
    for move in data["word"]:
        active = m.braid(active, [move])
        trajectory.append(mu(active))
    assert trajectory == data["terminal_mass_trajectory"] == [1, 0]
    at = data["captured_factor_index_zero_based"]
    assert active[at] == t[1] and masses(t)[1] > mu(t)
    assert all(h.degrees(t)[a] == 2 and h.degrees(active)[a] == 1
               for a in data["new_leaves"])
    smaller = h.delete(active, at, data["new_leaves"])
    assert m.cycle_lengths(m.product(smaller)) == [13]
    return {
        "degree": 15, "initial_minimizing_factor_zero_based": 6,
        "initial_mass": 1, "minimizing_outer_contacts": 0,
        "word": data["word"], "mass_trajectory": trajectory,
        "captured_factor_was_not_minimizing": True,
        "both_captured_sheets_were_initially_degree_two": True,
        "lower_full_cycle_degree": 13,
        "full_A15_scope": "Written odd-degree primitivity and the literal 3-cycle, then Jordan; no group closure claimed",
    }


def plateau_certificate_control():
    data = json.loads((HERE / "terminal_mass_hurwitz_descent_certificate.json").read_text())["plateau_witness"]
    t = tuple(m.perm(data["degree"], es) for es in data["tuple_edges"])
    assert edge_list(t) and masses(t) == tuple(data["factor_masses"])
    assert m.cycle_lengths(m.product(t)) == [15]
    assert not clean_contacts(t)
    assert (h.support(t[1]) & h.support(t[2])).bit_count() == 1
    assert m.power(m.compose(t[1], t[2]), 2) == m.perm(
        15, [data["three_cycle_from_second_pair_product_squared"]])
    profiles = []
    for at in range(6):
        period, _ = rotation_controls(t, at)
        profile = [mu(m.braid(t, [at + 1] * k)) for k in range(period)]
        assert all(v == 1 for v in profile)
        profiles.append(profile)
    assert profiles == data["initial_pair_rotation_mass_profiles"]
    assert sum(map(len, profiles)) == 24
    # Reconstruct both local conjugations independently of the stored endpoints.
    expected_flat = tuple(m.perm(15, es) for es in (
        ((1, 2), (13, 14)), ((6, 7), (10, 11)),
        ((3, 5), (9, 10)), ((3, 4), (11, 12)),
        ((8, 9), (12, 13)), ((0, 1), (2, 3)), ((5, 6), (7, 8))))
    flat = m.braid(t, [-3])
    assert flat == expected_flat and mu(flat) == 1
    assert masses(flat)[3] == 2 and not clean_contacts(flat)
    rotation_controls(flat, 3)
    active = m.braid(flat, [-4])
    assert active[4] == t[2]
    assert mu(active) == 0
    assert [mu(t), mu(flat), mu(active)] == data["terminal_mass_trajectory"]
    assert data["word"] == [-3, -4] and m.braid(t, data["word"]) == active
    assert all(h.degrees(t)[a] == 2 and h.degrees(active)[a] == 1
               for a in data["new_leaves"])
    smaller = h.delete(active, data["captured_factor_index_zero_based"], data["new_leaves"])
    assert m.cycle_lengths(m.product(smaller)) == [13]
    return {"degree": 15, "initial_mass": 1,
            "minimizing_outer_contacts": 0,
            "all_single_adjacent_pair_powers_have_mass_one": True,
            "pair_periods": list(map(len, profiles)),
            "complete_pair_power_cases": sum(map(len, profiles)),
            "word": data["word"], "mass_trajectory": data["terminal_mass_trajectory"],
            "captured_original_factor_mass": masses(t)[2],
            "captured_factor_mass_after_flat_step": masses(flat)[3],
            "new_leaves": data["new_leaves"], "lower_full_cycle_degree": 13,
            "scope": "Complete local non-descent for every power of one initially adjacent pair, not a global orbit obstruction"}


def larger_controls():
    answer = []
    for r in (4, 7, 12, 25):
        t = h.cyclic(r)
        assert mu(t) == 1
        contacts = clean_contacts(t)
        assert contacts
        i, j = contacts[0]
        length, first_drop, shrink = prune(t, i, j)
        # Reverse the branch order to exercise the other bringing convention.
        reversed_t = tuple(reversed(t))
        assert mu(reversed_t) == 1
        prune(reversed_t, r - 1 - i, r - 1 - j)
        rotation_controls(t, 0)
        answer.append({"r": r, "degree": 2 * r + 1,
                       "pruning_word_length": length,
                       "first_strict_drop_step": first_drop,
                       "terminal_component_shrink": shrink,
                       "both_bringing_orders_checked": True})
    return answer


def disjoint_pair_control():
    # Degree seven cannot contain a disjoint-support DT pair. This degree-nine
    # tree exercises the missing period-two case with nontrivial component weights.
    t = tuple(m.perm(9, es) for es in (
        ((0, 1), (2, 3)), ((4, 5), (6, 7)),
        ((1, 4), (6, 8)), ((3, 6), (5, 8))))
    assert m.cycle_lengths(m.product(t)) == [9]
    period, checks = rotation_controls(t, 0)
    assert period == 2 and checks == 4
    return {"degree": 9, "pair_period": period, "weighted_mass_checks": checks}


def invalid_controls():
    def rejected(action):
        try:
            action()
        except (AssertionError, StopIteration):
            return True
        raise AssertionError("Invalid control was accepted")
    bad_tree = (m.perm(7, [(0, 1), (2, 3)]),) * 3
    assert rejected(lambda: edge_list(bad_tree))
    good = h.cyclic(4)
    assert rejected(lambda: prune(good, 0, 1))  # Overlap three, not one.
    data = json.loads((HERE / "terminal_mass_hurwitz_descent_certificate.json").read_text())
    t = tuple(m.perm(15, es) for es in data["tuple_edges"])
    # The opposite sign captures the other original factor instead. It cannot
    # certify the specified factor/position in the retained word certificate.
    assert not removable(m.braid(t, [2]), 2)
    # Add an outside interaction at pivot 3; the one-move formula still holds,
    # but the claimed degree-two/exclusivity hypothesis must fail.
    extended = tuple(p + (15, 16) for p in t)
    extended += (m.perm(17, [(3, 15), (0, 16)]),)
    assert m.cycle_lengths(m.product(extended)) == [17]
    assert not capture_controls(extended, 1)
    assert not removable(m.braid(extended, [-2]), 2)
    return {"repeated_edge_tree_rejected": True,
            "wrong_overlap_for_pruning_rejected": True,
            "wrong_target_for_opposite_capture_sign_detected": True,
            "external_pivot_interaction_rejected": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "terminal_mass_hurwitz_descent_verification.json")
    args = parser.parse_args()
    controls = degree7_controls()
    certificate = certificate_control()
    plateau = plateau_certificate_control()
    larger = larger_controls()
    disjoint = disjoint_pair_control()
    invalid = invalid_controls()
    inputs = ["TERMINAL_MASS_HURWITZ_DESCENT.md",
              "terminal_mass_hurwitz_descent_certificate.json",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py",
              "alternating_hurwitz_source_ledger.json"]
    receipt = {
        "schema_version": 1, "kind": "terminal_mass_descent_controls", "result": "PASS",
        "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "degree7_mechanism_controls": controls,
        "degree15_changing_candidate_certificate": certificate,
        "degree15_local_plateau_certificate": plateau,
        "cyclic_larger_controls": larger, "invalid_controls": invalid,
        "disjoint_support_pair_control": disjoint,
        "scope": "Written cut-preserving bringing, strict outer-contact pruning, adjacent two-contact capture and five-weight pair rotation. A full-A15 input is flat under every power of each initial adjacent pair but extracts after a two-pair word. Finite controls do not prove global plateau escape or arbitrary extraction.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
