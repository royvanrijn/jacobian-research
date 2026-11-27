#!/usr/bin/env python3
"""Exact controls for complementary-twist compatibility; no all-degree extraction claim."""
from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import importlib.util
import json
from math import gcd, lcm
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("global compatibility controls require active assertions")

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "terminal_mass_helpers", HERE / "verify_terminal_mass_hurwitz_descent.py")
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
m, h = v.m, v.h


def support(g):
    return frozenset(i for i, j in enumerate(g) if i != j)


def front(t, at):
    """Conjugate the selected factor leftwards; leave crossed factors unchanged."""
    initial = v.factor_cut(t, at)
    sets = frozenset((initial[1],) + initial[2])
    word = list(range(-at, 0))
    active = t
    for move in word:
        active = m.braid(active, [move])
        cut = v.factor_cut(active, -move - 1)
        assert cut[0] == initial[0]
        assert frozenset((cut[1],) + cut[2]) == sets
    assert m.product(active) == m.product(t)
    return active, word


def cycle_coordinates(R, root, component):
    sequence = [root]
    for _ in range(len(component) - 1):
        sequence.append(R[sequence[-1]])
    assert set(sequence) == set(component) and R[sequence[-1]] == root
    return {a: i for i, a in enumerate(sequence)}


def phase_profile(t):
    """A is first. Test ALL complementary-twist phases by exact permutation powers."""
    A, rest = t[0], t[1:]
    mass, middle, terminal, outer = v.factor_cut(t, 0)
    roots = [next(iter(S & outer)) for S in terminal]
    R = m.product(rest)
    coordinates = [cycle_coordinates(R, root, S) for root, S in zip(roots, terminal)]
    b, c = sorted(support(A) & middle)
    middle_coordinates = cycle_coordinates(R, b, middle)
    period = lcm(len(middle), *(len(S) for S in terminal))
    phases, forbidden = [], []
    for k in range(period):
        power = m.power(R, k)
        # m.conjugate(x,g) = x g x^-1; a positive full twist uses R^-1.
        developed = (A,) + tuple(m.conjugate(m.inverse(power), B) for B in rest)
        cut = v.factor_cut(developed, 0)
        assert cut[0] == mass and frozenset(cut[2]) == frozenset(terminal)
        for j, B in enumerate(rest, 1):
            common = support(developed[j]) & support(A)
            S_hits = [power[root] in support(B) for root in roots]
            M_hit = bool({power[b], power[c]} & support(B))
            assert [bool(common & S) for S in terminal] == S_hits
            assert bool(common & middle) == M_hit
            if any(S_hits) and len(common) == 1:
                forbidden.append((k, j))
        phases.append(developed)
    return {
        "mass": mass, "sizes": [*(len(S) for S in terminal), len(middle)],
        "period": period, "forbidden": forbidden, "phases": phases,
        "coordinates": coordinates, "middle_coordinates": middle_coordinates,
        "roots": roots, "terminal": terminal, "middle": middle, "R": R,
    }


def check_necessary_conditions(t, profile):
    """Conclusions valid when every cyclic development has no outer one-contact."""
    assert not profile["forbidden"]
    S, U = profile["terminal"]
    middle = profile["middle"]
    s, u, ell = profile["sizes"]
    x = y = z = internal = 0
    for B in t[1:]:
        counts = [len(support(B) & V) for V in (S, U, middle)]
        assert counts not in ([4, 0, 0], [0, 4, 0])
        if counts == [2, 2, 0]:
            x += 1
            assert s == u
            left, right = profile["coordinates"]
            assert {left[a] for a in support(B) & S} == {
                right[a] for a in support(B) & U}
        elif counts == [2, 0, 2]:
            y += 1
            assert ell <= 2 * s and ell // gcd(s, ell) <= 4
        elif counts == [0, 2, 2]:
            z += 1
            assert ell <= 2 * u and ell // gcd(u, ell) <= 4
        else:
            assert counts == [0, 0, 4]
            internal += 1
    assert x + y == s - 1 and x + z == u - 1
    assert y + z + 2 * internal == ell - 1
    if s != u and min(s, u) > 1:
        assert sorted((s, u, ell)) in ([2, 3, 4], [3, 4, 6])
        assert abs(s - u) == 1 and ell == 2 * min(s, u) and internal == 0
    if min(s, u) == 1 and max(s, u) > 1:
        larger = max(s, u)
        assert (ell == larger or (ell == 2 * larger and larger % 2 == 0)
                or (2 * ell == 3 * larger and larger % 4 == 0)
                or (3 * ell == 4 * larger and larger % 6 == 0))
    if y + z == 0 and profile["mass"] > 0:
        # Rooted cyclic coordinates give an isomorphism of the two colored trees.
        assert s == u and s > 1
        degree = h.degrees(t)
        leaf = next(a for a in S if a != profile["roots"][0] and degree[a] == 1)
        image = next(a for a in U if profile["coordinates"][1][a]
                     == profile["coordinates"][0][leaf])
        factor = next(j for j, B in enumerate(t) if leaf in support(B))
        assert degree[image] == 1 and image in support(t[factor])
        assert v.removable(t, factor)
    return x, y, z, internal


def twist_control(t):
    R = m.product(t[1:])
    word = list(range(2, len(t))) * (len(t) - 1)
    developed = m.braid(t, word)
    expected = (t[0],) + tuple(m.conjugate(m.inverse(R), B) for B in t[1:])
    assert developed == expected and m.product(developed) == m.product(t)
    return len(word)


def degree7_controls():
    # Reuse the retained complete domain. Do not extend the degree census.
    ds = h.dt_class(7)
    D = set(ds)
    P = m.perm(7, [tuple(range(7))])
    domain = {(a, b, c) for a in ds for b in ds
              if (c := m.compose(m.inverse(m.compose(a, b)), P)) in D}
    assert len(domain) == 392
    seen, orbit_rows = set(), []
    for seed in sorted(domain):
        if seed in seen:
            continue
        orbit, queue = {seed}, deque([seed])
        while queue:
            t = queue.popleft()
            for j in range(2):
                for sign in (-1, 1):
                    new = m.hurwitz(t, j, sign)
                    assert new in domain
                    if new not in orbit:
                        orbit.add(new)
                        queue.append(new)
        seen |= orbit
        minimum = min(v.mu(t) for t in orbit)
        group_order = h.group_closure(seed)
        assert group_order == (168 if minimum == 1 else 2520)
        positive_tests = phase_tests = front_tests = 0
        sizes = Counter()
        for t in sorted(orbit):
            if v.mu(t) != minimum:
                continue
            for at, value in enumerate(v.masses(t)):
                if value != minimum:
                    continue
                developed, _ = front(t, at)
                front_tests += 1
                assert v.mu(developed) == minimum
                twist_control(developed)
                if minimum == 0:
                    continue
                profile = phase_profile(developed)
                assert not profile["forbidden"]
                x, y, z, internal = check_necessary_conditions(developed, profile)
                assert y + z > 0
                assert 3 * minimum >= 7 - 4
                assert all(v.mu(u) == minimum for u in profile["phases"])
                positive_tests += 1
                phase_tests += profile["period"]
                sizes[tuple(sorted(profile["sizes"][:2])) + (profile["sizes"][2],)] += 1
        orbit_rows.append({"orbit_size": len(orbit), "group_order": group_order,
                           "minimum_mass": minimum,
                           "front_and_twist_tests": front_tests,
                           "positive_minimizing_marked_inputs": positive_tests,
                           "positive_phase_tests": phase_tests,
                           "positive_cut_size_distribution": {
                               str(k): v for k, v in sorted(sizes.items())}})
    assert seen == set(domain)
    assert sorted((row["orbit_size"], row["minimum_mass"]) for row in orbit_rows) == [
        (49, 1), (49, 1), (147, 0), (147, 0)]
    return orbit_rows


def plateau_control():
    data = json.loads((HERE / "terminal_mass_hurwitz_descent_certificate.json").read_text())
    changing = tuple(m.perm(15, es) for es in data["tuple_edges"])
    changing_at = v.masses(changing).index(v.mu(changing))
    changing_front, changing_word = front(changing, changing_at)
    assert changing_word == [-6, -5, -4, -3, -2, -1]
    assert v.factor_cut(changing_front, 0)[0] == v.mu(changing) == 1
    twist_control(changing_front)
    t = tuple(m.perm(15, es) for es in data["plateau_witness"]["tuple_edges"])
    at = v.masses(t).index(v.mu(t))
    t, front_word = front(t, at)
    profile = phase_profile(t)
    assert profile["mass"] == 1 and profile["forbidden"]
    k, factor = profile["forbidden"][0]
    active = profile["phases"][k]
    assert v.mu(active) <= 1
    if v.mu(active) == 1:
        assert (0, factor) in v.clean_contacts(active)
        v.prune(active, 0, factor)
    return {"initial_mass": 1, "front_word": front_word,
            "separate_changing_candidate_front_word": changing_word,
            "cut_sizes_terminal_terminal_middle": profile["sizes"],
            "complete_phase_period": profile["period"],
            "first_forbidden_phase": k, "factor_index_zero_based": factor,
            "lower_mass_at_phase_or_after_outer_pruning": True}


def synchronized_controls():
    rows = []
    for size in (2, 3, 5, 9):
        # Two parallel path trees attached at their roots to a five-sheet D5 middle.
        S, U = list(range(size)), list(range(size, 2 * size))
        n = 2 * size + 5
        M = list(range(2 * size, n))
        A = m.perm(n, [(S[0], M[0]), (U[0], M[4])])
        paired = tuple(m.perm(n, [(S[i], S[i + 1]), (U[i], U[i + 1])])
                       for i in range(size - 1))
        central = (m.perm(n, [(M[0], M[1]), (M[2], M[3])]),
                   m.perm(n, [(M[1], M[2]), (M[3], M[4])]))
        t = (A,) + paired + central
        assert m.cycle_lengths(m.product(t)) == [n]
        profile = phase_profile(t)
        assert not profile["forbidden"]
        assert check_necessary_conditions(t, profile) == (size - 1, 0, 0, 2)
        assert v.mu(t) == 0
        rows.append({"terminal_size": size, "degree": n, "minimum_mass": 0})
    return rows


def negative_controls():
    # Equal terminal sizes alone do not synchronize their colored edges.
    t = tuple(m.perm(11, es) for es in (
        ((0, 6), (3, 10)), ((0, 1), (4, 5)), ((1, 2), (3, 4)),
        ((6, 7), (8, 9)), ((7, 8), (9, 10))))
    assert m.cycle_lengths(m.product(t)) == [11]
    p = phase_profile(t)
    assert p["sizes"] == [3, 3, 5] and p["forbidden"]
    # An entire DT in a terminal component is detected by a later cyclic phase.
    internal = tuple(m.perm(11, es) for es in (
        ((0, 6), (5, 10)), ((0, 1), (2, 3)), ((1, 2), (3, 4)),
        ((6, 7), (8, 9)), ((7, 8), (9, 10))))
    q = phase_profile(internal)
    assert q["sizes"] == [5, 1, 5] and q["forbidden"]
    assert any(len(support(B) & q["terminal"][0]) == 4 for B in internal[1:])
    return {"unsynchronized_equal_terminal_trees_detected": True,
            "DT_wholly_in_a_terminal_component_detected": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "terminal_mass_global_compatibility_verification.json")
    args = parser.parse_args()
    degree7 = degree7_controls()
    plateau = plateau_control()
    synchronized = synchronized_controls()
    negative = negative_controls()
    inputs = ["TERMINAL_MASS_HURWITZ_DESCENT.md",
              "terminal_mass_hurwitz_descent_certificate.json",
              "verify_terminal_mass_hurwitz_descent.py",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt = {
        "schema_version": 1, "kind": "terminal_mass_global_compatibility_controls",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "degree7_complete_orbit_controls": degree7, "degree15_plateau_control": plateau,
        "synchronized_terminal_tree_controls": synchronized,
        "negative_controls": negative,
        "scope": "Exact controls of cut-preserving front placement, complementary full twists, cyclic contact covers, paired-terminal synchronization, the positive orbit-minimum bound 3m>=N-4 and the terminal-size restrictions. General statements have written proofs. No higher-degree census or exclusion of all positive full-alternating orbit minima.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "degree7_orbits": len(degree7),
                      "degree15_first_forbidden_phase": plateau["first_forbidden_phase"],
                      "synchronized_family_controls": len(synchronized),
                      "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
