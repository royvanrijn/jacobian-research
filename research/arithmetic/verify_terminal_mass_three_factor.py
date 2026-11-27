#!/usr/bin/env python3
"""Exact seven-weight B3 orbits and inequalities; not arbitrary DT extraction."""
from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import importlib.util
import json
from itertools import combinations
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("three-factor checker requires active Python assertions")

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "terminal_mass_helpers", HERE / "verify_terminal_mass_hurwitz_descent.py")
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
m, h = v.m, v.h
MOVES = (1, -1, 2, -2)


def canonical_factor(edges):
    answer = tuple(sorted(tuple(sorted(e)) for e in edges))
    assert len(answer) == 2 and all(len(e) == 2 and e[0] != e[1] for e in answer)
    assert len({a for e in answer for a in e}) == 4
    return answer


def edge_move(state, move):
    """Hurwitz by endpoint substitution, independently of permutation braid()."""
    assert move in MOVES
    at = abs(move) - 1
    A, B = state[at:at + 2]
    conjugator, target = (B, A) if move > 0 else (A, B)
    image = {a: b for a, b in conjugator} | {b: a for a, b in conjugator}
    changed = canonical_factor(tuple(image.get(a, a) for a in e) for e in target)
    pair = (B, changed) if move > 0 else (changed, A)
    return state[:at] + pair + state[at + 2:]


def permutations(state, support_size):
    return tuple(m.perm(support_size, factor) for factor in state)


def complete_orbit(initial):
    """No cutoff, depth limit or assumed transitivity: explore to closure."""
    initial = tuple(canonical_factor(g) for g in initial)
    support = {a for g in initial for e in g for a in e}
    assert support == set(range(len(support))) and len(support) <= 12
    initial_product = m.product(permutations(initial, len(support)))
    words, queue = {initial: ()}, deque([initial])
    while queue:
        state = queue.popleft()
        ps = permutations(state, len(support))
        assert m.product(ps) == initial_product
        for move in MOVES:
            new = edge_move(state, move)
            assert permutations(new, len(support)) == m.braid(ps, [move])
            assert edge_move(new, -move) == state
            if new not in words:
                words[new] = words[state] + (move,)
                queue.append(new)
    assert all(edge_move(state, move) in words for state in words for move in MOVES)
    return words


def validate_weights(owners, weights):
    assert len(weights) == 7 and set(owners) == set(range(7))
    assert all(isinstance(w, int) and w > 0 for w in weights)
    multiplicities = Counter(owners)
    assert all(weights[i] >= multiplicities[i] for i in range(7))


def terminal_masks(state, owners):
    """A bit i means w_i occurs in this factor's terminal mass, before -2."""
    es = [tuple(owners[a] for a in edge) for g in state for edge in g]
    assert all(a != b for a, b in es)
    assert len(set(map(frozenset, es))) == 6
    assert len(v.components(7, es)) == 1
    answer = []
    for g in state:
        removed = {frozenset(owners[a] for a in edge) for edge in g}
        comp = v.components(7, [e for e in es if frozenset(e) not in removed])
        assert len(comp) == 3
        middle, = (S for S in comp if all(S & e for e in removed))
        answer.append(sum(1 << i for i in range(7) if i not in middle))
    return tuple(answer)


def mass(mask, weights):
    return sum(w for i, w in enumerate(weights) if mask & (1 << i)) - 2


def indices(mask):
    return [i for i in range(7) if mask & (1 << i)]


def analyze(initial, owners, weights):
    validate_weights(owners, weights)
    words = complete_orbit(initial)
    witnesses = {}
    for state, word in words.items():
        for j, mask in enumerate(terminal_masks(state, owners)):
            witnesses.setdefault(mask, (word, j))
    minimal = sorted(c for c in witnesses
                     if not any(d != c and d & c == d for d in witnesses))
    assert all(any(d & c == d for d in minimal) for c in witnesses)
    best = min(mass(c, weights) for c in minimal)
    orbit_bytes = json.dumps(sorted(words), separators=(",", ":")).encode()
    return {
        "support_size": len(owners), "orbit_states": len(words),
        "directed_generator_edges": 4 * len(words),
        "orbit_sha256": hashlib.sha256(orbit_bytes).hexdigest(),
        "all_distinct_terminal_forms": len(witnesses),
        "all_terminal_coefficient_sets": [indices(c) for c in sorted(witnesses)],
        "irredundant_for_positive_weights": [
            {"components": indices(c), "word": list(witnesses[c][0]),
             "factor_index_zero_based": witnesses[c][1]} for c in minimal],
        "weights": list(weights), "endpoint_label_counts": [owners.count(i) for i in range(7)],
        "initial_factor_masses": [mass(c, weights) for c in terminal_masks(initial, owners)],
        "minimum_candidate_mass": best,
    }, words


def compress(t, at):
    local = t[at:at + 3]
    assert len(local) == 3
    removed = {tuple(sorted(e)) for g in local for e in h.edges(g)}
    components = v.components(len(t[0]), [e for e in v.edge_list(t) if e not in removed])
    assert len(components) == 7
    support = sorted({a for g in local for e in h.edges(g) for a in e})
    relabel = {a: i for i, a in enumerate(support)}
    initial = tuple(canonical_factor(tuple(relabel[a] for a in e) for e in h.edges(g))
                    for g in local)
    owners = tuple(next(i for i, S in enumerate(components) if a in S) for a in support)
    return initial, owners, tuple(map(len, components)), support, components


def window_control(t, at):
    initial, owners, weights, support, comp = compress(t, at)
    report, orbit = analyze(initial, owners, weights)
    comparisons = 0
    for state in orbit:
        local = tuple(m.perm(len(t[0]), [tuple(support[a] for a in e) for e in g])
                      for g in state)
        changed = t[:at] + local + t[at + 3:]
        assert m.product(changed) == m.product(t)
        for j, mask in enumerate(terminal_masks(state, owners)):
            assert mass(mask, weights) == v.factor_cut(changed, at + j)[0]
            comparisons += 1
    report.update({"window_start_zero_based": at, "full_tree_mass_comparisons": comparisons,
                   "support_labels": support, "owners": list(owners),
                   "components": [sorted(S) for S in comp]})
    return report


def symbolic_path_control():
    # The seven-vertex path has colors A B C B A C. Distinct labels in
    # components 2 and 4 MUST NOT be identified when computing the orbit.
    initial = (((0, 1), (6, 7)), ((1, 2), (4, 5)), ((3, 4), (7, 8)))
    owners = (0, 1, 2, 2, 3, 4, 4, 5, 6)
    weights = (4, 1, 4, 1, 2, 1, 2)
    report, orbit = analyze(initial, owners, weights)
    pairs = sorted([(0, 3), (0, 5), (0, 6), (1, 3),
                    (1, 5), (1, 6), (3, 5), (3, 6)])
    assert report["orbit_states"] == 162 and report["all_distinct_terminal_forms"] == 27
    assert sorted(tuple(row["components"]) for row in report["irredundant_for_positive_weights"]) == pairs
    # Every other mask strictly contains one of these pairs; positivity alone
    # proves redundancy for ALL weights, not only the numerical control above.
    assert report["minimum_candidate_mass"] == 0
    initial_masks = terminal_masks(initial, owners)
    domination = []
    for j, c in enumerate(initial_masks):
        d = next(sum(1 << i for i in pair) for pair in pairs
                 if sum(1 << i for i in pair) & c == sum(1 << i for i in pair)
                 and sum(1 << i for i in pair) != c)
        domination.append({"initial_factor_index_zero_based": j,
                           "initial_components": indices(c),
                           "strictly_smaller_components": indices(d),
                           "positive_drop_components": indices(c ^ d)})
    active = edge_move(edge_move(initial, -1), -2)
    assert terminal_masks(active, owners)[2] == (1 << 1) + (1 << 5)
    # Deliberately replacing endpoint labels by component vertices changes
    # the action. This is a negative control for the tempting seven-sheet model.
    collapsed = tuple(canonical_factor(tuple(owners[a] for a in e) for e in g) for g in initial)
    collapsed_orbit = complete_orbit(collapsed)
    assert len(collapsed_orbit) == 147 != len(orbit)
    report.update({"initial_edges_on_support": initial, "owners": owners,
                   "uniform_strict_domination_of_original_factors": domination,
                   "word_capturing_components_1_5": [-1, -2],
                   "collapsed_endpoint_orbit_states": len(collapsed_orbit),
                   "uniform_weight_domain": "integer w_i >= endpoint_label_counts[i]; no upper bound"})
    return report


def owner_partition_controls():
    """All seven-component owner partitions for this fixed nine-label forest."""
    initial = (((0, 1), (6, 7)), ((1, 2), (4, 5)), ((3, 4), (7, 8)))
    orbit = complete_orbit(initial)
    partitions = set()
    # Any partition of nine labels into seven nonempty blocks is generated
    # by two pair identifications: a triple block or two double blocks.
    for e, f in combinations(list(combinations(range(9), 2)), 2):
        blocks = v.components(9, [e, f])
        owners = tuple(next(i for i, S in enumerate(blocks) if a in S) for a in range(9))
        try:
            terminal_masks(initial, owners)
        except AssertionError:
            continue
        partitions.add(owners)
    assert len(partitions) == 189
    distribution, counts, rows = Counter(), Counter(), []
    feasible_marked_templates = excluded_marked_templates = 0
    for owners in sorted(partitions):
        forms = {c for state in orbit for c in terminal_masks(state, owners)}
        minimal = {c for c in forms if not any(d != c and d & c == d for d in forms)}
        initial_masks = terminal_masks(initial, owners)
        allowed = tuple(j for j, c in enumerate(initial_masks) if c in minimal)
        assert len(forms) == 27 and len(minimal) in (8, 12)
        distribution[len(allowed)] += 1
        counts[len(minimal)] += 1
        for j, c in enumerate(initial_masks):
            if j not in allowed:
                assert any(d != c and d & c == d for d in minimal)
                excluded_marked_templates += 1
                continue
            # A minimal mask can actually minimize the ISOLATED linear
            # system at positive mass. This does not realize a global tuple.
            lower = tuple(owners.count(i) + 1 for i in range(7))
            large = 1 + sum(w for i, w in enumerate(lower) if c & (1 << i))
            weights = [lower[i] if c & (1 << i) else large for i in range(7)]
            outside = next(i for i in range(7) if not c & (1 << i))
            if sum(weights) % 2 == 0:
                weights[outside] += 1
            validate_weights(owners, weights)
            candidate = mass(c, weights)
            assert candidate > 0 and sum(weights) % 2 == 1
            assert all(mass(d, weights) >= candidate for d in forms)
            feasible_marked_templates += 1
        rows.append({"owners": list(owners), "initial_coefficient_bitmasks": list(initial_masks),
                     "minimal_coefficient_bitmasks": sorted(minimal),
                     "possible_initial_minimizer_indices": list(allowed)})
    assert distribution == {0: 44, 1: 111, 2: 33, 3: 1}
    assert counts == {8: 162, 12: 27}
    assert excluded_marked_templates == 387 and feasible_marked_templates == 180
    return {"scope": "All owner partitions for the fixed three-disjoint-two-edge-path support forest; not all triple forests",
            "owner_partitions": len(partitions), "common_local_orbit_states": len(orbit),
            "coefficient_bitmask_convention": "Bit i is the coefficient of w_i; every form has constant -2",
            "minimal_form_count_distribution": dict(sorted(counts.items())),
            "number_of_possible_initial_minimizers_distribution": dict(sorted(distribution.items())),
            "marked_templates_excluded_uniformly": excluded_marked_templates,
            "marked_templates_with_positive_isolated_linear_witness": feasible_marked_templates,
            "all_weight_witnesses_checked_with_odd_total": True, "templates": rows}


def controls():
    data = json.loads((HERE / "terminal_mass_hurwitz_descent_certificate.json").read_text())
    plateau = data["plateau_witness"]
    t = tuple(m.perm(15, es) for es in plateau["tuple_edges"])
    windows = [window_control(t, at) for at in range(5)]
    assert [row["orbit_states"] for row in windows] == [18, 144, 162, 18, 18]
    assert [row["minimum_candidate_mass"] for row in windows] == [1, 0, 0, 1, 1]
    assert sum(row["full_tree_mass_comparisons"] for row in windows) == 1080
    proper = h.cyclic(3)
    proper_report = window_control(proper, 0)
    assert proper_report["minimum_candidate_mass"] == 1
    assert h.group_closure(proper) == 168
    disjoint = tuple(m.perm(13, es) for es in (
        ((0, 1), (2, 3)), ((4, 5), (6, 7)), ((8, 9), (10, 11)),
        ((1, 2), (5, 6)), ((3, 4), (9, 10)), ((7, 8), (11, 12))))
    disjoint_report = window_control(disjoint, 0)
    assert disjoint_report["support_size"] == 12 and disjoint_report["orbit_states"] == 6
    return windows, proper_report, disjoint_report


def invalid_controls():
    def rejected(action):
        try:
            action()
        except (AssertionError, StopIteration):
            return True
        raise AssertionError("Invalid input accepted")
    assert rejected(lambda: canonical_factor(((0, 1), (1, 2))))
    assert rejected(lambda: validate_weights((0, 1, 2, 2, 3, 4, 4, 5, 6), (1,) * 7))
    state = (((0, 1), (6, 7)), ((1, 2), (4, 5)), ((3, 4), (7, 8)))
    assert rejected(lambda: terminal_masks(state, (0,) * 9))
    orbit = complete_orbit(state)
    corrupted = set(orbit) - {edge_move(state, 1)}
    assert not all(edge_move(u, k) in corrupted for u in corrupted for k in MOVES)
    return {"non_DT_rejected": True, "insufficient_component_weights_rejected": True,
            "quotient_loops_rejected": True, "incomplete_orbit_detected": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "terminal_mass_three_factor_verification.json")
    args = parser.parse_args()
    symbolic = symbolic_path_control()
    partitions = owner_partition_controls()
    windows, proper, disjoint = controls()
    invalid = invalid_controls()
    inputs = ["TERMINAL_MASS_HURWITZ_DESCENT.md",
              "terminal_mass_hurwitz_descent_certificate.json",
              "verify_terminal_mass_hurwitz_descent.py",
              "verify_paired_hurwitz_leaf_extraction.py",
              "verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt = {
        "schema_version": 1, "kind": "seven_weight_three_factor_orbits", "result": "PASS",
        "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "symbolic_nine_label_path_template": symbolic,
        "complete_owner_inventory_for_the_nine_label_support_forest": partitions,
        "degree15_all_initial_three_factor_windows": windows,
        "proper_order168_positive_local_minimum_control": proper,
        "twelve_distinct_support_labels_control": disjoint,
        "invalid_controls": invalid,
        "scope": "Complete local orbits and exact coefficient sets for the stated templates. The path template gives eight inequalities and uniform strict domination of each original factor. All 189 owner partitions for its nine-label support forest exclude 387 of 567 marked templates uniformly, leaving 180 with positive isolated linear witnesses. Seven-weight compression and global-minimum necessity have written proofs. Other support forests and exclusion of every positive global orbit minimum remain open.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "symbolic_orbit_states": symbolic["orbit_states"],
                      "symbolic_forms": symbolic["all_distinct_terminal_forms"],
                      "symbolic_minimal_forms": len(symbolic["irredundant_for_positive_weights"]),
                      "complete_nine_label_forest_owner_partitions": partitions["owner_partitions"],
                      "uniformly_excluded_marked_templates": partitions["marked_templates_excluded_uniformly"],
                      "positive_isolated_linear_marked_templates": partitions["marked_templates_with_positive_isolated_linear_witness"],
                      "window_orbit_states": [r["orbit_states"] for r in windows],
                      "full_tree_comparisons": sum(r["full_tree_mass_comparisons"]
                                                   for r in windows + [proper, disjoint]),
                      "proper_control_minimum": proper["minimum_candidate_mass"],
                      "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
