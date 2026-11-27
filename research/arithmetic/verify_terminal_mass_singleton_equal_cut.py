#!/usr/bin/env python3
"""Bounded exact controls of Lemma 15; the universal proof is in the note."""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
import heapq
import importlib.util
from itertools import permutations, product
import json
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("singleton-cut controls require active assertions")

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "terminal_mass_helpers", HERE / "verify_terminal_mass_hurwitz_descent.py")
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
m, h = v.m, v.h


def trees(labels):
    """All labeled trees, once each, by the Prüfer bijection."""
    labels = tuple(labels)
    for code in product(labels, repeat=len(labels) - 2):
        degree = dict.fromkeys(labels, 1)
        for a in code:
            degree[a] += 1
        leaves = [a for a in labels if degree[a] == 1]
        heapq.heapify(leaves)
        edges = []
        for a in code:
            leaf = heapq.heappop(leaves)
            edges.append(tuple(sorted((leaf, a))))
            degree[a] -= 1
            if degree[a] == 1:
                heapq.heappush(leaves, a)
        edges.append(tuple(sorted(leaves)))
        yield tuple(sorted(edges))


def root_side(edges, root, removed):
    """Separate graph traversal for the side sizes in inequality (35)."""
    reached, queue = {root}, [root]
    for a in queue:
        for x, y in edges:
            if (x, y) == removed:
                continue
            if x == a and y not in reached:
                reached.add(y)
                queue.append(y)
            if y == a and x not in reached:
                reached.add(x)
                queue.append(x)
    return reached


def path_from(edges, root):
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    assert len(adjacency[root]) == 1
    assert all(len(neighbors) <= 2 for neighbors in adjacency.values())
    path = [root]
    while len(path) < len(adjacency):
        new = [a for a in adjacency[path[-1]] if a not in path]
        assert len(new) == 1
        path.append(new[0])
    return path


def path_tuple(size, order=None, pairing=None):
    """a=0, M=(1,...,t), U=(t+1,...,2t)."""
    n = 2 * size + 1
    if order is None:
        order = range(size - 1)
    if pairing is None:
        pairing = list(range(size - 1))
    A = m.perm(n, [(0, 1), (size, size + 1)])
    Bs = [m.perm(n, [(i + 1, i + 2),
                     (size + 1 + pairing[i], size + 2 + pairing[i])])
          for i in range(size - 1)]
    return (A,) + tuple(Bs[i] for i in order)


def capture(t, middle_path, terminal_path):
    """Find a closest neighboring-index pair, commute it together, then capture."""
    size = len(middle_path)
    assert size >= 3 and len(terminal_path) == size
    Bs = [m.perm(len(t[0]), [(middle_path[i], middle_path[i + 1]),
                            (terminal_path[i], terminal_path[i + 1])])
          for i in range(size - 1)]
    assert set(t[1:]) == set(Bs) and len(set(Bs)) == size - 1
    positions = {B: t.index(B) for B in Bs}
    distance, i = min((abs(positions[Bs[i]] - positions[Bs[i + 1]]), i)
                      for i in range(size - 2))
    left, right = sorted((positions[Bs[i]], positions[Bs[i + 1]]))
    C, D = t[left], t[right]
    for E in t[left + 1:right]:
        assert m.compose(C, E) == m.compose(E, C)
        assert m.compose(D, E) == m.compose(E, D)
    shared = {middle_path[i + 1], terminal_path[i + 1]}
    degree = h.degrees(t)
    assert all(degree[a] == 2 for a in shared)
    assert all(t[0][a] == a for a in shared)
    active, word = t, []
    for at in range(left, right - 1):
        # Positive moves on commuting factors are literal swaps.
        assert active[at] == C
        crossed = active[at + 1]
        active = m.hurwitz(active, at, 1)
        assert active[at:at + 2] == (crossed, C)
        assert set(active) == set(t)
        word.append(at + 1)
    at = right - 1
    assert active[at:at + 2] == (C, D)
    # Both signs must create a removable factor; check via supports, not cuts.
    for sign, retained_at in ((1, at), (-1, at + 1)):
        endpoint = m.hurwitz(active, at, sign)
        assert v.removable(endpoint, retained_at)
        assert all(h.degrees(endpoint)[a] == 1 for a in shared)
        assert v.mu(endpoint) == 0
        assert m.product(endpoint) == m.product(t)
        keep = [a for a in range(len(t[0])) if a not in shared]
        relabel = {a: j for j, a in enumerate(keep)}
        lower = tuple(tuple(relabel[B[a]] for a in keep)
                      for j, B in enumerate(endpoint) if j != retained_at)
        assert m.cycle_lengths(m.product(lower)) == [len(keep)]
        v.edge_list(lower)
    word.append(at + 1)
    endpoint = m.braid(t, word)
    assert v.mu(endpoint) == 0
    assert len(word) == distance
    return word


def exhaustive_tree_controls():
    rows = []
    for size in range(2, 5):
        n = 2 * size + 1
        A = m.perm(n, [(0, 1), (size, size + 1)])
        middle_trees = list(trees(range(1, size + 1)))
        terminal_trees = list(trees(range(size + 1, 2 * size + 1)))
        assert len(set(middle_trees)) == len(set(terminal_trees)) == size ** (size - 2)
        tested = survivors = captured = 0
        for middle in middle_trees:
            for terminal in terminal_trees:
                for pairing in permutations(terminal):
                    t = (A,) + tuple(m.perm(n, [e, f]) for e, f in zip(middle, pairing))
                    assert m.cycle_lengths(m.product(t)) == [n]
                    assert v.factor_cut(t, 0)[0] == size - 1
                    middle_weights = []
                    for at, (e, f) in enumerate(zip(middle, pairing), 1):
                        Mc = root_side(middle, size, e)
                        Ud = root_side(terminal, size + 1, f)
                        weight = len(Mc) + len(Ud) + (1 in Mc)
                        assert len(v.factor_cut(t, at)[1]) == weight
                        middle_weights.append(weight)
                    tested += 1
                    if any(weight > size for weight in middle_weights):
                        assert v.mu(t) < size - 1
                        continue
                    survivors += 1
                    assert v.mu(t) == size - 1
                    M = path_from(middle, 1)
                    U = path_from(terminal, size + 1)
                    assert M[-1] == size
                    for e, f in zip(middle, pairing):
                        mi = next(i for i in range(size - 1) if set(e) == set(M[i:i + 2]))
                        ui = next(i for i in range(size - 1) if set(f) == set(U[i:i + 2]))
                        assert mi == ui
                    if size >= 3:
                        capture(t, M, U)
                        captured += 1
        rows.append({"t": size, "tree_pairings_tested": tested,
                     "all_factor_mass_survivors": survivors,
                     "survivors_captured_in_both_signs": captured})
    assert [row["tree_pairings_tested"] for row in rows] == [1, 18, 1536]
    assert [row["all_factor_mass_survivors"] for row in rows] == [1, 2, 12]
    return rows


def path_pairing_controls():
    rows = []
    for size in range(3, 9):
        tested = survivors = 0
        for pairing in permutations(range(size - 1)):
            t = path_tuple(size, pairing=pairing)
            weights = [len(v.factor_cut(t, at + 1)[1]) for at in range(size - 1)]
            assert weights == [size + pairing[i] - i for i in range(size - 1)]
            tested += 1
            if max(weights) <= size:
                survivors += 1
                assert pairing == tuple(range(size - 1))
        assert survivors == 1
        rows.append({"t": size, "pairings_tested": tested, "survivors": survivors})
    return rows


def order_controls():
    rows = []
    for size in range(3, 9):
        tested = longest = 0
        for order in permutations(range(size - 1)):
            t = path_tuple(size, order)
            word = capture(t, list(range(1, size + 1)), list(range(size + 1, 2 * size + 1)))
            tested += 1
            longest = max(longest, len(word))
        rows.append({"t": size, "orders_tested": tested, "longest_word": longest})
    larger = []
    for size in (16, 32):
        indices = list(range(size - 1))
        orders = [indices[::-1], indices[::2] + indices[1::2]]
        for order in orders:
            t = path_tuple(size, order)
            word = capture(t, list(range(1, size + 1)), list(range(size + 1, 2 * size + 1)))
            larger.append({"t": size, "order": order, "word": word})
    return rows, larger


def boundary_and_negative_controls():
    seed = path_tuple(2)
    orbit, queue = {seed}, deque([seed])
    while queue:
        active = queue.popleft()
        for sign in (-1, 1):
            endpoint = m.hurwitz(active, 0, sign)
            if endpoint not in orbit:
                orbit.add(endpoint)
                queue.append(endpoint)
    assert len(orbit) == 5 and all(v.mu(t) == 1 for t in orbit)
    assert h.group_closure(seed) == 10
    rejected = 0
    # Incorrect edge pairing violates an actual factor inequality.
    wrong_pairing = path_tuple(3, pairing=[1, 0])
    assert v.mu(wrong_pairing) < 2
    try:
        capture(wrong_pairing, [1, 2, 3], [4, 5, 6])
    except AssertionError:
        rejected += 1
    # The two-contact move requires the claimed degree-two contacts.
    wrong_paths = path_tuple(4)
    try:
        capture(wrong_paths, [1, 3, 2, 4], [5, 6, 7, 8])
    except AssertionError:
        rejected += 1
    # One-cut covers without all-factor minimality do not force (35).
    size = 7
    one_cut_only = (m.perm(15, [(0, 1), (2, 8)]),) + tuple(
        m.perm(15, [(i + 1, i + 2), (i + 8, i + 9)]) for i in range(size - 1))
    assert v.factor_cut(one_cut_only, 0)[0] == 6 and v.mu(one_cut_only) == 0
    try:
        capture(seed, [1, 2], [3, 4])
    except AssertionError:
        rejected += 1
    assert rejected == 3
    return {"t2_complete_orbit_size": len(orbit), "t2_orbit_minimum_mass": 1,
            "t2_generated_group_order": 10, "invalid_capture_inputs_rejected": rejected,
            "cut_mass_without_all_factor_minimality": 6,
            "same_tuple_global_mass": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "terminal_mass_singleton_equal_cut_verification.json")
    args = parser.parse_args()
    tree_rows = exhaustive_tree_controls()
    pair_rows = path_pairing_controls()
    order_rows, larger = order_controls()
    boundary = boundary_and_negative_controls()
    inputs = ["TERMINAL_MASS_HURWITZ_DESCENT.md", "verify_terminal_mass_hurwitz_descent.py",
              "verify_paired_hurwitz_leaf_extraction.py", "verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt = {
        "schema_version": 1, "kind": "terminal_mass_singleton_equal_cut_controls",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
        "all_labeled_tree_controls": tree_rows, "all_path_pairing_controls": pair_rows,
        "all_complementary_order_controls": order_rows, "larger_order_controls": larger,
        "boundary_and_negative_controls": boundary,
        "scope": "Internal bounded exact controls of Lemma 15. All labeled tree pairs and edge pairings for t=2..4, all path edge pairings and complementary orders for t=3..8, four larger literal orders, both signs of final capture, product preservation and lower full-cycle deletion. The unbounded exclusion follows from the written proof, not these finite enumerations. No all-degree extraction or mixed connectivity.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "PASS", "tree_pairings": sum(r["tree_pairings_tested"] for r in tree_rows),
                      "complementary_orders": sum(r["orders_tested"] for r in order_rows),
                      "t2_orbit_minimum": boundary["t2_orbit_minimum_mass"],
                      "receipt": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
