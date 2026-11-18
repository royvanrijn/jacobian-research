#!/usr/bin/env python3
"""Exact controls for the conditional pure leaf-pair extraction."""
from __future__ import annotations
import argparse
from collections import deque
from itertools import combinations
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
HELPER = HERE / "verify_unmarked_polynomial_hurwitz_reduction.py"
spec = importlib.util.spec_from_file_location("marking", HELPER)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
W = [-1, -2, -2, -1]

def cyclic(r):
    n = 2 * r + 1
    P = m.perm(n, [tuple(range(n))])
    d = m.perm(n, [(0, 2), (3, 4)])
    return tuple(m.conjugate(m.power(P, 2 * i), d) for i in range(r))

def edges(p):
    done, answer = set(), []
    for a in range(len(p)):
        if a in done:
            continue
        cycle, j = [], a
        while j not in done:
            done.add(j)
            cycle.append(j)
            j = p[j]
        if len(cycle) == 2:
            answer.append(tuple(cycle))
        elif len(cycle) == 3:
            answer.extend([(cycle[0], cycle[1]), (cycle[1], cycle[2])])
        else:
            assert len(cycle) == 1
    return answer

def degrees(t):
    degree = [0] * len(t[0])
    for p in t:
        for a, b in edges(p):
            degree[a] += 1
            degree[b] += 1
    return degree

def delete(t, factor_index, leaves):
    leaves = set(leaves)
    assert len(leaves) == 2
    for i, p in enumerate(t):
        if i != factor_index:
            assert all(p[a] == a for a in leaves)
    old = [a for a in range(len(t[0])) if a not in leaves]
    remap = {a: i for i, a in enumerate(old)}
    smaller = tuple(tuple(remap[p[a]] for a in old)
                    for i, p in enumerate(t) if i != factor_index)
    assert m.cycle_lengths(m.product(smaller)) == [len(old)]
    return smaller

def local_identity():
    n = 9
    t = tuple(m.perm(n, cs) for cs in [
        [(0, 2), (3, 4)], [(2, 4), (5, 6)], [(4, 6), (7, 8)]])
    expected = tuple(m.perm(n, cs) for cs in [
        [(2, 3), (4, 5)], [(0, 3), (5, 6)], [(3, 6), (7, 8)]])
    assert m.braid(t, W) == expected
    assert m.product(t) == m.product(expected)
    assert m.labels(3, W) == (0, 1, 2)
    return {"word": W, "strand_permutation": [0, 1, 2],
            "exact_triple_identity": True}

def one_cyclic_step(r):
    t = cyclic(r)
    moved = m.braid(t, W)
    assert moved[3:] == t[3:]
    degree = degrees(moved)
    assert degree[2] == degree[4] == 1
    smaller = delete(moved, 0, [2, 4])
    assert smaller == cyclic(r - 1)
    return smaller

def recursive_controls(r):
    count, active = 0, r
    while active > 4:
        one_cyclic_step(active)
        count += 1
        active -= 1
    assert count == r - 4
    return {"initial_r": r, "deletions_to_A9_core": count,
            "local_moves_on_successive_cores": 4 * count,
            "terminal_core_r": active, "terminal_degree": 9}

def external_controls():
    base = tuple(p + (11, 12) for p in cyclic(5))
    count = 0
    allowed = [a for a in range(11) if a not in (2, 4)]
    for a in allowed:
        for b in allowed:
            if a == b:
                continue
            t = base + (m.perm(13, [(a, 11), (b, 12)]),)
            assert m.cycle_lengths(m.product(t)) == [13]
            moved = m.braid(t, W)
            assert degrees(moved)[2] == degrees(moved)[4] == 1
            delete(moved, 0, [2, 4])
            count += 1
    bad = base + (m.perm(13, [(2, 11), (0, 12)]),)
    assert m.cycle_lengths(m.product(bad)) == [13]
    bad_moved = m.braid(bad, W)
    degree = degrees(bad_moved)
    assert [degree[a] for a in (2, 3, 4, 5)] == [2, 3, 1, 2]
    assert not all(degree[a] == 1 or degree[b] == 1
                   for a, b in edges(bad_moved[0]))
    # T remains the first marked branch while a later DT window is extracted.
    marked = (m.perm(13, [(1, 11, 12)]),) + base
    shifted = [-2, -3, -3, -2]
    assert m.labels(6, shifted) == tuple(range(6))
    moved = m.braid(marked, shifted)
    assert moved[0] == marked[0]
    assert m.cycle_lengths(m.product(marked)) == [13]
    delete(moved, 1, [2, 4])
    return {"all_allowed_ordered_anchor_pairs": count,
            "missing_external_support_condition_rejected": True,
            "counterexample_degrees_2_3_4_5": [2, 3, 1, 2],
            "marked_T_branch_retained": True}

def evacuation_word(moved_count):
    word = []
    for j in range(1, moved_count + 1):
        word.extend((j + 2, j + 1))
    word.extend(-j for j in range(1, moved_count + 1))
    return word

def evacuate(after_identity, moved_count):
    E, U1, U2 = after_identity[:3]
    etas = after_identity[3:3 + moved_count]
    rest = after_identity[3 + moved_count:]
    T = m.product(etas)
    output = m.braid(after_identity, evacuation_word(moved_count))
    expected = tuple(m.conjugate(E, eta) for eta in etas)
    expected += (E, m.conjugate(m.inverse(T), U1),
                    m.conjugate(m.inverse(T), U2))
    expected += rest
    assert output == expected
    assert m.product(output) == m.product(after_identity)
    assert len(evacuation_word(moved_count)) == 3 * moved_count
    F = {a for a in range(len(E)) if U1[a] == U2[a] == a}
    hypotheses = (all(eta[3] == 3 and eta[5] == 5 for eta in etas)
                  and T[2] in F and T[4] in F
                  and all(g[2] == 2 and g[4] == 4 for g in rest))
    if hypotheses:
        assert degrees(output)[2] == degrees(output)[4] == 1
        delete(output, moved_count, [2, 4])
    return output, hypotheses

def evacuation_controls():
    controls = []
    for count in (1, 2, 5, 20):
        n = 11 + 2 * count
        base = tuple(p + tuple(range(11, n)) for p in cyclic(5))
        etas = tuple(m.perm(n, [(2, 11 + 2 * j), (4, 12 + 2 * j)])
                     for j in range(count))
        original = base[:3] + etas + base[3:]
        assert m.cycle_lengths(m.product(original)) == [n]
        moved = m.braid(original, W)
        assert degrees(moved)[2] > 1 and degrees(moved)[4] > 1
        output, hypotheses = evacuate(moved, count)
        assert hypotheses
        controls.append({"external_factors": count, "degree": n,
                         "repair_moves": 3 * count,
                         "two_exclusive_leaf_ends_after_repair": True})
    base = tuple(p + (11, 12) for p in cyclic(5))
    original_bad = base + (m.perm(13, [(2, 11), (0, 12)]),)
    output, hypotheses = evacuate(m.braid(original_bad, W), 3)
    assert hypotheses
    # Paired pivot/anchor interaction violates only the anchor condition.
    paired = base[:3] + (m.perm(13, [(2, 11), (3, 12)]),) + base[3:]
    assert m.cycle_lengths(m.product(paired)) == [13]
    moved = m.braid(paired, W)
    E, U1, U2 = moved[:3]
    eta = moved[3]
    F = {a for a in range(13) if U1[a] == U2[a] == a}
    assert eta[2] in F and eta[4] in F and eta[3] != 3
    output, hypotheses = evacuate(moved, 1)
    assert not hypotheses
    assert output[0] == m.perm(13, [(3, 11), (2, 12)])
    deg = degrees(output)
    assert deg[2] > 1 and deg[3] > 1
    assert not all(deg[a] == 1 or deg[b] == 1 for a, b in edges(output[1]))
    return {"unbounded_interaction_family_controls": controls,
            "original_degree13_counterexample_repaired_with_block_length": 3,
            "paired_pivot_anchor_blocker_rejected": True,
            "noncommutative_product_conjugation_formula_checked": True}

def group_closure(generators):
    ident = tuple(range(len(generators[0])))
    seen, queue = {ident}, [ident]
    for g in queue:
        for s in generators:
            h = m.compose(g, s)
            if h not in seen:
                seen.add(h)
                queue.append(h)
    assert all(m.compose(g, s) in seen for g in seen for s in generators)
    return len(seen)

def support(p):
    return sum(1 << i for i, j in enumerate(p) if i != j)

def code_span(vectors):
    C = {0}
    for v in vectors:
        C |= {x ^ v for x in C}
    return frozenset(C)

def orthogonal(vectors):
    return all((x & y).bit_count() % 2 == 0 for x in vectors for y in vectors)

def dt_class(n):
    answer = set()
    for a, b, c, d in combinations(range(n), 4):
        for pairs in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            answer.add(m.perm(n, pairs))
    return sorted(answer)

def code_controls():
    ds = dt_class(7)
    checks = 0
    for a in ds:
        sa = support(a)
        for b in ds:
            sb = support(b)
            if (sa & sb).bit_count() % 2:
                continue
            changed = support(m.conjugate(a, b))
            assert changed in (sb, sa ^ sb)
            assert code_span([sa, sb]) == code_span([sa, changed])
            checks += 1
    proper = cyclic(3)
    C = code_span(map(support, proper))
    assert len(C) == 8 and sorted(x.bit_count() for x in C) == [0] + [4] * 7
    assert orthogonal(list(map(support, proper)))
    full = tuple(m.perm(7, pairs) for pairs in [
        [(0, 1), (3, 4)], [(1, 2), (4, 5)], [(2, 3), (5, 6)]])
    assert group_closure(full) == 2520
    assert not orthogonal(list(map(support, full)))
    return {"all_DT_elements_degree7": len(ds),
            "all_even_intersection_pair_updates_checked": checks,
            "proper_G3_code_dimension": 3,
            "proper_G3_code_weight_distribution": {"0": 1, "4": 7},
            "full_A7_seed_has_odd_overlap": True}

def bring_adjacent(t, i, j):
    if i > j:
        i, j = j, i
    left, right = t[i], t[j]
    word = list(range(j, i + 1, -1))
    moved = m.braid(t, word)
    assert moved[i] == left and moved[i + 1] == right
    return moved, word, i

def odd_overlap_gate(t):
    initial_product = m.product(t)
    chosen = next(((i, j) for i in range(len(t)) for j in range(i + 1, len(t))
                   if (support(t[i]) & support(t[j])).bit_count() % 2), None)
    assert chosen is not None, "The full-A_N/odd-overlap input hypothesis is absent"
    t, word, at = bring_adjacent(t, *chosen)
    overlap = (support(t[at]) & support(t[at + 1])).bit_count()
    rotation = 0
    if overlap == 3:
        core = {i for i in range(len(t[0]))
                if t[at][i] != i or t[at + 1][i] != i}
        assert len(core) == 5
        gi = next(i for i, g in enumerate(t) if i not in (at, at + 1)
                  and any(g[a] != a for a in core))
        contact = {a for a in core if t[gi][a] != a}
        assert 1 <= len(contact) <= 2
        if len(contact) == 1:
            first = next(i for i in (at, at + 1)
                         if (support(t[i]) & support(t[gi])).bit_count() == 1)
        else:
            before = t
            for k in range(5):
                trial = m.braid(before, [at + 1] * (2 * k))
                if (support(trial[at]) & support(trial[gi])).bit_count() == 1:
                    t = trial
                    local = [at + 1] * (2 * k)
                    rotation = len(local)
                    word += local
                    break
            else:
                raise AssertionError("Five-cycle conjugation did not sweep fixed vertices")
            first = at
        t, extra, at = bring_adjacent(t, first, gi)
        word += extra
    assert (support(t[at]) & support(t[at + 1])).bit_count() == 1
    assert m.product(t) == initial_product
    return t, word, at, rotation

def secondary_exchange(t, at):
    A, B = t[at:at + 2]
    common = support(A) & support(B)
    assert common.bit_count() == 1
    v = common.bit_length() - 1
    ae, be = edges(A), edges(B)
    primary_a = next(e for e in ae if v in e)
    primary_b = next(e for e in be if v in e)
    secondary_a = next(e for e in ae if v not in e)
    secondary_b = next(e for e in be if v not in e)
    n = len(A)
    expected = (m.perm(n, [primary_a, secondary_b]),
                m.perm(n, [primary_b, secondary_a]))
    result = m.braid(t, [at + 1] * 3)
    assert result[at:at + 2] == expected
    assert result[:at] == t[:at] and result[at + 2:] == t[at + 2:]
    old_edges = sorted(tuple(sorted(e)) for p in t for e in edges(p))
    new_edges = sorted(tuple(sorted(e)) for p in result for e in edges(p))
    assert old_edges == new_edges and m.product(result) == m.product(t)

def complete_degree7_gate():
    ds = dt_class(7)
    D = set(ds)
    P = m.perm(7, [tuple(range(7))])
    all_tuples = set()
    for a in ds:
        for b in ds:
            c = m.compose(m.inverse(m.compose(a, b)), P)
            if c in D:
                all_tuples.add((a, b, c))
    assert len(all_tuples) == 392
    unseen, orbit_data = set(all_tuples), []
    full_count, rotation_count, max_word = 0, 0, 0
    while unseen:
        seed = min(unseen)
        orbit, queue = {seed}, deque([seed])
        while queue:
            t = queue.popleft()
            for j in range(2):
                for sign in (-1, 1):
                    u = m.hurwitz(t, j, sign)
                    assert u in all_tuples
                    if u not in orbit:
                        orbit.add(u)
                        queue.append(u)
        unseen -= orbit
        size = group_closure(seed)
        if size == 168:
            C = code_span(map(support, seed))
            assert all(orthogonal(list(map(support, t)))
                       and code_span(map(support, t)) == C for t in orbit)
        else:
            assert size == 2520
            for t in orbit:
                final, word, at, rotation = odd_overlap_gate(t)
                assert m.braid(t, word) == final
                secondary_exchange(final, at)
                full_count += 1
                rotation_count += bool(rotation)
                max_word = max(max_word, len(word))
        orbit_data.append({"tuple_count": len(orbit), "group_order": size})
    assert sorted((x["tuple_count"], x["group_order"]) for x in orbit_data) == [
        (49, 168), (49, 168), (147, 2520), (147, 2520)]
    assert full_count == 294
    try:
        odd_overlap_gate(cyclic(3))
    except AssertionError:
        pass
    else:
        raise AssertionError("Proper even-overlap core was incorrectly accepted")
    return {"complete_fixed_product_domain": 392,
            "orbits": sorted(orbit_data, key=lambda x: (x["group_order"], x["tuple_count"])),
            "all_full_A7_inputs_gate_and_exchange_checked": full_count,
            "inputs_using_five_cycle_rotation": rotation_count,
            "maximum_gate_word_length_in_control": max_word,
            "proper_even_overlap_core_rejected": True}

def larger_gate_controls():
    answer = []
    for r in (4, 5, 6, 12, 25, 50):
        t = cyclic(r)
        final, word, at, rotation = odd_overlap_gate(t)
        assert m.braid(t, word) == final
        secondary_exchange(final, at)
        answer.append({"r": r, "degree": 2 * r + 1,
                       "gate_word_length": len(word),
                       "adjacent_overlap": 1, "secondary_exchange_checked": True})
    return answer

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path,
        default=HERE / "paired_hurwitz_leaf_extraction_verification.json")
    args = ap.parse_args()
    identity = local_identity()
    one_cyclic_step(4)
    controls = [recursive_controls(r) for r in (4, 5, 6, 7, 12, 25, 50)]
    external = external_controls()
    evacuation = evacuation_controls()
    code = code_controls()
    gate = complete_degree7_gate()
    larger_gate = larger_gate_controls()
    size = group_closure(cyclic(3))
    assert size == 168 and size < factorial(7) // 2
    paths = ["PAIRED_HURWITZ_LEAF_EXTRACTION.md",
             "UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md",
             "verify_unmarked_polynomial_hurwitz_reduction.py",
             "ALTERNATING_HURWITZ_CYCLE_FOLDING.md"]
    receipt = {
        "schema_version": 1, "kind": "paired_tree_extraction_and_overlap_gates",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in paths},
        "local_identity": identity, "cyclic_recursion_controls": controls,
        "external_support_controls": external,
        "block_evacuation_controls": evacuation,
        "binary_code_controls": code,
        "complete_degree7_overlap_gate": gate,
        "higher_degree_overlap_gate_controls": larger_gate,
        "complete_degree7_core_group_order": size,
        "full_A7_order_for_comparison": factorial(7) // 2,
        "scope": "Exact pure extraction, block-evacuation repair, cyclic recursion, even-overlap invariant and arbitrary-input overlap-one gate. No arbitrary-tree leaf extraction, extension connectivity or all-degree orbit classification is proved.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))

if __name__ == "__main__":
    main()

