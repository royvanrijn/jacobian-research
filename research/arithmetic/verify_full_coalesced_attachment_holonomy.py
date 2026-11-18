#!/usr/bin/env python3
"""Exact stabilizer words and the complete degree-seven V4 Schreier bound.

Uniform maximality in higher degree is proved in the companion note by
transpositions/three-cycles; the bounded controls do not substitute for that proof.
"""
from __future__ import annotations
import argparse
from collections import deque
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import verify_uniform_hurwitz_attachment as core

HERE = Path(__file__).resolve().parent
KINDS = ("double_transposition", "three_cycle")
DT_WORDS = ([2, 1, 1, 1, -2],
            [1, 2, 2, -1, -1, -2, -2, -1])
T9_WORDS = ([2, -1, -2],
            [-3, 2, -1, -2, 3],
            [-3, -3, -2, -1, 2, 3, 3])
T9_ACTIONS = (((3, 4), (5, 6)),
              ((3, 8), (4, 7)),
              ((4, 6), (7, 8)))

def inv_word(word):
    return [-c for c in reversed(word)]

def restriction(g, outside):
    index = {a: i for i, a in enumerate(outside)}
    return tuple(index[g[a]] for a in outside)

def closure(generators, degree):
    ident = tuple(range(degree))
    seen, queue = {ident}, deque([ident])
    while queue:
        g = queue.popleft()
        for h in generators:
            new = core.compose(h, g)
            if new not in seen:
                seen.add(new)
                queue.append(new)
    return seen

def pure(word, r):
    order = list(range(r))
    for code in word:
        if abs(code) > 1:
            j = abs(code) - 1
            order[j], order[j + 1] = order[j + 1], order[j]
    return order == list(range(r))

def even(g):
    return sum(g[i] > g[j] for i in range(len(g))
               for j in range(i + 1, len(g))) % 2 == 0

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def strong_cable(t, word, action):
    """Replay all attachments and ordered splittings, holding the split fixed."""
    n, delta = len(t[0]), t[0]
    outside = [a for a in range(n) if delta[a] == a]
    count = 0
    for a in outside:
        for p, q in core.splittings(delta):
            old = (core.permutation(n + 1, [(n, a), p]),
                   core.permutation(n + 1, [(n, a), q]))
            old += tuple(g + (n,) for g in t[1:])
            new = old
            for code in core.cable(word):
                new = core.hurwitz(new, abs(code) - 1, 1 if code > 0 else -1)
            expected = (core.permutation(n + 1, [(n, action[a]), p]),
                        core.permutation(n + 1, [(n, action[a]), q]))
            expected += old[2:]
            assert new == expected
            assert core.product(new) == core.product(old)
            count += 1
    return count

def checked_word(t, word, expected=None, cables=False):
    new, action = core.replay_codes(t, word)
    assert new == t
    assert pure(word, len(t))
    assert even(action)
    assert core.conjugate(action, t[0]) == t[0]
    if expected is not None:
        assert action == expected
    # Preserve each transposition of every ordered splitting, not just delta.
    n = len(t[0])
    for p, q in core.splittings(t[0]):
        assert core.conjugate(action, core.permutation(n, [p])) == core.permutation(n, [p])
        assert core.conjugate(action, core.permutation(n, [q])) == core.permutation(n, [q])
    return action, strong_cable(t, word, action) if cables else 0

def find_generator_word(generators, target):
    n = len(target)
    ident = tuple(range(n))
    paths, queue = {ident: []}, deque([ident])
    while queue:
        g = queue.popleft()
        if g == target:
            return paths[g]
        for i, h in enumerate(generators):
            new = core.compose(h, g)
            if new not in paths:
                paths[new] = paths[g] + [i]
                queue.append(new)
    raise AssertionError("target is not generated")

def outside_set(t):
    return [a for a in range(len(t[0])) if t[0][a] == a]

def dt_controls(r):
    t, c = core.base_tuple(KINDS[0], r), 6
    n, outside = len(t[0]), outside_set(t)
    assert len(core.cycles(core.product(t))) == 1
    checks, cables = 0, 0
    for word, expected_cycles in zip(DT_WORDS, core.SEEDS[KINDS[0]]["actions"]):
        _, k = checked_word(t, word, core.permutation(n, expected_cycles), r <= 6)
        checks += 1
        cables += k
    # On the old outside set, these are (5,6) and (2,6).
    bc = list(DT_WORDS[0])
    ac = list(DT_WORDS[1]) + bc + inv_word(DT_WORDS[1])
    for index in range(3, r):
        ext = core.extension_word(index)
        _, k = checked_word(t, ext, t[index], r <= 6)
        checks += 1
        cables += k
        x = 7 + 2 * (index - 3)
        for old_word, fresh in ((ac, x), (bc, x + 1)):
            word = ext + old_word + inv_word(ext)
            action, k = checked_word(t, word, cables=r <= 6)
            assert restriction(action, outside) == restriction(
                core.permutation(n, [(fresh, c)]), outside)
            checks += 1
            cables += k
    return {"r": r, "degree": n, "outside_size": len(outside),
            "pure_literal_loops_checked": checks,
            "ordered_split_cables_checked": cables,
            "explicit_new_star_transpositions": 2 * (r - 3)}

def mixed_seed():
    t = core.base_tuple(KINDS[1], 4)
    actions, cables = [], 0
    for word, expected in zip(T9_WORDS, T9_ACTIONS):
        action, k = checked_word(t, word, core.permutation(9, expected), True)
        assert all(action[i] == i for i in range(3))
        actions.append(action)
        cables += k
    outside = outside_set(t)
    group = closure([restriction(g, outside) for g in actions], 6)
    assert len(group) == 360
    target_words, generator_words = {}, {}
    for label in (3, 4):
        target = core.permutation(9, [(label, 5, 6)])
        indices = find_generator_word(actions, target)
        word = [code for i in indices for code in T9_WORDS[i]]
        checked_word(t, word, target)
        target_words[label] = word
        generator_words[str(label)] = indices
    return target_words, {
        "degree": 9, "loops": T9_WORDS, "full_actions": T9_ACTIONS,
        "outside_group_order": len(group),
        "ordered_split_cables_checked": cables,
        "old_three_cycle_generator_words_zero_based": generator_words,
        "old_three_cycle_compressed_words": {str(k): v for k, v in target_words.items()},
        "mutation_drop_third_loop_order": len(closure(
            [restriction(g, outside) for g in actions[:2]], 6)),
    }, group

def mixed_controls(r, target_words):
    assert r >= 4
    t = core.base_tuple(KINDS[1], r)
    n, outside = len(t[0]), outside_set(t)
    assert len(core.cycles(core.product(t))) == 1
    checks, cables = 0, 0
    for word, expected in zip(T9_WORDS, T9_ACTIONS):
        _, k = checked_word(t, word, core.permutation(n, expected), r <= 6)
        checks += 1
        cables += k
    for index in range(4, r):
        ext = core.extension_word(index)
        _, k = checked_word(t, ext, t[index], r <= 6)
        checks += 1
        cables += k
        x = 7 + 2 * (index - 3)
        for anchor, fresh in ((3, x), (4, x + 1)):
            word = ext + target_words[anchor] + inv_word(ext)
            action, k = checked_word(t, word, cables=r <= 6)
            assert action == core.permutation(n, [(fresh, 5, 6)])
            checks += 1
            cables += k
    return {"r": r, "degree": n, "outside_size": len(outside),
            "pure_literal_loops_checked": checks,
            "ordered_split_cables_checked": cables,
            "explicit_new_three_cycles": 2 * (r - 4)}

def degree_seven_schreier():
    t0 = core.base_tuple(KINDS[1], 3)
    outside, ident = outside_set(t0), tuple(range(7))
    v4 = {tuple(range(4)),
          (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)}
    # T_x transports the marked permutation from the seed to vertex x.
    transports, queue = {t0: ident}, deque([t0])
    generators, edges = set(), 0
    while queue:
        t = queue.popleft()
        for code in (1, -1, 2, -2):
            u, h = core.code_move(t, code)
            candidate = core.compose(h, transports[t])
            if u not in transports:
                transports[u] = candidate
                queue.append(u)
            g = core.compose(core.inverse(transports[u]), candidate)
            assert core.conjugate(g, t0[0]) == t0[0]
            assert even(g)
            restricted = restriction(g, outside)
            assert restricted in v4
            generators.add(restricted)
            edges += 1
    assert len(transports) == 196 and edges == 784
    assert closure(list(generators), 4) == v4
    # Independent complete fixed-product enumeration, not the BFS generator rule.
    ts = [core.permutation(7, [e]) for e in itertools.combinations(range(7), 2)]
    dts = {core.compose(a, b) for a in ts for b in ts
           if set(i for i in range(7) if a[i] != i).isdisjoint(
               i for i in range(7) if b[i] != i)}
    threes = {core.permutation(7, [c]) for c in itertools.permutations(range(7), 3)}
    assert len(dts) == 105 and len(threes) == 70
    P, all_tuples = core.product(t0), set()
    for first in threes:
        for second in dts:
            last = core.compose(core.inverse(core.compose(first, second)), P)
            if last in dts:
                all_tuples.add((first, second, last))
    assert len(all_tuples) == 196 and set(transports) == all_tuples
    # Both V4 generators also have pure loops fixing the marked support pointwise.
    pointwise_words = [list(core.SEEDS[KINDS[1]]["words"][0]),
                       list(core.SEEDS[KINDS[1]]["words"][1]) * 3]
    acts, cables = [], 0
    for word in pointwise_words:
        action, k = checked_word(t0, word, cables=True)
        assert action[:3] == (0, 1, 2)
        acts.append(restriction(action, outside))
        cables += k
    assert closure(acts, 4) == v4
    state_text = json.dumps(sorted(transports), separators=(",", ":")).encode()
    return {"literal_states": len(transports), "signed_edges_checked": edges,
            "independent_fixed_product_tuples": len(all_tuples),
            "complete_state_set_sha256": hashlib.sha256(state_text).hexdigest(),
            "all_schreier_transports_in_V4": True,
            "outside_group_order": 4,
            "outside_generators": sorted(generators),
            "pointwise_support_pure_words": pointwise_words,
            "ordered_split_cables_checked": cables}, v4

def ordered_orbit_counts(group, degree):
    ans = []
    for k in range(1, degree + 1):
        remaining = set(itertools.permutations(range(degree), k))
        count = 0
        while remaining:
            point = next(iter(remaining))
            orbit = {tuple(g[a] for a in point) for g in group}
            assert orbit <= remaining
            remaining -= orbit
            count += 1
        ans.append(count)
    return ans

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=HERE / "full_coalesced_attachment_holonomy_verification.json")
    args = parser.parse_args()
    started = time.perf_counter()
    target_words, t9, a6 = mixed_seed()
    assert t9["mutation_drop_third_loop_order"] < 360
    special, v4 = degree_seven_schreier()
    controls = {
        KINDS[0]: [dt_controls(r) for r in (3, 4, 5, 6, 12, 25)],
        KINDS[1]: [mixed_controls(r, target_words) for r in (4, 5, 6, 12, 25)],
    }
    s3 = set(itertools.permutations(range(3)))
    ordered = {"S3": ordered_orbit_counts(s3, 3),
               "V4": ordered_orbit_counts(v4, 4),
               "A6": ordered_orbit_counts(a6, 6)}
    assert ordered == {"S3": [1, 1, 1], "V4": [1, 3, 6, 6],
                       "A6": [1, 1, 1, 1, 2, 2]}
    # Negative controls catch the false full-A4 base and unrestricted-split claim.
    assert len(v4) != 12
    t = core.base_tuple(KINDS[1], 3)
    _, rotated = core.replay_codes(t, core.SEEDS[KINDS[1]]["words"][1])
    assert rotated[:3] != (0, 1, 2)
    assert core.replay_codes(core.base_tuple(KINDS[0], 4), [])[1] != core.base_tuple(KINDS[0], 4)[3]
    sources = [Path(__file__).resolve(), HERE / "FULL_COALESCED_ATTACHMENT_HOLONOMY.md",
               HERE / "verify_uniform_hurwitz_attachment.py",
               HERE / "ALTERNATING_HURWITZ_CYCLE_FOLDING.md"]
    receipt = {
        "schema_version": 1, "result": "PASS",
        "kind": "full_coalesced_attachment_holonomy_exact_certificate",
        "python": platform.python_version(),
        "source_sha256": {p.name: sha(p) for p in sources},
        "degree_seven_mixed_complete_schreier": special,
        "degree_nine_mixed_seed": t9,
        "uniform_word_controls": controls,
        "ordered_distinct_choice_orbit_counts": ordered,
        "elapsed_seconds": time.perf_counter() - started,
        "scope": "Literal pure stabilizer transports of the explicitly chosen coalesced bases: S_(2r-3) for DT r>=3, V4 for mixed r=3, A_(2r-2) for mixed r>=4. Pointwise-splitting lower bounds, complete196-state V4 upper certificate, and bounded literal cabling controls. Uniform maximality uses the written induction; no arbitrary polynomial-orbit connectivity or leaf-pair absorption claim.",
        "formal_verification": False, "external_human_review": False,
        "independent_external_replay": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))

if __name__ == "__main__":
    main()
