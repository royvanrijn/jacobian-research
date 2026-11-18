#!/usr/bin/env python3
"""Exact word certificate for uniform coalesced-polynomial attachment holonomy.

The four degree-seven seed loops and the degree-independent leaf-extension
proof are in ALTERNATING_HURWITZ_CYCLE_FOLDING.md, Corollary6.
No polynomial-space connectivity or higher-degree exhaustive census is claimed.
Only the Python standard library is used.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform

def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))

def inverse(p):
    ans = [0] * len(p)
    for i, j in enumerate(p):
        ans[j] = i
    return tuple(ans)

def permutation(n, cycles):
    p = list(range(n))
    for cycle in cycles:
        for i, a in enumerate(cycle):
            p[a] = cycle[(i + 1) % len(cycle)]
    return tuple(p)

def conjugate(g, x):
    return compose(compose(g, x), inverse(g))

def hurwitz(t, j, sign):
    a, b = t[j:j + 2]
    pair = (b, conjugate(inverse(b), a)) if sign == 1 else (conjugate(a, b), a)
    return t[:j] + pair + t[j + 2:]

def code_move(t, code):
    """Code1 is Q1^2; codej for j>=2 is Qj; negative codes invert."""
    j, sign = abs(code) - 1, 1 if code > 0 else -1
    ident = tuple(range(len(t[0])))
    if j == 0:
        a, b = t[:2]
        transport = b if sign == 1 else conjugate(a, b)
        return hurwitz(hurwitz(t, 0, sign), 0, sign), transport
    return hurwitz(t, j, sign), ident

def replay_codes(t, word):
    transport = tuple(range(len(t[0])))
    for code in word:
        t, step = code_move(t, code)
        transport = compose(step, transport)
    return t, transport

def cable(word):
    ans = []
    for code in word:
        if code == 1:
            ans.extend((2, 1, 1, 2))
        elif code == -1:
            ans.extend((-2, -1, -1, -2))
        else:
            ans.append((abs(code) + 1) * (1 if code > 0 else -1))
    return ans

def product(t):
    ans = tuple(range(len(t[0])))
    for p in t:
        ans = compose(ans, p)
    return ans

def cycles(p):
    seen, ans = set(), []
    for i in range(len(p)):
        if i in seen:
            continue
        cycle, j = [], i
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = p[j]
        ans.append(tuple(cycle))
    return ans

def edges(p):
    cs = cycles(p)
    assert all(len(c) in (1, 2) for c in cs)
    return {tuple(sorted(c)) for c in cs if len(c) == 2}

def edge_perm(n, e):
    return permutation(n, [e])

def splittings(delta):
    n = len(delta)
    candidates = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ans = []
    for e in candidates:
        p = edge_perm(n, e)
        q = compose(p, delta)
        if sorted(map(len, cycles(q))) == [1] * (n - 2) + [2]:
            ans.append((e, next(iter(edges(q)))))
    return ans

def attachment_orbits(outside, generators):
    parent = {a: a for a in outside}
    def root(a):
        while parent[a] != a:
            a = parent[a]
        return a
    for g in generators:
        assert {g[a] for a in outside} == set(outside)
        for a in outside:
            parent[root(g[a])] = root(a)
    return len({root(a) for a in outside})

SEEDS = {
    "double_transposition": {
        "cycles": [[(0, 1), (3, 4)], [(1, 2), (4, 5)], [(2, 3), (5, 6)]],
        "anchors": (2, 5),
        "words": [[2, 1, 1, 1, -2], [1, 2, 2, -1, -1, -2, -2, -1]],
        "actions": [[(3, 4), (5, 6)], [(0, 1), (2, 5)]],
    },
    "three_cycle": {
        "cycles": [[(0, 1, 2)], [(2, 3), (4, 5)], [(3, 4), (5, 6)]],
        "anchors": (3, 4),
        "words": [[2, 1, -2], [1, 1, 2, -1, -2, 1, -2, 1, 2]],
        "actions": [[(3, 4), (5, 6)], [(0, 1, 2), (3, 6), (4, 5)]],
    },
}

def base_tuple(kind, r):
    assert r >= 3
    data = SEEDS[kind]
    n = 2 * r + 1
    t = tuple(permutation(n, cs) for cs in data["cycles"])
    a, b = data["anchors"]
    for new in range(7, n, 2):
        t += (permutation(n, [(a, new), (b, new + 1)]),)
    return t

def extension_word(index):
    # index is zero-based; the chosen new factor occupies position index+1.
    assert index >= 3
    bring = list(range(index, 1, -1))
    return bring + [1] + [-x for x in reversed(bring)]

def check_cable(t, word, action):
    n = len(t[0])
    delta = t[0]
    outside = [a for a in range(n) if delta[a] == a]
    splitting_list = splittings(delta)
    expected_splittings = 2 if all(len(c) in (1, 2) for c in cycles(delta)) else 3
    assert len(splitting_list) == expected_splittings
    for a in outside:
        for p, q in splitting_list:
            first = permutation(n + 1, [(n, a), p])
            second = permutation(n + 1, [(n, a), q])
            original = (first, second) + tuple(x + (n,) for x in t[1:])
            moved = original
            for code in cable(word):
                moved = hurwitz(moved, abs(code) - 1, 1 if code > 0 else -1)
            assert moved[2:] == original[2:]
            assert compose(moved[0], moved[1]) == compose(original[0], original[1])
            assert edges(moved[0]) & edges(moved[1]) == {tuple(sorted((n, action[a])))}
            assert product(moved) == product(original)
    return len(outside) * len(splitting_list)

def check_case(kind, r, check_all_cables):
    t = base_tuple(kind, r)
    n, data = 2 * r + 1, SEEDS[kind]
    assert len(t) == r
    assert len(cycles(product(t))) == 1
    outside = [a for a in range(n) if t[0][a] == a]
    actions, counts = [], []
    for word, expected_cycles in zip(data["words"], data["actions"]):
        moved, action = replay_codes(t, word)
        assert moved == t
        assert action == permutation(n, expected_cycles)
        actions.append(action)
        if check_all_cables:
            counts.append(check_cable(t, word, action))
    for index in range(3, r):
        word = extension_word(index)
        moved, action = replay_codes(t, word)
        assert moved == t and action == t[index]
        actions.append(action)
        if check_all_cables:
            counts.append(check_cable(t, word, action))
    assert attachment_orbits(outside, actions) == 1
    return {
        "polynomial_degree": n, "factor_count": r,
        "outside_attachment_count": len(outside),
        "seed_stabilizers": len(data["words"]),
        "commuting_leaf_stabilizers": r - 3,
        "attachment_orbits": 1,
        "cabled_expansions_checked": sum(counts),
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
        default=Path(__file__).with_name("uniform_hurwitz_attachment_verification.json"))
    args = parser.parse_args()
    controls = (3, 4, 5, 6, 12, 25)
    cases = {kind: [check_case(kind, r, r <= 6) for r in controls] for kind in SEEDS}
    # Mutation controls: a claimed identity action cannot replace the second
    # seed transport, and the extension braid cannot be replaced by emptyword.
    for kind in SEEDS:
        t = base_tuple(kind, 3)
        _, action = replay_codes(t, SEEDS[kind]["words"][1])
        assert action != tuple(range(7))
        assert attachment_orbits([i for i in range(7) if t[0][i] == i],
            [permutation(7, SEEDS[kind]["actions"][0])]) > 1
    t = base_tuple("three_cycle", 4)
    assert replay_codes(t, [])[1] != t[3]
    receipt = {
        "schema_version": 1, "kind": "uniform_hurwitz_attachment_word_certificate",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seed_certificates": SEEDS, "deterministic_controls": cases,
        "scope": "Four exact degree7 stabilizer loops and their literal doubled-edge cables; leaf-extension controls. Uniformity follows from the written commuting-DT/tree proof. No complete higher-degree census or marked polynomial connectivity claim.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
if __name__ == "__main__":
    main()

