#!/usr/bin/env python3
"""Exact controls for the uniform unmarked-to-marked Hurwitz reduction."""
from __future__ import annotations
import argparse
import hashlib
from itertools import permutations
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

def power(p, k):
    if k < 0:
        return power(inverse(p), -k)
    ans = tuple(range(len(p)))
    while k:
        if k & 1:
            ans = compose(ans, p)
        p = compose(p, p)
        k >>= 1
    return ans

def perm(n, cycles):
    p = list(range(n))
    for cycle in cycles:
        for i, a in enumerate(cycle):
            p[a] = cycle[(i + 1) % len(cycle)]
    return tuple(p)

def conjugate(a, x):
    return compose(compose(a, x), inverse(a))

def cycle_lengths(p):
    done, lengths = set(), []
    for a in range(len(p)):
        if a in done:
            continue
        j, length = a, 0
        while j not in done:
            done.add(j)
            length += 1
            j = p[j]
        lengths.append(length)
    return sorted(lengths)

def product(t):
    ans = tuple(range(len(t[0])))
    for p in t:
        ans = compose(ans, p)
    return ans

def hurwitz(t, j, sign=1):
    a, b = t[j:j + 2]
    pair = (b, conjugate(inverse(b), a)) if sign == 1 else (conjugate(a, b), a)
    return t[:j] + pair + t[j + 2:]

def braid(t, word):
    for move in word:
        t = hurwitz(t, abs(move) - 1, 1 if move > 0 else -1)
    return t

def labels(r, word):
    lab = list(range(r))
    for move in word:
        j = abs(move) - 1
        lab[j], lab[j + 1] = lab[j + 1], lab[j]
    return tuple(lab)

def tuple_case(r):
    n = 2 * r + 1
    P = perm(n, [tuple(range(n))])
    c = power(P, -2)
    d = perm(n, [(0, 2), (3, 4)])
    t = tuple(conjugate(power(c, -i), d) for i in range(r))
    assert all(cycle_lengths(g) == [1] * (n - 4) + [2, 2] for g in t)
    h = compose(d, power(P, 2))
    assert cycle_lengths(h) == [1, r, r]
    assert power(c, r) == P
    assert product(t) == P
    cycling = list(range(1, r))
    cycled = braid(t, cycling)
    assert cycled == tuple(conjugate(power(P, 2), g) for g in t)
    assert cycled == t[1:] + (conjugate(inverse(P), t[0]),)
    full_twist = braid(t, cycling * r)
    assert full_twist == tuple(conjugate(inverse(P), g) for g in t)
    literal_word = cycling * n
    assert braid(t, literal_word) == t
    assert cycle_lengths(labels(r, literal_word)) == [r]
    assert labels(r, literal_word) == labels(r, cycling)
    three_cycle = None
    if r >= 4:
        e = compose(t[0], power(compose(t[0], t[1]), 2))
        assert e == perm(n, [(0, 4), (2, 3)])
        three_cycle = power(compose(e, t[2]), 2)
        assert three_cycle == perm(n, [(0, 6, 4)])
        assert cycled != t
    # Each meridian word in Lemma1 has identity strand permutation.
    for j in range(2, r + 1):
        bring = list(range(j - 1, 1, -1))
        word = bring + [1, 1] + [-x for x in reversed(bring)]
        assert labels(r, word) == tuple(range(r))
    corrupted = t[:-1] + (tuple(range(n)),)
    assert product(corrupted) != P
    assert cycle_lengths(labels(r, cycling * (n * r))) == [1] * r
    return {
        "r": r, "polynomial_degree": n,
        "tuple": [list(g) for g in t] if r <= 6 else None,
        "product_cycle_lengths": cycle_lengths(product(t)),
        "dP2_cycle_lengths": cycle_lengths(h),
        "literal_stabilizer_word_length": len(literal_word),
        "branch_label_orbits": 1,
        "three_cycle": list(three_cycle) if three_cycle is not None else None,
        "full_group_scope": "A_N by written primitivity/Jordan proof" if r >= 4 else "No full-A_N assertion",
    }

def coset_controls(r):
    # Abstract finite quotient control of the adjustment in Lemma3.
    cycle = perm(r, [tuple(range(r))])
    K = [power(cycle, j) for j in range(r)]
    count = 0
    for beta in permutations(range(r)):
        wanted = inverse(beta)[0]
        k = next(x for x in K if x[0] == wanted)
        assert compose(beta, k)[0] == 0
        count += 1
    return {"r": r, "all_Sr_elements_checked": count,
            "cyclic_stabilizer_transitive": True}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path,
        default=Path(__file__).with_name("unmarked_polynomial_hurwitz_reduction_verification.json"))
    args = ap.parse_args()
    checks = [tuple_case(r) for r in (2, 3, 4, 5, 6, 7, 12, 25, 50)]
    cosets = [coset_controls(r) for r in range(2, 8)]
    here = Path(__file__).resolve().parent
    inputs = ["UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md",
              "ALTERNATING_HURWITZ_CYCLE_FOLDING.md",
              "alternating_hurwitz_source_ledger.json"]
    report = {
        "schema_version": 1, "kind": "uniform_polynomial_marking_reduction",
        "result": "PASS", "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {x: hashlib.sha256((here / x).read_bytes()).hexdigest()
                         for x in inputs},
        "literal_tuple_and_word_controls": checks,
        "finite_coset_controls": cosets,
        "scope": "Uniform explicit DT tuple and literal cycling stabilizer; subgroup equality and all-r marking reduction use the written proofs. Finite controls do not prove either unmarked polynomial connectivity statement.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()

