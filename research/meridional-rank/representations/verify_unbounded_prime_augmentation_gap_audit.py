#!/usr/bin/env python3
"""Exact controls for the unbounded prime augmentation-gap audit.

Complete finite domain: all five-transposition factorizations of a fixed
six-cycle. Universal controls: free-group identities for protected tails.
Bounded controls: literal alternating tuples for 5 <= b <= 32.
The companion note supplies all quantified proofs and the primeness limits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
INPUT = HERE / "alternating_target_prime_families_inputs.json"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def cycle(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for a, b in zip(c, c[1:] + c[:1]):
            p[a] = b
    return tuple(p)


def perm_hurwitz(values, word):
    values = list(values)
    for letter in word:
        j = abs(letter) - 1
        a, b = values[j:j+2]
        if letter > 0:
            values[j:j+2] = [b, compose(compose(inverse(b), a), b)]
        else:
            values[j:j+2] = [compose(compose(a, b), inverse(a)), a]
    return values


def connected(edges, n):
    reached = {0}
    while True:
        new = reached | {b for a, b in edges if a in reached} | {a for a, b in edges if b in reached}
        if new == reached:
            return len(reached) == n
        reached = new


def reduce_word(values):
    out = []
    for a in values:
        if out and out[-1] == -a:
            out.pop()
        else:
            out.append(a)
    return tuple(out)


def word_mul(*values):
    return reduce_word(a for value in values for a in value)


def word_inv(value):
    return tuple(-a for a in value[::-1])


def free_hurwitz(values, word):
    values = list(values)
    for letter in word:
        j = abs(letter) - 1
        a, b = values[j:j+2]
        values[j:j+2] = (
            [b, word_mul(word_inv(b), a, b)] if letter > 0
            else [word_mul(a, b, word_inv(a)), a]
        )
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    data = json.loads(INPUT.read_text())
    beta = data["beta"]
    require(beta == [1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2],
            "literal braid binding")
    edges = list(combinations(range(6), 2))
    edge_index = {e: j for j, e in enumerate(edges)}
    transpositions = [cycle(6, e) for e in edges]
    transposition_index = {p: j for j, p in enumerate(transpositions)}
    identity = tuple(range(6))
    target = cycle(6, tuple(range(6)))
    all_tuples = set()
    # Every tuple is uniquely its first four factors plus the forced fifth.
    for first in product(range(15), repeat=4):
        current = identity
        for j in first:
            current = compose(current, transpositions[j])
        last = compose(inverse(current), target)
        if last in transposition_index:
            all_tuples.add(first + (transposition_index[last],))
    require(len(all_tuples) == 1296, "complete six-cycle factorization count")
    require(all(connected([edges[j] for j in t], 6) for t in all_tuples),
            "all enumerated tuples generate S6")

    def conjugated(edge, by):
        a, b = edge
        c, d = by
        def move(x):
            return d if x == c else c if x == d else x
        return tuple(sorted((move(a), move(b))))

    def act(value):
        value = [edges[j] for j in value]
        for letter in beta + [4,4,4]:
            j = abs(letter)-1
            a, b = value[j:j+2]
            value[j:j+2] = (
                [b, conjugated(a, b)] if letter > 0
                else [conjugated(b, a), a]
            )
        return tuple(edge_index[e] for e in value)

    seen = set()
    census = Counter()
    representatives = {}
    for initial in sorted(all_tuples):
        if initial in seen:
            continue
        current = initial
        orbit = []
        while current not in seen:
            seen.add(current)
            orbit.append(current)
            current = act(current)
            require(current in all_tuples, "Hurwitz preserves finite domain")
        require(current == initial, "complete permutation cycle")
        census[len(orbit)] += 1
        representatives.setdefault(len(orbit), initial)
    expected = {2:15, 6:2, 10:3, 24:1, 66:3, 72:1, 147:2, 180:1, 228:2}
    require(dict(census) == expected, "exact literal period census")
    require(len(seen) == 1296, "exhaustive cycle partition")
    require(sum(length * count for length, count in census.items()) == 1296,
            "cycle cardinalities")
    require(sum(d*c for d,c in census.items() if 2 % d == 0) == 30,
            "five known returning conjugacy classes times centralizer six")
    require(sum(d*c for d,c in census.items() if 3 % d == 0) == 0,
            "first viable exponent has no S6 transposition quotient")
    witness_edges = data["finite_symmetric_control"]["witness_transpositions"]
    witness = [cycle(6, tuple(e)) for e in witness_edges]
    require(connected(witness_edges, 6), "existing full S6 witness")
    require(perm_hurwitz(witness, beta+[4,4,4]) != witness, "witness not period one")
    require(perm_hurwitz(witness, (beta+[4,4,4])*2) == witness, "witness period two")
    require(not connected([[0,1],[0,2],[0,3],[2,4],[0,1]],6),
            "disconnected false full-rank witness rejected")

    # Universal group-word identity for a35 = sigma3 sigma4 sigma3^-1.
    free_seed = [(j,) for j in range(1,6)]
    outer_c = free_seed[3]
    framed_new = word_mul(outer_c, free_seed[4], word_inv(outer_c))
    result = free_hurwitz(free_seed, [3,4,-3])
    require(result[0] == free_seed[0] and result[3] == outer_c,
            "protected outer meridians fixed as free words")
    require(result[2] == framed_new, "first framed Hurwitz entry")
    require(word_mul(outer_c,result[4],word_inv(outer_c)) ==
            word_mul(word_inv(framed_new),free_seed[2],framed_new),
            "second framed Hurwitz entry")
    require(word_mul(free_seed[1],free_seed[2],framed_new) ==
            word_mul(result[1],result[2],outer_c,result[4],word_inv(outer_c)),
            "protected framed ordered product")
    bad = free_hurwitz(free_seed, [-3,4,3])
    require(bad[3] != outer_c, "wrong band conjugation rejected")

    # Every later ordinary Hurwitz generator commutes with framing by C.
    six = [(j,) for j in range(1,7)]
    result = free_hurwitz(six,[5])
    framed_before = [word_mul(six[3],v,word_inv(six[3])) for v in six[4:]]
    framed_after = [word_mul(result[3],v,word_inv(result[3])) for v in result[4:]]
    require(framed_after == [framed_before[1],
            word_mul(word_inv(framed_before[1]),framed_before[0],framed_before[1])],
            "later framed generators are ordinary Hurwitz moves")

    checked = []
    for b in range(5,33):
        n = 2*b
        seed = [cycle(n,(0,1),(2,3)),cycle(n,(1,4),(2,3)),
                cycle(n,(0,5),(1,6)),cycle(n,(1,2),(6,7))]
        seed += [cycle(n,(0,2*j),(2,2*j+1)) for j in range(4,b)]
        band3 = [3,4,4,4,-3]
        braid = beta + band3 + [j for j in range(5,b) for _ in range(3)]
        require(perm_hurwitz(seed,braid) == seed, "protected alternating tuple return")
        strand = list(range(b))
        for letter in braid:
            j = abs(letter)-1
            strand[j],strand[j+1] = strand[j+1],strand[j]
        current = 0
        visited = set()
        while current not in visited:
            visited.add(current)
            current = strand[current]
        require(len(visited) == b, "protected elementary chain is a knot")
        require(len(braid) == 17+5+3*(b-5), "elementary protected braid length")
        checked.append(b)

    report = {
        "schema_version": 1,
        "kind": "complete_finite_period_and_universal_free_group_certificate",
        "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "complete_domain": {
            "degree": 6, "factors": 5, "transpositions": 15,
            "fixed_product": list(target), "factorization_count": 1296,
            "literal_period_census": {str(k):v for k,v in sorted(census.items())},
            "period_representatives": {
                str(k):[list(edges[j])for j in value]
                for k,value in sorted(representatives.items())
            },
        },
        "conclusion_with_written_proof": {
            "full_S6_transposition_quotient_iff": "2 divides n or 147 divides n",
            "ar_C_equals_5_on_these_admissible_prime_powers": True,
            "first_unexcluded_power": {"b":5,"n":3,"crossings":60,"genus":28},
            "protected_tail_paired_span_bound": "min(r,q)<=b-1 and r+q<=2b-1",
            "elementary_protected_chain_is_connected_sum": "K_0 # trefoil^(b-4)",
        },
        "controls": {
            "protected_action": "exact reduced words in a free group",
            "literal_alternating_checks_b": checked,
            "known_S6_n2_classes_crosscheck": 5,
            "invalid_controls": 2,
        },
        "scope": (
            "The period classification exhausts every minimal S6 transposition tuple. "
            "It does not exclude other rank-five KCH representations at n=3. "
            "The protected-tail all-b rank theorem and composite identification "
            "are written proofs; bounded color controls do not prove primeness."
        ),
        "independent_replay": False, "formal_verification": False,
        "external_review": False, "result": "PASS",
    }
    rendered = json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(rendered)
    print(rendered,end="")


if __name__ == "__main__":
    main()
