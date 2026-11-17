#!/usr/bin/env python3
"""Exact conjugacy-restricted GL(3,2) quotient capacity.

No knot panel or historical search receipt is used.  Construct the group from
all invertible binary matrices, construct 336 explicit automorphisms, and
exhaust generating pairs and normalized involution triples.  The uniform
kernel-capacity implication is proved in the companion mathematical note.
Run from any directory; --output records a compact current-byte receipt.
"""

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import sys


def matrix_action(columns):
    """Action of a binary matrix on the seven nonzero column vectors."""
    return tuple(
        (columns[0] if vector & 1 else 0)
        ^ (columns[1] if vector & 2 else 0)
        ^ (columns[2] if vector & 4 else 0)
        for vector in range(1, 8)
    )


def group_data():
    actions = sorted(
        action
        for columns in product(range(1, 8), repeat=3)
        if len(set(action := matrix_action(columns))) == 7
    )
    assert len(actions) == len(set(actions)) == 168
    index = {action: i for i, action in enumerate(actions)}
    identity = index[tuple(range(1, 8))]
    multiplication = [
        [index[tuple(left[value - 1] for value in right)] for right in actions]
        for left in actions
    ]
    inverses = [
        next(j for j in range(168) if multiplication[i][j] == identity)
        for i in range(168)
    ]
    # Associativity also follows from literal function composition; check it
    # exhaustively so corrupted multiplication data cannot pass unnoticed.
    assert all(
        multiplication[multiplication[a][b]][c]
        == multiplication[a][multiplication[b][c]]
        for a in range(168) for b in range(168) for c in range(168)
    )
    assert all(multiplication[inverses[i]][i] == identity for i in range(168))
    assert [
        i for i in range(168)
        if all(multiplication[i][j] == multiplication[j][i] for j in range(168))
    ] == [identity]
    return actions, multiplication, inverses, identity


def automorphisms(actions, multiplication, inverses):
    index = {action: i for i, action in enumerate(actions)}
    theta = []
    for inverse in inverses:
        action = actions[inverse]
        columns = [action[(1 << j) - 1] for j in range(3)]
        transpose_columns = [
            sum(((columns[i] >> j) & 1) << i for i in range(3))
            for j in range(3)
        ]
        theta.append(index[matrix_action(transpose_columns)])
    assert sorted(theta) == list(range(168))
    assert all(theta[theta[i]] == i for i in range(168))
    assert all(
        theta[multiplication[i][j]] == multiplication[theta[i]][theta[j]]
        for i in range(168) for j in range(168)
    )
    inner = [
        tuple(multiplication[multiplication[g][h]][inverses[g]] for h in range(168))
        for g in range(168)
    ]
    explicit = inner + [tuple(a[theta[h]] for h in range(168)) for a in inner]
    assert len(set(explicit)) == 336
    # The inner maps and theta normalize one another, so these are a group.
    # Direct closure is unnecessary for the capacity upper bound: 336 distinct
    # freely acting automorphisms already suffice.
    return explicit


def generated_order(generators, multiplication, identity):
    seen, pending = {identity}, [identity]
    while pending:
        element = pending.pop()
        for generator in generators:
            value = multiplication[element][generator]
            if value not in seen:
                seen.add(value)
                pending.append(value)
    return len(seen)


def compute():
    actions, multiplication, inverses, identity = group_data()
    explicit = automorphisms(actions, multiplication, inverses)
    involutions = [
        i for i in range(168)
        if i != identity and multiplication[i][i] == identity
    ]
    assert len(involutions) == 21
    first = involutions[0]
    assert {a[first] for a in explicit} == set(involutions)
    stabilizer = [a for a in explicit if a[first] == first]
    assert len(stabilizer) == 16
    histogram = Counter(
        generated_order((first, b, c), multiplication, identity)
        for b, c in product(involutions, repeat=2)
    )
    expected = {2: 1, 4: 16, 6: 32, 8: 72, 24: 208, 168: 112}
    assert dict(histogram) == expected and sum(histogram.values()) == 21 ** 2
    marked, generating_orbits = set(), 0
    for b, c in product(involutions, repeat=2):
        if (b, c) in marked:
            continue
        orbit = {(a[b], a[c]) for a in stabilizer}
        assert not marked.intersection(orbit)
        marked.update(orbit)
        if generated_order((first, b, c), multiplication, identity) == 168:
            assert len(orbit) == 16
            generating_orbits += 1
    assert len(marked) == 21 ** 2 and generating_orbits == 7
    generating_pairs = sum(
        generated_order((a, b), multiplication, identity) == 168
        for a, b in product(range(168), repeat=2)
    )
    assert generating_pairs == 19152 == 57 * 336
    return {
        "status": "PASS",
        "group_order": 168,
        "nonidentity_involutions": 21,
        "explicit_distinct_automorphisms": 336,
        "normalized_triple_generated_order_histogram": {
            str(k): histogram[k] for k in sorted(histogram)
        },
        "generating_involution_triples": 21 * histogram[168],
        "generating_involution_triple_orbits_under_explicit_automorphisms": 7,
        "at_most_three_meridians_kernel_capacity_upper_bound": 7,
        "sharp_four_sufficient_distinct_kernel_count": 8,
        "arbitrary_generating_pairs": generating_pairs,
        "arbitrary_two_generator_kernel_capacity_upper_bound": 57,
        "claim_boundary": (
            "Exact finite group enumeration and elementary capacity gate; "
            "no knot panel tested, no novelty, formal verification or external review claimed."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = compute()
    source = Path(__file__).resolve()
    report["checker_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    report["python_version"] = sys.version.split()[0]
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
