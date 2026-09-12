#!/usr/bin/env python3
"""Explore degree-six transposition actions on the exact F2 r=8 cusp pair.

This is a numerical braid-continuation experiment, not a proof.  Exact
algebra in ``verify_f2_geometric_degree_six_stein_reduction.py`` isolates the
quartic cusp locus.  Here one real embedding is used to discover the finite
permutation relations that a future certified braid computation should test.
"""

from __future__ import annotations

import argparse
import cmath
from itertools import product
import math

import sympy as sp

import verify_f2_affine_target_k1_implicit_conductor as k1


def polynomial_roots(coefficients: list[complex]) -> list[complex]:
    variable = sp.symbols("root_variable")
    expression = sum(
        sp.N(coefficient, 40) * variable ** (len(coefficients) - index - 1)
        for index, coefficient in enumerate(coefficients)
    )
    return [complex(root) for root in sp.nroots(expression, n=35, maxsteps=200)]


def continue_cubic(
    previous: list[complex],
    target: complex,
    linear: complex,
    constant: complex,
) -> list[complex]:
    roots = []
    for initial in previous:
        root = initial
        for _ in range(20):
            value = root**3 + linear * root + constant - target
            derivative = 3 * root**2 + linear
            correction = value / derivative
            root -= correction
            if abs(correction) < 1e-13:
                break
        assert abs(root**3 + linear * root + constant - target) < 1e-8
        roots.append(root)
    assert min(abs(roots[i] - roots[j]) for i in range(3) for j in range(i)) > 1e-8
    return roots


def implicit_quintic(
    aa: sp.Expr,
    bb: sp.Expr,
    cc: sp.Expr,
    dd: sp.Expr,
    pp: sp.Expr,
    qq: sp.Expr,
) -> sp.Expr:
    return sp.expand(
        k1.expected_implicit_quintic().subs(
            {
                k1.a: aa,
                k1.b: bb,
                k1.c: cc,
                k1.d: dd,
                k1.P: pp,
                k1.Q: qq,
            },
            simultaneous=True,
        )
    )


def exact_witness(real_root_index: int = 1) -> dict[str, complex]:
    a = sp.symbols("a")
    quartic = (
        196000000 * a**4
        + 260940000 * a**3
        + 82362825 * a**2
        - 2390688 * a
        + 20736
    )
    real_roots = sorted(
        float(sp.re(root))
        for root in sp.nroots(quartic, n=35, maxsteps=200)
        if abs(float(sp.im(root))) < 1e-20
    )
    aa = real_roots[real_root_index]
    cc = -(
        51940000 * aa**3 + 65374350 * aa**2 + 15840099 * aa - 898128
    ) / 4361202
    dd = -(
        53410000 * aa**3 + 72989725 * aa**2 + 21185232 * aa - 301824
    ) / 7268670
    shift = -(
        7840000 * aa**3 + 22668000 * aa**2 + 25088409 * aa + 4472496
    ) / 5814936
    p0 = (
        42140000 * aa**3 + 55884050 * aa**2 + 19349013 * aa + 295344
    ) / 2422890

    A = aa + shift
    E = 5 * shift / 3
    C = (4 * shift + 3 * cc + 5 * p0) / 3
    D = (5 * A**2 - 15 * A * aa + 10 * aa**2 + 12 * p0 + 9 * dd) / 9
    q0 = (
        2 * A**2
        - 8 * A * aa
        + 6 * A * cc
        + 10 * A * p0
        + 6 * aa**2
        - 6 * aa * cc
        - 15 * aa * p0
    ) / 9

    P, Q = k1.P, k1.Q
    first = k1.expected_implicit_quintic()
    second = implicit_quintic(
        sp.Float(A, 35),
        sp.Integer(1),
        sp.Float(C, 35),
        sp.Float(D - E * A, 35),
        P - sp.Float(p0, 35),
        Q - sp.Float(q0, 35) - sp.Float(E, 35) * (P - sp.Float(p0, 35)),
    )
    difference = sp.Poly(sp.expand(second - first.subs({
        k1.a: sp.Float(aa, 35),
        k1.b: 1,
        k1.c: sp.Float(cc, 35),
        k1.d: sp.Float(dd, 35),
    })), P, Q)
    h10 = complex(difference.coeff_monomial(P))
    h00 = complex(difference.coeff_monomial(1))
    assert abs(h10) > 1e-8

    return {
        "a": aa,
        "c": cc,
        "d": dd,
        "A": A,
        "E": E,
        "C": C,
        "D": D,
        "p0": p0,
        "q0": q0,
        "dbar": D - E * A,
        "mutual_p": -h00 / h10,
    }


def discriminant_values(row: dict[str, complex]) -> list[tuple[complex, tuple[str, ...]]]:
    a, b, c, d = row["a"], row.get("b", 1), row["c"], row["d"]
    A, C, dbar, p0 = row["A"], row["C"], row["dbar"], row["p0"]

    p1 = lambda parameter: parameter**3 + a * parameter
    p2 = lambda parameter: parameter**3 + A * parameter + p0
    raw: list[tuple[str, complex]] = []
    for parameter in polynomial_roots([3, 0, a]):
        raw.append(("vertical-1", p1(parameter)))
    collision1 = [1, b, a, 2 * a * b - c, -(a**2 + d)]
    for parameter in polynomial_roots(collision1):
        raw.append(("collision-1", -parameter * (parameter**2 + a)))
    for parameter in polynomial_roots([3, 0, A]):
        raw.append(("vertical-2", p2(parameter)))
    collision2 = [1, b, A, 2 * A * b - C, -(A**2 + dbar)]
    for parameter in polynomial_roots(collision2):
        raw.append(("collision-2", p0 - parameter * (parameter**2 + A)))
    raw.append(("mutual", row["mutual_p"]))

    merged: list[list[object]] = []
    for label, value in raw:
        for item in merged:
            if abs(item[0] - value) < 1e-7:
                item[1].append(label)
                break
        else:
            merged.append([value, [label]])
    merged.sort(key=lambda item: (item[0].real, item[0].imag))
    return [(item[0], tuple(item[1])) for item in merged]


def segment(start: complex, end: complex, count: int) -> list[complex]:
    return [start + (end - start) * index / count for index in range(1, count + 1)]


def star_loops(
    critical: list[tuple[complex, tuple[str, ...]]],
    tail_steps: int = 180,
    circle_steps: int = 180,
) -> tuple[complex, list[tuple[list[complex], int, int]]]:
    radii = []
    for index, (point, _) in enumerate(critical):
        nearest = min(
            abs(point - other)
            for other_index, (other, _) in enumerate(critical)
            if other_index != index
        )
        radii.append(0.15 * nearest)

    candidates = (2 + 3j, 3 + 4j, -2 + 3j, 4 + 1.5j, -3 + 4j)
    for base in candidates:
        valid = True
        for index, (point, _) in enumerate(critical):
            edge = point - base
            for other_index, (other, _) in enumerate(critical):
                if other_index == index:
                    continue
                projection = (
                    (other - base).real * edge.real
                    + (other - base).imag * edge.imag
                ) / abs(edge) ** 2
                projection = max(0.0, min(1.0, projection))
                distance = abs(other - (base + projection * edge))
                if distance < 1.5 * radii[other_index]:
                    valid = False
        if valid:
            break
    else:
        raise AssertionError("no disjoint star system was found")

    loops = []
    for index, (center, _) in enumerate(critical):
        direction = (base - center) / abs(base - center)
        approach = center + radii[index] * direction
        angle = cmath.phase(direction)
        path = [base] + segment(base, approach, tail_steps)
        path.extend(
            center
            + radii[index] * cmath.exp(1j * (angle + 2 * math.pi * step / circle_steps))
            for step in range(1, circle_steps + 1)
        )
        path.extend(segment(approach, base, tail_steps))
        loops.append((path, tail_steps, circle_steps))
    return base, loops


def strand_trajectories(
    path: list[complex],
    row: dict[str, complex],
    scale: complex = 1e8 * cmath.exp(1j * math.pi / 4),
) -> list[list[complex]]:
    a, b, c, d = row["a"], row.get("b", 1), row["c"], row["d"]
    A, E, C, D, p0, q0 = (
        row["A"],
        row["E"],
        row["C"],
        row["D"],
        row["p0"],
        row["q0"],
    )
    first_roots = polynomial_roots([1, 0, a, -path[0]])
    second_roots = polynomial_roots([1, 0, A, p0 - path[0]])
    trajectories: list[list[complex]] = [[] for _ in range(6)]

    def q1(parameter: complex) -> complex:
        return scale * (
            parameter**5 + b * parameter**4 + c * parameter**2 + d * parameter
        )

    def q2(parameter: complex) -> complex:
        return scale * (
            parameter**5
            + b * parameter**4
            + E * parameter**3
            + C * parameter**2
            + D * parameter
            + q0
        )

    for index, parameter in enumerate(first_roots + second_roots):
        trajectories[index].append(q1(parameter) if index < 3 else q2(parameter))
    for target in path[1:]:
        first_roots = continue_cubic(first_roots, target, a, 0)
        second_roots = continue_cubic(second_roots, target, A, p0)
        for index, parameter in enumerate(first_roots + second_roots):
            trajectories[index].append(q1(parameter) if index < 3 else q2(parameter))
    return trajectories


def braid_from_piecewise(trajectories: list[list[complex]]) -> tuple[int, ...]:
    """Port the crossing-order core of Sage's braid_from_piecewise."""

    word: list[int] = []
    strand_count = len(trajectories)

    def sign(left: float, right: float) -> int:
        return 1 if left < right else -1 if left > right else 0

    for step in range(len(trajectories[0]) - 1):
        pairs = sorted(
            [
                (
                    [trajectory[step].real, trajectory[step].imag],
                    [trajectory[step + 1].real, trajectory[step + 1].imag],
                )
                for trajectory in trajectories
            ]
        )
        left = [pair[0] for pair in pairs]
        right = [pair[1] for pair in pairs]
        crossings = []
        for later, right_later in enumerate(right):
            for earlier in range(later):
                if right_later < right[earlier]:
                    crossing_time = (left[later][0] - left[earlier][0]) / (
                        (right[earlier][0] - right_later[0])
                        + (left[later][0] - left[earlier][0])
                    )
                    imag_earlier = (
                        left[earlier][1] * (1 - crossing_time)
                        + crossing_time * right[earlier][1]
                    )
                    imag_later = (
                        left[later][1] * (1 - crossing_time)
                        + crossing_time * right_later[1]
                    )
                    crossings.append(
                        [
                            crossing_time,
                            earlier,
                            later,
                            sign(imag_earlier, imag_later),
                        ]
                    )
        crossings.sort()
        permutation = list(range(strand_count))
        while crossings:
            same_time = [
                crossing for crossing in crossings if crossing[0] == crossings[0][0]
            ]
            crossings = crossings[len(same_time):]
            pending = [
                (permutation[crossing[2]] - permutation[crossing[1]], crossing)
                for crossing in same_time
            ]
            while pending:
                pending.sort(key=lambda item: item[0])
                _, crossing = pending.pop(0)
                _, earlier, later, crossing_sign = crossing
                word.append(
                    crossing_sign
                    * (min(permutation[earlier], permutation[later]) + 1)
                )
                permutation[earlier], permutation[later] = (
                    permutation[later],
                    permutation[earlier],
                )
                pending = [
                    (
                        permutation[item[1][2]] - permutation[item[1][1]],
                        item[1],
                    )
                    for item in pending
                ]
    return tuple(word)


def transposition_actions(
    words: list[tuple[int, ...]],
    cusp_tail: tuple[int, ...],
    cusp_local: tuple[int, ...],
) -> tuple[int, int, int]:
    identity = tuple(range(6))
    transpositions = []
    for left in range(6):
        for right in range(left + 1, 6):
            permutation = list(identity)
            permutation[left], permutation[right] = right, left
            transpositions.append(tuple(permutation))

    def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(left[right[index]] for index in range(6))

    conjugate = [[0] * 15 for _ in range(15)]
    lookup = {permutation: index for index, permutation in enumerate(transpositions)}
    for outer_index, outer in enumerate(transpositions):
        for inner_index, inner in enumerate(transpositions):
            conjugate[outer_index][inner_index] = lookup[
                compose(outer, compose(inner, outer))
            ]

    def act(values: tuple[int, ...], word: tuple[int, ...]) -> tuple[int, ...]:
        result = list(values)
        for letter in word:
            position = abs(letter) - 1
            left, right = result[position], result[position + 1]
            if letter > 0:
                result[position] = right
                result[position + 1] = conjugate[right][left]
            else:
                result[position] = conjugate[left][right]
                result[position + 1] = left
        return tuple(result)

    def transitive(values: tuple[int, ...]) -> bool:
        graph = [set() for _ in range(6)]
        for index in values:
            moved = [sheet for sheet in range(6) if transpositions[index][sheet] != sheet]
            graph[moved[0]].add(moved[1])
            graph[moved[1]].add(moved[0])
        reached, frontier = {0}, [0]
        while frontier:
            sheet = frontier.pop()
            for neighbor in graph[sheet]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    frontier.append(neighbor)
        return len(reached) == 6

    survivors = transitive_survivors = cusp_s3_survivors = 0
    examples: list[tuple[tuple[int, int], ...]] = []
    # Simultaneous conjugacy lets the first transposition be fixed.
    for tail in product(range(15), repeat=5):
        values = (0,) + tail
        if any(act(values, word) != values for word in words):
            continue
        survivors += 1
        examples.append(
            tuple(
                tuple(index for index in range(6) if transpositions[value][index] != index)
                for value in values
            )
        )
        if transitive(values):
            transitive_survivors += 1
        local_values = act(values, cusp_tail)
        assert act(local_values, cusp_local) == local_values
        generator = abs(cusp_local[0]) - 1
        left = transpositions[local_values[generator]]
        right = transpositions[local_values[generator + 1]]
        moved_left = {index for index in range(6) if left[index] != index}
        moved_right = {index for index in range(6) if right[index] != index}
        if left != right and len(moved_left & moved_right) == 1:
            cusp_s3_survivors += 1
    print("surviving transposition-edge tuples:")
    for example in examples:
        print(" ", example)
    return survivors, transitive_survivors, cusp_s3_survivors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--real-root-index",
        type=int,
        choices=(0, 1),
        default=1,
        help="choose one of the two real embeddings of the cusp quartic",
    )
    arguments = parser.parse_args()
    row = exact_witness(arguments.real_root_index)
    critical = discriminant_values(row)
    base, loops = star_loops(critical)
    print(f"witness a={row['a']:.15g}; base={base}; critical values={len(critical)}")
    words = []
    cusp_tail: tuple[int, ...] | None = None
    cusp_local: tuple[int, ...] | None = None
    for index, ((path, tail_steps, circle_steps), (_, labels)) in enumerate(
        zip(loops, critical, strict=True)
    ):
        trajectories = strand_trajectories(path, row)
        if index == 0:
            component_order = tuple(
                component
                for _, component in sorted(
                    [
                        (trajectory[0], 1 if position < 3 else 2)
                        for position, trajectory in enumerate(trajectories)
                    ],
                    key=lambda item: (item[0].real, item[0].imag),
                )
            )
            print("base-fiber component order:", component_order)
        word = braid_from_piecewise(trajectories)
        tail_word = braid_from_piecewise(
            [trajectory[: tail_steps + 1] for trajectory in trajectories]
        )
        local_word = braid_from_piecewise(
            [
                trajectory[tail_steps: tail_steps + circle_steps + 1]
                for trajectory in trajectories
            ]
        )
        words.append(word)
        print(index, labels, "word", word, "local", local_word)
        if labels == ("vertical-1", "collision-1"):
            cusp_tail, cusp_local = tail_word, local_word

    assert cusp_tail is not None and cusp_local is not None
    assert len(cusp_local) == 3
    assert len({abs(letter) for letter in cusp_local}) == 1
    counts = transposition_actions(words, cusp_tail, cusp_local)
    print(
        "labeled survivors with first transposition fixed: total=%d, "
        "transitive=%d, cusp-S3=%d" % counts
    )


if __name__ == "__main__":
    main()
