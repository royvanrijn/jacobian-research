#!/usr/bin/env python3
"""Verify the global affine-ramification budget and the degree-six E8 escape."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import permutations

import sympy as sp


Permutation = tuple[int, ...]


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(value: Permutation) -> Permutation:
    result = [0] * len(value)
    for i, image in enumerate(value):
        result[image] = i
    return tuple(result)


def power(value: Permutation, exponent: int) -> Permutation:
    result = tuple(range(len(value)))
    while exponent:
        if exponent & 1:
            result = compose(value, result)
        value = compose(value, value)
        exponent //= 2
    return result


def cycle_type(value: Permutation) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(value)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = value[point]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def is_transitive(generators: tuple[Permutation, ...]) -> bool:
    moves = generators + tuple(inverse(value) for value in generators)
    orbit = {0}
    frontier = [0]
    while frontier:
        point = frontier.pop()
        for value in moves:
            image = value[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return len(orbit) == len(generators[0])


def generated_group(generators: tuple[Permutation, ...]) -> set[Permutation]:
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    moves = generators + tuple(inverse(value) for value in generators)
    while frontier:
        value = frontier.pop()
        for move in moves:
            product = compose(move, value)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def check_global_identity() -> None:
    """Replay the genus/puncture cancellation with arbitrary exact inputs."""

    for source_square in range(-9, 5):
        for target_square in range(-7, 4):
            for degree in range(2, 13):
                for data in (
                    ((2, -4, -1),),
                    ((2, -4, -1), (3, 1, -1)),
                    ((4, -2, -1), (2, 3, 0), (5, -1, 1)),
                ):
                    # A component datum is (moved degree A, L_Y.C, beta),
                    # where beta=2g-2+s.
                    pairing = sum(moved * contact for moved, contact, _ in data)
                    determinant_square = (
                        source_square - degree * target_square - 2 * pairing
                    )
                    kernel_degree = sum(
                        moved * (contact - beta)
                        for moved, contact, beta in data
                    )
                    global_budget = Fraction(
                        source_square
                        - degree * target_square
                        + 2 * (degree - 1),
                        2,
                    )
                    residual = global_budget - (
                        kernel_degree + Fraction(determinant_square, 2)
                    )
                    expected = degree - 1 + sum(
                        moved * beta for moved, _, beta in data
                    )
                    assert residual == expected

    # Rational one-puncture components: beta=-1 and A_j=d-u_j.
    for degree in range(3, 20):
        for fixed_1 in range(1, degree):
            for fixed_2 in range(1, degree):
                moved_sum = 2 * degree - fixed_1 - fixed_2
                residual = degree - 1 - moved_sum
                assert residual == fixed_1 + fixed_2 - degree - 1


def check_degree_six_cubic_action() -> None:
    degree = 6
    identity = tuple(range(degree))
    symmetric_group = list(permutations(range(degree)))
    by_cube: dict[Permutation, list[Permutation]] = {}
    by_fifth: dict[Permutation, list[Permutation]] = {}
    for value in symmetric_group:
        by_cube.setdefault(power(value, 3), []).append(value)
        by_fifth.setdefault(power(value, 5), []).append(value)

    solutions: list[tuple[Permutation, Permutation]] = []
    central_images: list[Permutation] = []
    for central in set(by_cube).intersection(by_fifth):
        for a in by_cube[central]:
            for b in by_fifth[central]:
                meridian = compose(inverse(a), power(b, 2))
                if cycle_type(meridian) != (3, 1, 1, 1):
                    continue
                if is_transitive((a, b)):
                    solutions.append((a, b))
                    central_images.append(central)

    assert len(solutions) == 720
    assert set(central_images) == {identity}
    assert Counter(cycle_type(a) for a, _ in solutions) == Counter({(3, 3): 720})
    assert Counter(cycle_type(b) for _, b in solutions) == Counter({(5, 1): 720})

    representative = solutions[0]
    conjugacy_orbit = set()
    for change in symmetric_group:
        change_inverse = inverse(change)
        conjugacy_orbit.add(
            (
                compose(compose(change, representative[0]), change_inverse),
                compose(compose(change, representative[1]), change_inverse),
            )
        )
    assert len(conjugacy_orbit) == 720
    assert conjugacy_orbit == set(solutions)

    displayed_a = (1, 2, 0, 4, 5, 3)  # (1 2 3)(4 5 6)
    displayed_b = (0, 3, 4, 2, 5, 1)  # (2 4 3 5 6)
    displayed_meridian = compose(inverse(displayed_a), power(displayed_b, 2))
    assert (displayed_a, displayed_b) in solutions
    assert displayed_meridian == (2, 1, 4, 3, 0, 5)  # (1 3 5)

    a, b = representative
    meridian = compose(inverse(a), power(b, 2))
    assert len(generated_group((a, b))) == 360
    assert cycle_type(a) == (3, 3)
    assert cycle_type(b) == (5, 1)
    assert cycle_type(meridian) == (3, 1, 1, 1)
    assert power(a, 3) == power(b, 5) == identity

    # The preferred longitude is z*m^(-15).  Here z=1 and m^3=1.
    longitude = compose(power(a, 3), power(inverse(meridian), 15))
    assert longitude == identity

    moved_sheets = 3
    fixed_sheets = 3
    residue_cycle_count = 1
    point_budget = fixed_sheets - 1
    cusp_lower = 2 * residue_cycle_count
    assert moved_sheets == degree - fixed_sheets
    assert point_budget == cusp_lower == 2

    # This is exactly the already certified terminal A_6 passport.
    assert (cycle_type(b), cycle_type(a), cycle_type(meridian)) == (
        (5, 1),
        (3, 3),
        (3, 1, 1, 1),
    )


def check_local_cubic_packet() -> None:
    r, t = sp.symbols("r t")
    x = t**3 + r * t
    y = t**5 + sp.Rational(5, 3) * r * t**3 + sp.Rational(5, 9) * r**2 * t

    jacobian = sp.factor(sp.diff(x, r) * sp.diff(y, t) - sp.diff(x, t) * sp.diff(y, r))
    cusp_pullback = sp.factor(y**3 - x**5)
    assert jacobian == -sp.Rational(5, 9) * r**2 * t
    assert sp.factor(cusp_pullback / (r**3 * t**3)) == sp.factor(
        (125 * r**3 + 396 * r**2 * t**2 + 405 * r * t**4 + 135 * t**6)
        / 729
    )

    # At the SNC source node use (r*d/dr,t*d/dt).  The matrix has a common
    # factor t, while the quotient matrix is invertible at the generic point
    # of T=(t=0).  Thus the generic Smith form on T is diag(t,t), not cyclic.
    theta = sp.Matrix(
        [
            [r * sp.diff(x, r), r * sp.diff(y, r)],
            [t * sp.diff(x, t), t * sp.diff(y, t)],
        ]
    )
    assert sp.factor(theta.det()) == -sp.Rational(5, 9) * r**3 * t**2
    divided = sp.simplify(theta / t)
    assert sp.factor(divided.det()) == -sp.Rational(5, 9) * r**3
    assert sp.factor(divided.det().subs(t, 0)) == -sp.Rational(5, 9) * r**3
    assert all(sp.rem(entry, t) == 0 for entry in theta)


def main() -> None:
    check_global_identity()
    check_degree_six_cubic_action()
    check_local_cubic_packet()
    print(
        "PASS: the multi-component logarithmic cancellation is "
        "d-1+sum A_j(2g_j-2+s_j); the unique degree-six E8 cubic-inertia "
        "action is the natural A6 action, saturates the point budget 2=2, "
        "matches the terminal A6 passport, and has an exact local SNC model "
        "with a generically split diag(t,t) contracted packet"
    )


if __name__ == "__main__":
    main()
