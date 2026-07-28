#!/usr/bin/env python3
"""Screen length-two quadratic target shears for bordered flatness.

Let F=(P,B,C) be the normalized foundational cubic gauge and, for distinct
indices i,j with remaining index k, put

    A = F_i + Q(F_j,F_k),
    L = F_j + R(A,F_k),

where Q and R have positive total degree at most two.  There are six ordered
words (i,j).  This script exhaustively enumerates their ten coefficients over
a small prime field and rejects a tuple as soon as

    K(L) = -grad(L)^T adj(Hess(L)) grad(L)

is nonzero at one of the selected source points.

This is a bounded modular search, not a characteristic-zero proof.  Its role
is to identify possible strata for exact coefficient extraction.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


Vector = tuple[int, int, int]
Matrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]


@dataclass(frozen=True)
class Jet:
    value: int
    gradient: Vector
    hessian: Matrix


x, y, z = sp.symbols("x y z")
t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
FOUNDATIONAL_MAP = (
    sp.expand(t * q),
    sp.expand(y + 3 * x * q),
    sp.expand(x * (5 - 3 * t) - x**3 * z),
)
TARGET_NAMES = ("P", "B", "C")

# The first five points separate most tuples quickly.  The remaining points
# probe mixed signs and all three nonzero coordinates.
DEFAULT_POINTS = (
    (0, 0, 0),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 1),
    (-1, 1, 2),
    (2, -1, 1),
    (1, 2, -1),
    (-2, 1, 1),
    (1, -2, 2),
)


def add_vectors(left: Vector, right: Vector, prime: int) -> Vector:
    return tuple((left[index] + right[index]) % prime for index in range(3))  # type: ignore[return-value]


def add_matrices(left: Matrix, right: Matrix, prime: int) -> Matrix:
    return tuple(
        tuple((left[row][column] + right[row][column]) % prime for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def scale_vector(scalar: int, vector: Vector, prime: int) -> Vector:
    return tuple((scalar * entry) % prime for entry in vector)  # type: ignore[return-value]


def scale_matrix(scalar: int, matrix: Matrix, prime: int) -> Matrix:
    return tuple(
        tuple((scalar * entry) % prime for entry in row)
        for row in matrix
    )  # type: ignore[return-value]


def outer(left: Vector, right: Vector, prime: int) -> Matrix:
    return tuple(
        tuple((left[row] * right[column]) % prime for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def add_jets(left: Jet, right: Jet, prime: int) -> Jet:
    return Jet(
        (left.value + right.value) % prime,
        add_vectors(left.gradient, right.gradient, prime),
        add_matrices(left.hessian, right.hessian, prime),
    )


def quadratic_jet(
    coefficients: tuple[int, int, int, int, int],
    first: Jet,
    second: Jet,
    prime: int,
) -> Jet:
    """Compose a positive-degree quadratic polynomial with two scalar jets."""

    c10, c01, c20, c11, c02 = coefficients
    u = first.value
    v = second.value
    derivative_u = (c10 + 2 * c20 * u + c11 * v) % prime
    derivative_v = (c01 + c11 * u + 2 * c02 * v) % prime
    value = (
        c10 * u
        + c01 * v
        + c20 * u * u
        + c11 * u * v
        + c02 * v * v
    ) % prime
    gradient = add_vectors(
        scale_vector(derivative_u, first.gradient, prime),
        scale_vector(derivative_v, second.gradient, prime),
        prime,
    )
    hessian = add_matrices(
        scale_matrix(derivative_u, first.hessian, prime),
        scale_matrix(derivative_v, second.hessian, prime),
        prime,
    )
    hessian = add_matrices(
        hessian,
        scale_matrix(2 * c20, outer(first.gradient, first.gradient, prime), prime),
        prime,
    )
    hessian = add_matrices(
        hessian,
        scale_matrix(
            c11,
            add_matrices(
                outer(first.gradient, second.gradient, prime),
                outer(second.gradient, first.gradient, prime),
                prime,
            ),
            prime,
        ),
        prime,
    )
    hessian = add_matrices(
        hessian,
        scale_matrix(2 * c02, outer(second.gradient, second.gradient, prime), prime),
        prime,
    )
    return Jet(value, gradient, hessian)


def bordered_invariant(jet: Jet, prime: int) -> int:
    """Evaluate -g^T adj(H) g for a symmetric three-by-three Hessian."""

    gx, gy, gz = jet.gradient
    a, b, c = jet.hessian[0]
    _, d, e = jet.hessian[1]
    _, _, f = jet.hessian[2]
    value = (
        gx * gx * (d * f - e * e)
        + 2 * gx * gy * (c * e - b * f)
        + 2 * gx * gz * (b * e - c * d)
        + gy * gy * (a * f - c * c)
        + 2 * gy * gz * (b * c - a * e)
        + gz * gz * (a * d - b * b)
    )
    return (-value) % prime


def foundational_jets(point: tuple[int, int, int], prime: int) -> tuple[Jet, Jet, Jet]:
    substitutions = dict(zip((x, y, z), point, strict=True))
    result = []
    for polynomial in FOUNDATIONAL_MAP:
        gradient = tuple(
            int(sp.diff(polynomial, variable).subs(substitutions)) % prime
            for variable in (x, y, z)
        )
        hessian = tuple(
            tuple(
                int(sp.diff(polynomial, first, second).subs(substitutions)) % prime
                for second in (x, y, z)
            )
            for first in (x, y, z)
        )
        result.append(
            Jet(
                int(polynomial.subs(substitutions)) % prime,
                gradient,  # type: ignore[arg-type]
                hessian,  # type: ignore[arg-type]
            )
        )
    return tuple(result)  # type: ignore[return-value]


def word_value(
    coefficients: tuple[int, ...],
    indices: tuple[int, int, int],
    jets: tuple[Jet, Jet, Jet],
    prime: int,
) -> int:
    first_index, second_index, remaining_index = indices
    first_shear = add_jets(
        jets[first_index],
        quadratic_jet(
            coefficients[:5],  # type: ignore[arg-type]
            jets[second_index],
            jets[remaining_index],
            prime,
        ),
        prime,
    )
    retained = add_jets(
        jets[second_index],
        quadratic_jet(
            coefficients[5:],  # type: ignore[arg-type]
            first_shear,
            jets[remaining_index],
            prime,
        ),
        prime,
    )
    return bordered_invariant(retained, prime)


def all_words() -> Iterable[tuple[int, int, int]]:
    for first_index in range(3):
        for second_index in range(3):
            if first_index != second_index:
                remaining_index = 3 - first_index - second_index
                yield first_index, second_index, remaining_index


def word_name(indices: tuple[int, int, int]) -> str:
    first, second, remaining = (TARGET_NAMES[index] for index in indices)
    return f"A={first}+Q({second},{remaining}); L={second}+R(A,{remaining})"


def search(prime: int, points: tuple[tuple[int, int, int], ...]) -> dict[str, object]:
    point_jets = [foundational_jets(point, prime) for point in points]
    words: list[dict[str, object]] = []
    for indices in all_words():
        survivors = list(itertools.product(range(prime), repeat=10))
        counts = [len(survivors)]
        for jets in point_jets:
            survivors = [
                coefficients
                for coefficients in survivors
                if word_value(coefficients, indices, jets, prime) == 0
            ]
            counts.append(len(survivors))
            if not survivors:
                break
        words.append(
            {
                "indices": indices,
                "word": word_name(indices),
                "survivor_counts": counts,
                "survivors": survivors,
            }
        )
        print(f"{word_name(indices)}: {' -> '.join(map(str, counts))}")
    return {
        "status": "bounded modular experiment",
        "prime": prime,
        "coefficient_order": [
            "q10",
            "q01",
            "q20",
            "q11",
            "q02",
            "r10",
            "r01",
            "r20",
            "r11",
            "r02",
        ],
        "points": points,
        "words": words,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSON result path",
    )
    arguments = parser.parse_args()
    if arguments.prime < 2 or not sp.isprime(arguments.prime):
        parser.error("--prime must be prime")

    result = search(arguments.prime, DEFAULT_POINTS)
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
