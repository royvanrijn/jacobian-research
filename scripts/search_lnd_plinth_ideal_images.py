#!/usr/bin/env python3
"""Exact bounded LNED search across a nontrivial plinth divisor.

Work over Q[x,y,z] with the linear locally nilpotent derivation

    D = x*d/dy + y*d/dz.

The local slice y/x has plinth element x, and the invariant
2*x*z-y^2 records the normalization denominator.  For each homogeneous
test ideal, the script constructs D(I_n) exactly in every required degree,
using rational linear algebra.  It then tests sparse homogeneous seeds
through six pure powers and four fixed mixed multipliers.

This is candidate generation only.  A bounded prefix is not a proof of
pure membership, and a bounded mixed tail is not a counterexample.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def derivation(f: sp.Expr) -> sp.Expr:
    return sp.expand(x * sp.diff(f, y) + y * sp.diff(f, z))


def monomials_of_degree(degree: int) -> tuple[sp.Expr, ...]:
    if degree < 0:
        return ()
    return tuple(
        x**i * y**j * z ** (degree - i - j)
        for i in range(degree + 1)
        for j in range(degree - i + 1)
    )


def homogeneous_degree(f: sp.Expr) -> int:
    degrees = {
        sum(exponents)
        for exponents, coefficient in sp.Poly(f, *variables).terms()
        if coefficient
    }
    assert len(degrees) == 1
    return degrees.pop()


def coefficient_vector(f: sp.Expr, basis: tuple[sp.Expr, ...]) -> sp.Matrix:
    poly = sp.Poly(sp.expand(f), *variables, domain=sp.QQ)
    return sp.Matrix(
        [
            poly.coeff_monomial(monomial)
            for monomial in basis
        ]
    )


@dataclass(frozen=True)
class HomogeneousIdeal:
    name: str
    generators: tuple[sp.Expr, ...]


invariant = 2 * x * z - y**2
CHARTS = (
    HomogeneousIdeal("monomial-control", (x**2, x * y, y**2)),
    HomogeneousIdeal("invariant-double", (x**2, invariant)),
    HomogeneousIdeal("tilted-conductor", (x**2, y**2 + x * z)),
    HomogeneousIdeal("jet-coupling", (x**2, x * y, y**2 + x * z)),
    HomogeneousIdeal(
        "plinth-thickening",
        (x**3, x * y, y**2 + x * z),
    ),
)


@cache
def image_constraints(
    chart: HomogeneousIdeal, degree: int
) -> tuple[sp.Matrix, ...]:
    basis = monomials_of_degree(degree)
    source_spanning: list[sp.Expr] = []
    for generator in chart.generators:
        generator_degree = homogeneous_degree(generator)
        for multiplier in monomials_of_degree(degree - generator_degree):
            source_spanning.append(sp.expand(generator * multiplier))

    image_vectors = [
        coefficient_vector(derivation(source), basis)
        for source in source_spanning
    ]
    if image_vectors:
        image_matrix = sp.Matrix.hstack(*image_vectors)
    else:
        image_matrix = sp.zeros(len(basis), 0)
    return tuple(image_matrix.T.nullspace())


def in_image(f: sp.Expr, chart: HomogeneousIdeal) -> bool:
    degree = homogeneous_degree(f)
    vector = coefficient_vector(f, monomials_of_degree(degree))
    return all((constraint.T * vector)[0] == 0 for constraint in image_constraints(chart, degree))


def seed_family() -> tuple[sp.Expr, ...]:
    seeds: list[sp.Expr] = []
    for degree in (1, 2):
        monomials = monomials_of_degree(degree)
        seeds.extend(monomials)
        for left, right in combinations(monomials, 2):
            seeds.append(left + right)
            seeds.append(left - right)
    return tuple(seeds)


def main() -> None:
    pure_bound = 6
    mixed_start = 4
    seeds = seed_family()
    multipliers = (sp.Integer(1), x, y, z)

    assert derivation(x) == 0
    assert derivation(y) == x
    assert derivation(z) == y
    assert derivation(invariant) == 0
    assert derivation(derivation(derivation(z))) == 0

    print(f"SEARCH: D=x*d_y+y*d_z; {len(seeds)} sparse homogeneous seeds")
    total_survivors = 0
    total_tail_obstructions = 0

    for chart in CHARTS:
        survivors: list[sp.Expr] = []
        tail_obstructions: list[tuple[sp.Expr, sp.Expr]] = []

        for f in seeds:
            if all(in_image(sp.expand(f**m), chart) for m in range(1, pure_bound + 1)):
                survivors.append(f)
                for g in multipliers:
                    failures = [
                        not in_image(sp.expand(g * f**m), chart)
                        for m in range(mixed_start, pure_bound + 1)
                    ]
                    if all(failures):
                        tail_obstructions.append((f, g))

        total_survivors += len(survivors)
        total_tail_obstructions += len(tail_obstructions)
        print(
            "CHART:",
            chart.name,
            f"pure-prefix survivors={len(survivors)}",
            f"bounded-tail obstructions={len(tail_obstructions)}",
        )
        if chart.name != "monomial-control":
            for f in survivors:
                print("  SURVIVOR:", f"f={f}")
        for f, g in tail_obstructions[:5]:
            print("  CANDIDATE-ONLY:", f"f={f}", f"g={g}")

    print(f"SUMMARY: pure-prefix survivors={total_survivors}")
    print(f"SUMMARY: bounded-tail obstructions={total_tail_obstructions}")
    print("NOTE: bounded survival or failure is not an all-order certificate")


if __name__ == "__main__":
    main()
