#!/usr/bin/env python3
"""Direct searches for the all-k square/cube transfer basis.

This script never uses the refuted conductor-ribbon norm.  It works from the
exact division identity A^2 = S*Q + R and the coefficient equations of
9*Q^2 + 8*S*A*Q - 48*S^2*R - 64*A*R.

Default checks prove the relative monic presentations for k <= 4 and search
the maximally collided fibers through k=7.  Larger bounds are regressions,
not an all-k proof.
"""
from __future__ import annotations

import argparse
import itertools
import math
import time
from dataclasses import dataclass

import sympy as sp


Z = sp.Symbol("Z")


@dataclass(frozen=True)
class Presentation:
    k: int
    base: tuple[sp.Symbol, ...]
    transverse: tuple[sp.Symbol, ...]
    S: sp.Expr
    A: sp.Expr
    Q: sp.Expr
    R: sp.Expr
    equations: tuple[sp.Expr, ...]


def direct_presentation(k: int, base_values: tuple[sp.Expr, ...] | None = None) -> Presentation:
    """Return the exact U-eliminated presentation for a fixed k."""
    base = sp.symbols(f"s1:{k + 1}")
    transverse = sp.symbols(f"a0:{k}")
    if base_values is None:
        base_values = base
    if len(base_values) != k:
        raise ValueError("base_values must have length k")

    S = Z**k + sum(base_values[index - 1] * Z ** (k - index) for index in range(1, k + 1))
    A = sum(transverse[index] * Z**index for index in range(k))
    quotient, remainder = sp.div(sp.Poly(A**2, Z), sp.Poly(S, Z))
    Q = sp.expand(quotient.as_expr())
    R = sp.expand(remainder.as_expr())
    obstruction = sp.Poly(
        sp.expand(9 * Q**2 + 8 * S * A * Q - 48 * S**2 * R - 64 * A * R),
        Z,
    )
    equations = tuple(
        sp.expand(coefficient)
        for coefficient in obstruction.all_coeffs()
        if coefficient != 0
    )
    return Presentation(k, base, transverse, S, A, Q, R, equations)


def verify_elimination_identity(presentation: Presentation) -> None:
    """Check U^2-V^3 equals the compact obstruction divided by 64."""
    S, A, Q, R = (
        presentation.S,
        presentation.A,
        presentation.Q,
        presentation.R,
    )
    assert sp.expand(A**2 - S * Q - R) == 0
    U = S**3 + sp.Rational(3, 2) * S * A + sp.Rational(3, 8) * Q
    V = S**2 + A
    compact = 9 * Q**2 + 8 * S * A * Q - 48 * S**2 * R - 64 * A * R
    assert sp.expand(64 * (U**2 - V**3) - compact) == 0


def relative_groebner(presentation: Presentation) -> sp.GroebnerBasis:
    domain = sp.QQ.poly_ring(*presentation.base)
    return sp.groebner(
        presentation.equations,
        *presentation.transverse,
        order="grevlex",
        domain=domain,
    )


def collided_groebner(k: int) -> tuple[Presentation, sp.GroebnerBasis]:
    presentation = direct_presentation(k, (sp.Integer(0),) * k)
    basis = sp.groebner(
        presentation.equations,
        *presentation.transverse,
        order="grevlex",
        domain=sp.QQ,
    )
    return presentation, basis


def leading_exponents(basis: sp.GroebnerBasis) -> tuple[tuple[int, ...], ...]:
    return tuple(polynomial.LM(order=basis.order).exponents for polynomial in basis.polys)


def standard_exponents(basis: sp.GroebnerBasis) -> tuple[tuple[int, ...], ...]:
    """Enumerate a zero-dimensional monomial quotient by closure from 1."""
    variable_count = len(basis.gens)
    leading = leading_exponents(basis)
    zero = (0,) * variable_count
    seen = {zero}
    stack = [zero]
    while stack:
        exponent = stack.pop()
        for index in range(variable_count):
            neighbor = list(exponent)
            neighbor[index] += 1
            neighbor_tuple = tuple(neighbor)
            if neighbor_tuple in seen:
                continue
            if any(
                all(neighbor_tuple[position] >= monomial[position]
                    for position in range(variable_count))
                for monomial in leading
            ):
                continue
            seen.add(neighbor_tuple)
            stack.append(neighbor_tuple)
            if len(seen) > 1_000_000:
                raise RuntimeError("standard-monomial search exceeded safety bound")
    return tuple(sorted(seen, key=lambda exponent: (sum(exponent), exponent)))


def monomial(variables: tuple[sp.Symbol, ...], exponent: tuple[int, ...]) -> sp.Expr:
    return sp.prod(variable**power for variable, power in zip(variables, exponent))


def quotient_vectorizer(
    basis: sp.GroebnerBasis,
    standard: tuple[tuple[int, ...], ...],
):
    variables = basis.gens
    index = {exponent: position for position, exponent in enumerate(standard)}

    def vector(expression: sp.Expr) -> sp.Matrix:
        result = sp.zeros(len(standard), 1)
        remainder = sp.Poly(basis.reduce(expression)[1], *variables)
        for exponent, coefficient in remainder.terms():
            result[index[exponent]] = coefficient
        return result

    return vector


def maximal_ideal_filtration(
    basis: sp.GroebnerBasis,
    standard: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], list[sp.Matrix], list[sp.Matrix]]:
    """Return Hilbert function, multiplication matrices, and ideal powers."""
    variables = basis.gens
    standard_monomials = [monomial(variables, exponent) for exponent in standard]
    vector = quotient_vectorizer(basis, standard)
    multiplication = [
        sp.Matrix.hstack(*[vector(variable * item) for item in standard_monomials])
        for variable in variables
    ]
    current = sp.Matrix.hstack(*multiplication)
    powers: list[sp.Matrix] = []
    dimensions: list[int] = []
    while True:
        columns = current.columnspace()
        current = sp.Matrix.hstack(*columns) if columns else sp.zeros(len(standard), 0)
        powers.append(current)
        dimensions.append(current.rank())
        if not dimensions[-1]:
            break
        current = sp.Matrix.hstack(*[
            matrix * column
            for matrix in multiplication
            for column in columns
        ])
    hilbert = (len(standard) - dimensions[0],) + tuple(
        dimensions[index] - dimensions[index + 1]
        for index in range(len(dimensions) - 1)
    )
    return hilbert, multiplication, powers


# Filtration-compatible standard-monomial bases.  Variables a_i are the
# coefficients of Z^i in A, so a_(k-1) is the highest transverse coefficient.
CANDIDATE_EXPONENTS: dict[int, tuple[tuple[int, ...], ...]] = {
    1: ((0,), (1,)),
    2: ((0, 0), (0, 1), (1, 0), (0, 2)),
    3: (
        (0, 0, 0),
        (0, 0, 1), (0, 1, 0), (1, 0, 0),
        (0, 0, 2), (0, 1, 1), (1, 0, 1),
        (2, 0, 0),
    ),
}


def filtered_candidate_k4() -> tuple[sp.Expr, ...]:
    """Return the actual k=4 candidate; placeholders above are not used."""
    a0, a1, a2, a3 = sp.symbols("a0:4")
    return (
        sp.Integer(1),
        a3, a2, a1, a0,
        a3**2, a2 * a3, a2**2, a1 * a3, a0 * a3, a0 * a2,
        a0 * a1, a0**2, a1 * a3**2, a0 * a3**2,
        a0 * a1 * a3,
    )


def verify_relative_cases(maximum: int, print_relations: int | None) -> None:
    for k in range(1, maximum + 1):
        presentation = direct_presentation(k)
        verify_elimination_identity(presentation)
        started = time.monotonic()
        basis = relative_groebner(presentation)
        standard = standard_exponents(basis)
        assert all(polynomial.LC(order=basis.order) == 1 for polynomial in basis.polys)
        assert len(standard) == 2**k

        if k < 4:
            expected = CANDIDATE_EXPONENTS[k]
        elif k == 4:
            expected = tuple(
                sp.Poly(item, *presentation.transverse).monoms()[0]
                for item in filtered_candidate_k4()
            )
        else:
            expected = standard
        # The filtered list is a reordering of the relative standard basis.
        assert set(expected) == set(standard)

        special_presentation, special_basis = collided_groebner(k)
        special_standard = standard_exponents(special_basis)
        hilbert, _, powers = maximal_ideal_filtration(special_basis, special_standard)
        assert hilbert == tuple(math.comb(k, degree) for degree in range(k + 1))

        if k <= 4:
            candidate_expressions = [
                monomial(presentation.transverse, exponent) for exponent in expected
            ]
            vector = quotient_vectorizer(special_basis, special_standard)
            valuations = []
            for expression in candidate_expressions:
                value = vector(expression)
                valuation = 0
                for degree, power in enumerate(powers, 1):
                    if power.row_join(value).rank() == power.rank():
                        valuation = degree
                valuations.append(valuation)
            assert sorted(valuations) == sorted(
                degree
                for degree in range(k + 1)
                for _ in range(math.comb(k, degree))
            )

        elapsed = time.monotonic() - started
        print(
            f"PASS relative k={k}: {len(basis.polys)} monic relations, "
            f"rank={len(standard)}, filtration={hilbert}, seconds={elapsed:.2f}"
        )
        if print_relations == k:
            for index, polynomial in enumerate(basis.polys, 1):
                print(f"RELATION {index}: {polynomial.as_expr()}")


def verify_k2_missing_class_and_induction_obstruction() -> None:
    presentation, basis = collided_groebner(2)
    a0, a1 = presentation.transverse
    assert basis.reduce(a1**2)[1] == a1**2
    assert basis.reduce(a1**3)[1] == 0
    assert basis.reduce(a0 * a1)[1] == 0
    assert basis.reduce(a0**2)[1] == 0

    # Every square-zero element in the maximal ideal is b*a0+c*a1^2.
    # Its multiplication operator has rank at most one.  A free quadratic
    # extension of Q[eta]/eta^2 would make multiplication by eta have rank 2.
    standard = standard_exponents(basis)
    vector = quotient_vectorizer(basis, standard)
    standard_monomials = [monomial(presentation.transverse, exponent) for exponent in standard]
    for square_zero in (a0, a1**2, a0 + a1**2):
        assert basis.reduce(square_zero**2)[1] == 0
        multiplication = sp.Matrix.hstack(*[
            vector(square_zero * item) for item in standard_monomials
        ])
        assert multiplication.rank() <= 1
    print("PASS k=2: the indispensable class a1^2 survives")
    print("PASS k=1->2: no free monic quadratic-extension induction on collided fibers")


def search_collided(maximum: int, compute_filtration_through: int) -> None:
    for k in range(1, maximum + 1):
        started = time.monotonic()
        _, basis = collided_groebner(k)
        standard = standard_exponents(basis)
        assert len(standard) == 2**k
        message = (
            f"PASS collided k={k}: groebner={len(basis.polys)}, "
            f"length={len(standard)}"
        )
        if k <= compute_filtration_through:
            hilbert, _, _ = maximal_ideal_filtration(basis, standard)
            expected = tuple(math.comb(k, degree) for degree in range(k + 1))
            assert hilbert == expected
            message += f", filtration={hilbert}"
        print(message + f", seconds={time.monotonic() - started:.2f}")


def standard_count_from_leading(leading: tuple[tuple[int, ...], ...]) -> int:
    variable_count = len(leading[0])
    zero = (0,) * variable_count
    seen = {zero}
    stack = [zero]
    while stack:
        exponent = stack.pop()
        for index in range(variable_count):
            neighbor = list(exponent)
            neighbor[index] += 1
            neighbor = tuple(neighbor)
            if neighbor in seen or any(
                all(neighbor[position] >= monomial[position]
                    for position in range(variable_count))
                for monomial in leading
            ):
                continue
            seen.add(neighbor)
            stack.append(neighbor)
    return len(seen)


def search_artin_directions(k: int, order: int, selected_direction: int | None) -> None:
    """Test coordinate Artin directions S=Z^k+tZ^(k-j), t^order=0."""
    t = sp.Symbol("tau")
    directions = (selected_direction,) if selected_direction is not None else range(1, k + 1)
    for direction in directions:
        base_values = tuple(t if index == direction else sp.Integer(0)
                            for index in range(1, k + 1))
        presentation = direct_presentation(k, base_values)
        variables = (t,) + presentation.transverse
        basis = sp.groebner(
            presentation.equations + (t**order,),
            *variables,
            order="grevlex",
            domain=sp.QQ,
        )
        length = standard_count_from_leading(leading_exponents(basis))
        assert length == order * 2**k
        print(
            f"PASS Artin k={k}, order={order}, base direction={direction}: "
            f"length={length}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relative-max", type=int, default=4, choices=range(1, 6))
    parser.add_argument("--collision-max", type=int, default=7, choices=range(1, 10))
    parser.add_argument("--filtration-max", type=int, default=7, choices=range(0, 8))
    parser.add_argument("--artin-k", type=int, choices=range(1, 8))
    parser.add_argument("--artin-order", type=int, default=2, choices=(2, 3))
    parser.add_argument("--artin-direction", type=int)
    parser.add_argument("--print-relations", type=int, choices=range(1, 6))
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    verify_relative_cases(arguments.relative_max, arguments.print_relations)
    verify_k2_missing_class_and_induction_obstruction()
    search_collided(arguments.collision_max, arguments.filtration_max)
    if arguments.artin_k is not None:
        if arguments.artin_direction is not None and not (
            1 <= arguments.artin_direction <= arguments.artin_k
        ):
            raise ValueError("artin direction must lie between 1 and artin-k")
        search_artin_directions(
            arguments.artin_k,
            arguments.artin_order,
            arguments.artin_direction,
        )
    print("RESULT: no all-k proof and no counterexample in the tested range")


if __name__ == "__main__":
    main()
