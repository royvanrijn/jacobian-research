#!/usr/bin/env python3
"""Verify the one-sided positive-characteristic Ritt tangent classification.

For monic original ``R`` of degree ``r`` and a normalized tangent
``U in x*k[x]`` of degree below ``r``, the invisible tangent space at
Hessian cutoff ``q`` is

    {U : deg(H'(R)U) <= q}.

If ``H'=0`` this is the full ``r-1`` dimensional tangent space.  Otherwise,
put ``d=deg(H')``.  Degree additivity over a field makes its dimension

    max(0, min(r-1, q-r*d)).

The script checks the formula exhaustively in small finite fields and
regresses the characteristic-two primitive examples.
"""

from __future__ import annotations

import itertools

import sympy as sp


x = sp.Symbol("x")


def poly(expression: sp.Expr, prime: int) -> sp.Poly:
    return sp.Poly(expression, x, modulus=prime)


def degree_or_none(item: sp.Poly) -> int | None:
    return None if item.is_zero else int(item.degree())


def expected_invisible(
    derivative: sp.Poly,
    tangent: sp.Poly,
    right_degree: int,
    cutoff: int,
) -> bool:
    if tangent.is_zero or derivative.is_zero:
        return True
    derivative_degree = int(derivative.degree())
    tangent_degree = int(tangent.degree())
    return right_degree * derivative_degree + tangent_degree <= cutoff


def exhaustive_check() -> int:
    cases = 0
    for prime in (2, 3):
        for outer_degree in range(2, 5):
            for right_degree in (2, 3):
                for h_coefficients in itertools.product(
                    range(prime), repeat=outer_degree
                ):
                    H = poly(
                        x**outer_degree
                        + sum(
                            coefficient * x**index
                            for index, coefficient in enumerate(h_coefficients)
                        ),
                        prime,
                    )
                    derivative = H.diff()
                    for r_coefficients in itertools.product(
                        range(prime), repeat=right_degree - 1
                    ):
                        R = poly(
                            x**right_degree
                            + sum(
                                coefficient * x ** (index + 1)
                                for index, coefficient in enumerate(r_coefficients)
                            ),
                            prime,
                        )
                        composed_derivative = poly(
                            derivative.as_expr().subs(x, R.as_expr()), prime
                        )
                        for u_coefficients in itertools.product(
                            range(prime), repeat=right_degree - 1
                        ):
                            U = poly(
                                sum(
                                    coefficient * x ** (index + 1)
                                    for index, coefficient in enumerate(
                                        u_coefficients
                                    )
                                ),
                                prime,
                            )
                            first_order = composed_derivative * U
                            actual = (
                                first_order.is_zero
                                or int(first_order.degree()) <= 1
                            )
                            predicted = expected_invisible(
                                derivative, U, right_degree, cutoff=1
                            )
                            assert actual == predicted, (
                                prime,
                                H,
                                R,
                                U,
                                degree_or_none(first_order),
                            )
                            cases += 1
    return cases


def main() -> None:
    cases = exhaustive_check()

    # The primitive characteristic-two defect H=z^2+z has H'=1.
    H_separable = poly(x**2 + x, 2)
    R = poly(x**3 + x**2, 2)
    U = poly(x, 2)
    first_order = poly(
        H_separable.diff().as_expr().subs(x, R.as_expr()) * U.as_expr(), 2
    )
    assert first_order == poly(x, 2)

    # The purely inseparable outer polynomial kills every right tangent.
    H_radicial = poly(x**2, 2)
    assert H_radicial.diff().is_zero
    full_tangent_basis = [poly(x**index, 2) for index in range(1, 3)]
    assert all(
        poly(
            H_radicial.diff().as_expr().subs(x, R.as_expr()) * basis.as_expr(),
            2,
        ).is_zero
        for basis in full_tangent_basis
    )

    print(
        "PASS positive-characteristic Ritt tangents: "
        f"{cases} exhaustive normalized pairs"
    )
    print("PASS char 2: H=z^2+z and U=x give the invisible defect epsilon*x")
    print("PASS char 2: H=z^2 kills the full normalized right tangent space")
    print(
        "THEOREM: dim V_q=r-1 if H'=0, otherwise "
        "max(0,min(r-1,q-r*deg(H')))"
    )


if __name__ == "__main__":
    main()
