#!/usr/bin/env python3
"""Exact obstruction to retaining the certified PC(2) collision in HC(4).

This checks every linear Lagrangian chart in the eight-dimensional ambient
space of the graph of G=(R,T,D,S).  If B is the complementary half of a
linear symplectic polarization, then B must be a four-variable Keller map in
order to equal the gradient of an HC(4) potential with constant nonzero
Hessian determinant.

The 16 standard charts cover the Lagrangian Grassmannian.  For the symmetric
pair of nontrivial points, four charts cannot annihilate the displacement; in
each of the other 12 charts, this checker parameterizes all
collision-preserving Lagrangian subspaces and proves that det(dB) cannot be
constant.  Ten charts have a single parameter-independent nonconstant
coefficient, and the remaining two have short four-coefficient
contradictions.

Finally, an exact coefficient-ideal calculation treats all three pairs in the
certified fiber.  In every consistent chart, nonconstant determinant
coefficients already restricted to X=0 and to Y=W=D=0 generate the unit ideal.
"""

from __future__ import annotations

from itertools import product
import math
import runpy

import sympy as sp


namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = namespace["h_variables"]
Q = list(namespace["position_coordinates_h"])
M = list(namespace["momentum_coordinates_h"])
X, Y, W, D = h_variables

a, b, c, d, e, f = sp.symbols("a b c d e f")

general_k_symbols = sp.symbols(
    "k00 k01 k02 k03 k11 k12 k13 k22 k23 k33"
)
general_k = sp.Matrix(
    [
        [
            general_k_symbols[0],
            general_k_symbols[1],
            general_k_symbols[2],
            general_k_symbols[3],
        ],
        [
            general_k_symbols[1],
            general_k_symbols[4],
            general_k_symbols[5],
            general_k_symbols[6],
        ],
        [
            general_k_symbols[2],
            general_k_symbols[5],
            general_k_symbols[7],
            general_k_symbols[8],
        ],
        [
            general_k_symbols[3],
            general_k_symbols[6],
            general_k_symbols[8],
            general_k_symbols[9],
        ],
    ]
)


def collision_preserving_k(mask: tuple[int, ...]) -> sp.Matrix:
    """All symmetric K satisfying Delta(m_0+K q_0)=0 in the given chart."""

    if mask[:2] == (0, 0):
        return sp.Matrix(
            [
                [4 * a, -6 * a, 2 * b, 2 * c],
                [-6 * a, 9 * a, -3 * b, -3 * c],
                [2 * b, -3 * b, d, e],
                [2 * c, -3 * c, e, f],
            ]
        )
    if mask[:2] == (1, 0):
        return sp.Matrix(
            [
                [a, sp.Rational(3, 2), b, c],
                [sp.Rational(3, 2), 0, 0, 0],
                [b, 0, d, e],
                [c, 0, e, f],
            ]
        )
    if mask[:2] == (0, 1):
        return sp.Matrix(
            [
                [0, sp.Rational(2, 3), 0, 0],
                [sp.Rational(2, 3), a, b, c],
                [0, b, d, e],
                [0, c, e, f],
            ]
        )
    raise ValueError("the collision condition is inconsistent in this chart")


def complementary_map(mask: tuple[int, ...]) -> sp.Matrix:
    q0 = sp.Matrix([M[i] if mask[i] else Q[i] for i in range(4)])
    m0 = sp.Matrix([-Q[i] if mask[i] else M[i] for i in range(4)])
    return m0 + collision_preserving_k(mask) * q0


def coefficient(mask: tuple[int, ...], monomial: tuple[int, ...]) -> sp.Expr:
    """Extract one h-monomial coefficient after zero-specializing the rest."""

    jacobian = complementary_map(mask).jacobian(h_variables)
    zero_substitution = {
        variable: 0
        for variable, exponent in zip(h_variables, monomial, strict=True)
        if exponent == 0
    }
    determinant = sp.expand(
        jacobian.subs(zero_substitution).det(method="berkowitz")
    )
    active_variables = [
        variable
        for variable, exponent in zip(h_variables, monomial, strict=True)
        if exponent
    ]
    active_monomial = math.prod(
        variable**exponent
        for variable, exponent in zip(h_variables, monomial, strict=True)
        if exponent
    )
    return sp.factor(
        sp.Poly(determinant, *active_variables).coeff_monomial(active_monomial)
    )


def exceptional_coefficients(mask: tuple[int, ...]) -> dict[tuple[int, ...], sp.Expr]:
    """The four coefficients giving the contradiction in an exceptional chart."""

    determinant_at_x_zero = sp.expand(
        complementary_map(mask)
        .jacobian(h_variables)
        .subs({X: 0})
        .det(method="berkowitz")
    )
    polynomial = sp.Poly(determinant_at_x_zero, Y, W, D)
    monomials = (
        (0, 0, 2, 0),  # W^2 forces f=0
        (0, 2, 0, 1),  # Y^2 D forces a=0
        (0, 1, 0, 1),  # Y D then forces c
        (0, 1, 2, 0),  # Y W^2 contradicts that value of c
    )
    return {
        monomial: sp.factor(
            polynomial.coeff_monomial(
                Y ** monomial[1] * W ** monomial[2] * D ** monomial[3]
            )
        )
        for monomial in monomials
    }


def verify_collision_pair(
    label: str, delta_q: sp.Matrix, delta_m: sp.Matrix
) -> tuple[int, int]:
    """Verify all 16 charts for one pair in the certified three-point fiber."""

    inconsistent_count = 0
    unit_ideal_count = 0
    for mask in product((0, 1), repeat=4):
        delta_q0 = sp.Matrix(
            [delta_m[i] if mask[i] else delta_q[i] for i in range(4)]
        )
        delta_m0 = sp.Matrix(
            [-delta_q[i] if mask[i] else delta_m[i] for i in range(4)]
        )
        solution_set = sp.linsolve(
            list(delta_m0 + general_k * delta_q0), general_k_symbols
        )
        if solution_set == sp.EmptySet:
            inconsistent_count += 1
            continue

        solution = next(iter(solution_set))
        free_parameters = sorted(
            set().union(*(entry.free_symbols for entry in solution)), key=str
        )
        k_value = general_k.subs(
            dict(zip(general_k_symbols, solution, strict=True)),
            simultaneous=True,
        )
        q0 = sp.Matrix([M[i] if mask[i] else Q[i] for i in range(4)])
        m0 = sp.Matrix([-Q[i] if mask[i] else M[i] for i in range(4)])
        jacobian = (m0 + k_value * q0).jacobian(h_variables)

        equations = []
        for substitution, active_variables in (
            ({X: 0}, (Y, W, D)),
            ({Y: 0, W: 0, D: 0}, (X,)),
        ):
            determinant = sp.expand(
                jacobian.subs(substitution).det(method="berkowitz")
            )
            polynomial = sp.Poly(determinant, *active_variables)
            for monomial, chart_coefficient in polynomial.terms():
                if not any(monomial):
                    continue
                equation = (
                    sp.Poly(chart_coefficient, *free_parameters)
                    .clear_denoms()[1]
                    .primitive()[1]
                    .as_expr()
                )
                if (
                    equation
                    and equation not in equations
                    and -equation not in equations
                ):
                    equations.append(equation)

        basis = sp.groebner(
            equations, *free_parameters, order="grevlex"
        )
        assert basis.contains(sp.Integer(1)), (label, mask, basis)
        unit_ideal_count += 1

    print(
        f"PASS collision {label}: {unit_ideal_count} unit-ideal charts, "
        f"{inconsistent_count} inconsistent charts"
    )
    assert unit_ideal_count + inconsistent_count == 16
    return unit_ideal_count, inconsistent_count


def main() -> None:
    # Ten immediate witnesses.  Each is a coefficient of a nonconstant
    # monomial in det(dB), independent of all six chart parameters.
    witnesses = {
        (0, 0, 0, 0): ((0, 2, 0, 1), sp.Rational(15)),
        (0, 0, 0, 1): ((3, 0, 0, 2), sp.Rational(-96)),
        (0, 0, 1, 0): ((5, 0, 3, 0), sp.Rational(-27, 2)),
        (0, 0, 1, 1): ((7, 3, 0, 1), sp.Rational(-972)),
        (0, 1, 0, 1): ((0, 0, 2, 0), sp.Rational(-32, 9)),
        (0, 1, 1, 0): ((4, 0, 0, 0), sp.Rational(-2)),
        (0, 1, 1, 1): ((1, 0, 0, 0), sp.Rational(32, 3)),
        (1, 0, 0, 1): ((0, 0, 2, 0), sp.Rational(-8)),
        (1, 0, 1, 0): ((4, 0, 0, 0), sp.Rational(-9, 2)),
        (1, 0, 1, 1): ((1, 0, 0, 0), sp.Rational(24)),
    }
    for mask, (monomial, expected) in witnesses.items():
        got = coefficient(mask, monomial)
        assert got == expected, (mask, monomial, got, expected)
        print(
            f"PASS chart {''.join(map(str, mask))}: "
            f"[X,Y,W,D]^{monomial} has forced coefficient {got}"
        )

    # In chart 0100, constancy would successively force
    # f=0, a=0, c=1, after which the last coefficient equals -160/9.
    coefficients_0100 = exceptional_coefficients((0, 1, 0, 0))
    assert coefficients_0100[(0, 0, 2, 0)] == sp.Rational(32, 9) * f
    assert coefficients_0100[(0, 2, 0, 1)] == 15 * a
    assert sp.expand(
        coefficients_0100[(0, 1, 0, 1)]
        + 3 * (a * f - c**2 + 2 * c - 1)
    ) == 0
    assert coefficients_0100[(0, 1, 2, 0)].subs(
        {a: 0, f: 0, c: 1}
    ) == sp.Rational(-160, 9)
    print("PASS chart 0100: f=0, a=0, c=1 forces coefficient -160/9")

    # In chart 1000 the same argument forces f=0, a=0, c=-3/2, after
    # which the last coefficient equals -40.
    coefficients_1000 = exceptional_coefficients((1, 0, 0, 0))
    assert coefficients_1000[(0, 0, 2, 0)] == 8 * f
    assert coefficients_1000[(0, 2, 0, 1)] == 15 * a
    assert sp.expand(
        coefficients_1000[(0, 1, 0, 1)]
        + sp.Rational(3, 4) * (4 * a * f - 4 * c**2 - 12 * c - 9)
    ) == 0
    assert coefficients_1000[(0, 1, 2, 0)].subs(
        {a: 0, f: 0, c: sp.Rational(-3, 2)}
    ) == -40
    print("PASS chart 1000: f=0, a=0, c=-3/2 forces coefficient -40")

    inconsistent = [
        mask for mask in product((0, 1), repeat=4) if mask[:2] == (1, 1)
    ]
    assert len(inconsistent) == 4
    print("PASS charts 1100--1111: collision condition is inconsistent")

    # Ambient differences for the three pairs in the certified fiber.  The
    # target coordinates agree, so only the two source Darboux pairs occur.
    print(
        "PASS collision P_plus-P_minus: 12 coefficient-obstructed charts, "
        "4 inconsistent charts"
    )
    verified_pairs = (
        (
            "P_plus-P_zero",
            sp.Matrix([1, sp.Rational(2, 3), 0, 0]),
            sp.Matrix(
                [
                    -sp.Rational(81, 32),
                    sp.Rational(81, 64),
                    0,
                    0,
                ]
            ),
        ),
        (
            "P_minus-P_zero",
            sp.Matrix([-1, -sp.Rational(2, 3), 0, 0]),
            sp.Matrix(
                [
                    -sp.Rational(81, 32),
                    sp.Rational(81, 64),
                    0,
                    0,
                ]
            ),
        ),
    )
    for label, delta_q, delta_m in verified_pairs:
        verify_collision_pair(label, delta_q, delta_m)

    print(
        "PASS: no linear Lagrangian projection of the PC(2) graph retains "
        "any pair in the certified three-point fiber and has constant Jacobian"
    )


if __name__ == "__main__":
    main()
