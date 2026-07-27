#!/usr/bin/env python3
"""Prove that no linear Lagrangian projection of the PC(2) graph is Keller.

In each of the 16 standard Lagrangian charts, every linear complementary
projection has the form B=m_0+K q_0 with K symmetric.  Constancy of det(dB)
forces every nonconstant coefficient to vanish.  Exact coefficient ideals on
the three hyperplane slices X=0, Y=0, and W=0 already give:

* the unit ideal in 14 charts; and
* in charts 1110 and 1111, only K=diag(0,0,0,l), whose projection has
  identically zero determinant.

Thus no linear symplectic polarization can produce an HC(4) gradient,
regardless of which collision it might have.
"""

from __future__ import annotations

from itertools import product
import runpy

import sympy as sp


namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = namespace["h_variables"]
X, Y, W, D = h_variables
Q = list(namespace["position_coordinates_h"])
M = list(namespace["momentum_coordinates_h"])

a, b, c, d, e, f, g, iota, j, ell = sp.symbols(
    "a b c d e f g iota j ell"
)
parameters = (a, b, c, d, e, f, g, iota, j, ell)
K = sp.Matrix(
    [
        [a, b, c, d],
        [b, e, f, g],
        [c, f, iota, j],
        [d, g, j, ell],
    ]
)


def add_slice_equations(
    equations: list[sp.Expr],
    jacobian: sp.Matrix,
    substitution: dict[sp.Symbol, int],
    active_variables: tuple[sp.Symbol, ...],
) -> None:
    determinant = sp.expand(
        jacobian.subs(substitution).det(method="berkowitz")
    )
    polynomial = sp.Poly(determinant, *active_variables)
    for monomial, coefficient in polynomial.terms():
        if not any(monomial):
            continue
        equation = (
            sp.Poly(coefficient, *parameters)
            .clear_denoms()[1]
            .primitive()[1]
            .as_expr()
        )
        if equation and equation not in equations and -equation not in equations:
            equations.append(equation)


def main() -> None:
    exceptional = {(1, 1, 1, 0), (1, 1, 1, 1)}
    unit_charts = 0
    singular_charts = 0

    for mask in product((0, 1), repeat=4):
        q0 = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
        m0 = sp.Matrix(
            [-Q[index] if mask[index] else M[index] for index in range(4)]
        )
        complementary = m0 + K * q0
        jacobian = complementary.jacobian(h_variables)
        equations: list[sp.Expr] = []
        basis = None

        for substitution, active_variables in (
            ({X: 0}, (Y, W, D)),
            ({Y: 0}, (X, W, D)),
            ({W: 0}, (X, Y, D)),
        ):
            add_slice_equations(
                equations, jacobian, substitution, active_variables
            )
            basis = sp.groebner(equations, *parameters, order="grevlex")
            if basis.contains(sp.Integer(1)):
                break

        assert basis is not None
        label = "".join(map(str, mask))
        if mask not in exceptional:
            assert basis.contains(sp.Integer(1)), (mask, basis)
            unit_charts += 1
            print(f"PASS chart {label}: sliced coefficient ideal is the unit ideal")
            continue

        assert not basis.contains(sp.Integer(1))
        basis_expressions = {polynomial.as_expr() for polynomial in basis.polys}
        forced = {
            d**2,
            d * g,
            g**2,
            d * j,
            g * j,
            j**2,
            a,
            b,
            c,
            e,
            f,
            iota,
        }
        assert basis_expressions == forced, (mask, basis_expressions)

        singular_k = sp.zeros(4)
        singular_k[3, 3] = ell
        singular_projection = m0 + singular_k * q0
        singular_determinant = sp.expand(
            singular_projection.jacobian(h_variables).det(method="berkowitz")
        )
        assert singular_determinant == 0
        singular_charts += 1
        print(
            f"PASS chart {label}: constancy forces K=diag(0,0,0,ell), "
            "and det(dB)=0"
        )

    assert unit_charts == 14
    assert singular_charts == 2
    print(
        "PASS: no linear Lagrangian projection of the PC(2) graph has "
        "constant nonzero Jacobian"
    )


if __name__ == "__main__":
    main()
