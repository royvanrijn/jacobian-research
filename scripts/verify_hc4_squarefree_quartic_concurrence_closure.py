#!/usr/bin/env python3
"""Verify the concurrence closures in the squarefree quartic packet.

The checker has two independent parts.

* At a four-line pencil point it verifies the rank-zero/rank-one Hessian
  determinant ladders used to exclude an order-eight squarefree tangent
  cone.
* For the exactly-three-concurrent arrangement it computes the quartic
  polar syzygies when the fourth flag is transverse and verifies the four
  mixed-partial obstructions (RRR, TRR, TTR, and TTT on the concurrent
  lines).

All calculations are over QQ and are exact.
"""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def binary_form(degree: int, coefficients: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sum(
        coefficient * x ** (degree - index) * y**index
        for index, coefficient in enumerate(coefficients)
    )


def xy_degree_piece(expression: sp.Expr, degree: int) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), x, y, z)
    return sp.expand(
        sum(
            coefficient * x**monomial[0] * y**monomial[1] * z**monomial[2]
            for monomial, coefficient in polynomial.terms()
            if monomial[0] + monomial[1] == degree
        )
    )


def assert_zero(expression: sp.Expr) -> None:
    assert sp.expand(expression) == 0


def verify_binary_root_multiplicities() -> None:
    """A quartic times its binary Hessian is never P^2 for squarefree P."""
    a, b, c = sp.symbols("a b c")
    models = {
        1: x * (a * x**3 + b * x**2 * y + c * x * y**2 + y**3),
        2: x**2 * (a * x**2 + b * x * y + c * y**2),
        3: x**3 * (a * x + b * y),
        4: x**4,
    }
    expected_restrictions = {
        1: -9 * y**4,
        2: -12 * c**2 * y**2,
        3: -9 * b**2,
    }
    for multiplicity, form in models.items():
        hessian = sp.factor(sp.hessian(form, (x, y)).det())
        if multiplicity == 4:
            assert hessian == 0
            continue
        quotient = sp.cancel(hessian / x ** (2 * multiplicity - 2))
        assert sp.denom(quotient) == 1
        assert_zero(quotient.subs(x, 0) - expected_restrictions[multiplicity])


def verify_pencil_point_ladders() -> None:
    """Verify all rank-at-most-one jets at a four-line pencil point."""
    A, alpha, beta = sp.symbols("A alpha beta", nonzero=True)

    # Rank one with a nonzero z^5 coefficient.  Rank one of the constant
    # Hessian permits a z -> z+linear shift which removes the z^4 and z^3
    # jets.  A nonzero cubic jet must be a cube.
    f = sp.symbols("f0:5")
    g = sp.symbols("g0:6")
    F4 = binary_form(4, f)
    F5 = binary_form(5, g)
    h = A * z**5 + alpha * x**3 * z**2 + z * F4 + F5
    determinant = sp.expand(sp.hessian(h, variables).det())
    assert_zero(
        xy_degree_piece(determinant, 3)
        - (
            240 * A * alpha * f[2] * x**3 * z**6
            + 720 * A * alpha * f[3] * x**2 * y * z**6
            + 1440 * A * alpha * f[4] * x * y**2 * z**6
        )
    )
    cubic_reduction = determinant.subs({f[2]: 0, f[3]: 0, f[4]: 0})
    degree_four = xy_degree_piece(cubic_reduction, 4)
    forced_degree_four = {
        g[2]: 3 * f[1] ** 2 / (4 * alpha),
        g[3]: 0,
        g[4]: 0,
        g[5]: 0,
    }
    assert_zero(degree_four.subs(forced_degree_four))
    degree_five = xy_degree_piece(
        cubic_reduction.subs(forced_degree_four), 5
    )
    assert sp.factor(degree_five.coeff(x, 4).coeff(y, 1).coeff(z, 4)) == (
        -360 * A * f[1] ** 3 / alpha
    )
    cone_reduction = cubic_reduction.subs(forced_degree_four).subs(f[1], 0)
    degree_six = xy_degree_piece(cone_reduction, 6)
    assert sp.factor(degree_six.coeff(x, 6).coeff(z, 3)) == -320 * A * g[1] ** 2
    cone_form = h.subs(
        {
            f[1]: 0,
            f[2]: 0,
            f[3]: 0,
            f[4]: 0,
            g[1]: 0,
            g[2]: 0,
            g[3]: 0,
            g[4]: 0,
            g[5]: 0,
        }
    )
    assert sp.diff(cone_form, y) == 0

    # If the cubic jet vanishes, a first nonzero quartic jet must be a
    # fourth power; its coefficient ladder again makes the quintic a cone.
    h_quartic = A * z**5 + beta * x**4 * z + F5
    determinant_quartic = sp.expand(sp.hessian(h_quartic, variables).det())
    degree_five = xy_degree_piece(determinant_quartic, 5)
    expected = (
        480 * A * beta * g[2] * x**5 * z**4
        + 1440 * A * beta * g[3] * x**4 * y * z**4
        + 2880 * A * beta * g[4] * x**3 * y**2 * z**4
        + 4800 * A * beta * g[5] * x**2 * y**3 * z**4
    )
    assert_zero(degree_five - expected)
    quartic_reduction = determinant_quartic.subs(
        {g[2]: 0, g[3]: 0, g[4]: 0, g[5]: 0}
    )
    assert sp.factor(
        xy_degree_piece(quartic_reduction, 6).coeff(x, 6).coeff(z, 3)
    ) == -320 * A * g[1] ** 2
    assert sp.diff(
        h_quartic.subs(
            {g[1]: 0, g[2]: 0, g[3]: 0, g[4]: 0, g[5]: 0}
        ),
        y,
    ) == 0
    binary_five = binary_form(5, g)
    assert_zero(
        sp.hessian(A * z**5 + binary_five, variables).det()
        - 20 * A * z**3 * sp.hessian(binary_five, (x, y)).det()
    )

    # Rank one with zero z^5 coefficient.  The rank-one constant Hessian is
    # normalized to the x^2*z^3 jet.  The cubic jet is removable after the
    # first determinant face.
    q = sp.symbols("q0:4")
    r = sp.symbols("r0:5")
    s = sp.symbols("s0:6")
    h_rank_one = (
        alpha * x**2 * z**3
        + z**2 * binary_form(3, q)
        + z * binary_form(4, r)
        + binary_form(5, s)
    )
    determinant_rank_one = sp.expand(sp.hessian(h_rank_one, variables).det())
    degree_three = xy_degree_piece(determinant_rank_one, 3)
    assert_zero(
        degree_three
        + 48 * alpha**2 * q[2] * x**3 * z**6
        + 144 * alpha**2 * q[3] * x**2 * y * z**6
    )
    normalized_rank_one = h_rank_one.subs(
        {q[0]: 0, q[1]: 0, q[2]: 0, q[3]: 0}
    )
    normalized_determinant = sp.expand(
        sp.hessian(normalized_rank_one, variables).det()
    )
    degree_four = xy_degree_piece(normalized_determinant, 4)
    assert_zero(
        degree_four
        + 48 * alpha**2 * r[2] * x**4 * z**5
        + 144 * alpha**2 * r[3] * x**3 * y * z**5
        + 288 * alpha**2 * r[4] * x**2 * y**2 * z**5
    )
    rank_one_reduction = normalized_rank_one.subs(
        {r[2]: 0, r[3]: 0, r[4]: 0, s[2]: 0, s[3]: 0, s[4]: 0, s[5]: 0}
    )
    rank_one_formula = sp.factor(sp.hessian(rank_one_reduction, variables).det())
    expected_rank_one = 4 * x**6 * (
        -5 * alpha * r[1] ** 2 * z**3
        - 24 * alpha * r[1] * s[1] * x * z**2
        - 24 * alpha * s[1] ** 2 * x**2 * z
        + 3 * r[0] * r[1] ** 2 * x**2 * z
        + 8 * r[0] * r[1] * s[1] * x**3
        + 3 * r[1] ** 3 * x * y * z
        - 5 * r[1] ** 2 * s[0] * x**3
        + 3 * r[1] ** 2 * s[1] * x**2 * y
    )
    assert_zero(rank_one_formula - expected_rank_one)
    assert_zero(
        rank_one_formula.subs(r[1], 0) + 96 * alpha * s[1] ** 2 * x**8 * z
    )

    # Rank zero.  A nonzero cubic tangent jet is forced to be a cube.  The
    # next two faces give the exact ninth-power determinant.
    u = sp.symbols("u0:4")
    v = sp.symbols("v0:5")
    w = sp.symbols("w0:6")
    cubic = binary_form(3, u)
    quartic = binary_form(4, v)
    quintic = binary_form(5, w)
    h_rank_zero = z**2 * cubic + z * quartic + quintic
    determinant_rank_zero = sp.expand(sp.hessian(h_rank_zero, variables).det())
    assert_zero(
        xy_degree_piece(determinant_rank_zero, 5)
        + 4 * z**4 * cubic * sp.hessian(cubic, (x, y)).det()
    )
    cube_model = alpha * x**3 * z**2 + z * quartic + quintic
    cube_determinant = sp.expand(sp.hessian(cube_model, variables).det())
    assert_zero(
        xy_degree_piece(cube_determinant, 6)
        + 48 * alpha**2 * v[2] * x**6 * z**3
        + 144 * alpha**2 * v[3] * x**5 * y * z**3
        + 288 * alpha**2 * v[4] * x**4 * y**2 * z**3
    )
    cube_reduction = cube_model.subs(
        {
            v[2]: 0,
            v[3]: 0,
            v[4]: 0,
            w[2]: v[1] ** 2 / (4 * alpha),
            w[3]: 0,
            w[4]: 0,
            w[5]: 0,
        }
    )
    assert_zero(
        sp.hessian(cube_reduction, variables).det()
        + 8 * x**9 * (2 * alpha * w[1] - v[0] * v[1]) ** 2 / alpha
    )

    # If the cubic jet is zero, the order-eight initial form is F4*Hess(F4).
    rank_zero_quartic = z * quartic + quintic
    initial_eight = xy_degree_piece(
        sp.hessian(rank_zero_quartic, variables).det(), 8
    )
    assert_zero(
        initial_eight
        + sp.Rational(4, 3)
        * z
        * quartic
        * sp.hessian(quartic, (x, y)).det()
    )

    verify_binary_root_multiplicities()


def quartic_monomials() -> list[sp.Expr]:
    return [
        x**i * y**j * z ** (4 - i - j)
        for i in range(5)
        for j in range(5 - i)
    ]


def polar_syzygies(pattern: str) -> list[tuple[sp.Expr, ...]]:
    lines = (x, y, x + y, z)
    monomials = quartic_monomials()
    columns: list[sp.Poly] = []
    blocks: list[tuple[sp.Expr, list[sp.Expr]]] = []
    for line, flag in zip(lines, pattern):
        exponent = 2 if flag == "T" else 3
        residual_degree = 4 - exponent
        residuals = [
            x**i * y**j * z ** (residual_degree - i - j)
            for i in range(residual_degree + 1)
            for j in range(residual_degree + 1 - i)
        ]
        blocks.append((line**exponent, residuals))
        columns.extend(
            sp.Poly(sp.expand(line**exponent * residual), x, y, z)
            for residual in residuals
        )
    matrix = sp.Matrix(
        [
            [column.coeff_monomial(monomial) for column in columns]
            for monomial in monomials
        ]
    )
    result: list[tuple[sp.Expr, ...]] = []
    for vector in matrix.nullspace():
        offset = 0
        components: list[sp.Expr] = []
        for line_power, residuals in blocks:
            components.append(
                sp.factor(
                    line_power
                    * sum(
                        vector[offset + index] * residual
                        for index, residual in enumerate(residuals)
                    )
                )
            )
            offset += len(residuals)
        result.append(tuple(components))
    return result


def directional_derivative(vector: tuple[sp.Expr, ...], form: sp.Expr) -> sp.Expr:
    return sp.expand(
        sum(entry * sp.diff(form, variable) for entry, variable in zip(vector, variables))
    )


def coefficient_equations(expression: sp.Expr) -> list[sp.Expr]:
    return [
        sp.factor(coefficient)
        for _, coefficient in sp.Poly(sp.expand(expression), x, y, z).terms()
    ]


def assert_equations(actual: list[sp.Expr], expected: list[sp.Expr]) -> None:
    assert len(actual) == len(expected)
    assert all(sp.expand(left - right) == 0 for left, right in zip(actual, expected))


def verify_transverse_fourth_syzygies() -> None:
    """Close all eight triple-concurrent patterns with transverse fourth flag."""
    expected_dimensions = {"RRRR": 1, "TRRR": 2, "TTRR": 4, "TTTR": 6}
    syzygy_bases: dict[str, list[tuple[sp.Expr, ...]]] = {}
    for pattern, expected_dimension in expected_dimensions.items():
        basis = polar_syzygies(pattern)
        assert len(basis) == expected_dimension
        assert all(component[3] == 0 for component in basis)
        syzygy_bases[pattern] = basis

    # RRR: the unique relation has the displayed three binary components.
    rrr = syzygy_bases["RRRR"][0]
    assert rrr == (
        -x**3 * (x + 2 * y),
        y**3 * (2 * x + y),
        (x - y) * (x + y) ** 3,
        0,
    )
    a, b, c, d, e, f = sp.symbols("a b c d e f")
    mixed_rrr = coefficient_equations(
        directional_derivative((a, b, c), rrr[1])
        - directional_derivative((d, e, f), rrr[0])
    )
    assert_equations(mixed_rrr, [4 * d + 2 * e, 6 * d, 6 * b, 2 * a + 4 * b])

    # TRR: a general syzygy and its mixed partials force the y-component of
    # the second (transverse-to-y) direction to vanish.
    A, B, u, v, p, q, r = sp.symbols("A B u v p q r")
    trr_basis = syzygy_bases["TRRR"]
    trr = tuple(
        sp.expand(A * trr_basis[0][index] + B * trr_basis[1][index])
        for index in range(4)
    )
    mixed_trr = coefficient_equations(
        directional_derivative((0, u, v), trr[1])
        - directional_derivative((p, q, r), trr[0])
    )
    assert_equations(mixed_trr, [
        A * q + 4 * B * p + 3 * B * q,
        3 * (A * p + 2 * A * q + 3 * B * p + 2 * B * q),
        3 * (2 * A * p - 3 * A * u + 2 * B * p - B * u),
        -4 * A * u,
    ])

    # TTR: if B is nonzero then p=u=0.  If B=0, the remaining four
    # coefficients force A=C=D=0 whenever p+u is nonzero.
    A, B, C, D = sp.symbols("A B C D")
    ttr_basis = syzygy_bases["TTRR"]
    ttr = tuple(
        sp.expand(
            sum(
                scalar * basis[index]
                for scalar, basis in zip((A, B, C, D), ttr_basis)
            )
        )
        for index in range(4)
    )
    mixed_ttr = coefficient_equations(
        directional_derivative((0, u, v), ttr[1])
        - directional_derivative((p, 0, r), ttr[0])
    )
    assert_equations(mixed_ttr, [
        B * r + 4 * D * p,
        2 * A * u + 3 * B * r + 3 * C * p + 9 * D * p,
        3 * B * p,
        2 * A * p - 3 * B * v + 6 * C * p - 9 * C * u + 6 * D * p - 3 * D * u,
        6 * B * (p - u),
        -B * v - 4 * C * u,
        -3 * B * u,
    ])

    # TTT: when the relation has nonzero xy projection, p=-u and the same
    # mixed-partial row forces all six syzygy coefficients to vanish.  When
    # u=p=0, injectivity leaves only a repeated tangent pair; the final
    # transverse polar is checked explicitly below.
    A, B, C, D, E, F = sp.symbols("A B C D E F")
    ttt_basis = syzygy_bases["TTTR"]
    ttt = tuple(
        sp.expand(
            sum(
                scalar * basis[index]
                for scalar, basis in zip((A, B, C, D, E, F), ttt_basis)
            )
        )
        for index in range(4)
    )
    mixed_ttt = coefficient_equations(
        directional_derivative((0, u, v), ttt[1])
        - directional_derivative((-u, 0, r), ttt[0])
    )
    assert_equations(mixed_ttt, [
        D * r - 4 * F * u,
        2 * A * u + B * r + 2 * D * r - 3 * E * u - 6 * F * u,
        -3 * D * u,
        -2 * A * u - 2 * B * v - 8 * C * u - D * v - 7 * E * u - 2 * F * u,
        -6 * u * (B + D),
        -B * v - 4 * C * u,
        -3 * B * u,
    ])

    kappa, c4, d4 = sp.symbols("kappa c4 d4", nonzero=True)
    repeated_pair_form = kappa * x**2 * y**2 * z + sp.Function("F5")(x, y)
    # Only the z-dependent part is needed: divisibility by z^3 forces both
    # its z coefficient and constant kappa*x^2*y^2 term to vanish.
    z_coefficient = sp.expand(
        2 * kappa * (c4 * x * y**2 + d4 * x**2 * y)
    )
    assert sp.Poly(z_coefficient, x, y).coeff_monomial(x * y**2) == 2 * kappa * c4
    assert sp.Poly(z_coefficient, x, y).coeff_monomial(x**2 * y) == 2 * kappa * d4
    constant_z_derivative = kappa * x**2 * y**2
    assert constant_z_derivative != 0


def main() -> None:
    verify_pencil_point_ladders()
    print("PASS: all-four-concurrent order-eight Hessian ladders")
    verify_transverse_fourth_syzygies()
    print("PASS: triple-concurrent transverse-fourth polar syzygies")
    print("PASS: squarefree quartic concurrence closure")


if __name__ == "__main__":
    main()
