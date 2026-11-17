#!/usr/bin/env python3
"""Verify exact Meng--Yang polynomial-termination frontends.

The all-normal theorem makes the graph determinant recursively solvable in
the source-normal variable ``x``.  This checker supplies seven exact next
steps without promoting a bounded prefix to a termination theorem:

1. the pullback potential is quadratic in the graph and its determinant is
   semilinear in the highest normal derivative ``R_xx``;
2. the leading equation at ``x=infinity`` for a putative terminal
   coefficient ``R_n=U`` is a factored bordered-Hessian equation;
3. the pure weight-five ansatz ``R=(N/L)*y^5*f(x*y)`` has no polynomial
   constant-determinant member, by a dominant-balance coefficient with no
   integral resonance;
4. the terminal equation on the coupled five-function packet has exactly
   two leading-coefficient branches;
5. the complete 2,348-term weight-six equation has a two-chamber upper
   Newton envelope, one chamber being impossible in every degree;
6. the remaining genus-one resonance has no positive first-channel degree
   and reduces to one explicit exceptional degree ridge;
7. the coupled equation is affine algebraic in the ell-channel, with a
   factored coefficient and no ell derivatives.

Finally, for the collision-containing near miss already recorded in
``HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md``, the checker constructs the uniquely
forced normal coefficients through a requested finite order.  Nonvanishing
of that table is an exact bounded computation, not an all-order proof.
"""

from __future__ import annotations

import argparse
from itertools import permutations, product

import sympy as sp


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def truncate_in_x(expression: sp.Expr, x: sp.Symbol, order: int) -> sp.Expr:
    """Return the exact polynomial truncation modulo x**(order+1)."""

    return sp.Add(
        *(
            coefficient * x ** monomial[0]
            for monomial, coefficient in sp.Poly(
                sp.expand(expression), x
            ).terms()
            if monomial[0] <= order
        )
    )


def truncated_determinant(
    matrix: sp.Matrix, x: sp.Symbol, order: int
) -> sp.Expr:
    """Expand a 4-by-4 determinant inside K[x]/(x**(order+1))."""

    determinant = sp.Integer(0)
    for permutation in permutations(range(4)):
        product = sp.Integer(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            product = truncate_in_x(
                product * matrix[row, column], x, order
            )
        determinant = truncate_in_x(determinant + product, x, order)
    return determinant


def verify_semilinear_structure() -> None:
    x, y, p, q, r = sp.symbols("x y p q r")
    L, M, N = sp.symbols("L M N")
    u = 1 + x * y
    meng_a = u**3 * p + 3 * x * u**2 * q - x**3 * r
    meng_b = (
        y**2 * u * (4 + 3 * x * y) * p
        + (y + 3 * x * y**2 * (4 + 3 * x * y)) * q
        + (2 * x - 3 * x**2 * y) * r
    )
    potential = sp.expand(L * meng_a**2 + M * meng_a + N * meng_b)
    graph_free_part = sp.expand(potential.subs(r, 0))
    graph_linear_part = sp.expand(sp.diff(potential, r).subs(r, 0))

    assert sp.expand(
        potential
        - graph_free_part
        - graph_linear_part * r
        - L * x**6 * r**2
    ) == 0
    assert sp.factor(graph_linear_part - (
        -2 * L * x**3 * (u**3 * p + 3 * x * u**2 * q)
        - M * x**3
        + N * (2 * x - 3 * x**2 * y)
    )) == 0

    # For Phi=F+cR+L*x^6*R^2, R_xx occurs only in Phi_xx.  Its
    # coefficient is the graph derivative Phi_R=c+2*L*x^6*R, so the
    # Hessian determinant is semilinear in R_xx with the displayed
    # tangential cofactor.  The abstract matrix identity checks this without
    # expanding a generic differential polynomial.
    phi_r, r_xx = sp.symbols("phi_r r_xx")
    h00 = sp.symbols("h00")
    mixed = sp.symbols("h01 h02 h03")
    tangential_symbols = sp.symbols("h11 h12 h13 h22 h23 h33")
    tangential = sp.Matrix(
        (
            (tangential_symbols[0], tangential_symbols[1], tangential_symbols[2]),
            (tangential_symbols[1], tangential_symbols[3], tangential_symbols[4]),
            (tangential_symbols[2], tangential_symbols[4], tangential_symbols[5]),
        )
    )
    hessian = sp.zeros(4)
    hessian[0, 0] = h00 + phi_r * r_xx
    hessian[0, 1:] = sp.Matrix(1, 3, mixed)
    hessian[1:, 0] = sp.Matrix(3, 1, mixed)
    hessian[1:, 1:] = tangential
    coefficient = sp.factor(sp.diff(hessian.det(), r_xx))
    assert coefficient == phi_r * tangential.det()


def verify_terminal_equation() -> None:
    """Check the universal last-coefficient bordered-Hessian factor."""

    U, d = sp.symbols("U d")
    gradient = sp.Matrix(sp.symbols("g0:3"))
    h = sp.symbols("h00 h01 h02 h11 h12 h22")
    hessian = sp.Matrix(
        ((h[0], h[1], h[2]), (h[1], h[3], h[4]), (h[2], h[4], h[5]))
    )
    square_gradient = 2 * U * gradient
    square_hessian = 2 * (U * hessian + gradient * gradient.T)
    bordered = sp.zeros(4)
    bordered[0, 0] = d * (d - 1) * U**2
    bordered[0, 1:] = d * square_gradient.T
    bordered[1:, 0] = d * square_gradient
    bordered[1:, 1:] = square_hessian

    polar = (gradient.T * hessian.adjugate() * gradient)[0]
    expected = 8 * d * U**4 * (
        (d - 1) * U * hessian.det() - (d + 1) * polar
    )
    assert sp.factor(bordered.det() - expected) == 0

    # If U_m is homogeneous of tangential degree m>=2, Euler gives
    # polar=m*U_m*det(H)/(m-1).  The terminal bracket is therefore a
    # nonzero scalar times U_m*det(H) for d=2*n+6.
    m, n = sp.symbols("m n", integer=True, positive=True)
    homogeneous_bracket = sp.factor(
        (d - 1) - (d + 1) * m / (m - 1)
    )
    assert sp.factor(
        homogeneous_bracket + (d + 2 * m - 1) / (m - 1)
    ) == 0
    assert sp.expand((2 * n + 6) + 2 * m - 1) == 2 * n + 2 * m + 5


def verify_coupled_five_channel_terminal_gate() -> None:
    """Classify the terminal face of the coupled five-function packet."""

    n = sp.symbols("n", integer=True, positive=True)
    y, p, q, P, Q = sp.symbols("y p q P Q")
    A, B, C, D, E = sp.symbols("A B C D E")
    terminal_coefficient = (
        A * y ** (n + 5)
        + B * y ** (n + 3) * p
        + C * y ** (n + 3)
        + D * y ** (n + 2) * q
        + E * y ** (n + 1) * p**2
    )
    variables = (y, p, q)
    gradient = sp.Matrix(
        [sp.diff(terminal_coefficient, variable) for variable in variables]
    )
    hessian = sp.hessian(terminal_coefficient, variables)
    d = 2 * n + 6
    terminal_bracket = sp.expand(
        (d - 1) * terminal_coefficient * hessian.det()
        - (d + 1) * (gradient.T * hessian.adjugate() * gradient)[0]
    )
    radial_bracket = sp.factor(
        terminal_bracket.subs({p: P * y, q: Q * y})
    )
    expected = D**2 * y ** (4 * n + 6) * (
        (
            8 * A * E * n**2
            + 14 * A * E * n
            - 40 * A * E
            + 2 * B**2 * n
            + 7 * B**2
        )
        * y**2
        + 2 * B * E * P * (4 * n**2 + 15 * n + 8) * y
        + 2
        * E
        * (n + 2)
        * (4 * n + 11)
        * (C + D * Q + E * P**2)
    )
    assert sp.factor(radial_bracket - expected) == 0

    # If D is nonzero, the Q coefficient first forces E=0; the remaining
    # y**2 coefficient is then D**2*B**2*(2*n+7), so B=0 in
    # characteristic zero.  Both resulting branches vanish identically.
    q_coefficient = sp.factor(sp.diff(expected, Q))
    assert q_coefficient == (
        2 * D**3 * E * y ** (4 * n + 6) * (n + 2) * (4 * n + 11)
    )
    assert sp.factor(expected.subs(E, 0)) == (
        B**2 * D**2 * y ** (4 * n + 8) * (2 * n + 7)
    )
    assert expected.subs(D, 0) == 0
    assert expected.subs({B: 0, E: 0}) == 0
    assert sp.factor(terminal_bracket.subs(D, 0)) == 0
    assert sp.factor(terminal_bracket.subs({B: 0, E: 0})) == 0


def one_channel_extremal_equation() -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Construct the normalized pure weight-five determinant face."""

    x, y, p, q, z = sp.symbols("x y p q z")
    f = sp.Function("f")
    graph = y**5 * f(x * y)
    u = 1 + x * y
    meng_a = u**3 * p + 3 * x * u**2 * q - x**3 * graph
    meng_b = (
        y**2 * u * (4 + 3 * x * y) * p
        + (y + 3 * x * y**2 * (4 + 3 * x * y)) * q
        + (2 * x - 3 * x**2 * y) * graph
    )

    # Scaling R by N/L reduces the leading face to L=N=1; M has lower
    # weight and does not enter it.  Differentiation must precede p=q=0.
    hessian_on_axis = sp.hessian(
        meng_a**2 + meng_b, (x, y, p, q)
    ).subs({p: 0, q: 0})
    determinant = hessian_on_axis.det(method="domain-ge")
    extremal = sp.cancel(determinant.subs(x, z / y).doit() / y**6)
    f0, f1, f2 = sp.symbols("f0 f1 f2")
    extremal = sp.expand(
        extremal.xreplace(
            {
                f(z): f0,
                sp.diff(f(z), z): f1,
                sp.diff(f(z), z, 2): f2,
            }
        )
    )
    assert not extremal.has(y)
    assert sp.degree(extremal, f2) == 1

    expected_f2 = -2 * z * (z + 1) ** 6 * (2 * z**5 * f0 - 3 * z + 2) * (
        6 * z**6 * f0
        + 12 * z**5 * f0
        + 6 * z**4 * f0
        - 9 * z**2
        - 12 * z
        + 1
    ) ** 2
    assert sp.factor(sp.diff(extremal, f2) - expected_f2) == 0
    return extremal, (z, f0, f1, f2)


def verify_one_channel_dominant_balance() -> None:
    extremal, (z, f0, f1, f2) = one_channel_extremal_equation()
    degree, leading_coefficient = sp.symbols("degree leading_coefficient")
    sectors: dict[tuple[int, int], sp.Expr] = {}

    for (z_power, f_power, f1_power, f2_power), coefficient in sp.Poly(
        extremal, z, f0, f1, f2
    ).terms():
        function_degree = f_power + f1_power + f2_power
        intercept = z_power - f1_power - 2 * f2_power
        contribution = (
            coefficient
            * leading_coefficient**function_degree
            * degree**f1_power
            * (degree * (degree - 1)) ** f2_power
        )
        key = (function_degree, intercept)
        sectors[key] = sp.expand(sectors.get(key, 0) + contribution)

    expected = {
        4: (22, -144 * leading_coefficient**4 * (degree**2 - degree - 21)),
        3: (18, 216 * leading_coefficient**3 * (3 * degree**2 - 3 * degree - 64)),
        2: (14, -972 * leading_coefficient**2 * (degree**2 - degree - 23)),
        1: (10, 486 * leading_coefficient * (degree - 6) * (degree + 5)),
        0: (6, sp.Integer(2916)),
    }
    for function_degree, (intercept, coefficient) in expected.items():
        available = [
            current_intercept
            for current_degree, current_intercept in sectors
            if current_degree == function_degree
            and sectors[(current_degree, current_intercept)] != 0
        ]
        assert max(available) == intercept
        assert sp.factor(
            sectors[(function_degree, intercept)] - coefficient
        ) == 0

    # For every polynomial degree m>=0, the quartic sector strictly
    # dominates the other four sectors.  Its resonance polynomial has no
    # integral root because its discriminant is 85, not a square.
    assert sp.discriminant(degree**2 - degree - 21, degree) == 85
    assert sp.integer_nthroot(85, 2) == (9, False)


def coupled_five_channel_extremal_equation() -> tuple[
    sp.Expr, tuple[sp.Symbol, ...]
]:
    """Extract the complete weight-six five-function equation efficiently."""

    x, y, p, q, z, P, Q = sp.symbols("x y p q z P Q")
    L, M, N = sp.symbols("L M N")
    functions = [sp.Function(name) for name in "fghjl"]
    f, g, h, j, ell = functions
    graph = (
        y**5 * f(x * y)
        + y**3 * p * g(x * y)
        + y**3 * h(x * y)
        + y**2 * q * j(x * y)
        + y * p**2 * ell(x * y)
    )
    u = 1 + x * y
    meng_a = u**3 * p + 3 * x * u**2 * q - x**3 * graph
    meng_b = (
        y**2 * u * (4 + 3 * x * y) * p
        + (y + 3 * x * y**2 * (4 + 3 * x * y)) * q
        + (2 * x - 3 * x**2 * y) * graph
    )
    hessian = sp.hessian(
        sp.expand(L * meng_a**2 + M * meng_a + N * meng_b), (x, y, p, q)
    )
    hessian = hessian.subs({x: z / y, p: P * y, q: Q * y}).doit()

    derivative_symbols: list[sp.Symbol] = []
    replacements: dict[sp.Expr, sp.Symbol] = {}
    for function, name in zip(functions, "fghjl", strict=True):
        for derivative_order in range(3):
            atom = (
                function(z)
                if derivative_order == 0
                else sp.diff(function(z), z, derivative_order)
            )
            symbol = sp.Symbol(f"{name}{derivative_order}")
            derivative_symbols.append(symbol)
            replacements[atom] = symbol
    hessian = hessian.applyfunc(
        lambda entry: sp.expand(entry.xreplace(replacements))
    )

    def laurent_coefficients(entry: sp.Expr) -> dict[int, sp.Expr]:
        coefficients: dict[int, sp.Expr] = {}
        for term in sp.Add.make_args(entry):
            exponent = int(term.as_powers_dict().get(y, 0))
            coefficients[exponent] = (
                coefficients.get(exponent, 0) + term / y**exponent
            )
        return {
            exponent: sp.expand(coefficient)
            for exponent, coefficient in coefficients.items()
        }

    entries = [
        [laurent_coefficients(hessian[row, column]) for column in range(4)]
        for row in range(4)
    ]
    weight_six_parts: list[sp.Expr] = []
    for permutation in permutations(range(4)):
        coefficient_sets = [
            tuple(entries[row][permutation[row]]) for row in range(4)
        ]
        for exponents in product(*coefficient_sets):
            if sum(exponents) != 6:
                continue
            weight_six_parts.append(
                permutation_sign(permutation)
                * sp.prod(
                    entries[row][permutation[row]][exponents[row]]
                    for row in range(4)
                )
            )
    extremal = sp.expand(sp.Add(*weight_six_parts))

    # The slopes and the h-channel cancel from the complete face.  The
    # remaining 2,348 terms are retained for the Newton analysis rather than
    # printed as a formula.
    assert not extremal.has(P, Q)
    assert not extremal.has(M)
    assert not any(extremal.has(symbol) for symbol in derivative_symbols[6:9])
    assert len(sp.Add.make_args(extremal)) == 2348
    assert max(
        sum(monomial[1:])
        for monomial, _ in sp.Poly(
            extremal, z, *derivative_symbols
        ).terms()
    ) == 8
    return extremal, (z, L, N, *derivative_symbols)


def verify_coupled_five_channel_newton_gate(
    extremal: sp.Expr, data: tuple[sp.Symbol, ...]
) -> None:
    """Verify the upper Newton envelope and its two chamber equations."""

    z, L, N, *derivative_symbols = data
    a, b, d, e = sp.symbols("a b d e", integer=True, nonnegative=True)
    A, B, D, E = sp.symbols("A B D E")
    degrees: list[sp.Symbol | None] = [a, b, None, d, e]
    leading_coefficients: list[sp.Symbol | None] = [A, B, None, D, E]
    sectors: dict[tuple[int, int, int, int, int], sp.Expr] = {}

    for monomial, coefficient in sp.Poly(
        extremal, z, *derivative_symbols
    ).terms():
        intercept = monomial[0]
        channel_counts: list[int] = []
        contribution = coefficient
        for channel in range(5):
            derivative_powers = monomial[
                1 + 3 * channel : 1 + 3 * channel + 3
            ]
            channel_count = sum(derivative_powers)
            channel_counts.append(channel_count)
            if channel == 2:
                assert channel_count == 0
                continue
            degree = degrees[channel]
            leading_coefficient = leading_coefficients[channel]
            assert degree is not None and leading_coefficient is not None
            contribution *= leading_coefficient**channel_count
            for derivative_order in (1, 2):
                power = derivative_powers[derivative_order]
                contribution *= sp.prod(
                    degree - shift for shift in range(derivative_order)
                ) ** power
                intercept -= derivative_order * power
        key = (
            channel_counts[0],
            channel_counts[1],
            channel_counts[3],
            channel_counts[4],
            intercept,
        )
        sectors[key] = sectors.get(key, 0) + contribution

    sectors = {
        key: sp.factor(coefficient)
        for key, coefficient in sectors.items()
        if sp.factor(coefficient) != 0
    }
    assert len(sectors) == 396

    f_ell_key = (5, 0, 2, 1, 22)
    g_square_key = (4, 2, 2, 0, 22)

    # On a+e>=2*b the f**5*j**2*ell sector is on top; on the
    # complementary chamber the f**4*g**2*j**2 sector is on top.  These
    # elementary dual-cone checks certify the claim for all 396 sectors.
    for key in sectors:
        f_delta = (
            f_ell_key[0] - key[0],
            f_ell_key[1] - key[1],
            f_ell_key[2] - key[2],
            f_ell_key[3] - key[3],
            f_ell_key[4] - key[4],
        )
        f_lower = max(sp.Rational(0), -sp.Rational(f_delta[1], 2))
        f_upper = min(sp.Rational(f_delta[0]), sp.Rational(f_delta[3]))
        assert f_delta[2] >= 0 and f_delta[4] >= 0
        assert f_lower <= f_upper
        if key != f_ell_key:
            assert f_delta[2] > 0 or f_delta[4] > 0 or f_upper > 0

        g_delta = (
            g_square_key[0] - key[0],
            g_square_key[1] - key[1],
            g_square_key[2] - key[2],
            g_square_key[3] - key[3],
            g_square_key[4] - key[4],
        )
        g_lower = max(
            sp.Rational(0),
            -sp.Rational(g_delta[0]),
            -sp.Rational(g_delta[3]),
        )
        g_upper = sp.Rational(g_delta[1], 2)
        assert g_delta[2] >= 0 and g_delta[4] >= 0
        assert g_lower <= g_upper
        if key != g_square_key:
            assert g_delta[2] > 0 or g_delta[4] > 0 or g_upper > 0

    f_resonance = (
        2 * a**3
        - 4 * a**2 * d
        + 15 * a**2
        - 2 * a * d**2
        - 34 * a * d
        + 29 * a
        - 30 * d
        + 60
    )
    g_resonance = (
        a**2
        - 2 * a * b**2
        + 4 * a * b * d
        - 8 * a * b
        - 2 * a * d**2
        + 6 * a * d
        - a
        - b**2
        + 4 * b * d
        - 18 * b
        - 4 * d**2
        + 6 * d
        - 21
    )
    assert sp.factor(
        sectors[f_ell_key] + 32 * A**5 * D**2 * E * L**4 * f_resonance
    ) == 0
    # On the wall, replacing the leading ell coefficient E by its next
    # coefficient lowers the z-degree by one.  Because this sector is
    # affine in ell, the multiplier in that next equation is the derivative
    # below; it is nonzero once the integral f-resonances are removed.
    wall_next_ell_multiplier = sp.diff(sectors[f_ell_key], E)
    assert sp.factor(
        wall_next_ell_multiplier
        + 32 * A**5 * D**2 * L**4 * f_resonance
    ) == 0
    assert sp.factor(
        sectors[g_square_key]
        + 16 * A**4 * B**2 * D**2 * L**4 * g_resonance
    ) == 0

    # In the g**2 chamber, put t=2*b-a>0.  The discriminant of the
    # quadratic g_resonance(d) is the negative polynomial below, so this
    # chamber has no polynomial balance at any nonnegative degrees.
    t = sp.symbols("t", positive=True)
    discriminant = sp.factor(sp.discriminant(g_resonance, d))
    positive_polynomial = (
        a**3
        + 6 * a**2 * t
        + 22 * a**2
        + a * t**2
        + 44 * a * t
        + 116 * a
        + 60 * t
        + 150
    )
    assert sp.factor(
        discriminant.subs(b, (a + t) / 2) + 2 * positive_polynomial
    ) == 0
    assert sp.Poly(g_resonance, d).LC() == -2 * a - 4

    # The other strict chamber reduces to one explicit genus-one resonance.
    assert sp.factor(
        sp.discriminant(f_resonance, d)
        - 4 * (a + 3) * (a + 5) * (2 * a + 1) * (4 * a + 15)
    ) == 0
    assert sp.factor(f_resonance.subs(a, 0) + 30 * (d - 2)) == 0

    # If j has degree zero, its constant leading coefficient combines all
    # three tied j-sectors into the exact square (D-3)**2.
    f_degree_zero = sum(
        sectors[(5, 0, j_power, 1, 22)] for j_power in range(3)
    ).subs(d, 0)
    g_degree_zero = sum(
        sectors[(4, 2, j_power, 0, 22)] for j_power in range(3)
    ).subs(d, 0)
    f_zero_resonance = 2 * a**3 + 15 * a**2 + 29 * a + 60
    g_zero_resonance = (
        a**2 - 2 * a * b**2 - 8 * a * b - a - b**2 - 18 * b - 21
    )
    assert sp.factor(
        f_degree_zero
        + 32
        * A**5
        * E
        * L**4
        * f_zero_resonance
        * (D - 3) ** 2
    ) == 0
    assert sp.factor(
        g_degree_zero
        + 16
        * A**4
        * B**2
        * L**4
        * g_zero_resonance
        * (D - 3) ** 2
    ) == 0

    # Close every positive-a integral point of the f-resonance without an
    # elliptic-curve black box.  Modulo a, P_F=0 gives a|30*(d-2).
    assert sp.factor(
        f_resonance.subs(d, 1)
        - (2 * a**3 + 4 * a**2 + 7 * a * (a - 1) + 30)
    ) == 0
    assert sp.rem(f_resonance, a, domain=sp.ZZ[d]) == -30 * (d - 2)
    k = sp.symbols("k", integer=True, nonnegative=True)
    reduced_quadratic = (
        a**2 * (k**2 + 60 * k - 900)
        + a * (630 * k - 3150)
        + 450 * k
        + 21150
    )
    assert sp.factor(
        f_resonance.subs(d, 2 + k * a / 30)
        + a * reduced_quadratic / 450
    ) == 0
    assert sp.expand(
        (k**2 + 60 * k - 900)
        - (49 + (k - 13) * (k + 73))
    ) == 0
    assert 630 * k - 3150 == 630 * (k - 5)
    assert 450 * k + 21150 == 450 * (k + 47)

    reduced_discriminant = sp.factor(sp.discriminant(reduced_quadratic, a))
    assert sp.factor(
        reduced_discriminant
        - 900 * (51 - k) * (25 - k) * (75 - 2 * k)
    ) == 0
    for k_value in range(13):
        square_candidate = int(
            ((51 - k_value) * (25 - k_value) * (75 - 2 * k_value))
        )
        assert sp.integer_nthroot(square_candidate, 2)[1] is False

    # The only a=0 resonance is d=2.  Its next f**5*j*ell sector is
    # nonzero.  Comparing it with the g**2*j**2 sector leaves only the
    # degree ridge e=2*b+2, with the displayed amplitude equation.
    f_next_key = (5, 0, 1, 1, 22)
    f_next_resonance = (
        2 * a**3
        - 2 * a**2 * d
        + 15 * a**2
        - 17 * a * d
        + 29 * a
        - 15 * d
        + 60
    )
    assert sp.factor(
        sectors[f_next_key]
        - 192 * A**5 * D * E * L**4 * f_next_resonance
    ) == 0
    assert f_next_resonance.subs({a: 0, d: 2}) == 30
    exceptional_g_resonance = sp.factor(
        g_resonance.subs({a: 0, d: 2})
    )
    assert sp.factor(
        exceptional_g_resonance + (b + 5) ** 2
    ) == 0
    exceptional_tie = sp.factor(
        sectors[f_next_key].subs({a: 0, d: 2})
        + sectors[g_square_key].subs({a: 0, d: 2})
    )
    expected_tie = (
        16
        * L**4
        * A**4
        * D
        * (360 * A * E + B**2 * D * (b + 5) ** 2)
    )
    assert sp.factor(exceptional_tie - expected_tie) == 0


def verify_constant_j_exceptional_face(
    extremal: sp.Expr, data: tuple[sp.Symbol, ...]
) -> None:
    """Close both strict Newton chambers on the exceptional j=3 slice."""

    z, L, _N, *derivative_symbols = data
    exceptional = sp.expand(
        extremal.subs(
            {
                derivative_symbols[9]: 3,
                derivative_symbols[10]: 0,
                derivative_symbols[11]: 0,
            }
        )
    )
    assert len(sp.Add.make_args(exceptional)) == 512

    a, b, e = sp.symbols("a b e", integer=True, nonnegative=True)
    A, B, E = sp.symbols("A B E")
    channel_starts = (0, 3, 12)
    degrees = (a, b, e)
    leading_coefficients = (A, B, E)
    sectors: dict[tuple[int, int, int, int], sp.Expr] = {}
    for monomial, coefficient in sp.Poly(
        exceptional, z, *derivative_symbols
    ).terms():
        intercept = monomial[0]
        counts: list[int] = []
        contribution = coefficient
        for start, degree, leading_coefficient in zip(
            channel_starts, degrees, leading_coefficients, strict=True
        ):
            derivative_powers = monomial[1 + start : 1 + start + 3]
            count = sum(derivative_powers)
            counts.append(count)
            contribution *= leading_coefficient**count
            for derivative_order in (1, 2):
                power = derivative_powers[derivative_order]
                contribution *= sp.prod(
                    degree - shift for shift in range(derivative_order)
                ) ** power
                intercept -= derivative_order * power
        assert sum(monomial[7:13]) == 0
        key = (counts[0], counts[1], counts[2], intercept)
        sectors[key] = sectors.get(key, 0) + contribution
    sectors = {
        key: sp.factor(coefficient)
        for key, coefficient in sectors.items()
        if sp.factor(coefficient) != 0
    }
    assert len(sectors) == 121

    f_ell_key = (5, 0, 1, 20)
    g_square_key = (4, 2, 0, 20)
    for key in sectors:
        f_delta = (
            f_ell_key[0] - key[0],
            f_ell_key[1] - key[1],
            f_ell_key[2] - key[2],
            f_ell_key[3] - key[3],
        )
        f_lower = max(sp.Rational(0), -sp.Rational(f_delta[1], 2))
        f_upper = min(sp.Rational(f_delta[0]), sp.Rational(f_delta[2]))
        assert f_delta[3] >= 0 and f_lower <= f_upper
        if key != f_ell_key:
            assert f_delta[3] > 0 or f_upper > 0

        g_delta = (
            g_square_key[0] - key[0],
            g_square_key[1] - key[1],
            g_square_key[2] - key[2],
            g_square_key[3] - key[3],
        )
        g_lower = max(
            sp.Rational(0),
            -sp.Rational(g_delta[0]),
            -sp.Rational(g_delta[2]),
        )
        g_upper = sp.Rational(g_delta[1], 2)
        assert g_delta[3] >= 0 and g_lower <= g_upper
        if key != g_square_key:
            assert g_delta[3] > 0 or g_upper > 0

    f_exceptional = 2 * a**3 + 19 * a**2 + 61 * a + 90
    g_exceptional = (
        a**2
        - 2 * a * b**2
        - 12 * a * b
        - 9 * a
        - b**2
        - 22 * b
        - 31
    )
    assert sp.factor(
        sectors[f_ell_key] + 1152 * L**4 * A**5 * E * f_exceptional
    ) == 0
    assert sp.factor(
        sectors[g_square_key]
        + 576 * L**4 * A**4 * B**2 * g_exceptional
    ) == 0

    t = sp.symbols("t", positive=True)
    positive_polynomial = (
        2 * a**3
        + 4 * a**2 * t
        + 21 * a**2
        + 2 * a * t**2
        + 26 * a * t
        + 80 * a
        + t**2
        + 44 * t
        + 124
    )
    assert sp.factor(
        4 * g_exceptional.subs(b, (a + t) / 2) + positive_polynomial
    ) == 0


def verify_coupled_ell_elimination(
    extremal: sp.Expr, data: tuple[sp.Symbol, ...]
) -> None:
    """Verify exact algebraic elimination of the fifth channel."""

    z, L, N, *derivative_symbols = data
    ell0, ell1, ell2 = derivative_symbols[12:15]
    assert sp.diff(extremal, ell1) == 0
    assert sp.diff(extremal, ell2) == 0
    assert sp.Poly(extremal, ell0).degree() == 1

    coefficient = sp.expand(sp.diff(extremal, ell0))
    first_channel = derivative_symbols[0]
    graph_factor = 2 * L * first_channel * z**5 - 3 * N * z + 2 * N
    residual_factor = sp.cancel(coefficient / (2 * z * graph_factor))
    assert sp.denom(residual_factor) == 1
    assert len(sp.Add.make_args(sp.expand(residual_factor))) == 311
    assert not any(
        residual_factor.has(symbol)
        for symbol in (
            *derivative_symbols[3:9],
            *derivative_symbols[12:15],
        )
    )
    assert sp.factor(
        coefficient - 2 * z * graph_factor * residual_factor
    ) == 0

    remainder = sp.expand(extremal.subs(ell0, 0))
    assert len(sp.Add.make_args(remainder)) == 1895
    f0, f1, _f2 = derivative_symbols[0:3]
    g0 = derivative_symbols[3]
    j0 = derivative_symbols[9]
    expected_axis = N**3 * (
        80 * L * f0 * j0
        + 492 * L * f0
        - 8 * L * f1
        + 4 * N * g0**2
        - 64 * N * g0 * j0
        - 356 * N * g0
        + 256 * N * j0**2
        + 2848 * N * j0
        + 7921 * N
    )
    assert sp.factor(remainder.subs(z, 0) - expected_axis) == 0


def near_miss_prefix(maximum_order: int) -> list[dict[str, object]]:
    """Compute the unique exact normal prefix of the collision near miss."""

    x, y, p, q = sp.symbols("x y p q")
    u = 1 + x * y
    a0 = u**3 * p + 3 * x * u**2 * q
    b0 = (
        y**2 * u * (4 + 3 * x * y) * p
        + (y + 3 * x * y**2 * (4 + 3 * x * y)) * q
    )
    free_part = sp.expand(a0**2 + 13 * a0 + 2 * b0)
    linear_part = sp.expand(
        -2 * x**3 * a0 - 13 * x**3 + 4 * x - 6 * x**2 * y
    )
    target = sp.Rational(17165601, 25)
    graph = (
        sp.Rational(51, 100) * y**5
        - sp.Rational(47, 10) * y**3 * p
        - sp.Rational(123, 20) * y**2 * q
    )
    expected_extremal = {
        1: ((4, 0, 0), sp.Rational(1989, 40)),
        2: ((7, 0, 0), sp.Rational(459, 500)),
        3: ((8, 0, 0), sp.Rational(235467, 4000)),
        4: ((9, 0, 0), -sp.Rational(51583797, 20000)),
        5: ((10, 0, 0), sp.Rational(1094341514607, 10000000)),
    }
    expected_coupled_packet = {
        1: (
            -sp.Rational(2727, 200),
            sp.Rational(348, 5),
        ),
        2: (
            sp.Rational(243, 100),
            -sp.Rational(158457, 50),
        ),
        3: (
            -sp.Rational(10145133, 10000),
            sp.Rational(53188533, 400),
        ),
        4: (
            sp.Rational(2083052601, 50000),
            -sp.Rational(22284332193, 4000),
        ),
        5: (
            -sp.Rational(880031929203, 500000),
            sp.Rational(46567165753557, 200000),
        ),
    }
    rows: list[dict[str, object]] = []

    for normal_order in range(1, maximum_order + 1):
        determinant_order = normal_order - 1
        potential = truncate_in_x(
            free_part + linear_part * graph + x**6 * graph**2,
            x,
            determinant_order + 2,
        )
        hessian = sp.hessian(potential, (x, y, p, q)).applyfunc(
            lambda entry: truncate_in_x(entry, x, determinant_order)
        )
        determinant = truncated_determinant(hessian, x, determinant_order)
        residual = sp.Poly(determinant - target, x).coeff_monomial(
            x**determinant_order
        )
        coefficient = sp.cancel(
            residual / (32 * normal_order * (normal_order + 1))
        )
        coefficient_polynomial = sp.Poly(sp.expand(coefficient), y, p, q)
        total_degree = coefficient_polynomial.total_degree()
        top_terms = [
            (monomial, value)
            for monomial, value in coefficient_polynomial.terms()
            if sum(monomial) == total_degree
        ]
        assert coefficient != 0
        assert expected_extremal[normal_order] in top_terms

        # D_k and E_k are the z**k coefficients of the q- and p**2-
        # channels in the coupled packet.  They provide a compact regression
        # for the emerging five-function face.
        coupled_d = coefficient_polynomial.coeff_monomial(
            y ** (normal_order + 2) * q
        )
        coupled_e = coefficient_polynomial.coeff_monomial(
            y ** (normal_order + 1) * p**2
        )
        assert (coupled_d, coupled_e) == expected_coupled_packet[normal_order]

        # Apply the universal terminal equation to the complete coefficient,
        # including every term below the five displayed channels.  A nonzero
        # bracket proves that this exact prefix cannot terminate at k.
        variables = (y, p, q)
        gradient = sp.Matrix(
            [sp.diff(coefficient, variable) for variable in variables]
        )
        hessian = sp.hessian(coefficient, variables)
        d = 2 * normal_order + 6
        terminal_bracket = sp.expand(
            (d - 1) * coefficient * hessian.det()
            - (d + 1) * (gradient.T * hessian.adjugate() * gradient)[0]
        )
        assert terminal_bracket != 0
        rows.append(
            {
                "normal_order": normal_order,
                "coefficient_count": len(coefficient_polynomial.terms()),
                "total_degree": total_degree,
                "extremal_monomial": str(expected_extremal[normal_order][0]),
                "extremal_coefficient": str(expected_extremal[normal_order][1]),
                "coupled_D": str(coupled_d),
                "coupled_E": str(coupled_e),
                "terminal_gate": "nonzero",
            }
        )
        graph = truncate_in_x(
            graph + x**normal_order * coefficient, x, normal_order
        )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--maximum-normal-order",
        type=int,
        default=5,
        choices=range(1, 6),
        metavar="{1,2,3,4,5}",
        help="last exact normal coefficient computed for the near miss",
    )
    arguments = parser.parse_args()

    verify_semilinear_structure()
    verify_terminal_equation()
    verify_one_channel_dominant_balance()
    verify_coupled_five_channel_terminal_gate()
    coupled_extremal, coupled_data = coupled_five_channel_extremal_equation()
    verify_coupled_five_channel_newton_gate(coupled_extremal, coupled_data)
    verify_constant_j_exceptional_face(coupled_extremal, coupled_data)
    verify_coupled_ell_elimination(coupled_extremal, coupled_data)
    rows = near_miss_prefix(arguments.maximum_normal_order)

    print("PASS HC4MYPT1: every terminal coefficient satisfies the cone gate")
    print("PASS HC4MYPT2: the pure weight-five one-channel ansatz is empty")
    print("PASS HC4MYPT3: the coupled terminal packet has exactly two branches")
    print("PASS HC4MYPT4: the coupled Newton envelope has two exact chambers")
    print("PASS HC4MYPT5: the genus-one chamber reduces to one degree ridge")
    print("PASS HC4MYPT6: the ell-channel is algebraically unit-triangular")
    for row in rows:
        print(
            "NEAR_MISS_PREFIX "
            + " ".join(f"{key}={value}" for key, value in row.items())
        )
    print(
        "BOUNDED RESULT: the formal branch above the collision-containing "
        "plane-flat near miss is nonzero and fails the terminal equation "
        "at every normal order from 1 through "
        f"{arguments.maximum_normal_order}"
    )
    print("COLLISION SCOPE: higher normal corrections need not retain the marks")
    print("SCOPE: the bounded prefix is not an all-order nontermination proof")


if __name__ == "__main__":
    main()
