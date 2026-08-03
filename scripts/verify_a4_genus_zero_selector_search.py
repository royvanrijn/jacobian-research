#!/usr/bin/env python3
"""Verify the first bounded genus-zero selector search for the A4 core.

The bound is:

* ordinary total degree at most three in (a,b,T);
* selected-root degree at most two; and
* horizontal norm degree strictly below sixteen.

The script computes the complete valuation-filtered ansatz through degree
three.  It then proves that every exact selector in this ansatz has norm
degree at least sixteen.  The bound is sharp: explicit root-linear degree
three selectors have the required exceptional orders and irreducible
degree-sixteen norms.  Their genera are checked separately by
``verify_a4_genus_zero_selector_search.sing``.
"""

from fractions import Fraction

import sympy as sp


a, b, c, T, X, z = sp.symbols("a b c T X z")
x, y, u, k, t = sp.symbols("x y u k t")

rho = b**2 + 3 * b + 9
B = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
A = a**3 - b**3 - 9 * b**2 - 27 * b - 54
C = a**3 - b**3 + 27
P = (
    T**4
    - 6 * A * B * T**2
    - 8 * B**3 * T
    + B**2 * (9 * A**2 - 12 * C * B)
)


def monomials(total_degree, root_degree):
    """Return the ordinary-degree ansatz monomials in a stable order."""

    return tuple(
        a**a_power * b**b_power * T**root_power
        for root_power in range(root_degree + 1)
        for a_power in range(total_degree + 1)
        for b_power in range(total_degree + 1)
        if a_power + b_power + root_power <= total_degree
    )


def truncated_coefficients(expression, parameter, threshold, extras):
    """Return coefficient labels below a divisorial threshold.

    The cluster residue field is Q(c)/(c^2+27).  The returned labels retain
    the parameter order, the c exponent, and all generic divisor variables.
    """

    truncated = sp.series(
        expression,
        parameter,
        0,
        threshold,
    ).removeO()
    output = {}
    for order in range(threshold):
        coefficient = sp.expand(truncated).coeff(parameter, order)
        if coefficient == 0:
            continue
        numerator, denominator = sp.cancel(coefficient).as_numer_denom()
        assert not denominator.has(c, *extras)
        residue = sp.rem(
            sp.Poly(sp.expand(numerator), c),
            sp.Poly(c**2 + 27, c),
        ).as_expr()
        for exponent, value in sp.Poly(
            sp.expand(residue),
            c,
            *extras,
        ).terms():
            if value:
                output[(order, exponent)] = sp.cancel(value / denominator)
    return output


def valuation_truncations(expression):
    """Map an ansatz element to all seven truncated normalized branches."""

    result = {}

    # Write c=2b+3 and z=rho.  If z=s^m q, then
    #
    #   c=c_0-(2c_0/27)s^m q+O(s^(2m)),  c_0^2=-27.
    #
    # Only the displayed term can contribute below the required thresholds.

    # S1: a=x, z=xy, T=x(27-3c)y/2+O(x^2), require order >=2.
    c_series = c - sp.Rational(2, 27) * c * x * y
    substitutions = {
        a: x,
        b: (c_series - 3) / 2,
        T: x * (27 - 3 * c_series) * y / 2,
    }
    result["S1"] = truncated_coefficients(
        expression.subs(substitutions),
        x,
        2,
        (y,),
    )

    # S2: a=x, z=x^2y, T=x^2(27-3c)y/2+O(x^3), require >=3.
    c_series = c - sp.Rational(2, 27) * c * x**2 * y
    substitutions = {
        a: x,
        b: (c_series - 3) / 2,
        T: x**2 * (27 - 3 * c_series) * y / 2,
    }
    result["S2"] = truncated_coefficients(
        expression.subs(substitutions),
        x,
        3,
        (y,),
    )

    # Fs: a=u^2k, z=u^3k.  For G=T+a^3=u^3r, Hensel lifting of
    # the simple exceptional root gives
    #
    #   r=(27-3c)k/2 + 9(c+3)k^2 u^2/8 + O(u^3).
    #
    # This is exactly the jet needed to test order >=6.
    c_series = c - sp.Rational(2, 27) * c * u**3 * k
    simple_r = (
        (27 - 3 * c_series) * k / 2
        + sp.Rational(9, 8) * (c + 3) * k**2 * u**2
    )
    substitutions = {
        a: u**2 * k,
        b: (c_series - 3) / 2,
        T: u**3 * simple_r - u**6 * k**3,
    }
    result["Fs"] = truncated_coefficients(
        expression.subs(substitutions),
        u,
        6,
        (k,),
    )

    # The other four branches need no root jet below their thresholds.
    branch_data = {
        "Q": (x, 3, x, x**3 * y, x**3 * t, (y, t)),
        "R2": (x, 6, x**3, x**6 * y, x**6 * t, (y, t)),
        "Ft": (u, 3, u**2 * k, u**3 * k, u**3 * t, (k, t)),
        "R1": (x, 3, x**3, x**3 * y, x**3 * t, (y, t)),
    }
    for name, (
        parameter,
        threshold,
        a_value,
        z_value,
        t_value,
        extras,
    ) in branch_data.items():
        c_series = c - sp.Rational(2, 27) * c * z_value
        substitutions = {
            a: a_value,
            b: (c_series - 3) / 2,
            T: t_value,
        }
        result[name] = truncated_coefficients(
            expression.subs(substitutions),
            parameter,
            threshold,
            extras,
        )
    return result


def valuation_matrix(ansatz_monomials):
    """Build the exact linear map to all forbidden valuation jets."""

    columns = []
    labels = set()
    for monomial in ansatz_monomials:
        column = {}
        for branch, truncation in valuation_truncations(monomial).items():
            for label, coefficient in truncation.items():
                full_label = (branch, label)
                column[full_label] = coefficient
                labels.add(full_label)
        columns.append(column)
    labels = tuple(sorted(labels, key=str))
    return sp.Matrix(
        [
            [column.get(label, 0) for column in columns]
            for label in labels
        ]
    )


def coefficient_vector(expression, ansatz_monomials):
    """Write an expression in the declared monomial basis."""

    polynomial = sp.Poly(sp.expand(expression), a, b, T)
    return sp.Matrix(
        [
            polynomial.coeff_monomial(monomial)
            for monomial in ansatz_monomials
        ]
    )


def assert_kernel_basis(matrix, ansatz_monomials, expected_basis):
    """Check that the displayed expressions are exactly the matrix kernel."""

    expected = sp.Matrix.hstack(
        *(
            coefficient_vector(expression, ansatz_monomials)
            for expression in expected_basis
        )
    )
    assert matrix * expected == sp.zeros(matrix.rows, expected.cols)
    assert expected.rank() == len(expected_basis)
    assert matrix.rank() + expected.rank() == len(ansatz_monomials)


# ---------------------------------------------------------------------------
# 1. Complete degree-two and degree-three valuation spaces
# ---------------------------------------------------------------------------

q0 = a**3
q1 = (
    4 * b * T
    + 81 * a * b**2
    + 243 * a * b
    + 729 * a
    - 72 * b**3
    - 324 * b**2
    - 972 * b
    - 972
) / 4
q2 = (
    4 * b**2 * T
    + 36 * T
    - 243 * a * b**2
    - 729 * a * b
    - 2187 * a
    + 216 * b**3
    + 972 * b**2
    + 2916 * b
    + 2916
) / 4
q3 = (
    3 * a * T
    + 4 * T
    - 54 * a * b**2
    - 162 * a * b
    - 486 * a
    + 12 * b**3
    - 324
) / 3
q4 = (
    a * b * T
    - 8 * T
    + 27 * a * b**2
    + 81 * a * b
    + 243 * a
    - 24 * b**3
    + 648
)
q5 = a**2 * T

degree_two_linear = monomials(2, 1)
degree_two_quadratic = monomials(2, 2)
degree_three_linear = monomials(3, 1)
degree_three_quadratic = monomials(3, 2)

matrix_2_1 = valuation_matrix(degree_two_linear)
matrix_2_2 = valuation_matrix(degree_two_quadratic)
matrix_3_1 = valuation_matrix(degree_three_linear)
matrix_3_2 = valuation_matrix(degree_three_quadratic)

assert matrix_2_1.nullspace() == []
assert_kernel_basis(matrix_2_2, degree_two_quadratic, (T**2,))
assert_kernel_basis(
    matrix_3_1,
    degree_three_linear,
    (q0, q1, q2, q3, q4, q5),
)
assert_kernel_basis(
    matrix_3_2,
    degree_three_quadratic,
    (q0, q1, q2, q3, q4, q5, T**2, a * T**2, b * T**2),
)

print("PASS: no root-linear selector of total degree at most two reaches the valuation ideal")
print("PASS: the degree-two root-quadratic valuation space is exactly <T^2>")
print("PASS: the degree-three root-linear valuation space has the displayed dimension six")
print("PASS: allowing root degree two adds exactly <T^2,aT^2,bT^2>")


# ---------------------------------------------------------------------------
# 2. Exact branch vectors and the sharp norm-degree floor
# ---------------------------------------------------------------------------

local_relation = b**2 + 3 * b + 9 - z


def weighted_order(expression, a_weight, z_weight):
    """Return the exact coefficient-plane monomial order at the cluster."""

    polynomial = sp.Poly(
        sp.rem(sp.expand(expression), local_relation, b),
        a,
        z,
        b,
    )
    return min(
        a_weight * a_power + z_weight * z_power
        for (a_power, z_power, _), coefficient in polynomial.terms()
        if coefficient
    )


def newton_root_orders(coefficient_orders):
    """Read all root orders from a lower Newton polygon."""

    points = [
        (exponent, Fraction(order))
        for exponent, order in enumerate(coefficient_orders)
        if order is not None
    ]
    hull = []
    for point in points:
        while len(hull) >= 2:
            old_slope = Fraction(
                hull[-1][1] - hull[-2][1],
                hull[-1][0] - hull[-2][0],
            )
            new_slope = Fraction(
                point[1] - hull[-1][1],
                point[0] - hull[-1][0],
            )
            if new_slope <= old_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    orders = []
    for left, right in zip(hull, hull[1:]):
        slope = Fraction(
            right[1] - left[1],
            right[0] - left[0],
        )
        orders.extend([-slope] * (right[0] - left[0]))
    return tuple(orders)


def characteristic_orders(function, weights):
    """Return the four root orders of one root-algebra function."""

    characteristic = sp.Poly(
        sp.resultant(P, X - function, T),
        X,
    )
    coefficient_orders = [
        weighted_order(characteristic.nth(exponent), *weights)
        if characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    return newton_root_orders(coefficient_orders)


resolution_rays = {
    "E1": (1, 1),
    "E2": (1, 2),
    "E3": (1, 3),
    "F": (2, 3),
}
required_orders = {
    "E1": (Fraction(2), Fraction(1), Fraction(1), Fraction(1)),
    "E2": (Fraction(3), Fraction(2), Fraction(2), Fraction(2)),
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(6), Fraction(3), Fraction(3), Fraction(3)),
}
for selector in (q1, q2, q3, q4):
    for name, weights in resolution_rays.items():
        assert characteristic_orders(selector, weights) == required_orders[name]

# Under the pole-degree weights deg(a)=deg(b)=1 and deg(T)=3, the leading
# A4 polynomial is monic of weight twelve and has nonzero constant term.
# Therefore:
#
# * a nonzero T^2 coefficient of degree e has norm degree 24+4e;
# * without T^2, a T coefficient of degree two has norm degree 20;
# * a T coefficient of degree one has norm degree 16.
#
# In the computed kernel, the respective leading coefficients are
# <1,a,b>, <a^2,ab,b^2>, and <a,b>.  If all degree-one T coefficients
# vanish after removing the higher two cases, only q0=a^3 remains, and q0
# does not have the exact vector.  Thus sixteen is the sharp norm floor.
A_lead = a**3 - b**3
B_lead = a**3 - 3 * a * b**2 + 2 * b**3
C_lead = a**3 - b**3
constant_lead = sp.expand(
    B_lead**2 * (9 * A_lead**2 - 12 * C_lead * B_lead)
)
assert constant_lead != 0

for selector, expected_degree in ((q1, 16), (q2, 20), (q3, 16), (q4, 20)):
    norm = sp.resultant(P, sp.together(selector), T)
    assert sp.Poly(norm, a, b).total_degree() == expected_degree

for selector in (q1, q3):
    norm = sp.resultant(P, sp.together(selector), T)
    _, factors = sp.factor_list(norm, a, b)
    assert len(factors) == 1
    assert factors[0][1] == 1

assert characteristic_orders(q0, resolution_rays["E1"]) == (Fraction(3),) * 4
assert characteristic_orders(T**2, resolution_rays["E1"]) == (Fraction(2),) * 4

print("PASS: q1,q2,q3,q4 have exactly the required four exceptional vectors")
print("PASS: every exact selector in the bounded ansatz has norm degree at least sixteen")
print("PASS: q1 and q3 attain degree sixteen with irreducible norms")
print("THEOREM: below norm degree sixteen no total-degree-three root-quadratic selector can have the required valuations and genus zero")
print("NOTE: degree sixteen is sharp; its full two-parameter genus stratification remains open")


# ---------------------------------------------------------------------------
# 3. Two reducible directions and the one-jet affine-modification handoff
# ---------------------------------------------------------------------------

# The two additional root-linear directions q2 and q4 reveal two exact
# cancellations.  The first is only the old root and rho boundary.  The
# second is more informative: after adjoining the quotient by a, it gives a
# degree-two root-linear function which passes every normalized valuation
# condition except one coefficient on the simple F branch.
near_selector = sp.expand((b + 6) * T - 81 * rho)
assert sp.expand(3 * q1 + q2 - T * rho) == 0
assert sp.expand(6 * q3 + q4 - a * near_selector) == 0

near_truncations = valuation_truncations(near_selector)
assert near_truncations == {
    "S1": {},
    "S2": {},
    "Fs": {(5, (1, 2)): sp.Rational(27, 4)},
    "Q": {},
    "R2": {},
    "Ft": {},
    "R1": {},
}

near_expected_orders = {
    "E1": (Fraction(2), Fraction(1), Fraction(1), Fraction(1)),
    "E2": (Fraction(3), Fraction(2), Fraction(2), Fraction(2)),
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(5), Fraction(3), Fraction(3), Fraction(3)),
}
for name, weights in resolution_rays.items():
    assert characteristic_orders(near_selector, weights) == (
        near_expected_orders[name]
    )

near_norm = sp.resultant(P, near_selector, T)
near_unit, near_factors = sp.factor_list(near_norm, a, b)
assert near_unit != 0
assert len(near_factors) == 1
assert near_factors[0][1] == 1
assert sp.Poly(near_factors[0][0], a, b).total_degree() == 16

print("PASS: 3*q1+q2 is exactly T*rho")
print("PASS: (6*q3+q4)/a is polynomial after the a-modification")
print("NEAR MISS: that quotient has only the simple-F order-five coefficient 27*c*k^2/4")
print("PASS: its coefficient-plane norm is irreducible of degree 16")
