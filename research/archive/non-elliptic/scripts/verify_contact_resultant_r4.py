#!/usr/bin/env python3
"""Exact uniform certificate for the cancellation contact resultant at r=4.

The proof has three parts.  Small m are disposed of by exact endpoint gcds.
For m>=22 a Schur--Cohn inertia computation shows that the degree-eleven
endpoint eliminant has exactly one conjugate pair in the K-root disk.  A
rational Rouche cover localizes that pair, and a Bernstein certificate shows
that its linear-subresultant z-branch lies in the opposite half-plane from
y**m.

Floating point roots are used only to propose rational comparison
polynomials in the Rouche cover.  Every inequality that certifies the root
count is subsequently checked over QQ.
"""

from __future__ import annotations

import itertools
import math

import sympy as sp
from sympy.polys.domains import QQ, ZZ
from sympy.polys.rings import ring


m, n, y, z, u, x, t = sp.symbols("m n y z u x t")
D = 1 - y


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


def endpoint_moment(k: int) -> sp.Expr:
    return (beta(k) - y * z**k * endpoint_tail(k)) / D ** (k + 1)


# At a common K_{m,4},L_{m,4} root, M_4=0 and z is nonzero.  The two
# endpoint equations are therefore the following cubic and quartic in z.
M1, M2, M3, M4 = (endpoint_moment(k) for k in range(1, 5))
E4 = sp.cancel((z**3 - 4 * z**2 * M1 + 6 * z * M2 - 4 * M3) * D**4)
F4 = sp.cancel(M4 * D**5)
assert sp.degree(E4, z) == 3
assert sp.degree(F4, z) == 4

# The last two subresultants are A_m(y)z+B_m(y) and the z-resultant.
subresultants4 = sp.subresultants(E4, F4, z)
linear_subresultant4 = sp.Poly(subresultants4[-2], z)
endpoint_resultant4 = sp.factor(subresultants4[-1])
assert linear_subresultant4.degree() == 1

eliminant4 = sp.cancel(endpoint_resultant4 / (y - 1) ** 5)
eliminant4_numerator, eliminant4_denominator = sp.together(
    eliminant4
).as_numer_denom()
assert sp.degree(eliminant4_denominator, y) == 0
H4 = sp.Poly(eliminant4_numerator, y, domain=ZZ[m]).primitive()[1]
assert H4.degree() == 11
assert H4.eval(1) != 0
assert sp.cancel(endpoint_resultant4 - (y - 1) ** 5 * eliminant4) == 0


# For 1<=m<=21, imposing z=y^m directly leaves only the clearing-denominator
# factor (y-1)^5 in the gcd.  The endpoint y=1 is w=0 and K_{m,4}(0)=1/5.
excluded_endpoint = sp.Poly((y - 1) ** 5, y).monic()
for m_value in range(1, 22):
    substitutions = {m: m_value, z: y**m_value}
    e_value = sp.Poly(
        sp.together(E4.subs(substitutions)).as_numer_denom()[0], y
    )
    f_value = sp.Poly(
        sp.together(F4.subs(substitutions)).as_numer_denom()[0], y
    )
    assert sp.gcd(e_value, f_value).monic() == excluded_endpoint


# For m>=22, invert the comparison disk |y|<rho in the usual Schur--Cohn
# fashion.  The Hermitian form has inertia (9,2): nine H-roots are outside
# the disk and exactly two are inside.  Fraction-free Bareiss elimination
# obtains all eleven leading minors in one pass.
rho = m / (m + 1)
Q4_numerator = sp.together(
    u**11 * H4.as_expr().subs(y, rho / u)
).as_numer_denom()[0]
Q4 = sp.Poly(Q4_numerator, u, domain=ZZ[m]).primitive()[1]
assert Q4.degree() == 11

coefficients = Q4.all_coeffs()
size = Q4.degree()
forward = sp.zeros(size)
reverse = sp.zeros(size)
for row in range(size):
    for column in range(row + 1):
        forward[row, column] = coefficients[row - column]
        reverse[row, column] = coefficients[size - (row - column)]
schur_cohn = forward * forward.T - reverse * reverse.T

polynomial_ring, ring_m = ring("m", ZZ)
bareiss = [
    [polynomial_ring.from_expr(schur_cohn[i, j]) for j in range(size)]
    for i in range(size)
]
previous_pivot = polynomial_ring.one
leading_minors = []
for pivot_index in range(size):
    pivot = bareiss[pivot_index][pivot_index]
    leading_minors.append(pivot)
    if pivot_index == size - 1:
        break
    for row in range(pivot_index + 1, size):
        for column in range(pivot_index + 1, size):
            numerator = (
                pivot * bareiss[row][column]
                - bareiss[row][pivot_index] * bareiss[pivot_index][column]
            )
            bareiss[row][column] = numerator.exquo(previous_pivot)
    previous_pivot = pivot

expected_minor_degrees = [65, 128, 189, 248, 305, 360, 413, 464, 513, 560, 605]
expected_minor_signs = [1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1]
assert [minor.degree() for minor in leading_minors] == expected_minor_degrees
for minor, expected_sign in zip(leading_minors, expected_minor_signs):
    shifted = minor.compose({ring_m: ring_m + 22})
    assert all(expected_sign * coefficient > 0 for coefficient in shifted.values())

# Signs of the LDL pivots are Delta_k/Delta_{k-1}.  There are nine positive
# and two negative pivots, hence nine Q-roots inside and two outside |u|=1.
ldl_signs = [expected_minor_signs[0]] + [
    expected_minor_signs[index] * expected_minor_signs[index - 1]
    for index in range(1, size)
]
assert ldl_signs.count(1) == 9
assert ldl_signs.count(-1) == 2


# Scale the exceptional pair by x=m(y-1), t=1/m.  The factor t^14 removes
# the valuation at t=0 and produces a degree-eleven polynomial P(t,x).
P4_expression = sp.cancel(
    t**14 * H4.as_expr().subs({m: 1 / t, y: 1 + t * x})
)
P4 = sp.Poly(P4_expression, x)
assert P4.degree() == 11
assert sp.Poly(P4_expression, x, t).degree(t) == 24
leading_coefficient = P4.LC()
assert all(
    coefficient > 0
    for coefficient in sp.Poly(leading_coefficient, t).all_coeffs()
)

# The localization disk in the x-plane.
center_real = sp.Rational(-523, 200)
center_imag = sp.Rational(41, 8)
radius = sp.Rational(1, 2)
x_modulus_bound = sp.Rational(313, 50)
assert x_modulus_bound**2 > (
    sp.sqrt(center_real**2 + center_imag**2) + radius
) ** 2


def polynomial_absolute_bound(polynomial: sp.Poly, endpoint: sp.Rational) -> sp.Rational:
    """Coefficientwise upper bound for |p(s)| on 0<=s<=endpoint."""

    return sum(
        abs(coefficient) * endpoint ** monomial[0]
        for monomial, coefficient in polynomial.terms()
    )


# A global derivative bound for the coefficients of the monic P(t,x).
t_max = sp.Rational(1, 22)
leading_at_zero = leading_coefficient.subs(t, 0)
derivative_bound = sp.Rational(0)
for power in range(11):
    numerator = P4.coeff_monomial(x**power)
    derivative_numerator = sp.Poly(
        sp.expand(
            sp.diff(numerator, t) * leading_coefficient
            - numerator * sp.diff(leading_coefficient, t)
        ),
        t,
    )
    coefficient_bound = (
        polynomial_absolute_bound(derivative_numerator, t_max)
        / leading_at_zero**2
    )
    derivative_bound += coefficient_bound * x_modulus_bound**power


def rational_sqrt_bound(
    value: sp.Rational, *, upper: bool, digits: int = 14
) -> sp.Rational:
    """A rigorous decimal lower/upper bound for sqrt(value)."""

    numerator, denominator = map(int, sp.fraction(value))
    scale = 10**digits
    approximation = math.isqrt(numerator * scale * scale // denominator)
    while (approximation + 1) ** 2 * denominator <= numerator * scale * scale:
        approximation += 1
    while approximation**2 * denominator > numerator * scale * scale:
        approximation -= 1
    if upper:
        approximation += 1
    return sp.Rational(approximation, scale)


# Cover 0<=t<=1/22 by rational intervals of width at most 1/5000.  At each
# midpoint, numerical roots merely propose a rational split polynomial G.
# Exact distance and coefficient inequalities prove |p_t-G|<|G| on the
# circle, so Rouche gives exactly one P-root in the disk throughout the cell.
cover_step = sp.Rational(1, 5000)
root_rounding_denominator = 10**8
left_endpoint = sp.Rational(0)
cover_count = 0
minimum_margin = None
while left_endpoint < t_max:
    right_endpoint = min(left_endpoint + cover_step, t_max)
    midpoint = (left_endpoint + right_endpoint) / 2
    midpoint_polynomial = sp.Poly(P4_expression.subs(t, midpoint), x)
    midpoint_leading = midpoint_polynomial.LC()
    numerical_roots = [
        complex(root)
        for root in sp.nroots(midpoint_polynomial, n=25, maxsteps=100)
    ]

    real_root = min(numerical_roots, key=lambda root: abs(root.imag))
    rational_roots: list[sp.Expr] = [
        sp.Rational(
            round(real_root.real * root_rounding_denominator),
            root_rounding_denominator,
        )
    ]
    upper_roots = sorted(
        (root for root in numerical_roots if root.imag > 1e-7),
        key=lambda root: root.real,
    )
    assert len(upper_roots) == 5
    for root in upper_roots:
        root_real = sp.Rational(
            round(root.real * root_rounding_denominator),
            root_rounding_denominator,
        )
        root_imag = sp.Rational(
            round(root.imag * root_rounding_denominator),
            root_rounding_denominator,
        )
        rational_roots.extend(
            [root_real + sp.I * root_imag, root_real - sp.I * root_imag]
        )
    assert len(rational_roots) == 11

    comparison = sp.Poly(sp.prod(x - root for root in rational_roots), x)
    distance_factors = []
    local_roots = 0
    for root in rational_roots:
        squared_distance = (
            sp.re(root) - center_real
        ) ** 2 + (sp.im(root) - center_imag) ** 2
        if squared_distance < radius**2:
            local_roots += 1
            distance_factors.append(
                radius - rational_sqrt_bound(squared_distance, upper=True)
            )
        else:
            distance_factors.append(
                rational_sqrt_bound(squared_distance, upper=False) - radius
            )
    assert local_roots == 1
    assert all(factor > 0 for factor in distance_factors)
    comparison_lower_bound = sp.prod(distance_factors)

    midpoint_error = sum(
        abs(
            midpoint_polynomial.coeff_monomial(x**power) / midpoint_leading
            - comparison.coeff_monomial(x**power)
        )
        * x_modulus_bound**power
        for power in range(11)
    )
    parameter_error = derivative_bound * (right_endpoint - left_endpoint) / 2
    margin = comparison_lower_bound - midpoint_error - parameter_error
    assert margin > 0
    if minimum_margin is None or margin < minimum_margin:
        minimum_margin = margin

    left_endpoint = right_endpoint
    cover_count += 1

assert cover_count == 228
assert minimum_margin is not None and minimum_margin > 0


# The whole localization disk maps strictly inside |y|<m/(m+1).  Squaring
# the triangle-inequality bound leaves a cubic whose coefficients are positive
# after m=n+22.
disk_difference = sp.together(
    (m**2 / (m + 1) - radius) ** 2
    - ((m + center_real) ** 2 + center_imag**2)
).as_numer_denom()[0]
disk_certificate = sp.Poly(sp.expand(disk_difference.subs(m, n + 22)), n)
assert all(coefficient > 0 for coefficient in disk_certificate.all_coeffs())


# Every point in the upper disk has pi < m*arg(y) < 2*pi.  For the lower
# bound use atan(q)>q/(1+q^2) and pi<22/7.  For the upper bound compare the
# disk with the ray at angle 2*pi/m, use the alternating Taylor bounds for
# sin and cos, and 333/106<pi<355/113.  Both reduce to positive polynomials
# after m=n+22.
q_lower = sp.Rational(37, 8) / (m - sp.Rational(423, 200))
lower_angle_bound = m * q_lower / (1 + q_lower**2) - sp.Rational(22, 7)
lower_angle_numerator = sp.together(lower_angle_bound).as_numer_denom()[0]
lower_angle_certificate = sp.Poly(
    sp.expand(lower_angle_numerator.subs(m, n + 22)), n
)
assert all(coefficient > 0 for coefficient in lower_angle_certificate.all_coeffs())

pi_lower = sp.Rational(333, 106)
pi_upper = sp.Rational(355, 113)
center_abscissa = m + center_real
upper_ray_bound = (
    center_abscissa
    * (2 * pi_lower / m - sp.Rational(4, 3) * pi_upper**3 / m**3)
    - center_imag
    * (
        1
        - 2 * pi_lower**2 / m**2
        + sp.Rational(2, 3) * pi_upper**4 / m**4
    )
    - radius
)
upper_ray_numerator = sp.together(upper_ray_bound).as_numer_denom()[0]
upper_ray_certificate = sp.Poly(
    sp.expand(upper_ray_numerator.subs(m, n + 22)), n
)
assert all(coefficient > 0 for coefficient in upper_ray_certificate.all_coeffs())


# Recover z=-B/A from the linear subresultant.  On the entire rectangular
# box containing the disk, certify Im(z)>0 by a tensor-product Bernstein
# expansion.  If Z=N/Den, then S=Im(N*conj(Den))=|Den|^2 Im(Z).  The factor
# Q=-S/(4*b*(3*t+2)) has strictly negative Bernstein coefficients, so S>0;
# in particular Den cannot vanish.
a, b, cube_a, cube_b, cube_t = sp.symbols(
    "a b A B T", real=True
)
linear_a = linear_subresultant4.coeff_monomial(z)
linear_b = linear_subresultant4.coeff_monomial(1)
z_branch = sp.cancel(
    (-linear_b / linear_a).subs({m: 1 / t, y: 1 + t * x})
)
z_numerator, z_denominator = sp.together(z_branch).as_numer_denom()


def split_real_imaginary(polynomial: sp.Expr) -> tuple[sp.Poly, sp.Poly]:
    """Evaluate x=a+ib by even/odd binomial splitting over QQ."""

    real_part = sp.Rational(0)
    imaginary_part = sp.Rational(0)
    for (x_power, t_power), coefficient in sp.Poly(
        polynomial, x, t, domain=QQ
    ).terms():
        for b_power in range(x_power + 1):
            term = (
                coefficient
                * sp.binomial(x_power, b_power)
                * a ** (x_power - b_power)
                * b**b_power
                * t**t_power
            )
            if b_power % 2 == 0:
                real_part += (-1) ** (b_power // 2) * term
            else:
                imaginary_part += (-1) ** ((b_power - 1) // 2) * term
    return (
        sp.Poly(real_part, a, b, t, domain=QQ),
        sp.Poly(imaginary_part, a, b, t, domain=QQ),
    )


numerator_real, numerator_imaginary = split_real_imaginary(z_numerator)
denominator_real, denominator_imaginary = split_real_imaginary(z_denominator)
imaginary_cross_product = (
    numerator_imaginary * denominator_real
    - numerator_real * denominator_imaginary
)
positive_factor = sp.Poly(4 * b * (3 * t + 2), a, b, t, domain=QQ)
bernstein_base = -imaginary_cross_product.exquo(positive_factor)

unit_cube_expression = bernstein_base.as_expr().subs(
    {
        a: sp.Rational(-623, 200) + cube_a,
        b: sp.Rational(37, 8) + cube_b,
        t: cube_t / 22,
    }
)
unit_cube_polynomial = sp.Poly(
    unit_cube_expression, cube_a, cube_b, cube_t, domain=QQ
)
assert unit_cube_polynomial.degree_list() == (16, 16, 30)


def power_to_bernstein(polynomial: sp.Poly) -> dict[tuple[int, ...], object]:
    """Exact tensor-product power-to-Bernstein conversion on [0,1]^d."""

    degrees = polynomial.degree_list()
    shape = tuple(degree + 1 for degree in degrees)
    coefficients = {monomial: coefficient for monomial, coefficient in polynomial.terms()}
    for axis, degree in enumerate(degrees):
        transformed = {}
        other_ranges = [range(shape[index]) for index in range(len(shape)) if index != axis]
        for other_indices in itertools.product(*other_ranges):
            for bernstein_index in range(degree + 1):
                value = QQ.zero
                for power_index in range(bernstein_index + 1):
                    full_index = list(other_indices)
                    full_index.insert(axis, power_index)
                    value += QQ(
                        math.comb(bernstein_index, power_index),
                        math.comb(degree, power_index),
                    ) * coefficients.get(tuple(full_index), QQ.zero)
                full_index = list(other_indices)
                full_index.insert(axis, bernstein_index)
                transformed[tuple(full_index)] = value
        coefficients = transformed
    return coefficients


bernstein_coefficients = power_to_bernstein(unit_cube_polynomial)
assert len(bernstein_coefficients) == 17 * 17 * 31
assert all(coefficient < 0 for coefficient in bernstein_coefficients.values())


print("PASS contact resultant r=4: exact endpoint gcds for 1<=m<=21")
print("PASS contact resultant r=4: Schur--Cohn inertia (9,2) for every m>=22")
print(
    "PASS contact resultant r=4: 228-cell Rouche localization and argument separation"
)
print("PASS contact resultant: uniform nonvanishing for every m and r=4")
