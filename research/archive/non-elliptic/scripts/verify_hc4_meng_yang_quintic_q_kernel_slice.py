#!/usr/bin/env python3
"""Verify a cone-compatible quintic Meng--Yang graph obstruction.

This checker starts from the v2 Meng--Yang potential and the plane normal-jet
equation.  It treats the first minimal trace slice on the ``partial_q``
constant-kernel chart left by the projective-gradient/Gordan--Noether gate.

The plane equations can be solved rationally in one parameter ``rho``.  On
the q-free lower slice, two repair-immune transverse coefficients reduce to
coprime cubics.  After adjoining ``a*p*q``, three such coefficients generate
the unit ideal over Q.  The remaining degree-at-most-two terms are then
classified from the unsolved plane faces: the linear q term is invisible,
the q-squared term violates the top-kernel chart, and the exceptional y*q
branch ends in two coprime transverse polynomials.  The finite-field census
records empty admissible good-prime fibers and bad-prime collapses; it is not
evidence for an HC(4) example.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import permutations
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_meng_yang_quintic_q_kernel_slice.json"
)

x, y, p, q, r = sp.symbols("x y p q r")
trace_symbol, trace_y, trace_p, trace_q, normal_symbol = sp.symbols(
    "trace_symbol trace_y trace_p trace_q normal_symbol"
)


# Derive the v2 plane determinant directly from the five-variable potential.
u = 1 + x * y
A = u**3 * p + 3 * x * u**2 * q - x**3 * r
B = (
    y**2 * u * (4 + 3 * x * y) * p
    + (y + 3 * x * y**2 * (4 + 3 * x * y)) * q
    + (2 * x - 3 * x**2 * y) * r
)
ambient_potential = A**2 + 13 * A + 2 * B
ambient_hessian = sp.hessian(ambient_potential, (x, y, p, q, r)).subs(
    {x: 0, r: trace_symbol}
)
source_unit = sp.Matrix([1, 0, 0, 0])
trace_gradient = sp.Matrix(
    [normal_symbol, trace_y, trace_p, trace_q]
)
plane_graph_hessian = (
    ambient_hessian[:4, :4]
    + 4
    * (
        source_unit * trace_gradient.T
        + trace_gradient * source_unit.T
    )
)
plane_graph_determinant = sp.factor(plane_graph_hessian.det())


# Independent explicit form of the normal-free plane determinant.  Keeping
# this formula avoids repeatedly expanding the ambient 5-by-5 Hessian.
plane_forcing = sp.expand(
    96 * trace_symbol * y
    + 64 * trace_p**2
    - 1024 * trace_p * trace_q * y
    - 1152 * trace_p * p * y
    + 192 * trace_p * q
    - 5696 * trace_p * y**3
    - 8736 * trace_p * y
    - 512 * trace_q**2 * p
    + 4096 * trace_q**2 * y**2
    + 128 * trace_q * trace_y
    - 1344 * trace_q * p**2
    + 4416 * trace_q * p * y**2
    - 8736 * trace_q * p
    + 45568 * trace_q * y**4
    + 69888 * trace_q * y**2
    + 192 * trace_y * p
    + 768 * trace_y * y**2
    + 1248 * trace_y
    - 864 * p**3
    - 1104 * p**2 * y**2
    - 11232 * p**2
    + 96 * p * q * y
    + 40800 * p * y**4
    + 38688 * p * y**2
    - 36504 * p
    + 384 * q * y**3
    + 624 * q * y
    + 126736 * y**6
    + 388752 * y**4
    + 298116 * y**2
)
assert sp.expand(
    plane_graph_determinant - plane_forcing + 64 * normal_symbol
) == 0
assert plane_graph_hessian[1:, 1:].det() == -8


def forcing_for_trace(trace: sp.Expr) -> sp.Expr:
    return sp.expand(
        plane_forcing.subs(
            {
                trace_symbol: trace,
                trace_y: sp.diff(trace, y),
                trace_p: sp.diff(trace, p),
                trace_q: sp.diff(trace, q),
            }
        )
    )


def homogeneous_coefficient(
    expression: sp.Expr, exponents: tuple[int, int, int]
) -> sp.Expr:
    return sp.factor(
        sp.Poly(expression, y, p, q).coeff_monomial(
            y ** exponents[0] * p ** exponents[1] * q ** exponents[2]
        )
    )


# The base trace slice retains every q-free coefficient through degree two,
# the linear q coefficient, and ``mixed_pq*p*q``.  The two remaining
# quadratic q terms are classified separately from the unsolved plane faces.
# The degree-three coefficient h is arbitrary; e is solved by the high faces.
rho, V, kappa, d_coefficient, h, e_coefficient, mixed_pq = sp.symbols(
    "rho V kappa d_coefficient h e_coefficient mixed_pq"
)
linear_q, mixed_yq, quadratic_q = sp.symbols(
    "linear_q mixed_yq quadratic_q"
)
c00, c10, c01, c20, c11, c02 = sp.symbols(
    "c00 c10 c01 c20 c11 c02"
)
lower_trace = (
    c00
    + c10 * y
    + c01 * p
    + c20 * y**2
    + c11 * y * p
    + c02 * p**2
    + linear_q * q
    + mixed_pq * p * q
)
d_relation = (16 * rho + 89 - V) / 2
trace_before_solving = (
    kappa * y**5
    + d_relation * y**3 * p
    + rho * y**2 * q
    + h * y**3
    + e_coefficient * y * p**2
    + lower_trace
)
forcing_before_solving = forcing_for_trace(trace_before_solving)
Q_polynomial = 160 * rho**2 + 1968 * rho + 6021
P_polynomial = 8 * rho**2 + 99 * rho + 279

assert sp.factor(
    homogeneous_coefficient(forcing_before_solving, (6, 0, 0))
    - 16 * (V**2 + 2 * kappa * (20 * rho + 123))
) == 0
assert sp.factor(
    homogeneous_coefficient(forcing_before_solving, (4, 1, 0))
    - 16
    * (
        32 * V * mixed_pq
        + 40 * mixed_pq * kappa
        - 8 * V * e_coefficient
        - 12 * V * rho
        - 39 * V
        + 60 * kappa
        + Q_polynomial
    )
) == 0

# The quartic first normal jet is forcing_4/64.  On the partial_q kernel
# chart its y^3*q coefficient must vanish, giving
# (2*mixed_pq+3)*V=P(rho).
quartic_normal_q = homogeneous_coefficient(
    forcing_before_solving / 64, (3, 0, 1)
)
assert sp.factor(
    quartic_normal_q
    + ((2 * mixed_pq + 3) * V - P_polynomial) / 2
) == 0


# Audit the two omitted quadratic q terms before solving V,kappa,e.  The
# degree-five plane faces force mixed_yq*(4V+5kappa)=0 and
# quadratic_q*(4V+5kappa)=0.  More decisively on the partial_q top-kernel
# chart, the quartic normal jet must be q-independent, while its y^2*q^2
# coefficient is the square 256*quadratic_q^2.  Hence quadratic_q=0 in
# characteristic zero.  The y*q branch is handled after the generic ideal.
full_quadratic_trace = (
    trace_before_solving + mixed_yq * y * q + quadratic_q * q**2
)
full_quadratic_forcing = forcing_for_trace(full_quadratic_trace)
assert sp.factor(
    homogeneous_coefficient(full_quadratic_forcing, (5, 0, 0))
    - 128 * mixed_yq * (4 * V + 5 * kappa)
) == 0
assert sp.factor(
    homogeneous_coefficient(full_quadratic_forcing, (4, 0, 1))
    - 256 * quadratic_q * (4 * V + 5 * kappa)
) == 0
assert sp.factor(
    homogeneous_coefficient(full_quadratic_forcing / 64, (3, 0, 1))
    + (
        (2 * mixed_pq + 3) * V
        - 512 * quadratic_q * mixed_yq
        - P_polynomial
    )
    / 2
) == 0
assert sp.factor(
    homogeneous_coefficient(full_quadratic_forcing / 64, (2, 0, 2))
    - 256 * quadratic_q**2
) == 0


# Boundary audit before division.  If P=V=0, the plane equations force
# kappa=0 and Q=0, but P and Q are coprime.  The denominator
# 20*rho+123 cannot vanish on a surviving branch.  The separate
# 2*mixed_pq+3=0 boundary is checked after the generic transverse ideal.
boundary_resultant = sp.resultant(
    P_polynomial, Q_polynomial, rho
)
assert boundary_resultant == 16959456
assert sp.factor(
    P_polynomial.subs(rho, -sp.Rational(123, 20))
) == -sp.Rational(2727, 100)


kernel_denominator = 2 * mixed_pq + 3
V_solution = P_polynomial / kernel_denominator
kappa_solution = sp.factor(
    -V_solution**2 / (2 * (20 * rho + 123))
)
e_solution = sp.factor(
    (
        32 * V_solution * mixed_pq
        + 40 * mixed_pq * kappa_solution
        - 12 * V_solution * rho
        - 39 * V_solution
        + 60 * kappa_solution
        + Q_polynomial
    )
    / (8 * V_solution)
)
solution_substitution = {
    V: V_solution,
    kappa: kappa_solution,
    e_coefficient: e_solution,
}
trace = sp.factor(trace_before_solving.subs(solution_substitution))
forcing = sp.cancel(forcing_before_solving.subs(solution_substitution))
assert sp.Poly(sp.together(forcing).as_numer_denom()[0], y, p, q).total_degree() <= 4

target_constant = sp.Symbol("target_constant")
normal_jet = sp.cancel((forcing - target_constant) / 64)
assert sp.Poly(
    sp.together(normal_jet).as_numer_denom()[0], y, p, q
).total_degree() <= 4


# Exact x^3 jet of the graph pullback R=T+x*S.  A remaining x^2*U changes
# the first transverse determinant by -192*U; deg(U)<=3, so it cannot touch
# the three higher-degree witnesses extracted below.
a0 = p
a1 = 3 * y * p + 3 * q
a2 = 3 * y**2 * p + 6 * y * q
b0 = 4 * y**2 * p + y * q


def potential_jet_for(
    graph_trace: sp.Expr, graph_normal: sp.Expr
) -> sp.Expr:
    """Return the x^3 jet needed for the first transverse determinant."""

    a3 = y**3 * p + 3 * y**2 * q - graph_trace
    b1 = 7 * y**3 * p + 12 * y**2 * q + 2 * graph_trace
    b2 = (
        3 * y**4 * p
        + 9 * y**3 * q
        + 2 * graph_normal
        - 3 * y * graph_trace
    )
    b3 = -3 * y * graph_normal
    return (
        a0**2
        + 13 * a0
        + 2 * b0
        + x * (2 * a0 * a1 + 13 * a1 + 2 * b1)
        + x**2
        * (a1**2 + 2 * a0 * a2 + 13 * a2 + 2 * b2)
        + x**3
        * (2 * a0 * a3 + 2 * a1 * a2 + 13 * a3 + 2 * b3)
    )


potential_jet = potential_jet_for(trace, normal_jet)
assert -24 * 1 * 2**3 == -192


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def product_y_coefficient(
    factors: list[sp.Expr], target_degree: int
) -> sp.Expr:
    """Return one y coefficient without expanding the complete product."""

    coefficients = [sp.S.One] + [sp.S.Zero] * target_degree
    for factor in factors:
        factor_terms = {
            exponent[0]: coefficient
            for exponent, coefficient in sp.Poly(factor, y).terms()
            if exponent[0] <= target_degree
        }
        product = [sp.S.Zero] * (target_degree + 1)
        for left_degree, left_coefficient in enumerate(coefficients):
            for right_degree, right_coefficient in factor_terms.items():
                degree = left_degree + right_degree
                if degree <= target_degree:
                    product[degree] += left_coefficient * right_coefficient
        coefficients = product
    return coefficients[target_degree]


def transverse_witnesses(matrix: sp.Matrix) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Extract [x*y^7], [x*y^5*p], and [x*y^4*q] from det(H)."""

    hessian_zero = matrix.subs(x, 0)
    hessian_linear = matrix.diff(x).subs(x, 0)
    y7 = sp.S.Zero
    y5p = sp.S.Zero
    y4q = sp.S.Zero
    for permutation in permutations(range(4)):
        sign = permutation_sign(permutation)
        for linear_row in range(4):
            factors = [
                (
                    hessian_linear[row, permutation[row]]
                    if row == linear_row
                    else hessian_zero[row, permutation[row]]
                )
                for row in range(4)
            ]
            values = [
                factor.subs({p: 0, q: 0})
                for factor in factors
            ]
            y7 += sign * product_y_coefficient(values, 7)
            for differentiated in range(4):
                q_factors = list(values)
                q_factors[differentiated] = sp.diff(
                    factors[differentiated], q
                ).subs({p: 0, q: 0})
                y4q += sign * product_y_coefficient(q_factors, 4)
                p_factors = list(values)
                p_factors[differentiated] = sp.diff(
                    factors[differentiated], p
                ).subs({p: 0, q: 0})
                y5p += sign * product_y_coefficient(p_factors, 5)
    return sp.factor(y7), sp.factor(y5p), sp.factor(y4q)


jet_hessian = sp.hessian(potential_jet, (x, y, p, q))
y7_coefficient, y5p_coefficient, y4q_coefficient = transverse_witnesses(
    jet_hessian
)

# Normalize three repair-invariant coefficients with both linear_q*q and
# mixed_pq*p*q retained.  Their polynomial ideal over Q is the unit ideal;
# linear_q cancels from every equation.
denominator_linear = 20 * rho + 123
pq_equation_7 = sp.factor(
    y7_coefficient
    * kernel_denominator**2
    * denominator_linear
    / (24 * P_polynomial)
)
pq_equation_6 = sp.factor(
    -y5p_coefficient
    * kernel_denominator**2
    * denominator_linear
    * P_polynomial
    / 24
)
pq_equation_5 = sp.factor(-y4q_coefficient / 48)
for equation in (pq_equation_7, pq_equation_6, pq_equation_5):
    assert sp.denom(sp.cancel(equation)) == 1

pq_equation_7_expected = (
    1280 * mixed_pq * rho**3
    + 23232 * mixed_pq * rho**2
    + 156744 * mixed_pq * rho
    + 383022 * mixed_pq
    + 1968 * rho**3
    + 33498 * rho**2
    + 212733 * rho
    + 506736
)
assert sp.expand(pq_equation_7 - pq_equation_7_expected) == 0
pq_equation_5_expected = (
    64 * mixed_pq * rho**2
    + 768 * mixed_pq * rho
    + 3114 * mixed_pq
    - 16 * rho**3
    - 192 * rho**2
    - 801 * rho
    - 54
)
assert sp.expand(pq_equation_5 - pq_equation_5_expected) == 0
assert sp.Poly(pq_equation_7, mixed_pq).degree() == 1
assert sp.Poly(pq_equation_6, mixed_pq).degree() == 3
assert sp.Poly(pq_equation_5, mixed_pq).degree() == 1

cancelled_parameters = {
    h,
    c00,
    c10,
    c01,
    c20,
    c11,
    c02,
    linear_q,
    target_constant,
}
for equation in (pq_equation_7, pq_equation_6, pq_equation_5):
    assert not equation.free_symbols & cancelled_parameters

pq_resultant_5 = sp.factor(
    sp.resultant(pq_equation_7, pq_equation_5, mixed_pq)
)
pq_resultant_6 = sp.factor(
    sp.resultant(pq_equation_7, pq_equation_6, mixed_pq)
)
pq_resultant_gcd = sp.gcd(
    sp.Poly(pq_resultant_5, rho),
    sp.Poly(pq_resultant_6, rho),
).monic()
shared_resultant_quadratic = 32 * rho**2 + 384 * rho + 1557
assert pq_resultant_gcd == sp.Poly(
    shared_resultant_quadratic, rho
).monic()
pq_basis = sp.groebner(
    [pq_equation_7, pq_equation_6, pq_equation_5],
    mixed_pq,
    rho,
    order="lex",
)
assert [polynomial.as_expr() for polynomial in pq_basis.polys] == [1]


# The parameterization divided by kernel_denominator=2*mixed_pq+3.  On its
# zero chart the kernel equation first forces P(rho)=0.  The V=0 sub-branch
# is already excluded by resultant(P,Q).  For V nonzero, solve the two plane
# faces directly and extract one transverse coefficient from the original
# jet.  Modulo P it is 540*(25*rho+141), and the resultant is nonzero.
boundary_V = sp.Symbol("boundary_V", nonzero=True)
boundary_mixed_pq = -sp.Rational(3, 2)
boundary_kappa = -boundary_V**2 / (2 * denominator_linear)
boundary_e = (
    Q_polynomial / (8 * boundary_V) - (12 * rho + 87) / 8
)
boundary_d = (16 * rho + 89 - boundary_V) / 2
boundary_trace = (
    boundary_kappa * y**5
    + boundary_d * y**3 * p
    + rho * y**2 * q
    + h * y**3
    + boundary_e * y * p**2
    + boundary_mixed_pq * p * q
    + c00
    + c10 * y
    + c01 * p
    + c20 * y**2
    + c11 * y * p
    + c02 * p**2
    + linear_q * q
)
boundary_forcing = forcing_for_trace(boundary_trace)
boundary_normal = sp.together((boundary_forcing - target_constant) / 64)
boundary_potential_jet = potential_jet_for(
    boundary_trace, boundary_normal
)
boundary_y7, boundary_y5p, boundary_y4q = transverse_witnesses(
    sp.hessian(boundary_potential_jet, (x, y, p, q))
)
for equation in (boundary_y7, boundary_y5p, boundary_y4q):
    assert not equation.free_symbols & cancelled_parameters
boundary_y4q_expected = 48 * (
    80 * rho**3 + 1464 * rho**2 + 8937 * rho + 18117
)
assert sp.factor(boundary_y4q - boundary_y4q_expected) == 0
boundary_y4q_remainder = sp.rem(
    sp.Poly(boundary_y4q_expected, rho),
    sp.Poly(P_polynomial, rho),
).as_expr()
assert sp.factor(
    boundary_y4q_remainder - 540 * (25 * rho + 141)
) == 0
kernel_boundary_resultant = sp.resultant(
    P_polynomial, 25 * rho + 141, rho
)
assert kernel_boundary_resultant == -15552


# It remains to classify mixed_yq != 0 after quadratic_q=0.  The degree-five
# plane face gives 4*V+5*kappa=0.  Together with the degree-six face this has
# V=0 or V=8*(20*rho+123)/5.  The first case makes the kernel and second plane
# equations P=Q=0, already excluded by boundary_resultant.  On the second
# case, the kernel equation and the second plane face determine mixed_pq and
# e.  Two high transverse witnesses are independent of mixed_yq and of every
# other lower coefficient, and have nonzero resultant.
yq_V = 8 * denominator_linear / 5
yq_kappa = -32 * denominator_linear / 25
yq_mixed_pq = sp.factor(
    (40 * rho**2 + 15 * rho - 1557)
    / (16 * denominator_linear)
)
yq_e = sp.factor(
    (
        32 * yq_V * yq_mixed_pq
        + 40 * yq_mixed_pq * yq_kappa
        - 12 * yq_V * rho
        - 39 * yq_V
        + 60 * yq_kappa
        + Q_polynomial
    )
    / (8 * yq_V)
)
assert sp.factor(4 * yq_V + 5 * yq_kappa) == 0
assert sp.factor(
    yq_V**2 + 2 * yq_kappa * denominator_linear
) == 0
assert sp.factor(
    (2 * yq_mixed_pq + 3) * yq_V - P_polynomial
) == 0
yq_substitution = {
    V: yq_V,
    kappa: yq_kappa,
    mixed_pq: yq_mixed_pq,
    e_coefficient: yq_e,
}
yq_trace = sp.factor(
    trace_before_solving.subs(yq_substitution) + mixed_yq * y * q
)
yq_forcing = sp.cancel(forcing_for_trace(yq_trace))
assert sp.Poly(yq_forcing, y, p, q).total_degree() <= 4
yq_forcing_top = sum(
    coefficient * y**exponents[0] * p**exponents[1] * q**exponents[2]
    for exponents, coefficient in sp.Poly(yq_forcing, y, p, q).terms()
    if sum(exponents) == 4
)
assert sp.factor(sp.diff(yq_forcing_top, q)) == 0
yq_normal = sp.cancel((yq_forcing - target_constant) / 64)
yq_y7, yq_y5p, yq_y4q = transverse_witnesses(
    sp.hessian(potential_jet_for(yq_trace, yq_normal), (x, y, p, q))
)
yq_F = 160 * rho**2 + 1968 * rho + 5841
yq_G = (
    3968000 * rho**4
    + 98457600 * rho**3
    + 896972480 * rho**2
    + 3553345920 * rho
    + 5132433339
)
assert sp.factor(
    yq_y7 - 192 * denominator_linear * yq_F / 25
) == 0
assert sp.factor(
    yq_y5p - 3 * yq_G / (25 * denominator_linear)
) == 0
for equation in (yq_y7, yq_y5p, yq_y4q):
    assert not equation.free_symbols & (
        cancelled_parameters | {mixed_yq}
    )
yq_resultant = sp.resultant(yq_F, yq_G, rho)
assert yq_resultant == 986335129354383654912000
assert sp.gcd(sp.Poly(yq_F, rho), sp.Poly(yq_G, rho)) == sp.Poly(
    1, rho
)
yq_resultant_factorization = sp.factorint(yq_resultant)
assert yq_resultant_factorization == {
    2: 27,
    3: 8,
    5: 3,
    11: 1,
    24223: 1,
    33629: 1,
}


# On the zero-pq sub-slice, the first and third equations reduce to the two
# cubics from the original calculation.
A_cubic = 656 * rho**3 + 11166 * rho**2 + 70911 * rho + 168912
B_cubic = 16 * rho**3 + 192 * rho**2 + 801 * rho + 54
assert sp.factor(
    y7_coefficient.subs(mixed_pq, 0)
    - 8 * P_polynomial * A_cubic / denominator_linear
) == 0
assert sp.factor(
    y4q_coefficient.subs(mixed_pq, 0) - 48 * B_cubic
) == 0


terminal_resultant = sp.resultant(A_cubic, B_cubic, rho)
assert terminal_resultant == -108117004020524928
assert sp.gcd(
    sp.Poly(A_cubic, rho), sp.Poly(B_cubic, rho)
) == sp.Poly(1, rho)
terminal_factorization = sp.factorint(abs(terminal_resultant))
assert terminal_factorization == {
    2: 7,
    3: 17,
    11: 1,
    13: 1,
    53: 1,
    863: 1,
}


def common_roots(polynomials: tuple[sp.Expr, ...], prime: int) -> list[int]:
    return [
        value
        for value in range(prime)
        if all(
            int(polynomial.subs(rho, value)) % prime == 0
            for polynomial in polynomials
        )
    ]


good_prime_roots = {
    str(prime): common_roots((A_cubic, B_cubic), prime)
    for prime in (101, 103)
}
assert good_prime_roots == {"101": [], "103": []}
bad_prime_roots = {
    str(prime): common_roots((A_cubic, B_cubic), prime)
    for prime in (11, 13, 53, 863)
}
assert bad_prime_roots == {
    "11": [2],
    "13": [10],
    "53": [34],
    "863": [717],
}
for prime_text, roots in bad_prime_roots.items():
    prime = int(prime_text)
    assert all(
        int(P_polynomial.subs(rho, value)) % prime != 0
        and (20 * value + 123) % prime != 0
        for value in roots
    )

yq_good_prime_roots = {
    str(prime): common_roots((yq_F, yq_G), prime)
    for prime in (101, 103)
}
assert yq_good_prime_roots == {"101": [], "103": []}
yq_bad_prime_roots = {
    str(prime): common_roots((yq_F, yq_G), prime)
    for prime in (11, 24223, 33629)
}
assert yq_bad_prime_roots == {
    "11": [0],
    "24223": [4365],
    "33629": [30101],
}


pq_equation_polys = [
    sp.Poly(equation, mixed_pq, rho)
    for equation in (pq_equation_7, pq_equation_6, pq_equation_5)
]


def evaluate_mod_p(
    polynomial: sp.Poly, pq_value: int, rho_value: int, prime: int
) -> int:
    return sum(
        int(coefficient)
        * pow(pq_value, exponents[0], prime)
        * pow(rho_value, exponents[1], prime)
        for exponents, coefficient in polynomial.terms()
    ) % prime


def pq_common_points(prime: int, *, admissible_only: bool) -> list[list[int]]:
    points: list[list[int]] = []
    for rho_value in range(prime):
        if admissible_only and (
            int(P_polynomial.subs(rho, rho_value)) % prime == 0
            or int(denominator_linear.subs(rho, rho_value)) % prime == 0
        ):
            continue
        for pq_value in range(prime):
            if admissible_only and (2 * pq_value + 3) % prime == 0:
                continue
            if all(
                evaluate_mod_p(
                    polynomial, pq_value, rho_value, prime
                )
                == 0
                for polynomial in pq_equation_polys
            ):
                points.append([rho_value, pq_value])
    return points


pq_prime_counts = {}
for prime in (101, 103):
    raw_points = pq_common_points(prime, admissible_only=False)
    admissible_points = pq_common_points(prime, admissible_only=True)
    pq_prime_counts[str(prime)] = {
        "raw": len(raw_points),
        "admissible": len(admissible_points),
    }
assert pq_prime_counts == {
    "101": {"raw": 1, "admissible": 0},
    "103": {"raw": 0, "admissible": 0},
}


# Independent full-potential replay at exact rational parameter points.  This
# does not prove a symbolic identity by sampling; it checks that the x^3-jet
# extraction agrees with fresh Hessian calculations in the original
# five-variable formula on both transverse branches.
def full_potential_witnesses(
    graph_trace: sp.Expr, graph_normal: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    sample_graph = graph_trace + x * graph_normal
    sample_potential = sp.expand(ambient_potential.subs(r, sample_graph))
    sample_hessian = sp.hessian(sample_potential, (x, y, p, q))
    sample_hessian_zero = sample_hessian.subs(x, 0)
    sample_hessian_linear = sample_hessian.diff(x).subs(x, 0)
    sample_transverse = sp.cancel(
        sum(
            sample_hessian_zero.adjugate()[column, row]
            * sample_hessian_linear[row, column]
            for row in range(4)
            for column in range(4)
        )
    )
    sample_y7 = sp.factor(
        sp.Poly(
            sample_transverse.subs({p: 0, q: 0}), y
        ).coeff_monomial(y**7)
    )
    sample_y4q = sp.factor(
        sp.Poly(
            sp.diff(sample_transverse, q).subs({p: 0, q: 0}), y
        ).coeff_monomial(y**4)
    )
    sample_y5p = sp.factor(
        sp.Poly(
            sp.diff(sample_transverse, p).subs({p: 0, q: 0}), y
        ).coeff_monomial(y**5)
    )
    return sample_y7, sample_y5p, sample_y4q


sample_substitution = {
    rho: 0,
    h: 2,
    c00: 1,
    c10: 2,
    c01: 3,
    c20: 5,
    c11: 7,
    c02: 11,
    linear_q: 13,
    mixed_pq: 4,
    target_constant: 17,
}
sample_trace = sp.factor(trace.subs(sample_substitution))
sample_normal = sp.factor(normal_jet.subs(sample_substitution))
sample_y7, sample_y5p, sample_y4q = full_potential_witnesses(
    sample_trace, sample_normal
)
assert sp.factor(
    sample_y7 - y7_coefficient.subs(sample_substitution)
) == 0
assert sp.factor(
    sample_y5p - y5p_coefficient.subs(sample_substitution)
) == 0
assert sp.factor(
    sample_y4q - y4q_coefficient.subs(sample_substitution)
) == 0

yq_sample_substitution = {
    rho: 0,
    h: 2,
    c00: 1,
    c10: 2,
    c01: 3,
    c20: 5,
    c11: 7,
    c02: 11,
    linear_q: 13,
    mixed_yq: 19,
    target_constant: 17,
}
yq_sample_trace = sp.factor(yq_trace.subs(yq_sample_substitution))
yq_sample_normal = sp.factor(yq_normal.subs(yq_sample_substitution))
yq_sample_witnesses = full_potential_witnesses(
    yq_sample_trace, yq_sample_normal
)
for replayed, symbolic in zip(
    yq_sample_witnesses,
    (yq_y7, yq_y5p, yq_y4q),
):
    assert sp.factor(
        replayed - symbolic.subs(yq_sample_substitution)
    ) == 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = {
        "schema": "hc4-meng-yang-quintic-q-kernel-slice.v3",
        "status": "exact characteristic-zero exclusion; no HC4 candidate",
        "kernel_chart": "partial_q",
        "trace_slice": (
            "kappa*y^5+d*y^3*p+rho*y^2*q+h*y^3+e*y*p^2+"
            "T_le2(y,p,q)"
        ),
        "terminal_cubics": {
            "A": [656, 11166, 70911, 168912],
            "B": [16, 192, 801, 54],
        },
        "terminal_resultant": int(terminal_resultant),
        "terminal_resultant_factorization": {
            str(prime): exponent
            for prime, exponent in terminal_factorization.items()
        },
        "good_prime_common_roots": good_prime_roots,
        "bad_prime_common_roots": bad_prime_roots,
        "pq_extension": {
            "witness_bidegrees_in_a_rho": [[1, 3], [3, 7], [1, 3]],
            "resultant_gcd": "32*rho^2+384*rho+1557",
            "groebner_basis": ["1"],
            "finite_field_counts": pq_prime_counts,
            "kernel_denominator_boundary_resultant": int(
                kernel_boundary_resultant
            ),
            "linear_q_coefficient_cancels": True,
        },
        "complete_degree_two_extension": {
            "q_squared_forced_zero_by_quartic_normal_square": 256,
            "nonzero_yq_branch": {
                "witness_F": [160, 1968, 5841],
                "witness_G": [
                    3968000,
                    98457600,
                    896972480,
                    3553345920,
                    5132433339,
                ],
                "resultant": int(yq_resultant),
                "resultant_factorization": {
                    str(prime): exponent
                    for prime, exponent in yq_resultant_factorization.items()
                },
                "good_prime_common_roots": yq_good_prime_roots,
                "bad_prime_common_roots": yq_bad_prime_roots,
            },
        },
        "rational_reconstruction": (
            "empty: every degree-at-most-two lower trace branch is "
            "excluded over Q"
        ),
        "open": [
            "broader cubic and quartic trace terms on the partial_q chart",
            "the remaining finite and infinite constant-kernel charts",
            "lower transverse layers after a surviving high transverse ideal",
        ],
    }
    canonical = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(canonical.encode()).hexdigest()
    if args.output is not None:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(canonical)
        print(f"WROTE {output}")
        print(f"SHA256 {digest}")

    print("PASS: the plane determinant is unit-affine in the normal jet")
    print("PASS: the partial_q top-kernel plane equations are parameterized")
    print("PASS: the q and p*q lower terms leave a unit transverse ideal")
    print("PASS: the q^2 lower term violates the partial_q top-kernel chart")
    print("PASS: the exceptional y*q branch has coprime transverse witnesses")
    print("PASS: admissible good-prime fibers are empty")
    print("PASS: the original full potential independently replays both branches")


if __name__ == "__main__":
    main()
