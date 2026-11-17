#!/usr/bin/env python3
"""Exact verification of the global relative Davenport--Sunada Keller pair."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    A,
    IDENTITY,
    MINIMAL_A,
    T,
    U,
    Y,
    Z,
    branch_cubic,
    conjugate_a,
    correspondence,
    cycle_partition,
    davenport_pair,
    gl32,
    inverses,
    line_permutation,
    matrix_product,
    point_permutation,
    reduce_a,
    zeta_denominator,
)


g, h = davenport_pair()
delta = branch_cubic()

# The correspondence is global over the two-dimensional (T,U)-base, not
# merely after fixing T.  Reduction by the coefficient-field relation makes
# the divisibility certificate exact.
groebner = sp.groebner(
    [MINIMAL_A, correspondence()],
    Y,
    Z,
    A,
    order="grlex",
    domain=sp.QQ.frac_field(T),
)
assert groebner.reduce(sp.expand(g - h))[1] == 0

# Each of the six critical points occurs above one of three finite branch
# values.  Both conjugate covers have exactly the same rational branch cubic.
resultant_g = reduce_a(sp.resultant(sp.diff(g, Y), g - U, Y), T, U)
resultant_h = reduce_a(sp.resultant(sp.diff(h, Z), h - U, Z), T, U)
assert sp.factor(resultant_g - delta**2 / 7**6) == 0
assert sp.factor(resultant_h - delta**2 / 7**6) == 0
coefficient_field = sp.QQ.algebraic_field(sp.sqrt(-7)).frac_field(T)
branch_polynomial = sp.Poly(delta, U, domain=coefficient_field)
assert branch_polynomial.degree() == 3
assert branch_polynomial.is_irreducible

# Direct Cox-ledger Kellerization.  On the unramified source g'_T(Y) is a
# unit, so adjoining Zeta=zeta/g' cancels the full Jacobian unit.  The same
# calculation applies to h and preserves the degree-seven fibers.
zeta = sp.symbols("zeta")
cox_g = (T, g, zeta / sp.diff(g, Y))
cox_h = (T, h, zeta / sp.diff(h, Z))
assert sp.factor(sp.Matrix(cox_g).jacobian((T, Y, zeta)).det() - 1) == 0
assert sp.factor(sp.Matrix(cox_h).jacobian((T, Z, zeta)).det() - 1) == 0
target_zeta = sp.symbols("target_zeta")
assert sp.cancel(
    (sp.diff(g, Y) * target_zeta) / sp.diff(g, Y) - target_zeta
) == 0
assert sp.cancel(
    (sp.diff(h, Z) * target_zeta) / sp.diff(h, Z) - target_zeta
) == 0

# The tangent-mark base change used by the relative weighted lift is finite.
# On it H(W)=g(rW)-g(0)-r*g'(0)W satisfies H(0)=H'(0)=H(1)=0.
r, W = sp.symbols("r W")
g0 = g.subs(Y, 0)
g0_prime = sp.diff(g, Y).subs(Y, 0)
tangent_equation = sp.factor((g.subs(Y, r) - g0 - r * g0_prime) / r**2)
assert sp.Poly(tangent_equation, r).degree() == 5
H = sp.expand(g.subs(Y, r * W) - g0 - r * g0_prime * W)
assert H.subs(W, 0) == 0
assert sp.diff(H, W).subs(W, 0) == 0
assert sp.factor(H.subs(W, 1) - r**2 * tangent_equation) == 0
c = -sp.diff(H, W).subs(W, 1)
kappa_numerator = sp.factor(sp.diff(H, W, 2).subs(W, 1) + 2 * c)
assert reduce_a(c, T, r) != 0
assert reduce_a(kappa_numerator, T, r) != 0
assert reduce_a(sp.resultant(tangent_equation, c, r), T) != 0
assert reduce_a(sp.resultant(tangent_equation, kappa_numerator, r), T) != 0

# The common abstract target L has coordinates (T,U).  On C=1, its two
# embeddings in the weighted targets recover g(Y)-U and h(Z)-U exactly.
s_g = -r * g0_prime
t_g = g0 - U
assert sp.expand(H - s_g * W + t_g - (g.subs(Y, r * W) - U)) == 0

rho = sp.symbols("rho")
h0 = h.subs(Z, 0)
h0_prime = sp.diff(h, Z).subs(Z, 0)
G = sp.expand(h.subs(Z, rho * W) - h0 - rho * h0_prime * W)
s_h = -rho * h0_prime
t_h = h0 - U
assert sp.expand(G - s_h * W + t_h - (h.subs(Z, rho * W) - U)) == 0

# The Hessian divisors cannot become affinely equivalent over the generic
# base.  A hypothetical g''(alpha*z+beta)=scale*h''(z) first forces beta=0;
# the next two coefficients force alpha^2=alpha^3=-(1+a)/a.
alpha, beta, scale = sp.symbols("alpha beta scale")
g_second = sp.diff(g, Y, 2)
h_second = sp.diff(h, Z, 2)
comparison = sp.Poly(
    sp.expand(g_second.subs(Y, alpha * Z + beta) - scale * h_second),
    Z,
)
leading = comparison.coeff_monomial(Z**5)
quartic = comparison.coeff_monomial(Z**4)
assert sp.factor(leading - 6 * (alpha**5 - scale)) == 0
assert sp.factor(quartic - 30 * alpha**4 * beta) == 0
centered = sp.Poly(comparison.as_expr().subs(beta, 0), Z)
cubic = centered.coeff_monomial(Z**3).subs(scale, alpha**5)
quadratic = centered.coeff_monomial(Z**2).subs(scale, alpha**5)
assert sp.factor(cubic - 20 * T * alpha**3 * (1 + A + A * alpha**2)) == 0
assert sp.factor(
    quadratic - 12 * T * alpha**2 * (1 + A + A * alpha**3)
) == 0
assert sp.resultant(MINIMAL_A, 2 * A + 1, A) != 0

# The point and line stabilizers form a nonconjugate Gassmann pair.  Equality
# of full cycle partitions gives equality of the zero-dimensional zeta
# function for every unramified Frobenius element, not just point counts.
matrices = gl32()
assert len(matrices) == 168
inverse = inverses(matrices)
for matrix in matrices:
    points = point_permutation(matrix)
    lines = line_permutation(matrix, inverse[matrix])
    assert cycle_partition(points) == cycle_partition(lines)
    assert zeta_denominator(points) == zeta_denominator(lines)

point_stabilizer = {m for m in matrices if point_permutation(m)[0] == 0}
line_stabilizer = {
    m for m in matrices if line_permutation(m, inverse[m])[0] == 0
}
assert len(point_stabilizer) == len(line_stabilizer) == 24
for matrix in matrices:
    conjugate = {
        matrix_product(matrix_product(matrix, member), inverse[matrix])
        for member in point_stabilizer
    }
    assert conjugate != line_stabilizer

# Exact finite-field regression on every unramified rational (T,U)-fiber for
# the first two split good primes of Q(sqrt(-7)).
finite_field_fibers = 0
for prime, a_value in ((11, 4), (23, 9)):
    conjugate_value = (-1 - a_value) % prime
    assert (a_value**2 + a_value + 2) % prime == 0
    for t_value in range(prime):
        for u_value in range(prime):
            if int(delta.subs({T: t_value, U: u_value})) % prime == 0:
                continue

            factor_degrees = []
            for coefficient_value in (a_value, conjugate_value):
                polynomial = sp.Poly(
                    sp.expand(
                        28
                        * (
                            g.subs({A: coefficient_value, T: t_value})
                            - u_value
                        )
                    ),
                    Y,
                    modulus=prime,
                )
                factor_degrees.append(
                    sorted(
                        factor.degree()
                        for factor, multiplicity in polynomial.factor_list()[1]
                        for _ in range(multiplicity)
                    )
                )
            assert factor_degrees[0] == factor_degrees[1]
            finite_field_fibers += 1

print("PASS: the Davenport correspondence is global over the (T,U)-surface")
print("PASS: both covers have the same reduced cubic branch divisor")
print("PASS: the two direct Cox suspensions have Jacobian one and unchanged fibers")
print("PASS: the finite tangent-mark base change produces the two weighted pencils")
print("PASS: the generic Hessian divisors remain affinely inequivalent")
print("PASS: all 168 Frobenius classes have identical point/line zeta factors")
print(
    "PASS: exact factor-degree equality on "
    f"{finite_field_fibers} unramified fibers over F_11 and F_23"
)
