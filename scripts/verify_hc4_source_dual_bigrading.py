#!/usr/bin/env python3
"""Verify the 2+2 HC4 bigrading and weighted-cone identities.

The script checks the finite algebra used in HC4_SOURCE_DUAL_BIGRADING.md:

* the quartic and sextic source/dual monomial ledgers;
* covariance of Hessian faces under weighted diagonal scaling;
* the cotangent block determinant on the dual-linear JC(2) locus;
* the successive Schur-cone recursion through order eight; and
* the rotating-cone model (x*t+y*m)^2.

The synchronization lemma itself is an all-order identity over a filtered
domain.  The order-eight calculation is an independent symbolic regression
of its displayed recursion, not a bounded substitute for the proof.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


def weak_compositions(total: int, length: int):
    """Yield weak compositions of total into length parts."""
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first, *tail)


# 1. With source variables (x,y) and dual variables (t,m), the number of
# monomials of bidegree (a,b) is (a+1)(b+1).
quartic_ledger = {
    (source_degree, 4 - source_degree):
    (source_degree + 1) * (5 - source_degree)
    for source_degree in range(5)
}
sextic_ledger = {
    (source_degree, 6 - source_degree):
    (source_degree + 1) * (7 - source_degree)
    for source_degree in range(7)
}
assert quartic_ledger == {
    (4, 0): 5,
    (3, 1): 8,
    (2, 2): 9,
    (1, 3): 8,
    (0, 4): 5,
}
assert sextic_ledger == {
    (6, 0): 7,
    (5, 1): 12,
    (4, 2): 15,
    (3, 3): 16,
    (2, 4): 15,
    (1, 5): 12,
    (0, 6): 7,
}
assert sum(quartic_ledger.values()) == 35
assert sum(sextic_ledger.values()) == 84
assert quartic_ledger[(4, 0)] + quartic_ledger[(3, 1)] == 13
assert sextic_ledger[(6, 0)] + sextic_ledger[(5, 1)] == 19


# 2. Weighted Hessian covariance.  For a monomial z^alpha of rho-weight d,
# D_rho Hess(z^alpha)(lambda^rho z) D_rho = lambda^d Hess(z^alpha)(z).
# Check every quartic and sextic monomial for several nontrivial weights.
x, y, t, m, lam = sp.symbols("x y t m lambda")
variables = (x, y, t, m)
weights_to_check = (
    (0, 0, 1, 1),  # dual degree
    (1, 1, 0, 0),  # source degree
    (0, 1, 0, 0),  # the rotating-cone chart
    (1, 2, 3, 5),  # generic positive weight
)

for total_degree in (4, 6):
    for exponents in weak_compositions(total_degree, 4):
        monomial = sp.prod(
            variable**exponent
            for variable, exponent in zip(variables, exponents, strict=True)
        )
        hessian = sp.hessian(monomial, variables)
        for weights in weights_to_check:
            diagonal = sp.diag(*(lam**weight for weight in weights))
            substitution = {
                variable: lam**weight * variable
                for variable, weight in zip(variables, weights, strict=True)
            }
            weighted_degree = sum(
                exponent * weight
                for exponent, weight in zip(exponents, weights, strict=True)
            )
            transformed = diagonal * hessian.subs(substitution) * diagonal
            assert all(
                sp.expand(entry) == 0
                for entry in transformed - lam**weighted_degree * hessian
            )


# 3. On the dual-linear locus
# psi=t*F(x,y)+m*G(x,y)+H(x,y), so the Hessian has a cotangent block and
# determinant Jac(F,G)^2, independently of Hess(H).
Fx, Fy, Gx, Gy = sp.symbols("Fx Fy Gx Gy")
uxx, uxy, uyy = sp.symbols("uxx uxy uyy")
dual_linear_hessian = sp.Matrix(
    [
        [0, 0, Fx, Fy],
        [0, 0, Gx, Gy],
        [Fx, Gx, uxx, uxy],
        [Fy, Gy, uxy, uyy],
    ]
)
plane_jacobian = Fx * Gy - Fy * Gx
assert sp.expand(dual_linear_hessian.det() - plane_jacobian**2) == 0


# 4. Successive-cone recursion.  Work modulo epsilon^(N+1).  Given a unit
# a and arbitrary b, recursively solve b=a*r.  The determinant identity
# a*c-b^2=0 then forces c=a*r^2.  This is the coefficient-by-coefficient
# Schur synchronization used in the proof.
order = 8
epsilon = sp.symbols("epsilon")
a_coefficients = sp.symbols(f"a0:{order + 1}")
b_coefficients = sp.symbols(f"b0:{order + 1}")
r_coefficients: list[sp.Expr] = []

for degree in range(order + 1):
    convolution = sum(
        a_coefficients[index] * r_coefficients[degree - index]
        for index in range(1, degree + 1)
    )
    r_coefficients.append(
        sp.cancel((b_coefficients[degree] - convolution) / a_coefficients[0])
    )

a_series = sum(
    coefficient * epsilon**degree
    for degree, coefficient in enumerate(a_coefficients)
)
b_series = sum(
    coefficient * epsilon**degree
    for degree, coefficient in enumerate(b_coefficients)
)
r_series = sum(
    coefficient * epsilon**degree
    for degree, coefficient in enumerate(r_coefficients)
)

def truncate(expression: sp.Expr) -> sp.Expr:
    return sp.series(expression, epsilon, 0, order + 1).removeO().expand()


assert truncate(b_series - a_series * r_series) == 0
c_series = truncate(a_series * r_series**2)
assert truncate(a_series * c_series - b_series**2) == 0

# Conversely the exact Schur identity identifies the determinant defect.
c = sp.symbols("c")
a = sp.symbols("a", nonzero=True)
b = sp.symbols("b")
r = b / a
schur_matrix = sp.Matrix([[a, b], [b, c]])
cone_matrix = a * sp.Matrix([[1, r], [r, r**2]])
assert sp.simplify(schur_matrix.det() - a * (c - a * r**2)) == 0
assert cone_matrix * sp.Matrix([-r, 1]) == sp.zeros(2, 1)


# 5. The near-miss is not an exception to synchronization.  Under the
# source weight wt(y)=1, wt(x)=wt(t)=wt(m)=0 its dual Hessian is the
# synchronized rotating cone with slope epsilon*y/x.
rotating_quartic = (x * t + y * m) ** 2
dual_hessian = sp.hessian(rotating_quartic, (t, m))
weighted_dual_hessian = dual_hessian.subs({y: epsilon * y})
rotating_slope = epsilon * y / x
expected_rotating_cone = 2 * x**2 * sp.Matrix(
    [
        [1, rotating_slope],
        [rotating_slope, rotating_slope**2],
    ]
)
assert all(
    sp.cancel(entry) == 0
    for entry in weighted_dual_hessian - expected_rotating_cone
)
assert sp.expand(weighted_dual_hessian.det()) == 0
assert weighted_dual_hessian * sp.Matrix(
    [-epsilon * y, x]
) == sp.zeros(2, 1)


print("quartic bigrading:", quartic_ledger)
print("sextic bigrading:", sextic_ledger)
print("dual-linear dimensions: quartic=13 sextic=19")
print("weighted Hessian covariance: 119 monomials x 4 weights")
print("dual-linear cotangent determinant: Jac(F,G)^2")
print(f"successive Schur-cone recursion: exact through order {order}")
print("rotating cone: slope epsilon*y/x, kernel (-epsilon*y,x)")
print("all HC4 source/dual bigrading checks passed")
