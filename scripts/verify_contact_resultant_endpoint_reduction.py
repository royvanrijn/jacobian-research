#!/usr/bin/env python3
"""Exact endpoint-moment reduction for the cancellation contact resultant.

This checker separates proved identities from the remaining conjecture.
It verifies the all-(m,r) endpoint formulas on a bounded exact grid and proves
symbolically, with m left as a parameter, that Res(K_{m,r},L_{m,r}) is
nonzero for the three complete columns r=1, r=2, and r=3.  The r=3 proof is
a Schur--Cohn separation certificate for the fixed-degree endpoint eliminant.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.boundary import (  # noqa: E402
    cancellation_branch_polynomial,
    cancellation_critical_coefficient,
)


w, y, z, u = sp.symbols("w y z u")


def endpoint_moment(m: int, k: int) -> sp.Expr:
    """M_k after y=1-w and z=y^m, with z retained independently."""

    beta = sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))
    tail = y * z**k * sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )
    return sp.factor((beta - tail) / (1 - y) ** (k + 1))


# The endpoint formula and the triangular binomial identity are general; this
# finite loop is an exact regression, not their all-parameter proof.
for m_value in range(1, 7):
    for r_value in range(1, 7):
        substitution = {y: 1 - w, z: (1 - w) ** m_value}
        moments = [endpoint_moment(m_value, k).subs(substitution) for k in range(r_value + 1)]
        K = cancellation_branch_polynomial(m_value, r_value, w)
        L = cancellation_critical_coefficient(m_value, r_value, w)
        assert sp.factor(moments[r_value] - K) == 0
        triangular = sum(
            (-1) ** k
            * sp.binomial(r_value, k)
            * (1 - w) ** (m_value * (r_value - k))
            * moments[k]
            for k in range(r_value + 1)
        )
        assert sp.factor(triangular - L) == 0


# Uniform r=1 proof.  Here L_{m,1}+K_{m,1}=(1-w)^m and
# K_{m,1}(1)=Beta(2,m+1)=1/((m+1)(m+2)).
m = sp.symbols("m", integer=True, positive=True)
K1_at_one = sp.factorial(1) / sp.prod(m + j for j in range(1, 3))
assert sp.factor(K1_at_one - 1 / ((m + 1) * (m + 2))) == 0
for m_value in range(1, 12):
    K1 = cancellation_branch_polynomial(m_value, 1, w)
    L1 = cancellation_critical_coefficient(m_value, 1, w)
    assert sp.factor(K1 + L1 - (1 - w) ** m_value) == 0
    resultant = sp.factor(sp.resultant(K1, L1, w))
    assert resultant == sp.Rational(1, ((m_value + 1) * (m_value + 2)) ** m_value)


# Uniform r=2 proof.  Put M_k=K_{m,k}, y=1-w and z=y^m.  At a common
# K_{m,2},L_{m,2} root, z is nonzero and L_2=z^2-2zM_1+M_2, hence 2M_1=z.
# The first equation is z*A=2.  Combining it with M_2=0 eliminates z and
# gives the displayed factored polynomial in y.
M1 = endpoint_moment(m, 1)
M2 = endpoint_moment(m, 2)
A = m * (m + 1) * y**2 - 2 * m * (m + 2) * y + (m + 1) * (m + 2)
assert sp.factor(
    (2 * M1 - z) * (m + 1) * (m + 2) * (1 - y) ** 2 + (z * A - 2)
) == 0

beta2 = 2 / ((2 * m + 1) * (2 * m + 2) * (2 * m + 3))
tail2 = 1 / (2 * m + 1) - 2 * y / (2 * m + 2) + y**2 / (2 * m + 3)
assert sp.factor(M2 - (beta2 - y * z**2 * tail2) / (1 - y) ** 3) == 0

# Substitute z=2/A into M_2=0 and clear the nonzero denominators.
eliminant = sp.factor(
    2 * y * tail2 * (2 * m + 1) * (2 * m + 2) * (2 * m + 3) - A**2
)
expected = -(
    (m + 1) ** 2
    * (y - 1) ** 3
    * (m**2 * y - (m + 2) ** 2)
)
assert sp.factor(eliminant - expected) == 0

# y=1 is w=0, where K_{m,2}(0)=1/3.  At the other candidate,
# y=(m+2)^2/m^2>1, both y^m>1 and A>2, contradicting z*A=2.
y_candidate = (m + 2) ** 2 / m**2
A_candidate = sp.factor(A.subs(y, y_candidate))
assert sp.factor(
    A_candidate - 2 * (m + 2) * (5 * m**2 + 10 * m + 4) / m**3
) == 0

for m_value in range(1, 12):
    K2 = cancellation_branch_polynomial(m_value, 2, w)
    L2 = cancellation_critical_coefficient(m_value, 2, w)
    assert sp.resultant(K2, L2, w) != 0


# Uniform r=3 proof.  At a common root M_3=0 and L_3=0.  Since z is
# nonzero, the triangular identity reduces L_3=0 to
#
#     z^2 - 3*z*M_1 + 3*M_2 = 0.
#
# Clearing endpoint denominators gives a quadratic E=a*z^2+b*z+c and the
# binomial cubic F=d*z^3+e.  Their compact resultant factors off the excluded
# endpoint (y-1)^3 and leaves a genuine sextic H_m(y).
M3 = endpoint_moment(m, 3)
D = 1 - y
beta1 = 1 / ((m + 1) * (m + 2))
beta2 = 2 / ((2 * m + 1) * (2 * m + 2) * (2 * m + 3))
beta3 = 6 / sp.prod(3 * m + j for j in range(1, 5))


def endpoint_tail(k: int) -> sp.Expr:
    return sp.factor(
        sum(
            (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
            for j in range(k + 1)
        )
    )


T1 = endpoint_tail(1)
T2 = endpoint_tail(2)
T3 = endpoint_tail(3)
a = sp.factor(D**3 + 3 * y * D * T1 - 3 * y * T2)
b = sp.factor(-3 * beta1 * D)
c = sp.factor(3 * beta2)
d = sp.factor(-y * T3)
e = sp.factor(beta3)

E3 = sp.cancel((z**2 - 3 * z * M1 + 3 * M2) * D**3)
F3 = sp.cancel(M3 * D**4)
assert sp.factor(E3 - (a * z**2 + b * z + c)) == 0
assert sp.factor(F3 - (d * z**3 + e)) == 0

# Universal identity for Res_z(a*z^2+b*z+c, d*z^3+e).
aa, bb, cc, dd, ee = sp.symbols("aa bb cc dd ee")
compact_resultant = (
    aa**3 * ee**2
    + 3 * aa * bb * cc * dd * ee
    - bb**3 * dd * ee
    + cc**3 * dd**2
)
assert sp.factor(
    sp.resultant(aa * z**2 + bb * z + cc, dd * z**3 + ee, z)
    - compact_resultant
) == 0

endpoint_resultant3 = sp.factor(
    compact_resultant.subs({aa: a, bb: b, cc: c, dd: d, ee: e})
)
eliminant3 = sp.cancel(endpoint_resultant3 / (y - 1) ** 3)
eliminant3_numerator, eliminant3_denominator = sp.together(
    eliminant3
).as_numer_denom()
assert sp.degree(eliminant3_denominator, y) == 0
H3 = sp.Poly(eliminant3_numerator, y, domain=sp.ZZ[m])
assert H3.degree() == 6
assert H3.eval(1) != 0
assert sp.cancel(endpoint_resultant3 - (y - 1) ** 3 * eliminant3) == 0

# For every r, the endpoint transform is the negative-binomial section
#
#   K_{m,r}(1-y) = beta_r * sum_{j=0}^{mr} binomial(j+r,r) y^j.
#
# Its consecutive coefficient ratios are j/(j+r), so Enestrom--Kakeya puts
# every K_{m,r} root in |y| <= m/(m+1).  The finite loop is an exact
# regression of the general coefficient identity in the r=3 column.
for m_value in range(1, 9):
    K3 = cancellation_branch_polynomial(m_value, 3, w)
    beta3_value = sp.factorial(3) / sp.prod(
        3 * m_value + j for j in range(1, 5)
    )
    negative_binomial_section = sum(
        sp.binomial(j + 3, 3) * y**j for j in range(3 * m_value + 1)
    )
    assert sp.factor(
        K3.subs(w, 1 - y) - beta3_value * negative_binomial_section
    ) == 0

# To prove that every H_m root lies strictly outside that disk, put
# rho=m/(m+1) and Q_m(u)=u^6 H_m(rho/u).  The roots of Q are rho/y.  The
# Schur--Cohn Hermitian matrix is positive definite exactly when all Q roots
# lie in the open unit disk.  Sylvester's criterion reduces this to its six
# leading principal minors.  After a harmless common integral scaling, every
# coefficient of every minor is strictly positive in m.
rho = m / (m + 1)
Q3_numerator = sp.together(
    u**6 * H3.as_expr().subs(y, rho / u)
).as_numer_denom()[0]
Q3 = sp.Poly(Q3_numerator, u, domain=sp.ZZ[m]).primitive()[1]
assert Q3.degree() == 6

coefficients = Q3.all_coeffs()
size = Q3.degree()
forward = sp.zeros(size)
reverse = sp.zeros(size)
for row in range(size):
    for column in range(row + 1):
        forward[row, column] = coefficients[row - column]
        reverse[row, column] = coefficients[size - (row - column)]

schur_cohn = forward * forward.T - reverse * reverse.T
expected_minor_degrees = [29, 56, 81, 104, 125, 144]
minor_degrees: list[int] = []
for principal_size in range(1, size + 1):
    principal_minor = sp.factor(
        schur_cohn[:principal_size, :principal_size].det(method="domain-ge")
    )
    minor_polynomial = sp.Poly(principal_minor, m, domain=sp.ZZ)
    minor_degrees.append(minor_polynomial.degree())
    assert all(coefficient > 0 for coefficient in minor_polynomial.all_coeffs())

assert minor_degrees == expected_minor_degrees


print("PASS contact endpoint reduction: general moment identities on 1<=m,r<=6")
print("PASS contact resultant: uniform nonvanishing for every m and r=1")
print("PASS contact resultant: uniform nonvanishing for every m and r=2")
print("PASS contact resultant: Schur--Cohn separation for every m and r=3")
print("SCOPE: this checker covers r<=3; the separate r=4 checker completes the next column")
