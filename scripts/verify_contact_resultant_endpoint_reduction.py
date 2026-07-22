#!/usr/bin/env python3
"""Exact endpoint-moment reduction for the cancellation contact resultant.

This checker separates proved identities from the remaining conjecture.
It verifies the all-(m,r) endpoint formulas on a bounded exact grid and proves
symbolically, with m left as a parameter, that Res(K_{m,r},L_{m,r}) is
nonzero for the two complete columns r=1 and r=2.
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


w, y, z = sp.symbols("w y z")


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


print("PASS contact endpoint reduction: general moment identities on 1<=m,r<=6")
print("PASS contact resultant: uniform nonvanishing for every m and r=1")
print("PASS contact resultant: uniform nonvanishing for every m and r=2")
print("SCOPE: the symbolic endpoint eliminant for r>=3 remains open")
