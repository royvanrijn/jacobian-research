#!/usr/bin/env python3
"""Exact fixed-r=5 reduction for the cancellation contact resultant.

This checker constructs the bounded-degree endpoint equations without forming
the degree-5*m polynomials in w.  It uses the sparse quartic--binomial-quintic
resultant and subresultant formulas to obtain the primitive endpoint eliminant
H_5(m,y) and the rational reconstruction z=-B/A.

The final uniform separation certificate is developed below this algebraic
core.  No floating-point computation is used for the asserted identities or
the exact small-m gcd checks.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


m, y, z = sp.symbols("m y z")
D = 1 - y


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


# At a common K_{m,5},L_{m,5} root, M_5=0 and z != 0.  Dividing the
# triangular identity by z gives a quartic E_5; M_5 gives the binomial
# quintic F_5.  Writing the coefficients directly preserves their sparsity.
a = (
    D**5
    + 5 * y * D**3 * endpoint_tail(1)
    - 10 * y * D**2 * endpoint_tail(2)
    + 10 * y * D * endpoint_tail(3)
    - 5 * y * endpoint_tail(4)
)
b = -5 * beta(1) * D**3
c = 10 * beta(2) * D**2
d = -10 * beta(3) * D
e = 5 * beta(4)
p = -y * endpoint_tail(5)
q = beta(5)

E5 = a * z**4 + b * z**3 + c * z**2 + d * z + e
F5 = p * z**5 + q

assert sp.degree(E5, z) == 4
assert sp.degree(F5, z) == 5


# Clear parameter denominators once.  Scaling either endpoint equation by a
# nonzero function of m changes neither its roots nor the reconstructed z.
parameter_field = QQ.frac_field(m)
_, integral_E5 = sp.Poly(E5, z, y, domain=parameter_field).clear_denoms(
    convert=True
)
_, integral_F5 = sp.Poly(F5, z, y, domain=parameter_field).clear_denoms(
    convert=True
)

aa, bb, cc, dd, ee = sp.Poly(integral_E5.as_expr(), z).all_coeffs()
pp = sp.Poly(integral_F5.as_expr(), z).coeff_monomial(z**5)
qq = sp.Poly(integral_F5.as_expr(), z).coeff_monomial(1)


# Universal sparse resultant for a quartic and p*z^5+q.  It is also the
# constant member of their subresultant recurrence.  Keeping this 26-term
# norm formula explicit avoids a prohibitively expensive generic PRS over
# QQ(m,y).
def sparse_resultant(
    a0: object,
    b0: object,
    c0: object,
    d0: object,
    e0: object,
    p0: object,
    q0: object,
) -> object:
    return (
        a0**5 * q0**4
        + 5 * a0**3 * b0 * e0 * p0 * q0**3
        + 5 * a0**3 * c0 * d0 * p0 * q0**3
        - 5 * a0**2 * b0**2 * d0 * p0 * q0**3
        - 5 * a0**2 * b0 * c0**2 * p0 * q0**3
        + 5 * a0**2 * c0 * e0**2 * p0**2 * q0**2
        + 5 * a0**2 * d0**2 * e0 * p0**2 * q0**2
        + 5 * a0 * b0**3 * c0 * p0 * q0**3
        + 5 * a0 * b0**2 * e0**2 * p0**2 * q0**2
        - 5 * a0 * b0 * c0 * d0 * e0 * p0**2 * q0**2
        - 5 * a0 * b0 * d0**3 * p0**2 * q0**2
        - 5 * a0 * c0**3 * e0 * p0**2 * q0**2
        + 5 * a0 * c0**2 * d0**2 * p0**2 * q0**2
        + 5 * a0 * d0 * e0**3 * p0**3 * q0
        - b0**5 * p0 * q0**3
        - 5 * b0**3 * d0 * e0 * p0**2 * q0**2
        + 5 * b0**2 * c0**2 * e0 * p0**2 * q0**2
        + 5 * b0**2 * c0 * d0**2 * p0**2 * q0**2
        - 5 * b0 * c0**3 * d0 * p0**2 * q0**2
        + 5 * b0 * c0 * e0**3 * p0**3 * q0
        - 5 * b0 * d0**2 * e0**2 * p0**3 * q0
        + c0**5 * p0**2 * q0**2
        - 5 * c0**2 * d0 * e0**2 * p0**3 * q0
        + 5 * c0 * d0**3 * e0 * p0**3 * q0
        - d0**5 * p0**3 * q0
        + e0**5 * p0**4
    )


# Verify the compact formula once over the universal coefficient ring.
A0, B0, C0, D0, E0, P0, Q0 = sp.symbols("A0 B0 C0 D0 E0 P0 Q0")
assert sp.expand(
    sp.resultant(
        A0 * z**4 + B0 * z**3 + C0 * z**2 + D0 * z + E0,
        P0 * z**5 + Q0,
        z,
    )
    - sparse_resultant(A0, B0, C0, D0, E0, P0, Q0)
) == 0


# Polynomial-ring arithmetic is much faster here than expanding a SymPy
# expression with 49 parameter degrees.  The endpoint factor is exactly
# (y-1)^5 and the remaining primitive eliminant has degree 20 in y.
coefficient_ring, ring_m, ring_y = ring("m,y", QQ)
ring_coefficients = [
    coefficient_ring.from_expr(value) for value in (aa, bb, cc, dd, ee, pp, qq)
]
endpoint_resultant5 = sparse_resultant(*ring_coefficients)
for _ in range(5):
    endpoint_resultant5 = endpoint_resultant5.exquo(ring_y - 1)
assert endpoint_resultant5.rem(ring_y - 1) != 0

eliminant_with_content = sp.Poly(
    endpoint_resultant5.as_expr(), y, domain=QQ.poly_ring(m)
)
eliminant_content, H5 = eliminant_with_content.primitive()
assert sp.factor(eliminant_content - (m + 1)) == 0
assert H5.degree() == 20
assert sp.degree(H5.as_expr(), m) == 49
assert H5.eval(1) != 0


# The penultimate subresultant is A_5(m,y) z+B_5(m,y).  These are the two
# coefficient formulas supplied by the same sparse PRS.
linear_A5 = (
    -aa**3 * cc * pp * qq**2
    + aa**2 * bb**2 * pp * qq**2
    - 2 * aa**2 * dd * ee * pp**2 * qq
    + aa * bb * cc * ee * pp**2 * qq
    + 3 * aa * bb * dd**2 * pp**2 * qq
    - 2 * aa * cc**2 * dd * pp**2 * qq
    - aa * ee**3 * pp**3
    + bb**3 * ee * pp**2 * qq
    - 2 * bb**2 * cc * dd * pp**2 * qq
    + bb * cc**3 * pp**2 * qq
    + 2 * bb * dd * ee**2 * pp**3
    + cc**2 * ee**2 * pp**3
    - 3 * cc * dd**2 * ee * pp**3
    + dd**4 * pp**3
)
linear_B5 = (
    aa**3 * dd * pp * qq**2
    - 2 * aa**2 * bb * cc * pp * qq**2
    + aa**2 * ee**2 * pp**2 * qq
    + aa * bb**3 * pp * qq**2
    - aa * bb * dd * ee * pp**2 * qq
    - 3 * aa * cc**2 * ee * pp**2 * qq
    + 2 * aa * cc * dd**2 * pp**2 * qq
    + 2 * bb**2 * cc * ee * pp**2 * qq
    + bb**2 * dd**2 * pp**2 * qq
    - 3 * bb * cc**2 * dd * pp**2 * qq
    + bb * ee**3 * pp**3
    + cc**4 * pp**2 * qq
    - 2 * cc * dd * ee**2 * pp**3
    + dd**3 * ee * pp**3
)

universal_subresultants = sp.subresultants(
    A0 * z**4 + B0 * z**3 + C0 * z**2 + D0 * z + E0,
    P0 * z**5 + Q0,
    z,
)
assert [sp.degree(member, z) for member in universal_subresultants] == [
    5,
    4,
    3,
    2,
    1,
    0,
]
assert sp.expand(
    universal_subresultants[-2].subs(
        {A0: aa, B0: bb, C0: cc, D0: dd, E0: ee, P0: pp, Q0: qq}
    )
    - (linear_A5 * z + linear_B5)
) == 0


# Exact finite front end.  Directly imposing z=y^m leaves only the excluded
# clearing-denominator root y=1.  The threshold is chosen so that the uniform
# tail certificate can use t=1/m <= 1/20.
excluded_endpoint = sp.Poly((y - 1) ** 5, y).monic()
for m_value in range(1, 20):
    substitutions = {m: m_value, z: y**m_value}
    e_value = sp.Poly(integral_E5.as_expr().subs(substitutions), y)
    f_value = sp.Poly(integral_F5.as_expr().subs(substitutions), y)
    assert sp.gcd(e_value, f_value).monic() == excluded_endpoint


print("PASS contact resultant r=5: sparse eliminant has bidegree (49,20)")
print("PASS contact resultant r=5: explicit linear subresultant reconstructs z")
print("PASS contact resultant r=5: exact endpoint gcds for 1<=m<=19")
print("SCOPE: the uniform m>=20 localization and modulus certificate remains")
