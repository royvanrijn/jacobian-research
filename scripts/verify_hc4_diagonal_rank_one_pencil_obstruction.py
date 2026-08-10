#!/usr/bin/env python3
"""Verify the rank-one pencil obstruction on the diagonal HC4 quintic top.

Use constant-kernel coordinates ``(x,y,z,t)`` and normalize

    Hess(h5)|_(x,y,z) = diag(x^3,y^3,z^3).

The quintic Schur face then has

    D_t h4 = (a*x^3+b*y^3+c*z^3)/3,
    D_t^2 h3 = a^2*x+b^2*y+c^2*z.

For arbitrary compatible lower pieces, this checker extracts six coefficients
of ``ell^T adj(Hess(psi)) ell``.  They prove that a nonzero Schur channel
forces the projective constant-null-covector scheme to be empty.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def weak_compositions(length: int, degree: int) -> tuple[tuple[int, ...], ...]:
    if length == 1:
        return ((degree,),)
    return tuple(
        (first, *tail)
        for first in range(degree + 1)
        for tail in weak_compositions(length - 1, degree - first)
    )


def generic_homogeneous(
    variables: Sequence[sp.Symbol], degree: int, prefix: str
) -> sp.Expr:
    return sp.Add(
        *(
            sp.Symbol(prefix + "_" + "".join(map(str, exponents)))
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(variables, exponents, strict=True)
            )
            for exponents in weak_compositions(len(variables), degree)
        )
    )


x, y, z, t = variables = sp.symbols("x y z t")
a, b, c = sp.symbols("a b c")
p, q, r = sp.symbols("p q r")

# The diagonal Schur pair solves the degree-ten determinant face exactly.
C = sp.diag(x**3, y**3, z**3)
d = sp.Matrix([a * x**2, b * y**2, c * z**2])
f = a**2 * x + b**2 * y + c**2 * z
Delta = C.det()
assert sp.expand(Delta * f - (d.T * C.adjugate() * d)[0]) == 0

# The degree-nine metric face is tau^2*Delta, so any constant null covector
# has zero t-component.  This face is independent of all lower terms.
tau = sp.symbols("tau")
top_hessian = sp.diag(x**3, y**3, z**3, 0)
top_covector = sp.Matrix([p, q, r, tau])
top_metric = sp.expand(
    (top_covector.T * top_hessian.adjugate() * top_covector)[0]
)
assert top_metric == tau**2 * x**3 * y**3 * z**3

# Retain every lower coefficient compatible with the first Schur face.  Their
# names are deliberately generic: the decisive six metric coefficients below
# must not depend on any of them.
r4 = generic_homogeneous((x, y, z), 4, "r4")
g2 = generic_homogeneous((x, y, z), 2, "g2")
r3 = generic_homogeneous((x, y, z), 3, "r3")
q2 = generic_homogeneous(variables, 2, "q2")

h5 = (x**5 + y**5 + z**5) / 20
s3 = (a * x**3 + b * y**3 + c * z**3) / 3
psi = sp.expand(h5 + t * s3 + r4 + t**2 * f / 2 + t * g2 + r3 + q2)
H = sp.hessian(psi, variables)
active_covector = sp.Matrix([p, q, r, 0])
metric = sp.Poly(
    sp.expand((active_covector.T * H.adjugate(method="domain-ge") * active_covector)[0]),
    *variables,
)

# The degree-seven channel coefficients first impose ap=bq=cr=0.
assert sp.factor(metric.coeff_monomial(x * y**3 * z**3)) == a**2 * p**2
assert sp.factor(metric.coeff_monomial(x**3 * y * z**3)) == b**2 * q**2
assert sp.factor(metric.coeff_monomial(x**3 * y**3 * z)) == c**2 * r**2

# The next immutable t^2 coefficients synchronize the zero-channel directions.
assert sp.expand(
    metric.coeff_monomial(x**3 * t**2) + (b**2 * r - c**2 * q) ** 2
) == 0
assert sp.expand(
    metric.coeff_monomial(y**3 * t**2) + (a**2 * r - c**2 * p) ** 2
) == 0
assert sp.expand(
    metric.coeff_monomial(z**3 * t**2) + (a**2 * q - b**2 * p) ** 2
) == 0

# Over a characteristic-zero field, if (say) a is nonzero, the first channel
# gives p=0 and the last two displayed squares give r=q=0.  The b- and c-
# charts are identical by symmetry.  Hence a nonzero (p,q,r) requires
# a=b=c=0, which is the already-closed aligned branch.

print("PASS: diagonal quintic Schur pair solves the determinant face")
print("PASS: the top metric face forces the null covector to be active")
print("PASS: three degree-seven coefficients impose ap=bq=cr=0")
print("PASS: three immutable t^2 squares kill every zero-channel direction")
print("PASS: a nonaligned diagonal Schur packet has empty rank-one null scheme")
print("SCOPE: higher-rank and nonlinear pencil directions remain open")
