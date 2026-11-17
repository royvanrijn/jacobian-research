#!/usr/bin/env python3
"""Exclude all degree-at-most-three triangular target shears for HC(4).

Use the foundational determinant-minus-two quadratic-gauge map F=(P,B,C).
For each target permutation (U,V,W), consider the polynomial target
coordinate

    L = W + H(U,V),

where H is an arbitrary polynomial of total degree at most three and its
constant term is omitted as irrelevant.  A double-dual Schur descent would
leave a four-variable potential

    psi(x,y,z,s) = s*L(F(x,y,z)) + h(x,y,z).

The coefficient of s^2 in det Hess(psi) is the bordered-Hessian invariant

    K(L) = -grad(L)^T adj(Hess(L)) grad(L).

This checker proves K(L) cannot vanish in any of the three triangular
orientations.  It uses a sparse polynomial ring and a short triangular list
of exact spatial coefficients.  No full four-variable determinant is needed
because every candidate fails this necessary first gate.

Scope: constant target permutations followed by one triangular shear of
degree at most three, on the foundational normalized cubic gauge.  General
nontriangular target automorphisms, preliminary nonlinear changes of the
other two target coordinates, degree at least four, and nonlinear
source--dual symplectic transformations are not tested.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import PolyElement, ring


# Every admissible cubic gauge is linearly source--target equivalent to the
# foundational representative used below.  If alpha=g3/g1 and beta=g2/g1,
# use z'=alpha*z and (P',B',C')=(alpha*P,B-2*beta*P,C).
source_x, source_y, source_z = sp.symbols("source_x source_y source_z")
alpha, beta = sp.symbols("alpha beta", nonzero=True)
source_t = 1 + source_x * source_y
source_q = (
    source_t**2 * source_z
    + source_y**2 * (1 + 3 * source_t) / alpha
)
source_P = source_t * source_q
source_B = (
    source_y
    + 3 * alpha * source_x * source_q
    + 2 * beta * source_t * source_q
)
source_C = (
    source_x * (5 - 3 * source_t)
    - alpha * source_x**3 * source_z
)
normalized_z = alpha * source_z
normalized_q = (
    source_t**2 * normalized_z
    + source_y**2 * (1 + 3 * source_t)
)
assert sp.expand(normalized_q - alpha * source_q) == 0
assert sp.expand(source_t * normalized_q - alpha * source_P) == 0
assert sp.expand(
    source_y + 3 * source_x * normalized_q
    - (source_B - 2 * beta * source_P)
) == 0
assert sp.expand(
    source_x * (5 - 3 * source_t)
    - source_x**3 * normalized_z
    - source_C
) == 0


polynomial_ring, *ring_generators = ring(
    "x,y,z," + ",".join(f"a{index}" for index in range(9)),
    QQ,
)
x, y, z, *parameters = ring_generators

t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
P = t * q
B = y + 3 * x * q
C = x * (5 - 3 * t) - x**3 * z
mapping = (P, B, C)


def shear_monomials(first: PolyElement, second: PolyElement) -> tuple[PolyElement, ...]:
    """Monomials of positive total degree at most three in two variables."""

    return (
        first,
        second,
        first**2,
        first * second,
        second**2,
        first**3,
        first**2 * second,
        first * second**2,
        second**3,
    )


def bordered_invariant(
    remain: int,
    zero_source: tuple[int, ...] = (),
) -> PolyElement:
    """Return K(W+H(U,V)) for one target orientation."""

    other = [index for index in range(3) if index != remain]
    first = mapping[other[0]]
    second = mapping[other[1]]
    retained = mapping[remain]
    polynomial = retained + sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            parameters,
            shear_monomials(first, second),
            strict=True,
        )
    )

    gradient = [polynomial.diff(index) for index in range(3)]
    hessian = [
        [gradient[row].diff(column) for column in range(3)]
        for row in range(3)
    ]
    if zero_source:
        for source_index in zero_source:
            source_variable = ring_generators[source_index]
            gradient = [
                entry.compose(source_variable, polynomial_ring.zero)
                for entry in gradient
            ]
            hessian = [
                [
                    entry.compose(source_variable, polynomial_ring.zero)
                    for entry in row
                ]
                for row in hessian
            ]
    h11, h12, h13 = hessian[0]
    _, h22, h23 = hessian[1]
    _, _, h33 = hessian[2]

    adj11 = h22 * h33 - h23**2
    adj12 = h13 * h23 - h12 * h33
    adj13 = h12 * h23 - h13 * h22
    adj22 = h11 * h33 - h13**2
    adj23 = h12 * h13 - h11 * h23
    adj33 = h11 * h22 - h12**2

    return -(
        gradient[0] ** 2 * adj11
        + 2 * gradient[0] * gradient[1] * adj12
        + 2 * gradient[0] * gradient[2] * adj13
        + gradient[1] ** 2 * adj22
        + 2 * gradient[1] * gradient[2] * adj23
        + gradient[2] ** 2 * adj33
    )


def spatial_coefficient(
    polynomial: PolyElement,
    exponents: tuple[int, int, int],
    zero_parameters: tuple[int, ...] = (),
) -> PolyElement:
    """Extract one spatial coefficient after optional parameter zeros."""

    result = polynomial_ring.zero
    for monomial, coefficient in polynomial.terms():
        if (
            monomial[:3] == exponents
            and all(monomial[3 + index] == 0 for index in zero_parameters)
        ):
            result += polynomial_ring.term_new(
                (0, 0, 0) + monomial[3:],
                coefficient,
            )
    return result


# Orientation P+H(B,C).  Four high layers first kill the cubic coefficients:
#
#   a5=B^3, a6=B^2 C, a7=B C^2, a8=C^3.
#
# Afterward two quadratic-layer coefficients simultaneously force a2=0 and
# 12*a2=1.  Thus no lower coefficients need to be considered.
K_P = bordered_invariant(0)
assert spatial_coefficient(K_P, (34, 22, 10)) == (
    -1_549_681_956 * parameters[5] ** 4
)
assert spatial_coefficient(K_P, (32, 0, 8)) == (
    729 * parameters[8] ** 4
)
assert spatial_coefficient(K_P, (34, 14, 10), (5, 8)) == (
    -11_337_408 * parameters[6] ** 4
)
assert spatial_coefficient(K_P, (34, 6, 10), (5, 8)) == (
    -61_236 * parameters[7] ** 4
)
assert spatial_coefficient(K_P, (8, 0, 6), (5, 6, 7, 8)) == (
    -1_259_712 * parameters[2] ** 4
)
assert spatial_coefficient(K_P, (0, 0, 2), (5, 6, 7, 8)) == (
    -9 * (12 * parameters[2] - 1)
)


# Orientation B+H(P,C).  The same four-step top-layer argument removes the
# cubic part.  The remaining displayed coefficients successively force
#
#   a2=a4=a3=a0=a1=0,
#
# after which the constant spatial coefficient is 9.
K_B = bordered_invariant(1)
assert spatial_coefficient(K_B, (34, 34, 10)) == (
    -5_103 * parameters[5] ** 4
)
assert spatial_coefficient(K_B, (32, 0, 8)) == (
    729 * parameters[8] ** 4
)
assert spatial_coefficient(K_B, (34, 22, 10), (5, 8)) == (
    -2_916 * parameters[6] ** 4
)
assert spatial_coefficient(K_B, (34, 10, 10), (5, 8)) == (
    -1_215 * parameters[7] ** 4
)
assert spatial_coefficient(K_B, (0, 0, 6), (5, 6, 7, 8)) == (
    144 * parameters[2] ** 4
)
assert spatial_coefficient(K_B, (20, 0, 4), (5, 6, 7, 8)) == (
    144 * parameters[4] ** 4
)
assert spatial_coefficient(K_B, (12, 0, 6), (2, 4, 5, 6, 7, 8)) == (
    -216 * parameters[3] ** 4
)
assert spatial_coefficient(K_B, (0, 0, 2), (2, 4, 5, 6, 7, 8)) == (
    9 * parameters[0] ** 4
)
assert spatial_coefficient(K_B, (8, 0, 0), (2, 4, 5, 6, 7, 8)) == (
    9 * parameters[1] ** 4
)
assert spatial_coefficient(K_B, (0, 0, 0), tuple(range(9))) == 9


# Orientation C+H(P,B) is excluded before any parameter comparison: the
# x^8 coefficient is the scalar 9 for the full cubic shear.
K_C = bordered_invariant(2, (1, 2))
assert spatial_coefficient(K_C, (8, 0, 0)) == 9


print("PASS: P+H(B,C) has no bordered-flat cubic-or-lower triangular shear")
print("PASS: B+H(P,C) has no bordered-flat cubic-or-lower triangular shear")
print("PASS: C+H(P,B) retains the unavoidable bordered coefficient 9")
