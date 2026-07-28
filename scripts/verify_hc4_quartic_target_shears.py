#!/usr/bin/env python3
"""Exclude degree-at-most-four triangular target shears for HC(4).

Use the normalized foundational map F=(P,B,C).  For each target permutation
(U,V,W), consider

    L = W + H(U,V),

where H is arbitrary of positive total degree at most four.  The existing
cubic checker excludes the subfamily of degree at most three.  This checker
proves that the five homogeneous quartic coefficients must vanish in the
P- and B-retained orientations, using exact extreme spatial coefficients of

    K(L) = -grad(L)^T adj(Hess(L)) grad(L).

In the C-retained orientation, [x^8]K=9 remains independent of every
coefficient through degree four.  Thus the result is exact over every
characteristic-zero coefficient field.

The two large sparse invariants are constructed sequentially so they are
not retained in memory at the same time.
"""

from __future__ import annotations

import gc

from sympy.polys.domains import QQ
from sympy.polys.rings import PolyElement, ring


parameter_names = (
    "a10",
    "a01",
    "a20",
    "a11",
    "a02",
    "a30",
    "a21",
    "a12",
    "a03",
    "a40",
    "a31",
    "a22",
    "a13",
    "a04",
)
polynomial_ring, *generators = ring(
    "x,y,z," + ",".join(parameter_names),
    QQ,
)
x, y, z, *parameters = generators

t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
P = t * q
B = y + 3 * x * q
C = x * (5 - 3 * t) - x**3 * z
mapping = (P, B, C)


def shear_monomials(
    first: PolyElement,
    second: PolyElement,
) -> tuple[PolyElement, ...]:
    """Positive monomials of total degree at most four in two variables."""

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
        first**4,
        first**3 * second,
        first**2 * second**2,
        first * second**3,
        second**4,
    )


def bordered_invariant(
    remain: int,
    zero_source: tuple[int, ...] = (),
) -> PolyElement:
    other = [index for index in range(3) if index != remain]
    first = mapping[other[0]]
    second = mapping[other[1]]
    retained = mapping[remain]
    polynomial = retained + sum(
        coefficient_value * monomial
        for coefficient_value, monomial in zip(
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
    for source_index in zero_source:
        source_variable = generators[source_index]
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
    gx, gy, gz = gradient
    return -(
        gx**2 * (h22 * h33 - h23**2)
        + 2 * gx * gy * (h13 * h23 - h12 * h33)
        + 2 * gx * gz * (h12 * h23 - h13 * h22)
        + gy**2 * (h11 * h33 - h13**2)
        + 2 * gy * gz * (h12 * h13 - h11 * h23)
        + gz**2 * (h11 * h22 - h12**2)
    )


def spatial_coefficient(
    polynomial: PolyElement,
    exponents: tuple[int, int, int],
    zero_parameters: tuple[int, ...] = (),
) -> PolyElement:
    result = polynomial_ring.zero
    for monomial, coefficient_value in polynomial.terms():
        if (
            monomial[:3] == exponents
            and all(monomial[3 + index] == 0 for index in zero_parameters)
        ):
            result += polynomial_ring.term_new(
                (0, 0, 0) + monomial[3:],
                coefficient_value,
            )
    return result


# Parameter indices 9,...,13 are a40,a31,a22,a13,a04.
#
# P+H(B,C): the five displayed extreme coefficients successively force
# a40=a04=a31=a13=a22=0.  Lower-degree shear terms are present in this full
# 14-parameter calculation and do not contaminate any displayed coefficient.
K_P = bordered_invariant(0)
assert spatial_coefficient(K_P, (46, 30, 14)) == (
    -396_718_580_736 * parameters[9] ** 4
)
assert spatial_coefficient(K_P, (44, 0, 12), (9,)) == (
    2_304 * parameters[13] ** 4
)
assert spatial_coefficient(K_P, (46, 22, 14), (9, 13)) == (
    -3_367_210_176 * parameters[10] ** 4
)
assert spatial_coefficient(K_P, (46, 6, 14), (9, 13, 10)) == (
    -139_968 * parameters[12] ** 4
)
assert spatial_coefficient(K_P, (46, 14, 14), (9, 13, 10, 12)) == (
    -25_194_240 * parameters[11] ** 4
)
del K_P
gc.collect()


# B+H(P,C): the same endpoint-to-middle elimination removes its quartic
# layer.  Again the calculation contains every lower-degree coefficient.
K_B = bordered_invariant(1)
assert spatial_coefficient(K_B, (46, 46, 14)) == (
    -16_128 * parameters[9] ** 4
)
assert spatial_coefficient(K_B, (44, 0, 12), (9,)) == (
    2_304 * parameters[13] ** 4
)
assert spatial_coefficient(K_B, (46, 34, 14), (9, 13)) == (
    -10_800 * parameters[10] ** 4
)
assert spatial_coefficient(K_B, (46, 10, 14), (9, 13, 10)) == (
    -2_736 * parameters[12] ** 4
)
assert spatial_coefficient(K_B, (46, 22, 14), (9, 13, 10, 12)) == (
    -6_336 * parameters[11] ** 4
)
del K_B
gc.collect()


# C+H(P,B): on the x-axis, every quartic contribution disappears before
# the coefficient x^8.  The unavoidable scalar from the retained C remains.
K_C_axis = bordered_invariant(2, (1, 2))
assert spatial_coefficient(K_C_axis, (8, 0, 0)) == 9


print("PASS: P+H(B,C) has no nonzero homogeneous quartic layer")
print("PASS: B+H(P,C) has no nonzero homogeneous quartic layer")
print("PASS: C+H(P,B) retains the unavoidable bordered coefficient 9")
print("PASS: no quartic-or-lower triangular target shear has K(L)=0")
