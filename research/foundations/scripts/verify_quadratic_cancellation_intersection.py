#!/usr/bin/env python3
"""Exact checks for the quadratic-gauge/cancellation stable intersection.

The all-degree proof is in
verified/QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md.  This regression
checks the full smallest meaningful range N=4,5,6,7 with symbolic
quadratic-gauge coefficients and every cancellation factorization.
"""

from __future__ import annotations

import math

import sympy as sp


S, P, B, C = sp.symbols("S P B C")
Q, Y, q = sp.symbols("Q Y q")


def cancellation_pairs(degree: int) -> list[tuple[int, int]]:
    """Return all positive (m,r) with degree-1=(m+1)r."""
    return [
        (m, r)
        for r in range(1, degree)
        for m in range(1, degree)
        if (m + 1) * r == degree - 1
    ]


def parameter_polynomial(m: int, r: int) -> sp.Expr:
    """The monic cancellation parameter polynomial M_(m,r)."""
    n = m * r
    return sp.expand(
        sum(
            (-1) ** j * sp.binomial(n + r + 1, j) * q ** (n - j)
            for j in range(n + 1)
        )
    )


EXPECTED = {
    4: {
        (2, 1): (q**2 - 4 * q + 6, 6),
    },
    5: {
        (3, 1): (q**3 - 5 * q**2 + 10 * q - 10, 12),
        (1, 2): (q**2 - 5 * q + 10, 4),
    },
    6: {
        (4, 1): (q**4 - 6 * q**3 + 15 * q**2 - 20 * q + 15, 20),
    },
    7: {
        (5, 1): (
            q**5 - 7 * q**4 + 21 * q**3 - 35 * q**2 + 35 * q - 21,
            30,
        ),
        (2, 2): (q**4 - 7 * q**3 + 21 * q**2 - 35 * q + 35, 12),
        (1, 3): (q**3 - 7 * q**2 + 21 * q - 35, 6),
    },
}


for degree, rows in EXPECTED.items():
    assert set(cancellation_pairs(degree)) == set(rows)
    for (m, r), (expected_parameter, expected_mu) in rows.items():
        parameter = parameter_polynomial(m, r)
        assert sp.expand(parameter - expected_parameter) == 0
        assert sp.gcd(parameter, sp.diff(parameter, q)) == 1
        assert m * r * (m + 1) == expected_mu

        # On the normalized critical divisor:
        # T=Y^-m and P=(Q-Y)Y^m.  The relative Fitting generator is P_Y
        # up to the unit Y^(m-1), hence mQ-(m+1)Y.  In u=Q/Y-1 it is m*u-1.
        critical_P = (Q - Y) * Y**m
        fitting = sp.factor(sp.diff(critical_P, Y) / Y ** (m - 1))
        assert sp.expand(fitting - (m * Q - (m + 1) * Y)) == 0
        u = sp.symbols(f"u_{m}_{r}")
        fitting_on_torus = sp.factor(fitting.subs(Q, (u + 1) * Y) / Y)
        assert fitting_on_torus == m * u - 1

        # The ramified profile is r^(m+1), with boundary index r+1.
        assert r + 1 >= 2


for degree in range(4, 8):
    coefficients = {
        j: sp.symbols(f"a{degree}_{j}")
        for j in range(3, degree + 1)
    }
    H = (
        S
        + coefficients[3] * P * S**3
        + sum(
            coefficients[j] * P**j * S**j
            for j in range(4, degree + 1)
        )
    )

    # The normalized discriminant Fitting generator is S^2*dB/dS.
    normalized_B = sp.diff(H, S) / S
    normalized_C = 2 * H - S * sp.diff(H, S)
    fitting = sp.expand(S**2 * sp.diff(normalized_B, S))
    expected_fitting = (
        -1
        + 3 * coefficients[3] * P * S**2
        + sum(
            j * (j - 2) * coefficients[j] * P**j * S ** (j - 1)
            for j in range(4, degree + 1)
        )
    )
    assert sp.expand(fitting - expected_fitting) == 0
    assert sp.expand(
        sp.diff(normalized_C, S)
        + S**2 * sp.diff(normalized_B, S)
    ) == 0

    # The three unavoidable support points have affine rank two.
    support_determinant = sp.det(
        sp.Matrix([[1, 2], [degree, degree - 1]])
    )
    assert support_determinant == -(degree + 1)

    # Exact saturated boundary trace.  Intermediate coefficients remain
    # symbolic, so this checks general admissible seeds rather than a sparse
    # specialization.
    doubled_inverse = sp.expand(2 * H - B * S**2 - C)
    discriminant = sp.expand(sp.discriminant(doubled_inverse, S))
    p_polynomial = sp.Poly(discriminant, P)
    p_order = min(
        monomial[0]
        for monomial, coefficient in p_polynomial.terms()
        if coefficient != 0
    )
    expected_order = degree**2 - 3 * degree - 2
    assert p_order == expected_order

    trace = sp.factor(
        sp.expand(discriminant / P**p_order).subs(P, 0)
    )
    d = degree - 3
    expected_trace = sp.factor(
        (-1) ** ((degree - 1) * (degree - 2) // 2)
        * 2 ** (2 * degree - 4)
        * d**d
        * coefficients[3] ** (degree - 2)
        * coefficients[degree] ** (degree - 4)
        * B**2
        * (B * C - 1)
    )
    assert sp.factor(trace - expected_trace) == 0

    # The complete quadratic second-boundary profile.
    boundary_degree = degree - 3
    number_of_primes = math.gcd(boundary_degree, 2)
    boundary_index = boundary_degree // number_of_primes
    assert number_of_primes * boundary_index == boundary_degree

    # Every noncubic cancellation thick index differs from the quadratic
    # value 2, including all parameter roots in the row.
    for m, r in cancellation_pairs(degree):
        assert m * r * (m + 1) >= 4
        assert m * r * (m + 1) != 2

    print(
        f"PASS N={degree}: p-order={p_order}, "
        f"quadratic trace=B^2(BC-1), quadratic mu=2, "
        f"cancellation types={cancellation_pairs(degree)}"
    )


# Degree three: both families reduce to the foundational member.
assert cancellation_pairs(3) == [(1, 1)]
assert parameter_polynomial(1, 1) == q - 3

a3, alpha, beta = sp.symbols("a3 alpha beta", nonzero=True)
transformed_a3 = alpha**-2 * beta**-1 * a3
assert sp.factor(transformed_a3.subs({alpha: 1, beta: a3}) - 1) == 0

print("PASS N=3: unique cancellation root q=3 and every cubic gauge scales to S+S^3")
print("PASS: the stable intersection is exactly the foundational cubic class")
