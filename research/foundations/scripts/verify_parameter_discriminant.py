#!/usr/bin/env python3
"""Exact regressions for the closed cancellation construction parameter-discriminant formula."""
from __future__ import annotations

import sympy as sp

from master_cancellation import (
    even_square_discriminant_family,
    parameter_discriminant,
    parameter_discriminant_is_square,
    parameter_polynomial,
)


q, x, u = sp.symbols("q x u")
square_pairs: set[tuple[int, int]] = set()
checked = 0

for m in range(1, 31):
    for r in range(1, 31):
        n = m * r
        if n > 30:
            continue

        polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
        reciprocal = sum(sp.binomial(j + r, r) * x**j for j in range(n + 1))
        geometric = sum(x**j for j in range((m + 1) * r + 1))
        normalized_derivative = sp.diff(geometric, x, r) / sp.factorial(r)
        assert sp.expand(reciprocal - normalized_derivative) == 0
        assert sp.cancel(
            x**n * polynomial.as_expr().subs(q, 1 - 1 / x)
            - (-1) ** n * reciprocal
        ) == 0

        leading = sp.binomial(n + r, r)
        differential_identity = (
            (1 - x) * sp.diff(reciprocal, x)
            - (r + 1) * reciprocal
            + (n + r + 1) * leading * x**n
        )
        assert sp.expand(differential_identity) == 0

        actual = sp.Integer(sp.discriminant(polynomial.as_expr(), q))
        expected = parameter_discriminant(m, r)
        assert actual == expected

        actual_square = actual > 0 and sp.integer_nthroot(int(actual), 2)[1]
        criterion_square = parameter_discriminant_is_square(m, r)
        assert bool(actual_square) == criterion_square
        if criterion_square:
            square_pairs.add((m, r))
        checked += 1

assert square_pairs == {
    (1, 1),
    (4, 3),
    (2, 8),
    (16, 1),
    (17, 1),
    (1, 24),
    (12, 2),
}
print(f"PASS: closed parameter-discriminant formula for all {checked} pairs with mr<=30")
print("PASS: parity square criterion and complete square list through degree 30")

# Complete parametrization of the even-degree square locus.  If
# r+1=d*a^2 with d squarefree, then the criterion is equivalent to
# n+r+1=d*b^2.  The explicit family takes b=a+2*r*k.
for r in range(1, 51):
    factors = sp.factorint(r + 1)
    squarefree = sp.prod(
        prime for prime, exponent in factors.items() if exponent % 2
    )
    a, exact = sp.integer_nthroot((r + 1) // squarefree, 2)
    assert exact

    for k in range(1, 6):
        m = even_square_discriminant_family(r, k)
        b = a + 2 * r * k
        assert m == squarefree * (b * b - a * a) // r
        assert m * r % 4 == 0
        assert m * r + r + 1 == squarefree * b * b
        assert parameter_discriminant_is_square(m, r)

    for m in range(1, 51):
        n = m * r
        if n % 4:
            continue
        total = n + r + 1
        represented = total % squarefree == 0 and sp.integer_nthroot(
            total // squarefree, 2
        )[1]
        assert parameter_discriminant_is_square(m, r) == bool(represented)

print("PASS: complete even-degree parametrization and infinite family for each r<=50")

# Complete odd-degree square locus in the r=1 row.  It is exactly
# m=2*a^2-1 with a odd.
for a in range(1, 32, 2):
    m = 2 * a * a - 1
    assert m % 4 == 1
    assert parameter_discriminant_is_square(m, 1)

for m in range(1, 2001, 2):
    represented = False
    if m % 4 == 1:
        a, exact = sp.integer_nthroot((m + 1) // 2, 2)
        represented = bool(exact and a % 2 == 1)
    assert parameter_discriminant_is_square(m, 1) == represented

print("PASS: geometric-derivative identity and complete odd r=1 square family")

# For fixed odd r>=3, the odd-degree square condition is an integral-point
# problem on a squarefree hyperelliptic curve of degree r.  Clear the
# denominator r! to obtain an integral model without changing its roots.
for r in range(3, 22, 2):
    curve_rhs = (r + 1) * sp.factorial(r) * sp.prod(
        r * (u + 1) - index for index in range(r)
    )
    curve = sp.Poly(curve_rhs, u, domain=sp.ZZ)
    assert curve.degree() == r
    assert sp.gcd(curve, curve.diff()).degree() == 0

print("PASS: odd fixed-r square locus has squarefree hyperelliptic model")
