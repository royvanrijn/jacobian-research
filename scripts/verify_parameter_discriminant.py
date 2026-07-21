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


q, x = sp.symbols("q x")
square_pairs: set[tuple[int, int]] = set()
checked = 0

for m in range(1, 31):
    for r in range(1, 31):
        n = m * r
        if n > 30:
            continue

        polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.ZZ)
        reciprocal = sum(sp.binomial(j + r, r) * x**j for j in range(n + 1))
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
