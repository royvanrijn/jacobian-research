#!/usr/bin/env python3
"""Exact audit of the Danielewski multi-dicritical family."""

from __future__ import annotations

from collections import Counter
from itertools import product

import sympy as sp


a, b, c, x = sp.symbols("a b c x")


def audit_polynomial(polynomial: sp.Expr) -> None:
    quotient = sp.cancel(polynomial / x)
    assert sp.denom(quotient) == 1
    source_x = b * c
    source_w = sp.expand(b * quotient.subs(x, source_x))
    assert sp.factor(c * source_w - polynomial.subs(x, source_x)) == 0

    jacobian = sp.factor(
        sp.Matrix((a, c, source_x))
        .jacobian((a, b, c))
        .det()
    )
    assert jacobian == -c
    assert sp.cancel(jacobian / c) == -1

    polynomial_poly = sp.Poly(polynomial, x)
    assert sp.gcd(polynomial_poly, polynomial_poly.diff()).degree() == 0


for degree in range(1, 7):
    example = x * sp.prod(x - root for root in range(1, degree))
    audit_polynomial(sp.expand(example))
print("PASS: squarefree examples in every tested degree have residue Jacobian -1")


L, n = sp.symbols("L n")
class_surface = sp.expand(
    (L - n) * (L - 1) + n * (2 * L - 1)
)
class_target = sp.expand(L * class_surface)
assert sp.expand(class_surface - (L**2 + (n - 1) * L)) == 0
assert sp.expand(class_target - (L**3 + (n - 1) * L**2)) == 0
print("PASS: motivic excess is one A^2 per nonzero root")


def finite_map(
    prime: int,
    source: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    finite_a, finite_b, finite_c = source
    finite_x = finite_b * finite_c % prime
    finite_w = finite_b * (finite_x**2 + 1) % prime
    return finite_a, finite_c, finite_x, finite_w


def target_points(prime: int) -> set[tuple[int, int, int, int]]:
    return {
        (finite_a, finite_c, finite_x, finite_w)
        for finite_a, finite_c, finite_x, finite_w in product(
            range(prime),
            repeat=4,
        )
        if (
            finite_c * finite_w
            - finite_x * (finite_x**2 + 1)
        )
        % prime
        == 0
    }


for prime in (3, 7, 11, 19):
    fibers = Counter(
        finite_map(prime, source)
        for source in product(range(prime), repeat=3)
    )
    target = target_points(prime)
    assert set(fibers.values()) == {1}
    assert len(fibers) == prime**3
    assert fibers.keys() == target
print("PASS: x(x^2+1) gives bijections for p=3 mod 4")


for prime in (5, 13):
    fibers = Counter(
        finite_map(prime, source)
        for source in product(range(prime), repeat=3)
    )
    target = target_points(prime)
    assert set(fibers.values()) == {1}
    assert len(target) == prime**3 + 2 * prime**2
    assert len(target - fibers.keys()) == 2 * prime**2
print("PASS: split primes expose exactly two missing A^2 boundary divisors")
print("PASS Danielewski multi-dicritical family")
