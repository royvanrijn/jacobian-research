"""Symmetric-group fixed-point laws and finite-field pencil diagnostics."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import math

import sympy as sp


def derangement_count(n: int) -> int:
    """Return the number of fixed-point-free permutations of ``n`` letters."""
    if n < 0:
        raise ValueError("n must be nonnegative")
    if n == 0:
        return 1
    if n == 1:
        return 0
    previous, current = 1, 0
    for degree in range(2, n + 1):
        previous, current = current, (degree - 1) * (current + previous)
    return current


def fixed_point_count(n: int, fixed: int) -> int:
    """Count permutations in S_n with exactly ``fixed`` fixed points."""
    if not 0 <= fixed <= n:
        return 0
    return math.comb(n, fixed) * derangement_count(n - fixed)


def fixed_point_distribution(n: int) -> dict[int, Fraction]:
    """Return the exact fixed-point probability law in the natural S_n action."""
    if n < 1:
        raise ValueError("n must be positive")
    order = math.factorial(n)
    return {
        fixed: Fraction(fixed_point_count(n, fixed), order)
        for fixed in range(n + 1)
        if fixed_point_count(n, fixed)
    }


def _coefficients_mod(polynomial, variable, prime):
    coefficients = {}
    for (exponent,), coefficient in sp.Poly(polynomial, variable).terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        if denominator % prime == 0:
            raise ValueError(f"coefficient denominator is not a unit modulo {prime}")
        reduced = numerator * pow(denominator, -1, prime) % prime
        if reduced:
            coefficients[exponent] = reduced
    return coefficients


def _evaluate_mod(coefficients, value, prime):
    return sum(
        coefficient * pow(value, exponent, prime)
        for exponent, coefficient in coefficients.items()
    ) % prime


def pencil_simple_root_histogram(primitive, variable, prime: int) -> dict[int, int]:
    """Enumerate simple F_p-roots of H(W)-sW+t over all ``(s,t)``."""
    if not sp.isprime(prime):
        raise ValueError("the diagnostic enumerator currently expects a prime field")
    H = sp.Poly(primitive, variable)
    coefficients = _coefficients_mod(H.as_expr(), variable, prime)
    derivative = _coefficients_mod(sp.diff(H.as_expr(), variable), variable, prime)
    if not coefficients or max(coefficients) != H.degree():
        raise ValueError(f"inverse degree drops modulo {prime}")

    histogram = Counter()
    for slope in range(prime):
        for intercept in range(prime):
            simple_roots = 0
            for root in range(prime):
                value = (
                    _evaluate_mod(coefficients, root, prime)
                    - slope * root
                    + intercept
                ) % prime
                derivative_value = (
                    _evaluate_mod(derivative, root, prime) - slope
                ) % prime
                if value == 0 and derivative_value != 0:
                    simple_roots += 1
            histogram[simple_roots] += 1
    return dict(sorted(histogram.items()))
