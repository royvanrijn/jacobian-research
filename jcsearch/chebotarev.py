"""Symmetric-group fixed-point laws and finite-field pencil diagnostics."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import math

import sympy as sp


def _rational_integer_data(expressions) -> tuple[int, list[sp.Expr]]:
    """Clear rational denominators and return their lcm and integral forms."""
    denominator = 1
    expanded = []
    for expression in expressions:
        expression = sp.expand(sp.sympify(expression))
        expanded.append(expression)
        coefficients = (sp.Poly(expression).coeffs()
                        if expression.free_symbols else [expression])
        for coefficient in coefficients:
            denominator = math.lcm(
                denominator, abs(int(sp.denom(sp.cancel(coefficient))))
            )
    return denominator, [sp.expand(denominator * expression) for expression in expanded]


def rational_good_reduction_certificate(
    primitive, variable, *, c=sp.Integer(1), b0=sp.Integer(1), a0=sp.Integer(0)
):
    """Return an explicit rational-seed bad-prime certificate.

    Every prime not dividing ``bad_integer`` preserves the degree, the
    weighted Jacobian unit, the squarefree multiplicity profile of ``H``, and
    the characteristic hypotheses used in the uniform S_n monodromy proof.
    """
    H = sp.Poly(sp.expand(primitive), variable, domain=sp.QQ)
    n = H.degree()
    denominator, integral = _rational_integer_data((H.as_expr(), c, b0, a0))
    integral_H, integral_c, integral_b0, _ = integral

    # The squarefree decomposition is defined over Q.  Its discriminants and
    # pairwise resultants are precisely the collision certificate for the
    # primitive-root multiplicity profile.
    _, squarefree = sp.sqf_list(H.as_expr(), variable)
    boundary_resultant = sp.Integer(1)
    factors = [sp.Poly(factor, variable, domain=sp.QQ) for factor, _ in squarefree]
    for factor in factors:
        if factor.degree() > 1:
            boundary_resultant *= sp.discriminant(factor.as_expr(), variable)
    for left in range(len(factors)):
        for right in range(left + 1, len(factors)):
            boundary_resultant *= sp.resultant(
                factors[left].as_expr(), factors[right].as_expr(), variable
            )
    boundary_resultant = sp.cancel(boundary_resultant)
    boundary_num, boundary_den = map(int, sp.fraction(boundary_resultant))

    leading = int(sp.Poly(integral_H, variable).LC())
    c_integer = int(integral_c)
    b_integer = int(integral_b0)
    bad_integer = abs(
        denominator
        * math.factorial(n)
        * leading
        * c_integer
        * b_integer
        * boundary_num
        * boundary_den
    )
    if bad_integer == 0:
        raise ValueError("the seed/model data do not define a nonzero certificate")

    slope, intercept = sp.symbols("s t")
    pencil = H.as_expr() - slope * variable + intercept
    discriminant = sp.primitive(sp.together(sp.discriminant(pencil, variable)))[1]
    return {
        "degree": n,
        "denominator": denominator,
        "boundary_resultant": sp.factor(boundary_resultant),
        "bad_integer": bad_integer,
        "discriminant": sp.factor(discriminant),
        "discriminant_degree": sp.Poly(discriminant, slope, intercept).total_degree(),
    }


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


def cycle_type_centralizer_size(parts) -> int:
    """Return ``z_lambda`` for a partition giving a symmetric-group cycle type."""
    partition = tuple(int(part) for part in parts)
    if not partition or any(part < 1 for part in partition):
        raise ValueError("cycle-type parts must be positive")
    multiplicities = Counter(partition)
    return math.prod(
        length**multiplicity * math.factorial(multiplicity)
        for length, multiplicity in multiplicities.items()
    )


def cycle_type_probability(parts) -> Fraction:
    """Return the proportion of permutations having the prescribed cycle type."""
    return Fraction(1, cycle_type_centralizer_size(parts))


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
