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


def integer_partitions(total: int, largest: int | None = None):
    """Yield partitions of ``total`` as nonincreasing tuples."""
    if total < 0:
        raise ValueError("total must be nonnegative")
    if total == 0:
        yield ()
        return
    if largest is None or largest > total:
        largest = total
    for part in range(largest, 0, -1):
        for tail in integer_partitions(total - part, part):
            yield (part,) + tail


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


def _pencil_mod_prime(primitive, variable, slope: int, intercept: int, prime: int):
    """Construct ``H-sW+t`` over ``F_prime``, including rational inputs."""
    if not sp.isprime(prime):
        raise ValueError("the diagnostic routines currently expect a prime field")
    H = sp.Poly(primitive, variable)
    coefficients = _coefficients_mod(H.as_expr(), variable, prime)
    if not coefficients or max(coefficients) != H.degree():
        raise ValueError(f"inverse degree drops modulo {prime}")
    coefficients[1] = (coefficients.get(1, 0) - int(slope)) % prime
    coefficients[0] = (coefficients.get(0, 0) + int(intercept)) % prime
    expression = sum(
        coefficient * variable**exponent
        for exponent, coefficient in coefficients.items()
    )
    return sp.Poly(expression, variable, modulus=prime)


def factorization_type_mod_prime(
    primitive, variable, slope: int, intercept: int, prime: int
) -> tuple[int, ...] | None:
    """Return the squarefree factor-degree partition of ``H-sW+t``.

    ``None`` records a discriminant point (at least one repeated factor).
    """
    polynomial = _pencil_mod_prime(primitive, variable, slope, intercept, prime)
    factors = polynomial.factor_list()[1]
    if any(multiplicity != 1 for _, multiplicity in factors):
        return None
    return tuple(sorted((factor.degree() for factor, _ in factors), reverse=True))


def pencil_factorization_histogram(primitive, variable, prime: int):
    """Count every squarefree cycle type and the ramified locus over ``F_p``."""
    histogram = Counter()
    for slope in range(prime):
        for intercept in range(prime):
            cycle_type = factorization_type_mod_prime(
                primitive, variable, slope, intercept, prime
            )
            histogram[cycle_type] += 1
    return dict(histogram)


def balanced_residue(value: int, prime: int) -> int:
    """Return the least-absolute integer representative modulo ``prime``."""
    residue = int(value) % prime
    return residue - prime if residue > prime // 2 else residue


def find_pencil_factorization_witness(
    primitive, variable, prime: int, cycle_type
):
    """Find a deterministic finite-field witness for a prescribed cycle type.

    The returned integer lift preserves the factorization fingerprint modulo
    ``prime``.  It does not in general preserve factorization over Q.
    """
    H = sp.Poly(primitive, variable)
    requested = tuple(sorted((int(part) for part in cycle_type), reverse=True))
    if not requested or any(part < 1 for part in requested):
        raise ValueError("cycle_type must be a nonempty partition")
    if sum(requested) != H.degree():
        raise ValueError(f"cycle_type must partition the inverse degree {H.degree()}")

    for slope in range(prime):
        for intercept in range(prime):
            if factorization_type_mod_prime(
                H.as_expr(), variable, slope, intercept, prime
            ) != requested:
                continue
            polynomial = _pencil_mod_prime(
                H.as_expr(), variable, slope, intercept, prime
            )
            unit, factor_data = polynomial.factor_list()
            factors = tuple(
                sp.Poly(factor, variable, modulus=prime).as_expr()
                for factor, _ in factor_data
            )
            lifted_slope = balanced_residue(slope, prime)
            lifted_intercept = balanced_residue(intercept, prime)
            return {
                "prime": prime,
                "cycle_type": requested,
                "slope": slope,
                "intercept": intercept,
                "factorization_unit": int(unit) % prime,
                "factors": factors,
                "lifted_slope": lifted_slope,
                "lifted_intercept": lifted_intercept,
                "lifted_polynomial": sp.expand(
                    H.as_expr() - lifted_slope * variable + lifted_intercept
                ),
            }
    return None


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
