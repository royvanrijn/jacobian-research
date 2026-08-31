"""Exact local-discriminant fingerprints and rational parameter recovery.

The input family is an integral one-parameter Weierstrass model together with
its discriminant polynomial.  At every selected prime this module solves the
homogenized discriminant congruence on both charts of ``P^1(Q_p)``, compresses
complete root sets to maximal balls, combines one ball per prime by CRT, and
Gauss-reduces the resulting homogeneous congruence lattice.

Repeated discriminant valuations generally determine balls rather than a
unique residue.  Consequently the bounded lattice enumeration is an
experiment.  An exact match is reported only after comparing the candidate
and target j-invariants by integer arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from typing import Any, Iterable, Sequence

from .crt_lattice import (
    crt,
    gauss_reduce_linear_congruence_lattice,
    hensel_lift_roots,
    p_adic_valuation,
)


Polynomial = tuple[int, ...]
Vector = tuple[int, int]


def _trim(coefficients: Sequence[int]) -> Polynomial:
    result = list(coefficients)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def polynomial_add(*polynomials: Sequence[int]) -> Polynomial:
    degree = max((len(polynomial) for polynomial in polynomials), default=1)
    return _trim(
        [
            sum(
                polynomial[index] if index < len(polynomial) else 0
                for polynomial in polynomials
            )
            for index in range(degree)
        ]
    )


def polynomial_scale(polynomial: Sequence[int], scalar: int) -> Polynomial:
    return _trim([scalar * coefficient for coefficient in polynomial])


def polynomial_multiply(*polynomials: Sequence[int]) -> Polynomial:
    result: Polynomial = (1,)
    for polynomial in polynomials:
        next_result = [0] * (len(result) + len(polynomial) - 1)
        for left_degree, left in enumerate(result):
            for right_degree, right in enumerate(polynomial):
                next_result[left_degree + right_degree] += left * right
        result = _trim(next_result)
    return result


def weierstrass_invariant_polynomials(
    model: dict[str, Sequence[int]],
) -> dict[str, Polynomial]:
    """Derive ``c4``, ``c6`` and ``Delta`` from polynomial a-invariants."""

    required = ("a1", "a2", "a3", "a4", "a6")
    missing = [name for name in required if name not in model]
    if missing:
        raise ValueError(f"family model omits {missing}")
    a1, a2, a3, a4, a6 = (
        _trim(tuple(int(value) for value in model[name])) for name in required
    )
    b2 = polynomial_add(polynomial_multiply(a1, a1), polynomial_scale(a2, 4))
    b4 = polynomial_add(polynomial_scale(a4, 2), polynomial_multiply(a1, a3))
    b6 = polynomial_add(polynomial_multiply(a3, a3), polynomial_scale(a6, 4))
    b8 = polynomial_add(
        polynomial_multiply(a1, a1, a6),
        polynomial_scale(polynomial_multiply(a2, a6), 4),
        polynomial_scale(polynomial_multiply(a1, a3, a4), -1),
        polynomial_multiply(a2, a3, a3),
        polynomial_scale(polynomial_multiply(a4, a4), -1),
    )
    c4 = polynomial_add(
        polynomial_multiply(b2, b2), polynomial_scale(b4, -24)
    )
    c6 = polynomial_add(
        polynomial_scale(polynomial_multiply(b2, b2, b2), -1),
        polynomial_scale(polynomial_multiply(b2, b4), 36),
        polynomial_scale(b6, -216),
    )
    discriminant = polynomial_add(
        polynomial_scale(polynomial_multiply(b2, b2, b8), -1),
        polynomial_scale(polynomial_multiply(b4, b4, b4), -8),
        polynomial_scale(polynomial_multiply(b6, b6), -27),
        polynomial_scale(polynomial_multiply(b2, b4, b6), 9),
    )
    return {"c4": c4, "c6": c6, "discriminant": discriminant}


def weierstrass_invariants(coefficients: Sequence[int]) -> dict[str, int]:
    if len(coefficients) != 5:
        raise ValueError("five integral Weierstrass coefficients are required")
    constant_model = {
        name: (int(value),)
        for name, value in zip(("a1", "a2", "a3", "a4", "a6"), coefficients)
    }
    result = weierstrass_invariant_polynomials(constant_model)
    return {name: values[0] for name, values in result.items()}


def homogeneous_value(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> int:
    """Evaluate ``b^degree*f(a/b)`` without rational arithmetic."""

    coefficients = _trim(coefficients)
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def small_prime_valuations(value: int, maximum_prime: int) -> tuple[tuple[int, int], ...]:
    """Return all nonzero valuations at primes through ``maximum_prime``."""

    if value == 0:
        raise ValueError("zero has no finite prime factorization")
    if maximum_prime < 2:
        return ()
    primes: list[int] = []
    for candidate in range(2, maximum_prime + 1):
        if all(candidate % prime for prime in primes if prime <= isqrt(candidate)):
            primes.append(candidate)
    return tuple(
        (prime, p_adic_valuation(value, prime))
        for prime in primes
        if value % prime == 0
    )


def factor_over_known_primes(
    value: int, primes: Iterable[int]
) -> tuple[tuple[tuple[int, int], ...], int]:
    """Divide by a declared prime support and return factors plus cofactor."""

    if value == 0:
        raise ValueError("zero has no finite factorization")
    remainder = abs(value)
    factors = []
    for prime in sorted(set(int(prime) for prime in primes)):
        if prime < 2:
            raise ValueError("declared factors must be at least two")
        if remainder % prime:
            continue
        exponent = 0
        while remainder % prime == 0:
            remainder //= prime
            exponent += 1
        factors.append((prime, exponent))
    return tuple(factors), remainder


def select_repeated_prime_constraints(
    value: int,
    *,
    maximum_prime: int,
    minimum_valuation: int = 2,
    excluded_primes: Iterable[int] = (),
) -> tuple[tuple[int, int], ...]:
    """Select unusually repeated small-prime powers by a declared threshold."""

    if minimum_valuation < 1:
        raise ValueError("minimum_valuation must be positive")
    excluded = {int(prime) for prime in excluded_primes}
    return tuple(
        (prime, exponent)
        for prime, exponent in small_prime_valuations(value, maximum_prime)
        if exponent >= minimum_valuation and prime not in excluded
    )


@dataclass(frozen=True, order=True)
class ProjectiveRootBall:
    """A local projective class forcing a target discriminant valuation."""

    prime: int
    target_exponent: int
    chart: str
    exponent: int
    residue: int

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    @property
    def row(self) -> tuple[int, int]:
        if self.chart == "affine":
            return 1, -self.residue
        if self.chart == "infinity":
            return -self.residue, 1
        raise ValueError(f"unknown projective chart {self.chart!r}")

    @property
    def label(self) -> str:
        variable = "a/b" if self.chart == "affine" else "b/a"
        return (
            f"p={self.prime} {variable}={self.residue} "
            f"(mod {self.prime}^{self.exponent})"
        )

    def matches(self, numerator: int, denominator: int) -> bool:
        if self.chart == "affine":
            return (
                denominator % self.prime != 0
                and (numerator - self.residue * denominator) % self.modulus == 0
            )
        return (
            numerator % self.prime != 0
            and (denominator - self.residue * numerator) % self.modulus == 0
        )


def _compress_roots(
    roots: Iterable[int], prime: int, target_exponent: int
) -> tuple[tuple[int, int], ...]:
    balls: set[tuple[int, int]] = {
        (target_exponent, root) for root in roots
    }
    for child_exponent in range(target_exponent, 0, -1):
        child_modulus = prime**child_exponent
        parent_modulus = child_modulus // prime
        parents = {
            residue % parent_modulus
            for exponent, residue in balls
            if exponent == child_exponent
        }
        for parent in parents:
            children = {
                (child_exponent, parent + digit * parent_modulus)
                for digit in range(prime)
            }
            if children <= balls:
                balls.difference_update(children)
                balls.add((child_exponent - 1, parent))
    return tuple(sorted(balls))


def discover_projective_root_balls(
    coefficients: Sequence[int],
    prime: int,
    exponent: int,
    *,
    maximum_roots: int = 100_000,
) -> tuple[tuple[ProjectiveRootBall, ...], dict[str, Any]]:
    """Solve a homogeneous polynomial congruence on both projective charts."""

    coefficients = _trim(tuple(int(value) for value in coefficients))
    affine_roots = tuple(
        hensel_lift_roots(
            coefficients,
            prime,
            exponent,
            maximum_roots=maximum_roots,
        )
    )
    reciprocal_roots = hensel_lift_roots(
        tuple(reversed(coefficients)),
        prime,
        exponent,
        maximum_roots=maximum_roots,
    )
    # The affine chart already contains all b-units.  Retaining only b/a
    # divisible by p makes the two charts a disjoint partition.
    infinity_roots = tuple(root for root in reciprocal_roots if root % prime == 0)
    affine_balls = tuple(
        ProjectiveRootBall(prime, exponent, "affine", ball_exponent, residue)
        for ball_exponent, residue in _compress_roots(
            affine_roots, prime, exponent
        )
    )
    infinity_balls = tuple(
        ProjectiveRootBall(prime, exponent, "infinity", ball_exponent, residue)
        for ball_exponent, residue in _compress_roots(
            infinity_roots, prime, exponent
        )
    )
    balls = tuple(sorted((*affine_balls, *infinity_balls)))
    profile = {
        "prime": prime,
        "target_valuation": exponent,
        "affine_target_root_count": len(affine_roots),
        "infinity_target_root_count": len(infinity_roots),
        "maximal_balls": [
            {
                "chart": ball.chart,
                "residue": ball.residue,
                "exponent": ball.exponent,
                "modulus": ball.modulus,
                "leaf_count": prime ** (exponent - ball.exponent),
            }
            for ball in balls
        ],
    }
    return balls, profile


def select_projective_ball(
    balls: Sequence[ProjectiveRootBall], branch: dict[str, Any]
) -> ProjectiveRootBall:
    """Select and verify a declared branch from a complete local root set."""

    if not balls:
        raise ValueError("cannot select from an empty projective root set")
    chart = str(branch["chart"])
    exponent = int(branch["exponent"])
    if exponent < 0 or exponent > balls[0].target_exponent:
        raise ValueError("branch precision must lie between zero and the target")
    residue = int(branch["residue"]) % (balls[0].prime**exponent)
    candidates = [
        ball
        for ball in balls
        if ball.chart == chart
        and exponent >= ball.exponent
        and residue % ball.modulus == ball.residue
    ]
    if not candidates:
        raise ValueError(
            f"declared branch {chart}:{residue} mod p^{exponent} is not a root"
        )
    prime = balls[0].prime
    target_exponent = balls[0].target_exponent
    return ProjectiveRootBall(
        prime, target_exponent, chart, exponent, residue
    )


def combine_projective_rows(
    choices: Sequence[ProjectiveRootBall],
) -> tuple[int, int, int]:
    coefficient_a, coefficient_b, modulus = 0, 0, 1
    seen_primes: set[int] = set()
    for choice in choices:
        if choice.prime in seen_primes:
            raise ValueError("use at most one independent constraint per prime")
        seen_primes.add(choice.prime)
        row_a, row_b = choice.row
        coefficient_a, next_modulus = crt(
            ((coefficient_a, modulus), (row_a, choice.modulus))
        )
        coefficient_b, check_modulus = crt(
            ((coefficient_b, modulus), (row_b, choice.modulus))
        )
        if next_modulus != check_modulus:
            raise AssertionError("CRT row moduli diverged")
        modulus = next_modulus
    if gcd(gcd(coefficient_a, coefficient_b), modulus) != 1:
        raise AssertionError("the combined projective row is not primitive")
    return coefficient_a, coefficient_b, modulus


def _normalize_vector(vector: Vector, modulus: int) -> Vector | None:
    numerator, denominator = vector
    if denominator == 0:
        return None
    common = gcd(abs(numerator), abs(denominator))
    if gcd(common, modulus) != 1:
        return None
    numerator //= common
    denominator //= common
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if gcd(abs(numerator), denominator) != 1:
        raise AssertionError("rational normalization failed")
    return numerator, denominator


def enumerate_projective_lattice(
    basis: tuple[Vector, Vector],
    modulus: int,
    *,
    coefficient_radius: int,
    height_cap: int | None = None,
) -> Iterable[Vector]:
    """Enumerate a declared reduced-basis box of primitive rational pairs."""

    if coefficient_radius < 0:
        raise ValueError("coefficient_radius must be nonnegative")
    first, second = basis
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            pair = _normalize_vector(
                (
                    left * first[0] + right * second[0],
                    left * first[1] + right * second[1],
                ),
                modulus,
            )
            if pair is None:
                continue
            if height_cap is not None and max(abs(pair[0]), pair[1]) > height_cap:
                continue
            yield pair


def exact_j_match(
    numerator: int,
    denominator: int,
    *,
    family_c4: Sequence[int],
    family_discriminant: Sequence[int],
    target_c4: int,
    target_discriminant: int,
) -> bool:
    """Compare j-invariants using a homogeneous integer identity."""

    c4 = homogeneous_value(family_c4, numerator, denominator)
    discriminant = homogeneous_value(
        family_discriminant, numerator, denominator
    )
    if discriminant == 0:
        return False
    denominator_exponent = (
        len(family_discriminant) - 1 - 3 * (len(family_c4) - 1)
    )
    left = c4**3 * target_discriminant
    right = target_c4**3 * discriminant
    if denominator_exponent >= 0:
        left *= denominator**denominator_exponent
    else:
        right *= denominator ** (-denominator_exponent)
    return left == right


def recover_parameters(
    discriminant_coefficients: Sequence[int],
    constraint_groups: Sequence[Sequence[ProjectiveRootBall]],
    *,
    coefficient_radius: int,
    height_cap: int | None,
    weights: Vector = (1, 1),
    family_c4: Sequence[int] | None = None,
    target_c4: int | None = None,
    target_discriminant: int | None = None,
    candidate_limit: int = 20,
) -> dict[str, Any]:
    """Enumerate a complete declared CRT/Gauss box and verify exact matches."""

    if candidate_limit < 0:
        raise ValueError("candidate_limit must be nonnegative")
    if not constraint_groups or any(not group for group in constraint_groups):
        raise ValueError("every selected prime needs at least one root ball")
    candidates: dict[Vector, tuple[ProjectiveRootBall, ...]] = {}
    class_count = 0
    vector_count = 0
    for choices in product(*constraint_groups):
        class_count += 1
        coefficient_a, coefficient_b, modulus = combine_projective_rows(choices)
        basis = gauss_reduce_linear_congruence_lattice(
            coefficient_a,
            coefficient_b,
            modulus,
            weights=weights,
        )
        for pair in enumerate_projective_lattice(
            basis,
            modulus,
            coefficient_radius=coefficient_radius,
            height_cap=height_cap,
        ):
            vector_count += 1
            if not all(choice.matches(*pair) for choice in choices):
                raise AssertionError("a lattice vector left its projective balls")
            candidates.setdefault(pair, tuple(choices))

    nonsingular_candidates = {
        pair: choices
        for pair, choices in candidates.items()
        if homogeneous_value(discriminant_coefficients, *pair) != 0
    }
    singular_count = len(candidates) - len(nonsingular_candidates)
    candidates = nonsingular_candidates
    ordered = sorted(
        candidates,
        key=lambda pair: (
            max(abs(pair[0]), pair[1]),
            pair[0] * pair[0] + pair[1] * pair[1],
            abs(pair[0]),
            pair[1],
            pair[0],
        ),
    )
    exact_matches: list[Vector] = []
    if family_c4 is not None:
        if target_c4 is None or target_discriminant is None:
            raise ValueError("target c4 and discriminant are required for j matching")
        exact_matches = [
            pair
            for pair in ordered
            if exact_j_match(
                *pair,
                family_c4=family_c4,
                family_discriminant=discriminant_coefficients,
                target_c4=target_c4,
                target_discriminant=target_discriminant,
            )
        ]

    def record(pair: Vector) -> dict[str, Any]:
        numerator, denominator = pair
        homogeneous = homogeneous_value(
            discriminant_coefficients, numerator, denominator
        )
        choices = candidates[pair]
        return {
            "numerator": numerator,
            "denominator": denominator,
            "parameter": str(Fraction(numerator, denominator)),
            "height": max(abs(numerator), denominator),
            "homogeneous_discriminant": str(homogeneous),
            "forced_valuations": {
                str(choice.prime): p_adic_valuation(homogeneous, choice.prime)
                for choice in choices
            },
            "local_branches": [choice.label for choice in choices],
        }

    return {
        "scope": (
            "complete enumeration of the declared projective branches and "
            "reduced-basis coefficient box"
        ),
        "projective_class_count": class_count,
        "normalized_vectors_in_declared_boxes": vector_count,
        "unique_bounded_parameters": len(ordered),
        "singular_parameters_removed": singular_count,
        "shortest_candidates": [record(pair) for pair in ordered[:candidate_limit]],
        "exact_j_matches": [record(pair) for pair in exact_matches],
    }
