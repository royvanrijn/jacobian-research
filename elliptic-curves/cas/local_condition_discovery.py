#!/usr/bin/env python3
"""Exact discovery of power-efficient local conditions for Fermigier's H(T).

This module constructs the invariant polynomials of the normalized
Fermigier--Mestre family without a symbolic-algebra dependency, enumerates
roots of ``H`` modulo declared prime powers, compresses them to maximal
p-adic balls, and classifies the resulting conditions at primes ``p >= 5``.

The clean classifications are uniform on the whole ball.  A clean
multiplicative ball has a unit minimal ``c4`` and a uniform split symbol.  A
clean additive ball also carries an exact obstruction to another integral
minimalizing scale.  Balls for which those statements cannot be proved at
the current precision are explicitly marked mixed and may be refined.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from math import log
from typing import Any, Sequence

from ek_k3 import legendre_symbol
from fermigier_mestre import (
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
)
from multiple_root_lifting import (
    RootBall,
    affine_variable_coefficients,
    all_roots_mod_prime_power,
    fixed_divisor_valuation,
)


Polynomial = tuple[int, ...]


def _trim(coefficients: list[int]) -> Polynomial:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def polynomial_add(*polynomials: Sequence[int]) -> Polynomial:
    if not polynomials:
        return (0,)
    degree = max(len(polynomial) for polynomial in polynomials)
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
    answer = (1,)
    for polynomial in polynomials:
        product = [0] * (len(answer) + len(polynomial) - 1)
        for left_degree, left in enumerate(answer):
            for right_degree, right in enumerate(polynomial):
                product[left_degree + right_degree] += left * right
        answer = _trim(product)
    return answer


def polynomial_value(polynomial: Sequence[int], value: int) -> int:
    answer = 0
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


def fermigier_invariant_polynomials() -> dict[str, Polynomial]:
    """Return exact power-basis polynomials for I, J, a4, a6, c4, and c6."""

    quartic_a = (1_149_050, 0, 1)
    quartic_b = (-30 * 68_377_393, 0, -30 * 62)
    quartic_c = (1_195_214_262_641, 0, 1_718_550, 0, -2)
    quartic_d = (-30 * 8_594_794_400_346, 0, -30 * 21_690_305, 0, 30 * 62)
    quartic_e = (
        18_103_855_887_324_900,
        0,
        102_302_344_648,
        0,
        -879_500,
        0,
        1,
    )
    invariant_i = polynomial_add(
        polynomial_scale(polynomial_multiply(quartic_a, quartic_e), 12),
        polynomial_scale(polynomial_multiply(quartic_b, quartic_d), -3),
        polynomial_multiply(quartic_c, quartic_c),
    )
    invariant_j = polynomial_add(
        polynomial_scale(
            polynomial_multiply(quartic_a, quartic_c, quartic_e), 72
        ),
        polynomial_scale(
            polynomial_multiply(quartic_b, quartic_c, quartic_d), 9
        ),
        polynomial_scale(
            polynomial_multiply(quartic_a, quartic_d, quartic_d), -27
        ),
        polynomial_scale(
            polynomial_multiply(quartic_b, quartic_b, quartic_e), -27
        ),
        polynomial_scale(
            polynomial_multiply(quartic_c, quartic_c, quartic_c), -2
        ),
    )
    return {
        "I": invariant_i,
        "J": invariant_j,
        "a4": polynomial_scale(invariant_i, -27),
        "a6": polynomial_scale(invariant_j, -27),
        "c4": polynomial_scale(invariant_i, 1296),
        "c6": polynomial_scale(invariant_j, 23328),
    }


INVARIANT_POLYNOMIALS = fermigier_invariant_polynomials()


def integer_valuation(integer: int, prime: int) -> int:
    if integer == 0:
        return 10**9
    value = abs(integer)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def coefficient_content_valuation(
    polynomial: Sequence[int], prime: int
) -> int:
    return min(
        integer_valuation(coefficient, prime)
        for coefficient in polynomial
        if coefficient
    )


def uniform_exact_valuation(
    polynomial: Sequence[int], prime: int
) -> int | None:
    """Prove a uniform exact valuation using coefficient content and F_p.

    If the content-normalized polynomial has no root in ``F_p``, every
    integer value has exactly the returned valuation.  ``None`` means only
    that this inexpensive uniform proof did not succeed.
    """

    content = coefficient_content_valuation(polynomial, prime)
    divisor = prime**content
    normalized = tuple(coefficient // divisor for coefficient in polynomial)
    if all(polynomial_value(normalized, residue) % prime for residue in range(prime)):
        return content
    return None


@dataclass(frozen=True)
class ClassifiedBall:
    prime: int
    exponent: int
    residue: int
    forced_h_valuation: int
    h_coefficient_content_valuation: int
    presented_model_scaling: int
    forced_minimal_discriminant_valuation: int
    reduction: str
    split_symbols_seen: tuple[bool, ...]
    proof: str

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    @property
    def density(self) -> Fraction:
        return Fraction(1, self.modulus)


def classify_ball(ball: RootBall) -> ClassifiedBall:
    """Classify reduction uniformly on one exact p-adic ball when possible."""

    prime = ball.prime
    modulus = ball.modulus
    transformed = {
        "H": affine_variable_coefficients(
            DISCRIMINANT_FACTOR_COEFFICIENTS,
            ball.residue,
            modulus,
        ),
        **{
            name: affine_variable_coefficients(polynomial, ball.residue, modulus)
            for name, polynomial in INVARIANT_POLYNOMIALS.items()
            if name in ("a4", "a6", "c4", "c6")
        },
    }
    forced_h = fixed_divisor_valuation(transformed["H"], prime)
    contents = {
        name: coefficient_content_valuation(polynomial, prime)
        for name, polynomial in transformed.items()
    }
    scaling = min(
        contents["a4"] // 4,
        contents["a6"] // 6,
        forced_h // 12,
    )
    minimal_delta = forced_h - 12 * scaling
    if minimal_delta < 1:
        return ClassifiedBall(
            prime,
            ball.exponent,
            ball.residue,
            forced_h,
            contents["H"],
            scaling,
            minimal_delta,
            "not-forced-bad",
            (),
            "the ball does not force positive minimal-discriminant valuation",
        )

    a4_divisor = prime ** (4 * scaling)
    a6_divisor = prime ** (6 * scaling)
    a4_scaled = tuple(coefficient // a4_divisor for coefficient in transformed["a4"])
    a6_scaled = tuple(coefficient // a6_divisor for coefficient in transformed["a6"])
    c4_divisor = prime ** (4 * scaling)
    c4_scaled = tuple(coefficient // c4_divisor for coefficient in transformed["c4"])
    c4_values = [
        polynomial_value(c4_scaled, digit) % prime for digit in range(prime)
    ]

    if all(c4_values):
        split_symbols: set[bool] = set()
        for digit in range(prime):
            coefficient_a = polynomial_value(a4_scaled, digit) % prime
            coefficient_b = polynomial_value(a6_scaled, digit) % prime
            if coefficient_a == 0:
                raise AssertionError("unit c4 unexpectedly gave a4=0 modulo p")
            double_root = (
                -3 * coefficient_b * pow(2 * coefficient_a, -1, prime)
            ) % prime
            tangent_symbol = legendre_symbol(3 * double_root, prime)
            if tangent_symbol == 0:
                raise AssertionError("a multiplicative fiber became cuspidal")
            split_symbols.add(tangent_symbol == 1)
        if split_symbols == {True}:
            reduction = "split multiplicative"
        elif split_symbols == {False}:
            reduction = "nonsplit multiplicative"
        else:
            reduction = "mixed multiplicative split symbol"
        return ClassifiedBall(
            prime,
            ball.exponent,
            ball.residue,
            forced_h,
            contents["H"],
            scaling,
            minimal_delta,
            reduction,
            tuple(sorted(split_symbols)),
            "minimal c4 is a unit on every next p-adic digit",
        )

    exact_a4 = uniform_exact_valuation(transformed["a4"], prime)
    exact_a6 = uniform_exact_valuation(transformed["a6"], prime)
    exact_h = uniform_exact_valuation(transformed["H"], prime)
    minimality_obstructions = (
        exact_a4 is not None and exact_a4 - 4 * scaling < 4,
        exact_a6 is not None and exact_a6 - 6 * scaling < 6,
        exact_h is not None and exact_h - 12 * scaling < 12,
    )
    if all(value == 0 for value in c4_values) and any(minimality_obstructions):
        return ClassifiedBall(
            prime,
            ball.exponent,
            ball.residue,
            forced_h,
            contents["H"],
            scaling,
            minimal_delta,
            "additive",
            (),
            (
                "minimal c4 is divisible by p throughout and a uniform exact "
                "invariant valuation prevents another integral scaling"
            ),
        )

    return ClassifiedBall(
        prime,
        ball.exponent,
        ball.residue,
        forced_h,
        contents["H"],
        scaling,
        minimal_delta,
        "mixed or unresolved",
        (),
        "the current ball contains more than one provable local behavior",
    )


def _contains(outer: RootBall, inner: RootBall) -> bool:
    return (
        outer.prime == inner.prime
        and outer.exponent <= inner.exponent
        and inner.residue % outer.modulus == outer.residue
    )


def power_efficient_root_balls(
    prime: int,
    lift_exponent: int,
    *,
    max_roots: int,
) -> tuple[tuple[RootBall, ...], list[dict[str, Any]]]:
    """Discover coarsest balls whose forced H-valuation exceeds their cost exponent."""

    profiles: list[dict[str, Any]] = []
    candidates: set[RootBall] = set()
    for exponent in range(2, lift_exponent + 1):
        result = all_roots_mod_prime_power(
            DISCRIMINANT_FACTOR_COEFFICIENTS,
            prime,
            exponent,
            max_roots=max_roots,
            verify=True,
        )
        balls = result.maximal_balls()
        profiles.append(
            {
                "target_exponent": exponent,
                "root_count": len(result.roots),
                "level_counts": list(result.level_counts),
                "compressed_ball_count": len(balls),
                "candidate_digits_checked": result.candidate_digits_checked,
            }
        )
        for ball in balls:
            forced = fixed_divisor_valuation(
                affine_variable_coefficients(
                    DISCRIMINANT_FACTOR_COEFFICIENTS,
                    ball.residue,
                    ball.modulus,
                ),
                prime,
            )
            if forced > ball.exponent:
                candidates.add(ball)

    coarsest = tuple(
        sorted(
            (
                ball
                for ball in candidates
                if not any(
                    other != ball and _contains(other, ball)
                    for other in candidates
                )
            ),
            key=lambda ball: (ball.exponent, ball.residue),
        )
    )
    return coarsest, profiles


def refine_and_classify(
    ball: RootBall,
    *,
    maximum_exponent: int,
) -> tuple[ClassifiedBall, ...]:
    """Refine a mixed ball digit-by-digit, retaining exact exhaustive coverage."""

    classified = classify_ball(ball)
    if (
        classified.reduction not in (
            "mixed or unresolved",
            "mixed multiplicative split symbol",
        )
        or ball.exponent >= maximum_exponent
    ):
        return (classified,)
    children: list[ClassifiedBall] = []
    for digit in range(ball.prime):
        child = RootBall(
            ball.prime,
            ball.exponent + 1,
            ball.residue + digit * ball.modulus,
        )
        children.extend(
            refine_and_classify(child, maximum_exponent=maximum_exponent)
        )
    return tuple(children)


@dataclass(frozen=True)
class ConditionGroup:
    prime: int
    exponent: int
    residues: tuple[int, ...]
    forced_h_valuation: int
    presented_model_scaling: int
    forced_minimal_discriminant_valuation: int
    unconditional_h_valuation: int
    reduction: str
    reciprocal_density: Fraction
    log_congruence_cost: float
    radical_saving_log: float
    efficiency: float

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    def serializable(self) -> dict[str, Any]:
        answer = asdict(self)
        answer["modulus"] = self.modulus
        answer["reciprocal_density"] = str(self.reciprocal_density)
        return answer


def group_classified_balls(
    balls: Sequence[ClassifiedBall],
) -> tuple[ConditionGroup, ...]:
    """Union balls with identical exact local guarantees at one prime."""

    buckets: dict[tuple[Any, ...], list[int]] = {}
    for ball in balls:
        key = (
            ball.prime,
            ball.exponent,
            ball.forced_h_valuation,
            ball.presented_model_scaling,
            ball.forced_minimal_discriminant_valuation,
            ball.reduction,
        )
        buckets.setdefault(key, []).append(ball.residue)
    groups: list[ConditionGroup] = []
    for key, residues in buckets.items():
        prime, exponent, forced_h, scaling, minimal_delta, reduction = key
        canonical_residues = tuple(sorted(set(residues)))
        reciprocal_density = Fraction(prime**exponent, len(canonical_residues))
        log_cost = log(float(reciprocal_density))
        unconditional_h = fixed_divisor_valuation(
            DISCRIMINANT_FACTOR_COEFFICIENTS,
            prime,
        )
        # Multiplicative reduction pays one copy of p in the conductor.  The
        # remainder measures deliberately manufactured radical saving.  Any
        # unconditional fixed divisor is baseline rather than a discovery
        # benefit (notably v_5(H)>=2 for every specialization).
        radical_saving = (
            max(0, minimal_delta - max(1, unconditional_h)) * log(prime)
        )
        efficiency = radical_saving / log_cost if log_cost > 0 else 0.0
        groups.append(
            ConditionGroup(
                prime,
                exponent,
                canonical_residues,
                forced_h,
                scaling,
                minimal_delta,
                unconditional_h,
                reduction,
                reciprocal_density,
                log_cost,
                radical_saving,
                efficiency,
            )
        )
    return tuple(
        sorted(
            groups,
            key=lambda group: (
                -group.efficiency,
                group.prime,
                group.exponent,
                group.residues,
            ),
        )
    )


def discover_prime(
    prime: int,
    *,
    lift_exponent: int,
    classification_exponent: int,
    max_roots: int,
) -> dict[str, Any]:
    efficient_balls, lift_profiles = power_efficient_root_balls(
        prime,
        lift_exponent,
        max_roots=max_roots,
    )
    classified: list[ClassifiedBall] = []
    for ball in efficient_balls:
        classified.extend(
            refine_and_classify(
                ball,
                maximum_exponent=max(classification_exponent, ball.exponent),
            )
        )
    groups = group_classified_balls(classified)
    return {
        "prime": prime,
        "lift_profiles": lift_profiles,
        "power_efficient_balls_before_classification_refinement": [
            asdict(ball) | {"modulus": ball.modulus} for ball in efficient_balls
        ],
        "classified_balls": [
            asdict(ball) | {"modulus": ball.modulus} for ball in classified
        ],
        "condition_groups": [group.serializable() for group in groups],
        "_groups": groups,
    }


def select_condition_groups(
    groups: Sequence[ConditionGroup],
    *,
    count: int,
    maximum_crt_classes: int,
) -> tuple[ConditionGroup, ...]:
    """Greedily select the best clean split group at distinct primes."""

    eligible = sorted(
        (
            group
            for group in groups
            if group.reduction == "split multiplicative"
            and group.forced_minimal_discriminant_valuation >= 2
            and group.reciprocal_density > 1
        ),
        key=lambda group: (
            -group.efficiency,
            group.prime,
            group.exponent,
            group.residues,
        ),
    )
    selected: list[ConditionGroup] = []
    selected_primes: set[int] = set()
    class_count = 1
    for group in eligible:
        if group.prime in selected_primes:
            continue
        proposed = class_count * len(group.residues)
        if proposed > maximum_crt_classes:
            continue
        selected.append(group)
        selected_primes.add(group.prime)
        class_count = proposed
        if len(selected) == count:
            break
    if len(selected) < count:
        raise ValueError("not enough eligible groups satisfy the CRT class cap")
    return tuple(selected)


def validate_invariant_polynomials() -> None:
    """Cross-check the constructed polynomials against the family implementation."""

    for parameter in (-17, 0, 23):
        invariant_i, invariant_j = FermigierMestreFamily.binary_invariants(
            Fraction(parameter)
        )
        if polynomial_value(INVARIANT_POLYNOMIALS["I"], parameter) != invariant_i:
            raise AssertionError("the constructed I(T) polynomial is incorrect")
        if polynomial_value(INVARIANT_POLYNOMIALS["J"], parameter) != invariant_j:
            raise AssertionError("the constructed J(T) polynomial is incorrect")
