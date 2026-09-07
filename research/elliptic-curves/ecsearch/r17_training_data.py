"""Leakage-controlled cheap features for the published R17 family.

This module builds discovery data, not rank certificates.  It deliberately
does not know the ranks or points on the four published high-rank fibres.  The
parameters themselves appear only in :data:`EMBARGOED_PARAMETERS`, so callers
can prove that they were excluded from score development.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, log, sqrt
from pathlib import Path
from random import Random
from typing import Iterable, Sequence


EMBARGOED_PARAMETERS = frozenset(
    {
        (-2, 377),
        (-308, 251),
        (2456, 135),
        (-9529, 5471),
    }
)


def normalized_parameter(a: int, b: int) -> tuple[int, int]:
    if a == 0 and b == 0:
        raise ValueError("(0:0) is not a projective parameter")
    common = gcd(abs(a), abs(b))
    a //= common
    b //= common
    if b < 0 or (b == 0 and a < 0):
        a, b = -a, -b
    return a, b


def parameter_text(parameter: tuple[int, int]) -> str:
    a, b = parameter
    return "infinity" if b == 0 else f"{a}/{b}"


def deterministic_sample(
    *, count: int, height: int, seed: int, excluded: Iterable[tuple[int, int]] = ()
) -> list[tuple[int, int]]:
    """Sample canonical finite parameters uniformly from a bounded box.

    Proposals are uniform in ``[-height,height] x [1,height]`` and are accepted
    exactly when primitive, new, and not excluded.  Conditional on acceptance,
    every finite primitive projective parameter in the box has equal weight.
    """

    if count < 1 or height < 1:
        raise ValueError("count and height must be positive")
    excluded_set = {normalized_parameter(*parameter) for parameter in excluded}
    available_upper_bound = (2 * height + 1) * height - len(excluded_set)
    if count > available_upper_bound:
        raise ValueError("requested sample is larger than the parameter box")
    generator = Random(seed)
    selected: set[tuple[int, int]] = set()
    while len(selected) < count:
        parameter = (generator.randint(-height, height), generator.randint(1, height))
        if gcd(abs(parameter[0]), parameter[1]) != 1:
            continue
        if parameter in excluded_set:
            continue
        selected.add(parameter)
    return sorted(selected, key=lambda item: (max(abs(item[0]), item[1]), item[1], item[0]))


def split_bucket(parameter: tuple[int, int], salt: str) -> str:
    digest = sha256(f"{salt}|{parameter[0]}/{parameter[1]}".encode()).digest()
    value = int.from_bytes(digest[:8], "big") % 100
    if value < 70:
        return "train"
    if value < 85:
        return "validation"
    return "internal_test"


def projective_index(parameter: tuple[int, int], prime: int) -> int:
    a, b = parameter
    if b % prime == 0:
        return prime
    return a % prime * pow(b % prime, -1, prime) % prime


@dataclass(frozen=True)
class CheapPrimeSymbol:
    contribution: float
    good: bool
    two_quotient_dimension: int
    three_quotient_dimension: int


@dataclass(frozen=True)
class CheapPrimeTable:
    rows: tuple[CheapPrimeSymbol, ...]
    good_mean: float
    good_standard_deviation: float

    def __getitem__(self, index: int) -> CheapPrimeSymbol:
        return self.rows[index]


def _legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    answer = pow(value, (prime - 1) // 2, prime)
    return -1 if answer == prime - 1 else 1


def _torsion_quotient_dimensions(coefficient_a: int, coefficient_b: int, prime: int) -> tuple[int, int]:
    roots_two = sum(
        (x**3 + coefficient_a * x + coefficient_b) % prime == 0
        for x in range(prime)
    )
    two_dimension = {0: 0, 1: 1, 3: 2}[roots_two]
    roots_three = 0
    for x in range(prime):
        division_value = (
            3 * x**4
            + 6 * coefficient_a * x**2
            + 12 * coefficient_b * x
            - coefficient_a**2
        ) % prime
        if division_value != 0:
            continue
        curve_value = (x**3 + coefficient_a * x + coefficient_b) % prime
        if _legendre(curve_value, prime) == 1:
            roots_three += 1
    three_torsion_size = 1 + 2 * roots_three
    three_dimension = {1: 0, 3: 1, 9: 2}[three_torsion_size]
    return two_dimension, three_dimension


def build_cheap_prime_tables(model, primes: Sequence[int]) -> dict[int, CheapPrimeTable]:
    """Precompute Nagao and finite-quotient structure on ``P^1(F_p)``."""

    # Imported lazily so the module remains easy to unit test in isolation.
    from search_h92_q12o5867_rootless_nagao import reduced_coefficients, residue_table

    tables: dict[int, CheapPrimeTable] = {}
    for prime in primes:
        nagao_table = residue_table(model, prime)
        a_coefficients, b_coefficients = reduced_coefficients(model, prime)
        rows = []
        for index, symbol in enumerate(nagao_table):
            if not symbol.good_reduction:
                rows.append(CheapPrimeSymbol(0.0, False, -1, -1))
                continue
            numerator, denominator = ((index, 1) if index < prime else (1, 0))
            coefficient_a = sum(
                coefficient
                * pow(numerator, power, prime)
                * pow(denominator, 8 - power, prime)
                for power, coefficient in enumerate(a_coefficients)
            ) % prime
            coefficient_b = sum(
                coefficient
                * pow(numerator, power, prime)
                * pow(denominator, 12 - power, prime)
                for power, coefficient in enumerate(b_coefficients)
            ) % prime
            dimension_two, dimension_three = _torsion_quotient_dimensions(
                coefficient_a, coefficient_b, prime
            )
            rows.append(
                CheapPrimeSymbol(
                    symbol.contribution_units / 10**12,
                    True,
                    dimension_two,
                    dimension_three,
                )
            )
        good_values = [row.contribution for row in rows if row.good]
        mean = sum(good_values) / len(good_values)
        variance = sum((value - mean) ** 2 for value in good_values) / len(good_values)
        tables[prime] = CheapPrimeTable(tuple(rows), mean, sqrt(variance))
    return tables


def nagao_features(
    parameter: tuple[int, int],
    tables: dict[int, CheapPrimeTable],
    prime_blocks: Sequence[Sequence[int]],
) -> dict[str, object]:
    block_signals: list[float] = []
    good_count = 0
    bad_count = 0
    for block in prime_blocks:
        normalized_sum = 0.0
        for prime in block:
            table = tables[prime]
            symbol = table[projective_index(parameter, prime)]
            if symbol.good:
                normalized_sum += (
                    symbol.contribution - table.good_mean
                ) / table.good_standard_deviation
                good_count += 1
            else:
                bad_count += 1
        block_signals.append(normalized_sum / sqrt(len(block)))
    return {
        "standardized_block_signals": block_signals,
        "worst_block_signal": min(block_signals),
        "mean_block_signal": sum(block_signals) / len(block_signals),
        "good_prime_count": good_count,
        "bad_prime_count": bad_count,
    }


def quotient_code_features(
    parameter: tuple[int, int],
    tables: dict[int, CheapPrimeTable],
    primes: Sequence[int],
) -> dict[str, object]:
    code = []
    for prime in primes:
        symbol = tables[prime][projective_index(parameter, prime)]
        code.append([symbol.two_quotient_dimension, symbol.three_quotient_dimension])
    return {
        "primes": list(primes),
        "local_E_mod_2_and_mod_3_dimensions": code,
        "code": ",".join(f"{two}/{three}" for two, three in code),
    }


def homogeneous_integer_value(coefficients: Sequence[int], parameter: tuple[int, int], degree: int) -> int:
    a, b = parameter
    return sum(
        coefficient * a**power * b ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def _valuation(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def conductor_proxy_features(
    parameter: tuple[int, int],
    a_coefficients: Sequence[int],
    b_coefficients: Sequence[int],
    primes: Sequence[int],
) -> dict[str, object]:
    """Return a partial local ``log|Delta|-log(N)`` saving proxy.

    For ``p>=5`` the short model is first scaled to a p-minimal short model.
    The conductor exponent is then 0 for good, 1 for multiplicative, and 2 for
    additive reduction.  The omitted primes and the wild primes 2 and 3 make
    this a ranking proxy, never an exact global conductor.
    """

    coefficient_a = homogeneous_integer_value(a_coefficients, parameter, 8)
    coefficient_b = homogeneous_integer_value(b_coefficients, parameter, 12)
    discriminant_core = 4 * coefficient_a**3 + 27 * coefficient_b**2
    if discriminant_core == 0:
        raise ValueError("singular rational fibre")
    local_rows = []
    saving = 0.0
    scaling_log = 0.0
    for prime in primes:
        valuation_a = _valuation(coefficient_a, prime)
        valuation_b = _valuation(coefficient_b, prime)
        scale = min(valuation_a // 4, valuation_b // 6)
        locally_minimal_a = coefficient_a // prime ** (4 * scale)
        valuation_delta = _valuation(discriminant_core, prime) - 12 * scale
        if prime == 2:
            valuation_delta += 4
        conductor_exponent = None
        if prime >= 5:
            if valuation_delta == 0:
                conductor_exponent = 0
            elif locally_minimal_a % prime != 0:
                conductor_exponent = 1
            else:
                conductor_exponent = 2
            saving += max(valuation_delta - conductor_exponent, 0) * log(prime)
        scaling_log += 12 * scale * log(prime)
        local_rows.append([prime, valuation_delta, conductor_exponent, scale])
    log_discriminant_after_known_scaling = log(16 * abs(discriminant_core)) - scaling_log
    return {
        "local_rows_p_vdelta_f_scale": local_rows,
        "known_prime_log_discriminant_saving": saving,
        "log_discriminant_after_known_scaling": log_discriminant_after_known_scaling,
        "quality_proxy": saving / max(log_discriminant_after_known_scaling, 1.0),
        "boundary": "partial p>=5 local saving; not an exact conductor",
    }


def normalize_quadratic(
    coefficients: Sequence[str],
) -> tuple[tuple[int, int, int], int]:
    """Return a primitive polynomial and the exact rational scalar squareclass."""

    fractions = [Fraction(value) for value in coefficients]
    denominator = 1
    for value in fractions:
        denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in fractions]
    common = gcd(gcd(abs(integers[0]), abs(integers[1])), abs(integers[2]))
    integers = [value // common for value in integers]
    scalar = Fraction(common, denominator)
    if next(value for value in reversed(integers) if value) < 0:
        integers = [-value for value in integers]
        scalar = -scalar
    scalar_squareclass = scalar.numerator * scalar.denominator
    return tuple(integers), scalar_squareclass  # type: ignore[return-value]


def split_quadratic_indices(
    parameter: tuple[int, int],
    quadratics: Sequence[tuple[tuple[int, int, int], int]],
) -> tuple[list[int], list[int]]:
    """Return the nonzero split and ramified indices in an exact cover list.

    A normalized entry represents ``scalar * (q0 + q1*t + q2*t^2)`` modulo
    rational squares.  At ``t=a/b`` its squareclass is therefore represented
    by the integer ``scalar_squareclass * (q0*b^2+q1*a*b+q2*a^2)``.  Testing
    that integer avoids constructing a :class:`Fraction` for every atlas row.
    Zero values are reported separately because they are ramified
    specializations, not two-branch split covers.
    """

    a, b = normalized_parameter(*parameter)
    split: list[int] = []
    ramified: list[int] = []
    for index, ((q0, q1, q2), scalar_squareclass) in enumerate(quadratics):
        value = scalar_squareclass * (q0 * b * b + q1 * a * b + q2 * a * a)
        if value == 0:
            ramified.append(index)
        elif value > 0 and isqrt(value) ** 2 == value:
            split.append(index)
    return split, ramified


def select_cover_panel(
    bisections: Sequence[dict], count: int
) -> list[tuple[str, tuple[tuple[int, int, int], int]]]:
    if count < 1 or count > len(bisections):
        raise ValueError("invalid cover panel size")
    ordered = sorted(
        bisections,
        key=lambda row: sha256(str(row["label"]).encode()).digest(),
    )[:count]
    return [
        (
            str(row["label"]),
            normalize_quadratic(row["residual_chord"]["q_coefficients"]),
        )
        for row in ordered
    ]


def cover_character_tables(
    panel: Sequence[tuple[str, tuple[tuple[int, int, int], int]]],
    primes: Sequence[int],
) -> dict[int, tuple[tuple[int, ...], ...]]:
    tables = {}
    for prime in primes:
        rows = []
        for index in range(prime + 1):
            a, b = ((index, 1) if index < prime else (1, 0))
            signs = tuple(
                _legendre(
                    scalar * (q0 * b * b + q1 * a * b + q2 * a * a),
                    prime,
                )
                for _label, ((q0, q1, q2), scalar) in panel
            )
            rows.append(signs)
        tables[prime] = tuple(rows)
    return tables


def cover_diversity_features(
    parameter: tuple[int, int],
    tables: dict[int, tuple[tuple[int, ...], ...]],
    primes: Sequence[int],
) -> dict[str, object]:
    per_prime = [tables[prime][projective_index(parameter, prime)] for prime in primes]
    patterns = list(zip(*per_prime))
    frequencies = Counter(patterns)
    entropy = -sum(
        count / len(patterns) * log(count / len(patterns)) for count in frequencies.values()
    )
    return {
        "primes": list(primes),
        "distinct_character_patterns": len(frequencies),
        "character_pattern_entropy": entropy,
        "branch_zero_count": sum(value == 0 for pattern in patterns for value in pattern),
        "boundary": "modular character diversity on a fixed cover panel; not a split-cover count",
    }


def development_lane_memberships(records: list[dict], lane_size: int, salt: str) -> dict[str, list[int]]:
    if not 1 <= lane_size <= len(records):
        raise ValueError("lane_size must lie between one and the population size")

    def nagao_key(index: int) -> tuple:
        row = records[index]["features"]["level1_nagao"]
        return (row["worst_block_signal"], row["mean_block_signal"], -records[index]["height"])

    split_quotas = {
        "train": round(0.70 * lane_size),
        "validation": round(0.15 * lane_size),
    }
    split_quotas["internal_test"] = lane_size - sum(split_quotas.values())
    split_indices = {
        split: [index for index, row in enumerate(records) if row["split"] == split]
        for split in split_quotas
    }
    if any(len(split_indices[split]) < quota for split, quota in split_quotas.items()):
        raise ValueError("a data split is too small for the requested lane quota")

    def stratified_rank(key, *, reverse: bool) -> list[int]:
        return [
            index
            for split, quota in split_quotas.items()
            for index in sorted(split_indices[split], key=key, reverse=reverse)[:quota]
        ]

    lanes = {
        "top_nagao": stratified_rank(nagao_key, reverse=True),
        "high_conductor_quality_proxy": stratified_rank(
            key=lambda index: (
                records[index]["features"]["level2_conductor_proxy"]["quality_proxy"],
                nagao_key(index),
            ),
            reverse=True,
        ),
        "unusual_quotient_code": stratified_rank(
            key=lambda index: (
                records[index]["features"]["level2_quotient_code"]["rarity"],
                nagao_key(index),
            ),
            reverse=True,
        ),
        "high_cover_diversity": stratified_rank(
            key=lambda index: (
                records[index]["features"]["level2_cover_diversity"]["distinct_character_patterns"],
                records[index]["features"]["level2_cover_diversity"]["character_pattern_entropy"],
                nagao_key(index),
            ),
            reverse=True,
        ),
    }
    lanes["random_controls"] = stratified_rank(
        key=lambda index: sha256(
            f"{salt}|random-control|{records[index]['parameter']}".encode()
        ).digest(),
        reverse=False,
    )
    return lanes
