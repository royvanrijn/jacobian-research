#!/usr/bin/env python3
"""Third-generation rank-gain search in Nagao's rank-13 family.

The first two rank-gain searches exhaustively screened a Farey grid through
denominator 12 and mutations about ``u=42,118`` through denominator 64.  This
experiment is disjoint from both populations.  It searches three new,
deterministic sources of rational parameters:

* a Pell-like shell close to ``u^2=23550``, where a large-denominator ``u``
  can nevertheless give a relatively small base parameter ``T``;
* local CRT classes whose reductions at 29, 41, and 43 agree with the most
  favourable residues of the four exactly certified rank-17 specializations;
* high-denominator tangent mutations about all four exact rank-17 labels.

Before point searching, candidates are selected by a declared multiobjective
union: rank-17 residue affinity, a short-prime Nagao statistic, small-prime
discriminant-power savings, and the height of ``T``.  Later selection uses
only exact quartic-point yield (with all eighteen classified linear sections
removed) and precision-stable numerical height rank.  The final stage adds
eighteen section-centred charts, so large global ``x`` height is not the only
geometry searched.

Every subprocess is a foreground process group with a strict timeout.  Point
membership and bounded enumerations are exact.  Height ranks are numerical.
If a stable rank at least 18 appears, the script immediately attempts an
unconditional finite-reduction independence certificate; it does not promote
the numerical rank by itself.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, gcd, isqrt, log
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair
from ek_k3 import rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_BASE_CHANGE_CONSTANT,
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from nagao_linear_sections import omitted_companion_sections
from nagao_rank13_local import BASE_CHANGED_DISCRIMINANT, polynomial_value_mod
from pari_bridge import pari_version
from search_extra_points import parse_point_vector
from search_nagao_rank13_rank_gain import (
    ParameterCandidate,
    PointScreen,
    _screen_record,
    batch_point_screen,
    canonical_positive_u,
    conductor_probe,
    evaluate_height_rank,
    quartic_gp_polynomial,
    run_gp_capped,
)
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXACT_RANK17_PARAMETERS = (Q(135, 2), Q(471, 11), Q(42), Q(74))
LABEL_PRIMES = (29, 41, 43)
LABEL_RESIDUE_COUNTS = {29: 3, 41: 2, 43: 2}
PREDICTOR_PRIMES = (29, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
POWER_PRIMES = (7, 11, 13, 17, 19, 23, 31)


@dataclass(frozen=True)
class CandidateFeatures:
    parameter_u: Fraction
    origins: tuple[str, ...]
    parameter_t: Fraction
    parameter_t_height: int
    rank17_residue_affinity: float
    short_prime_score: float
    discriminant_power_savings: float
    combined_objective: float

    @property
    def candidate(self) -> ParameterCandidate:
        return ParameterCandidate(self.parameter_u, self.origins)

    def compact_json(self) -> dict[str, Any]:
        return {
            "parameter_u": rational_to_string(self.parameter_u),
            "parameter_t": rational_to_string(self.parameter_t),
            "origins": list(self.origins),
            "parameter_t_height": self.parameter_t_height,
            "rank17_residue_affinity": self.rank17_residue_affinity,
            "short_prime_score": self.short_prime_score,
            "discriminant_power_savings": self.discriminant_power_savings,
            "combined_objective": self.combined_objective,
        }


def _fraction_mod(value: Fraction, prime: int) -> int | None:
    value = Q(value)
    if value.denominator % prime == 0:
        return None
    return value.numerator * pow(value.denominator, -1, prime) % prime


def _local_trace(parameter_u_residue: int, prime: int) -> int | None:
    """Return ``a_p`` for a nonzero ``u mod p``, or ``None`` if bad."""

    if not 0 < parameter_u_residue < prime:
        return None
    coefficients = rank13_base_changed_short_jacobian_coefficients(
        Q(parameter_u_residue)
    )
    coefficient_a_q, coefficient_b_q = coefficients[3:]
    coefficient_a = _fraction_mod(coefficient_a_q, prime)
    coefficient_b = _fraction_mod(coefficient_b_q, prime)
    if coefficient_a is None or coefficient_b is None:
        return None
    discriminant = -16 * (
        4 * coefficient_a**3 + 27 * coefficient_b**2
    )
    if discriminant % prime == 0:
        return None
    character_sum = 0
    for abscissa in range(prime):
        rhs = (abscissa**3 + coefficient_a * abscissa + coefficient_b) % prime
        if rhs:
            character_sum += 1 if pow(rhs, (prime - 1) // 2, prime) == 1 else -1
    return -character_sum


def local_trace_tables(
    primes: Sequence[int] = PREDICTOR_PRIMES,
) -> dict[int, dict[int, int | None]]:
    return {
        prime: {
            residue: _local_trace(residue, prime)
            for residue in range(1, prime)
        }
        for prime in primes
    }


def rank17_label_residues(
    trace_tables: dict[int, dict[int, int | None]],
) -> dict[int, tuple[int, ...]]:
    """Select favourable residues actually occupied by exact rank-17 labels."""

    answer: dict[int, tuple[int, ...]] = {}
    for prime in LABEL_PRIMES:
        labelled = set()
        for parameter_u in EXACT_RANK17_PARAMETERS:
            residue = _fraction_mod(parameter_u, prime)
            if residue not in (None, 0) and trace_tables[prime][residue] is not None:
                labelled.add(residue)
        ordered = sorted(
            labelled,
            key=lambda residue: (trace_tables[prime][residue], residue),
        )
        retained = tuple(ordered[: LABEL_RESIDUE_COUNTS[prime]])
        if len(retained) != LABEL_RESIDUE_COUNTS[prime]:
            raise AssertionError(f"too few usable exact-rank-17 residues at p={prime}")
        answer[prime] = retained
    return answer


def label_crt_classes(
    residue_groups: dict[int, tuple[int, ...]],
) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    classes = []
    for residues in product(*(residue_groups[prime] for prime in LABEL_PRIMES)):
        combined_residue, modulus = 0, 1
        for prime, residue in zip(LABEL_PRIMES, residues):
            combined_residue, modulus = crt_pair(
                combined_residue, modulus, residue, prime
            )
        classes.append((combined_residue, modulus, tuple(residues)))
    return tuple(classes)


def load_prior_population(paths: Sequence[Path]) -> set[Fraction]:
    prior: set[Fraction] = set()
    for path in paths:
        data = json.loads(path.read_text())
        records = data["initial_point_screen"]["records"]
        prior.update(Q(record["parameter_u"]) for record in records)
    return prior


def load_exact_rank17_labels(path: Path) -> tuple[Fraction, ...]:
    data = json.loads(path.read_text())
    labels = []
    for record in data["certificates"]:
        if (
            int(
                record["finite_reduction_certificate"][
                    "certified_algebraic_rank_lower_bound"
                ]
            )
            != 17
        ):
            raise AssertionError("an exact frontier label is not rank at least 17")
        labels.append(Q(record["parameter_u"]))
    if tuple(labels) != EXACT_RANK17_PARAMETERS:
        raise AssertionError("the exact rank-17 label set changed")
    return tuple(labels)


def generate_generation3_population(
    residue_groups: dict[int, tuple[int, ...]],
    prior_population: set[Fraction],
    *,
    pell_denominator_max: int = 1024,
    pell_offset_radius: int = 8,
    crt_denominator_max: int = 512,
    mutation_denominator_max: int = 384,
    mutation_radius: int = 8,
) -> tuple[ParameterCandidate, ...]:
    """Return the symmetry-normalized population, disjoint from generations 1/2."""

    if min(
        pell_denominator_max,
        crt_denominator_max,
        mutation_denominator_max,
    ) < 65:
        raise ValueError("generation 3 starts strictly above denominator 64")
    origins: dict[Fraction, set[str]] = {}

    def add(parameter_u: Fraction, origin: str) -> None:
        if parameter_u <= 0:
            return
        parameter_u = canonical_positive_u(parameter_u)
        if parameter_u.denominator <= 64 or parameter_u in prior_population:
            return
        origins.setdefault(parameter_u, set()).add(origin)

    constant = RANK13_BASE_CHANGE_CONSTANT
    for denominator in range(65, pell_denominator_max + 1):
        boundary = isqrt(constant * denominator**2)
        for offset in range(-pell_offset_radius, pell_offset_radius + 1):
            add(Q(boundary + offset, denominator), "pell-small-T-shell")

    for residue, modulus, _ in label_crt_classes(residue_groups):
        for denominator in range(65, crt_denominator_max + 1):
            if gcd(denominator, modulus) != 1:
                continue
            maximum_numerator = isqrt(constant * denominator**2 - 1)
            base_numerator = residue * denominator % modulus
            minimum_step = (1 - base_numerator + modulus - 1) // modulus
            maximum_step = (maximum_numerator - base_numerator) // modulus
            for step in range(minimum_step, maximum_step + 1):
                add(
                    Q(base_numerator + step * modulus, denominator),
                    "rank17-local-crt",
                )

    for center in EXACT_RANK17_PARAMETERS:
        for denominator in range(65, mutation_denominator_max + 1):
            for numerator_delta in range(-mutation_radius, mutation_radius + 1):
                if numerator_delta:
                    add(
                        center + Q(numerator_delta, denominator),
                        f"rank17-tangent-{rational_to_string(center)}",
                    )

    population = tuple(
        ParameterCandidate(parameter_u, tuple(sorted(labels)))
        for parameter_u, labels in sorted(origins.items())
    )
    if not population or any(
        candidate.parameter_u in prior_population
        or candidate.parameter_u.denominator <= 64
        for candidate in population
    ):
        raise AssertionError("generation 3 overlaps a prior population")
    return population


def _discriminant_valuation_proxy(
    parameter_u: Fraction, prime: int, *, exponent_cap: int = 6
) -> int:
    if parameter_u.denominator % prime == 0:
        return 0
    modulus = prime**exponent_cap
    residue = (
        parameter_u.numerator
        * pow(parameter_u.denominator, -1, modulus)
        % modulus
    )
    value = polynomial_value_mod(BASE_CHANGED_DISCRIMINANT, residue, modulus)
    if value == 0:
        return exponent_cap
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def candidate_features(
    candidate: ParameterCandidate,
    trace_tables: dict[int, dict[int, int | None]],
    residue_groups: dict[int, tuple[int, ...]],
) -> CandidateFeatures:
    parameter_u = candidate.parameter_u
    affinity = 0.0
    short_score = 0.0
    for prime in PREDICTOR_PRIMES:
        residue = _fraction_mod(parameter_u, prime)
        if residue in (None, 0):
            continue
        trace = trace_tables[prime][residue]
        if trace is not None:
            short_score += (2 - trace) / (prime + 1 - trace) * log(prime)
        if prime in residue_groups and residue in residue_groups[prime]:
            affinity += log(prime / len(residue_groups[prime]))
    power_savings = sum(
        max(0, _discriminant_valuation_proxy(parameter_u, prime) - 1)
        * log(prime)
        for prime in POWER_PRIMES
    )
    parameter_t = candidate.parameter_t
    parameter_t_height = max(abs(parameter_t.numerator), parameter_t.denominator)
    combined = (
        2.0 * affinity
        + short_score
        + 0.35 * power_savings
        - 0.4 * log(max(2, parameter_t_height))
    )
    return CandidateFeatures(
        parameter_u=parameter_u,
        origins=candidate.origins,
        parameter_t=parameter_t,
        parameter_t_height=parameter_t_height,
        rank17_residue_affinity=affinity,
        short_prime_score=short_score,
        discriminant_power_savings=power_savings,
        combined_objective=combined,
    )


def multiobjective_prefilter(
    features: Sequence[CandidateFeatures],
    *,
    combined_keep: int,
    axis_keep: int,
    origin_keep: int,
) -> tuple[ParameterCandidate, ...]:
    if min(combined_keep, axis_keep, origin_keep) <= 0:
        raise ValueError("prefilter keep counts must be positive")
    selected: set[Fraction] = set()

    def retain(ordered: Iterable[CandidateFeatures], count: int) -> None:
        selected.update(item.parameter_u for item in list(ordered)[:count])

    retain(
        sorted(features, key=lambda item: (-item.combined_objective, item.parameter_u)),
        combined_keep,
    )
    retain(
        sorted(features, key=lambda item: (-item.short_prime_score, item.parameter_u)),
        axis_keep,
    )
    retain(
        sorted(
            features,
            key=lambda item: (-item.discriminant_power_savings, item.parameter_u),
        ),
        axis_keep,
    )
    retain(
        sorted(features, key=lambda item: (item.parameter_t_height, item.parameter_u)),
        axis_keep,
    )
    origins = sorted({origin for item in features for origin in item.origins})
    for origin in origins:
        retain(
            sorted(
                (item for item in features if origin in item.origins),
                key=lambda item: (-item.combined_objective, item.parameter_u),
            ),
            origin_keep,
        )
    by_u = {item.parameter_u: item.candidate for item in features}
    return tuple(by_u[parameter_u] for parameter_u in sorted(selected))


def _feature_sort_map(
    features: Sequence[CandidateFeatures],
) -> dict[Fraction, CandidateFeatures]:
    return {item.parameter_u: item for item in features}


def select_point_survivors_generation3(
    screens: Sequence[PointScreen],
    features: Sequence[CandidateFeatures],
    *,
    point_keep: int,
    objective_keep: int,
    origin_keep: int,
) -> tuple[ParameterCandidate, ...]:
    by_feature = _feature_sort_map(features)
    ordered = sorted(
        screens,
        key=lambda screen: (
            -int(screen.record["unexpected_nonzero_quartic_x_values"]),
            -int(screen.record["distinct_quartic_x_values"]),
            -by_feature[screen.candidate.parameter_u].combined_objective,
            screen.candidate.parameter_u,
        ),
    )
    selected = {
        screen.candidate.parameter_u for screen in ordered[:point_keep]
    }
    objective_order = sorted(
        screens,
        key=lambda screen: (
            -by_feature[screen.candidate.parameter_u].combined_objective,
            screen.candidate.parameter_u,
        ),
    )
    selected.update(
        screen.candidate.parameter_u for screen in objective_order[:objective_keep]
    )
    origins = sorted(
        {origin for screen in screens for origin in screen.candidate.origins}
    )
    for origin in origins:
        origin_order = [
            screen for screen in ordered if origin in screen.candidate.origins
        ]
        selected.update(
            screen.candidate.parameter_u for screen in origin_order[:origin_keep]
        )
    by_u = {screen.candidate.parameter_u: screen.candidate for screen in screens}
    return tuple(by_u[parameter_u] for parameter_u in sorted(selected))


def shifted_polynomial_coefficients(
    coefficients: Sequence[Fraction], center: Fraction
) -> tuple[Fraction, ...]:
    """Coefficients of ``Q(w+center)`` in ascending powers of ``w``."""

    center = Q(center)
    return tuple(
        sum(
            Q(coefficients[source_degree])
            * comb(source_degree, target_degree)
            * center ** (source_degree - target_degree)
            for source_degree in range(target_degree, len(coefficients))
        )
        for target_degree in range(len(coefficients))
    )


def section_centered_screen(
    candidate: ParameterCandidate,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[PointScreen, dict[str, Any]]:
    """Search all 18 classified linear-section charts for one candidate."""

    parameter_t = candidate.parameter_t
    coefficients = primitive_quartic_coefficients(RANK13_CONSTRUCTION, parameter_t)
    centers = [
        point[0] for point in rank13_known_quartic_points(parameter_t)
    ] + [
        section.point(parameter_t)[0] for section in omitted_companion_sections()
    ]
    if len(centers) != 18 or len(set(centers)) != 18:
        raise AssertionError("the classified linear-section chart set changed")
    commands: list[str] = []
    for index, center in enumerate(centers):
        shifted = shifted_polynomial_coefficients(coefficients, center)
        commands.extend(
            (
                f"Q={quartic_gp_polynomial(shifted)};",
                f"R=hyperellratpoints(Q,{height_bound});",
                f'print("ROW|{index}|",R);',
            )
        )
    commands.append("quit")
    output, wall_seconds = run_gp_capped(
        "\n".join(commands) + "\n",
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    points: list[tuple[Fraction, Fraction]] = []
    chart_signed_counts = []
    for line in output.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, index_text, vector_text = line.split("|", 2)
        index = int(index_text)
        chart_points = parse_point_vector(vector_text)
        chart_signed_counts.append(len(chart_points))
        points.extend((x_value + centers[index], z_value) for x_value, z_value in chart_points)
    if len(chart_signed_counts) != len(centers):
        raise AssertionError("PARI omitted one or more section-centred charts")
    distinct_points = tuple(dict.fromkeys(points))
    screen = _screen_record(candidate, distinct_points, height_bound=height_bound)
    return screen, {
        "chart_count": len(centers),
        "chart_naive_height_bound": height_bound,
        "signed_hits_across_charts_before_deduplication": sum(chart_signed_counts),
        "distinct_signed_points_after_global_mapping": len(distinct_points),
        "wall_seconds": wall_seconds,
    }


def merge_screens(
    left: PointScreen, right: PointScreen, *, height_bound: int
) -> PointScreen:
    if left.candidate != right.candidate:
        raise ValueError("cannot merge point screens from different candidates")
    points = tuple(dict.fromkeys((*left.raw_points, *right.raw_points)))
    return _screen_record(left.candidate, points, height_bound=height_bound)


def screen_summary(screens: Sequence[PointScreen]) -> dict[str, Any]:
    rows = [
        "|".join(
            (
                rational_to_string(screen.candidate.parameter_u),
                str(screen.record["signed_points_found"]),
                str(screen.record["distinct_quartic_x_values"]),
                str(screen.record["unexpected_nonzero_quartic_x_values"]),
                str(screen.record["unexpected_point_sha256"]),
            )
        )
        for screen in screens
    ]
    histogram = Counter(
        int(screen.record["unexpected_nonzero_quartic_x_values"])
        for screen in screens
    )
    positive = [
        dict(screen.record)
        for screen in screens
        if int(screen.record["unexpected_nonzero_quartic_x_values"]) > 0
    ]
    return {
        "candidate_count": len(screens),
        "screen_rows_sha256": hashlib.sha256("\n".join(rows).encode()).hexdigest(),
        "unexpected_x_count_histogram": {
            str(count): frequency for count, frequency in sorted(histogram.items())
        },
        "positive_yield_records": positive,
    }


def _height_record_sort_key(
    record: dict[str, Any], feature_by_u: dict[Fraction, CandidateFeatures]
) -> tuple[int, int, int, float, Fraction]:
    rank = record.get("stable_pool_numerical_rank")
    gain = record.get("stable_numerical_rank_gain")
    parameter_u = Q(record["parameter_u"])
    return (
        -(int(rank) if rank is not None else -1),
        -(int(gain) if gain is not None else -1),
        -int(record["unexpected_nonzero_quartic_x_values"]),
        -feature_by_u[parameter_u].combined_objective,
        parameter_u,
    )


def _serialize_signatures(signatures: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": signature.prime,
            "group_order": signature.group_order,
            "doubled_subgroup_order": signature.doubled_subgroup_order,
            "quotient_dimension": signature.quotient_dimension,
            "rows": [list(row) for row in signature.rows],
        }
        for signature in signatures
    ]


def exact_independence_attempt(
    parameter_u: Fraction,
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    numerical_rank: int,
    saturation_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    """Try direct mod-2 certification, then one capped saturation replay."""

    coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)
    selected = tuple(points[:numerical_rank])
    if len(selected) < 18 or any(
        not point_on_short_curve(coefficients, point) for point in selected
    ):
        return {"status": "not_attempted", "reason": "fewer than 18 exact points"}

    def reduction_attempt(basis: Sequence[tuple[Fraction, Fraction]]) -> dict[str, Any]:
        signatures = find_mod2_reduction_certificate(
            coefficients, basis, prime_bound=certificate_prime_bound
        )
        binary_rank = combined_mod2_rank(signatures, len(basis))
        answer = {
            "point_count": len(basis),
            "point_sha256": point_digest(basis),
            "certificate_primes": [signature.prime for signature in signatures],
            "combined_exact_rank_over_F2": binary_rank,
            "signatures": _serialize_signatures(signatures),
        }
        if binary_rank == len(basis):
            answer["two_torsion_certificate_prime"] = (
                find_two_torsion_certificate_prime(coefficients, prime_bound=200)
            )
            answer["certified_algebraic_rank_lower_bound"] = len(basis)
        return answer

    direct = reduction_attempt(selected)
    if direct["combined_exact_rank_over_F2"] == len(selected):
        return {"status": "certified_directly", "direct_reduction": direct}

    try:
        saturated, saturation = saturate_exact_basis(
            coefficients,
            selected,
            prime_bound=20,
            timeout=saturation_timeout,
            stack_bytes=stack_bytes,
        )
    except Exception as error:
        return {
            "status": "direct_rank_deficient_saturation_failed",
            "direct_reduction": direct,
            "saturation_error": str(error)[:1000],
        }
    saturated_reduction = reduction_attempt(saturated)
    status = (
        "certified_after_saturation"
        if saturated_reduction["combined_exact_rank_over_F2"] == len(saturated)
        else "bounded_certificate_search_rank_deficient"
    )
    return {
        "status": status,
        "direct_reduction": direct,
        "small_prime_saturation": saturation,
        "saturated_reduction": saturated_reduction,
    }


def _parse_precisions(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(values) < 2 or values != tuple(sorted(set(values))) or values[0] < 32:
        raise argparse.ArgumentTypeError("provide increasing precisions >=32")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prior-input",
        type=Path,
        action="append",
        default=None,
        help="generation-1/2 artifact; repeat twice (defaults are pinned)",
    )
    parser.add_argument(
        "--rank17-input",
        type=Path,
        default=generated / "elliptic_nagao_rank17_frontier_certificate.json",
    )
    parser.add_argument("--prefilter-combined-keep", type=int, default=2048)
    parser.add_argument("--prefilter-axis-keep", type=int, default=1024)
    parser.add_argument("--prefilter-origin-keep", type=int, default=384)
    parser.add_argument("--initial-height", type=int, default=2_000)
    parser.add_argument("--initial-point-keep", type=int, default=192)
    parser.add_argument("--initial-objective-keep", type=int, default=64)
    parser.add_argument("--initial-origin-keep", type=int, default=8)
    parser.add_argument("--rank-box-height", type=int, default=20_000)
    parser.add_argument("--height-evaluation-keep", type=int, default=96)
    parser.add_argument("--final-keep", type=int, default=12)
    parser.add_argument("--final-height", type=int, default=250_000)
    parser.add_argument("--section-chart-height", type=int, default=20_000)
    parser.add_argument("--escalation-keep", type=int, default=2)
    parser.add_argument("--escalation-height", type=int, default=1_000_000)
    parser.add_argument("--precisions", type=_parse_precisions, default=(72, 120))
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--batch-timeout", type=float, default=25.0)
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--final-timeout", type=float, default=60.0)
    parser.add_argument("--escalation-timeout", type=float, default=90.0)
    parser.add_argument("--chart-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--conductor-timeout", type=float, default=25.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank13_generation3.json",
    )
    return parser


def _default_prior_paths(root: Path) -> tuple[Path, Path]:
    generated = root / "artifacts" / "generated-results"
    return (
        generated / "elliptic_nagao_rank13_rank_gain_search.json",
        generated / "elliptic_nagao_rank13_rank_gain_mutations.json",
    )


def _validate_args(args: argparse.Namespace) -> None:
    integers = (
        args.prefilter_combined_keep,
        args.prefilter_axis_keep,
        args.prefilter_origin_keep,
        args.initial_height,
        args.initial_point_keep,
        args.initial_objective_keep,
        args.initial_origin_keep,
        args.rank_box_height,
        args.height_evaluation_keep,
        args.final_keep,
        args.final_height,
        args.section_chart_height,
        args.escalation_keep,
        args.escalation_height,
        args.batch_size,
        args.certificate_prime_bound,
    )
    if any(value <= 0 for value in integers):
        raise SystemExit("all count, height, and prime bounds must be positive")
    if not args.initial_height < args.rank_box_height < args.final_height < args.escalation_height:
        raise SystemExit("uniform search heights must be strictly increasing")
    if not args.final_keep <= args.height_evaluation_keep:
        raise SystemExit("final keep exceeds the height-evaluation count")
    timeouts = (
        args.batch_timeout,
        args.height_timeout,
        args.final_timeout,
        args.escalation_timeout,
        args.chart_timeout,
        args.saturation_timeout,
        args.conductor_timeout,
    )
    if any(timeout <= 0 or timeout > 120 for timeout in timeouts):
        raise SystemExit("every subprocess timeout must lie in (0,120]")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack bound is too small")


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    root = Path(__file__).resolve().parents[2]
    prior_paths = tuple(args.prior_input or _default_prior_paths(root))
    if len(prior_paths) != 2:
        raise SystemExit("exactly two prior population artifacts are required")
    load_exact_rank17_labels(args.rank17_input)
    prior_population = load_prior_population(prior_paths)
    trace_tables = local_trace_tables()
    residue_groups = rank17_label_residues(trace_tables)
    population = generate_generation3_population(
        residue_groups, prior_population
    )
    features = tuple(
        candidate_features(candidate, trace_tables, residue_groups)
        for candidate in population
    )
    feature_by_u = _feature_sort_map(features)
    population_text = "\n".join(
        f"{rational_to_string(item.parameter_u)}|{','.join(item.origins)}"
        for item in features
    )
    print(f"generation-3 population={len(population)}", flush=True)

    prefiltered = multiobjective_prefilter(
        features,
        combined_keep=args.prefilter_combined_keep,
        axis_keep=args.prefilter_axis_keep,
        origin_keep=args.prefilter_origin_keep,
    )
    print(f"multiobjective prefilter={len(prefiltered)}", flush=True)
    initial_screens = batch_point_screen(
        prefiltered,
        height_bound=args.initial_height,
        batch_size=args.batch_size,
        timeout_per_batch=args.batch_timeout,
        stack_bytes=args.stack_bytes,
    )
    initial_survivors = select_point_survivors_generation3(
        initial_screens,
        features,
        point_keep=min(args.initial_point_keep, len(initial_screens)),
        objective_keep=min(args.initial_objective_keep, len(initial_screens)),
        origin_keep=args.initial_origin_keep,
    )
    print(f"initial point survivors={len(initial_survivors)}", flush=True)

    rank_box_screens = batch_point_screen(
        initial_survivors,
        height_bound=args.rank_box_height,
        batch_size=args.batch_size,
        timeout_per_batch=args.batch_timeout,
        stack_bytes=args.stack_bytes,
    )
    height_candidates = select_point_survivors_generation3(
        rank_box_screens,
        features,
        point_keep=min(args.height_evaluation_keep, len(rank_box_screens)),
        objective_keep=min(16, len(rank_box_screens)),
        origin_keep=4,
    )
    rank_screen_by_u = {
        screen.candidate.parameter_u: screen for screen in rank_box_screens
    }
    height_records: list[dict[str, Any]] = []
    for index, candidate in enumerate(height_candidates, 1):
        record, _ = evaluate_height_rank(
            rank_screen_by_u[candidate.parameter_u],
            precisions=args.precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=False,
        )
        height_records.append(record)
        print(
            f"height {index}/{len(height_candidates)} u={candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']} "
            f"extra={record['unexpected_nonzero_quartic_x_values']}",
            flush=True,
        )

    height_records.sort(key=lambda record: _height_record_sort_key(record, feature_by_u))
    final_parameters = {
        Q(record["parameter_u"])
        for record in height_records[: min(args.final_keep, len(height_records))]
    }
    # Preserve one best candidate from each genuinely different geometry.
    for origin in sorted({origin for item in features for origin in item.origins}):
        match = next(
            (
                record
                for record in height_records
                if origin in feature_by_u[Q(record["parameter_u"])].origins
            ),
            None,
        )
        if match is not None:
            final_parameters.add(Q(match["parameter_u"]))
    final_candidates = tuple(
        rank_screen_by_u[parameter_u].candidate
        for parameter_u in sorted(final_parameters)
    )
    print(f"final candidates={len(final_candidates)}", flush=True)

    final_uniform_screens = batch_point_screen(
        final_candidates,
        height_bound=args.final_height,
        batch_size=1,
        timeout_per_batch=args.final_timeout,
        stack_bytes=args.stack_bytes,
    )
    final_records: list[dict[str, Any]] = []
    final_selected: dict[Fraction, tuple[tuple[Fraction, Fraction], ...]] = {}
    final_screen_by_u: dict[Fraction, PointScreen] = {}
    final_chart_metadata: dict[Fraction, dict[str, Any]] = {}
    for uniform_screen in final_uniform_screens:
        chart_screen, chart_metadata = section_centered_screen(
            uniform_screen.candidate,
            height_bound=args.section_chart_height,
            timeout=args.chart_timeout,
            stack_bytes=args.stack_bytes,
        )
        merged = merge_screens(
            uniform_screen,
            chart_screen,
            height_bound=args.final_height,
        )
        record, selected = evaluate_height_rank(
            merged,
            precisions=args.precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=True,
        )
        record["section_centered_charts"] = chart_metadata
        record["conductor_probe"] = conductor_probe(
            rank13_base_changed_short_jacobian_coefficients(
                uniform_screen.candidate.parameter_u
            ),
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        final_records.append(record)
        final_selected[uniform_screen.candidate.parameter_u] = selected
        final_screen_by_u[uniform_screen.candidate.parameter_u] = merged
        final_chart_metadata[uniform_screen.candidate.parameter_u] = chart_metadata
        print(
            f"final u={uniform_screen.candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']} "
            f"extra={record['unexpected_nonzero_quartic_x_values']} "
            f"logN={record['conductor_probe'].get('log_conductor')}",
            flush=True,
        )

    final_records.sort(key=lambda record: _height_record_sort_key(record, feature_by_u))
    escalation_parameters = tuple(
        Q(record["parameter_u"])
        for record in final_records[: min(args.escalation_keep, len(final_records))]
    )
    escalation_candidates = tuple(
        final_screen_by_u[parameter].candidate for parameter in escalation_parameters
    )
    escalation_uniform = batch_point_screen(
        escalation_candidates,
        height_bound=args.escalation_height,
        batch_size=1,
        timeout_per_batch=args.escalation_timeout,
        stack_bytes=args.stack_bytes,
    )
    final_by_u = {Q(record["parameter_u"]): record for record in final_records}
    escalation_records: list[dict[str, Any]] = []
    escalation_selected: dict[Fraction, tuple[tuple[Fraction, Fraction], ...]] = {}
    for uniform_screen in escalation_uniform:
        merged = merge_screens(
            uniform_screen,
            final_screen_by_u[uniform_screen.candidate.parameter_u],
            height_bound=args.escalation_height,
        )
        record, selected = evaluate_height_rank(
            merged,
            precisions=args.precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=True,
        )
        record["section_centered_charts"] = final_chart_metadata[
            uniform_screen.candidate.parameter_u
        ]
        record["conductor_probe"] = final_by_u[
            uniform_screen.candidate.parameter_u
        ]["conductor_probe"]
        escalation_records.append(record)
        escalation_selected[uniform_screen.candidate.parameter_u] = selected
        print(
            f"escalated u={uniform_screen.candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']} "
            f"extra={record['unexpected_nonzero_quartic_x_values']}",
            flush=True,
        )

    escalation_records.sort(
        key=lambda record: _height_record_sort_key(record, feature_by_u)
    )
    escalated_u = {Q(record["parameter_u"]) for record in escalation_records}
    frontier_records = escalation_records + [
        record for record in final_records if Q(record["parameter_u"]) not in escalated_u
    ]
    frontier_records.sort(key=lambda record: _height_record_sort_key(record, feature_by_u))

    exact_attempts = []
    for record in frontier_records:
        numerical_rank = record.get("stable_pool_numerical_rank")
        if numerical_rank is None or int(numerical_rank) < 18:
            continue
        parameter_u = Q(record["parameter_u"])
        selected = escalation_selected.get(parameter_u, final_selected[parameter_u])
        attempt = exact_independence_attempt(
            parameter_u,
            selected,
            numerical_rank=int(numerical_rank),
            saturation_timeout=args.saturation_timeout,
            stack_bytes=args.stack_bytes,
            certificate_prime_bound=args.certificate_prime_bound,
        )
        exact_attempts.append(
            {"parameter_u": rational_to_string(parameter_u), **attempt}
        )
        print(f"exact attempt u={parameter_u}: {attempt['status']}", flush=True)

    stable_ranks = [
        int(record["stable_pool_numerical_rank"])
        for record in frontier_records
        if record.get("stable_pool_numerical_rank") is not None
    ]
    maximum_rank = max(stable_ranks, default=0)
    exact_lower_bound = max(
        (
            int(
                attempt.get("direct_reduction", {}).get(
                    "certified_algebraic_rank_lower_bound", 0
                )
            )
            for attempt in exact_attempts
        ),
        default=0,
    )
    exact_lower_bound = max(
        exact_lower_bound,
        max(
            (
                int(
                    attempt.get("saturated_reduction", {}).get(
                        "certified_algebraic_rank_lower_bound", 0
                    )
                )
                for attempt in exact_attempts
            ),
            default=0,
        ),
    )

    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    input_hashes = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (*prior_paths, args.rank17_input)
    }
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded third-generation exact-point search with numerical height "
            "triage; only finite-reduction results are exact rank certificates"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "inputs": input_hashes,
        "population": {
            "count": len(population),
            "sha256": hashlib.sha256(population_text.encode()).hexdigest(),
            "minimum_reduced_denominator": min(
                candidate.parameter_u.denominator for candidate in population
            ),
            "maximum_reduced_denominator": max(
                candidate.parameter_u.denominator for candidate in population
            ),
            "disjoint_from_both_prior_populations": True,
            "prior_population_union_count": len(prior_population),
            "origin_counts": dict(
                sorted(
                    Counter(
                        origin
                        for candidate in population
                        for origin in candidate.origins
                    ).items()
                )
            ),
            "construction": {
                "pell_small_T_shell": {
                    "denominators": "65..1024",
                    "offsets_from_floor_sqrt_23550_d": "-8..8",
                },
                "rank17_local_crt": {
                    "primes": list(LABEL_PRIMES),
                    "residue_groups": {
                        str(prime): list(residues)
                        for prime, residues in residue_groups.items()
                    },
                    "class_count": len(label_crt_classes(residue_groups)),
                    "denominators": "65..512",
                },
                "rank17_tangent_mutations": {
                    "centers": [str(value) for value in EXACT_RANK17_PARAMETERS],
                    "denominators": "65..384",
                    "nonzero_numerator_deltas": "-8..8",
                },
            },
        },
        "predictor": {
            "exact_rank17_labels": [str(value) for value in EXACT_RANK17_PARAMETERS],
            "label_primes": list(LABEL_PRIMES),
            "short_prime_score_primes": list(PREDICTOR_PRIMES),
            "discriminant_power_primes": list(POWER_PRIMES),
            "combined_objective": (
                "2*rank17_residue_affinity + short_prime_score + "
                "0.35*discriminant_power_savings - 0.4*log(height(T))"
            ),
            "prefilter_count": len(prefiltered),
            "prefilter_parameter_u": [
                rational_to_string(candidate.parameter_u) for candidate in prefiltered
            ],
            "top_feature_records": [
                item.compact_json()
                for item in sorted(
                    features,
                    key=lambda item: (-item.combined_objective, item.parameter_u),
                )[:100]
            ],
        },
        "initial_point_screen": {
            "height_bound": args.initial_height,
            **screen_summary(initial_screens),
            "survivor_parameter_u": [
                rational_to_string(candidate.parameter_u)
                for candidate in initial_survivors
            ],
        },
        "rank_box_point_screen": {
            "height_bound": args.rank_box_height,
            **screen_summary(rank_box_screens),
            "height_evaluated_parameter_u": [
                rational_to_string(candidate.parameter_u)
                for candidate in height_candidates
            ],
        },
        "rank_box_height_evaluation": {
            "precisions": list(args.precisions),
            "records": height_records,
        },
        "final_box": {
            "uniform_height_bound": args.final_height,
            "section_centered_chart_height_bound": args.section_chart_height,
            "records": final_records,
        },
        "escalation_box": {
            "uniform_height_bound": args.escalation_height,
            "selection_count": args.escalation_keep,
            "records": escalation_records,
        },
        "exact_independence_attempts": exact_attempts,
        "summary": {
            "maximum_stable_numerical_rank": maximum_rank,
            "maximum_exact_rank_lower_bound_newly_certified": exact_lower_bound,
            "stable_rank_at_least_18_lead_count": sum(rank >= 18 for rank in stable_ranks),
            "exact_rank_at_least_18_certificate_count": sum(
                attempt["status"] in ("certified_directly", "certified_after_saturation")
                for attempt in exact_attempts
            ),
            "frontier": frontier_records,
            "target_hit": any(
                int(
                    attempt.get("direct_reduction", {}).get(
                        "certified_algebraic_rank_lower_bound", 0
                    )
                )
                >= 21
                or int(
                    attempt.get("saturated_reduction", {}).get(
                        "certified_algebraic_rank_lower_bound", 0
                    )
                )
                >= 21
                for attempt in exact_attempts
            ),
            "interpretation": (
                "negative bounded searches do not bound Mordell-Weil rank; "
                "stable height ranks remain numerical evidence"
            ),
        },
        "parameters": {
            "prefilter_combined_keep": args.prefilter_combined_keep,
            "prefilter_axis_keep": args.prefilter_axis_keep,
            "prefilter_origin_keep": args.prefilter_origin_keep,
            "initial_point_keep": args.initial_point_keep,
            "initial_objective_keep": args.initial_objective_keep,
            "initial_origin_keep": args.initial_origin_keep,
            "height_evaluation_keep": args.height_evaluation_keep,
            "final_keep": args.final_keep,
            "batch_size": args.batch_size,
            "batch_timeout_seconds": args.batch_timeout,
            "height_timeout_seconds": args.height_timeout,
            "final_timeout_seconds": args.final_timeout,
            "escalation_timeout_seconds": args.escalation_timeout,
            "chart_timeout_seconds": args.chart_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "conductor_timeout_seconds": args.conductor_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}", flush=True)
    print(
        f"maximum stable numerical rank={maximum_rank}; "
        f"new exact lower bound={exact_lower_bound}",
        flush=True,
    )


if __name__ == "__main__":
    main()
