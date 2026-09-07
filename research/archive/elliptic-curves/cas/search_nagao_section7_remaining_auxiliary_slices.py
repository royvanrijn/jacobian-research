#!/usr/bin/env python3
"""Search the eleven remaining minimum-intercept Section-7 genus-one slices.

The pinned H=200000 accidental-slice population supplies rational points on
sixteen auxiliary genus-one quartics through T0=5081/47.  Five positive-point-
cardinality or already-closed slices are owned by other search lanes.  This
script handles exactly the complementary eleven slices.

For each slice, the pointed-quartic map sends the chosen T0 point to the
auxiliary identity.  Thus mapped H=200000 points are differences from T0.  We
double those differences, select a numerically height-independent subset at
two precisions, and exhaust all nonzero {-1,0,1} coefficient vectors, with the
declared cap 3^10=59049 per slice.  Every inverse image is checked exactly on
the original quartic.  The pinned H population, all generic-section images,
singular fibres, and duplicates are excluded before conductor-proxy work.

This is a bounded constructive computation.  A stable numerical auxiliary
height rank is not an exact Mordell-Weil rank certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import quartic_value, short_jacobian_coefficients
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_ROOTS,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank20_t5081_neighborhood import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from search_nagao_rank21_accidental_slices import (
    PINNED_IDENTITY_HEIGHT_200000_PARAMETERS,
    Slice,
    build_slices,
    replay_pinned_identity_parameters,
    select_minimum_intercept_priority_slices,
)
from search_nagao_section7_auxiliary_jacobians import (
    AuxiliarySlice,
    SECTION7_QUARTIC_COEFFICIENT_POLYNOMIALS,
    SliceSpecification,
    T0,
    canonical_parameter,
    make_auxiliary_slice,
    polynomial_value,
    projective_height,
    weierstrass_add,
    weierstrass_multiply,
    weierstrass_negate,
)
from triage_nagao_rank13_finalists import height_matrix_replay, stable_height_rank


Q = Fraction
PINNED_ACCIDENTAL_ARTIFACT_SHA256 = (
    "125a6b0df7941099547039302b6f1878b5009dcde774328527952699877b1670"
)
OPEN_SLICE_IDS = (
    "a02_sp01",
    "a03_sp01",
    "a05_sp01",
    "a06_sp01",
    "a07_sp01",
    "a11_sm01",
    "a12_sp01",
    "a13_sp01",
    "a14_sp01",
    "a15_sp01",
    "a16_sp01",
)
COORDINATION_EXCLUDED_SLICE_IDS = (
    "a01_sp01",
    "a04_sp01",
    "a08_sm01",
    "a09_sm01",
    "a10_sp01",
)
MAX_BASIS_DIMENSION = 10
MAX_VECTOR_POPULATION_INCLUDING_ZERO = 59_049
PROXY_THRESHOLD = 190.0
PROXY_TRIAL_BOUND = 2_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")
MAX_EXACT_CONDUCTORS = 256
DEFAULT_INPUT = Path(
    "artifacts/generated-results/elliptic_nagao_rank21_accidental_slices.json"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_section7_remaining_auxiliary_slices.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_remaining_auxiliary_slices.py"
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
    }


def parameter_digest(parameters: Iterable[Fraction]) -> str:
    payload = "\n".join(
        rational_to_string(value)
        for value in sorted(
            {Q(value) for value in parameters},
            key=lambda value: (
                projective_height(value),
                value.numerator,
                value.denominator,
            ),
        )
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def ternary_vectors(dimension: int) -> tuple[tuple[int, ...], ...]:
    if not 0 <= dimension <= MAX_BASIS_DIMENSION:
        raise ValueError("the ternary search dimension is outside [0,10]")
    if 3**dimension > MAX_VECTOR_POPULATION_INCLUDING_ZERO:
        raise ValueError("the declared ternary vector cap was exceeded")
    return tuple(
        vector
        for vector in itertools.product((-1, 0, 1), repeat=dimension)
        if any(vector)
    )


@dataclass(frozen=True)
class PreparedSlice:
    source: Slice
    auxiliary: AuxiliarySlice
    known_points: tuple[tuple[Fraction, Fraction], ...]


def load_open_slices(
    input_path: Path,
) -> tuple[tuple[PreparedSlice, ...], dict[str, Any]]:
    source_hash = file_sha256(input_path)
    if source_hash != PINNED_ACCIDENTAL_ARTIFACT_SHA256:
        raise ValueError(
            "the pinned accidental-slice artifact changed; inspect and repin it"
        )
    artifact = json.loads(input_path.read_text(encoding="utf-8"))
    accidental = tuple(
        (Q(record["x"]), Q(record["y"]))
        for record in artifact["decontamination_at_T0"]["accidental_points"]
    )
    priority = select_minimum_intercept_priority_slices(build_slices(accidental))
    replayed, replay_metadata = replay_pinned_identity_parameters(priority)
    by_identifier = {
        item.identifier: (item, points)
        for item, points in zip(priority, replayed, strict=True)
    }
    expected = set(OPEN_SLICE_IDS) | set(COORDINATION_EXCLUDED_SLICE_IDS)
    if set(by_identifier) != expected:
        raise AssertionError("the sixteen priority-slice identifiers changed")

    prepared = []
    association_counts = {}
    for identifier in OPEN_SLICE_IDS:
        item, normalized_points = by_identifier[identifier]
        specification = SliceSpecification(
            label=identifier,
            slope=Q(item.slope),
            intercept=item.intercept,
            base_x=item.source_point[0],
            base_y=item.source_point[1],
            parent_accidental_label=identifier,
        )
        auxiliary = make_auxiliary_slice(specification)
        known_points = []
        for parameter, normalized_ordinate in sorted(normalized_points):
            original_ordinate = item.normalized.original_ordinate(
                parameter, normalized_ordinate
            )
            if original_ordinate**2 != auxiliary.quartic_value(parameter):
                raise AssertionError("a pinned normalized point lost its raw ordinate")
            known_points.append((parameter, original_ordinate))
        prepared.append(PreparedSlice(item, auxiliary, tuple(known_points)))
        association_counts[identifier] = {
            "parameters_including_T0": len({point[0] for point in known_points}),
            "signed_affine_points_including_T0": len(known_points),
        }
    return tuple(prepared), {
        "input_sha256": source_hash,
        "pinned_identity_height": replay_metadata["height"],
        "pinned_global_parameter_count_excluding_T0": replay_metadata[
            "signless_parameter_count_excluding_T0"
        ],
        "open_slice_association_counts": association_counts,
        "all_associations_replayed_by_exact_square_tests": True,
    }


def doubled_difference_population(
    prepared: PreparedSlice,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
]:
    """Map P-T0 and return distinct doubled images with their exact sources."""

    auxiliary = prepared.auxiliary
    unique: dict[tuple[Fraction, Fraction], list[dict[str, str]]] = {}
    for parameter, ordinate in prepared.known_points:
        if parameter == T0:
            continue
        mapped_difference = auxiliary.forward((parameter, ordinate))
        doubled = weierstrass_multiply(
            auxiliary.weierstrass_coefficients, mapped_difference, 2
        )
        if doubled is None:
            continue
        unique.setdefault(doubled, []).append(
            {
                "T": rational_to_string(parameter),
                "quartic_ordinate": rational_to_string(ordinate),
            }
        )
    points = tuple(unique)
    records = tuple(
        {
            "doubled_point": point_record(point),
            "source_known_points": unique[point],
        }
        for point in points
    )
    return points, records


def select_basis(
    prepared: PreparedSlice,
    doubled_points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
    tuple[int, ...],
]:
    if not doubled_points:
        return (), (), ()
    runs = height_matrix_replay(
        prepared.auxiliary.weierstrass_coefficients,
        doubled_points,
        precisions=(72, 120),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    rank = stable_height_rank(runs)
    indices = tuple(int(value) for value in runs[-1]["subset_indices_one_based"])
    if rank != len(indices):
        raise AssertionError("the stable height rank and subset length differ")
    selected_indices = indices[:MAX_BASIS_DIMENSION]
    selected = tuple(doubled_points[index - 1] for index in selected_indices)
    return selected, runs, selected_indices


def prior_parameters() -> set[Fraction]:
    return {
        canonical_parameter(value)
        for value in (*PINNED_IDENTITY_HEIGHT_200000_PARAMETERS, T0)
    }


def explicit_section7_quartic(parameter: Fraction) -> tuple[Fraction, ...]:
    """Evaluate the repository's explicit e,d,c,b,a coefficient polynomials."""

    return tuple(
        polynomial_value(coefficients, Q(parameter))
        for coefficients in SECTION7_QUARTIC_COEFFICIENT_POLYNOMIALS
    )


def generic_abscissae(parameter: Fraction) -> tuple[tuple[str, Fraction], ...]:
    parameter = Q(parameter)
    records = [
        (f"visible-{index:02d}", Q(root) + slope * parameter)
        for index, (root, slope) in enumerate(
            itertools.product(sorted(Q(root) for root in SECTION7_ROOTS), (Q(-1), Q(1)))
        )
    ]
    records.extend(
        (
            f"linear-{section.label}",
            section.slope * parameter + section.intercept,
        )
        for section in SECTION7_LINEAR_COMPANION_SECTIONS
    )
    records.extend(
        (
            f"quadratic-{section.label}",
            section.quadratic_coefficient * parameter**2
            + section.linear_coefficient * parameter
            + section.constant_coefficient,
        )
        for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
    )
    if len(records) != 21:
        raise AssertionError("the generic abscissa count changed")
    return tuple(records)


def quartic_g_coefficients(
    coefficients: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    e_value, d_value, c_value, b_value, a_value = (
        Q(value) for value in coefficients
    )
    return (
        d_value**2 / 16 - c_value * e_value / 6,
        c_value * d_value / 12 - b_value * e_value / 2,
        c_value**2 / 12 - b_value * d_value / 8 - a_value * e_value,
        b_value * c_value / 12 - a_value * d_value / 2,
        b_value**2 / 16 - a_value * c_value / 6,
    )


def generic_labels_for_point(
    parameter: Fraction,
    point: tuple[Fraction, Fraction],
) -> tuple[str, ...]:
    """Compare exact quartic x and covariant-Jacobian X with all 21 sections."""

    parameter = Q(parameter)
    x_value, ordinate = (Q(value) for value in point)
    coefficients = explicit_section7_quartic(parameter)
    if ordinate**2 != quartic_value(coefficients, x_value):
        raise AssertionError("the generated point missed the explicit quartic")
    g_coefficients = quartic_g_coefficients(coefficients)
    generated_jacobian_x = (
        36 * polynomial_value(g_coefficients, x_value) / ordinate**2
        if ordinate
        else None
    )
    labels = []
    for label, generic_x in generic_abscissae(parameter):
        generic_y_squared = quartic_value(coefficients, generic_x)
        if generic_x == x_value:
            labels.append(f"quartic-x:{label}")
        if (
            generated_jacobian_x is not None
            and generic_y_squared
            and 36 * polynomial_value(g_coefficients, generic_x)
            / generic_y_squared
            == generated_jacobian_x
        ):
            labels.append(f"jacobian-sign-pair:{label}")
    if generated_jacobian_x is None:
        labels.append("ramification-ordinate-zero")
    return tuple(sorted(set(labels)))


def ternary_group_elements(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[int, ...], tuple[Fraction, Fraction] | None], ...]:
    """Build the ternary cube recursively, reusing every partial group sum."""

    states: list[tuple[tuple[int, ...], tuple[Fraction, Fraction] | None]] = [
        ((), None)
    ]
    for basis_point in basis:
        negative = weierstrass_negate(coefficients, basis_point)
        next_states = []
        for vector, point in states:
            next_states.extend(
                (
                    (
                        vector + (-1,),
                        weierstrass_add(coefficients, point, negative),
                    ),
                    (vector + (0,), point),
                    (
                        vector + (1,),
                        weierstrass_add(coefficients, point, basis_point),
                    ),
                )
            )
        states = next_states
    return tuple((vector, point) for vector, point in states if any(vector))


def generate_slice_candidates(
    prepared: PreparedSlice,
    basis: Sequence[tuple[Fraction, Fraction]],
) -> tuple[dict[Fraction, list[dict[str, Any]]], dict[str, int]]:
    auxiliary = prepared.auxiliary
    prior = prior_parameters()
    candidates: dict[Fraction, list[dict[str, Any]]] = {}
    counts = {
        "coefficient_vectors": 0,
        "identity_or_torsion_relations": 0,
        "exceptional_inverse_images": 0,
        "zero_parameters": 0,
        "prior_H200000_or_T0_parameters": 0,
        "generic_section_images": 0,
        "singular_section7_parameters": 0,
        "new_parameters_at_height_at_most_200000": 0,
        "accepted_source_images": 0,
    }
    elements = ternary_group_elements(auxiliary.weierstrass_coefficients, basis)
    if len(elements) != len(ternary_vectors(len(basis))):
        raise AssertionError("the recursive ternary population changed size")
    for vector, point in elements:
        counts["coefficient_vectors"] += 1
        if point is None:
            counts["identity_or_torsion_relations"] += 1
            continue
        inverse = auxiliary.inverse(point)
        if inverse is None:
            counts["exceptional_inverse_images"] += 1
            continue
        signed_parameter, ordinate = inverse
        if ordinate**2 != auxiliary.quartic_value(signed_parameter):
            raise AssertionError("a generated inverse missed the raw slice quartic")
        parameter = canonical_parameter(signed_parameter)
        if parameter == 0:
            counts["zero_parameters"] += 1
            continue
        if parameter in prior:
            counts["prior_H200000_or_T0_parameters"] += 1
            continue
        x_value = (
            prepared.source.slope * signed_parameter
            + prepared.source.intercept
        )
        quartic_point = x_value, ordinate
        labels = generic_labels_for_point(signed_parameter, quartic_point)
        if labels:
            counts["generic_section_images"] += 1
            continue
        try:
            homogenized_discriminant(parameter)
        except ValueError:
            counts["singular_section7_parameters"] += 1
            continue
        if projective_height(parameter) <= 200_000:
            counts["new_parameters_at_height_at_most_200000"] += 1
        counts["accepted_source_images"] += 1
        candidates.setdefault(parameter, []).append(
            {
                "slice": prepared.source.identifier,
                "coefficient_vector": list(vector),
                "signed_T": rational_to_string(signed_parameter),
                "forced_quartic_point": point_record(quartic_point),
            }
        )
    counts["unique_accepted_parameters"] = len(candidates)
    return candidates, counts


def proxy_records(
    candidates: dict[Fraction, list[dict[str, Any]]],
    *,
    trial_bound: int,
) -> tuple[dict[str, Any], ...]:
    records = []
    for parameter, sources in candidates.items():
        proxy = conductor_radical_proxy(
            parameter, trial_prime_bound=trial_bound
        )
        records.append(
            {
                "constructor_parameter_T": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "radical_proxy": proxy,
                "source_count": len(sources),
                "sources": sources,
            }
        )
    return tuple(
        sorted(
            records,
            key=lambda record: (
                record["radical_proxy"]["log_radical_upper_proxy"],
                record["projective_height"],
                Q(record["constructor_parameter_T"]),
            ),
        )
    )


def exact_conductors(
    records: Sequence[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    selected = tuple(
        record
        for record in records
        if record["radical_proxy"]["log_radical_upper_proxy"] < PROXY_THRESHOLD
    )
    if len(selected) > MAX_EXACT_CONDUCTORS:
        raise RuntimeError("the proxy<190 population exceeded its exact-conductor cap")
    completed = []
    failures = []
    for index, record in enumerate(selected, start=1):
        parameter = Q(record["constructor_parameter_T"])
        try:
            result = minimal_curve_data(
                short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError) as error:
            failures.append(
                {
                    "constructor_parameter_T": rational_to_string(parameter),
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "one_attempt_no_retry": True,
                    "error": str(error)[:500],
                }
            )
            continue
        completed.append(
            {
                **record,
                "conductor": str(result["conductor"]),
                "log_conductor": result["log_conductor"],
                "root_number": result["root_number"],
                "minimal_discriminant": str(result["minimal_discriminant"]),
                "minimal_model": [str(value) for value in result["minimal_model"]],
                "below_strict_182_72_target": (
                    Decimal(result["log_conductor"]) < TARGET_LOG_CONDUCTOR
                ),
                "status": "completed",
            }
        )
        print(
            f"conductor {index}/{len(selected)} T={parameter} "
            f"lnN={result['log_conductor']}",
            flush=True,
        )
    return tuple(completed), tuple(failures)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--proxy-trial-bound", type=int, default=PROXY_TRIAL_BOUND)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    # Exact group-law parameters can have more than 4300 decimal digits.  The
    # proxy routine records the exact cofactor digit count, so Python's display
    # guard must not truncate or reject those integers.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    if not 0 < args.height_timeout <= 60:
        raise SystemExit("--height-timeout must be in (0,60]")
    if not 0 < args.conductor_timeout <= 60:
        raise SystemExit("--conductor-timeout must be in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    if not 2 <= args.proxy_trial_bound <= 10_000:
        raise SystemExit("--proxy-trial-bound must be in [2,10000]")

    prepared_slices, source_record = load_open_slices(args.input)
    all_candidates: dict[Fraction, list[dict[str, Any]]] = {}
    slice_records = []
    total_vectors = 0
    for prepared in prepared_slices:
        doubled_points, doubled_records = doubled_difference_population(prepared)
        basis, height_runs, selected_indices = select_basis(
            prepared,
            doubled_points,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        vectors = ternary_vectors(len(basis))
        candidates, counts = generate_slice_candidates(prepared, basis)
        total_vectors += len(vectors)
        for parameter, sources in candidates.items():
            all_candidates.setdefault(parameter, []).extend(sources)
        slice_records.append(
            {
                "id": prepared.source.identifier,
                "slope": prepared.source.slope,
                "intercept": rational_to_string(prepared.source.intercept),
                "source_T0_point": point_record(prepared.source.source_point),
                "known_H200000_signed_point_count_including_T0": len(
                    prepared.known_points
                ),
                "known_parameters_including_T0": sorted(
                    {
                        rational_to_string(point[0])
                        for point in prepared.known_points
                    }
                ),
                "mapped_distinct_doubled_difference_count": len(doubled_points),
                "doubled_difference_population": list(doubled_records),
                "height_replay": list(height_runs),
                "stable_numerical_height_rank": len(selected_indices),
                "selected_indices_one_based": list(selected_indices),
                "selected_doubled_basis": [point_record(point) for point in basis],
                "ternary_vector_count_excluding_zero": len(vectors),
                "generation_counts": counts,
            }
        )
        print(
            f"slice {prepared.source.identifier} directions={len(basis)} "
            f"vectors={len(vectors)} unique_T={len(candidates)}",
            flush=True,
        )

    proxies = proxy_records(
        all_candidates, trial_bound=args.proxy_trial_bound
    )
    proxy_below = tuple(
        record
        for record in proxies
        if record["radical_proxy"]["log_radical_upper_proxy"] < PROXY_THRESHOLD
    )
    print(
        "PROXY_CHECKPOINT "
        + json.dumps(
            {
                "unique_parameters": len(all_candidates),
                "candidate_sha256": parameter_digest(all_candidates),
                "proxy_below_190": len(proxy_below),
                "minimum_proxy": (
                    proxies[0]["radical_proxy"]["log_radical_upper_proxy"]
                    if proxies
                    else None
                ),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    conductor_records, conductor_failures = exact_conductors(
        proxies,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    subtarget = tuple(
        record
        for record in conductor_records
        if record["below_strict_182_72_target"]
    )
    artifact = {
        "schema_version": 1,
        "status": "bounded remaining auxiliary-slice search complete",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "scope": {
            "open_slice_ids": list(OPEN_SLICE_IDS),
            "coordination_excluded_slice_ids": list(
                COORDINATION_EXCLUDED_SLICE_IDS
            ),
            "reason_for_coordination_exclusion": (
                "a01,a04,a08,a09 are owned by the concurrent ellfromeqn lane; "
                "a10 was already closed by the root lane"
            ),
            "no_other_slices_searched": True,
        },
        "source_population": source_record,
        "method": {
            "constructor_T0": rational_to_string(T0),
            "pointed_map_identity": (
                "the chosen T0 point maps to O, so every mapped known point "
                "is its exact difference from T0"
            ),
            "searched_subgroup": "twice the mapped known-point difference lattice",
            "height_precisions": [72, 120],
            "height_independence_status": (
                "stable numerical direction selection, not an exact rank proof"
            ),
            "coefficient_alphabet": [-1, 0, 1],
            "maximum_basis_dimension": MAX_BASIS_DIMENSION,
            "maximum_vectors_including_zero_per_slice": (
                MAX_VECTOR_POPULATION_INCLUDING_ZERO
            ),
            "actual_total_nonzero_vectors": total_vectors,
            "inverse_map": (
                "exact explicit pointed-quartic birational inverse; no mapX "
                "factor ambiguity remains"
            ),
            "prior_parameter_exclusion_count_after_T_to_abs_T": len(
                prior_parameters()
            ),
            "generic_filter": (
                "exact comparison with all 21 generic quartic abscissas and "
                "their short-Jacobian sign pairs at the signed parameter"
            ),
        },
        "slices": slice_records,
        "generation": {
            "unique_accepted_parameters": len(all_candidates),
            "candidate_parameter_sha256": parameter_digest(all_candidates),
            "all_candidate_records": list(proxies),
        },
        "proxy_filter": {
            "trial_prime_bound": args.proxy_trial_bound,
            "strict_log_radical_upper_proxy_threshold": PROXY_THRESHOLD,
            "population_count": len(proxies),
            "below_threshold_count": len(proxy_below),
            "below_threshold_records": list(proxy_below),
            "minimum_log_radical_upper_proxy": (
                proxies[0]["radical_proxy"]["log_radical_upper_proxy"]
                if proxies
                else None
            ),
        },
        "exact_conductors": {
            "attempted": len(proxy_below),
            "completed": len(conductor_records),
            "failures": list(conductor_failures),
            "records": list(conductor_records),
            "sub_182_72_count": len(subtarget),
            "sub_182_72_records": list(subtarget),
        },
        "outcome": {
            "subtarget_parameters_for_H50000_H250000_triage": [
                record["constructor_parameter_T"] for record in subtarget
            ],
            "rank21_certified": False,
            "breakthrough_curve_found": False,
            "negative_result_scope": (
                "only the declared doubled sublattices and ternary coefficient "
                "boxes on the eleven named slices"
            ),
        },
        "bounded_execution": {
            "height_timeout_seconds_each": args.height_timeout,
            "conductor_timeout_seconds_each": args.conductor_timeout,
            "maximum_exact_conductors": MAX_EXACT_CONDUCTORS,
            "one_attempt_no_retry": True,
            "no_detached_processes": True,
            "python_integer_string_digit_limit_disabled_for_exact_artifact": True,
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {args.output}: candidates={len(all_candidates)} "
        f"proxy<190={len(proxy_below)} subtarget={len(subtarget)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
