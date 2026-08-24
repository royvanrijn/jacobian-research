#!/usr/bin/env python3
"""Complete the exact max-root-100 Mestre T=1,...,8 fiber screen.

The frozen max-root-100 census contains 235 generically nonsingular,
nonreflection six-root Mestre families, but its rank-blind specialization
tranche used only 40 of the 191 families beyond diameter 50.  This standalone
continuation closes the entire 235 by 8 integer panel.  Exact admissibility is
recomputed for every proposed fiber.  Earlier conductor, visible-height, and
H=5000 results are reused only at their respective, separately hashed phase
boundaries; every genuinely missing phase is run once.

For every admissible fiber the script checks the twelve displayed quartic
points and their exact Jacobian images.  A two-precision canonical-height
matrix chooses a numerical independent subset, which is then certified by
finite reductions whenever the fixed prime bound suffices.  Exact PARI
minimal-model/conductor computations precede the complete missing H=5000
point-search phase.  After that population is closed, a deterministic tranche
of the strongest signals is escalated through H=50000, H=250000, and
H=1000000.  Every numerical signal of rank at least 14 immediately receives
an exact finite-reduction attempt; numerical ranks themselves are never
reported as proofs.

All external work uses capped foreground process groups.  A timeout is a final
record for that declared call: this driver has no retry path and leaves no
worker or enumerator process behind.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable, Sequence

from mestre_root_tuples import SixRootMestreConstruction
from extend_nagao_u42_frontier import saturate_exact_basis
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    TARGET_LOG_CONDUCTOR,
    bounded_quartic_points,
    canonical_signless_points,
    capped_minimal_curve_data,
    classify_nonreflection,
    finite_reduction_attempt,
    height_matrix_replay,
    numerical_subset,
    point_digest,
    point_record,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    sha256_file,
    tuple_digest,
    verify_enumerator_records,
)
from search_mestre_root_tuple_scale_max100 import (
    EXPECTED_MAX100_COUNTS,
    EXPECTED_MAX100_NONSINGULAR_SHA256,
    STACK_BYTES,
    compiled_enumeration_max100,
    stable_json_digest,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

PARAMETERS = tuple(range(1, 9))
VISIBLE_CERTIFICATE_PRIME_BOUND = 500
SIGNAL_CERTIFICATE_PRIME_BOUND = 2_000
H5000_HEIGHT = 5_000
H5000_MAPPING_CAP = 128
ESCALATION_STAGES = (
    ("H50000", 50_000, 16, 256),
    ("H250000", 250_000, 4, 384),
    ("H1000000", 1_000_000, 1, 512),
)

FROZEN = {
    "compiled_source_sha256": (
        "31650333800698201819eddc91bf228089824bca026c629c9360683324a69eb5"
    ),
    "max50_driver_sha256": (
        "5e7228b95ae995019fbc50b9f7667de41e06a86b4490f0feacff5702bb5cc174"
    ),
    "max100_driver_sha256": (
        "34677c38be30aa15e99b3239a6d487a51c158fa33326826d37ceead310555600"
    ),
    "max50_artifact_sha256": (
        "fd2dccb1fd08aad70857df7ca19df77bd521e2be017b98f5579a748fd26cfc14"
    ),
    "max100_artifact_sha256": (
        "63dcd39555ad8b39c7b584a16663164bf73e6c6c59906b6a230bfa9b9f65a3bb"
    ),
}


def fiber_identifier(roots: Sequence[int], parameter: int) -> str:
    return "r" + "_".join(map(str, roots)) + f"_t{parameter}"


def identifier_digest(identifiers: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(sorted(identifiers)).encode()).hexdigest()


def compact_height_record(
    status: str,
    *,
    quartic_points: Sequence[tuple[Fraction, Fraction]],
    jacobian_points: Sequence[tuple[Fraction, Fraction]],
    height: Sequence[dict[str, Any]],
    subset: Sequence[tuple[Fraction, Fraction]],
    certificate: dict[str, Any],
) -> dict[str, Any]:
    return {
        "status": status,
        "quartic_point_count": len(quartic_points),
        "distinct_quartic_abscissas": len({point[0] for point in quartic_points}),
        "quartic_point_sha256": point_digest(quartic_points),
        "jacobian_point_count": len(jacobian_points),
        "jacobian_point_sha256": point_digest(jacobian_points),
        "exact_quartic_and_jacobian_membership_checked": True,
        "height_matrix_runs": list(height),
        "stable_numerical_rank": int(height[-1]["numerical_rank"]),
        "numerical_subset": [point_record(point) for point in subset],
        "finite_reduction_of_numerical_subset": certificate,
        "numerical_rank_is_not_an_independence_certificate": True,
    }


def compute_visible_phase(
    construction: SixRootMestreConstruction,
    parameter: int,
    *,
    height_timeout: float,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...] | None]:
    parameter_q = Q(parameter)
    coefficients = construction.primitive_jacobian_coefficients(parameter_q)
    try:
        quartic_points = primitive_visible_points(construction, parameter_q)
        jacobian_points = tuple(
            quartic_point_to_jacobian(construction, parameter_q, point)
            for point in quartic_points
        )
        if (
            len(quartic_points) != 12
            or len({point[0] for point in quartic_points}) != 12
            or len(jacobian_points) != 12
        ):
            raise AssertionError("an admissible fiber lost a displayed point")
        height = height_matrix_replay(
            coefficients,
            jacobian_points,
            precisions=(72, 120),
            timeout=height_timeout,
            stack_bytes=STACK_BYTES,
        )
        subset = numerical_subset(jacobian_points, height)
        certificate = finite_reduction_attempt(
            coefficients,
            subset,
            prime_bound=VISIBLE_CERTIFICATE_PRIME_BOUND,
        )
        if certificate["status"] == "certified" and certificate.get(
            "two_torsion_certificate_prime"
        ) is None:
            raise AssertionError("a visible certificate omitted its no-2-torsion prime")
        return (
            compact_height_record(
                "completed exact visible membership, height triage, and finite reductions",
                quartic_points=quartic_points,
                jacobian_points=jacobian_points,
                height=height,
                subset=subset,
                certificate=certificate,
            ),
            subset,
        )
    except CappedProcessTimeout:
        return (
            {
                "status": "timeout-no-retry",
                "timeout_seconds": height_timeout,
            },
            None,
        )
    except Exception as error:
        return ({"status": "error-no-retry", "error": str(error)[:1000]}, None)


def imported_visible_record(
    record: dict[str, Any], coefficients: Sequence[Fraction]
) -> dict[str, Any]:
    phase = record["visible_point_triage"]
    if not phase["status"].startswith("completed"):
        raise AssertionError("a frozen visible exclusion was not complete")
    subset = tuple(
        (Q(point["x"]), Q(point["y"])) for point in phase["numerical_subset"]
    )
    certificate = finite_reduction_attempt(
        coefficients,
        subset,
        prime_bound=VISIBLE_CERTIFICATE_PRIME_BOUND,
    )
    if certificate["status"] == "certified" and certificate.get(
        "two_torsion_certificate_prime"
    ) is None:
        raise AssertionError("an imported visible certificate omitted no 2-torsion")
    return {
        "status": "reused complete frozen max-root-50 visible phase",
        "source_status": phase["status"],
        "quartic_point_count": phase["quartic_point_count"],
        "distinct_quartic_abscissas": phase["distinct_quartic_abscissas"],
        "quartic_point_sha256": phase["quartic_point_sha256"],
        "jacobian_point_count": phase["jacobian_point_count"],
        "jacobian_point_sha256": phase["jacobian_point_sha256"],
        "exact_quartic_and_jacobian_membership_checked": phase[
            "exact_quartic_and_jacobian_membership_checked"
        ],
        "height_matrix_runs": phase["height_matrix_runs"],
        "stable_numerical_rank": phase["stable_numerical_rank"],
        "numerical_subset": phase["numerical_subset"],
        "finite_reduction_of_numerical_subset": certificate,
        "numerical_rank_is_not_an_independence_certificate": True,
    }


def normalized_conductor(
    phase: dict[str, Any], *, source: str
) -> dict[str, Any]:
    status = phase["status"]
    if not status.startswith("completed"):
        raise AssertionError("a frozen conductor exclusion was not complete")
    return {
        "status": f"reused complete {source} conductor phase",
        "source_status": status,
        "minimal_model": phase["minimal_model"],
        "conductor": phase["conductor"],
        "log_conductor": phase["log_conductor"],
        "minimal_discriminant": phase["minimal_discriminant"],
        "root_number": phase["root_number"],
        "below_strict_log_conductor_target_numerically": phase[
            "below_strict_log_conductor_target_numerically"
        ],
    }


def normalized_h5000_from_max50(record: dict[str, Any]) -> dict[str, Any]:
    search = record["point_search"]
    triage = record["augmented_rank_triage"]
    if not search["status"].startswith("complete"):
        raise AssertionError("a frozen max-root-50 H5000 exclusion was incomplete")
    return {
        "status": "reused complete frozen max-root-50 H5000 point/height phase",
        "source_status": search["status"],
        "height_bound": search["height_bound"],
        "signed_points_returned": search["signed_points_returned"],
        "distinct_nonzero_ordinate_abscissas": search[
            "distinct_nonzero_ordinate_abscissas"
        ],
        "abscissas_retained_for_mapping": search["abscissas_retained_for_mapping"],
        "mapping_cap": search["mapping_cap"],
        "mapping_truncated": search["mapping_truncated"],
        "pool_point_count_modulo_inverse": triage["pool_point_count_modulo_inverse"],
        "pool_point_sha256": triage["pool_point_sha256"],
        "height_matrix_runs": triage["height_matrix_runs"],
        "stable_numerical_rank": triage["stable_numerical_rank"],
        "numerical_subset": triage["numerical_subset"],
        "numerical_rank_is_not_an_independence_certificate": True,
    }


def normalized_h5000_from_max100(record: dict[str, Any]) -> dict[str, Any]:
    phase = record["point_triage"]
    if not phase["status"].startswith("completed"):
        raise AssertionError("a frozen max-root-100 H5000 exclusion was incomplete")
    result = dict(phase)
    result["source_status"] = result["status"]
    result["status"] = "reused complete frozen max-root-100 H5000 point/height phase"
    return result


def bounded_rank_search(
    construction: SixRootMestreConstruction,
    parameter: int,
    *,
    height_bound: int,
    mapping_cap: int,
    point_timeout: float,
    height_timeout: float,
    saturation_timeout: float,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...] | None]:
    parameter_q = Q(parameter)
    coefficients = construction.primitive_jacobian_coefficients(parameter_q)
    visible_quartic = primitive_visible_points(construction, parameter_q)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in visible_quartic
    )
    try:
        quartic_coefficients = construction.primitive_quartic_coefficients(parameter_q)
        raw = bounded_quartic_points(
            quartic_coefficients,
            height_bound=height_bound,
            timeout=point_timeout,
            stack_bytes=STACK_BYTES,
        )
        signless = canonical_signless_points(raw)
        retained = signless[:mapping_cap]
        if any(
            point[1] ** 2 != quartic_value(quartic_coefficients, point[0])
            for point in retained
        ):
            raise AssertionError("a bounded search returned a point off the quartic")
        searched_jacobian = tuple(
            quartic_point_to_jacobian(construction, parameter_q, point)
            for point in retained
        )
        pool_by_x = {point[0]: point for point in visible_jacobian}
        for point in searched_jacobian:
            pool_by_x.setdefault(point[0], point)
        pool = tuple(pool_by_x.values())
        height = height_matrix_replay(
            coefficients,
            pool,
            precisions=(72, 120),
            timeout=height_timeout,
            stack_bytes=STACK_BYTES,
        )
        stable_rank = int(height[-1]["numerical_rank"])
        subset = numerical_subset(pool, height)
        result: dict[str, Any] = {
            "status": "completed exact bounded point checks and numerical height triage",
            "height_bound": height_bound,
            "signed_points_returned": len(raw),
            "distinct_nonzero_ordinate_abscissas": len(signless),
            "abscissas_retained_for_mapping": len(retained),
            "mapping_cap": mapping_cap,
            "mapping_truncated": len(signless) > len(retained),
            "visible_quartic_point_count": len(visible_quartic),
            "visible_jacobian_point_count": len(visible_jacobian),
            "pool_point_count_modulo_inverse": len(pool),
            "pool_point_sha256": point_digest(pool),
            "height_matrix_runs": list(height),
            "stable_numerical_rank": stable_rank,
            "numerical_subset": [point_record(point) for point in subset],
            "numerical_rank_is_not_an_independence_certificate": True,
        }
        if stable_rank >= 14:
            certificate_basis = subset
            try:
                saturated_basis, saturation = saturate_exact_basis(
                    coefficients,
                    subset,
                    prime_bound=50,
                    timeout=saturation_timeout,
                    stack_bytes=STACK_BYTES,
                )
                if len(saturated_basis) != len(subset):
                    raise AssertionError("small-prime saturation changed basis length")
                certificate_basis = saturated_basis
                result["small_prime_saturation"] = saturation
            except CappedProcessTimeout:
                result["small_prime_saturation"] = {
                    "status": "timeout-no-retry",
                    "timeout_seconds": saturation_timeout,
                }
            except Exception as error:
                result["small_prime_saturation"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            certificate = finite_reduction_attempt(
                coefficients,
                certificate_basis,
                prime_bound=SIGNAL_CERTIFICATE_PRIME_BOUND,
            )
            if certificate["status"] == "certified" and certificate.get(
                "two_torsion_certificate_prime"
            ) is None:
                raise AssertionError("a high-rank certificate omitted no 2-torsion")
            result["finite_reduction_attempt"] = certificate
        else:
            result["finite_reduction_attempt"] = {
                "status": "not triggered",
                "trigger_stable_numerical_rank": 14,
            }
        return result, subset
    except CappedProcessTimeout:
        return (
            {
                "status": "timeout-no-retry",
                "height_bound": height_bound,
                "point_timeout_seconds": point_timeout,
                "height_timeout_seconds": height_timeout,
            },
            None,
        )
    except Exception as error:
        return (
            {
                "status": "error-no-retry",
                "height_bound": height_bound,
                "error": str(error)[:1000],
            },
            None,
        )


def phase_rank(phase: dict[str, Any]) -> int:
    return int(phase.get("stable_numerical_rank", -1))


def escalation_key(record: dict[str, Any], prior_rank: int) -> tuple[Any, ...]:
    phase = record["phase"]
    rank = phase_rank(phase)
    conductor = record["conductor"]
    root_number = int(conductor["root_number"])
    expected_parity = 1 if root_number == -1 else 0
    parity_mismatch = rank >= 0 and rank % 2 != expected_parity
    return (
        -rank,
        -int(parity_mismatch),
        -(rank - prior_rank),
        Decimal(conductor["log_conductor"]),
        record["identifier"],
    )


def result_digest(artifact: dict[str, Any]) -> str:
    records = artifact["fiber_records"]
    compact = []
    for record in records:
        certificate = record["H5000_phase"].get("finite_reduction_attempt", {})
        compact.append(
            [
                record["identifier"],
                record["visible_phase"]["status"],
                record["visible_phase"].get("stable_numerical_rank"),
                record["conductor_phase"]["status"],
                record["conductor_phase"].get("conductor"),
                record["H5000_phase"]["status"],
                record["H5000_phase"].get("stable_numerical_rank"),
                certificate.get("certified_algebraic_rank_lower_bound"),
                certificate.get("point_sha256"),
                certificate.get("certificate_primes"),
                certificate.get("two_torsion_certificate_prime"),
            ]
        )
    escalations = {
        stage: [
            [
                record["identifier"],
                record["phase"]["status"],
                record["phase"].get("stable_numerical_rank"),
                record["phase"].get("pool_point_sha256"),
                record["phase"].get("finite_reduction_attempt", {}).get(
                    "certified_algebraic_rank_lower_bound"
                ),
                record["phase"].get("finite_reduction_attempt", {}).get(
                    "point_sha256"
                ),
            ]
            for record in values
        ]
        for stage, values in artifact["escalation_records"].items()
    }
    return stable_json_digest(
        {
            "scope": artifact["scope"],
            "exclusions": artifact["exact_phase_exclusions"],
            "population": artifact["population"],
            "records": compact,
            "escalations": escalations,
            "certified_frontier": artifact["certified_frontier"],
            "hits": artifact["target"]["hits"],
        }
    )


def load_frozen_inputs(root: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cas = root / "elliptic-curves" / "cas"
    generated = root / "artifacts" / "generated-results"
    paths = {
        "compiled_source_sha256": cas / "enumerate_mestre_root_tuples_scale.cpp",
        "max50_driver_sha256": cas / "search_mestre_root_tuple_scale.py",
        "max100_driver_sha256": cas / "search_mestre_root_tuple_scale_max100.py",
        "max50_artifact_sha256": generated / "elliptic_mestre_root_tuple_scale.json",
        "max100_artifact_sha256": generated
        / "elliptic_mestre_root_tuple_scale_max100.json",
    }
    observed = {name: sha256_file(path) for name, path in paths.items()}
    if observed != FROZEN:
        raise AssertionError("a frozen complete-panel input changed")
    max50 = json.loads(paths["max50_artifact_sha256"].read_text())
    max100 = json.loads(paths["max100_artifact_sha256"].read_text())
    return observed, max50, max100


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--enumeration-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=12.0)
    parser.add_argument("--conductor-timeout", type=float, default=10.0)
    parser.add_argument("--point-timeout", type=float, default=15.0)
    parser.add_argument("--escalation-timeout", type=float, default=25.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_mestre_root_tuple_scale_max100_complete.json",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    caps = (
        args.compile_timeout,
        args.enumeration_timeout,
        args.height_timeout,
        args.conductor_timeout,
        args.point_timeout,
        args.escalation_timeout,
    )
    if min(caps) <= 0 or max(caps) > 30:
        raise SystemExit("all subprocess caps must lie in (0,30]")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the complete max-root-100 artifact")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    root = Path(__file__).resolve().parents[2]
    frozen, max50, max100 = load_frozen_inputs(root)
    source = root / "elliptic-curves" / "cas" / "enumerate_mestre_root_tuples_scale.cpp"
    total_started = time.monotonic()

    enumeration, enumeration_timings = compiled_enumeration_max100(
        source,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    if (
        (
            enumeration.normalized_count,
            enumeration.obstruction_count,
            enumeration.reflection_count,
            enumeration.nonreflection_count,
        )
        != EXPECTED_MAX100_COUNTS
    ):
        raise AssertionError("the frozen exact census counts changed")
    verify_enumerator_records(enumeration)
    families, singular, witnesses = classify_nonreflection(enumeration)
    if (
        len(families) != 235
        or len(singular) != 542
        or tuple_digest(families) != EXPECTED_MAX100_NONSINGULAR_SHA256
        or set(witnesses.values()) != {1}
    ):
        raise AssertionError("the frozen exact nonsingular population changed")

    max50_records = {
        record["identifier"]: record
        for record in max50["specialization_screen"][
            "conductor_first_fiber_records"
        ]
    }
    max100_records = {
        record["identifier"]: record
        for record in max100["specialization_screen"]["conductor_records"]
    }
    max50_h5000 = {
        record["identifier"]: record
        for record in max50["specialization_screen"]["point_search_finalists"]
    }
    max100_h5000 = {
        record["identifier"]: record
        for record in max100["specialization_screen"]["h5000_records"]
    }
    if set(max50_records) & set(max100_records) or set(max50_h5000) & set(max100_h5000):
        raise AssertionError("the frozen diameter strata unexpectedly overlap")

    constructions: dict[tuple[int, ...], SixRootMestreConstruction] = {}
    fiber_records: list[dict[str, Any]] = []
    runtime: dict[str, tuple[tuple[int, ...], int]] = {}
    inadmissible: list[dict[str, Any]] = []
    for roots in families:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        constructions[roots] = construction
        for parameter in PARAMETERS:
            parameter_q = Q(parameter)
            discriminant = construction.quartic_discriminant(parameter_q)
            degeneracy = construction.visible_point_degeneracy(parameter_q)
            identifier = fiber_identifier(roots, parameter)
            if (
                discriminant == 0
                or degeneracy.collision_loss
                or degeneracy.zero_ordinates
            ):
                inadmissible.append(
                    {
                        "identifier": identifier,
                        "roots": list(roots),
                        "parameter": parameter,
                        "singular_quartic": discriminant == 0,
                        "visible_collision_loss": degeneracy.collision_loss,
                        "visible_zero_ordinates": degeneracy.zero_ordinates,
                    }
                )
                continue
            runtime[identifier] = (roots, parameter)
            fiber_records.append(
                {
                    "identifier": identifier,
                    "roots": list(roots),
                    "parameter": parameter,
                    "admissibility": {
                        "exact_nonzero_quartic_discriminant": True,
                        "visible_collision_loss": 0,
                        "visible_zero_ordinates": 0,
                    },
                }
            )
    fiber_records.sort(key=lambda record: record["identifier"])
    admissible_ids = set(runtime)
    if not set(max50_records) <= admissible_ids or not set(max100_records) <= admissible_ids:
        raise AssertionError("a frozen conductor exclusion left the exact population")
    if not set(max50_h5000) <= admissible_ids or not set(max100_h5000) <= admissible_ids:
        raise AssertionError("a frozen H5000 exclusion left the exact population")
    if not set(max50_h5000) <= set(max50_records):
        raise AssertionError("a frozen max-root-50 H5000 record lacks a conductor")
    if not set(max100_h5000) <= set(max100_records):
        raise AssertionError("a frozen max-root-100 H5000 record lacks a conductor")

    # Phase 1: exact visible membership, height rank, then finite reductions.
    visible_started = time.monotonic()
    for index, record in enumerate(fiber_records, 1):
        identifier = record["identifier"]
        roots, parameter = runtime[identifier]
        if identifier in max50_records:
            coefficients = constructions[roots].primitive_jacobian_coefficients(
                Q(parameter)
            )
            record["visible_phase"] = imported_visible_record(
                max50_records[identifier], coefficients
            )
        else:
            phase, _ = compute_visible_phase(
                constructions[roots],
                parameter,
                height_timeout=args.height_timeout,
            )
            record["visible_phase"] = phase
        if index % 100 == 0:
            print(f"visible {index}/{len(fiber_records)}", flush=True)
    visible_wall_seconds = time.monotonic() - visible_started

    # Phase 2: close every conductor before launching one new bounded search.
    conductor_started = time.monotonic()
    for index, record in enumerate(fiber_records, 1):
        identifier = record["identifier"]
        roots, parameter = runtime[identifier]
        if identifier in max50_records:
            phase = normalized_conductor(
                max50_records[identifier]["conductor_phase"],
                source="frozen max-root-50",
            )
        elif identifier in max100_records:
            phase = normalized_conductor(
                max100_records[identifier]["conductor_phase"],
                source="frozen max-root-100",
            )
        else:
            coefficients = constructions[roots].primitive_jacobian_coefficients(
                Q(parameter)
            )
            try:
                computed = capped_minimal_curve_data(
                    coefficients,
                    timeout=args.conductor_timeout,
                    stack_bytes=STACK_BYTES,
                )
                phase = {
                    "status": "completed fresh exact PARI minimal-model/conductor computation",
                    **computed,
                    "below_strict_log_conductor_target_numerically": Decimal(
                        computed["log_conductor"]
                    )
                    < TARGET_LOG_CONDUCTOR,
                }
            except CappedProcessTimeout:
                phase = {
                    "status": "timeout-no-retry",
                    "timeout_seconds": args.conductor_timeout,
                }
            except Exception as error:
                phase = {"status": "error-no-retry", "error": str(error)[:1000]}
        record["conductor_phase"] = phase
        if index % 100 == 0:
            print(f"conductor {index}/{len(fiber_records)}", flush=True)
    conductor_wall_seconds = time.monotonic() - conductor_started

    # Phase 3: every H=5000 phase not already frozen is run exactly once.
    h5000_started = time.monotonic()
    fresh_h5000 = 0
    for index, record in enumerate(fiber_records, 1):
        identifier = record["identifier"]
        roots, parameter = runtime[identifier]
        if identifier in max50_h5000:
            phase = normalized_h5000_from_max50(max50_h5000[identifier])
        elif identifier in max100_h5000:
            phase = normalized_h5000_from_max100(max100_h5000[identifier])
        else:
            fresh_h5000 += 1
            phase, _ = bounded_rank_search(
                constructions[roots],
                parameter,
                height_bound=H5000_HEIGHT,
                mapping_cap=H5000_MAPPING_CAP,
                point_timeout=args.point_timeout,
                height_timeout=args.height_timeout,
                saturation_timeout=args.escalation_timeout,
            )
        record["H5000_phase"] = phase
        if index % 50 == 0:
            print(
                f"H5000 {index}/{len(fiber_records)} fresh={fresh_h5000}",
                flush=True,
            )
    h5000_wall_seconds = time.monotonic() - h5000_started

    completed_records = [
        record
        for record in fiber_records
        if phase_rank(record["H5000_phase"]) >= 0
        and "log_conductor" in record["conductor_phase"]
    ]
    by_identifier = {record["identifier"]: record for record in fiber_records}
    escalation_records: dict[str, list[dict[str, Any]]] = {}
    prior_ranks = {
        record["identifier"]: phase_rank(record["H5000_phase"])
        for record in completed_records
    }
    previous = [
        {
            "identifier": record["identifier"],
            "roots": record["roots"],
            "parameter": record["parameter"],
            "conductor": record["conductor_phase"],
            "phase": record["H5000_phase"],
        }
        for record in completed_records
    ]
    escalation_started = time.monotonic()
    for stage, height_bound, keep, mapping_cap in ESCALATION_STAGES:
        ranked = sorted(
            previous,
            key=lambda item: escalation_key(
                item, prior_ranks.get(item["identifier"], -1)
            ),
        )
        selected = ranked[: min(keep, len(ranked))]
        stage_records: list[dict[str, Any]] = []
        for selection_position, item in enumerate(selected, 1):
            identifier = item["identifier"]
            record = by_identifier[identifier]
            roots = tuple(record["roots"])
            parameter = int(record["parameter"])
            phase, _ = bounded_rank_search(
                constructions[roots],
                parameter,
                height_bound=height_bound,
                mapping_cap=mapping_cap,
                point_timeout=args.escalation_timeout,
                height_timeout=args.height_timeout,
                saturation_timeout=args.escalation_timeout,
            )
            stage_records.append(
                {
                    "selection_position": selection_position,
                    "identifier": identifier,
                    "roots": list(roots),
                    "parameter": parameter,
                    "input_stage": "H5000" if stage == "H50000" else previous_stage,
                    "input_stable_numerical_rank": phase_rank(item["phase"]),
                    "selection_rule": (
                        "descending input stable numerical rank, parity mismatch, "
                        "rank gain over H5000, then log conductor and identifier"
                    ),
                    "conductor": record["conductor_phase"],
                    "phase": phase,
                }
            )
            print(
                f"{stage} {selection_position}/{len(selected)} {identifier} "
                f"rank={phase_rank(phase)}",
                flush=True,
            )
        escalation_records[stage] = stage_records
        previous = stage_records
        previous_stage = stage
    escalation_wall_seconds = time.monotonic() - escalation_started

    all_phases: list[tuple[str, str, dict[str, Any]]] = []
    for record in fiber_records:
        all_phases.append((record["identifier"], "H5000", record["H5000_phase"]))
    for stage, records in escalation_records.items():
        for record in records:
            all_phases.append((record["identifier"], stage, record["phase"]))
    numerical_leaders = sorted(
        (
            {
                "identifier": identifier,
                "stage": stage,
                "stable_numerical_rank": phase_rank(phase),
                "finite_reduction_attempt": phase.get("finite_reduction_attempt"),
            }
            for identifier, stage, phase in all_phases
            if phase_rank(phase) >= 0
        ),
        key=lambda record: (
            -record["stable_numerical_rank"],
            record["identifier"],
            record["stage"],
        ),
    )[:25]
    hits = []
    certified_by_identifier: dict[str, dict[str, Any]] = {}
    for identifier, stage, phase in all_phases:
        certificate = phase.get("finite_reduction_attempt", {})
        rank = certificate.get("certified_algebraic_rank_lower_bound")
        if rank is None:
            continue
        conductor = by_identifier[identifier]["conductor_phase"]
        prior = certified_by_identifier.get(identifier)
        if prior is None or rank > prior["certified_algebraic_rank_lower_bound"]:
            roots, parameter = runtime[identifier]
            coefficients = constructions[roots].primitive_jacobian_coefficients(
                Q(parameter)
            )
            certified_by_identifier[identifier] = {
                "identifier": identifier,
                "roots": list(roots),
                "parameter": parameter,
                "primitive_jacobian_coefficients": [
                    str(value) for value in coefficients
                ],
                "stage": stage,
                "certified_algebraic_rank_lower_bound": rank,
                "basis_point_sha256": certificate["point_sha256"],
                "certificate_primes": certificate["certificate_primes"],
                "combined_exact_rank_over_F2": certificate[
                    "combined_exact_rank_over_F2"
                ],
                "two_torsion_certificate_prime": certificate[
                    "two_torsion_certificate_prime"
                ],
                "conductor": conductor,
            }
        if rank >= 30 or (
            rank >= 21
            and "log_conductor" in conductor
            and Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
        ):
            hits.append(
                {
                    "identifier": identifier,
                    "stage": stage,
                    "certified_algebraic_rank_lower_bound": rank,
                    "conductor": conductor["conductor"],
                    "log_conductor": conductor["log_conductor"],
                }
            )

    certified_frontier = sorted(
        certified_by_identifier.values(),
        key=lambda record: (
            -record["certified_algebraic_rank_lower_bound"],
            Decimal(record["conductor"]["log_conductor"]),
            record["identifier"],
        ),
    )

    visible_statuses = Counter(record["visible_phase"]["status"] for record in fiber_records)
    conductor_statuses = Counter(
        record["conductor_phase"]["status"] for record in fiber_records
    )
    h5000_statuses = Counter(record["H5000_phase"]["status"] for record in fiber_records)
    h5000_ranks = Counter(
        str(phase_rank(record["H5000_phase"]))
        for record in fiber_records
        if phase_rank(record["H5000_phase"]) >= 0
    )
    visible_mod2_ranks = Counter(
        str(
            record["visible_phase"]["finite_reduction_of_numerical_subset"].get(
                "combined_exact_rank_over_F2", -1
            )
        )
        for record in fiber_records
    )
    visible_certificate_statuses = Counter(
        record["visible_phase"]["finite_reduction_of_numerical_subset"]["status"]
        for record in fiber_records
    )
    subtarget = sum(
        bool(record["conductor_phase"].get("below_strict_log_conductor_target_numerically"))
        for record in fiber_records
    )
    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "complete exact 235-family T=1..8 max-root-100 Mestre fiber screen; "
            "bounded height ranks are triage unless a finite-reduction certificate is attached"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": hits,
        },
        "scope": {
            "max_root": 100,
            "complete_family_count": len(families),
            "family_sha256": tuple_digest(families),
            "families": [list(roots) for roots in families],
            "integer_parameters": list(PARAMETERS),
            "proposed_fiber_count": len(families) * len(PARAMETERS),
            "proposed_fiber_identifier_sha256": identifier_digest(
                fiber_identifier(roots, parameter)
                for roots in families
                for parameter in PARAMETERS
            ),
        },
        "exact_phase_exclusions": {
            "exclusions_are_phase_specific": True,
            "frozen_max50_visible_and_conductor_count": len(max50_records),
            "frozen_max50_visible_and_conductor_identifier_sha256": identifier_digest(
                max50_records
            ),
            "frozen_max100_conductor_only_count": len(max100_records),
            "frozen_max100_conductor_identifier_sha256": identifier_digest(max100_records),
            "frozen_max50_H5000_count": len(max50_h5000),
            "frozen_max50_H5000_identifier_sha256": identifier_digest(max50_h5000),
            "frozen_max100_H5000_count": len(max100_h5000),
            "frozen_max100_H5000_identifier_sha256": identifier_digest(max100_h5000),
            "visible_fresh_count": len(fiber_records) - len(max50_records),
            "conductor_fresh_count": len(fiber_records)
            - len(max50_records)
            - len(max100_records),
            "H5000_fresh_count": fresh_h5000,
            "no_phase_retry": True,
        },
        "population": {
            "admissible_fiber_count": len(fiber_records),
            "admissible_fiber_identifier_sha256": identifier_digest(admissible_ids),
            "inadmissible_fiber_count": len(inadmissible),
            "inadmissible_fiber_identifier_sha256": identifier_digest(
                record["identifier"] for record in inadmissible
            ),
            "visible_status_histogram": dict(sorted(visible_statuses.items())),
            "visible_finite_reduction_status_histogram": dict(
                sorted(visible_certificate_statuses.items())
            ),
            "visible_combined_exact_mod2_rank_histogram": dict(
                sorted(visible_mod2_ranks.items(), key=lambda item: int(item[0]))
            ),
            "conductor_status_histogram": dict(sorted(conductor_statuses.items())),
            "H5000_status_histogram": dict(sorted(h5000_statuses.items())),
            "H5000_stable_numerical_rank_histogram": dict(
                sorted(h5000_ranks.items(), key=lambda item: int(item[0]))
            ),
            "maximum_H5000_stable_numerical_rank": max(
                (phase_rank(record["H5000_phase"]) for record in fiber_records),
                default=-1,
            ),
            "maximum_any_stage_stable_numerical_rank": max(
                (record["stable_numerical_rank"] for record in numerical_leaders),
                default=-1,
            ),
            "maximum_certified_algebraic_rank_lower_bound": max(
                (
                    record["certified_algebraic_rank_lower_bound"]
                    for record in certified_frontier
                ),
                default=-1,
            ),
            "distinct_certified_signal_count": len(certified_frontier),
            "subtarget_conductor_count": subtarget,
        },
        "inadmissible_fibers": inadmissible,
        "fiber_records": fiber_records,
        "escalation_protocol": {
            "stages": [
                {
                    "name": name,
                    "height_bound": height,
                    "keep": keep,
                    "mapping_cap": cap,
                }
                for name, height, keep, cap in ESCALATION_STAGES
            ],
            "selection_after_complete_H5000_population": True,
            "signal_certificate_trigger_stable_numerical_rank": 14,
            "signal_certificate_prime_bound": SIGNAL_CERTIFICATE_PRIME_BOUND,
        },
        "escalation_records": escalation_records,
        "numerical_leaders": numerical_leaders,
        "certified_frontier": certified_frontier,
        "parameters": {
            "compile_timeout_seconds": args.compile_timeout,
            "enumeration_timeout_seconds": args.enumeration_timeout,
            "height_timeout_seconds": args.height_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "point_timeout_seconds": args.point_timeout,
            "escalation_timeout_seconds": args.escalation_timeout,
            "visible_certificate_prime_bound": VISIBLE_CERTIFICATE_PRIME_BOUND,
            "H5000_height_bound": H5000_HEIGHT,
            "H5000_mapping_cap": H5000_MAPPING_CAP,
            "stack_bytes": STACK_BYTES,
        },
        "frozen_inputs": {**frozen, "all_frozen_files_read_only": True},
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "timings": {
            **enumeration_timings,
            "visible_phase_wall_seconds": visible_wall_seconds,
            "conductor_phase_wall_seconds": conductor_wall_seconds,
            "H5000_phase_wall_seconds": h5000_wall_seconds,
            "escalation_phase_wall_seconds": escalation_wall_seconds,
            "total_wall_seconds": time.monotonic() - total_started,
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "reproducing_command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/search_mestre_root_tuple_scale_max100_complete.py"
            ),
            "all_external_processes_foreground_and_capped": True,
            "whole_process_group_killed_and_reaped_on_timeout": True,
            "no_retries": True,
            "owned_processes_remaining": 0,
        },
    }
    artifact["result_sha256"] = result_digest(artifact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(args.output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "families": len(families),
                "admissible": len(fiber_records),
                "fresh_H5000": fresh_h5000,
                "maximum_H5000_rank": artifact["population"][
                    "maximum_H5000_stable_numerical_rank"
                ],
                "maximum_any_stage_rank": artifact["population"][
                    "maximum_any_stage_stable_numerical_rank"
                ],
                "target_hits": hits,
                "result_sha256": artifact["result_sha256"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
