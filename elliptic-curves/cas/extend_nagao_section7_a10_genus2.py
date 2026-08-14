#!/usr/bin/env python3
"""Final bounded extension of the productive Section-7 ``a10_sp08`` slice.

The completed 240-slice tranche found two new small-conductor parameters on
the single slice ``a10_sp08``.  This standalone extension leaves that artifact
untouched and runs exactly seven one-shot auxiliary searches: the original
normalized sextic at height 1,000,000 and three determinant-one charts around
each of ``T=163`` and ``T=1049/10``, also at height 1,000,000.

All mapped points are verified exactly before canonicalizing ``T``.  The
previous parameter populations and the two center parameters are excluded,
and the remaining points are compared with all twenty-one generic sections
using exact quartic and covariant-Jacobian abscissas.  Only genuinely new
proxy-plausible parameters receive exact conductor calls.  Only exact
sub-182.72 conductors receive the declared specialization searches at heights
50,000 and 250,000 and two-precision numerical height-rank screens.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import signless_quartic_points
from search_nagao_rank20_t5081_neighborhood import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from search_nagao_rank21_accidental_slices import (
    homogenized_transform,
    map_transformed_point,
)
from search_nagao_rank21_t956_skew import search_original_quartic
from search_nagao_section7_accidental_genus2_slices import (
    ROOT,
    candidate_records,
    file_sha256,
    load_genus_two_slices,
    load_prior_parameters,
    parameter_stream_sha256,
    point_record,
    projective_height,
)
from search_nagao_section7_remaining_auxiliary_slices import (
    generic_abscissae,
    generic_labels_for_point,
)
from search_nagao_u42_skew_height import centered_unimodular_matrix
from triage_nagao_rank13_finalists import height_matrix_replay, stable_height_rank


Q = Fraction
SLICE_ID = "a10_sp08"
DIRECT_HEIGHT = 1_000_000
CHART_HEIGHT = 1_000_000
CHART_CENTERS = (Q(163), Q(1049, 10))
CHART_SHIFTS = (-1, 0, 1)
DECLARED_AUXILIARY_CALL_CAP = 7
RANK_TRIAGE_HEIGHTS = (50_000, 250_000)
MAX_PROXY_PLAUSIBLE_PARAMETERS = 16
MAX_RANK_TRIAGE_PARAMETERS = 16
PROXY_THRESHOLD = 190.0
PROXY_TRIAL_BOUND = 2_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")

PREVIOUS_ARTIFACT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_section7_accidental_genus2_slices.json"
)
PREVIOUS_ARTIFACT_SHA256 = (
    "112fb09a12aca5982311a64d161449887074eb13dad52b5b2cb752d93cf7c320"
)
PREVIOUS_SCRIPT_RELATIVE = Path(
    "elliptic-curves/cas/search_nagao_section7_accidental_genus2_slices.py"
)
PREVIOUS_SCRIPT_SHA256 = (
    "abb63f29b1e95eaa23e983e4f8003e682eb2cc4c1b3158ec93aa1f2a694d1e09"
)
SCRIPT_RELATIVE = Path("elliptic-curves/cas/extend_nagao_section7_a10_genus2.py")
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_section7_a10_genus2_extension.json"
)
DEFAULT_OUTPUT = ROOT / OUTPUT_RELATIVE
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/extend_nagao_section7_a10_genus2.py"
)


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def load_extension_inputs() -> tuple[Any, frozenset[Fraction], dict[str, Any]]:
    previous_artifact = ROOT / PREVIOUS_ARTIFACT_RELATIVE
    previous_script = ROOT / PREVIOUS_SCRIPT_RELATIVE
    if file_sha256(previous_artifact) != PREVIOUS_ARTIFACT_SHA256:
        raise RuntimeError("the completed genus-two tranche artifact changed")
    if file_sha256(previous_script) != PREVIOUS_SCRIPT_SHA256:
        raise RuntimeError("the completed genus-two tranche script changed")
    slices, metadata = load_genus_two_slices()
    item = next((value for value in slices if value.identifier == SLICE_ID), None)
    if item is None:
        raise AssertionError("the productive a10_sp08 slice disappeared")
    prior, prior_records = load_prior_parameters()
    payload = json.loads(previous_artifact.read_text(encoding="utf-8"))
    previous_candidates = tuple(
        abs(Q(record["constructor_parameter_T"]))
        for record in payload["new_candidate_population"][
            "records_sorted_by_radical_proxy"
        ]
    )
    if set(previous_candidates) != set(CHART_CENTERS):
        raise AssertionError("the two productive genus-two parameters changed")
    extended_prior = frozenset((*prior, *previous_candidates))
    if len(extended_prior) != 7_631:
        raise AssertionError("the exact extension prior union changed")
    return item, extended_prior, {
        "previous_artifact": str(PREVIOUS_ARTIFACT_RELATIVE),
        "previous_artifact_sha256": PREVIOUS_ARTIFACT_SHA256,
        "previous_script": str(PREVIOUS_SCRIPT_RELATIVE),
        "previous_script_sha256": PREVIOUS_SCRIPT_SHA256,
        "previous_candidate_parameters": [
            rational_to_string(value) for value in previous_candidates
        ],
        "base_prior_artifact_records": list(prior_records),
        "canonical_prior_union_count": len(extended_prior),
        "canonical_prior_union_sha256": parameter_stream_sha256(extended_prior),
        "slice_manifest_sha256": metadata["manifest_sha256"],
    }


def search_plans(item: Any) -> tuple[dict[str, Any], ...]:
    original = tuple(Q(value) for value in item.normalized.normalized_coefficients)
    plans = [
        {
            "id": "direct_H1000000",
            "kind": "direct",
            "height": DIRECT_HEIGHT,
            "matrix": None,
            "center": None,
            "shift": None,
            "polynomial": original,
        }
    ]
    for center in CHART_CENTERS:
        center_label = rational_to_string(center).replace("/", "_")
        for shift in CHART_SHIFTS:
            matrix = centered_unimodular_matrix(center, shift)
            if matrix[0] * matrix[3] - matrix[1] * matrix[2] != 1:
                raise AssertionError("a declared chart lost determinant one")
            if Q(matrix[1], matrix[3]) != center:
                raise AssertionError("a declared chart lost its centre")
            plans.append(
                {
                    "id": f"center_{center_label}_shift_{shift:+d}",
                    "kind": "centered_unimodular",
                    "height": CHART_HEIGHT,
                    "matrix": matrix,
                    "center": center,
                    "shift": shift,
                    "polynomial": homogenized_transform(
                        original, matrix, total_degree=6
                    ),
                }
            )
    if len(plans) != DECLARED_AUXILIARY_CALL_CAP:
        raise AssertionError("the declared seven-plan tranche changed")
    return tuple(plans)


def polynomial_value(
    coefficients: Sequence[Fraction], parameter: Fraction
) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(parameter) + Q(coefficient)
    return answer


def mapped_plan_points(
    item: Any,
    plan: dict[str, Any],
    raw_points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int]:
    mapped = []
    poles = 0
    for raw_point in signless_quartic_points(raw_points):
        local_parameter, local_ordinate = map(Q, raw_point)
        if local_ordinate**2 != polynomial_value(
            plan["polynomial"], local_parameter
        ):
            raise AssertionError("PARI returned a point off a transformed sextic")
        if plan["matrix"] is None:
            point = local_parameter, local_ordinate
        else:
            point = map_transformed_point(
                (local_parameter, local_ordinate),
                plan["matrix"],
                total_degree=6,
            )
            if point is None:
                poles += 1
                continue
        if point[1] ** 2 != item.normalized.normalized_value(point[0]):
            raise AssertionError("a mapped chart point missed the original sextic")
        mapped.append(point)
    return tuple(mapped), poles


def classify_extension_point(
    item: Any,
    normalized_point: tuple[Fraction, Fraction],
    prior: frozenset[Fraction],
    *,
    plan_id: str,
) -> dict[str, Any]:
    signed_parameter, normalized_ordinate = map(Q, normalized_point)
    original_ordinate = item.normalized.original_ordinate(
        signed_parameter, normalized_ordinate
    )
    quartic_x = item.x_value(signed_parameter)
    if original_ordinate**2 != quartic_value(
        section7_primitive_quartic_coefficients(signed_parameter), quartic_x
    ):
        raise AssertionError("an extension point missed the signed Section-7 quartic")
    parameter = abs(signed_parameter)
    if original_ordinate**2 != quartic_value(
        section7_primitive_quartic_coefficients(parameter), quartic_x
    ):
        raise AssertionError("canonicalizing T changed the Section-7 quartic")
    singular = False
    try:
        homogenized_discriminant(parameter)
    except ValueError:
        singular = True
    labels: set[str] = set()
    if not singular:
        for test_parameter in {signed_parameter, parameter}:
            labels.update(
                generic_labels_for_point(
                    test_parameter, (quartic_x, original_ordinate)
                )
            )
    generic_labels = tuple(
        sorted(
            label
            for label in labels
            if label.startswith("quartic-x:")
            or label.startswith("jacobian-sign-pair:")
        )
    )
    ramification = "ramification-ordinate-zero" in labels
    prior_match = parameter in prior
    zero = parameter == 0
    accepted = not (singular or generic_labels or prior_match or zero)
    if accepted:
        classification = "new_forced_non_generic_parameter"
    elif generic_labels:
        classification = "known_generic_section_intersection"
    elif prior_match:
        classification = "prior_parameter_population"
    elif singular:
        classification = "singular_parameter"
    else:
        classification = "zero_parameter"
    return {
        "slice": item.identifier,
        "plan": plan_id,
        "signed_constructor_parameter_T": rational_to_string(signed_parameter),
        "constructor_parameter_T": rational_to_string(parameter),
        "projective_height": projective_height(parameter),
        "normalized_point": point_record((signed_parameter, normalized_ordinate)),
        "forced_quartic_point_on_canonical_fiber": point_record(
            (quartic_x, original_ordinate)
        ),
        "classification": classification,
        "accepted_new_parameter": accepted,
        "prior_population_match": prior_match,
        "singular": singular,
        "ramification_ordinate_zero": ramification,
        "all_21_generic_sections_checked_exactly": not singular,
        "generic_labels": list(generic_labels),
    }


def exact_conductors(
    proxy_records: Sequence[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    selected = tuple(
        record
        for record in proxy_records
        if record["radical_proxy"]["log_radical_upper_proxy"] < PROXY_THRESHOLD
    )
    if len(selected) > MAX_PROXY_PLAUSIBLE_PARAMETERS:
        raise RuntimeError("proxy-plausible extension population exceeds its cap")
    completed = []
    failures = []
    for index, record in enumerate(selected, start=1):
        parameter = Q(record["constructor_parameter_T"])
        try:
            data = minimal_curve_data(
                short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError, FileNotFoundError) as error:
            failures.append(
                {
                    "constructor_parameter_T": rational_to_string(parameter),
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "error": str(error)[:500],
                    "one_attempt_no_retry": True,
                }
            )
            continue
        completed.append(
            {
                **record,
                "status": "completed",
                "minimal_model": [str(value) for value in data["minimal_model"]],
                "conductor": str(data["conductor"]),
                "log_conductor": data["log_conductor"],
                "minimal_discriminant": str(data["minimal_discriminant"]),
                "root_number": data["root_number"],
                "below_strict_182_72_target": (
                    Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
                ),
            }
        )
        print(
            f"conductor {index}/{len(selected)} T={parameter} "
            f"lnN={data['log_conductor']}",
            flush=True,
        )
    return tuple(completed), tuple(failures)


def generic_seed_points(parameter: Fraction) -> tuple[tuple[Fraction, Fraction], ...]:
    quartic = section7_primitive_quartic_coefficients(parameter)
    by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for _, x_value in generic_abscissae(parameter):
        root = rational_square_root(quartic_value(quartic, x_value))
        if root is None:
            raise AssertionError("a declared generic abscissa lost its rational ordinate")
        by_x.setdefault(x_value, (x_value, root))
    return tuple(by_x.values())


def rank_triage_at_height(
    exact_record: dict[str, Any],
    *,
    height: int,
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    parameter = Q(exact_record["constructor_parameter_T"])
    quartic = section7_primitive_quartic_coefficients(parameter)
    raw, search = search_original_quartic(
        quartic,
        str(height),
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    if search["status"] != "completed":
        return {
            "height": height,
            "status": search["status"],
            "quartic_search": search,
            "rank_claim": False,
        }
    searched = signless_quartic_points(raw)
    generic = generic_seed_points(parameter)
    forced = tuple(
        (
            Q(source["forced_quartic_point_on_canonical_fiber"]["x"]),
            Q(source["forced_quartic_point_on_canonical_fiber"]["y"]),
        )
        for source in exact_record["sources"]
    )
    quartic_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in generic + searched + forced:
        if point[1] ** 2 != quartic_value(quartic, point[0]):
            raise AssertionError("a rank-triage point missed the specialized quartic")
        quartic_by_x.setdefault(point[0], point)
    coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
    jacobian_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    exact_rows = []
    for point in quartic_by_x.values():
        image = quartic_point_to_short_jacobian(
            SECTION7_CONSTRUCTION, parameter, point
        )
        jacobian_by_x.setdefault(image[0], image)
        exact_rows.append(
            {"quartic": point_record(point), "jacobian": point_record(image)}
        )
    jacobian_points = tuple(jacobian_by_x.values())
    try:
        runs = height_matrix_replay(
            coefficients,
            jacobian_points,
            precisions=(72, 120),
            timeout=height_timeout,
            stack_bytes=stack_bytes,
        )
        rank = stable_height_rank(runs)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError) as error:
        return {
            "height": height,
            "status": "height_error",
            "error": str(error)[:500],
            "quartic_search": search,
            "rank_claim": False,
        }
    return {
        "height": height,
        "status": "completed",
        "quartic_search": search,
        "signed_quartic_point_count": len(raw),
        "signless_quartic_point_count": len(searched),
        "declared_generic_abscissa_count": len(generic),
        "forced_source_point_count": len(forced),
        "distinct_quartic_abscissa_count": len(quartic_by_x),
        "distinct_jacobian_sign_pair_count": len(jacobian_points),
        "exact_point_rows": exact_rows,
        "height_replay": list(runs),
        "stable_numerical_rank": rank,
        "numerical_only_not_a_rank_certificate": True,
        "rank_at_least_21_triggered": rank >= 21,
        "rank_claim": False,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--auxiliary-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--rank-search-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--proxy-trial-bound", type=int, default=PROXY_TRIAL_BOUND)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "auxiliary_timeout",
        "conductor_timeout",
        "rank_search_timeout",
        "height_timeout",
    ):
        if not 0 < float(getattr(args, name)) <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must lie in (0,60]")
    if not 64_000_000 <= args.stack_bytes <= 1_000_000_000:
        raise SystemExit("--stack-bytes must lie in [64MB,1GB]")
    if not 2 <= args.proxy_trial_bound <= 100_000:
        raise SystemExit("--proxy-trial-bound must lie in [2,100000]")

    started = time.monotonic()
    item, prior, input_metadata = load_extension_inputs()
    plans = search_plans(item)
    run_records = []
    returned_records = []
    sources_by_parameter: dict[Fraction, list[dict[str, Any]]] = {}
    for index, plan in enumerate(plans, start=1):
        raw, process = search_original_quartic(
            plan["polynomial"],
            str(plan["height"]),
            timeout=args.auxiliary_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped, poles = mapped_plan_points(item, plan, raw)
        run_records.append(
            {
                "id": plan["id"],
                "kind": plan["kind"],
                "height": plan["height"],
                "matrix": None if plan["matrix"] is None else list(plan["matrix"]),
                "center": (
                    None
                    if plan["center"] is None
                    else rational_to_string(plan["center"])
                ),
                "shift": plan["shift"],
                "transformed_sextic_coefficients_ascending": [
                    rational_to_string(value) for value in plan["polynomial"]
                ],
                **process,
                "mapped_signless_point_count": len(mapped),
                "chart_pole_count": poles,
                "one_call_no_retry": True,
            }
        )
        for point in mapped:
            record = classify_extension_point(
                item, point, prior, plan_id=plan["id"]
            )
            returned_records.append(record)
            if record["accepted_new_parameter"]:
                sources_by_parameter.setdefault(
                    Q(record["constructor_parameter_T"]), []
                ).append(record)
        print(
            f"auxiliary {index}/{len(plans)} plan={plan['id']} "
            f"mapped={len(mapped)} candidates={len(sources_by_parameter)}",
            flush=True,
        )

    proxies = candidate_records(
        sources_by_parameter, trial_bound=args.proxy_trial_bound
    )
    exact, conductor_failures = exact_conductors(
        proxies,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    subtarget = tuple(
        record for record in exact if record["below_strict_182_72_target"]
    )
    if len(subtarget) > MAX_RANK_TRIAGE_PARAMETERS:
        raise RuntimeError("subtarget extension population exceeds the rank-triage cap")
    rank_records = []
    for record in subtarget:
        parameter = Q(record["constructor_parameter_T"])
        tiers = []
        for height in RANK_TRIAGE_HEIGHTS:
            tier = rank_triage_at_height(
                record,
                height=height,
                search_timeout=args.rank_search_timeout,
                height_timeout=args.height_timeout,
                stack_bytes=args.stack_bytes,
            )
            tiers.append(tier)
            print(
                f"rank T={parameter} H={height} "
                f"stable={tier.get('stable_numerical_rank')}",
                flush=True,
            )
        rank_records.append(
            {
                "constructor_parameter_T": rational_to_string(parameter),
                "tiers": tiers,
                "maximum_stable_numerical_rank": max(
                    (
                        int(tier.get("stable_numerical_rank", -1))
                        for tier in tiers
                    ),
                    default=None,
                ),
                "rank_claim": False,
            }
        )

    status_counts = Counter(record["status"] for record in run_records)
    classification_counts = Counter(
        record["classification"] for record in returned_records
    )
    plan_manifest = [
        {
            "id": plan["id"],
            "kind": plan["kind"],
            "height": plan["height"],
            "matrix": None if plan["matrix"] is None else list(plan["matrix"]),
            "center": (
                None
                if plan["center"] is None
                else rational_to_string(plan["center"])
            ),
            "shift": plan["shift"],
            "polynomial": [rational_to_string(value) for value in plan["polynomial"]],
        }
        for plan in plans
    ]
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_section7_a10_genus2_final_extension",
        "status": "bounded_final_a10_extension_complete",
        "claim_scope": {
            "exact": (
                "seven declared sextic searches and mappings, prior and generic "
                "decontamination, and completed exact conductor outputs"
            ),
            "bounded": (
                "the one direct and six centered H=1000000 boxes only, followed "
                "by H=50000/H=250000 rank triage only on exact subtarget conductors"
            ),
            "rank_certificate": False,
            "numerical_rank_warning": (
                "stable two-precision height ranks are triage evidence only"
            ),
        },
        "inputs": input_metadata,
        "slice": {
            "id": item.identifier,
            "slope": item.slope,
            "intercept": rational_to_string(item.intercept),
            "normalized_sextic_coefficients_ascending": list(
                item.normalized.normalized_coefficients
            ),
        },
        "search_budget": {
            "direct_height": DIRECT_HEIGHT,
            "chart_height": CHART_HEIGHT,
            "chart_centers": [rational_to_string(value) for value in CHART_CENTERS],
            "chart_shifts": list(CHART_SHIFTS),
            "declared_auxiliary_call_cap": DECLARED_AUXILIARY_CALL_CAP,
            "actual_auxiliary_call_count": len(plans),
            "per_call_timeout_seconds": args.auxiliary_timeout,
            "one_call_per_plan_no_retry": True,
            "stack_bytes": args.stack_bytes,
            "plan_manifest_sha256": hashlib.sha256(
                json.dumps(plan_manifest, separators=(",", ":")).encode("utf-8")
            ).hexdigest(),
            "plans": plan_manifest,
        },
        "search_runs": {
            "status_counts": dict(status_counts),
            "records": run_records,
        },
        "returned_point_population": {
            "mapped_signless_incidence_count": len(returned_records),
            "classification_counts": dict(classification_counts),
            "records": returned_records,
            "all_nonsingular_points_checked_against_all_21_generic_sections": True,
        },
        "new_candidate_population": {
            "unique_parameter_count": len(sources_by_parameter),
            "parameter_stream_sha256": parameter_stream_sha256(
                sources_by_parameter
            ),
            "records_sorted_by_radical_proxy": list(proxies),
        },
        "proxy_filter": {
            "strict_log_radical_upper_proxy_threshold": PROXY_THRESHOLD,
            "trial_prime_bound": args.proxy_trial_bound,
            "plausible_parameter_cap": MAX_PROXY_PLAUSIBLE_PARAMETERS,
            "below_threshold_count": sum(
                record["radical_proxy"]["log_radical_upper_proxy"]
                < PROXY_THRESHOLD
                for record in proxies
            ),
        },
        "exact_conductors": {
            "attempted": len(exact) + len(conductor_failures),
            "completed": len(exact),
            "failures": list(conductor_failures),
            "records": list(exact),
            "sub_182_72_count": len(subtarget),
        },
        "rank_triage": {
            "gate": "genuinely new exact conductor with log N strictly below 182.72",
            "heights": list(RANK_TRIAGE_HEIGHTS),
            "parameter_cap": MAX_RANK_TRIAGE_PARAMETERS,
            "parameter_count": len(rank_records),
            "records": rank_records,
            "maximum_stable_numerical_rank": max(
                (
                    int(record["maximum_stable_numerical_rank"])
                    for record in rank_records
                ),
                default=None,
            ),
            "rank_certificate_claimed": False,
        },
        "outcome": {
            "rank21_certified": False,
            "rank30_certified": False,
            "breakthrough_curve_found": False,
            "new_subtarget_conductor_parameters": [
                record["constructor_parameter_T"] for record in subtarget
            ],
            "bounded_lane_stopped_after_this_tranche": True,
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(shlex.quote(value) for value in sys.argv),
            "script_sha256": file_sha256(ROOT / SCRIPT_RELATIVE),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "wall_seconds": time.monotonic() - started,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        f"wrote {args.output} runs={len(plans)} "
        f"new={len(sources_by_parameter)} subtarget={len(subtarget)} "
        f"max_rank={artifact['rank_triage']['maximum_stable_numerical_rank']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
