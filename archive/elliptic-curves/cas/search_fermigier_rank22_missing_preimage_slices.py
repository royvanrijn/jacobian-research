#!/usr/bin/env python3
"""Search the five published Fermigier preimage directions missed at H=10^6.

The exact quartic preimages of published points P14, P15, P20, P21, and P22
are absent from the 14-source record-fiber catalog used by the first
constructive slice tranche.  Each point defines both exact genus-one slices

``x = +T+n`` and ``x = -T+n``.

This script searches those ten quartics once at H=200000.  It removes exact
generic-section collisions, intersections already represented by the prior
14-source slice catalog, the record fiber, and a pinned 78-fiber Fermigier
prior-search set.  Every surviving parameter receives a bounded exact
conductor computation.  Only completed specializations with log(N)<182.72
receive an H=50000 specialized-quartic rank screen; only stable numerical rank
at least 16 can trigger a single H=10^6 escalation.  There are no retries.

Numerical height rank remains triage.  An exact target claim is emitted only
after the imported finite-reduction certificate routine succeeds.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any, Iterable, Sequence

import sympy as sp

from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_rank22_accidental_slices import (
    TARGET_LOG_CONDUCTOR,
    T0,
    build_slices,
    canonical_signless_points,
    conductor_probe,
    finite_reduction_attempt,
    point_digest,
    point_record,
    projective_height,
    quartic_group_pullback,
    search_slice,
    select_reconstruction_convention,
    specialized_quartic_screen,
)


Q = Fraction
MISSING_LABELS = ("P14", "P15", "P20", "P21", "P22")
LEGACY_SLICE_PARAMETERS = tuple(
    Q(value) for value in ("19033/135", "22253/114", "31331/104", "38633/138")
)
PRIOR_FILENAMES = (
    "elliptic_fermigier_batch_rank_triage.json",
    "elliptic_fermigier_rank22_accidental_slices.json",
    "elliptic_fermigier_record_residue_deep_tranche.json",
)
EXPECTED_CURRENT15_SHA256 = (
    "e8410fbcba4491165fd114e86cff11a1eabc1373f74cebda8dbc856bbcf0045f"
)
EXPECTED_BASELINE30_SHA256 = (
    "fb76dd2e8726aff4d6c2742750189ee8a794e0e8939b7ede40bd62d69fd4cbbd"
)
EXPECTED_DEEP48_SHA256 = (
    "d452648e2a645f8d80e9bda4070e6950cbbda525a2a4159ed9937ff4e2940df9"
)
EXPECTED_EXTENDED78_SHA256 = (
    "285f8dca1ab945c3d116f8f286ad50b77eaff9bff2f979048ac6ae683c14c867"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_rank22_missing_preimage_slices.json"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter_digest(parameters: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for parameter in sorted(set(Q(value) for value in parameters)):
        digest.update((rational_to_string(parameter) + "\n").encode())
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def missing_preimages() -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    _, _, reconstruction, _ = select_reconstruction_convention()
    by_label = {
        f"P{index}": point
        for index, point in enumerate(reconstruction, start=1)
    }
    answer = tuple((label, by_label[label]) for label in MISSING_LABELS)
    expected_x = {
        "P14": Q(1185892, 741),
        "P15": Q(1161926, 975),
        "P20": Q(-1717790, 6747),
        "P21": Q(1756342, 8229),
        "P22": Q(1038854, 12831),
    }
    if any(point[0] != expected_x[label] for label, point in answer):
        raise AssertionError("a missing published preimage changed")
    return answer


def load_prior_state(
    generated_dir: Path,
) -> tuple[set[Fraction], tuple[tuple[str, tuple[Fraction, Fraction]], ...], dict[str, Any]]:
    paths = {name: generated_dir / name for name in PRIOR_FILENAMES}
    batch = load_json(paths[PRIOR_FILENAMES[0]])
    current = load_json(paths[PRIOR_FILENAMES[1]])
    deep = load_json(paths[PRIOR_FILENAMES[2]])

    batch_parameters = {abs(Q(record["t"])) for record in batch["results"]}
    current_parameters = {
        abs(Q(record["parameter_t"]))
        for record in current["candidate_conductor_screen"]
    }
    current_parameters.add(abs(Q(current["record_parameter_normalized_T"])))
    deep_parameters = {abs(Q(record["t"])) for record in deep["candidates"]}
    legacy_parameters = {abs(parameter) for parameter in LEGACY_SLICE_PARAMETERS}
    if (len(batch_parameters), len(current_parameters), len(deep_parameters)) != (
        11,
        15,
        48,
    ):
        raise AssertionError("the pinned Fermigier prior-search counts changed")
    baseline = batch_parameters | current_parameters | legacy_parameters
    exclusions = baseline | deep_parameters
    if len(baseline) != 30 or len(exclusions) != 78:
        raise AssertionError("the pinned Fermigier prior exclusions changed")
    digests = (
        parameter_digest(current_parameters),
        parameter_digest(baseline),
        parameter_digest(deep_parameters),
        parameter_digest(exclusions),
    )
    if digests != (
        EXPECTED_CURRENT15_SHA256,
        EXPECTED_BASELINE30_SHA256,
        EXPECTED_DEEP48_SHA256,
        EXPECTED_EXTENDED78_SHA256,
    ):
        raise AssertionError("the pinned Fermigier prior parameter digests changed")

    priority_sources = tuple(
        (
            str(record["label"]),
            (Q(record["quartic_x"]), Q(record["quartic_z"])),
        )
        for record in current["record_fiber_height_1000000_replay"][
            "priority_sources"
        ]
    )
    if len(priority_sources) != 14:
        raise AssertionError("the prior source catalog no longer has 14 points")
    manifest = {
        "canonicalization": "T -> abs(T), because the family is exactly even",
        "batch_rank_triage_parameters": len(batch_parameters),
        "current_slice_parameters_including_record": len(current_parameters),
        "legacy_superseded_slice_parameters": [
            rational_to_string(value) for value in LEGACY_SLICE_PARAMETERS
        ],
        "baseline_Fermigier_only_parameters": len(baseline),
        "deep_tranche_parameters": len(deep_parameters),
        "total_unique_prior_Fermigier_parameters": len(exclusions),
        "current15_parameter_sha256": digests[0],
        "baseline30_parameter_sha256": digests[1],
        "deep48_parameter_sha256": digests[2],
        "extended78_parameter_sha256": digests[3],
        "source_sha256": {
            name: sha256_file(path) for name, path in paths.items()
        },
    }
    return exclusions, priority_sources, manifest


def prior_source_slice_ids(
    raw_parameter: Fraction,
    x_value: Fraction,
    priority_sources: Sequence[tuple[str, tuple[Fraction, Fraction]]],
) -> tuple[str, ...]:
    matches = []
    for label, point in priority_sources:
        for slope in (-1, 1):
            intercept = point[0] - slope * T0
            if slope * raw_parameter + intercept == x_value:
                suffix = "m1" if slope == -1 else "p1"
                matches.append(f"{label.lower()}_{suffix}")
    return tuple(sorted(matches))


def classify_hit(
    raw_parameter: Fraction,
    x_value: Fraction,
    *,
    prior_parameters: set[Fraction],
    priority_sources: Sequence[tuple[str, tuple[Fraction, Fraction]]],
) -> tuple[str, tuple[str, ...]]:
    normalized = abs(Q(raw_parameter))
    if normalized == 0:
        return "zero parameter", ()
    if normalized == abs(T0):
        return "record source fiber", ()
    if FermigierMestreFamily.discriminant_factor(normalized) == 0:
        return "singular fiber", ()
    generic_x = {
        point[0]
        for signed_parameter in (normalized, -normalized)
        for point in FermigierMestreFamily.known_quartic_points(signed_parameter)
    }
    if x_value in generic_x:
        return "generic-or-T-sign-conjugate collision", ()
    prior_slice_ids = prior_source_slice_ids(
        Q(raw_parameter), Q(x_value), priority_sources
    )
    if prior_slice_ids:
        return "prior 14-source slice incidence", prior_slice_ids
    if normalized in prior_parameters:
        return "previously searched Fermigier fiber", ()
    return "genuinely new fiber incidence", ()


def aggregate_new_hits(
    slices: Sequence[Any],
    search_results: Sequence[tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]],
    *,
    prior_parameters: set[Fraction],
    priority_sources: Sequence[tuple[str, tuple[Fraction, Fraction]]],
) -> tuple[dict[Fraction, list[dict[str, Any]]], list[dict[str, Any]]]:
    by_parameter: dict[Fraction, dict[Fraction, dict[str, Any]]] = defaultdict(dict)
    all_incidences = []
    for slice_data, (raw_points, _) in zip(slices, search_results, strict=True):
        for raw_parameter, ordinate in canonical_signless_points(raw_points):
            normalized = abs(raw_parameter)
            x_value = slice_data.x_value(raw_parameter)
            if ordinate**2 != FermigierMestreFamily.quartic_value(
                normalized, x_value
            ):
                raise AssertionError("a missing-direction slice point missed its fiber")
            classification, prior_slice_ids = classify_hit(
                raw_parameter,
                x_value,
                prior_parameters=prior_parameters,
                priority_sources=priority_sources,
            )
            incidence = {
                "slice_id": slice_data.identifier,
                "raw_parameter_t": rational_to_string(raw_parameter),
                "normalized_parameter_t": rational_to_string(normalized),
                "quartic_x": rational_to_string(x_value),
                "quartic_z": rational_to_string(ordinate),
                "classification": classification,
                "prior_source_slice_ids": list(prior_slice_ids),
            }
            all_incidences.append(incidence)
            if classification != "genuinely new fiber incidence":
                continue
            existing = by_parameter[normalized].get(x_value)
            if existing is None:
                by_parameter[normalized][x_value] = {
                    "ordinate": abs(ordinate),
                    "slice_ids": {slice_data.identifier},
                    "raw_parameters": {raw_parameter},
                }
            else:
                existing["slice_ids"].add(slice_data.identifier)
                existing["raw_parameters"].add(raw_parameter)

    records: dict[Fraction, list[dict[str, Any]]] = {}
    for parameter, by_x in by_parameter.items():
        point_records = []
        for x_value, source in sorted(by_x.items()):
            point = (x_value, Q(source["ordinate"]))
            pullback = quartic_group_pullback(parameter, point)
            if pullback is None:
                raise AssertionError("a genuinely new point pulled back to zero")
            point_records.append(
                {
                    "quartic_x": rational_to_string(x_value),
                    "quartic_z": rational_to_string(point[1]),
                    "slice_ids": sorted(source["slice_ids"]),
                    "raw_parameters": [
                        rational_to_string(value)
                        for value in sorted(source["raw_parameters"])
                    ],
                    "basepoint_group_pullback": point_record(pullback),
                    "exact_group_pullback_checked": True,
                }
            )
        records[parameter] = point_records
    return records, all_incidences


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        precisions = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(precisions) < 2 or precisions != tuple(sorted(set(precisions))):
        raise argparse.ArgumentTypeError("provide increasing distinct precisions")
    return precisions


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-height", type=int, default=200_000)
    parser.add_argument("--slice-timeout", type=float, default=10.0)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--specialization-height", type=int, default=50_000)
    parser.add_argument("--specialization-timeout", type=float, default=15.0)
    parser.add_argument("--escalation-height", type=int, default=1_000_000)
    parser.add_argument("--escalation-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--height-precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--generated-dir",
        type=Path,
        default=root / "artifacts" / "generated-results",
    )
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if (
        args.slice_height != 200_000
        or args.specialization_height != 50_000
        or args.escalation_height != 1_000_000
    ):
        raise SystemExit("the canonical H bounds are 200000, 50000, and 1000000")
    timeouts = (
        args.slice_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.escalation_timeout,
        args.height_timeout,
        args.saturation_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("every subprocess timeout must lie in (0,60]")
    if args.stack_bytes < 8_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid stack or certificate prime bound")

    sources = missing_preimages()
    slices = build_slices(sources)
    if len(slices) != 10:
        raise AssertionError("the missing-direction construction needs ten slices")
    prior_parameters, priority_sources, prior_manifest = load_prior_state(
        args.generated_dir
    )

    search_results = []
    slice_records = []
    for slice_data in slices:
        points, search = search_slice(
            slice_data,
            height_bound=args.slice_height,
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
        )
        search_results.append((points, search))
        slice_records.append(
            {
                "slice_id": slice_data.identifier,
                "source_label": slice_data.accidental_label,
                "source_quartic_point": {
                    "x": rational_to_string(slice_data.source_point[0]),
                    "z": rational_to_string(slice_data.source_point[1]),
                },
                "slope": slice_data.slope,
                "intercept": rational_to_string(slice_data.intercept),
                "auxiliary_coefficients_ascending": [
                    rational_to_string(value) for value in slice_data.coefficients
                ],
                "auxiliary_degree": len(slice_data.coefficients) - 1,
                "auxiliary_genus": 1,
                "source_parameter_replay": True,
                "search": search,
            }
        )
        print(
            f"slice={slice_data.identifier} status={search['status']} "
            f"parameters={search.get('distinct_parameter_count', 0)}",
            flush=True,
        )

    new_hits, all_incidences = aggregate_new_hits(
        slices,
        search_results,
        prior_parameters=prior_parameters,
        priority_sources=priority_sources,
    )
    candidate_records = []
    target_hits = []
    alternative_hits = []
    for parameter in sorted(
        new_hits,
        key=lambda value: (projective_height(value), value),
    ):
        conductor = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        candidate = {
            "parameter_t": rational_to_string(parameter),
            "projective_height": projective_height(parameter),
            "forced_points": new_hits[parameter],
            "conductor_probe": conductor,
            "specialized_quartic_height_50000_screen": None,
            "height_1000000_escalation": None,
            "finite_reduction_certificate": None,
        }
        if (
            conductor.get("status") == "completed"
            and conductor["below_strict_log_conductor_target"]
        ):
            try:
                screen, selected = specialized_quartic_screen(
                    parameter,
                    height_bound=args.specialization_height,
                    search_timeout=args.specialization_timeout,
                    height_timeout=args.height_timeout,
                    precisions=args.height_precisions,
                    stack_bytes=args.stack_bytes,
                )
                candidate["specialized_quartic_height_50000_screen"] = screen
                rank = int(screen["height_rank"]["stable_numerical_rank"])
                if rank >= 16:
                    escalation, escalation_selected = specialized_quartic_screen(
                        parameter,
                        height_bound=args.escalation_height,
                        search_timeout=args.escalation_timeout,
                        height_timeout=args.height_timeout,
                        precisions=args.height_precisions,
                        stack_bytes=args.stack_bytes,
                    )
                    candidate["height_1000000_escalation"] = escalation
                    escalation_rank = int(
                        escalation["height_rank"]["stable_numerical_rank"]
                    )
                    selected = escalation_selected
                    certificate = None
                    if escalation_rank >= 21:
                        certificate = finite_reduction_attempt(
                            FermigierMestreFamily.coefficients(parameter),
                            selected,
                            saturation_timeout=args.saturation_timeout,
                            stack_bytes=args.stack_bytes,
                            certificate_prime_bound=args.certificate_prime_bound,
                        )
                        candidate["finite_reduction_certificate"] = certificate
                        if (
                            certificate["status"] == "certified"
                            and int(
                                certificate[
                                    "certified_algebraic_rank_lower_bound"
                                ]
                            )
                            >= 21
                        ):
                            candidate["exact_log_conductor_bound"] = (
                                exact_log_conductor_certificate(
                                    int(conductor["conductor"])
                                )
                            )
                            target_hits.append(rational_to_string(parameter))
                    if (
                        escalation_rank >= 30
                        and certificate is not None
                        and certificate["status"] == "certified"
                        and int(
                            certificate["certified_algebraic_rank_lower_bound"]
                        )
                        >= 30
                    ):
                        alternative_hits.append(rational_to_string(parameter))
            except (
                subprocess.TimeoutExpired,
                RuntimeError,
                AssertionError,
                ValueError,
            ) as error:
                candidate["specialization_error"] = str(error)[:1000]
        candidate_records.append(candidate)
        print(
            f"T={parameter} conductor={conductor['status']} "
            f"logN={conductor.get('log_conductor')} "
            f"rank={((candidate['specialized_quartic_height_50000_screen'] or {}).get('height_rank') or {}).get('stable_numerical_rank')}",
            flush=True,
        )

    classification_counts: dict[str, int] = defaultdict(int)
    for incidence in all_incidences:
        classification_counts[incidence["classification"]] += 1
    completed_conductors = [
        record
        for record in candidate_records
        if record["conductor_probe"]["status"] == "completed"
    ]
    below_target = [
        record
        for record in completed_conductors
        if record["conductor_probe"]["below_strict_log_conductor_target"]
    ]
    screens = [
        record["specialized_quartic_height_50000_screen"]
        for record in candidate_records
        if record["specialized_quartic_height_50000_screen"] is not None
    ]
    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "exact five-preimage construction plus bounded ten-slice, conductor, "
            "and specialization triage; no target claim without finite reductions"
        ),
        "record_parameter_normalized_T": rational_to_string(T0),
        "missing_published_preimages": [
            {
                "label": label,
                "quartic_x": rational_to_string(point[0]),
                "quartic_z": rational_to_string(point[1]),
                "exact_quartic_membership_checked": (
                    point[1] ** 2
                    == FermigierMestreFamily.quartic_value(T0, point[0])
                ),
                "basepoint_group_pullback": point_record(
                    quartic_group_pullback(T0, point)
                ),
            }
            for label, point in sources
        ],
        "prior_decontamination": prior_manifest,
        "slices": slice_records,
        "all_signless_parameter_incidences": all_incidences,
        "candidate_conductor_screen": candidate_records,
        "summary": {
            "missing_published_preimages": len(sources),
            "genus_one_slices": len(slices),
            "slice_searches_completed": sum(
                record["search"]["status"] == "completed"
                for record in slice_records
            ),
            "slice_search_timeouts": sum(
                record["search"]["status"] == "timeout"
                for record in slice_records
            ),
            "incidence_classification_counts": dict(
                sorted(classification_counts.items())
            ),
            "genuinely_new_parameters": len(candidate_records),
            "new_parameters": [
                record["parameter_t"] for record in candidate_records
            ],
            "conductor_calls_completed": len(completed_conductors),
            "completed_conductors_below_strict_target": len(below_target),
            "specialized_height_50000_screens": len(screens),
            "specialized_height_50000_ranks": {
                screen["parameter_t"]: screen["height_rank"][
                    "stable_numerical_rank"
                ]
                for screen in screens
            },
            "height_1000000_escalations": sum(
                record["height_1000000_escalation"] is not None
                for record in candidate_records
            ),
            "target_hits": target_hits,
            "alternative_rank_hits": alternative_hits,
        },
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "low_conductor_certified_hits": target_hits,
            "alternative_rank_certified_hits": alternative_hits,
        },
        "parameters": {
            "slice_height_bound": args.slice_height,
            "slice_timeout_seconds_per_process": args.slice_timeout,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "specialization_height_bound": args.specialization_height,
            "specialization_timeout_seconds_per_process": (
                args.specialization_timeout
            ),
            "escalation_minimum_stable_numerical_rank": 16,
            "escalation_height_bound": args.escalation_height,
            "escalation_timeout_seconds_per_process": args.escalation_timeout,
            "height_precisions": list(args.height_precisions),
            "height_timeout_seconds_per_candidate": args.height_timeout,
            "saturation_timeout_seconds_per_candidate": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "pari_stack_bytes": args.stack_bytes,
            "retries": 0,
            "temporary_directory_dependency": None,
            "generated_dir": str(args.generated_dir),
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"sources={len(sources)} slices={len(slices)} "
        f"new={len(candidate_records)} below={len(below_target)} "
        f"screens={len(screens)} escalations={artifact['summary']['height_1000000_escalations']}"
    )


if __name__ == "__main__":
    main()
