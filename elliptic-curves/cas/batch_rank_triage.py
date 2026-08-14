#!/usr/bin/env python3
"""Two-stage bounded rank triage for Fermigier specializations.

Every candidate starts with the same modest ``hyperellratpoints`` bound.  The
exact Jacobian images of generic sections 2--13 are always injected as a
twelve-point baseline, independent of whether their quartic abscissas lie in
the bounded search.  Newly enumerated quartic points are mapped exactly and a
PARI height matrix is replayed at two precisions.  At most a configured number
of candidates then advance to the larger bound, ordered by numerical rank gain
and new-abscissa yield.

This is rank triage, not a rank certificate.  Exact curve membership is
verified; height-matrix rank remains explicitly numerical evidence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import time
from typing import Any, Iterable

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_extra_points import (
    exact_point_record,
    height_replay,
    parse_precisions,
    search_quartic_points,
    signless_quartic_points,
)


@dataclass(frozen=True)
class TriageCandidate:
    identifier: str
    parameter: Fraction
    provenance: str
    selection_metadata: str
    known_log_conductor: str | None = None
    known_global_root_number: int | None = None


CANDIDATES = (
    TriageCandidate(
        "primary-1666-9",
        Fraction(1666, 9),
        "record residue-class scan",
        "top score at numerical prime cutoff 500",
        "128.959882907388234115945927583952120187586146367640134572484",
    ),
    TriageCandidate(
        "staged-1547-492",
        Fraction(1547, 492),
        "leakage-free staged record-class rescore",
        "top finalist at numerical prime cutoffs 10000 and 100000",
    ),
    TriageCandidate(
        "staged-833-3464",
        Fraction(833, 3464),
        "leakage-free staged record-class rescore",
        "second finalist at numerical prime cutoff 100000",
    ),
    TriageCandidate(
        "score-714-617",
        Fraction(714, 617),
        "record residue-class score comparison",
        "best retained candidate at numerical prime cutoff 2000",
    ),
    TriageCandidate(
        "score-3332-2571",
        Fraction(3332, 2571),
        "record residue-class score comparison",
        "third candidate at numerical prime cutoff 500",
    ),
    TriageCandidate(
        "score-2975-1493",
        Fraction(2975, 1493),
        "record residue-class score comparison",
        "fourth candidate at numerical prime cutoff 500; sign normalized",
    ),
    TriageCandidate(
        "low-conductor-644-87",
        Fraction(644, 87),
        "expanded multiple-root CRT scan",
        "verified low-conductor candidate; height rank 4 and short-score rank 376",
        "153.964400023010213083969770128430120323448601004558877669472",
        -1,
    ),
    TriageCandidate(
        "low-conductor-154-103",
        Fraction(154, 103),
        "automatic local-condition discovery",
        "lowest-height conductor-only candidate forcing p=7,11,13,17,19",
        "162.234032455648408902235970522085591159150828298194308496774",
        1,
    ),
    TriageCandidate(
        "low-conductor-847-184",
        Fraction(847, 184),
        "height-first multiple-root frontier screen",
        "below-target conductor candidate from the height-first frontier",
        "162.852739805513698244891454752320221986022685402317389213062",
        1,
    ),
    TriageCandidate(
        "low-conductor-70-223",
        Fraction(70, 223),
        "multiple-root CRT scan",
        "lowest-height conductor-engineered control with a low short-prime score",
        "165.979271710943517295247782698287404041467326409598207233727",
    ),
    TriageCandidate(
        "low-conductor-1057-218",
        Fraction(1057, 218),
        "height-first multiple-root frontier screen",
        "below-target conductor candidate from the height-first frontier",
        "169.606323764358951772312077715131200953848832572805084059891",
        1,
    ),
)
DEFAULT_CANDIDATE_IDS = tuple(candidate.identifier for candidate in CANDIDATES)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/batch_rank_triage.py"
)


def candidate_lookup() -> dict[str, TriageCandidate]:
    return {candidate.identifier: candidate for candidate in CANDIDATES}


def parse_candidate_ids(value: str) -> tuple[str, ...]:
    identifiers = tuple(item for item in value.split(",") if item)
    known = candidate_lookup()
    if not identifiers:
        raise argparse.ArgumentTypeError("at least one candidate id is required")
    if len(set(identifiers)) != len(identifiers):
        raise argparse.ArgumentTypeError("candidate ids must be unique")
    unknown = [identifier for identifier in identifiers if identifier not in known]
    if unknown:
        raise argparse.ArgumentTypeError(f"unknown candidate ids: {unknown}")
    return identifiers


def stable_rank(height_runs: Iterable[dict[str, Any]]) -> int:
    ranks = {run["numerical_rank"] for run in height_runs}
    subsets = {tuple(run["subset_indices_one_based"]) for run in height_runs}
    if len(ranks) != 1 or len(subsets) != 1:
        raise AssertionError("numerical rank or selected subset changed with precision")
    return next(iter(ranks))


def unique_new_points(
    parameter: Fraction,
    bounded_points: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]], ...]:
    """Return signless nonvisible quartic points and distinct Jacobian images."""

    visible_x = {
        point[0] for point in FermigierMestreFamily.known_quartic_points(parameter)
    }
    seed_jacobian_x = {
        point[0]
        for point in FermigierMestreFamily.known_jacobian_points(parameter)[1:]
    }
    answers = []
    seen_jacobian_x = set(seed_jacobian_x)
    for quartic_point in signless_quartic_points(bounded_points):
        if quartic_point[0] in visible_x:
            continue
        jacobian_point = FermigierMestreFamily.quartic_point_to_jacobian(
            parameter, quartic_point
        )
        # On the short model, equal x-coordinates are a point/negative pair.
        if jacobian_point[0] in seen_jacobian_x:
            continue
        seen_jacobian_x.add(jacobian_point[0])
        answers.append((quartic_point, jacobian_point))
    return tuple(answers)


def explicit_subset(
    parameter: Fraction,
    source_points: tuple[tuple[Fraction, Fraction], ...],
    jacobian_points: tuple[tuple[Fraction, Fraction], ...],
    indices: list[int],
    seed_count: int,
) -> list[dict[str, Any]]:
    records = []
    for index in indices:
        record: dict[str, Any] = exact_point_record(
            parameter, source_points[index - 1], jacobian_points[index - 1]
        )
        record["pool_index_one_based"] = index
        record["provenance"] = (
            "exact generic-section image" if index <= seed_count else "bounded search"
        )
        records.append(record)
    return records


def run_stage(
    candidate: TriageCandidate,
    *,
    height_bound: int,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
    seed_height_runs: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    parameter = candidate.parameter
    started = time.monotonic()
    bounded_points, search_wall, pari_milliseconds = search_quartic_points(
        parameter,
        height_bound=height_bound,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    signless = signless_quartic_points(bounded_points)
    visible_x = {
        point[0] for point in FermigierMestreFamily.known_quartic_points(parameter)
    }
    new_points = unique_new_points(parameter, bounded_points)
    seed_source = FermigierMestreFamily.known_quartic_points(parameter)[1:]
    seed_jacobian = FermigierMestreFamily.known_jacobian_points(parameter)[1:]
    source_pool = seed_source + tuple(point[0] for point in new_points)
    jacobian_pool = seed_jacobian + tuple(point[1] for point in new_points)
    height_started = time.monotonic()
    height_runs = height_replay(
        parameter,
        jacobian_pool,
        precisions=precisions,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    height_wall = time.monotonic() - height_started
    rank = stable_rank(height_runs)
    seed_rank = stable_rank(seed_height_runs)
    indices = height_runs[-1]["subset_indices_one_based"]
    return {
        "quartic_height_bound": height_bound,
        "search_wall_seconds": search_wall,
        "pari_reported_search_milliseconds": pari_milliseconds,
        "height_replay_wall_seconds": height_wall,
        "total_stage_wall_seconds": time.monotonic() - started,
        "signed_quartic_points_found": len(bounded_points),
        "distinct_quartic_x_values": len(signless),
        "visible_section_x_values_found_within_bound": sum(
            point[0] in visible_x for point in signless
        ),
        "new_x_values_beyond_visible_sections": sum(
            point[0] not in visible_x for point in signless
        ),
        "unique_new_jacobian_images": len(new_points),
        "exact_seed_count": len(seed_jacobian),
        "stable_seed_numerical_rank": seed_rank,
        "stable_pool_numerical_rank": rank,
        "numerical_rank_gain_over_seed": rank - seed_rank,
        "height_matrix_runs": list(height_runs),
        "explicit_numerical_basis": explicit_subset(
            parameter, source_pool, jacobian_pool, indices, len(seed_jacobian)
        ),
        "status": (
            "bounded enumeration and exact membership with numerical height-rank "
            "evidence; no exact independence or rank claim"
        ),
    }


def escalation_key(result: dict[str, Any]) -> tuple[int, int, int, str]:
    stage = result["stages"][0]
    return (
        -stage["stable_pool_numerical_rank"],
        -stage["numerical_rank_gain_over_seed"],
        -stage["new_x_values_beyond_visible_sections"],
        result["candidate_id"],
    )


def select_escalations(
    results: list[dict[str, Any]],
    *,
    minimum_rank_gain: int,
    minimum_new_x: int,
    limit: int,
) -> tuple[str, ...]:
    eligible = [
        result
        for result in results
        if result["stages"][0]["numerical_rank_gain_over_seed"] >= minimum_rank_gain
        or result["stages"][0]["new_x_values_beyond_visible_sections"] >= minimum_new_x
    ]
    eligible.sort(key=escalation_key)
    return tuple(result["candidate_id"] for result in eligible[:limit])


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-ids", type=parse_candidate_ids, default=DEFAULT_CANDIDATE_IDS
    )
    parser.add_argument("--stage-one-height", type=int, default=50_000)
    parser.add_argument("--stage-two-height", type=int, default=1_000_000)
    parser.add_argument("--minimum-rank-gain", type=int, default=3)
    parser.add_argument("--minimum-new-x", type=int, default=20)
    parser.add_argument("--max-escalations", type=int, default=2)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_batch_rank_triage.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.stage_one_height < args.stage_two_height:
        raise SystemExit("require 1 <= stage-one height < stage-two height")
    if args.minimum_rank_gain < 0 or args.minimum_new_x < 0:
        raise SystemExit("escalation thresholds must be nonnegative")
    if args.max_escalations < 0:
        raise SystemExit("--max-escalations must be nonnegative")
    if args.timeout <= 0 or args.stack_bytes < 8_000_000:
        raise SystemExit("PARI timeout and stack bounds must be positive")

    lookup = candidate_lookup()
    selected = [lookup[identifier] for identifier in args.candidate_ids]
    results: list[dict[str, Any]] = []
    experiment_started = time.monotonic()
    for candidate in selected:
        parameter = candidate.parameter
        seed_jacobian = FermigierMestreFamily.known_jacobian_points(parameter)[1:]
        seed_height_runs = height_replay(
            parameter,
            seed_jacobian,
            precisions=args.precisions,
            timeout=args.timeout,
            stack_bytes=args.stack_bytes,
        )
        if stable_rank(seed_height_runs) != 12:
            raise AssertionError(f"{candidate.identifier}: seed was not numerical rank 12")
        first_stage = run_stage(
            candidate,
            height_bound=args.stage_one_height,
            precisions=args.precisions,
            timeout=args.timeout,
            stack_bytes=args.stack_bytes,
            seed_height_runs=seed_height_runs,
        )
        results.append(
            {
                "candidate_id": candidate.identifier,
                "t": rational_to_string(parameter),
                "provenance": candidate.provenance,
                "selection_metadata": candidate.selection_metadata,
                "known_log_conductor": candidate.known_log_conductor,
                "known_global_root_number": candidate.known_global_root_number,
                "seed_height_matrix_runs": list(seed_height_runs),
                "stages": [first_stage],
            }
        )

    escalated_ids = select_escalations(
        results,
        minimum_rank_gain=args.minimum_rank_gain,
        minimum_new_x=args.minimum_new_x,
        limit=args.max_escalations,
    )
    by_id = {result["candidate_id"]: result for result in results}
    for identifier in escalated_ids:
        result = by_id[identifier]
        candidate = lookup[identifier]
        second_stage = run_stage(
            candidate,
            height_bound=args.stage_two_height,
            precisions=args.precisions,
            timeout=args.timeout,
            stack_bytes=args.stack_bytes,
            seed_height_runs=tuple(result["seed_height_matrix_runs"]),
        )
        result["stages"].append(second_stage)

    results.sort(
        key=lambda result: (
            -result["stages"][-1]["stable_pool_numerical_rank"],
            -result["stages"][-1]["new_x_values_beyond_visible_sections"],
            result["candidate_id"],
        )
    )
    maximum_numerical_rank = max(
        result["stages"][-1]["stable_pool_numerical_rank"] for result in results
    )
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded two-stage experiment; rational point membership is exact, "
            "while every height rank is numerical evidence and no rank is claimed"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": "182.72",
            "hit": False,
            "reason": (
                f"maximum stable numerical rank was {maximum_numerical_rank}; "
                "no exact independence certificate was generated"
            ),
        },
        "selection_protocol": {
            "stage_one": (
                "all pinned candidates searched at the same quartic-height bound"
            ),
            "baseline": "exact Jacobian images of generic sections 2--13",
            "escalation_order": (
                "descending stable numerical rank, rank gain, and new-x yield; "
                "then candidate id"
            ),
            "escalation_eligibility": (
                f"rank gain >= {args.minimum_rank_gain} or new x-values >= "
                f"{args.minimum_new_x}"
            ),
            "escalated_candidate_ids": list(escalated_ids),
        },
        "results": results,
        "parameters": {
            "candidate_ids": list(args.candidate_ids),
            "stage_one_quartic_height": args.stage_one_height,
            "stage_two_quartic_height": args.stage_two_height,
            "minimum_rank_gain": args.minimum_rank_gain,
            "minimum_new_x": args.minimum_new_x,
            "max_escalations": args.max_escalations,
            "height_precisions": list(args.precisions),
            "timeout_seconds_per_pari_call": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "experiment_wall_seconds": time.monotonic() - experiment_started,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "sources": {
            "family": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
            "candidate_score_artifacts": [
                "artifacts/generated-results/elliptic_fermigier_score_cutoffs.json",
                "artifacts/generated-results/elliptic_fermigier_record_rescore_h5000.json",
                "artifacts/generated-results/elliptic_fermigier_multiple_root_crt.json",
                "artifacts/generated-results/elliptic_fermigier_discovered_local_conditions.json",
            ],
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"escalated={','.join(escalated_ids) or 'none'}")
    for result in results:
        final_stage = result["stages"][-1]
        print(
            f"{result['candidate_id']} T={result['t']} "
            f"bound={final_stage['quartic_height_bound']} "
            f"new_x={final_stage['new_x_values_beyond_visible_sections']} "
            f"numerical_rank={final_stage['stable_pool_numerical_rank']}"
        )


if __name__ == "__main__":
    main()
