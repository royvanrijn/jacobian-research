#!/usr/bin/env python3
"""Exhaust the unsearched exact weight-three Fermigier record directions.

The predecessor artifact pins all 6,160 global-sign-quotiented l1=3 vectors
and the 399 directions already searched.  This script regenerates that exact
population without using its modular score, removes those 399 identifiers,
and searches the remaining 5,761 directions in deterministic projective-
height order.  Both slopes are searched once at parameter height 50,000.

Each completed two-slope direction is flushed to a new append-only JSONL
stream.  The slice phase has a declared 1,800-second wall cap checked only
between directions; if it fires, the completed deterministic prefix is pinned
and the remainder stays explicitly open.  No subprocess is retried.

Every genuinely new parameter receives exact conductor first.  Every
subtarget conductor receives forced-point, H=50,000, and H=250,000 height-rank
triage.  H=1,000,000 is attempted only when H=250,000 strictly raises the
stable numerical rank over H=50,000.  Stable rank at least 21 triggers exact
saturation and finite-reduction certification immediately.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Sequence

import sympy as sp

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    extract_parameter_values,
    generic_abscissas,
    published_preimage_digest,
    rational_digest,
    sha256_file,
)
from search_fermigier_published_pair_fiber_products_h50000 import (
    EXPECTED_H50000_RESULT_SHA256,
    h50000_result_digest,
)
from search_fermigier_rank22_auxiliary_orbits import prior_parameter_manifest
from search_fermigier_rank22_accidental_slices import (
    T0,
    canonical_signless_points,
    conductor_probe,
    finite_reduction_attempt,
    generic_group_seed_points,
    point_record,
    quartic_group_pullback,
    search_specialized_quartic,
)
from search_fermigier_rank22_record_group_directions import (
    AUXILIARY_ORBIT_ARTIFACT,
    EXPECTED_PUBLISHED_BASIS_POINT_SHA256,
    EXPECTED_TRANSPORT_SOURCE_SHA256,
    H50000_PAIR_ARTIFACT,
    RecordQuarticAuxiliary,
    aggregate_candidates,
    direction_digest,
    known_record_abscissas,
    load_transport_source,
    projective_height,
    search_direction,
    short_combination,
    slice_result_digest,
    transport_source_digest,
)
from search_fermigier_rank22_record_group_triples import (
    EXPECTED_PREVIOUS_ARTIFACT_SHA256,
    EXPECTED_PREVIOUS_SLICE_RESULT_SHA256,
    generate_triple_directions,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    stable_height_rank,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

TRIPLE_ARTIFACT = "elliptic_fermigier_rank22_record_group_triples.json"
EXPECTED_TRIPLE_ARTIFACT_SHA256 = (
    "2803b1fa276c80eccceac5ce83215f8678d9fb771abdccb8e5043a9962b1ed36"
)
EXPECTED_FULL_VECTOR_DIRECTION_SHA256 = (
    "76bf3cfc5de64af12865209492bf45616470a3e0f5d19e5527551989e1fb24d3"
)
EXPECTED_FULL_NEW_DIRECTION_SHA256 = (
    "4673b556e0a60943d86b23c7f293bfc3a9952f6acd25ef5e152297e84a106455"
)
EXPECTED_SEARCHED_SELECTION_SHA256 = (
    "f22867994c80eecd7f33b48b3bf1d788cf128a9d86d77352187cbc6a9df9815e"
)
EXPECTED_FULL_DIRECTION_COUNT = 6_160
EXPECTED_PREVIOUSLY_SEARCHED_COUNT = 399
EXPECTED_REMAINDER_COUNT = 5_761
EXPECTED_PREDECESSOR_PRIOR_PARAMETER_COUNT = 1_239
EXPECTED_PREDECESSOR_PRIOR_PARAMETER_SHA256 = (
    "9482e61650aa8bb1fd45c3765e5db92c1474090faee8d831e0d73cee4fc864c4"
)
EXPECTED_PRIOR_PARAMETER_COUNT = 1_240
EXPECTED_PRIOR_PARAMETER_SHA256 = (
    "ca5be373dbe7934de0f7ee680203d6ae799f57fb55f712e73ad3777947b20cdd"
)
EXPECTED_AUXILIARY_PARAMETER_COUNT = 666
EXPECTED_AUXILIARY_PARAMETER_SHA256 = (
    "22e296780827d722bed88acb678323349b8f68abf7b298603901dbfdb49a8be1"
)
BOUNDED_EXECUTION_SCRIPT_SHA256 = (
    "bb11ad0b5460d1b10518cf2d8620da255447a3b3f616e7c184a27f73703f3c71"
)
BOUNDED_EXECUTION_AUXILIARY_ARTIFACT_SHA256 = (
    "1008336232ac65bb2bace6ff7008ffb18c1b491c6f0295e16fda14949d5b94d6"
)
BOUNDED_EXECUTION_ARTIFACT_SHA256 = (
    "f8aad4a0ed562c71ed6f6f8a15b8d229dac2d2cb60601d501aeeffd49c3df891"
)
EXPECTED_STREAM_SHA256 = (
    "66e69019bb53b28310bc1a4fa0d40989fe2fd95e677b53bca55feaafa2f3b5de"
)
EXPECTED_STREAM_RESULT_SHA256 = (
    "6818781aed81ca1bfd6822e3f276190d79b42da48b6762f3be9d5638500cb37e"
)
DEFAULT_SLICE_WALL_CAP_SECONDS = 1_800.0
SLICE_HEIGHT = 50_000


def remainder_digest(records: Sequence[dict[str, Any]]) -> str:
    return direction_digest(records)


def load_terminal_prior_parameters(
    artifact_directory: Path, triple_artifact: dict[str, Any]
) -> tuple[set[Fraction], dict[str, Any]]:
    base, base_record = prior_parameter_manifest(artifact_directory)
    auxiliary_path = artifact_directory / AUXILIARY_ORBIT_ARTIFACT
    auxiliary_values = extract_parameter_values(json.loads(auxiliary_path.read_text()))
    auxiliary_digest = rational_digest(sorted(auxiliary_values))
    if (
        len(auxiliary_values) != EXPECTED_AUXILIARY_PARAMETER_COUNT
        or auxiliary_digest != EXPECTED_AUXILIARY_PARAMETER_SHA256
    ):
        raise AssertionError("the exact auxiliary-orbit parameter set changed")

    pair_path = artifact_directory / H50000_PAIR_ARTIFACT
    pair_artifact = json.loads(pair_path.read_text())
    pair_result_digest = h50000_result_digest(pair_artifact["pair_searches"])
    if pair_result_digest != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 pair-product result changed")
    pair_seen = {
        abs(Q(incidence["canonical_parameter_t"]))
        for row in pair_artifact["pair_searches"]
        for incidence in row["search"]["incidences"]
    }
    parameters = base | auxiliary_values | pair_seen
    predecessor_digest = rational_digest(sorted(parameters))
    if (
        len(parameters) != EXPECTED_PREDECESSOR_PRIOR_PARAMETER_COUNT
        or predecessor_digest != EXPECTED_PREDECESSOR_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact predecessor Fermigier prior set changed")
    prior_triple_parameters = {
        abs(Q(candidate["parameter_t"]))
        for candidate in triple_artifact["candidates"]
    }
    parameters |= prior_triple_parameters
    digest = rational_digest(sorted(parameters))
    if (
        len(parameters) != EXPECTED_PRIOR_PARAMETER_COUNT
        or digest != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact terminal Fermigier prior set changed")
    return parameters, {
        "predecessor_parameter_manifest": {
            "base_manifest_parameter_count": len(base),
            "base_manifest_parameter_sha256": base_record["parameter_sha256"],
            "base_manifest": base_record,
            "auxiliary_orbit_artifact": AUXILIARY_ORBIT_ARTIFACT,
            "auxiliary_orbit_artifact_sha256_observed": sha256_file(auxiliary_path),
            "auxiliary_orbit_artifact_sha256_used_as_replay_gate": False,
            "auxiliary_orbit_extracted_parameter_count": len(auxiliary_values),
            "auxiliary_orbit_extracted_parameter_sha256": auxiliary_digest,
            "H50000_pair_artifact": H50000_PAIR_ARTIFACT,
            "H50000_exact_pair_result_sha256": pair_result_digest,
            "H50000_seen_parameters": [
                rational_to_string(value) for value in sorted(pair_seen)
            ],
            "terminal_prior_parameter_count": (
                EXPECTED_PREDECESSOR_PRIOR_PARAMETER_COUNT
            ),
            "terminal_prior_parameter_sha256": predecessor_digest,
        },
        "previous_weight3_candidate_parameters": [
            rational_to_string(value) for value in sorted(prior_triple_parameters)
        ],
        "terminal_prior_parameter_count": len(parameters),
        "terminal_prior_parameter_sha256": digest,
    }


def forced_pool(
    candidate: dict[str, Any],
) -> tuple[tuple[Fraction, Fraction], ...]:
    parameter = Q(candidate["parameter_t"])
    pool = list(generic_group_seed_points(parameter))
    seen_x = {point[0] for point in pool}
    for record in candidate["forced_points"]:
        point = (
            Q(record["basepoint_group_pullback"]["jacobian_x"]),
            Q(record["basepoint_group_pullback"]["jacobian_y"]),
        )
        if point[0] not in seen_x:
            seen_x.add(point[0])
            pool.append(point)
    return tuple(pool)


def exact_height_pool(
    candidate: dict[str, Any],
    *,
    height_bound: int,
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
    precisions: tuple[int, ...],
    saturation_timeout: float,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    parameter = Q(candidate["parameter_t"])
    coefficients = FermigierMestreFamily.coefficients(parameter)
    base_pool = forced_pool(candidate)
    raw_points, search = search_specialized_quartic(
        parameter,
        height_bound=height_bound,
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    record: dict[str, Any] = {
        "height_bound": height_bound,
        "specialized_quartic_search": search,
        "forced_pool_point_count": len(base_pool),
    }
    if search["status"] != "completed":
        return record

    generic_x = generic_abscissas(parameter)
    forced_x = {Q(point["quartic_x"]) for point in candidate["forced_points"]}
    pool = list(base_pool)
    seen_x = {point[0] for point in pool}
    new_points = []
    for quartic_point in canonical_signless_points(raw_points):
        if quartic_point[0] in generic_x or quartic_point[0] in forced_x:
            continue
        pullback = quartic_group_pullback(parameter, quartic_point)
        if pullback is None or pullback[0] in seen_x:
            continue
        seen_x.add(pullback[0])
        pool.append(pullback)
        new_points.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "basepoint_group_pullback": point_record(pullback),
            }
        )
    exact_pool = tuple(pool)
    runs = height_matrix_replay(
        coefficients,
        exact_pool,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    rank = stable_height_rank(runs)
    indices = tuple(runs[-1]["subset_indices_one_based"])
    selected = tuple(exact_pool[index - 1] for index in indices)
    record.update(
        {
            "signed_quartic_point_count": len(raw_points),
            "signless_quartic_abscissa_count": len(
                canonical_signless_points(raw_points)
            ),
            "new_search_group_pullbacks": new_points,
            "exact_pool_point_count": len(exact_pool),
            "stable_numerical_rank": rank,
            "height_runs": list(runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_point_sha256": point_digest(selected),
            "selected_points": [point_record(point) for point in selected],
            "scope_warning": "stable height rank is numerical triage evidence only",
        }
    )
    if rank >= 21:
        record["finite_reduction_attempt"] = finite_reduction_attempt(
            coefficients,
            selected,
            saturation_timeout=saturation_timeout,
            stack_bytes=stack_bytes,
            certificate_prime_bound=certificate_prime_bound,
        )
    return record


def staged_rank_triage(
    candidate: dict[str, Any],
    *,
    h50_timeout: float,
    h250_timeout: float,
    h1m_timeout: float,
    height_timeout: float,
    stack_bytes: int,
    precisions: tuple[int, ...],
    saturation_timeout: float,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    parameter = Q(candidate["parameter_t"])
    coefficients = FermigierMestreFamily.coefficients(parameter)
    base_pool = forced_pool(candidate)
    forced_runs = height_matrix_replay(
        coefficients,
        base_pool,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    record: dict[str, Any] = {
        "forced_pool_point_count": len(base_pool),
        "forced_pool_stable_numerical_rank": stable_height_rank(forced_runs),
        "forced_pool_height_runs": list(forced_runs),
        "stages": [],
        "H1000000_escalation_rule": (
            "attempt only if completed H250000 stable rank strictly exceeds "
            "completed H50000 stable rank"
        ),
    }
    h50 = exact_height_pool(
        candidate,
        height_bound=50_000,
        search_timeout=h50_timeout,
        height_timeout=height_timeout,
        stack_bytes=stack_bytes,
        precisions=precisions,
        saturation_timeout=saturation_timeout,
        certificate_prime_bound=certificate_prime_bound,
    )
    record["stages"].append(h50)
    h250 = exact_height_pool(
        candidate,
        height_bound=250_000,
        search_timeout=h250_timeout,
        height_timeout=height_timeout,
        stack_bytes=stack_bytes,
        precisions=precisions,
        saturation_timeout=saturation_timeout,
        certificate_prime_bound=certificate_prime_bound,
    )
    record["stages"].append(h250)
    h50_rank = h50.get("stable_numerical_rank")
    h250_rank = h250.get("stable_numerical_rank")
    meaningful_gain = (
        h50_rank is not None
        and h250_rank is not None
        and int(h250_rank) > int(h50_rank)
    )
    record["meaningful_rank_gain_at_H250000"] = meaningful_gain
    if meaningful_gain:
        h1m = exact_height_pool(
            candidate,
            height_bound=1_000_000,
            search_timeout=h1m_timeout,
            height_timeout=height_timeout,
            stack_bytes=stack_bytes,
            precisions=precisions,
            saturation_timeout=saturation_timeout,
            certificate_prime_bound=certificate_prime_bound,
        )
        record["stages"].append(h1m)
    ranks = [
        stage["stable_numerical_rank"]
        for stage in record["stages"]
        if "stable_numerical_rank" in stage
    ]
    record["maximum_stable_numerical_rank"] = max(
        [record["forced_pool_stable_numerical_rank"], *ranks]
    )
    return record


def stream_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stream_result_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            for search_row in row["slice_searches"]:
                search = search_row["search"]
                digest.update(
                    (
                        f"{search_row['slice_id']}|"
                        f"{search_row['quartic_polynomial_sha256']}|"
                        f"{search['status']}|{search.get('signed_point_count')}|"
                        f"{search_row['record_T0_calibration_count']}|"
                        f"{len(search_row['incidences'])}\n"
                    ).encode()
                )
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-height", type=int, default=SLICE_HEIGHT)
    parser.add_argument(
        "--slice-wall-cap", type=float, default=DEFAULT_SLICE_WALL_CAP_SECONDS
    )
    parser.add_argument("--slice-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--h50-timeout", type=float, default=30.0)
    parser.add_argument("--h250-timeout", type=float, default=45.0)
    parser.add_argument("--h1m-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated
        / "elliptic_fermigier_rank22_record_group_triples_remainder.json",
    )
    parser.add_argument(
        "--stream-output",
        type=Path,
        default=generated
        / "elliptic_fermigier_rank22_record_group_triples_remainder_stream.jsonl",
    )
    parser.add_argument(
        "--finalize-provenance-metadata",
        action="store_true",
        help=(
            "normalize the completed bounded artifact to stable extracted-parameter "
            "provenance without rerunning any slice search"
        ),
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != SLICE_HEIGHT:
        raise SystemExit("this remainder tranche is pinned at slice H=50000")
    if not 60 <= args.slice_wall_cap <= 3_600:
        raise SystemExit("the slice wall cap must lie in [60,3600] seconds")
    caps = (
        args.slice_timeout,
        args.conductor_timeout,
        args.h50_timeout,
        args.h250_timeout,
        args.h1m_timeout,
        args.height_timeout,
        args.saturation_timeout,
    )
    if min(caps) <= 0 or max(caps) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")
    if args.finalize_provenance_metadata:
        if sha256_file(args.output) != BOUNDED_EXECUTION_ARTIFACT_SHA256:
            raise AssertionError("the raw bounded-execution artifact changed")
        if sha256_file(args.stream_output) != EXPECTED_STREAM_SHA256:
            raise AssertionError("the completed append-only stream changed")
        if stream_result_digest(args.stream_output) != EXPECTED_STREAM_RESULT_SHA256:
            raise AssertionError("the exact completed slice result changed")
        artifact = json.loads(args.output.read_text())
        if (
            artifact.get("script_sha256") != BOUNDED_EXECUTION_SCRIPT_SHA256
            or artifact.get("outcome", {}).get("completed_direction_count")
            != EXPECTED_REMAINDER_COUNT
            or artifact.get("outcome", {}).get("slice_calls_completed")
            != 2 * EXPECTED_REMAINDER_COUNT
            or artifact.get("outcome", {}).get("genuinely_new_forced_fibres") != 0
        ):
            raise AssertionError("the bounded execution summary changed")
        root = Path(__file__).resolve().parents[2]
        artifact_directory = root / "artifacts" / "generated-results"
        triple_path = artifact_directory / TRIPLE_ARTIFACT
        if sha256_file(triple_path) != EXPECTED_TRIPLE_ARTIFACT_SHA256:
            raise AssertionError("the selected triple artifact changed")
        prior_parameters, prior_record = load_terminal_prior_parameters(
            artifact_directory, json.loads(triple_path.read_text())
        )
        if (
            len(prior_parameters) != EXPECTED_PRIOR_PARAMETER_COUNT
            or rational_digest(sorted(prior_parameters))
            != EXPECTED_PRIOR_PARAMETER_SHA256
        ):
            raise AssertionError("stable prior replay failed")
        current_script_sha = sha256_file(Path(__file__).resolve())
        current_auxiliary_sha = prior_record["predecessor_parameter_manifest"][
            "auxiliary_orbit_artifact_sha256_observed"
        ]
        artifact["prior_parameter_decontamination"] = prior_record
        artifact["execution_provenance"] = {
            "bounded_execution_script_sha256": BOUNDED_EXECUTION_SCRIPT_SHA256,
            "bounded_execution_artifact_sha256_before_metadata_normalization": (
                BOUNDED_EXECUTION_ARTIFACT_SHA256
            ),
            "bounded_execution_auxiliary_orbit_artifact_sha256": (
                BOUNDED_EXECUTION_AUXILIARY_ARTIFACT_SHA256
            ),
            "stable_replay_script_sha256": current_script_sha,
            "stable_replay_auxiliary_orbit_artifact_sha256_observed": (
                current_auxiliary_sha
            ),
            "stable_replay_auxiliary_parameter_count": (
                EXPECTED_AUXILIARY_PARAMETER_COUNT
            ),
            "stable_replay_auxiliary_parameter_sha256": (
                EXPECTED_AUXILIARY_PARAMETER_SHA256
            ),
            "metadata_normalization_only": True,
            "bounded_slice_searches_rerun_during_normalization": 0,
            "normalization_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
        }
        artifact["script_sha256"] = current_script_sha
        artifact["metadata_normalized_at_utc"] = datetime.now(
            timezone.utc
        ).isoformat()
        write_artifact(args.output, artifact)
        return
    if args.output.exists() or args.stream_output.exists():
        raise SystemExit("refusing to overwrite a remainder artifact or stream")

    root = Path(__file__).resolve().parents[2]
    artifact_directory = root / "artifacts" / "generated-results"
    primary_path = artifact_directory / PRIMARY_ARTIFACT
    primary = json.loads(primary_path.read_text())
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published preimages changed")
    triple_path = artifact_directory / TRIPLE_ARTIFACT
    if sha256_file(triple_path) != EXPECTED_TRIPLE_ARTIFACT_SHA256:
        raise AssertionError("the selected triple artifact changed")
    triple_artifact = json.loads(triple_path.read_text())
    full = triple_artifact["full_triple_population"]
    selection = triple_artifact["rank_blind_selection"]
    if (
        full["full_vector_direction_sha256"]
        != EXPECTED_FULL_VECTOR_DIRECTION_SHA256
        or full["genuinely_new_direction_sha256"]
        != EXPECTED_FULL_NEW_DIRECTION_SHA256
        or selection["selected_direction_sha256"]
        != EXPECTED_SEARCHED_SELECTION_SHA256
        or selection["selected_direction_count"] != EXPECTED_PREVIOUSLY_SEARCHED_COUNT
    ):
        raise AssertionError("the pinned full/selected triple digests changed")

    auxiliary = RecordQuarticAuxiliary.construct()
    auxiliary_basis, short_basis, _ = load_transport_source(primary, auxiliary)
    if (
        transport_source_digest(primary["published_point_preimages"])
        != EXPECTED_TRANSPORT_SOURCE_SHA256
        or point_digest(short_basis) != EXPECTED_PUBLISHED_BASIS_POINT_SHA256
    ):
        raise AssertionError("the transported rank22 basis changed")
    known_x, _ = known_record_abscissas(primary)
    previous_path = artifact_directory / "elliptic_fermigier_rank22_record_group_directions.json"
    if sha256_file(previous_path) != EXPECTED_PREVIOUS_ARTIFACT_SHA256:
        raise AssertionError("the l1<=2 direction artifact changed")
    previous = json.loads(previous_path.read_text())
    if previous["outcome"]["exact_slice_result_sha256"] != EXPECTED_PREVIOUS_SLICE_RESULT_SHA256:
        raise AssertionError("the l1<=2 exact slice result changed")
    weight2_x = {Q(row["quartic_x"]) for row in previous["direction_searches"]}
    directions, regenerated = generate_triple_directions(
        auxiliary, auxiliary_basis, short_basis, known_x | weight2_x
    )
    if (
        len(directions) != EXPECTED_FULL_DIRECTION_COUNT
        or regenerated["full_vector_direction_sha256"]
        != EXPECTED_FULL_VECTOR_DIRECTION_SHA256
        or regenerated["genuinely_new_direction_sha256"]
        != EXPECTED_FULL_NEW_DIRECTION_SHA256
    ):
        raise AssertionError("the full exact triple population failed replay")
    selected_ids = {
        record["direction_id"] for record in triple_artifact["selected_directions"]
    }
    if len(selected_ids) != EXPECTED_PREVIOUSLY_SEARCHED_COUNT:
        raise AssertionError("the searched triple identifier set changed")
    remainder = [
        record for record in directions if record["direction_id"] not in selected_ids
    ]
    if len(remainder) != EXPECTED_REMAINDER_COUNT:
        raise AssertionError("the exact unsearched triple remainder changed")
    prior_parameters, prior_record = load_terminal_prior_parameters(
        artifact_directory, triple_artifact
    )

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "in-progress deterministic weight3 remainder exhaustion",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "hit": False,
        },
        "source": {
            "primary_artifact": PRIMARY_ARTIFACT,
            "published_accidental_preimage_sha256": published_preimage_digest(primary),
            "transport_source_sha256": transport_source_digest(
                primary["published_point_preimages"]
            ),
            "certified_published_rank22_basis_point_sha256": point_digest(short_basis),
            "selected_triple_artifact": TRIPLE_ARTIFACT,
            "selected_triple_artifact_sha256": EXPECTED_TRIPLE_ARTIFACT_SHA256,
            "full_vector_direction_sha256": regenerated[
                "full_vector_direction_sha256"
            ],
            "full_new_direction_sha256": regenerated[
                "genuinely_new_direction_sha256"
            ],
            "previously_searched_selection_sha256": selection[
                "selected_direction_sha256"
            ],
        },
        "population": {
            "full_exact_weight3_direction_count": len(directions),
            "previously_searched_direction_count": len(selected_ids),
            "declared_remainder_direction_count": len(remainder),
            "declared_remainder_slice_call_count": 2 * len(remainder),
            "remainder_direction_sha256": remainder_digest(remainder),
            "order": "ascending projective height, rational x, direction id",
            "score_used_for_pruning": False,
        },
        "prior_parameter_decontamination": prior_record,
        "parameters": {
            "slice_height": args.slice_height,
            "slice_wall_cap_seconds": args.slice_wall_cap,
            "wall_cap_checked_only_between_two_slope_directions": True,
            "slice_timeout_seconds": args.slice_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "H50000_timeout_seconds": args.h50_timeout,
            "H250000_timeout_seconds": args.h250_timeout,
            "H1000000_timeout_seconds": args.h1m_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": [72, 120],
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
            "append_only_checkpoint_after_each_direction": True,
        },
        "stream": {
            "path": str(args.stream_output),
            "format": "one canonical JSON direction record per line",
            "created_exclusive": True,
        },
        "candidates": [],
        "execution": {
            "phase": "slice-search-in-progress",
            "directions_completed": 0,
            "slice_calls_completed_or_attempted": 0,
            "wall_cap_triggered": False,
            "owned_processes_remaining": None,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    write_artifact(args.output, artifact)

    qualifying_all: list[dict[str, Any]] = []
    classifications = Counter()
    completed = 0
    slice_completed = 0
    slice_failed = 0
    calibrated = 0
    started = time.monotonic()
    args.stream_output.parent.mkdir(parents=True, exist_ok=True)
    with args.stream_output.open("x") as stream:
        for direction in remainder:
            if time.monotonic() - started >= args.slice_wall_cap:
                artifact["execution"]["wall_cap_triggered"] = True
                break
            direction_row, qualifying = search_direction(
                direction,
                height_bound=args.slice_height,
                timeout=args.slice_timeout,
                stack_bytes=args.stack_bytes,
                prior_parameters=prior_parameters,
            )
            stream.write(json.dumps(direction_row, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
            qualifying_all.extend(qualifying)
            completed += 1
            for row in direction_row["slice_searches"]:
                if row["search"]["status"] == "completed":
                    slice_completed += 1
                else:
                    slice_failed += 1
                calibrated += row["record_T0_calibration_count"] == 1
                classifications.update(
                    incidence["classification"] for incidence in row["incidences"]
                )
            artifact["execution"].update(
                {
                    "directions_completed": completed,
                    "slice_calls_completed_or_attempted": completed * 2,
                    "last_direction_id": direction["direction_id"],
                    "wall_seconds_so_far": time.monotonic() - started,
                }
            )
            if completed % 25 == 0:
                write_artifact(args.output, artifact)
                print(
                    f"remainder {completed}/{len(remainder)} "
                    f"qualifying_incidences={len(qualifying_all)}",
                    flush=True,
                )

    slice_phase_wall_seconds = time.monotonic() - started
    candidates = aggregate_candidates(qualifying_all)
    artifact["candidates"] = candidates
    artifact["execution"]["phase"] = "conductor-first"
    write_artifact(args.output, artifact)
    for candidate in candidates:
        parameter = Q(candidate["parameter_t"])
        candidate["conductor_probe"] = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        write_artifact(args.output, artifact)
        if candidate["conductor_probe"].get("below_strict_log_conductor_target"):
            try:
                candidate["staged_rank_triage"] = staged_rank_triage(
                    candidate,
                    h50_timeout=args.h50_timeout,
                    h250_timeout=args.h250_timeout,
                    h1m_timeout=args.h1m_timeout,
                    height_timeout=args.height_timeout,
                    stack_bytes=args.stack_bytes,
                    precisions=(72, 120),
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except subprocess.TimeoutExpired as error:
                candidate["staged_rank_triage"] = {
                    "status": "timeout-no-retry",
                    "error": str(error)[:1000],
                }
            except (RuntimeError, AssertionError, ValueError) as error:
                candidate["staged_rank_triage"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            write_artifact(args.output, artifact)

    remainder_count = len(remainder) - completed
    stream_sha = stream_sha256(args.stream_output)
    exact_result_sha = stream_result_digest(args.stream_output)
    rank_records = [
        candidate["staged_rank_triage"]
        for candidate in candidates
        if "staged_rank_triage" in candidate
        and "maximum_stable_numerical_rank" in candidate["staged_rank_triage"]
    ]
    artifact["stream"].update(
        {
            "completed_direction_line_count": completed,
            "sha256": stream_sha,
            "exact_slice_result_sha256": exact_result_sha,
        }
    )
    artifact["outcome"] = {
        "declared_remainder_direction_count": len(remainder),
        "completed_direction_count": completed,
        "open_remainder_direction_count": remainder_count,
        "open_remainder_claimed_negative": False,
        "full_remainder_exhausted": remainder_count == 0,
        "slice_calls_attempted": completed * 2,
        "slice_calls_completed": slice_completed,
        "slice_calls_timed_out_or_errored": slice_failed,
        "record_T0_calibrated_slices": calibrated,
        "incidence_classification_counts": dict(sorted(classifications.items())),
        "genuinely_new_forced_fibres": len(candidates),
        "completed_conductors": sum(
            candidate.get("conductor_probe", {}).get("status") == "completed"
            for candidate in candidates
        ),
        "subtarget_conductors": sum(
            candidate.get("conductor_probe", {}).get(
                "below_strict_log_conductor_target"
            )
            is True
            for candidate in candidates
        ),
        "rank_triage_count": len(rank_records),
        "maximum_stable_numerical_rank": max(
            (
                record["maximum_stable_numerical_rank"]
                for record in rank_records
            ),
            default=None,
        ),
        "slice_phase_wall_seconds": slice_phase_wall_seconds,
    }
    artifact["target"]["hit"] = any(
        stage.get("finite_reduction_attempt", {}).get(
            "certified_algebraic_rank_lower_bound", 0
        )
        >= 21
        and candidate.get("conductor_probe", {}).get(
            "below_strict_log_conductor_target"
        )
        for candidate in candidates
        for stage in candidate.get("staged_rank_triage", {}).get("stages", [])
    )
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no completed remainder direction produced a certified subtarget rank21 fibre"
        )
    artifact["status"] = (
        "completed full exact weight3 remainder exhaustion"
        if remainder_count == 0
        else "safe-wall-cap checkpoint of deterministic exact weight3 remainder prefix"
    )
    artifact["execution"].update(
        {
            "phase": "complete",
            "wall_seconds": time.monotonic() - started,
            "owned_processes_remaining": 0,
        }
    )
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()
