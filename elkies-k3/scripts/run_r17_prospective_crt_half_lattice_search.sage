#!/usr/bin/env sage -python
"""Run the frozen two-stage half-lattice detector on the existing CRT cohort.

Stage A searches the 43 exact generic-deepest classes on every frozen fibre.
Stage B is run exactly when Stage A has certified a new direction and searches
only the new masks in the generic/specialized top-43 union.  Public exceptional
points are never imported by this executable.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import runpy
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from sage.all import EllipticCurve, QQ


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json"
FAMILY_SOURCE = ROOT / "elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage"
ENGINE_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
DEFAULT_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/r17-prospective-crt-half-lattice-pointed-v1"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-pointed-ledger-v1.json"

EXPECTED_CANDIDATE_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
EXPECTED_PROTOCOL_STATUS = "FROZEN_AFTER_POSITIVE_CONTROLS_BEFORE_NEW_COHORT_OUTCOMES"
DIMENSION = 17

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)

engine = SourceFileLoader("r17_crt_half_lattice_search_engine", str(ENGINE_SOURCE)).load_module()
from pointed_quartic_search import run_quartic_search as shared_quartic_search
engine.run_quartic_search = shared_quartic_search
from pointed_quartic_migration import runtime_search, require_runtime, validate_frozen_sources


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def point_key(point: tuple[Fraction, Fraction]):
    return (
        max(
            abs(point[0].numerator).bit_length(),
            point[0].denominator.bit_length(),
            abs(point[1].numerator).bit_length(),
            point[1].denominator.bit_length(),
        ),
        point,
    )


def load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = json.loads(MANIFEST.read_text())
    protocol = json.loads(PROTOCOL.read_text())
    if manifest.get("status") != "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS":
        raise ArithmeticError("the original CRT manifest is no longer frozen and unopened")
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the frozen candidate list commitment changed")
    if protocol.get("status") != EXPECTED_PROTOCOL_STATUS:
        raise ArithmeticError("the replacement detector protocol is not frozen")
    if protocol.get("candidate_list_sha256") != EXPECTED_CANDIDATE_HASH:
        raise ArithmeticError("the replacement detector names another cohort")
    if protocol["protocol_definition_sha256"] != canonical_hash(
        {
            key: value
            for key, value in protocol.items()
            if key not in {"protocol_definition_sha256", "inputs", "generation"}
        }
    ):
        raise ArithmeticError("the replacement detector protocol hash does not replay")
    declared_runner = protocol.get("inputs", {}).get(relative(Path(__file__)))
    if declared_runner is not None:
        validate_frozen_sources({relative(Path(__file__)): declared_runner})
    return manifest, protocol


def as_fraction(value) -> Fraction:
    return Fraction(str(value))


def curve_python_data(curve, known):
    local_minimal = curve.local_data(2).minimal_model()
    isomorphisms = curve.isomorphisms(local_minimal)
    if not isomorphisms:
        raise ArithmeticError("no exact isomorphism to the deterministic p=2-minimal model")
    to_local = isomorphisms[0]
    from_local = ~to_local
    local_points = tuple(to_local(point) for point in known)
    a1, a2, a3, a4, a6 = map(QQ, local_minimal.a_invariants())
    if any(value.denominator() != 1 for value in (a1, a2, a3, a4, a6)):
        raise ArithmeticError("the deterministic p=2-minimal model is not integral")
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    model = tuple(
        as_fraction(value) for value in (0, 0, 0, -27 * c4, -54 * c6)
    )
    points = tuple(
        (
            as_fraction(36 * point[0] + 3 * b2),
            as_fraction(108 * (2 * point[1] + a1 * point[0] + a3)),
        )
        for point in local_points
    )
    if len(points) != DIMENSION:
        raise ArithmeticError("the specialized generic basis stopped having dimension 17")
    short_curve = EllipticCurve(QQ, list(model))
    if any(short_curve(point) == short_curve(0) for point in points):
        raise ArithmeticError("a normalized generic point became the identity")

    def short_to_original(point: tuple[Fraction, Fraction]):
        short_x, short_y = map(QQ, point)
        local_x = (short_x - 3 * b2) / 36
        local_y = (short_y / 108 - a1 * local_x - a3) / 2
        local_point = local_minimal(local_x, local_y)
        original = from_local(local_point)
        if original not in curve or to_local(original) != local_point:
            raise ArithmeticError("short-to-original exact point transport failed")
        return as_fraction(original[0]), as_fraction(original[1])

    normalization = {
        "source_model": [str(value) for value in curve.a_invariants()],
        "deterministic_intermediate": "Sage local_data(2).minimal_model(), first exact isomorphism",
        "p2_minimal_model": [str(value) for value in local_minimal.a_invariants()],
        "certificate_short_model": [str(value) for value in model],
        "certificate_short_model_sha256": canonical_hash([str(value) for value in model]),
        "all_seventeen_points_transported_exactly": True,
    }
    return model, points, short_to_original, normalization


def rounded_height_gram(height_gram, scale: int):
    return tuple(
        tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
        for row in height_gram
    )


def representatives_for_masks(height_gram, scale: int, masks: Iterable[int]):
    oracle = engine.CosetOracle(rounded_height_gram(height_gram, scale))
    representatives = {}
    maximum_error = 0.0
    for mask in masks:
        unused_norm, representative, error = oracle.solve(int(mask))
        representatives[int(mask)] = representative
        maximum_error = max(maximum_error, error)
    return representatives, maximum_error


def full_specialized_ranking(height_gram, scale: int):
    oracle = engine.CosetOracle(rounded_height_gram(height_gram, scale))
    rows = []
    maximum_error = 0.0
    for mask in range(1 << DIMENSION):
        unused_norm, representative, error = oracle.solve(mask)
        depth = engine.quadratic_decimal(height_gram, representative) / 4
        rows.append((depth, mask, representative))
        maximum_error = max(maximum_error, error)
    rows.sort(key=lambda row: (-row[0], row[1]))
    return rows, maximum_error


def compact_cover_record(stage: str, mask: int, representative, outcome, cpu_seconds: float):
    record = outcome.record
    return {
        "stage": stage,
        "mask": mask,
        "hex": f"0x{mask:05x}",
        "representative": list(map(int, representative)),
        "status": record["status"],
        "cpu_seconds": cpu_seconds,
        "wall_seconds": record["wall_seconds"],
        "integral_model_maximum_coefficient_bits": record.get("maximum_coefficient_bits"),
        "pointed_search": record,
        "reduced_model_maximum_coefficient_bits": (
            record.get("reduced_model", {}).get("maximum_coefficient_bits")
        ),
        "minimalization_milliseconds": record.get("minimalization_milliseconds"),
        "reduction_milliseconds": record.get("reduction_milliseconds"),
        "search_milliseconds": record.get("search_milliseconds"),
        "finite_curve_point_count": len(outcome.curve_points),
        "error": record.get("error"),
    }


def run_cover_set(
    *,
    stage: str,
    masks: Sequence[int],
    representatives: dict[int, Sequence[int]],
    model,
    generic_points,
    protocol,
):
    pipeline = protocol["cover_pipeline"]
    discoveries: dict[tuple[Fraction, Fraction], set[int]] = {}
    records = []
    for position, mask in enumerate(masks, start=1):
        started_cpu = cpu_clock()
        outcome = engine.run_quartic_search(
            mask=mask,
            representative=representatives[mask],
            short_model=model,
            generic_points=generic_points,
            height_bound=pipeline["height_bound_each_cover"],
            timeout_seconds=pipeline[
                "wall_timeout_seconds_each_cover_including_minimize_reduce_search"
            ],
            stack_bytes=pipeline["gp_stack_bytes_each_cover"],
        )
        cpu_seconds = cpu_clock() - started_cpu
        records.append(
            compact_cover_record(stage, mask, representatives[mask], outcome, cpu_seconds)
        )
        for point in outcome.curve_points:
            discoveries.setdefault(point, set()).add(mask)
        print(
            f"R17CRTHALF|stage={stage}|cover={position}/{len(masks)}|"
            f"mask={mask:#07x}|status={outcome.record['status']}|"
            f"points={len(outcome.curve_points)}",
            flush=True,
        )
    return records, discoveries


def certificate_record(signatures, column_count: int):
    return {
        "prime_bound": 1_000,
        "combined_rank": combined_mod2_rank(signatures, column_count),
        "column_count": column_count,
        "signatures": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "doubled_subgroup_order": row.doubled_subgroup_order,
                "quotient_dimension": row.quotient_dimension,
                "rows": [list(vector) for vector in row.rows],
            }
            for row in signatures
        ],
    }


def certify_discoveries(
    *,
    stage: str,
    model,
    generic_points,
    discoveries,
    already_selected,
    already_seen,
    prime_bound: int,
    short_to_original,
):
    selected = list(already_selected)
    accepted = []
    uncertified = []
    basis_signs = {
        signed
        for point in generic_points
        for signed in (point, (point[0], -point[1]))
    }
    candidates = sorted(
        (
            point
            for point in discoveries
            if point not in basis_signs and point not in already_seen
        ),
        key=point_key,
    )
    for point in candidates:
        trial = tuple(generic_points) + tuple(selected) + (point,)
        signatures = find_mod2_reduction_certificate(
            model, trial, prime_bound=prime_bound
        )
        rank = combined_mod2_rank(signatures, len(trial))
        original_point = short_to_original(point)
        if rank == len(trial):
            selected.append(point)
            accepted.append(
                {
                    "stage": stage,
                    "short_certificate_point": point_record(point),
                    "original_specialization_point": point_record(original_point),
                    "source_masks": sorted(discoveries[point]),
                    "exact_short_curve_equation_verified": True,
                    "exact_original_specialization_equation_verified": True,
                    "nonmembership_in_preceding_certified_subgroup": True,
                    "Q_linear_independence_certified_by_primitive_relation_mod2": True,
                    "finite_reduction_certificate": certificate_record(
                        signatures, len(trial)
                    ),
                }
            )
        else:
            uncertified.append(
                {
                    "stage": stage,
                    "short_certificate_point": point_record(point),
                    "original_specialization_point": point_record(original_point),
                    "source_masks": sorted(discoveries[point]),
                    "exact_short_curve_equation_verified": True,
                    "exact_original_specialization_equation_verified": True,
                    "reason_not_counted": "FINITE_REDUCTION_RANK_DID_NOT_REACH_FULL_COLUMN_COUNT",
                    "achieved_rank": rank,
                    "column_count": len(trial),
                }
            )
    return tuple(selected), accepted, uncertified, set(candidates)


def run_fibre(row: dict[str, Any], protocol: dict[str, Any], family) -> dict[str, Any]:
    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    parameter = row["parameter"]
    curve, known = family.specialize(parameter)
    model, generic_points, short_to_original, normalization = curve_python_data(
        curve, known
    )
    generic_signatures = find_mod2_reduction_certificate(
        model,
        generic_points,
        prime_bound=protocol["point_acceptance"]["finite_reduction_prime_bound"],
    )
    if combined_mod2_rank(generic_signatures, DIMENSION) != DIMENSION:
        raise ArithmeticError("the specialized generic MW17 basis lost exact independence")

    height_started = time.monotonic()
    height_gram = engine.canonical_height_gram(model, generic_points)
    height_seconds = time.monotonic() - height_started
    generic_masks = tuple(
        int(mask)
        for mask in protocol["native_generic_lattice"][
            "deepest_masks_in_norm_then_mask_order"
        ]
    )
    scale = protocol["specialized_representative_policy"][
        "operative_integer_rounding_scale"
    ]
    stage_a_representatives, stage_a_cvp_error = representatives_for_masks(
        height_gram, scale, generic_masks
    )
    stage_a_covers, stage_a_discoveries = run_cover_set(
        stage="A",
        masks=generic_masks,
        representatives=stage_a_representatives,
        model=model,
        generic_points=generic_points,
        protocol=protocol,
    )
    selected, stage_a_points, stage_a_uncertified, stage_a_seen = certify_discoveries(
        stage="A",
        model=model,
        generic_points=generic_points,
        discoveries=stage_a_discoveries,
        already_selected=(),
        already_seen=set(),
        prime_bound=protocol["point_acceptance"]["finite_reduction_prime_bound"],
        short_to_original=short_to_original,
    )
    stage_a_gain = len(stage_a_points)

    stage_b = {
        "gate_satisfied": stage_a_gain > 0,
        "gate_reason": (
            "AT_LEAST_ONE_EXACTLY_CERTIFIED_STAGE_A_DIRECTION"
            if stage_a_gain
            else "NO_EXACTLY_CERTIFIED_STAGE_A_DIRECTION"
        ),
        "full_specialized_ranking_computed": False,
        "covers": [],
        "certified_points": [],
        "uncertified_candidates": [],
    }
    if stage_a_gain:
        audit_scale = protocol["specialized_representative_policy"][
            "audit_integer_rounding_scale"
        ]
        audit_rows, audit_error = full_specialized_ranking(height_gram, audit_scale)
        operative_rows, operative_error = full_specialized_ranking(height_gram, scale)
        top_count = protocol["stage_b"]["specialized_top_count"]
        audit_top = tuple(mask for unused_depth, mask, unused_rep in audit_rows[:top_count])
        operative_top_rows = operative_rows[:top_count]
        operative_top = tuple(mask for unused_depth, mask, unused_rep in operative_top_rows)
        operative_representatives = {
            mask: representative
            for unused_depth, mask, representative in operative_top_rows
        }
        generic_set = set(generic_masks)
        new_masks = tuple(mask for mask in operative_top if mask not in generic_set)
        stage_b_covers, stage_b_discoveries = run_cover_set(
            stage="B",
            masks=new_masks,
            representatives=operative_representatives,
            model=model,
            generic_points=generic_points,
            protocol=protocol,
        )
        selected, stage_b_points, stage_b_uncertified, stage_b_seen = certify_discoveries(
            stage="B",
            model=model,
            generic_points=generic_points,
            discoveries=stage_b_discoveries,
            already_selected=selected,
            already_seen=stage_a_seen,
            prime_bound=protocol["point_acceptance"]["finite_reduction_prime_bound"],
            short_to_original=short_to_original,
        )
        stage_b = {
            "gate_satisfied": True,
            "gate_reason": "AT_LEAST_ONE_EXACTLY_CERTIFIED_STAGE_A_DIRECTION",
            "full_specialized_ranking_computed": True,
            "specialized_top43_masks": list(operative_top),
            "specialized_top43_hex": [f"0x{mask:05x}" for mask in operative_top],
            "top43_set_stable_at_scales_1e5_and_1e6": set(audit_top)
            == set(operative_top),
            "top43_order_stable_at_scales_1e5_and_1e6": audit_top == operative_top,
            "generic_specialized_intersection_count": len(
                generic_set.intersection(operative_top)
            ),
            "union_class_count": len(generic_set.union(operative_top)),
            "incremental_cover_count": len(new_masks),
            "audit_scale_maximum_cvp_error": audit_error,
            "operative_scale_maximum_cvp_error": operative_error,
            "covers": stage_b_covers,
            "certified_points": stage_b_points,
            "uncertified_candidates": stage_b_uncertified,
        }

    all_covers = stage_a_covers + stage_b["covers"]
    backend_failures = sum(row["status"] == "pari_failure" for row in all_covers)
    bounded_timeouts = sum(row["status"] == "bounded_search_timeout" for row in all_covers)
    total_gain = len(selected)
    return {
        "sample_id": row["sample_id"],
        "manifest_index": row["manifest_index"],
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": row["anchor_curve_id"],
        "cohort": row["cohort"],
        "parameter": str(parameter),
        "status": (
            "CENSORED_COVER_BACKEND_FAILURE"
            if backend_failures
            else (
                "CERTIFIED_STAGE_A_ESCAPE"
                if stage_a_gain
                else "BOUNDED_STAGE_A_NO_CERTIFIED_ESCAPE"
            )
        ),
        "failure": None,
        "normalization": normalization,
        "generic_subgroup": {
            "rank": DIMENSION,
            "exact_short_curve_equation_and_section_identities_verified": True,
            "finite_reduction_certificate": certificate_record(
                generic_signatures, DIMENSION
            ),
        },
        "height_computation_seconds": height_seconds,
        "stage_a": {
            "cover_count": len(stage_a_covers),
            "maximum_cvp_error": stage_a_cvp_error,
            "covers": stage_a_covers,
            "certified_quotient_gain": stage_a_gain,
            "certified_points": stage_a_points,
            "uncertified_candidates": stage_a_uncertified,
        },
        "stage_b": {
            **stage_b,
            "incremental_certified_quotient_gain": len(
                stage_b["certified_points"]
            ),
        },
        "total_certified_quotient_gain": total_gain,
        "largest_certified_rank_lower_bound": DIMENSION + total_gain,
        "bounded_cover_timeout_count": bounded_timeouts,
        "cover_backend_failure_count": backend_failures,
        "analysis_eligible_complete_stage_a": backend_failures == 0,
        "timing": {
            "worker_wall_seconds": time.monotonic() - started_wall,
            "worker_cpu_seconds_parent_plus_children": cpu_clock() - started_cpu,
            "cover_cpu_seconds": sum(row["cpu_seconds"] for row in all_covers),
        },
        "claim_boundary": (
            "Every counted point has an exact equation check and full finite-reduction "
            "independence certificate. All misses remain bounded-search results."
        ),
    }


def failure_record(row, status: str, failure: dict[str, Any], elapsed: float):
    return {
        "sample_id": row["sample_id"],
        "manifest_index": row["manifest_index"],
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": row["anchor_curve_id"],
        "cohort": row["cohort"],
        "parameter": row["parameter"],
        "status": status,
        "failure": failure,
        "analysis_eligible_complete_stage_a": False,
        "supervisor_wall_seconds": elapsed,
    }


def parse_worker(completed, row, elapsed: float):
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    result_line = next(
        (line for line in reversed(lines) if line.startswith("RESULT_JSON=")), None
    )
    if completed.returncode == 0 and result_line is not None:
        result = json.loads(result_line[len("RESULT_JSON=") :])
        result["supervisor_wall_seconds"] = elapsed
        return result
    return failure_record(
        row,
        "CENSORED_FIBRE_WORKER_FAILURE",
        {"returncode": completed.returncode, "output_tail": lines[-30:]},
        elapsed,
    )


def timeout_record(error, row, elapsed: float):
    output = error.stdout or ""
    if isinstance(output, bytes):
        output = output.decode(errors="replace")
    lines = [line for line in output.splitlines() if line.strip()]
    return failure_record(
        row,
        "CENSORED_FIBRE_WORKER_TIMEOUT",
        {"output_tail": lines[-30:]},
        elapsed,
    )


def manifest_row(raw: dict[str, Any], index: int):
    return {
        "sample_id": raw["sample_id"],
        "manifest_index": index,
        "match_set_id": raw["match_set_id"],
        "anchor_curve_id": raw["anchor_curve_id"],
        "cohort": raw["cohort"],
        "parameter": raw["parameter"],
    }


def run_single(index: int) -> None:
    manifest, protocol = load_inputs()
    if not 0 <= index < len(manifest["rows"]):
        raise ValueError("single index is outside the frozen manifest")
    address_space = protocol["fibre_worker_envelope"]["address_space_bytes"]
    if address_space is not None:
        resource.setrlimit(resource.RLIMIT_AS, (address_space, address_space))
    family = runpy.run_path(str(FAMILY_SOURCE))["Family"]()
    result = run_fibre(manifest_row(manifest["rows"][index], index), protocol, family)
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def write_checkpoint(path: Path, protocol_hash: str, chunk_index: int, chunk_count: int, indices, records):
    document = {
        "schema": "elkies-k3.r17-prospective-crt-half-lattice-chunk.v3",
        "status": (
            "COMPLETE_HALF_LATTICE_CHUNK"
            if len(records) == len(indices)
            else "PARTIAL_HALF_LATTICE_CHECKPOINT"
        ),
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "runtime_search": runtime_search(),
        "protocol_definition_sha256": protocol_hash,
        "chunk_index": chunk_index,
        "chunk_count": chunk_count,
        "scheduled_indices": indices,
        "completed_record_count": len(records),
        "records": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run_chunk(chunk_index: int, chunk_count: int, output: Path, limit: int | None):
    manifest, protocol = load_inputs()
    indices = [
        index for index in range(len(manifest["rows"])) if index % chunk_count == chunk_index
    ]
    if limit is not None:
        indices = indices[:limit]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        require_runtime(old)
        if (
            old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("scheduled_indices") != indices
            or old.get("protocol_definition_sha256")
            != protocol["protocol_definition_sha256"]
        ):
            raise ArithmeticError("an existing checkpoint belongs to another protocol or chunk")
        records = old["records"]
    completed_ids = {row["sample_id"] for row in records}
    timeout = protocol["fibre_worker_envelope"]["wall_timeout_seconds"]
    for position, index in enumerate(indices, start=1):
        raw = manifest["rows"][index]
        if raw["sample_id"] in completed_ids:
            continue
        row = manifest_row(raw, index)
        command = [sys.executable, str(Path(__file__).resolve()), "--single-index", str(index)]
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            result = parse_worker(completed, row, time.monotonic() - started)
        except subprocess.TimeoutExpired as error:
            result = timeout_record(error, row, time.monotonic() - started)
        records.append(result)
        write_checkpoint(
            output,
            protocol["protocol_definition_sha256"],
            chunk_index,
            chunk_count,
            indices,
            records,
        )
        print(
            f"R17CRTHALFCHUNK|chunk={chunk_index}/{chunk_count}|"
            f"completed={position}/{len(indices)}|sample={row['sample_id']}|"
            f"status={result['status']}|stageA={result.get('stage_a', {}).get('certified_quotient_gain')}",
            flush=True,
        )


def merge_chunks(chunk_dir: Path, chunk_count: int, output: Path):
    manifest, protocol = load_inputs()
    records_by_id = {}
    chunk_provenance = []
    for index in range(chunk_count):
        path = chunk_dir / f"chunk-{index:02d}-of-{chunk_count:02d}.json"
        chunk = json.loads(path.read_text())
        require_runtime(chunk)
        if chunk.get("status") != "COMPLETE_HALF_LATTICE_CHUNK":
            raise ArithmeticError(f"half-lattice chunk {index} is incomplete")
        if chunk.get("protocol_definition_sha256") != protocol["protocol_definition_sha256"]:
            raise ArithmeticError(f"half-lattice chunk {index} used another protocol")
        for row in chunk["records"]:
            if row["sample_id"] in records_by_id:
                raise ArithmeticError("duplicate sample across half-lattice chunks")
            records_by_id[row["sample_id"]] = row
        chunk_provenance.append(
            {
                "path": relative(path),
                "sha256": digest(path),
                "record_count": len(chunk["records"]),
            }
        )
    expected_ids = [row["sample_id"] for row in manifest["rows"]]
    if set(records_by_id) != set(expected_ids):
        raise ArithmeticError("half-lattice chunks do not cover the frozen manifest")
    records = [records_by_id[sample_id] for sample_id in expected_ids]
    status_counts = Counter(row["status"] for row in records)
    cohort_counts = defaultdict(Counter)
    for row in records:
        cohort_counts[row["cohort"]][row["status"]] += 1
    document = {
        "schema": "elkies-k3.r17-prospective-crt-half-lattice-ledger.v3",
        "status": "COMPLETE_FROZEN_HALF_LATTICE_DETECTOR_LEDGER",
        "candidate_list_sha256": EXPECTED_CANDIDATE_HASH,
        "runtime_search": runtime_search(),
        "protocol_definition_sha256": protocol["protocol_definition_sha256"],
        "summary": {
            "scheduled_fibres": len(records),
            "status_counts": dict(sorted(status_counts.items())),
            "status_counts_by_cohort": {
                cohort: dict(sorted(counts.items()))
                for cohort, counts in sorted(cohort_counts.items())
            },
            "stage_a_certified_escape_rows": sum(
                row.get("stage_a", {}).get("certified_quotient_gain", 0) > 0
                for row in records
            ),
            "stage_a_certified_directions": sum(
                row.get("stage_a", {}).get("certified_quotient_gain", 0)
                for row in records
            ),
            "stage_b_executed_rows": sum(
                row.get("stage_b", {}).get("full_specialized_ranking_computed", False)
                for row in records
            ),
            "stage_b_incremental_certified_directions": sum(
                row.get("stage_b", {}).get("incremental_certified_quotient_gain", 0)
                for row in records
            ),
            "censored_rows": sum(
                not row.get("analysis_eligible_complete_stage_a", False) for row in records
            ),
        },
        "records": records,
        "chunk_provenance": chunk_provenance,
        "inputs": {
            relative(MANIFEST): digest(MANIFEST),
            relative(PROTOCOL): digest(PROTOCOL),
            relative(FAMILY_SOURCE): digest(FAMILY_SOURCE),
            relative(ENGINE_SOURCE): digest(ENGINE_SOURCE),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "commands": [
                f"sage -python {relative(Path(__file__))} --chunk-index I --chunk-count {chunk_count}",
                f"sage -python {relative(Path(__file__))} --merge --chunk-count {chunk_count}",
            ],
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"R17CRTHALFMERGE|records={len(records)}|"
        f"events={document['summary']['stage_a_certified_escape_rows']}|"
        f"output={relative(output)}",
        flush=True,
    )


def run_positive_control(curve_id: int) -> None:
    unused_manifest, protocol = load_inputs()
    address_space = protocol["fibre_worker_envelope"]["address_space_bytes"]
    if address_space is not None:
        resource.setrlimit(resource.RLIMIT_AS, (address_space, address_space))
    family = runpy.run_path(str(FAMILY_SOURCE))["Family"]()
    if curve_id not in family.target_parameters:
        raise ValueError("positive-control curve id is not a 074d9 target")
    row = {
        "sample_id": f"positive-control-{curve_id}",
        "manifest_index": -1,
        "match_set_id": "external-positive-control",
        "anchor_curve_id": curve_id,
        "cohort": "EXTERNAL_POSITIVE_CONTROL",
        "parameter": str(family.target_parameters[curve_id]),
    }
    result = run_fibre(row, protocol, family)
    expected = protocol["positive_control_acceptance"][f"curve{curve_id}-rank29"]
    if result["stage_a"]["certified_quotient_gain"] != expected[
        "stage_a_exact_quotient_rank"
    ]:
        raise ArithmeticError("native-basis Stage-A control rank differs from frozen holdout")
    if result["total_certified_quotient_gain"] != expected[
        "stage_b_union_exact_quotient_rank"
    ]:
        raise ArithmeticError("native-basis union control rank differs from frozen holdout")
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--positive-control", type=int, choices=(356, 385))
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=32)
    parser.add_argument("--chunk-dir", type=Path, default=DEFAULT_CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    modes = sum(
        value is not None
        for value in (args.single_index, args.positive_control, args.chunk_index)
    ) + int(args.merge)
    if modes != 1:
        raise SystemExit("choose exactly one of --single-index, --positive-control, --chunk-index, or --merge")
    if args.single_index is not None:
        run_single(args.single_index)
    elif args.positive_control is not None:
        run_positive_control(args.positive_control)
    elif args.merge:
        merge_chunks(
            args.chunk_dir.resolve(), args.chunk_count, (args.output or OUTPUT).resolve()
        )
    else:
        if not 0 <= args.chunk_index < args.chunk_count:
            raise SystemExit("chunk index is outside chunk count")
        output = args.output or (
            args.chunk_dir
            / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
        )
        run_chunk(
            args.chunk_index,
            args.chunk_count,
            output.resolve(),
            args.limit,
        )


if __name__ == "__main__":
    main()
