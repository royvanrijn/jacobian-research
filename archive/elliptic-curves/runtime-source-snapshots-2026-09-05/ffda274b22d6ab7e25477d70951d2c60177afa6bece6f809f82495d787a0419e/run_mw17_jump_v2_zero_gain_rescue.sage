#!/usr/bin/env sage-python
"""Run the prospective MW17-jump-v2 zero-gain rescue arm.

Only hash-assigned fibres with a clean exact zero from the original 43-chart
phase are eligible.  The rescue searches the next 301 generic half-classes in
seven batches of 43.  Once a batch certifies the first quotient direction, the
remaining slots switch to the existing adaptive quotient policy.  Thus every
treated fibre retains the original 344-chart total cap.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time
from typing import Any

from sage.all import Matrix, ZZ


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-mw17-jump-v2-zero-gain-rescue-arm-v1.json"
)
CAMPAIGN = (
    ROOT / "artifacts/generated-results/elkies-k3-mw17-jump-v2-campaign-v1.json"
)
BASE_RUNNER = ROOT / "elliptic-curves/cas/run_mw17_jump_v2.sage"
PRODUCTION_GATES = ROOT / "elliptic-curves/cas/production_search_gates.py"
BASE_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/mw17-jump-v2-pointed-v1"
CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/mw17-jump-v2-zero-gain-rescue-pointed-v1"
LEDGER = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-mw17-jump-v2-zero-gain-rescue-pointed-ledger-v1.json"
)
STOP_SENTINEL = CHUNK_DIR / "STOP_GAIN15.json"

EXPECTED_STATUS = "FROZEN_PROSPECTIVE_ZERO_GAIN_RESCUE_BEFORE_RESCUE_SEARCH"
GENERIC_DIMENSION = 17
BASE_INITIAL_CHARTS = 43
RESCUE_CHARTS = 301
RESCUE_BATCH_SIZE = 43
TOTAL_CHART_CAP = 344
HEIGHT_BOUND = 100_000
CHART_TIMEOUT_SECONDS = 15.0
STACK_BYTES = 1_000_000_000
RELATION_CHUNK_SIZE = 64
RELATION_TIMEOUT_SECONDS = 180.0


from pointed_quartic_migration import validate_frozen_sources, runtime_search, require_runtime

def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def atomic_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def load_base_module():
    return SourceFileLoader("mw17_jump_v2_rescue_base", str(BASE_RUNNER)).load_module()


def load_gate_module():
    return SourceFileLoader(
        "mw17_jump_v2_rescue_production_gates", str(PRODUCTION_GATES)
    ).load_module()


def load_protocol() -> tuple[dict[str, Any], dict[str, Any]]:
    protocol = json.loads(PROTOCOL.read_text())
    campaign = json.loads(CAMPAIGN.read_text())
    if protocol.get("status") != EXPECTED_STATUS:
        raise ArithmeticError("the zero-gain rescue protocol is not frozen")
    source = protocol.get("source_campaign", {})
    if (
        source.get("sha256") != digest(CAMPAIGN)
        or source.get("candidate_list_sha256")
        != campaign.get("candidate_list_sha256")
        or source.get("candidate_count") != campaign.get("candidate_count")
    ):
        raise ArithmeticError("the rescue protocol names another source campaign")
    definition = {
        key: value
        for key, value in protocol.items()
        if key != "protocol_definition_sha256"
    }
    if protocol.get("protocol_definition_sha256") != canonical_hash(definition):
        raise ArithmeticError("the rescue protocol definition does not replay")
    validate_frozen_sources(protocol["implementation_hashes"])
    if [row["sample_id"] for row in protocol["assignments"]] != [
        row["sample_id"] for row in campaign["rows"]
    ]:
        raise ArithmeticError("the rescue assignment no longer covers the source rows in order")
    detector = protocol["rescue_detector"]
    if (
        detector["base_initial_chart_count"],
        detector["additional_budget_chart_count"],
        detector["generic_rescue_batch_count"],
        detector["generic_rescue_charts_per_batch"],
        detector["height_bound_each_chart"],
        detector["wall_timeout_seconds_each_chart"],
        detector["gp_stack_bytes_each_chart"],
        detector["relation_chunk_size"],
        detector["relation_timeout_seconds"],
    ) != (
        BASE_INITIAL_CHARTS,
        RESCUE_CHARTS,
        7,
        RESCUE_BATCH_SIZE,
        HEIGHT_BOUND,
        CHART_TIMEOUT_SECONDS,
        STACK_BYTES,
        RELATION_CHUNK_SIZE,
        RELATION_TIMEOUT_SECONDS,
    ):
        raise ArithmeticError("the rescue executable budget differs from the protocol")
    return protocol, campaign


def base_records(chunk_dir: Path) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    campaign_hash = digest(CAMPAIGN)
    for path in sorted(chunk_dir.glob("chunk-*-of-*.json")):
        chunk = json.loads(path.read_text())
        require_runtime(chunk)
        if (
            chunk.get("schema") != "elkies-k3.mw17-jump-v2-chunk.v1"
            or chunk.get("campaign_sha256") != campaign_hash
        ):
            continue
        for record in chunk.get("records", []):
            index = int(record["campaign_index"])
            if index in records:
                raise ArithmeticError("duplicate base campaign index across checkpoints")
            records[index] = record
    return records


def clean_zero_eligible(
    assignment: dict[str, Any], base_record: dict[str, Any] | None
) -> bool:
    if not assignment["assigned_to_rescue_arm"] or base_record is None:
        return False
    return (
        base_record.get("sample_id") == assignment["sample_id"]
        and base_record.get("status") == "PASS_EXACT_CERTIFIED_QUOTIENT_GAIN"
        and base_record.get("actual_certified_quotient_rank_gain") == 0
        and base_record.get("attempted_chart_count") == BASE_INITIAL_CHARTS
        and base_record.get("bounded_cover_timeout_count") == 0
        and base_record.get("cover_backend_failure_count") == 0
        and base_record.get("initial", {}).get(
            "initial_43_chart_certified_gain_diagnostic_only"
        )
        == 0
        and base_record.get("adaptive", {}).get("status")
        == "NOT_APPLICABLE_ZERO_DISCOVERED_QUOTIENT"
        and len(base_record.get("initial", {}).get("cover_records", []))
        == BASE_INITIAL_CHARTS
        and all(
            cover.get("search", {}).get("status") == "bounded_search_complete"
            for cover in base_record.get("initial", {}).get("cover_records", [])
        )
    )


def rank_generic_rescue_charts(ladder, legacy, model, generic, generic_rows):
    selected = generic_rows[BASE_INITIAL_CHARTS:TOTAL_CHART_CAP]
    if len(selected) != RESCUE_CHARTS:
        raise ArithmeticError("the disjoint generic rescue pool is not 301 classes")
    gram, asymmetry = legacy.canonical_height_gram(model, generic)
    runs = {}
    for scale in (ladder.AUDIT_SCALE, ladder.OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(gram, scale))
        rows = []
        maximum_error = 0.0
        for generic_norm, mask, generic_representative in selected:
            residue = tuple(
                (mask >> index) & 1 for index in range(GENERIC_DIMENSION)
            )
            unused_norm, representative, error = oracle.solve(residue)
            maximum_error = max(maximum_error, error)
            depth = legacy.quadratic_decimal(gram, representative) / 4
            rows.append(
                (depth, mask, representative, generic_representative, generic_norm)
            )
        rows.sort(key=lambda row: (-row[0], row[1]))
        runs[scale] = rows, maximum_error
    operative = runs[ladder.OPERATIVE_SCALE][0]
    audit = runs[ladder.AUDIT_SCALE][0]
    audit_map = {row[1]: row[2] for row in audit}
    return gram, operative, {
        "fixed_generic_class_ranks_one_based": [44, 344],
        "selected_mask_sha256": canonical_hash(
            [row[1] for row in generic_rows[BASE_INITIAL_CHARTS:TOTAL_CHART_CAP]]
        ),
        "canonical_height_maximum_asymmetry": str(asymmetry),
        "operative_rounding_scale": ladder.OPERATIVE_SCALE,
        "audit_rounding_scale": ladder.AUDIT_SCALE,
        "maximum_cvp_distance_error": {
            str(scale): error for scale, (unused_rows, error) in runs.items()
        },
        "representative_disagreement_count": sum(
            audit_map[row[1]] != row[2] for row in operative
        ),
        "priority_order_identical_between_scales": [row[1] for row in audit]
        == [row[1] for row in operative],
    }


class DetectorArgs:
    relation_chunk_size = RELATION_CHUNK_SIZE
    relation_timeout_seconds = RELATION_TIMEOUT_SECONDS
    stack_bytes = STACK_BYTES


def run_fibre(
    row: dict[str, Any],
    assignment: dict[str, Any],
    base_record: dict[str, Any],
    base,
) -> dict[str, Any]:
    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    families = base.Families()
    ladder, legacy, engine, chart_policy = base.load_modules()
    curve, known, generic_gram, specialization = families.specialize(row)
    (
        model,
        generic,
        search_model,
        search_generic,
        search_to_certificate,
        certificate_to_search,
        normalization,
    ) = base.normalize_curve(curve, known, row["family"])
    signatures = legacy.find_mod2_reduction_certificate(
        model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND
    )
    if legacy.combined_mod2_rank(signatures, GENERIC_DIMENSION) != GENERIC_DIMENSION:
        raise ArithmeticError("the specialized generic MW17 lost exact independence")
    if base_record["generic_subgroup"]["finite_reduction_independence_certificate"] != base.signature_record(
        legacy, signatures, GENERIC_DIMENSION
    ):
        raise ArithmeticError("the base zero record no longer matches the exact specialization")

    generic_rows, generic_cvp_error = base.complete_generic_census(
        legacy, generic_gram
    )
    rescue_gram, rescue_ranked, rescue_ranking = rank_generic_rescue_charts(
        ladder, legacy, model, generic, generic_rows
    )
    rescue_ids = [f"mask:{entry[1]:05x}" for entry in rescue_ranked]
    rescue_universe = f"mw17-jump-v2-rescue:{row['sample_id']}:generic-ranks44-344"
    rescue_certificate = chart_policy.bind_ordering(
        basis_records=[legacy.point_record(point) for point in generic],
        height_gram_rows=rescue_gram,
        generic_coordinate_rows=ladder.identity_rows(GENERIC_DIMENSION),
        quotient_coordinate_rows=[],
        chart_universe_id=rescue_universe,
        ordered_chart_ids=rescue_ids,
        heuristics=["disjoint_next_generic_half_classes", "specialized_depth"],
    )
    chart_policy.validate_ordering(
        rescue_certificate,
        basis_records=[legacy.point_record(point) for point in generic],
        height_gram_rows=rescue_gram,
        generic_coordinate_rows=ladder.identity_rows(GENERIC_DIMENSION),
        quotient_coordinate_rows=[],
        chart_universe_id=rescue_universe,
        ordered_chart_ids=rescue_ids,
    )

    search_budget = load_gate_module().production_gate_record(
        target_rank=32,
        search_limits={
            "additional_chart_count": RESCUE_CHARTS,
            "height_bound_each_chart": HEIGHT_BOUND,
            "wall_timeout_seconds_each_chart": CHART_TIMEOUT_SECONDS,
            "worker_wall_timeout_seconds": 7_200,
            "worker_address_space_bytes": 24_000_000_000,
        },
        scheduling_information=[
            {
                "kind": "residual_selmer_state",
                "status": "NOT_REQUIRED_FOR_BOUNDED_RESCUE",
            }
        ],
    )
    if not search_budget["search_budget_gate"]["bounded_search_authorized"]:
        raise ArithmeticError("the independent production budget gate is not open")

    discoveries = {}
    searched_keys = {
        cover["base_point_key"]
        for cover in base_record["initial"]["cover_records"]
    }
    basis = generic
    rescue_covers = []
    batch_classifications = []
    first_escape_batch = None
    for offset in range(0, RESCUE_CHARTS, RESCUE_BATCH_SIZE):
        batch = rescue_ranked[offset:offset + RESCUE_BATCH_SIZE]
        for within_batch, (
            depth,
            mask,
            representative,
            generic_representative,
            generic_norm,
        ) in enumerate(batch, 1):
            priority = offset + within_batch
            base_point = legacy.exact_linear_combination(
                model[3], generic, representative
            )
            if base_point is None:
                raise ArithmeticError("a generic rescue class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                raise ArithmeticError("the generic rescue repeated an earlier chart")
            outcome = base.run_quartic_search_raw(
                engine,
                mask=mask,
                representative=representative,
                short_model=search_model,
                generic_points=search_generic,
                height_bound=HEIGHT_BOUND,
                timeout_seconds=CHART_TIMEOUT_SECONDS,
                stack_bytes=STACK_BYTES,
            )
            searched_keys.add(base_key)
            source = f"rescue:priority:{priority}:mask:{mask:#07x}"
            for point in outcome.curve_points:
                certificate_point = search_to_certificate(point)
                discoveries.setdefault(
                    legacy.canonical_point(certificate_point), set()
                ).add(source)
            search_record = ladder.compact_search_record(outcome.record)
            search_record["error"] = outcome.record.get("error")
            rescue_covers.append(
                {
                    "priority": priority,
                    "batch": offset // RESCUE_BATCH_SIZE + 1,
                    "mask": mask,
                    "exact_generic_norm": generic_norm,
                    "generic_representative": list(generic_representative),
                    "specialized_representative": list(representative),
                    "specialized_depth": str(depth),
                    "base_point_key": base_key,
                    "search": search_record,
                }
            )
        basis, classification = legacy.classify_discovered_group(
            model=model,
            basis=generic,
            discoveries=discoveries,
            relation_chunk_size=RELATION_CHUNK_SIZE,
            relation_timeout_seconds=RELATION_TIMEOUT_SECONDS,
            stack_bytes=STACK_BYTES,
        )
        if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            raise ArithmeticError("rescue discoveries could not be exactly classified")
        gain = len(basis) - GENERIC_DIMENSION
        batch_classifications.append(
            {
                "batch": offset // RESCUE_BATCH_SIZE + 1,
                "cumulative_chart_count": len(rescue_covers),
                "exact_certified_quotient_gain": gain,
                "classification": classification,
            }
        )
        if gain:
            first_escape_batch = offset // RESCUE_BATCH_SIZE + 1
            break

    rescue_gain = len(basis) - GENERIC_DIMENSION
    remaining_slots = RESCUE_CHARTS - len(rescue_covers)
    adaptive = {
        "status": (
            "NOT_APPLICABLE_RESCUE_FOUND_NO_QUOTIENT"
            if not rescue_gain
            else "NO_UNUSED_BUDGET_AFTER_RESCUE_ESCAPE"
        ),
        "budget_chart_count": remaining_slots,
        "cover_records": [],
        "exact_incremental_quotient_rank_recovered": 0,
    }
    if rescue_gain and remaining_slots:
        (
            adaptive_gram,
            adaptive_ranked,
            generic_coordinates,
            complement,
            adaptive_ranking,
        ) = ladder.rank_adaptive_pool(
            legacy, model, basis, generic, generic_rows, DetectorArgs()
        )
        adaptive_ranked = adaptive_ranked[:remaining_slots]
        adaptive_ids = [
            f"gmask:{entry[1]:05x}:qword:{entry[2]:x}"
            for entry in adaptive_ranked
        ]
        adaptive_universe = (
            f"mw17-jump-v2-rescue:{row['sample_id']}:adaptive-after-batch-{first_escape_batch}"
        )
        adaptive_certificate = chart_policy.bind_ordering(
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=adaptive_universe,
            ordered_chart_ids=adaptive_ids,
            heuristics=["legacy_half_lattice_depth", "quotient_hamming_weight"],
        )
        chart_policy.validate_ordering(
            adaptive_certificate,
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=adaptive_universe,
            ordered_chart_ids=adaptive_ids,
        )
        adaptive_covers = []
        search_basis = tuple(certificate_to_search(point) for point in basis)
        for priority, (
            depth,
            generic_mask,
            quotient_word,
            residue,
            representative,
            generic_norm,
        ) in enumerate(adaptive_ranked, 1):
            base_point = legacy.exact_linear_combination(
                model[3], basis, representative
            )
            if base_point is None:
                raise ArithmeticError("an adaptive rescue class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                raise ArithmeticError("the rescue adaptive phase repeated a chart")
            mask = sum(int(bit) << index for index, bit in enumerate(residue))
            outcome = base.run_quartic_search_raw(
                engine,
                mask=mask,
                representative=representative,
                short_model=search_model,
                generic_points=search_basis,
                height_bound=HEIGHT_BOUND,
                timeout_seconds=CHART_TIMEOUT_SECONDS,
                stack_bytes=STACK_BYTES,
            )
            searched_keys.add(base_key)
            source = (
                f"adaptive-after-rescue:priority:{priority}:"
                f"gmask:{generic_mask:#07x}:qword:{quotient_word:#x}"
            )
            for point in outcome.curve_points:
                certificate_point = search_to_certificate(point)
                discoveries.setdefault(
                    legacy.canonical_point(certificate_point), set()
                ).add(source)
            search_record = ladder.compact_search_record(outcome.record)
            search_record["error"] = outcome.record.get("error")
            adaptive_covers.append(
                {
                    "priority": priority,
                    "generic_mask": generic_mask,
                    "exact_generic_norm": generic_norm,
                    "quotient_word": quotient_word,
                    "current_basis_residue": list(residue),
                    "representative": list(representative),
                    "canonical_depth": str(depth),
                    "base_point_key": base_key,
                    "search": search_record,
                }
            )
        pre_adaptive_gain = len(basis) - GENERIC_DIMENSION
        basis, final_classification = legacy.classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=RELATION_CHUNK_SIZE,
            relation_timeout_seconds=RELATION_TIMEOUT_SECONDS,
            stack_bytes=STACK_BYTES,
        )
        if final_classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            raise ArithmeticError("adaptive rescue discoveries were not exactly classified")
        adaptive = {
            "status": "PASS_EXACTLY_CLASSIFIED",
            "budget_chart_count": remaining_slots,
            "ranking": adaptive_ranking,
            "ordering_certificate": adaptive_certificate,
            "cover_records": adaptive_covers,
            "discovered_group_classification": final_classification,
            "exact_incremental_quotient_rank_recovered": (
                len(basis) - GENERIC_DIMENSION - pre_adaptive_gain
            ),
        }

    all_new_covers = rescue_covers + adaptive["cover_records"]
    if len(all_new_covers) != RESCUE_CHARTS:
        raise ArithmeticError("the rescue did not consume exactly its 301-chart allowance")
    backend_failures = sum(
        entry["search"]["status"] == "pari_failure" for entry in all_new_covers
    )
    timeouts = sum(
        entry["search"]["status"] == "bounded_search_timeout"
        for entry in all_new_covers
    )
    final_gain = len(basis) - GENERIC_DIMENSION
    exact_score = final_gain if final_gain or not (backend_failures or timeouts) else None
    final_basis_sha256 = canonical_hash(
        [legacy.point_record(point) for point in basis]
    )
    lower_bound = load_gate_module().certified_point_lower_bound_record(
        certified_independent_rank=GENERIC_DIMENSION + final_gain,
        target_rank=32,
        curve_equations_verified=True,
        independence_evidence=f"basis-sha256:{final_basis_sha256}",
    )
    result = {
        "sample_id": row["sample_id"],
        "campaign_index": row["campaign_index"],
        "source_population": row["source_population"],
        "family": row["family"],
        "parameter": row["parameter"],
        "assignment": assignment,
        "status": (
            "CENSORED_RESCUE_BACKEND_FAILURE"
            if backend_failures
            else "CENSORED_RESCUE_TIMEOUT"
            if timeouts
            else "PASS_EXACT_CERTIFIED_RESCUE_GAIN"
        ),
        "base_zero_measurement_sha256": canonical_hash(base_record),
        "specialization": specialization,
        "normalization": normalization,
        "production_gates": search_budget,
        "generic_census": {
            "complete_class_count": 1 << GENERIC_DIMENSION,
            "maximum_cvp_distance_error": generic_cvp_error,
            "base_mask_sha256": canonical_hash(
                [entry[1] for entry in generic_rows[:BASE_INITIAL_CHARTS]]
            ),
            "rescue_mask_sha256": rescue_ranking["selected_mask_sha256"],
            "base_and_rescue_masks_disjoint": True,
        },
        "generic_rescue": {
            "status": "PASS_BATCHED_EXACT_CLASSIFICATION",
            "ranking": rescue_ranking,
            "ordering_certificate": rescue_certificate,
            "first_escape_batch": first_escape_batch,
            "cover_records": rescue_covers,
            "batch_classifications": batch_classifications,
            "exact_quotient_rank_after_rescue_batches": rescue_gain,
        },
        "adaptive_after_rescue": adaptive,
        "actual_certified_quotient_rank_gain": exact_score,
        "certified_rank_lower_bound": (
            GENERIC_DIMENSION + final_gain if exact_score is not None else None
        ),
        "point_lower_bound_certificate": lower_bound,
        "base_chart_count": BASE_INITIAL_CHARTS,
        "additional_attempted_chart_count": len(all_new_covers),
        "total_attempted_chart_count_including_base": (
            BASE_INITIAL_CHARTS + len(all_new_covers)
        ),
        "bounded_cover_timeout_count": timeouts,
        "cover_backend_failure_count": backend_failures,
        "final_basis": [legacy.point_record(point) for point in basis],
        "timing": {
            "worker_wall_seconds": time.monotonic() - started_wall,
            "worker_cpu_seconds_parent_plus_children": cpu_clock() - started_cpu,
        },
        "claim_boundary": (
            "The exact score is a certified point-based rank lower bound. A rescue "
            "miss or censored chart is not a rank or Selmer upper bound."
        ),
    }
    if result["total_attempted_chart_count_including_base"] != TOTAL_CHART_CAP:
        raise ArithmeticError("the rescue exceeded or underspent its fixed total cap")
    if final_gain >= 15:
        atomic_write(
            STOP_SENTINEL,
            {
                "schema": "elkies-k3.mw17-jump-v2-zero-gain-rescue-stop.v1",
                "status": "CERTIFIED_GAIN_AT_LEAST_15_GLOBAL_STOP",
                "sample_id": row["sample_id"],
                "campaign_index": row["campaign_index"],
                "actual_certified_quotient_rank_gain": final_gain,
                "certified_rank_lower_bound": GENERIC_DIMENSION + final_gain,
                "protocol_sha256": digest(PROTOCOL),
                "runtime_search": runtime_search(),
            },
        )
    return result


def failure_record(row, status: str, failure: Any, elapsed: float):
    return {
        "sample_id": row["sample_id"],
        "campaign_index": row["campaign_index"],
        "source_population": row["source_population"],
        "family": row["family"],
        "parameter": row["parameter"],
        "status": status,
        "failure": failure,
        "actual_certified_quotient_rank_gain": None,
        "certified_rank_lower_bound": None,
        "supervisor_wall_seconds": elapsed,
    }


def run_single(index: int, base_chunk_dir: Path) -> None:
    protocol, campaign = load_protocol()
    if not 0 <= index < len(campaign["rows"]):
        raise ValueError("single index is outside the frozen campaign")
    assignment = protocol["assignments"][index]
    base_record = base_records(base_chunk_dir).get(index)
    if not clean_zero_eligible(assignment, base_record):
        raise ArithmeticError("the requested fibre is not an assigned clean-zero rescue")
    memory = protocol["resource_gate"]["worker_address_space_bytes"]
    resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    result = run_fibre(
        campaign["rows"][index], assignment, base_record, load_base_module()
    )
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def write_checkpoint(
    path: Path,
    chunk_index: int,
    chunk_count: int,
    assigned_indices: list[int],
    records: list[dict[str, Any]],
) -> None:
    atomic_write(
        path,
        {
            "schema": "elkies-k3.mw17-jump-v2-zero-gain-rescue-chunk.v1",
            "status": "ELIGIBILITY_DEPENDENT_CHECKPOINT",
            "protocol_sha256": digest(PROTOCOL),
            "runtime_search": runtime_search(),
            "source_campaign_sha256": digest(CAMPAIGN),
            "chunk_index": chunk_index,
            "chunk_count": chunk_count,
            "hash_assigned_indices": assigned_indices,
            "completed_record_count": len(records),
            "records": records,
        },
    )


def run_chunk(
    chunk_index: int,
    chunk_count: int,
    output: Path,
    base_chunk_dir: Path,
    max_new: int | None,
) -> None:
    protocol, campaign = load_protocol()
    if not 0 <= chunk_index < chunk_count:
        raise ValueError("chunk index must lie in [0, chunk count)")
    indices = [
        row["campaign_index"]
        for row in protocol["assignments"]
        if row["assigned_to_rescue_arm"]
        and row["campaign_index"] % chunk_count == chunk_index
    ]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        require_runtime(old)
        if (
            old.get("protocol_sha256") != digest(PROTOCOL)
            or old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("hash_assigned_indices") != indices
        ):
            raise ArithmeticError("an existing rescue checkpoint belongs elsewhere")
        records = old["records"]
    completed_ids = {record["sample_id"] for record in records}
    current_base = base_records(base_chunk_dir)
    new_count = 0
    timeout = protocol["resource_gate"]["worker_wall_timeout_seconds"]
    for position, index in enumerate(indices, 1):
        row = campaign["rows"][index]
        if row["sample_id"] in completed_ids:
            continue
        if not clean_zero_eligible(protocol["assignments"][index], current_base.get(index)):
            continue
        if STOP_SENTINEL.exists():
            print(
                f"MW17JUMPV2RESCUECHUNK|chunk={chunk_index}/{chunk_count}|status=STOP_GAIN15",
                flush=True,
            )
            break
        if max_new is not None and new_count >= max_new:
            break
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single-index",
            str(index),
            "--base-chunk-dir",
            str(base_chunk_dir),
        ]
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
            lines = [line for line in completed.stdout.splitlines() if line.strip()]
            result_line = next(
                (line for line in reversed(lines) if line.startswith("RESULT_JSON=")),
                None,
            )
            if completed.returncode == 0 and result_line is not None:
                result = json.loads(result_line[len("RESULT_JSON="):])
                result["supervisor_wall_seconds"] = time.monotonic() - started
            else:
                result = failure_record(
                    row,
                    "CENSORED_RESCUE_WORKER_FAILURE",
                    {"returncode": completed.returncode, "output_tail": lines[-40:]},
                    time.monotonic() - started,
                )
        except subprocess.TimeoutExpired as error:
            output_text = error.stdout or ""
            if isinstance(output_text, bytes):
                output_text = output_text.decode(errors="replace")
            result = failure_record(
                row,
                "CENSORED_RESCUE_WORKER_TIMEOUT",
                {"output_tail": output_text.splitlines()[-40:]},
                time.monotonic() - started,
            )
        records.append(result)
        completed_ids.add(row["sample_id"])
        new_count += 1
        write_checkpoint(output, chunk_index, chunk_count, indices, records)
        print(
            f"MW17JUMPV2RESCUECHUNK|chunk={chunk_index}/{chunk_count}"
            f"|position={position}/{len(indices)}|sample={row['sample_id']}"
            f"|status={result['status']}"
            f"|gain={result.get('actual_certified_quotient_rank_gain')}",
            flush=True,
        )
    write_checkpoint(output, chunk_index, chunk_count, indices, records)


def merge_chunks(chunk_dir: Path, chunk_count: int, output: Path) -> None:
    protocol, _campaign = load_protocol()
    records_by_index = {}
    provenance = []
    for chunk_index in range(chunk_count):
        path = chunk_dir / f"chunk-{chunk_index:02d}-of-{chunk_count:02d}.json"
        if not path.exists():
            continue
        chunk = json.loads(path.read_text())
        require_runtime(chunk)
        if (
            chunk.get("schema")
            != "elkies-k3.mw17-jump-v2-zero-gain-rescue-chunk.v1"
            or chunk.get("protocol_sha256") != digest(PROTOCOL)
            or chunk.get("chunk_index") != chunk_index
            or chunk.get("chunk_count") != chunk_count
        ):
            raise ArithmeticError(f"rescue chunk {chunk_index} is malformed")
        for record in chunk["records"]:
            index = int(record["campaign_index"])
            if index in records_by_index:
                raise ArithmeticError("duplicate rescue index across chunks")
            records_by_index[index] = record
        provenance.append(
            {
                "path": relative(path),
                "sha256": digest(path),
                "record_count": len(chunk["records"]),
            }
        )
    scored = [
        record
        for record in records_by_index.values()
        if record.get("actual_certified_quotient_rank_gain") is not None
    ]
    leaderboard = sorted(
        scored,
        key=lambda record: (
            -record["actual_certified_quotient_rank_gain"],
            record["campaign_index"],
        ),
    )
    measurements = [
        {
            "leaderboard_position": position,
            "campaign_index": record["campaign_index"],
            "sample_id": record["sample_id"],
            "source_population": record["source_population"],
            "family": record["family"],
            "parameter": record["parameter"],
            "first_escape_batch": record["generic_rescue"]["first_escape_batch"],
            "actual_certified_quotient_rank_gain": record[
                "actual_certified_quotient_rank_gain"
            ],
            "certified_rank_lower_bound": record["certified_rank_lower_bound"],
            "detailed_chunk_record_sha256": canonical_hash(record),
        }
        for position, record in enumerate(leaderboard, 1)
    ]
    document = {
        "schema": "elkies-k3.mw17-jump-v2-zero-gain-rescue-ledger.v1",
        "status": "PARTIAL_ELIGIBILITY_DEPENDENT_RESCUE_LEDGER",
        "protocol_sha256": digest(PROTOCOL),
        "runtime_search": runtime_search(),
        "source_campaign_sha256": digest(CAMPAIGN),
        "hash_assigned_candidate_count": protocol["assignment"][
            "assigned_candidate_count"
        ],
        "completed_rescue_worker_count": len(records_by_index),
        "measured_certified_gain_count": len(measurements),
        "status_counts": dict(
            sorted(Counter(row["status"] for row in records_by_index.values()).items())
        ),
        "measurements_ranked_only_by_actual_certified_quotient_rank_gain": measurements,
        "gain_at_least_15_found": any(
            row["actual_certified_quotient_rank_gain"] >= 15
            for row in measurements
        ),
        "chunk_provenance": provenance,
        "claim_boundary": protocol["claim_boundary"],
        "generation": {
            "script": relative(Path(__file__)),
            "python": platform.python_version(),
            "command": (
                f"sage -python {relative(Path(__file__))} --merge "
                f"--chunk-count {chunk_count}"
            ),
        },
    }
    atomic_write(output, document)
    print(
        f"MW17JUMPV2RESCUEMERGE|rescued={len(records_by_index)}"
        f"|measured={len(measurements)}|output={relative(output)}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=4)
    parser.add_argument("--base-chunk-dir", type=Path, default=BASE_CHUNK_DIR)
    parser.add_argument("--chunk-dir", type=Path, default=CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-new", type=int)
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP is required")
    if args.check:
        protocol, campaign = load_protocol()
        load_base_module().load_campaign()
        load_gate_module()
        print(
            "MW17JUMPV2RESCUE|status=PASS_PREFLIGHT"
            f"|candidates={campaign['candidate_count']}"
            f"|assigned={protocol['assignment']['assigned_candidate_count']}"
        )
        return
    if args.single_index is not None:
        run_single(args.single_index, args.base_chunk_dir.resolve())
        return
    if args.merge:
        merge_chunks(
            args.chunk_dir.resolve(), args.chunk_count, (args.output or LEDGER).resolve()
        )
        return
    if args.chunk_index is None:
        raise SystemExit("choose --single-index, --chunk-index, --merge, or --check")
    output = args.output or args.chunk_dir / (
        f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
    )
    run_chunk(
        args.chunk_index,
        args.chunk_count,
        output.resolve(),
        args.base_chunk_dir.resolve(),
        args.max_new,
    )


if __name__ == "__main__":
    main()
