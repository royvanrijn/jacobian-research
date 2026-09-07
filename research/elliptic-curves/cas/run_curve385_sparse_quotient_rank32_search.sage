#!/usr/bin/env sage -python
"""Checkpointed sparse quotient-mask search for rank at least 32 on curve 385.

The input is the frozen blind M29 ledger.  Each lattice state searches complete
predeclared sparse quotient-word stages.  Any exactly certified rank or finite-
index enlargement restarts at weight one in the enlarged lattice.  A miss at a
sparse stage is never reported as a rank upper bound.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import platform
import shutil
import sys
import time
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local"
BLIND = ART / "curve385_iterated_half_lattice_blind_v1.json"
PROFILE = ART / "curve385_quotient_weight_profile_v1.json"
PROTOCOL = ART / "curve385_sparse_quotient_rank32_protocol_v1.json"
POLICY_SOURCE = CAS / "curve385_sparse_quotient_policy.py"
LEGACY_SOURCE = CAS / "run_curve385_iterated_half_lattice_search.sage"
OUTPUT = ROOT / "artifacts/local/elliptic-curves/pointed-quartic-search/campaigns/run_curve385_sparse_quotient_rank32_search.json"

EXPECTED_PROTOCOL_DEFINITION_HASH = "5723679da2907e036095f90376cdabde457a4f7ba5bc284ad4a4ca3edea1aa37"
TARGET_RANK = 32
AUDIT_SCALE = 100_000
OPERATIVE_SCALE = 1_000_000

sys.path.insert(0, str(CAS))
from curve385_sparse_quotient_policy import (  # noqa: E402
    canonical_hash as policy_hash,
    stage_plan,
    validate_stage_plan,
)

legacy = SourceFileLoader("curve385_sparse_legacy", str(LEGACY_SOURCE)).load_module()

Point = tuple[Fraction, Fraction]


from pointed_quartic_migration import runtime_search, require_runtime, validate_frozen_sources

def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(payload))


def read_discoveries(rows: Iterable[dict[str, Any]]) -> dict[Point, set[str]]:
    answer: dict[Point, set[str]] = {}
    for row in rows:
        point = legacy.canonical_point(legacy.read_point(row["point"]))
        answer.setdefault(point, set()).update(map(str, row["sources"]))
    return answer


def load_protocol() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text())
    if protocol.get("status") != "FROZEN_BEFORE_SPARSE_RANK32_POINT_SEARCH":
        raise ArithmeticError("the sparse rank-32 protocol is not frozen")
    definition_hash = policy_hash(protocol["protocol_definition"])
    if definition_hash != protocol.get("protocol_definition_hash"):
        raise ArithmeticError("the sparse protocol definition hash is invalid")
    if definition_hash != EXPECTED_PROTOCOL_DEFINITION_HASH:
        raise ArithmeticError("the sparse protocol definition differs from the runner")
    for path in (BLIND, PROFILE, POLICY_SOURCE, LEGACY_SOURCE, Path(__file__)):
        expected = protocol["input_hashes"].get(relative(path))
        if expected is None:
            raise ArithmeticError(f"protocol omitted input: {relative(path)}")
        validate_frozen_sources({relative(path): expected})
    for bit_count_text, plan in protocol["stage_plans_by_quotient_bit_count"].items():
        validate_stage_plan(plan, int(bit_count_text), protocol["old_class_count"])
    return protocol


def initial_state(protocol: dict[str, Any], args) -> tuple[
    dict[str, Any], tuple[Point, ...], tuple[Point, ...], tuple[int, ...],
    dict[Point, set[str]], set[str]
]:
    blind = json.loads(BLIND.read_text())
    if blind.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("the source blind ledger is not frozen at M29")
    basis = tuple(legacy.read_point(row) for row in blind["current_basis"])
    generic = tuple(legacy.read_point(row) for row in blind["curve"]["generic_points"])
    old_masks = tuple(map(int, blind["old_deep43"]["masks"]))
    discoveries = read_discoveries(blind["discoveries"])
    searched_keys = set(map(str, blind["searched_base_point_keys"]))
    if len(basis) != protocol["starting_rank"] or len(generic) != 17:
        raise ArithmeticError("the frozen sparse-search starting lattice changed")
    if legacy.canonical_hash([legacy.point_record(row) for row in basis]) != protocol[
        "starting_basis_sha256"
    ]:
        raise ArithmeticError("the frozen M29 basis hash changed")
    payload = {
        "schema": "elliptic-curves.curve385-sparse-quotient-rank32-search.v1",
        "status": "PARTIAL_CHECKPOINT",
        "runtime_search": runtime_search(),
        "protocol": {
            "path": relative(PROTOCOL),
            "whole_file_sha256": digest(PROTOCOL),
            "definition_sha256": protocol["protocol_definition_hash"],
        },
        "configuration": {
            "height_bound_each_quartic": args.height_bound,
            "wall_timeout_seconds_each_quartic": args.timeout_seconds,
            "gp_stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "maximum_stage_each_lattice_state": args.max_stage,
            "maximum_lattice_states": args.max_lattice_states,
            "checkpoint_every_completed_searches": args.checkpoint_every,
            "retries": 0,
        },
        "curve": blind["curve"],
        "old_deep43": blind["old_deep43"],
        "source_blind_ledger_sha256": digest(BLIND),
        "current_basis": [legacy.point_record(row) for row in basis],
        "discoveries": legacy.discovery_records(discoveries),
        "searched_base_point_keys": sorted(searched_keys),
        "lattice_states": [],
        "generation": {
            "python": platform.python_version(),
            "sage": str(sys.modules.get("sage")),
            "pari": str(legacy.pari("default(parisizemax)")),
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    return payload, basis, generic, old_masks, discoveries, searched_keys


def resumed_state(protocol: dict[str, Any], args) -> tuple[
    dict[str, Any], tuple[Point, ...], tuple[Point, ...], tuple[int, ...],
    dict[Point, set[str]], set[str]
]:
    payload = json.loads(args.output.read_text())
    require_runtime(payload)
    if payload.get("schema") != "elliptic-curves.curve385-sparse-quotient-rank32-search.v1":
        raise ArithmeticError("the checkpoint has the wrong schema")
    if payload.get("status") != "PARTIAL_CHECKPOINT":
        raise ArithmeticError(f"checkpoint is terminal: {payload.get('status')}")
    if payload["protocol"]["whole_file_sha256"] != digest(PROTOCOL):
        raise ArithmeticError("the checkpoint uses a different protocol file")
    expected_configuration = {
        "height_bound_each_quartic": args.height_bound,
        "wall_timeout_seconds_each_quartic": args.timeout_seconds,
        "gp_stack_bytes_each_quartic": args.stack_bytes,
        "relation_chunk_size": args.relation_chunk_size,
        "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
        "maximum_stage_each_lattice_state": args.max_stage,
        "maximum_lattice_states": args.max_lattice_states,
        "checkpoint_every_completed_searches": args.checkpoint_every,
        "retries": 0,
    }
    if payload["configuration"] != expected_configuration:
        raise ArithmeticError("resume arguments differ from the checkpoint")
    basis = tuple(legacy.read_point(row) for row in payload["current_basis"])
    generic = tuple(legacy.read_point(row) for row in payload["curve"]["generic_points"])
    old_masks = tuple(map(int, payload["old_deep43"]["masks"]))
    discoveries = read_discoveries(payload["discoveries"])
    searched_keys = set(map(str, payload["searched_base_point_keys"]))
    return payload, basis, generic, old_masks, discoveries, searched_keys


def quotient_residue(
    old_mask: int,
    physical_word: int,
    generic_coordinates: Sequence[Sequence[int]],
    complement: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    residue = (0,) * len(generic_coordinates[0])
    for index, row in enumerate(generic_coordinates):
        if (old_mask >> index) & 1:
            residue = legacy.add_mod2(residue, row)
    for index, row in enumerate(complement):
        if (physical_word >> index) & 1:
            residue = legacy.add_mod2(residue, row)
    return residue


def rank_sparse_lifts(
    model,
    basis: Sequence[Point],
    generic: Sequence[Point],
    old_masks: Sequence[int],
    physical_words: Sequence[int],
    args,
):
    height_started = time.monotonic()
    gram, asymmetry = legacy.canonical_height_gram(model, basis)
    height_seconds = time.monotonic() - height_started
    args.generic_points = generic
    generic_coordinates, complement = legacy.gf2_lift_data(model, basis, generic, args)
    quotient_bits = len(complement)
    if any(not 0 < word < (1 << quotient_bits) for word in physical_words):
        raise ArithmeticError("a sparse physical word is outside the quotient")
    residues = [
        (
            old_mask,
            word,
            quotient_residue(old_mask, word, generic_coordinates, complement),
        )
        for old_mask in old_masks
        for word in physical_words
    ]
    if len({row[2] for row in residues}) != len(residues):
        raise ArithmeticError("sparse lifted parity classes collided")
    runs = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(gram, scale))
        rows = []
        maximum_error = 0.0
        for old_mask, word, residue in residues:
            unused_norm, representative, error = oracle.solve(residue)
            depth = legacy.quadratic_decimal(gram, representative) / Decimal(4)
            rows.append((depth, old_mask, word, residue, representative))
            maximum_error = max(maximum_error, error)
        rows.sort(key=lambda row: (-row[0], row[1], row[2]))
        runs[scale] = rows, maximum_error
    operative = runs[OPERATIVE_SCALE][0]
    audit = runs[AUDIT_SCALE][0]
    audit_representatives = {(row[1], row[2]): row[4] for row in audit}
    ranking = {
        "canonical_height_seconds": height_seconds,
        "canonical_height_maximum_asymmetry": str(asymmetry),
        "operative_rounding_scale": OPERATIVE_SCALE,
        "audit_rounding_scale": AUDIT_SCALE,
        "maximum_cvp_distance_error": {
            str(scale): error for scale, (unused_rows, error) in runs.items()
        },
        "representative_disagreement_count": sum(
            audit_representatives[(row[1], row[2])] != row[4] for row in operative
        ),
        "priority_order_identical_between_scales": [
            (row[1], row[2]) for row in audit
        ] == [(row[1], row[2]) for row in operative],
        "generic_coordinate_rows_in_current_basis": [
            list(row) for row in generic_coordinates
        ],
        "quotient_complement_rows_mod2": [list(row) for row in complement],
        "quotient_bit_count": quotient_bits,
        "physical_word_count": len(physical_words),
        "ranked_lift_count": len(operative),
        "ranked_lifts_sha256": legacy.canonical_hash(
            [
                {
                    "old_mask": row[1],
                    "physical_word": row[2],
                    "residue": list(row[3]),
                    "representative": list(row[4]),
                    "depth": str(row[0]),
                }
                for row in operative
            ]
        ),
    }
    return operative, ranking


def checkpoint_state(
    args,
    payload: dict[str, Any],
    basis: Sequence[Point],
    discoveries: dict[Point, set[str]],
    searched_keys: set[str],
) -> None:
    payload["current_basis"] = [legacy.point_record(row) for row in basis]
    payload["discoveries"] = legacy.discovery_records(discoveries)
    payload["searched_base_point_keys"] = sorted(searched_keys)
    write_payload(args.output, payload)


def stage_counts(stage: dict[str, Any]) -> None:
    records = stage["cover_records"]
    stage["searched_new_chart_count"] = len(records)
    stage["bounded_complete_count"] = sum(
        row["search"]["status"] == "bounded_search_complete" for row in records
    )
    stage["timeout_count"] = sum(
        row["search"]["status"] == "bounded_search_timeout" for row in records
    )
    stage["pari_failure_count"] = sum(
        row["search"]["status"] == "pari_failure" for row in records
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--max-stage", type=int, default=2, choices=range(1, 7))
    parser.add_argument("--max-lattice-states", type=int, default=4)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if args.max_lattice_states < 1 or args.checkpoint_every < 1:
        raise SystemExit("invalid checkpoint/restart budget")
    protocol = load_protocol()
    if args.plan_only:
        for bit_count, plan in protocol["stage_plans_by_quotient_bit_count"].items():
            counts = ",".join(
                f"{row['index']}:{row['new_chart_count']}" for row in plan[: args.max_stage]
            )
            print(f"C385SPARSE|plan|bits={bit_count}|stages={counts}")
        return
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")
    if args.resume:
        if not args.output.exists():
            raise SystemExit(f"resume checkpoint does not exist: {args.output}")
        state = resumed_state(protocol, args)
    else:
        if args.output.exists():
            raise SystemExit(f"refusing to overwrite existing search ledger: {args.output}")
        state = initial_state(protocol, args)
    payload, basis, generic, old_masks, discoveries, searched_keys = state
    model = tuple(Fraction(value) for value in payload["curve"]["short_model"])
    checkpoint_state(args, payload, basis, discoveries, searched_keys)

    while True:
        if len(basis) >= TARGET_RANK:
            payload["status"] = "PASS_RANK_AT_LEAST_32_DISCOVERY"
            payload["certified_rank_lower_bound"] = len(basis)
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            print(f"C385SPARSE|status=RANK32|rank={len(basis)}")
            return

        active_state = payload["lattice_states"][-1] if payload["lattice_states"] else None
        if active_state is None or active_state["status"] == "GROUP_CHANGED_RESTART":
            if len(payload["lattice_states"]) >= args.max_lattice_states:
                payload["status"] = "STOPPED_AT_DECLARED_LATTICE_STATE_LIMIT"
                payload["stop"] = {
                    "basis_rank": len(basis),
                    "lattice_state_count": len(payload["lattice_states"]),
                }
                checkpoint_state(args, payload, basis, discoveries, searched_keys)
                return
            quotient_bits = len(basis) - len(generic)
            plans = protocol["stage_plans_by_quotient_bit_count"].get(str(quotient_bits))
            if plans is None:
                raise ArithmeticError("the current quotient dimension is outside the protocol")
            active_state = {
                "index": len(payload["lattice_states"]) + 1,
                "status": "ACTIVE",
                "basis_rank": len(basis),
                "basis": [legacy.point_record(row) for row in basis],
                "basis_sha256": legacy.canonical_hash(
                    [legacy.point_record(row) for row in basis]
                ),
                "quotient_bit_count": quotient_bits,
                "stages": [],
            }
            payload["lattice_states"].append(active_state)
            checkpoint_state(args, payload, basis, discoveries, searched_keys)

        plans = protocol["stage_plans_by_quotient_bit_count"][
            str(active_state["quotient_bit_count"])
        ]
        active_stage = active_state["stages"][-1] if active_state["stages"] else None
        ranked = None
        if active_stage is not None and active_stage["status"] == "SEARCHING":
            stage_index = active_stage["index"]
        else:
            stage_index = 1 if active_stage is None else active_stage["index"] + 1
            if stage_index > args.max_stage:
                active_state["status"] = "DECLARED_STAGE_LIMIT_REACHED_WITHOUT_GROWTH"
                payload["status"] = "STOPPED_AFTER_DECLARED_SPARSE_STAGE_LIMIT"
                payload["stop"] = {
                    "basis_rank": len(basis),
                    "maximum_stage": args.max_stage,
                    "no_rank_upper_bound_claimed": True,
                }
                checkpoint_state(args, payload, basis, discoveries, searched_keys)
                print(
                    f"C385SPARSE|status=STAGE_LIMIT|rank={len(basis)}|stage={args.max_stage}"
                )
                return
            plan = plans[stage_index - 1]
            ranked, ranking = rank_sparse_lifts(
                model,
                basis,
                generic,
                old_masks,
                tuple(map(int, plan["new_physical_words"])),
                args,
            )
            if len(ranked) != plan["new_chart_count"]:
                raise ArithmeticError("ranked sparse lift count differs from protocol")
            active_stage = {
                "index": stage_index,
                "id": plan["id"],
                "status": "SEARCHING",
                "policy": plan,
                "ranking": ranking,
                "cover_records": [],
                "unchanged_previously_searched_chart_count": 0,
                "unchanged_previously_searched_chart_keys": [],
                "wall_started_unix": time.time(),
                "cpu_seconds_completed_searches": 0.0,
            }
            active_state["stages"].append(active_stage)
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            print(
                f"C385SPARSE|state={active_state['index']}|rank={len(basis)}|"
                f"stage={stage_index}|planned={len(ranked)}|status=START",
                flush=True,
            )

        plan = plans[stage_index - 1]
        if ranked is None:
            ranked, ranking = rank_sparse_lifts(
                model,
                basis,
                generic,
                old_masks,
                tuple(map(int, plan["new_physical_words"])),
                args,
            )
        if ranking["ranked_lifts_sha256"] != active_stage["ranking"]["ranked_lifts_sha256"]:
            raise ArithmeticError("sparse priority order changed while resuming")
        completed_pairs = {
            (int(row["old_mask"]), int(row["physical_quotient_word"]))
            for row in active_stage["cover_records"]
        }
        skipped_keys = set(active_stage["unchanged_previously_searched_chart_keys"])
        completed_since_checkpoint = 0
        for priority, (depth, old_mask, word, residue, representative) in enumerate(ranked, 1):
            if (old_mask, word) in completed_pairs:
                continue
            base_point = legacy.exact_linear_combination(model[3], basis, representative)
            if base_point is None:
                raise ArithmeticError("a sparse nonzero class produced the point at infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                if base_key not in skipped_keys:
                    active_stage["unchanged_previously_searched_chart_count"] += 1
                    active_stage["unchanged_previously_searched_chart_keys"].append(base_key)
                    skipped_keys.add(base_key)
                continue
            started_cpu = legacy.cpu_clock()
            outcome = legacy.engine.run_quartic_search(
                mask=sum(int(bit) << index for index, bit in enumerate(residue)),
                representative=representative,
                short_model=model,
                generic_points=basis,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            cpu_seconds = legacy.cpu_clock() - started_cpu
            searched_keys.add(base_key)
            source = (
                f"sparse:state:{active_state['index']}:stage:{stage_index}:"
                f"old:{old_mask}:q:{word}"
            )
            for point in outcome.curve_points:
                point = legacy.canonical_point(point)
                discoveries.setdefault(point, set()).add(source)
            active_stage["cover_records"].append(
                {
                    "priority": priority,
                    "old_mask": old_mask,
                    "old_hex": f"0x{old_mask:05x}",
                    "physical_quotient_word": word,
                    "physical_quotient_word_binary": (
                        f"{word:0{active_state['quotient_bit_count']}b}"
                    ),
                    "current_basis_residue": list(residue),
                    "canonical_depth": str(depth),
                    "representative": list(representative),
                    "base_point_key": base_key,
                    "cpu_seconds": cpu_seconds,
                    "search": outcome.record,
                }
            )
            active_stage["cpu_seconds_completed_searches"] += cpu_seconds
            completed_since_checkpoint += 1
            if completed_since_checkpoint >= args.checkpoint_every:
                checkpoint_state(args, payload, basis, discoveries, searched_keys)
                completed_since_checkpoint = 0
            print(
                f"C385SPARSE|state={active_state['index']}|stage={stage_index}|"
                f"priority={priority}/{len(ranked)}|old={old_mask:#07x}|q={word:#x}|"
                f"status={outcome.record['status']}|points={len(outcome.curve_points)}",
                flush=True,
            )

        stage_counts(active_stage)
        old_basis_hash = legacy.canonical_hash([legacy.point_record(row) for row in basis])
        basis, saturation = legacy.classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=args.relation_chunk_size,
            relation_timeout_seconds=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        active_stage["discovered_group_saturation"] = saturation
        active_stage["basis_rank_after"] = len(basis)
        active_stage["basis_after"] = [legacy.point_record(row) for row in basis]
        active_stage["basis_after_sha256"] = legacy.canonical_hash(
            [legacy.point_record(row) for row in basis]
        )
        active_stage["group_changed"] = active_stage["basis_after_sha256"] != old_basis_hash
        active_stage["new_independent_direction_count"] = sum(
            event["type"] == "NEW_Q_INDEPENDENT_DIRECTION"
            for event in saturation["events"]
        )
        active_stage["finite_index_saturation_event_count"] = sum(
            event["type"] == "FINITE_INDEX_SATURATION_INSIDE_DISCOVERED_GROUP"
            for event in saturation["events"]
        )
        active_stage["wall_seconds"] = time.time() - active_stage["wall_started_unix"]
        if saturation["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            active_stage["status"] = "UNKNOWN_UNCLASSIFIED_DISCOVERIES"
            payload["status"] = "STOPPED_FAIL_CLOSED_UNCLASSIFIED_DISCOVERIES"
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            return
        if active_stage["timeout_count"] or active_stage["pari_failure_count"]:
            active_stage["status"] = "INCOMPLETE_SEARCH_STAGE"
            payload["status"] = "STOPPED_AT_INCOMPLETE_SPARSE_STAGE"
            payload["stop"] = {
                "basis_rank": len(basis),
                "timeouts": active_stage["timeout_count"],
                "pari_failures": active_stage["pari_failure_count"],
                "no_retry_and_no_rank_upper_bound_claimed": True,
            }
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            return
        active_stage["status"] = "CLASSIFIED"
        if len(basis) >= TARGET_RANK:
            active_state["status"] = "TARGET_REACHED"
            payload["status"] = "PASS_RANK_AT_LEAST_32_DISCOVERY"
            payload["certified_rank_lower_bound"] = len(basis)
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            print(f"C385SPARSE|status=RANK32|rank={len(basis)}", flush=True)
            return
        if active_stage["group_changed"]:
            active_state["status"] = "GROUP_CHANGED_RESTART"
            checkpoint_state(args, payload, basis, discoveries, searched_keys)
            print(
                f"C385SPARSE|state={active_state['index']}|stage={stage_index}|"
                f"status=GROW|rank={active_state['basis_rank']}->{len(basis)}",
                flush=True,
            )
            continue
        checkpoint_state(args, payload, basis, discoveries, searched_keys)
        print(
            f"C385SPARSE|state={active_state['index']}|stage={stage_index}|"
            f"status=NO_GROWTH|rank={len(basis)}",
            flush=True,
        )


if __name__ == "__main__":
    main()
