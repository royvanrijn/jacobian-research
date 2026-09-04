#!/usr/bin/env sage-python
"""Blind MW16 half-lattice and adaptive-quotient search on curve 398.

The sole mathematical input is the redacted curve-398 fixture: a short curve,
sixteen specialized generic A1/MW16 sections, and their exact generic height
Gram.  This runner never imports the public 30-point list, the recovered
fibration artifact, its public embedding, or any held-out point coordinates.

It completely enumerates M/2M for the generic half-integral MW16 lattice,
searches all deepest classes, classifies every returned point by exact group
law, and then applies the curve-385 adaptive quotient lift rule to every
nonzero quotient word within the declared bound.  Search misses remain bounded.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import platform
import resource
import shutil
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_curve398_mw16_blind_input_v1.json"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/curve398_mw16_adaptive_half_lattice_blind_v1.json"
EXPECTED_INPUT_SHA256 = "b5c41d27cdf5ce3707de29eec1febceec15cf599b7a5110f75be3dd0594fc0fc"
GENERIC_DIMENSION = 16
EXPECTED_DEEPEST_COUNT = 12
EXPECTED_DEEPEST_TWICE_NORM = 23
OPERATIVE_SCALE = 1_000_000
AUDIT_SCALE = 100_000


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def write_payload(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--max-quotient-bits", type=int, default=5)
    parser.add_argument("--max-planned-lifts", type=int, default=372)
    parser.add_argument("--max-iterations", type=int, default=4)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if args.max_quotient_bits < 1 or args.max_planned_lifts < EXPECTED_DEEPEST_COUNT:
        raise SystemExit("invalid adaptive-lift bound")
    if args.max_iterations < 1 or args.checkpoint_every < 1:
        raise SystemExit("invalid iteration/checkpoint bound")
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")
    if digest(INPUT) != EXPECTED_INPUT_SHA256:
        raise ArithmeticError("redacted curve-398 blind input changed")

    blind = json.loads(INPUT.read_text())
    if blind.get("status") != "PASS_REDACTED_GENERIC_MW16_INPUT":
        raise ArithmeticError("redacted input is not certified")
    if blind["redaction"] != {
        "contains_public_embedding_coordinates": False,
        "held_out_point_count": 14,
        "permitted_information": "curve equation, sixteen specialized generic sections, exact generic MW16 height Gram",
        "public_rank30_fixture_loaded_by_search": False,
    }:
        raise ArithmeticError("blindness boundary changed")

    legacy = SourceFileLoader("curve398_adaptive_legacy", str(LEGACY)).load_module()
    # The reusable curve-385 routines parameterize their dimensions by module
    # constants.  Override them before any classification or quotient lift.
    legacy.GENERIC_DIMENSION = GENERIC_DIMENSION
    legacy.OLD_CLASS_COUNT = EXPECTED_DEEPEST_COUNT
    model = tuple(Fraction(value) for value in blind["short_model"])
    generic = tuple(
        (Fraction(record["x"]), Fraction(record["y"]))
        for record in blind["generic_points"]
    )
    generic_gram = tuple(
        tuple(Fraction(value) for value in row)
        for row in blind["generic_height_gram"]
    )
    if len(generic) != GENERIC_DIMENSION or any(len(row) != GENERIC_DIMENSION for row in generic_gram):
        raise ArithmeticError("redacted generic input has the wrong dimension")

    signatures = legacy.find_mod2_reduction_certificate(
        model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND
    )
    if legacy.combined_mod2_rank(signatures, GENERIC_DIMENSION) != GENERIC_DIMENSION:
        raise ArithmeticError("generic MW16 specialization lost independence")

    # The exact generic Gram is half-integral.  CosetOracle sees 2G, so its
    # integer norm N represents actual half-lattice depth N/8.
    twice_gram = tuple(
        tuple(int(2 * value) for value in row) for row in generic_gram
    )
    if any(Fraction(twice_gram[i][j], 2) != generic_gram[i][j] for i in range(16) for j in range(16)):
        raise ArithmeticError("generic MW16 Gram is not half-integral")
    generic_oracle = legacy.CosetOracle(twice_gram)
    generic_rows = []
    histogram = Counter()
    for mask in range(1 << GENERIC_DIMENSION):
        residue = tuple((mask >> index) & 1 for index in range(GENERIC_DIMENSION))
        norm, representative, error = generic_oracle.solve(residue)
        if error > 1.0e-6:
            raise ArithmeticError("generic half-lattice CVP exact recomputation failed")
        histogram[norm] += 1
        generic_rows.append((norm, mask, representative))
    maximum_norm = max(row[0] for row in generic_rows)
    deepest = [row for row in generic_rows if row[0] == maximum_norm]
    if maximum_norm != EXPECTED_DEEPEST_TWICE_NORM or len(deepest) != EXPECTED_DEEPEST_COUNT:
        raise ArithmeticError("complete generic MW16 half-lattice census changed")
    old_masks = tuple(row[1] for row in deepest)

    specialized_gram, asymmetry = legacy.canonical_height_gram(model, generic)
    specialized_runs = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(specialized_gram, scale))
        ranked = []
        maximum_error = 0.0
        for unused_norm, mask, generic_representative in deepest:
            residue = tuple((mask >> index) & 1 for index in range(GENERIC_DIMENSION))
            unused_scaled_norm, representative, error = oracle.solve(residue)
            depth = legacy.quadratic_decimal(specialized_gram, representative) / 4
            ranked.append((depth, mask, representative, generic_representative))
            maximum_error = max(maximum_error, error)
        ranked.sort(key=lambda row: (-row[0], row[1]))
        specialized_runs[scale] = (ranked, maximum_error)
    initial_ranked = specialized_runs[OPERATIVE_SCALE][0]
    audit_ranked = specialized_runs[AUDIT_SCALE][0]
    audit_map = {row[1]: row[2] for row in audit_ranked}

    payload = {
        "schema": "elliptic-curves.curve398-mw16-adaptive-half-lattice-blind.v1",
        "status": "PARTIAL_CHECKPOINT",
        "blindness_boundary": {
            "sole_data_input": relative(INPUT),
            "public_rank30_fixture_loaded": False,
            "public_embedding_loaded": False,
            "held_out_point_coordinates_loaded": False,
            "held_out_point_count_known_only": 14,
            "forbidden_truth_artifact": "artifacts/generated-results/elliptic-curves/icarm_curve398_hidden_a1_mw16_v1.json",
        },
        "curve": {
            "label": blind["curve_label"],
            "short_model": [str(value) for value in model],
            "generic_rank": GENERIC_DIMENSION,
            "generic_points": [legacy.point_record(point) for point in generic],
            "generic_mod2_independence_rank": GENERIC_DIMENSION,
        },
        "generic_half_lattice": {
            "complete_class_count": 1 << GENERIC_DIMENSION,
            "twice_norm_histogram": {str(key): value for key, value in sorted(histogram.items())},
            "maximum_twice_norm": maximum_norm,
            "maximum_depth": str(Fraction(maximum_norm, 8)),
            "deepest_class_count": len(deepest),
            "deepest_masks": list(old_masks),
            "deepest_hex": [f"0x{mask:04x}" for mask in old_masks],
            "generic_representatives": [list(row[2]) for row in deepest],
        },
        "specialized_ranking": {
            "canonical_height_maximum_asymmetry": str(asymmetry),
            "operative_rounding_scale": OPERATIVE_SCALE,
            "audit_rounding_scale": AUDIT_SCALE,
            "maximum_cvp_distance_error": {
                str(scale): error for scale, (unused_rows, error) in specialized_runs.items()
            },
            "representative_disagreement_count": sum(
                audit_map[row[1]] != row[2] for row in initial_ranked
            ),
            "priority_order_identical_between_scales": [row[1] for row in initial_ranked] == [row[1] for row in audit_ranked],
        },
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "maximum_quotient_bits": args.max_quotient_bits,
            "maximum_nonzero_lifts_per_iteration": args.max_planned_lifts,
            "maximum_iterations": args.max_iterations,
            "checkpoint_every_completed_searches": args.checkpoint_every,
        },
        "initial_search": {"status": "SEARCHING", "cover_records": []},
        "iterations": [],
        "discoveries": [],
        "current_basis": [legacy.point_record(point) for point in generic],
        "searched_base_point_keys": [],
        "input_hashes": {
            relative(INPUT): digest(INPUT),
            relative(LEGACY): digest(LEGACY),
            relative(ENGINE): digest(ENGINE),
            relative(Path(__file__)): digest(Path(__file__)),
        },
        "software": {"python": platform.python_version()},
        "claim_boundary": [
            "No public complement coordinates or public-coordinate relations are read by this runner.",
            "Every accepted point and relation is replayed by exact rational group law.",
            "The generic half-lattice census is complete and its depths use the exact half-integral MW16 form.",
            "Specialized height ordering is high-precision numerical evidence checked at two rounded scales.",
            "Every point-search miss is bounded by the declared height and wall timeout.",
            "The blind run proves only its discovered subgroup; comparison with the held-out public subgroup is separate.",
        ],
    }
    write_payload(args.output, payload)

    discoveries = {}
    searched_keys = set()
    initial_started_wall = time.monotonic()
    initial_started_cpu = cpu_clock()
    for priority, (depth, mask, representative, generic_representative) in enumerate(initial_ranked, 1):
        base_point = legacy.exact_linear_combination(model[3], generic, representative)
        if base_point is None:
            raise ArithmeticError("a deepest half-class produced infinity")
        base_key = legacy.point_identifier(base_point)
        outcome = legacy.engine.run_quartic_search(
            mask=mask,
            representative=representative,
            short_model=model,
            generic_points=generic,
            height_bound=args.height_bound,
            timeout_seconds=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        searched_keys.add(base_key)
        source = f"initial:priority:{priority}:mask:{mask:#06x}"
        for point in outcome.curve_points:
            discoveries.setdefault(legacy.canonical_point(point), set()).add(source)
        payload["initial_search"]["cover_records"].append(
            {
                "priority": priority,
                "mask": mask,
                "hex": f"0x{mask:04x}",
                "generic_representative": list(generic_representative),
                "specialized_representative": list(representative),
                "specialized_depth": str(depth),
                "base_point_key": base_key,
                "search": outcome.record,
            }
        )
        payload["discoveries"] = legacy.discovery_records(discoveries)
        payload["searched_base_point_keys"] = sorted(searched_keys)
        write_payload(args.output, payload)
        print(
            f"C398MW16|initial={priority}/{len(initial_ranked)}|mask={mask:#06x}|"
            f"status={outcome.record['status']}|points={len(outcome.curve_points)}",
            flush=True,
        )

    basis, classification = legacy.classify_discovered_group(
        model=model,
        basis=generic,
        discoveries=discoveries,
        relation_chunk_size=args.relation_chunk_size,
        relation_timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    payload["initial_search"].update(
        {
            "status": "CLASSIFIED" if classification["status"] == "PASS_BASIS_EQUALS_DISCOVERED_GROUP" else classification["status"],
            "basis_rank_before": GENERIC_DIMENSION,
            "basis_rank_after": len(basis),
            "discovered_group_saturation": classification,
            "wall_seconds": time.monotonic() - initial_started_wall,
            "cpu_seconds": cpu_clock() - initial_started_cpu,
        }
    )
    if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        payload["status"] = "STOPPED_FAIL_CLOSED_UNCLASSIFIED_DISCOVERIES"
        write_payload(args.output, payload)
        return
    payload["current_basis"] = [legacy.point_record(point) for point in basis]
    write_payload(args.output, payload)
    print(f"C398MW16|transition=M16->M{len(basis)}|classification=PASS", flush=True)

    args.generic_points = generic
    for iteration_index in range(1, args.max_iterations + 1):
        rank_before = len(basis)
        quotient_bits = rank_before - GENERIC_DIMENSION
        if quotient_bits == 0:
            payload["status"] = "PASS_STABLE_DISCOVERED_SUBGROUP"
            payload["stable_rank"] = rank_before
            payload["stable_after_phase"] = "initial"
            write_payload(args.output, payload)
            return
        planned_count = EXPECTED_DEEPEST_COUNT * ((1 << quotient_bits) - 1)
        if quotient_bits > args.max_quotient_bits or planned_count > args.max_planned_lifts:
            payload["status"] = "STOPPED_AT_DECLARED_LIFT_LIMIT"
            payload["stop"] = {
                "basis_rank": rank_before,
                "quotient_bits": quotient_bits,
                "next_nonzero_lift_count": planned_count,
            }
            write_payload(args.output, payload)
            print(f"C398MW16|status=LIMIT|rank={rank_before}|bits={quotient_bits}|planned={planned_count}", flush=True)
            return

        ranked, ranking = legacy.rank_lifts(model, basis, old_masks, args)
        if len(ranked) != planned_count:
            raise ArithmeticError("adaptive lift census has the wrong size")
        iteration = {
            "iteration": iteration_index,
            "status": "SEARCHING",
            "basis_rank_before": rank_before,
            "quotient_bits": quotient_bits,
            "planned_nonzero_lifts": planned_count,
            "ranking": ranking,
            "cover_records": [],
            "unchanged_previously_searched_chart_count": 0,
        }
        payload["iterations"].append(iteration)
        write_payload(args.output, payload)
        print(f"C398MW16|iteration={iteration_index}|rank={rank_before}|bits={quotient_bits}|lifts={planned_count}|status=START", flush=True)
        started_wall = time.monotonic()
        started_cpu = cpu_clock()
        searched_this_iteration = 0
        for priority, (depth, old_mask, quotient_word, residue, representative) in enumerate(ranked, 1):
            base_point = legacy.exact_linear_combination(model[3], basis, representative)
            if base_point is None:
                raise ArithmeticError("adaptive half-class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                iteration["unchanged_previously_searched_chart_count"] += 1
                continue
            outcome = legacy.engine.run_quartic_search(
                mask=sum(int(bit) << index for index, bit in enumerate(residue)),
                representative=representative,
                short_model=model,
                generic_points=basis,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            searched_keys.add(base_key)
            source = f"iteration:{iteration_index}:priority:{priority}"
            for point in outcome.curve_points:
                discoveries.setdefault(legacy.canonical_point(point), set()).add(source)
            iteration["cover_records"].append(
                {
                    "priority": priority,
                    "old_mask": old_mask,
                    "old_hex": f"0x{old_mask:04x}",
                    "quotient_word": quotient_word,
                    "current_basis_residue": list(residue),
                    "canonical_depth": str(depth),
                    "representative": list(representative),
                    "base_point_key": base_key,
                    "search": outcome.record,
                }
            )
            searched_this_iteration += 1
            if searched_this_iteration % args.checkpoint_every == 0:
                payload["discoveries"] = legacy.discovery_records(discoveries)
                payload["searched_base_point_keys"] = sorted(searched_keys)
                write_payload(args.output, payload)
            print(
                f"C398MW16|iteration={iteration_index}|priority={priority}/{planned_count}|"
                f"old={old_mask:#06x}|qword={quotient_word}|status={outcome.record['status']}|"
                f"points={len(outcome.curve_points)}",
                flush=True,
            )

        old_basis_hash = legacy.canonical_hash([legacy.point_record(point) for point in basis])
        basis, classification = legacy.classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=args.relation_chunk_size,
            relation_timeout_seconds=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        new_basis_hash = legacy.canonical_hash([legacy.point_record(point) for point in basis])
        iteration.update(
            {
                "status": "CLASSIFIED" if classification["status"] == "PASS_BASIS_EQUALS_DISCOVERED_GROUP" else classification["status"],
                "searched_new_chart_count": len(iteration["cover_records"]),
                "discovered_group_saturation": classification,
                "basis_rank_after": len(basis),
                "group_changed": old_basis_hash != new_basis_hash,
                "wall_seconds": time.monotonic() - started_wall,
                "cpu_seconds": cpu_clock() - started_cpu,
            }
        )
        payload["current_basis"] = [legacy.point_record(point) for point in basis]
        payload["discoveries"] = legacy.discovery_records(discoveries)
        payload["searched_base_point_keys"] = sorted(searched_keys)
        if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            payload["status"] = "STOPPED_FAIL_CLOSED_UNCLASSIFIED_DISCOVERIES"
            write_payload(args.output, payload)
            return
        if old_basis_hash == new_basis_hash:
            payload["status"] = "PASS_STABLE_DISCOVERED_SUBGROUP"
            payload["stable_rank"] = len(basis)
            payload["stable_after_iteration"] = iteration_index
            write_payload(args.output, payload)
            print(f"C398MW16|iteration={iteration_index}|status=STABLE|rank={len(basis)}", flush=True)
            return
        write_payload(args.output, payload)
        print(f"C398MW16|iteration={iteration_index}|status=GROW|rank={rank_before}->{len(basis)}", flush=True)

    payload["status"] = "STOPPED_AT_DECLARED_ITERATION_LIMIT"
    payload["stop"] = {"basis_rank": len(basis)}
    write_payload(args.output, payload)


if __name__ == "__main__":
    main()
