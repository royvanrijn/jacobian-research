#!/usr/bin/env sage-python
"""Run one complement-blind adaptive MW16 calibration on curve 400.

The frozen input contains the target curve, sixteen specialized generic
sections, and the exact generic MW16 height Gram, but no exceptional points or
target rank.  The initial wave is the complete maximum-depth M/2M stratum.
If that exact wave recovers a five-dimensional quotient, this runner searches
all 4*(2^5-1)=124 nonzero quotient lifts and classifies the resulting group by
exact rational arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import platform
import resource
import shutil
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
LADDER = ROOT / "elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_curve400_adaptive_calibration_v1.json"
PARENT_ID = "curve400-p53042"
GENERIC_DIMENSION = 16
EXPECTED_INITIAL_QUOTIENT_BITS = 5
EXPECTED_DEEPEST_COUNT = 4
EXPECTED_ADAPTIVE_COUNT = EXPECTED_DEEPEST_COUNT * (
    (1 << EXPECTED_INITIAL_QUOTIENT_BITS) - 1
)


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


def read_discoveries(records):
    discoveries = {}
    for record in records:
        point_record = record["point"]
        point = (
            Fraction(point_record["x"]),
            Fraction(point_record["y"]),
        )
        discoveries[point] = set(record["sources"])
    return discoveries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--parent-id", default=PARENT_ID)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()
    if args.parent_id != PARENT_ID:
        raise SystemExit(f"this frozen calibration requires --parent-id {PARENT_ID}")
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if args.checkpoint_every <= 0 or shutil.which("gp") is None:
        raise SystemExit("invalid checkpoint bound or missing PARI/GP")

    blind = json.loads(args.input.read_text())
    if blind.get("status") != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("complement-blind MW16 input is not passing")
    parent = next(
        (row for row in blind["parents"] if row["parent_id"] == args.parent_id),
        None,
    )
    if parent is None or int(parent["curve_id"]) != 400:
        raise ArithmeticError("frozen curve-400 parent is absent")
    if any(
        blind["blindness_boundary"][key]
        for key in (
            "public_point_lists_loaded",
            "public_complement_coordinates_loaded",
            "target_rank_lower_bounds_loaded",
        )
    ):
        raise ArithmeticError("input crossed the complement-blind boundary")

    ladder = SourceFileLoader("mw16_c400_ladder", str(LADDER)).load_module()
    legacy = SourceFileLoader("mw16_c400_legacy", str(LEGACY)).load_module()
    legacy.GENERIC_DIMENSION = GENERIC_DIMENSION
    model = tuple(Fraction(value) for value in parent["target_short_model"])
    generic = tuple(
        (Fraction(row["x"]), Fraction(row["y"]))
        for row in parent["specialized_generic_points"]
    )

    initial = ladder.run_parent(parent, legacy, args)
    if initial["status"] != "PASS_COMPLETE_INITIAL_HALF_LATTICE_WAVE":
        raise ArithmeticError("initial exact half-lattice wave did not pass")
    if initial["generic_half_lattice"]["deepest_class_count"] != EXPECTED_DEEPEST_COUNT:
        raise ArithmeticError("curve-400 deepest-stratum size changed")
    if initial["exact_quotient_rank_recovered"] != EXPECTED_INITIAL_QUOTIENT_BITS:
        raise ArithmeticError("curve-400 initial calibration no longer recovers five directions")

    discoveries = read_discoveries(initial["discoveries"])
    basis, replay = legacy.classify_discovered_group(
        model=model,
        basis=generic,
        discoveries=discoveries,
        relation_chunk_size=args.relation_chunk_size,
        relation_timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if replay["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP" or len(basis) != 21:
        raise ArithmeticError("initial discovered M21 group did not replay")

    old_masks = tuple(initial["generic_half_lattice"]["deepest_masks"])
    legacy.OLD_CLASS_COUNT = len(old_masks)
    args.generic_points = generic
    ranked, ranking = legacy.rank_lifts(model, basis, old_masks, args)
    if len(ranked) != EXPECTED_ADAPTIVE_COUNT:
        raise ArithmeticError("adaptive curve-400 lift count changed")

    payload = {
        "schema": "elliptic-curves.icarm-mw16-curve400-adaptive-calibration.v1",
        "status": "SEARCHING",
        "blindness_boundary": {
            "sole_mathematical_input": relative(args.input),
            "parent_id": args.parent_id,
            "public_point_lists_loaded": False,
            "public_complement_coordinates_loaded": False,
            "target_rank_lower_bound_loaded": False,
        },
        "initial": initial,
        "adaptive": {
            "basis_rank_before": len(basis),
            "quotient_bits": len(basis) - GENERIC_DIMENSION,
            "complete_nonzero_lift_count": len(ranked),
            "ranking": ranking,
            "cover_records": [],
            "status": "SEARCHING",
        },
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "complete_adaptive_nonzero_lift_count": EXPECTED_ADAPTIVE_COUNT,
        },
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, LADDER, LEGACY, ENGINE, Path(__file__))
        },
        "software": {"python": platform.python_version()},
        "claim_boundary": [
            "The calibration is blind to the public exceptional complement and target rank.",
            "All four initial and all 124 adaptive quotient charts receive the same declared bounded search.",
            "Only exact rational points, exact group law, and finite-reduction independence certificates contribute rank.",
            "A timeout or completed miss has no point-absence, rank-upper-bound, saturation, covering, or Selmer meaning.",
            "Comparison with the held-out atlas jump occurs only in a separate verification layer.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/run_icarm_mw16_curve400_adaptive_calibration.sage"
        ),
    }
    write_payload(args.output, payload)

    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    adaptive = payload["adaptive"]
    for priority, (depth, old_mask, quotient_word, residue, representative) in enumerate(ranked, 1):
        base_point = legacy.exact_linear_combination(model[3], basis, representative)
        if base_point is None:
            raise ArithmeticError("adaptive half-class produced infinity")
        outcome = legacy.engine.run_quartic_search(
            mask=sum(int(bit) << index for index, bit in enumerate(residue)),
            representative=representative,
            short_model=model,
            generic_points=basis,
            height_bound=args.height_bound,
            timeout_seconds=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        source = f"adaptive:priority:{priority}:old:{old_mask:#06x}:qword:{quotient_word}"
        for point in outcome.curve_points:
            discoveries.setdefault(legacy.canonical_point(point), set()).add(source)
        adaptive["cover_records"].append(
            {
                "priority": priority,
                "old_mask": old_mask,
                "old_hex": f"0x{old_mask:04x}",
                "quotient_word": quotient_word,
                "current_basis_residue": list(residue),
                "canonical_depth": str(depth),
                "representative": list(representative),
                "base_point_key": legacy.point_identifier(base_point),
                "search": outcome.record,
            }
        )
        if priority % args.checkpoint_every == 0:
            payload["discoveries"] = legacy.discovery_records(discoveries)
            write_payload(args.output, payload)
        print(
            f"MW16C400ADAPT|chart={priority}/{len(ranked)}|old={old_mask:#06x}|"
            f"qword={quotient_word}|status={outcome.record['status']}|"
            f"points={len(outcome.curve_points)}",
            flush=True,
        )

    final_basis, classification = legacy.classify_discovered_group(
        model=model,
        basis=basis,
        discoveries=discoveries,
        relation_chunk_size=args.relation_chunk_size,
        relation_timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        payload["status"] = "STOPPED_FAIL_CLOSED_UNCLASSIFIED_DISCOVERIES"
        write_payload(args.output, payload)
        return
    adaptive.update(
        {
            "status": "CLASSIFIED",
            "basis_rank_after": len(final_basis),
            "incremental_exact_quotient_rank_recovered": len(final_basis) - len(basis),
            "discovered_group_saturation": classification,
            "wall_seconds": time.monotonic() - started_wall,
            "cpu_seconds": cpu_clock() - started_cpu,
        }
    )
    payload["discoveries"] = legacy.discovery_records(discoveries)
    payload["current_basis"] = [legacy.point_record(point) for point in final_basis]
    payload["exact_quotient_rank_recovered_total"] = len(final_basis) - GENERIC_DIMENSION
    payload["status"] = "PASS_COMPLETE_CURVE400_ADAPTIVE_CALIBRATION"
    write_payload(args.output, payload)
    print(
        f"MW16C400ADAPT|rank={GENERIC_DIMENSION}->{len(final_basis)}|"
        f"output={relative(args.output)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
