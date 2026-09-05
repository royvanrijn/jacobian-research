#!/usr/bin/env sage-python
"""Equal-budget blind adaptive half-lattice recovery on refreshed R17 fibres.

The only mathematical data input is the MW17-only redacted fixture.  For each
fibre this runner searches the same 43 generic-deepest charts.  If at least one
quotient direction is recovered, it constructs a deterministic mod-2
complement in the discovered basis, binds a new chart order to that exact
lattice state, and searches a 301-chart adaptive pool.  Public points P18 and
beyond, displayed ranks, and displayed jump labels are never loaded.

Every reported score is rank(discovered group)-17, certified before the public
complement is opened.  Search misses retain only their bounded meaning.
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
import resource
import shutil
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
INPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v1.json"
LEGACY = CAS / "run_curve385_iterated_half_lattice_search.sage"
ENGINE_SOURCE = CAS / "half_lattice_fake_descent_replay.sage"
POLICY_SOURCE = CAS / "half_lattice_chart_policy.py"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_universal_pointed_v1.json"

GENERIC_DIMENSION = 17
INITIAL_CHARTS = 43
ADAPTIVE_CHARTS = 301
TOTAL_CHART_CAP = INITIAL_CHARTS + ADAPTIVE_CHARTS
OPERATIVE_SCALE = 1_000_000
AUDIT_SCALE = 100_000

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]


from pointed_quartic_migration import validate_frozen_sources, runtime_search, require_runtime

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


def identity_rows(dimension: int):
    return [
        [int(row == column) for column in range(dimension)]
        for row in range(dimension)
    ]


def add_mod2(left, right):
    return tuple((int(a) + int(b)) & 1 for a, b in zip(left, right))


def base_residue(mask, generic_coordinates, dimension):
    residue = (0,) * dimension
    for index, vector in enumerate(generic_coordinates):
        if (mask >> index) & 1:
            residue = add_mod2(residue, vector)
    return residue


def quotient_words(bit_count: int):
    return sorted(
        range(1, 1 << bit_count),
        key=lambda word: (word.bit_count(), word),
    )


def lift_residue(mask, quotient_word, generic_coordinates, complement, dimension):
    residue = base_residue(mask, generic_coordinates, dimension)
    for index, vector in enumerate(complement):
        if (quotient_word >> index) & 1:
            residue = add_mod2(residue, vector)
    return residue


def compact_search_record(record):
    return {
        "status": record["status"],
        "pointed_search": record if record.get("backend") == "pointed_quartic_search_v1" else None,
        "backend": record.get("backend"),
        "wall_seconds": record.get("wall_seconds"),
        "search_milliseconds": record.get("search_milliseconds"),
        "integral_model_maximum_coefficient_bits": record.get(
            "integral_model_maximum_coefficient_bits"
        ),
        "reduced_model_maximum_coefficient_bits": (
            record.get("reduced_model", {}).get("maximum_coefficient_bits")
        ),
        "finite_curve_point_count": len(record.get("finite_curve_points", [])),
        "curve_points": record.get("curve_points", []),
        "pari_error": record.get("pari_error"),
    }


def complete_generic_census(legacy, gram):
    oracle = legacy.CosetOracle(gram)
    rows = []
    maximum_error = 0.0
    for mask in range(1 << GENERIC_DIMENSION):
        residue = tuple((mask >> index) & 1 for index in range(GENERIC_DIMENSION))
        norm, representative, error = oracle.solve(residue)
        maximum_error = max(maximum_error, error)
        rows.append((norm, mask, representative))
    rows.sort(key=lambda row: (-row[0], row[1]))
    if len(rows) != 1 << GENERIC_DIMENSION:
        raise ArithmeticError("the complete generic half-lattice census is incomplete")
    if rows[0][0] != 12 or sum(row[0] == 12 for row in rows) != INITIAL_CHARTS:
        raise ArithmeticError("the determinant-948 deepest-hole census changed")
    return rows, maximum_error


def rank_initial_charts(legacy, model, generic, generic_rows):
    gram, asymmetry = legacy.canonical_height_gram(model, generic)
    selected = generic_rows[:INITIAL_CHARTS]
    runs = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(gram, scale))
        rows = []
        maximum_error = 0.0
        for generic_norm, mask, generic_representative in selected:
            residue = tuple((mask >> index) & 1 for index in range(GENERIC_DIMENSION))
            unused_norm, representative, error = oracle.solve(residue)
            maximum_error = max(maximum_error, error)
            depth = legacy.quadratic_decimal(gram, representative) / 4
            rows.append((depth, mask, representative, generic_representative, generic_norm))
        rows.sort(key=lambda row: (-row[0], row[1]))
        runs[scale] = rows, maximum_error
    operative = runs[OPERATIVE_SCALE][0]
    audit = runs[AUDIT_SCALE][0]
    audit_map = {row[1]: row[2] for row in audit}
    return gram, operative, {
        "canonical_height_maximum_asymmetry": str(asymmetry),
        "operative_rounding_scale": OPERATIVE_SCALE,
        "audit_rounding_scale": AUDIT_SCALE,
        "maximum_cvp_distance_error": {
            str(scale): error for scale, (unused_rows, error) in runs.items()
        },
        "representative_disagreement_count": sum(
            audit_map[row[1]] != row[2] for row in operative
        ),
        "priority_order_identical_between_scales": [row[1] for row in audit]
        == [row[1] for row in operative],
    }


def rank_adaptive_pool(legacy, model, basis, generic, generic_rows, args):
    class LiftArgs:
        pass

    lift_args = LiftArgs()
    lift_args.generic_points = generic
    lift_args.relation_chunk_size = args.relation_chunk_size
    lift_args.relation_timeout_seconds = args.relation_timeout_seconds
    lift_args.stack_bytes = args.stack_bytes
    generic_coordinates, complement = legacy.gf2_lift_data(
        model, basis, generic, lift_args
    )
    quotient_bit_count = len(complement)
    if quotient_bit_count <= 0:
        raise ArithmeticError("adaptive ranking requires a nonzero quotient")
    words = quotient_words(quotient_bit_count)
    pool = []
    for index, (generic_norm, mask, unused_generic_representative) in enumerate(
        generic_rows[:ADAPTIVE_CHARTS]
    ):
        quotient_word = words[index % len(words)]
        residue = lift_residue(
            mask,
            quotient_word,
            generic_coordinates,
            complement,
            len(basis),
        )
        pool.append((generic_norm, mask, quotient_word, residue))
    if len({row[3] for row in pool}) != ADAPTIVE_CHARTS:
        raise ArithmeticError("the deterministic adaptive pool contains collisions")

    gram, asymmetry = legacy.canonical_height_gram(model, basis)
    runs = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(gram, scale))
        rows = []
        maximum_error = 0.0
        for generic_norm, mask, quotient_word, residue in pool:
            unused_norm, representative, error = oracle.solve(residue)
            maximum_error = max(maximum_error, error)
            depth = legacy.quadratic_decimal(gram, representative) / 4
            rows.append(
                (depth, mask, quotient_word, residue, representative, generic_norm)
            )
        rows.sort(key=lambda row: (-row[0], row[1], row[2]))
        runs[scale] = rows, maximum_error
    operative = runs[OPERATIVE_SCALE][0]
    audit = runs[AUDIT_SCALE][0]
    audit_map = {(row[1], row[2]): row[4] for row in audit}
    return gram, operative, generic_coordinates, complement, {
        "canonical_height_maximum_asymmetry": str(asymmetry),
        "operative_rounding_scale": OPERATIVE_SCALE,
        "audit_rounding_scale": AUDIT_SCALE,
        "maximum_cvp_distance_error": {
            str(scale): error for scale, (unused_rows, error) in runs.items()
        },
        "representative_disagreement_count": sum(
            audit_map[(row[1], row[2])] != row[4] for row in operative
        ),
        "priority_order_identical_between_scales": [
            (row[1], row[2]) for row in audit
        ]
        == [(row[1], row[2]) for row in operative],
        "quotient_bit_count": quotient_bit_count,
        "generic_coordinate_rows_in_current_basis": [
            list(row) for row in generic_coordinates
        ],
        "quotient_complement_rows_mod2": [list(row) for row in complement],
        "pool_construction": (
            "first 301 exact-generic half-classes in decreasing generic depth; "
            "paired cyclically with nonzero quotient words ordered by Hamming "
            "weight then integer value"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--protocol", type=Path, default=PROTOCOL)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()
    if args.height_bound != 100_000 or args.timeout_seconds != 15.0:
        raise SystemExit("the frozen ladder requires height 100000 and timeout 15s")
    if args.stack_bytes != 1_000_000_000:
        raise SystemExit("the frozen ladder requires a 1GB PARI stack per search")
    if args.relation_chunk_size != 64 or args.relation_timeout_seconds != 180.0:
        raise SystemExit("the frozen relation-classification budget changed")
    if args.checkpoint_every < 1 or shutil.which("gp") is None:
        raise SystemExit("invalid checkpoint setting or missing PARI/GP")

    blind_input = json.loads(args.input.read_text())
    protocol = json.loads(args.protocol.read_text())
    if blind_input.get("status") != "FROZEN_MW17_ONLY_NO_PUBLIC_COMPLEMENT":
        raise ArithmeticError("the MW17-only ladder input is not frozen")
    if blind_input["redaction"]["contains_displayed_complement_coordinates"]:
        raise ArithmeticError("the blind input contains public complement coordinates")
    if blind_input["redaction"]["contains_displayed_rank_or_jump"]:
        raise ArithmeticError("the blind input contains displayed truth labels")
    if protocol.get("status") != "FROZEN_BEFORE_BLIND_RECOVERY":
        raise ArithmeticError("the jump-ladder protocol is not frozen")
    if protocol["blind_input_sha256"] != digest(args.input):
        raise ArithmeticError("the protocol names another blind input")
    validate_frozen_sources(protocol["implementation_hashes"])
    if protocol["search_policy"]["total_chart_cap_per_fibre"] != TOTAL_CHART_CAP:
        raise ArithmeticError("the protocol chart budget changed")

    legacy = SourceFileLoader("r17_jump_ladder_legacy", str(LEGACY)).load_module()
    engine = SourceFileLoader("r17_jump_ladder_engine", str(ENGINE_SOURCE)).load_module()
    from pointed_quartic_search import run_quartic_search as shared_quartic_search
    engine.run_quartic_search = shared_quartic_search
    chart_policy = SourceFileLoader("r17_jump_ladder_policy", str(POLICY_SOURCE)).load_module()
    legacy.GENERIC_DIMENSION = GENERIC_DIMENSION
    legacy.OLD_CLASS_COUNT = INITIAL_CHARTS

    payload = {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-blind.v1",
        "status": "PARTIAL_CHECKPOINT",
        "blindness_boundary": {
            "sole_mathematical_data_input": relative(args.input),
            "public_complement_loaded": False,
            "displayed_rank_or_jump_loaded": False,
            "truth_artifact_loaded": False,
            "response_field": "exact_quotient_rank_recovered_before_public_complement",
        },
        "protocol_sha256": digest(args.protocol),
        "runtime_search": runtime_search(),
        "declared_budget": protocol["search_policy"],
        "results": [],
        "input_hashes": {
            relative(args.input): digest(args.input),
            relative(args.protocol): digest(args.protocol),
            relative(LEGACY): digest(LEGACY),
            relative(ENGINE_SOURCE): digest(ENGINE_SOURCE),
            relative(POLICY_SOURCE): digest(POLICY_SOURCE),
            relative(Path(__file__).resolve()): digest(Path(__file__).resolve()),
        },
        "generation": {
            "python": platform.python_version(),
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
            ),
        },
        "claim_boundary": [
            "Only exactly certified independent directions change the recovered quotient rank.",
            "Discovered-group saturation is relative to all returned blind points, not the unknown full E(Q).",
            "Every miss is bounded by the fixed height and timeout budgets.",
            "No public point P18 or beyond, displayed rank, or displayed jump is available to this runner.",
        ],
    }
    write_payload(args.output, payload)

    for case_index, case in enumerate(blind_input["cases"], 1):
        curve_id = int(case["curve_id"])
        model = tuple(Fraction(value) for value in case["short_model"])
        generic = tuple(
            (Fraction(row["x"]), Fraction(row["y"]))
            for row in case["generic_points"]
        )
        generic_gram = tuple(
            tuple(int(value) for value in row)
            for row in case["generic_height_gram"]
        )
        if len(generic) != GENERIC_DIMENSION:
            raise ArithmeticError(f"curve {curve_id}: generic dimension changed")
        signatures = legacy.find_mod2_reduction_certificate(
            model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND
        )
        if legacy.combined_mod2_rank(signatures, GENERIC_DIMENSION) != GENERIC_DIMENSION:
            raise ArithmeticError(f"curve {curve_id}: generic MW17 lost independence")

        started_wall = time.monotonic()
        started_cpu = cpu_clock()
        generic_rows, generic_cvp_error = complete_generic_census(legacy, generic_gram)
        initial_gram, initial_ranked, initial_ranking = rank_initial_charts(
            legacy, model, generic, generic_rows
        )
        initial_chart_ids = [f"mask:{row[1]:05x}" for row in initial_ranked]
        initial_certificate = chart_policy.bind_ordering(
            basis_records=[legacy.point_record(point) for point in generic],
            height_gram_rows=initial_gram,
            generic_coordinate_rows=identity_rows(GENERIC_DIMENSION),
            quotient_coordinate_rows=[],
            chart_universe_id=(
                f"r17-refresh-jump-ladder-v1:curve-{curve_id}:generic-deepest43"
            ),
            ordered_chart_ids=initial_chart_ids,
            heuristics=["legacy_half_lattice_depth"],
        )
        chart_policy.validate_ordering(
            initial_certificate,
            basis_records=[legacy.point_record(point) for point in generic],
            height_gram_rows=initial_gram,
            generic_coordinate_rows=identity_rows(GENERIC_DIMENSION),
            quotient_coordinate_rows=[],
            chart_universe_id=(
                f"r17-refresh-jump-ladder-v1:curve-{curve_id}:generic-deepest43"
            ),
            ordered_chart_ids=initial_chart_ids,
        )
        result = {
            "curve_id": curve_id,
            "representative_class": case["representative_class"],
            "native_chart": case["native_chart"],
            "status": "SEARCHING_INITIAL",
            "generic_census": {
                "complete_class_count": 1 << GENERIC_DIMENSION,
                "maximum_exact_generic_norm": generic_rows[0][0],
                "deepest_class_count": sum(
                    row[0] == generic_rows[0][0] for row in generic_rows
                ),
                "maximum_cvp_distance_error": generic_cvp_error,
                "top301_masks_sha256": legacy.canonical_hash(
                    [row[1] for row in generic_rows[:ADAPTIVE_CHARTS]]
                ),
            },
            "initial": {
                "status": "SEARCHING",
                "ranking": initial_ranking,
                "ordering_certificate": initial_certificate,
                "cover_records": [],
            },
            "adaptive": None,
        }
        payload["results"].append(result)
        write_payload(args.output, payload)
        print(
            f"R17JUMPLADDER|case={case_index}/{len(blind_input['cases'])}|"
            f"curve={curve_id}|stage=initial|charts={INITIAL_CHARTS}",
            flush=True,
        )

        discoveries = {}
        searched_keys = set()
        for priority, (depth, mask, representative, generic_representative, generic_norm) in enumerate(
            initial_ranked, 1
        ):
            base_point = legacy.exact_linear_combination(model[3], generic, representative)
            if base_point is None:
                raise ArithmeticError("an initial half-class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                raise ArithmeticError("an initial pointed chart was duplicated")
            outcome = engine.run_quartic_search(
                mask=mask,
                representative=representative,
                short_model=model,
                generic_points=generic,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            searched_keys.add(base_key)
            source = f"initial:priority:{priority}:mask:{mask:#07x}"
            for point in outcome.curve_points:
                point = legacy.canonical_point(point)
                discoveries.setdefault(point, set()).add(source)
            result["initial"]["cover_records"].append(
                {
                    "priority": priority,
                    "mask": mask,
                    "hex": f"0x{mask:05x}",
                    "exact_generic_norm": generic_norm,
                    "generic_representative": list(generic_representative),
                    "specialized_representative": list(representative),
                    "specialized_depth": str(depth),
                    "base_point_key": base_key,
                    "search": compact_search_record(outcome.record),
                }
            )
            if priority % args.checkpoint_every == 0:
                write_payload(args.output, payload)

        basis, initial_classification = legacy.classify_discovered_group(
            model=model,
            basis=generic,
            discoveries=discoveries,
            relation_chunk_size=args.relation_chunk_size,
            relation_timeout_seconds=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        if initial_classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            result["status"] = "UNKNOWN_UNCLASSIFIED_INITIAL_DISCOVERIES"
            result["initial"]["status"] = result["status"]
            payload["status"] = "STOPPED_FAIL_CLOSED"
            write_payload(args.output, payload)
            return
        initial_gain = len(basis) - GENERIC_DIMENSION
        result["initial"].update(
            {
                "status": "CLASSIFIED",
                "discovered_group_classification": initial_classification,
                "exact_quotient_rank_recovered": initial_gain,
            }
        )
        print(
            f"R17JUMPLADDER|curve={curve_id}|stage=initial-classified|gain={initial_gain}",
            flush=True,
        )

        if initial_gain == 0:
            result.update(
                {
                    "status": "PASS_COMPLETE_STRUCTURAL_STOP_NO_ADAPTIVE_QUOTIENT",
                    "adaptive": {
                        "status": "NOT_APPLICABLE_ZERO_DISCOVERED_QUOTIENT",
                        "searched_chart_count": 0,
                    },
                    "attempted_chart_count": INITIAL_CHARTS,
                    "exact_quotient_rank_recovered_before_public_complement": 0,
                    "final_basis": [legacy.point_record(point) for point in basis],
                    "wall_seconds": time.monotonic() - started_wall,
                    "cpu_seconds": cpu_clock() - started_cpu,
                }
            )
            write_payload(args.output, payload)
            continue

        adaptive_gram, adaptive_ranked, generic_coordinates, complement, adaptive_ranking = rank_adaptive_pool(
            legacy, model, basis, generic, generic_rows, args
        )
        adaptive_chart_ids = [
            f"gmask:{row[1]:05x}:qword:{row[2]:x}" for row in adaptive_ranked
        ]
        adaptive_certificate = chart_policy.bind_ordering(
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=(
                f"r17-refresh-jump-ladder-v1:curve-{curve_id}:adaptive301"
            ),
            ordered_chart_ids=adaptive_chart_ids,
            heuristics=[
                "legacy_half_lattice_depth",
                "quotient_hamming_weight",
            ],
        )
        chart_policy.validate_ordering(
            adaptive_certificate,
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=(
                f"r17-refresh-jump-ladder-v1:curve-{curve_id}:adaptive301"
            ),
            ordered_chart_ids=adaptive_chart_ids,
        )
        result["status"] = "SEARCHING_ADAPTIVE"
        result["adaptive"] = {
            "status": "SEARCHING",
            "basis_rank_before": len(basis),
            "basis_before": [legacy.point_record(point) for point in basis],
            "ranking": adaptive_ranking,
            "ordering_certificate": adaptive_certificate,
            "cover_records": [],
        }
        write_payload(args.output, payload)
        print(
            f"R17JUMPLADDER|curve={curve_id}|stage=adaptive|"
            f"bits={initial_gain}|charts={ADAPTIVE_CHARTS}",
            flush=True,
        )

        for priority, (depth, generic_mask, quotient_word, residue, representative, generic_norm) in enumerate(
            adaptive_ranked, 1
        ):
            base_point = legacy.exact_linear_combination(model[3], basis, representative)
            if base_point is None:
                raise ArithmeticError("an adaptive half-class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                raise ArithmeticError("an adaptive pointed chart was already searched")
            mask = sum(int(bit) << index for index, bit in enumerate(residue))
            outcome = engine.run_quartic_search(
                mask=mask,
                representative=representative,
                short_model=model,
                generic_points=basis,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            searched_keys.add(base_key)
            source = f"adaptive:priority:{priority}:gmask:{generic_mask:#07x}:qword:{quotient_word:#x}"
            for point in outcome.curve_points:
                point = legacy.canonical_point(point)
                discoveries.setdefault(point, set()).add(source)
            result["adaptive"]["cover_records"].append(
                {
                    "priority": priority,
                    "generic_mask": generic_mask,
                    "generic_hex": f"0x{generic_mask:05x}",
                    "exact_generic_norm": generic_norm,
                    "quotient_word": quotient_word,
                    "quotient_word_binary": f"{quotient_word:0{initial_gain}b}",
                    "current_basis_residue": list(residue),
                    "representative": list(representative),
                    "canonical_depth": str(depth),
                    "base_point_key": base_key,
                    "search": compact_search_record(outcome.record),
                }
            )
            if priority % args.checkpoint_every == 0:
                write_payload(args.output, payload)

        basis, final_classification = legacy.classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=args.relation_chunk_size,
            relation_timeout_seconds=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        if final_classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            result["status"] = "UNKNOWN_UNCLASSIFIED_ADAPTIVE_DISCOVERIES"
            result["adaptive"]["status"] = result["status"]
            payload["status"] = "STOPPED_FAIL_CLOSED"
            write_payload(args.output, payload)
            return
        final_gain = len(basis) - GENERIC_DIMENSION
        result["adaptive"].update(
            {
                "status": "CLASSIFIED",
                "discovered_group_classification": final_classification,
                "exact_incremental_quotient_rank_recovered": final_gain - initial_gain,
                "basis_rank_after": len(basis),
            }
        )
        result.update(
            {
                "status": "PASS_COMPLETE_EQUAL_BUDGET_BLIND_RECOVERY",
                "attempted_chart_count": TOTAL_CHART_CAP,
                "bounded_complete_chart_count": sum(
                    row["search"]["status"] == "bounded_search_complete"
                    for phase in (result["initial"], result["adaptive"])
                    for row in phase["cover_records"]
                ),
                "timeout_chart_count": sum(
                    row["search"]["status"] == "bounded_search_timeout"
                    for phase in (result["initial"], result["adaptive"])
                    for row in phase["cover_records"]
                ),
                "pari_failure_chart_count": sum(
                    row["search"]["status"] == "pari_failure"
                    for phase in (result["initial"], result["adaptive"])
                    for row in phase["cover_records"]
                ),
                "exact_quotient_rank_recovered_before_public_complement": final_gain,
                "final_basis": [legacy.point_record(point) for point in basis],
                "wall_seconds": time.monotonic() - started_wall,
                "cpu_seconds": cpu_clock() - started_cpu,
            }
        )
        write_payload(args.output, payload)
        print(
            f"R17JUMPLADDER|curve={curve_id}|status=PASS|"
            f"response={final_gain}|charts={TOTAL_CHART_CAP}",
            flush=True,
        )

    if len(payload["results"]) != blind_input["case_count"]:
        raise ArithmeticError("the blind ladder did not cover every eligible case")
    payload["status"] = "PASS_COMPLETE_BLIND_RECOVERY_BEFORE_PUBLIC_COMPLEMENT"
    payload["response"] = [
        {
            "curve_id": row["curve_id"],
            "exact_quotient_rank_recovered_before_public_complement": row[
                "exact_quotient_rank_recovered_before_public_complement"
            ],
        }
        for row in payload["results"]
    ]
    write_payload(args.output, payload)
    print(
        "R17JUMPLADDER|status=PASS|response="
        + ",".join(
            f"{row['curve_id']}:{row['exact_quotient_rank_recovered_before_public_complement']}"
            for row in payload["response"]
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
