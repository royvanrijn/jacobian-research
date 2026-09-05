#!/usr/bin/env sage-python
"""Run the complement-blind initial half-lattice wave on nine MW16 parents.

Each parent is searched only through the complete deepest stratum of its exact
generic MW16 half-lattice.  Returned points are checked by exact
group law and the discovered group is classified by the same relation engine
used for the curve-398 and curve-385 adaptive replays.  Public point lists and
the five target jump labels are never loaded.

This is a bounded detector calibration.  A miss has no rank, saturation,
covering, or Selmer meaning.
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
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/local/elliptic-curves/pointed-quartic-search/campaigns/run_icarm_mw16_parent_ladder_blind.json"
GENERIC_DIMENSION = 16
EXPECTED_CLASS_COUNT = 1 << GENERIC_DIMENSION
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


def run_parent(parent, legacy, args):
    model = tuple(Fraction(value) for value in parent["target_short_model"])
    generic = tuple(
        (Fraction(record["x"]), Fraction(record["y"]))
        for record in parent["specialized_generic_points"]
    )
    generic_gram = tuple(
        tuple(Fraction(value) for value in row)
        for row in parent["generic_height_gram"]
    )
    if len(generic) != GENERIC_DIMENSION or any(
        len(row) != GENERIC_DIMENSION for row in generic_gram
    ):
        raise ArithmeticError("a parent input has the wrong MW dimension")
    signatures = legacy.find_mod2_reduction_certificate(
        model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND
    )
    if legacy.combined_mod2_rank(signatures, GENERIC_DIMENSION) != GENERIC_DIMENSION:
        raise ArithmeticError("specialized generic MW16 lost independence")

    twice_gram = tuple(
        tuple(int(2 * value) for value in row) for row in generic_gram
    )
    if any(
        Fraction(twice_gram[i][j], 2) != generic_gram[i][j]
        for i in range(GENERIC_DIMENSION)
        for j in range(GENERIC_DIMENSION)
    ):
        raise ArithmeticError("generic MW16 Gram is not half-integral")
    generic_oracle = legacy.CosetOracle(twice_gram)
    generic_rows = []
    histogram = Counter()
    for mask in range(EXPECTED_CLASS_COUNT):
        residue = tuple(
            (mask >> index) & 1 for index in range(GENERIC_DIMENSION)
        )
        norm, representative, error = generic_oracle.solve(residue)
        if error > 1.0e-6:
            raise ArithmeticError("generic half-lattice CVP recomputation failed")
        histogram[norm] += 1
        generic_rows.append((norm, mask, representative))
    maximum_norm = max(row[0] for row in generic_rows)
    deepest = [row for row in generic_rows if row[0] == maximum_norm]
    if maximum_norm != EXPECTED_DEEPEST_TWICE_NORM or not deepest:
        raise ArithmeticError("complete generic MW16 half-lattice census changed")

    if getattr(args, "exact_generic_order_only", False):
        # When the entire maximum-depth stratum is searched, a costly
        # specialized canonical-height order cannot change coverage.  Use the
        # exact generic representatives and mask order for prospective fibres;
        # minimization/canonical heights are deferred to positive candidates.
        initial_ranked = sorted(
            (
                Fraction(unused_norm, 8),
                mask,
                representative,
                representative,
            )
            for unused_norm, mask, representative in deepest
        )
        specialized_ranking = {
            "mode": "exact_generic_maximum_depth_then_mask",
            "specialized_canonical_height_computed": False,
            "coverage_depends_on_within_stratum_order": False,
        }
    else:
        specialized_gram, asymmetry = legacy.canonical_height_gram(model, generic)
        specialized_runs = {}
        for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
            oracle = legacy.CosetOracle(legacy.rounded_gram(specialized_gram, scale))
            ranked = []
            maximum_error = 0.0
            for unused_norm, mask, generic_representative in deepest:
                residue = tuple(
                    (mask >> index) & 1 for index in range(GENERIC_DIMENSION)
                )
                unused_scaled_norm, representative, error = oracle.solve(residue)
                depth = legacy.quadratic_decimal(specialized_gram, representative) / 4
                ranked.append(
                    (depth, mask, representative, generic_representative)
                )
                maximum_error = max(maximum_error, error)
            ranked.sort(key=lambda row: (-row[0], row[1]))
            specialized_runs[scale] = (ranked, maximum_error)
        initial_ranked = specialized_runs[OPERATIVE_SCALE][0]
        audit_ranked = specialized_runs[AUDIT_SCALE][0]
        audit_map = {row[1]: row[2] for row in audit_ranked}
        specialized_ranking = {
            "mode": "specialized_canonical_height",
            "specialized_canonical_height_computed": True,
            "canonical_height_maximum_asymmetry": str(asymmetry),
            "operative_rounding_scale": OPERATIVE_SCALE,
            "audit_rounding_scale": AUDIT_SCALE,
            "representative_disagreement_count": sum(
                audit_map[row[1]] != row[2] for row in initial_ranked
            ),
            "priority_order_identical_between_scales": [
                row[1] for row in initial_ranked
            ]
            == [row[1] for row in audit_ranked],
            "maximum_cvp_distance_error": {
                str(scale): error
                for scale, (unused_rows, error) in specialized_runs.items()
            },
        }

    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    discoveries = {}
    covers = []
    for priority, (
        depth,
        mask,
        representative,
        generic_representative,
    ) in enumerate(initial_ranked, 1):
        outcome = legacy.engine.run_quartic_search(
            mask=mask,
            representative=representative,
            short_model=model,
            generic_points=generic,
            height_bound=args.height_bound,
            timeout_seconds=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        source = f"initial:priority:{priority}:mask:{mask:#06x}"
        for point in outcome.curve_points:
            discoveries.setdefault(legacy.canonical_point(point), set()).add(
                source
            )
        covers.append(
            {
                "priority": priority,
                "mask": mask,
                "hex": f"0x{mask:04x}",
                "generic_representative": list(generic_representative),
                "specialized_representative": list(representative),
                "specialized_depth": str(depth),
                "search": outcome.record,
            }
        )
        print(
            f"MW16LADDERBLIND|parent={parent['parent_id']}|"
            f"chart={priority}/{len(initial_ranked)}|mask={mask:#06x}|"
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
    if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("returned points could not be classified exactly")
    quotient_rank = len(basis) - GENERIC_DIMENSION
    return {
        "parent_id": parent["parent_id"],
        "curve_id": int(parent["curve_id"]),
        "priority_rank": int(parent["priority_rank"]),
        "status": "PASS_COMPLETE_INITIAL_HALF_LATTICE_WAVE",
        "generic_rank": GENERIC_DIMENSION,
        "generic_mod2_independence_rank": GENERIC_DIMENSION,
        "generic_half_lattice": {
            "complete_class_count": EXPECTED_CLASS_COUNT,
            "twice_norm_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "maximum_twice_norm": maximum_norm,
            "maximum_depth": str(Fraction(maximum_norm, 8)),
            "deepest_class_count": len(deepest),
            "deepest_masks": [row[1] for row in deepest],
        },
        "specialized_ranking": specialized_ranking,
        "cover_records": covers,
        "discoveries": legacy.discovery_records(discoveries),
        "discovered_group_saturation": classification,
        "basis_rank_after": len(basis),
        "exact_quotient_rank_recovered": quotient_rank,
        "wall_seconds": time.monotonic() - started_wall,
        "cpu_seconds": cpu_clock() - started_cpu,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--parent-ids", default="all")
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")

    blind = json.loads(args.input.read_text())
    if blind.get("status") != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("MW16 ladder input status changed")
    if blind["observation_unit"]["curve_count"] != 5 or blind["observation_unit"]["presentation_count"] != 9:
        raise ArithmeticError("five-curve/nine-presentation partition changed")
    if blind["blindness_boundary"]["public_point_lists_loaded"] or blind["blindness_boundary"]["public_complement_coordinates_loaded"] or blind["blindness_boundary"]["target_rank_lower_bounds_loaded"]:
        raise ArithmeticError("the search input is not complement-blind")
    available = {row["parent_id"]: row for row in blind["parents"]}
    if args.parent_ids == "all":
        selected_ids = sorted(
            available,
            key=lambda key: (
                int(available[key]["curve_id"]),
                int(available[key]["priority_rank"]),
            ),
        )
    else:
        selected_ids = [value.strip() for value in args.parent_ids.split(",") if value.strip()]
        if len(selected_ids) != len(set(selected_ids)) or any(
            value not in available for value in selected_ids
        ):
            raise SystemExit("--parent-ids contains an unknown or repeated parent")

    legacy = SourceFileLoader("mw16_ladder_legacy", str(LEGACY)).load_module()
    legacy.GENERIC_DIMENSION = GENERIC_DIMENSION
    payload = {
        "schema": "elliptic-curves.icarm-mw16-parent-ladder-blind.v1",
        "status": "SEARCHING",
        "observation_unit": blind["observation_unit"],
        "blindness_boundary": {
            "sole_data_input": relative(args.input),
            "public_point_lists_loaded": False,
            "public_complement_coordinates_loaded": False,
            "target_jump_labels_loaded": False,
        },
        "declared_budget": {
            "parent_ids": selected_ids,
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "complete_generic_half_lattice_classes_each_parent": EXPECTED_CLASS_COUNT,
            "searched_chart_rule": (
                "every class in the exact maximum-depth stratum; the count is "
                "recomputed separately for each parent presentation"
            ),
            "adaptive_quotient_lifts": 0,
        },
        "parents": [],
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, LEGACY, ENGINE, Path(__file__))
        },
        "software": {"python": platform.python_version()},
        "claim_boundary": [
            "The search is blind to public complements and target jump labels.",
            "Each parent receives the same complete-maximum-depth initial MW16 rule; the exact stratum sizes may differ.",
            "The nine parent responses are nested within five target-curve observations and are not nine independent observations.",
            "Only exact returned points and exact relation classification contribute to recovered quotient rank.",
            "Every miss is bounded and has no point-absence, rank, saturation, covering, or Selmer meaning.",
            "No adaptive lift, unrestricted point search, or expensive continuation is authorized by this artifact.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage"
        ),
    }
    for parent_id in selected_ids:
        result = run_parent(available[parent_id], legacy, args)
        payload["parents"].append(result)
        write_payload(args.output, payload)
        print(
            f"MW16LADDERBLIND|parent={parent_id}|"
            f"quotient_rank={result['exact_quotient_rank_recovered']}|"
            "status=PASS_COMPLETE_INITIAL_HALF_LATTICE_WAVE",
            flush=True,
        )
    payload["status"] = "PASS_COMPLETE_NINE_PARENT_INITIAL_HALF_LATTICE_LADDER"
    write_payload(args.output, payload)
    print(
        f"MW16LADDERBLIND|parents={len(payload['parents'])}|"
        f"output={relative(args.output)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
