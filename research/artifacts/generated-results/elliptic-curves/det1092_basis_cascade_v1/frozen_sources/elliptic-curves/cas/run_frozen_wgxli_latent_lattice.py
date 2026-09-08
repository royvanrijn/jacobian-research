#!/usr/bin/env python3
"""Run the frozen latent-lattice method on the five hash-pinned wgxli fibres.

The algorithm is defined by ``LATENT-LATTICE-WGXLI-FROZEN-2026-09-01-v1``.
This adapter may parse target records and serialize results, but it must not
change a frozen search or scoring parameter after target inspection.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from analyze_icarm_wgxli_rank17_lineages import (  # noqa: E402
    TARGET_CURVE_IDS,
    TARGET_SOURCES,
    fetch_bytes,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    beam_subspace_scan,
    build_relation_complex,
    enumerate_short_vectors,
    row_embedding_smith_invariant_factors,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
FREEZE = ARTIFACTS / "latent_lattice_target_method_freeze_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_wgxli_frozen_dimension_v1.json"
EXPECTED_FREEZE_SHA256 = "ef6f8b7be7a14095efa7529fb795d237e06465ba1cec023dcb4845287609c9f4"
EXPECTED_TAG = "LATENT-LATTICE-WGXLI-FROZEN-2026-09-01-v1"
OUTPUT_SCHEMA = "elliptic-curves.latent-lattice-wgxli-frozen-dimension.v1"
ADAPTER_PATH = Path(__file__).resolve()
HEIGHT_BOUNDS = tuple(range(20, 82, 2))
MINIMUM_RAYS = 1_800
MAXIMUM_RAYS = 100_000
PERSISTENCE_FRACTIONS = ((3, 4), (7, 8), (1, 1))
DIMENSIONS = tuple(range(10, 21))
POOL = 300
BEAM_WIDTH = 8
BRANCH_WIDTH = 80
HEIGHT_DIGITS = 80


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_freeze() -> dict[str, object]:
    observed = digest(FREEZE)
    if observed != EXPECTED_FREEZE_SHA256:
        raise SystemExit(f"frozen manifest changed: {observed}")
    document = json.loads(FREEZE.read_text())
    if document["algorithm_tag"] != EXPECTED_TAG:
        raise SystemExit("frozen algorithm tag changed")
    expected = document["cloud_protocol"]
    if expected["height_candidate_bounds"] != list(HEIGHT_BOUNDS):
        raise SystemExit("adapter height bounds differ from freeze")
    dimension = document["dimension_protocol"]
    if dimension["candidate_dimensions"] != list(DIMENSIONS):
        raise SystemExit("adapter dimensions differ from freeze")
    if dimension["beam_subspace_scan"] != {
        "pool": POOL,
        "beam_width": BEAM_WIDTH,
        "branch_width": BRANCH_WIDTH,
        "seed_rule": "10000 + ICARM record number",
    }:
        raise SystemExit("adapter beam bounds differ from freeze")
    return document


def point_pair(raw) -> tuple[Fraction, Fraction]:
    return Fraction(raw[0]), Fraction(raw[1])


def dimension_persistence(score_tables):
    persistence = {}
    for dimension in DIMENSIONS:
        neighbours = tuple(
            candidate
            for candidate in (dimension - 1, dimension, dimension + 1)
            if candidate in DIMENSIONS
        )
        persistence[dimension] = min(
            float(table[candidate]["integrality_likelihood_ratio"])
            for table in score_tables
            for candidate in neighbours
        )
    selected = max(
        DIMENSIONS,
        key=lambda dimension: (
            persistence[dimension],
            min(
                float(table[dimension]["integrality_likelihood_ratio"])
                for table in score_tables
            ),
            -dimension,
        ),
    )
    return selected, persistence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    freeze = validate_freeze()
    partial_output = args.output.with_suffix(args.output.suffix + ".partial")
    if args.check:
        if not args.output.exists():
            raise SystemExit("frozen wgxli dimension artifact is missing")
        document = json.loads(args.output.read_text())
        if document.get("algorithm_tag") != EXPECTED_TAG:
            raise SystemExit("dimension artifact has the wrong algorithm tag")
        if document.get("freeze_sha256") != EXPECTED_FREEZE_SHA256:
            raise SystemExit("dimension artifact has the wrong freeze hash")
        if document.get("target_source_sha256") != {
            str(key): value for key, value in TARGET_SOURCES.items()
        }:
            raise SystemExit("dimension artifact target hashes changed")
        print(f"WGXLIDIMENSION|check=PASS|sha256={digest(args.output)}")
        return

    targets = {
        curve_id: json.loads(
            fetch_bytes(
                f"https://elliptic-rank.icarm.cloud/curve/{curve_id}.json",
                TARGET_SOURCES[curve_id],
            )
        )
        for curve_id in TARGET_CURVE_IDS
    }
    results = []
    started = time.monotonic()
    for curve_id in TARGET_CURVE_IDS:
        record = targets[curve_id]
        curve = EllipticCurve(tuple(Fraction(value) for value in record["ainvs"]))
        points = tuple(point_pair(point) for point in record["points"])
        chosen_bound = None
        count_ledger = []
        for bound in HEIGHT_BOUNDS:
            short = enumerate_short_vectors(
                curve,
                points,
                height_bound=bound,
                digits=HEIGHT_DIGITS,
                maximum_lines=MAXIMUM_RAYS,
                materialize_points=False,
            )
            count_ledger.append({"height_bound": bound, "primitive_unoriented_ray_count": len(short)})
            if len(short) >= MAXIMUM_RAYS:
                raise ArithmeticError(f"curve {curve_id} reached the frozen ray cap")
            if len(short) >= MINIMUM_RAYS:
                chosen_bound = bound
                break
        if chosen_bound is None:
            raise ArithmeticError(f"curve {curve_id} did not reach the frozen minimum ray count")
        short = enumerate_short_vectors(
            curve,
            points,
            height_bound=chosen_bound,
            digits=HEIGHT_DIGITS,
            maximum_lines=MAXIMUM_RAYS,
            materialize_points=True,
        )
        cutoff_results = []
        score_tables = []
        for numerator, denominator in PERSISTENCE_FRACTIONS:
            cutoff = numerator * len(short) // denominator
            population = short[:cutoff]
            complex_ = build_relation_complex(
                tuple(item.coordinates for item in population),
                tuple(
                    {
                        "integral": bool(item.arithmetic["integral"]),
                        "total_bits": int(item.arithmetic["total_bits"]),
                    }
                    for item in population
                ),
            )
            candidates = beam_subspace_scan(
                population,
                complex_,
                dimensions=DIMENSIONS,
                pool=POOL,
                beam_width=BEAM_WIDTH,
                branch_width=BRANCH_WIDTH,
                seed=10_000 + curve_id,
            )
            by_dimension = {candidate.dimension: candidate.to_record() for candidate in candidates}
            missing_dimensions = sorted(set(DIMENSIONS) - set(by_dimension))
            for candidate in candidates:
                dimension = candidate.dimension
                smith = row_embedding_smith_invariant_factors(candidate.primitive_basis_rows)
                by_dimension[dimension]["smith_invariant_factors"] = list(smith)
                by_dimension[dimension]["primitive_in_public_subgroup"] = all(value == 1 for value in smith)
            score_tables.append(by_dimension)
            cutoff_results.append(
                {
                    "fraction": f"{numerator}/{denominator}",
                    "ray_count": len(population),
                    "ternary_relation_count": len(complex_.ternary_relations),
                    "scaled_relation_count": len(complex_.scaled_relations),
                    "missing_dimensions": missing_dimensions,
                    "dimension_candidates": [by_dimension[dimension] for dimension in sorted(by_dimension)],
                }
            )
        all_dimensions_present = all(set(table) == set(DIMENSIONS) for table in score_tables)
        if all_dimensions_present:
            selected_dimension, persistence = dimension_persistence(score_tables)
        else:
            selected_dimension, persistence = None, {}
        results.append(
            {
                "curve_id": curve_id,
                "ambient_rank": len(points),
                "rank_lower_bound": int(record["rank_lower_bound"]),
                "adaptive_height_count_ledger": count_ledger,
                "selected_height_bound": chosen_bound,
                "full_ray_count": len(short),
                "cutoffs": cutoff_results,
                "selected_dimension": selected_dimension,
                "persistence_scores": [
                    {
                        "dimension": dimension,
                        "cross_cutoff_three_level_bottleneck_integrality_llr": f"{persistence[dimension]:.17g}",
                    }
                    for dimension in sorted(persistence)
                ],
                "selected_embedding_matrix_rows_by_cutoff": [
                    score_tables[index][selected_dimension]["primitive_basis_rows"]
                    for index in range(len(score_tables))
                ] if selected_dimension is not None else [],
            }
        )
        partial_output.write_text(
            json.dumps(
                {
                    "algorithm_tag": EXPECTED_TAG,
                    "freeze_sha256": EXPECTED_FREEZE_SHA256,
                    "completed_fibres": results,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(
            f"WGXLIDIMENSIONPROGRESS|curve={curve_id}|ambient={len(points)}|bound={chosen_bound}|"
            f"rays={len(short)}|selected_k={selected_dimension if selected_dimension is not None else 'FAIL'}|"
            f"missing={','.join(map(str, sorted(set(DIMENSIONS)-set.intersection(*(set(table) for table in score_tables))))) or 'none'}|"
            f"seconds={time.monotonic()-started:.3f}",
            flush=True,
        )

    selected_dimensions = [record["selected_dimension"] for record in results]
    counts = {dimension: selected_dimensions.count(dimension) for dimension in DIMENSIONS}
    recurring = [dimension for dimension, count in counts.items() if count >= 4]
    status = (
        "PASS_FROZEN_DIMENSION_RECURRENCE"
        if len(recurring) == 1
        else "FAIL_FROZEN_DIMENSION_RECURRENCE"
    )
    payload = {
        "schema": OUTPUT_SCHEMA,
        "algorithm_tag": EXPECTED_TAG,
        "freeze_sha256": EXPECTED_FREEZE_SHA256,
        "status": status,
        "target_source_sha256": {str(key): value for key, value in TARGET_SOURCES.items()},
        "selected_dimensions": selected_dimensions,
        "dimensions_recurring_in_at_least_four_fibres": recurring,
        "fibres": results,
        "proof_boundary": (
            "Public source hashes, point membership used by height enumeration, primitive coordinates, "
            "relations, saturation indices, and replayed cloud membership are exact within the frozen "
            "bounds. Canonical heights and integrality-LLR persistence are numerical/statistical. The "
            "beam scan is bounded and does not exhaust all primitive subspaces."
        ),
        "inputs": {
            str(FREEZE.relative_to(ROOT)): digest(FREEZE),
            str(ADAPTER_PATH.relative_to(ROOT)): digest(ADAPTER_PATH),
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"WGXLIDIMENSION|status={status}|selected={','.join('FAIL' if value is None else str(value) for value in selected_dimensions)}|"
        f"recurring={','.join(map(str, recurring)) or 'none'}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
