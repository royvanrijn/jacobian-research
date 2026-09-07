#!/usr/bin/env python3
"""Measure exact coefficient-relation consensus on the labelled R17 controls.

This is a supervised Phase-0 signal audit, not a blind selector.  It uses the
published control embeddings to place retained rays in corresponding generic
coordinates, then asks how much exact additive structure survives
leave-one-fibre-out cutoff variation.  No wgxli module is imported.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib
import json
from math import lcm
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    build_relation_complex,
    canonical_rational_unoriented,
    enumerate_short_vectors,
    exact_span_mask,
    rational_rank,
    row_basis_coordinates,
)
from calibrate_finite_aware_latent_lattice import (  # noqa: E402
    FERMIGIER_RANK20_SOURCE,
    CURVE282_SOURCE,
    public_controls,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
HEIGHT_BOUNDS = {25: 40.0, 26: 43.0, 27: 52.0, 28: 60.0}
FERMIGIER_HEIGHT_BOUNDS = {
    "ICARM_245_Fermigier_negative_control": 28.0,
    "ICARM_282_Fermigier_sibling": 36.0,
    "Fermigier_u_28917_over_20_sibling": 140.0,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def consensus_records(coefficient_sets, *, support_levels, training_support):
    """Return exact common-ray cores and leave-one-out visibility records."""

    counts = Counter(
        vector for population in coefficient_sets.values() for vector in population
    )
    support_records = []
    for support in support_levels:
        vectors = tuple(
            sorted(vector for vector, count in counts.items() if count >= support)
        )
        complex_ = build_relation_complex(vectors)
        support_records.append(
            {
                "minimum_fibre_support": support,
                "ray_count": len(vectors),
                "rank": rational_rank(vectors),
                "relation_complex": complex_.to_record(include_relations=True),
            }
        )

    labels = tuple(coefficient_sets)
    leave_one_out = []
    for held_index, held_label in enumerate(labels):
        training = [
            coefficient_sets[label]
            for index, label in enumerate(labels)
            if index != held_index
        ]
        training_counts = Counter(vector for population in training for vector in population)
        training_core = {
            vector for vector, count in training_counts.items()
            if count >= training_support
        }
        held_visible = training_core & coefficient_sets[held_label]
        leave_one_out.append(
            {
                "held_out_label": held_label,
                "training_support_threshold": training_support,
                "training_core_ray_count": len(training_core),
                "training_core_rank": rational_rank(tuple(training_core)),
                "held_out_visible_ray_count": len(held_visible),
                "held_out_visible_rank": rational_rank(tuple(held_visible)),
                "held_out_visible_relation_digest": build_relation_complex(
                    tuple(held_visible)
                ).canonical_digest,
                "training_core_vectors": [
                    list(vector) for vector in sorted(training_core)
                ],
                "held_out_visible_vectors": [
                    list(vector) for vector in sorted(held_visible)
                ],
            }
        )
    return support_records, leave_one_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    truth_document = json.loads(TRUTH.read_text())
    truth_records = {
        record["label"]: record for record in truth_document["positive_controls"]
    }
    coefficient_sets = {}
    controls = []
    for rank in range(25, 29):
        label = f"rank_at_least_{rank}"
        truth = truth_records[label]
        basis = tuple(
            tuple(map(int, row)) for row in truth["embedding_matrix_columns"]
        )
        if truth["smith_invariant_factors"] != [1] * 17:
            raise ArithmeticError("R17 calibration embedding is not primitive")
        module = importlib.import_module(f"elkies_rank{rank}")
        records = enumerate_short_vectors(
            EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS)),
            tuple(module.POINTS),
            height_bound=HEIGHT_BOUNDS[rank],
            digits=80,
            maximum_lines=100_000,
            materialize_points=False,
        )
        vectors = tuple(record.coordinates for record in records)
        mask = exact_span_mask(vectors, basis)
        coordinates = row_basis_coordinates(
            [vectors[index] for index in np.flatnonzero(mask)], basis
        )
        coefficient_sets[label] = set(coordinates)
        complex_ = build_relation_complex(coordinates)
        controls.append(
            {
                "label": label,
                "height_bound": HEIGHT_BOUNDS[rank],
                "ambient_short_vector_lines": len(records),
                "ambient_short_vector_coordinate_rows": [
                    " ".join(map(str, vector)) for vector in vectors
                ],
                "ambient_short_vector_heights": [
                    record.canonical_height for record in records
                ],
                "truth_subgroup_lines": len(coordinates),
                "truth_subgroup_line_rank": rational_rank(coordinates),
                "truth_subgroup_ternary_relations": len(
                    complex_.ternary_relations
                ),
                "truth_subgroup_relation_digest": complex_.canonical_digest,
                "truth_subgroup_coefficient_vectors": [
                    list(vector) for vector in sorted(coordinates)
                ],
            }
        )

    support_records, leave_one_out = consensus_records(
        coefficient_sets, support_levels=(4, 3, 2), training_support=2
    )
    # Preserve the established R17 field names for downstream control tools.
    for record in leave_one_out:
        record["training_two_of_three_ray_count"] = record.pop(
            "training_core_ray_count"
        )
        record["training_two_of_three_rank"] = record.pop("training_core_rank")
        record["training_two_of_three_vectors"] = record.pop(
            "training_core_vectors"
        )

    raw_truth_bases = {
        record["label"]: tuple(
            tuple(map(int, row)) for row in record["embedding_matrix_columns"]
        )
        for record in truth_document["fermigier_family_controls"]
    }
    all_controls = public_controls()
    fermigier_sets = {}
    fermigier_controls = []
    for label, bound in FERMIGIER_HEIGHT_BOUNDS.items():
        model, points = all_controls[label]
        curve = EllipticCurve(tuple(model))
        records = enumerate_short_vectors(
            curve,
            tuple(points),
            height_bound=bound,
            digits=80,
            maximum_lines=100_000,
            materialize_points=False,
        )
        vectors = tuple(record.coordinates for record in records)
        basis = raw_truth_bases[label]
        mask = exact_span_mask(vectors, basis)
        rational_coordinates = row_basis_coordinates(
            [vectors[index] for index in np.flatnonzero(mask)],
            basis,
            require_integral=False,
        )
        coordinate_denominators = []
        coordinates = []
        for vector in rational_coordinates:
            denominator = 1
            for value in vector:
                denominator = lcm(denominator, Fraction(value).denominator)
            coordinate_denominators.append(denominator)
            coordinates.append(canonical_rational_unoriented(vector))
        coordinates = tuple(coordinates)
        fermigier_sets[label] = set(coordinates)
        complex_ = build_relation_complex(coordinates)
        fermigier_controls.append(
            {
                "label": label,
                "height_bound": bound,
                "ambient_rank": len(points),
                "ambient_short_vector_lines": len(records),
                "ambient_short_vector_coordinate_rows": [
                    " ".join(map(str, vector)) for vector in vectors
                ],
                "ambient_short_vector_heights": [
                    record.canonical_height for record in records
                ],
                "primitive_truth_subgroup_lines": len(coordinates),
                "reference_subgroup_coordinate_denominator_histogram": [
                    {"denominator": denominator, "count": count}
                    for denominator, count in sorted(
                        Counter(coordinate_denominators).items()
                    )
                ],
                "primitive_truth_subgroup_line_rank": rational_rank(coordinates),
                "primitive_truth_subgroup_ternary_relations": len(
                    complex_.ternary_relations
                ),
                "primitive_truth_subgroup_relation_digest": complex_.canonical_digest,
                "primitive_truth_subgroup_coefficient_vectors": [
                    list(vector) for vector in sorted(coordinates)
                ],
            }
        )
    sibling_labels = (
        "ICARM_282_Fermigier_sibling",
        "Fermigier_u_28917_over_20_sibling",
    )
    sibling_vectors = tuple(
        sorted(fermigier_sets[sibling_labels[0]] & fermigier_sets[sibling_labels[1]])
    )
    sibling_complex = build_relation_complex(sibling_vectors)
    fermigier_sibling_pair = {
        "labels": list(sibling_labels),
        "common_rational_ray_count": len(sibling_vectors),
        "common_rational_ray_rank": rational_rank(sibling_vectors),
        "relation_complex": sibling_complex.to_record(include_relations=True),
    }
    fermigier245 = next(
        record
        for record in fermigier_controls
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )

    passed_r17 = all(record["held_out_visible_rank"] == 17 for record in leave_one_out)
    passed_fermigier = (
        fermigier_sibling_pair["common_rational_ray_rank"] == 12
        and fermigier245["primitive_truth_subgroup_line_rank"] == 12
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-relation-consensus.v1",
        "status": (
            "PASS_CONTROL_EXACT_RELATION_SIGNALS"
            if passed_r17 and passed_fermigier
            else "FAIL_TWO_FAMILY_RELATION_SIGNAL"
        ),
        "scope": "Supervised Phase-0 R17 controls only; no wgxli target is loaded",
        "role": (
            "Exact recoverability benchmark for a future basis-free hypergraph matcher; "
            "published embeddings align the controls in this artifact and therefore "
            "this artifact is not itself a blind recovery."
        ),
        "controls": controls,
        "coefficient_support_cores": support_records,
        "leave_one_fibre_out": leave_one_out,
        "fermigier_controls": fermigier_controls,
        "fermigier_sibling_pair": fermigier_sibling_pair,
        "fermigier_245_negative_truth_audit": fermigier245,
        "gate_decision": (
            "CLOSED. Full-rank exact held-out R17 signal and rank-12 Fermigier "
            "rational-ray signals exist, but the ICARM 245 truth audit is supervised "
            "and a basis-free unequal-cloud matcher has not recovered it blindly."
        ),
        "proof_boundary": (
            "Integer and rational-ray coordinate conversion, denominator clearing, "
            "primitivity checks, coefficient-set intersections, ranks, and additive "
            "relations are exact. Canonical-height cutoffs are numerical. The R17 "
            "and Fermigier reference bases are supplied by the withheld control "
            "artifact, not discovered by this audit. ICARM 245 is a distinct "
            "Fermigier--Mestre family from the aligned 282/u sibling pair and is not "
            "asserted to share their labelled rays."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                TRUTH,
                CURVE282_SOURCE,
                FERMIGIER_RANK20_SOURCE,
                ELLIPTIC / "cas/icarm_curve245.py",
                ELLIPTIC / "cas/calibrate_finite_aware_latent_lattice.py",
                *(ELLIPTIC / "cas" / f"elkies_rank{rank}.py" for rank in range(25, 29)),
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice relation-consensus artifact is stale")
        print(
            f"LATENTRELATION|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTRELATION|status={payload['status']}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
