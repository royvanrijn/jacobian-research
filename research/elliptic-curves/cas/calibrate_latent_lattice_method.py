#!/usr/bin/env python3
"""Calibrate the bounded latent-lattice selector before any target search.

The selector never receives a generic section or embedding.  Exact withheld
embeddings are loaded only after selection to score the returned rational
subspace.  Failure closes the gate on the unknown wgxli cluster.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from icarm_curve245 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as CURVE245_MODEL,
    POINTS as CURVE245_POINTS,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    beam_subspace_scan,
    build_relation_complex,
    enumerate_short_vectors,
    exact_span_mask,
    finite_quotient_block,
    multiplicative_component_block,
    rational_rank,
    relation_seeded_subspace_scan,
)
from latent_lattice.subspace import integrality_llr  # noqa: E402


TRUTH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_truth_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_v2.json"
)
POSITIVE_BOUNDS = {25: 40.0, 26: 43.0, 27: 52.0, 28: 60.0}
DIMENSIONS = tuple(range(10, 21))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def overlap_dimension(left, right) -> int:
    return len(left) + len(right) - rational_rank(tuple(left) + tuple(right))


def positive_controls(truth: dict[str, object]) -> list[dict[str, object]]:
    answers = []
    for rank, bound in POSITIVE_BOUNDS.items():
        module = importlib.import_module(f"elkies_rank{rank}")
        curve = EllipticCurve(module.GENERAL_WEIERSTRASS_COEFFICIENTS)
        records = enumerate_short_vectors(
            curve,
            module.POINTS,
            height_bound=bound,
            digits=60,
            maximum_lines=10_000,
            materialize_points=False,
        )
        complex_ = build_relation_complex([record.coordinates for record in records])
        candidate = beam_subspace_scan(
            records,
            complex_,
            dimensions=(17,),
            pool=300,
            beam_width=8,
            branch_width=80,
            seed=1700 + rank,
        )[0]
        withheld = next(
            record
            for record in truth["positive_controls"]
            if record["label"] == f"rank_at_least_{rank}"
        )
        true_basis = tuple(map(tuple, withheld["embedding_matrix_columns"]))
        intersection = overlap_dimension(candidate.primitive_basis_rows, true_basis)
        answers.append(
            {
                "label": f"rank_at_least_{rank}",
                "ambient_rank": rank,
                "height_bound": bound,
                "short_vector_lines": len(records),
                "complete_relation_edges": len(complex_.ternary_relations),
                "relation_invariant_digest": complex_.canonical_digest,
                "selected_dimension": 17,
                "selected_support": candidate.support,
                "selected_relation_mass": candidate.relation_mass,
                "selected_primitive_basis_rows": [
                    list(row) for row in candidate.primitive_basis_rows
                ],
                "withheld_truth_intersection_dimension": intersection,
                "withheld_truth_recovered_exactly": intersection == 17,
            }
        )
    return answers


def negative_control(truth: dict[str, object]) -> dict[str, object]:
    curve = EllipticCurve(CURVE245_MODEL)
    records = enumerate_short_vectors(
        curve,
        CURVE245_POINTS,
        height_bound=28.0,
        digits=80,
        maximum_lines=100_000,
        materialize_points=True,
    )
    metadata = [
        {
            "integral": bool(record.arithmetic["integral"]),
            "x_denominator_bits": int(record.arithmetic["x_denominator_bits"]),
            "y_denominator_bits": int(record.arithmetic["y_denominator_bits"]),
        }
        for record in records
    ]
    complex_ = build_relation_complex(
        [record.coordinates for record in records], metadata
    )
    scan = beam_subspace_scan(
        records,
        complex_,
        dimensions=DIMENSIONS,
        pool=300,
        beam_width=8,
        branch_width=80,
        seed=24512,
    )
    truth_rows = truth["negative_control"][
        "primitive_closure_embedding_matrix_rows"
    ]
    true_basis = tuple(map(tuple, zip(*truth_rows)))
    vectors = [record.coordinates for record in records]
    integral = np.asarray(
        [bool(record.arithmetic["integral"]) for record in records], dtype=bool
    )
    degree_by_vector = dict(zip(complex_.vertices, complex_.additive_degrees))
    degrees = np.asarray([degree_by_vector[record.coordinates] for record in records])
    true_mask = exact_span_mask(vectors, true_basis)
    true_stats = {
        "dimension": 12,
        "support": int(true_mask.sum()),
        "integral_support": int(np.sum(true_mask & integral)),
        "integrality_likelihood_ratio": integrality_llr(true_mask, integral),
        "relation_mass": int(np.sum(degrees[true_mask])),
    }
    scan_records = []
    for candidate in scan:
        scan_records.append(
            {
                **candidate.to_record(),
                "withheld_truth_intersection_dimension": overlap_dimension(
                    candidate.primitive_basis_rows, true_basis
                ),
            }
        )
    selected = max(scan, key=lambda candidate: candidate.integrality_llr)

    relation_seeded = relation_seeded_subspace_scan(
        records,
        complex_,
        dimension=12,
        seed_edges=3_000,
        beam_width=48,
        branch_width=100,
        candidates=1,
    )[0]

    finite_blocks = []
    for relation_prime, reduction_primes in (
        (2, (11, 23, 29, 41)),
        (3, (11, 29, 41)),
    ):
        for reduction_prime in reduction_primes:
            block = finite_quotient_block(
                curve,
                CURVE245_POINTS,
                reduction_prime,
                relation_prime,
            )
            classes = Counter(
                block.vector_class(record.coordinates) for record in records
            )
            finite_blocks.append(
                {
                    **block.to_record(),
                    "short_vector_class_histogram": [
                        {"class": list(key), "count": value}
                        for key, value in sorted(classes.items())
                    ],
                }
            )

    component_blocks = []
    for prime, order, split in (
        (2, 17, True),
        (5, 4, True),
        (13, 4, True),
        (19, 5, True),
        (37, 2, False),
    ):
        block = multiplicative_component_block(
            curve,
            CURVE245_POINTS,
            prime=prime,
            fibre_order=order,
            split=split,
        )
        classes = Counter(block.vector_class(record.coordinates) for record in records)
        component_blocks.append(
            {
                **block.to_record(),
                "short_vector_class_histogram": [
                    {"class": key, "count": value}
                    for key, value in sorted(classes.items())
                ],
            }
        )

    return {
        "label": "ICARM_245_Fermigier_negative_control",
        "ambient_rank": 20,
        "height_bound": 28.0,
        "short_vector_lines": len(records),
        "complete_relation_edges": len(complex_.ternary_relations),
        "relation_invariant_digest": complex_.canonical_digest,
        "dimension_scan": scan_records,
        "dimension_selected_by_max_integrality_llr": selected.dimension,
        "selected_truth_intersection_dimension": overlap_dimension(
            selected.primitive_basis_rows, true_basis
        ),
        "selected_truth_recovered_exactly": False,
        "relation_seeded_rank12": {
            **relation_seeded.to_record(),
            "withheld_truth_intersection_dimension": overlap_dimension(
                relation_seeded.primitive_basis_rows, true_basis
            ),
        },
        "withheld_truth_stats": true_stats,
        "finite_good_reduction_quotient_blocks": finite_blocks,
        "multiplicative_component_blocks": component_blocks,
        "primitive_vector_divisible_by_2_or_3_in_displayed_Z_basis": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    truth = json.loads(TRUTH.read_text())
    positives = positive_controls(truth)
    negative = negative_control(truth)
    exact_positive_count = sum(
        record["withheld_truth_recovered_exactly"] for record in positives
    )
    payload = {
        "schema": "elliptic-curves.latent-lattice-calibration.v2",
        "status": "FAIL_CALIBRATION_TARGET_GATE_CLOSED",
        "algorithm": {
            "short_vector_population": (
                "complete primitive unoriented enumeration through each declared "
                "canonical-height bound in the full displayed subgroup"
            ),
            "relation_structure": (
                "complete unoriented unit-content ternary additive hypergraph; "
                "non-unit a+/-b=m*c relations retain multiplier and source/target "
                "roles; coordinate-free four-round color-refinement invariant"
            ),
            "selector": (
                "bounded exact-saturation beam search with pool=300, beam=8, "
                "branch=80; negative-control second pass uses 3000 relation seeds, "
                "beam=48, branch=100"
            ),
            "truth_policy": "withheld until a blind candidate has been selected",
        },
        "truth_artifact": {
            "path": str(TRUTH.relative_to(ROOT)),
            "sha256": digest(TRUTH),
        },
        "positive_controls": positives,
        "positive_controls_recovered_exactly": exact_positive_count,
        "positive_controls_total": len(positives),
        "negative_control": negative,
        "gate_decision": (
            "Calibration requires recovery on several R17 controls and recovery "
            "of the Fermigier rank-12 space. Those conditions fail, so no search "
            "on curves 351,356,376,377,385 is authorized by this calibrated method."
        ),
        "proof_boundary": (
            "Enumerated vectors, additive relations, finite quotient codes, component "
            "codes, exact span intersections, saturations, and supplied point identities "
            "are exact within the declared bounds. Canonical heights and selector scores "
            "are numerical. Failure to recover withheld controls proves only that this "
            "bounded method is inadequate; it is not a nonexistence theorem for a common "
            "generic lattice."
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(f"PASS|{args.output}|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTCAL|status={payload['status']}|positive={exact_positive_count}/4|"
        f"negative_intersection={negative['selected_truth_intersection_dimension']}/12|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
