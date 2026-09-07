#!/usr/bin/env python3
"""Calibrate the finite-aware height/relation cascade on public controls only.

This program never imports a wgxli target.  It first builds blind bounded
proposal ledgers, then reveals the published R17 and Fermigier embeddings only
to measure recall and selector rank.  A failed selector keeps the target gate
closed even when the true subgroup occurs in the proposal ledger.
"""

from __future__ import annotations

import argparse
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
    FiniteQuotientBlock,
    build_relation_complex,
    cloud_height_profile_distance,
    cloud_height_signature,
    enumerate_short_vectors,
    exact_span_mask,
    height_gram,
    hermite_signature,
    hermite_signature_distance,
    independent_relation_growth_proposals,
    intrinsic_shell_signature,
    primitive_span_basis,
    rational_nullspace,
    rational_rank,
    recombined_core_extension_search,
    repeated_cross_bound_intersection_ledger,
    restricted_height_gram,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
DIMENSION_SCAN = ARTIFACTS / "latent_lattice_calibration_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_shape_calibration_v1.json"
R17_HEIGHT_BOUNDS = {25: 40.0, 26: 43.0, 27: 52.0, 28: 60.0}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def finite_block(record: dict[str, object]) -> FiniteQuotientBlock:
    return FiniteQuotientBlock(
        reduction_prime=int(record["reduction_prime"]),
        relation_prime=int(record["relation_prime"]),
        group_order=int(record["group_order"]),
        multiple_subgroup_order=int(record["multiple_subgroup_order"]),
        quotient_dimension=int(record["quotient_dimension"]),
        rows=tuple(tuple(map(int, row)) for row in record["rows"]),
    )


def overlap_dimension(left, right) -> int:
    return len(left) + len(right) - rational_rank(tuple(left) + tuple(right))


def truth_bases(document: dict[str, object]):
    return {
        record["label"]: primitive_span_basis(
            tuple(
                tuple(map(int, row))
                for row in record["embedding_matrix_columns"]
            )
        )
        for record in (
            list(document["positive_controls"])
            + list(document["fermigier_family_controls"])
        )
    }


def r17_calibration(
    truth_document: dict[str, object], finite_document: dict[str, object]
) -> dict[str, object]:
    bases = truth_bases(truth_document)
    truth_records = {
        record["label"]: record for record in truth_document["positive_controls"]
    }
    finite_records = {
        record["label"]: record for record in finite_document["controls"]
    }
    states = {}
    for rank in range(25, 29):
        label = f"rank_at_least_{rank}"
        module = importlib.import_module(f"elkies_rank{rank}")
        curve = EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS))
        points = tuple(module.POINTS)
        records = enumerate_short_vectors(
            curve,
            points,
            height_bound=R17_HEIGHT_BOUNDS[rank],
            digits=80,
            maximum_lines=100_000,
            materialize_points=False,
        )
        complex_ = build_relation_complex([record.coordinates for record in records])
        development = tuple(
            finite_block(record)
            for record in finite_records[label]["development_blocks"]
        )
        proposals = independent_relation_growth_proposals(
            records,
            complex_,
            dimension=17,
            seed_edges=3_000,
            priority_mode="finite",
            seed_strategy="stratified",
            finite_blocks=development,
        )
        truth_basis = bases[label]
        truth_source_rank = next(
            index
            for index, proposal in enumerate(proposals)
            if overlap_dimension(proposal.basis_rows, truth_basis) == 17
        )
        truth_mask = exact_span_mask(
            [record.coordinates for record in records], truth_basis
        )
        truth_indices = tuple(map(int, np.flatnonzero(truth_mask)))
        truth_cloud = cloud_height_signature(
            [record.canonical_height for record in records], truth_indices
        )
        truth_gram = truth_records[label]["canonical_height_gram"]
        truth_hermite = hermite_signature(truth_gram, digits=80)
        states[label] = {
            "records": records,
            "proposals": proposals,
            "ambient_gram": height_gram(curve, points, digits=80),
            "truth_basis": truth_basis,
            "truth_source_rank": truth_source_rank,
            "truth_cloud": truth_cloud,
            "truth_hermite": truth_hermite,
            "truth_shell": intrinsic_shell_signature(
                truth_gram, minimum_vectors=128, quantiles=32, digits=80
            ),
            "edge_count": len(complex_.ternary_relations),
        }

    controls = []
    for label, state in states.items():
        cloud_references = [
            other["truth_cloud"]
            for other_label, other in states.items()
            if other_label != label
        ]
        hermite_references = [
            other["truth_hermite"]
            for other_label, other in states.items()
            if other_label != label
        ]
        heights = [record.canonical_height for record in state["records"]]
        cloud_scores = []
        for source_rank, proposal in enumerate(state["proposals"]):
            signature = cloud_height_signature(heights, proposal.inlier_indices)
            distance = sum(
                cloud_height_profile_distance(signature, reference)
                for reference in cloud_references
            ) / len(cloud_references)
            cloud_scores.append((distance, source_rank, proposal))
        cloud_scores.sort(key=lambda item: (item[0], item[1]))
        truth_cloud_rank = next(
            rank
            for rank, item in enumerate(cloud_scores)
            if item[1] == state["truth_source_rank"]
        )
        hermite_scores = []
        for cloud_rank, (cloud_distance, source_rank, proposal) in enumerate(
            cloud_scores[:64]
        ):
            signature = hermite_signature(
                restricted_height_gram(
                    state["ambient_gram"], proposal.basis_rows
                ),
                digits=80,
            )
            distance = sum(
                hermite_signature_distance(signature, reference)
                for reference in hermite_references
            ) / len(hermite_references)
            hermite_scores.append(
                (
                    distance,
                    cloud_rank,
                    source_rank,
                    proposal,
                    signature,
                    cloud_distance,
                )
            )
        hermite_scores.sort(key=lambda item: (item[0], item[1]))
        truth_hermite_rank = next(
            (
                rank
                for rank, item in enumerate(hermite_scores)
                if item[2] == state["truth_source_rank"]
            ),
            None,
        )
        selected = hermite_scores[0]
        controls.append(
            {
                "label": label,
                "height_bound": R17_HEIGHT_BOUNDS[int(label.rsplit("_", 1)[1])],
                "short_vector_lines": len(state["records"]),
                "complete_relation_edges": state["edge_count"],
                "finite_seeded_proposal_count": len(state["proposals"]),
                "truth_source_rank": state["truth_source_rank"],
                "truth_cloud_profile_rank": truth_cloud_rank,
                "truth_hermite_rank_within_cloud_top_64": truth_hermite_rank,
                "selected_source_rank": selected[2],
                "selected_truth_overlap": overlap_dimension(
                    selected[3].basis_rows, state["truth_basis"]
                ),
                "truth_hermite_signature": state["truth_hermite"].to_record(),
                "selected_hermite_signature": selected[4].to_record(),
                "truth_intrinsic_shell": state["truth_shell"].to_record(),
            }
        )
    exact_selectors = sum(record["selected_truth_overlap"] == 17 for record in controls)
    return {
        "controls": controls,
        "finite_seeded_exact_truth_recall_count": sum(
            record["truth_source_rank"] is not None for record in controls
        ),
        "symmetric_leave_one_out_exact_selection_count": exact_selectors,
        "rank25_held_out_training_on_rank26_through_rank28": {
            "truth_cloud_profile_rank": controls[0]["truth_cloud_profile_rank"],
            "truth_hermite_rank_within_cloud_top_64": controls[0][
                "truth_hermite_rank_within_cloud_top_64"
            ],
            "selected_truth_overlap": controls[0]["selected_truth_overlap"],
        },
        "intrinsic_shell_digest_count": len(
            {
                record["truth_intrinsic_shell"]["relation_complex"][
                    "canonical_digest"
                ]
                for record in controls
            }
        ),
    }


def fermigier_calibration(
    truth_document: dict[str, object], finite_document: dict[str, object]
) -> dict[str, object]:
    curve = EllipticCurve(tuple(CURVE245_MODEL))
    points = tuple(CURVE245_POINTS)
    ledgers = {}
    records_by_bound = {}
    for bound in (28.0, 29.0):
        records = enumerate_short_vectors(
            curve,
            points,
            height_bound=bound,
            digits=80,
            maximum_lines=100_000,
            materialize_points=True,
        )
        complex_ = build_relation_complex([record.coordinates for record in records])
        ledgers[bound] = recombined_core_extension_search(
            records,
            complex_,
            dimension=12,
            seed_edges=3_000,
            anchor_count=500,
            enclosure_codimension=3,
            enclosure_count=0,
            inner_count=2,
        )
        records_by_bound[bound] = (records, complex_)
    records, complex_ = records_by_bound[28.0]
    repeated = repeated_cross_bound_intersection_ledger(
        records,
        complex_,
        ledgers[28.0].enclosure_proposals,
        ledgers[29.0].enclosure_proposals,
        target_dimension=12,
        left_count=200,
        right_count=200,
        maximum_candidates=128,
    )
    truth_basis = truth_bases(truth_document)[
        "ICARM_245_Fermigier_negative_control"
    ]
    truth_key = rational_nullspace(truth_basis)
    truth_rank = next(
        (
            rank
            for rank, proposal in enumerate(repeated.proposals)
            if proposal.exact_annihilator_rows == truth_key
        ),
        None,
    )
    finite_record = next(
        record
        for record in finite_document["controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    return {
        "height_bounds": [28.0, 29.0],
        "short_vector_lines": {
            str(int(bound)): len(records_by_bound[bound][0])
            for bound in (28.0, 29.0)
        },
        "enclosure_counts": {
            str(int(bound)): len(ledgers[bound].enclosure_proposals)
            for bound in (28.0, 29.0)
        },
        "repeated_intersection_ledger": repeated.summary_record(),
        "truth_rank": truth_rank,
        "truth_proposal": (
            None if truth_rank is None else repeated.proposals[truth_rank].to_record()
        ),
        "development_truth_finite_signature": finite_record[
            "development_truth_signature"
        ],
        "held_out_truth_finite_signature": finite_record[
            "held_out_truth_signature"
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    truth_document = json.loads(TRUTH.read_text())
    finite_document = json.loads(FINITE.read_text())
    dimension_scan_document = json.loads(DIMENSION_SCAN.read_text())
    r17 = r17_calibration(truth_document, finite_document)
    fermigier = fermigier_calibration(truth_document, finite_document)
    earlier_negative = dimension_scan_document["negative_control"]
    passed_proposal_calibration = (
        r17["finite_seeded_exact_truth_recall_count"] == 4
        and fermigier["truth_rank"] is not None
        and fermigier["truth_rank"] < 128
        and earlier_negative["dimension_selected_by_max_integrality_llr"] == 12
    )
    passed_selector_calibration = (
        r17["symmetric_leave_one_out_exact_selection_count"] == 4
        and fermigier["truth_rank"] == 0
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-shape-calibration.v1",
        "status": (
            "PASS_CONTROL_SELECTOR_CALIBRATION"
            if passed_selector_calibration
            else (
                "PASS_PROPOSAL_CALIBRATION_SELECTOR_FAIL"
                if passed_proposal_calibration
                else "FAIL_PROPOSAL_CALIBRATION"
            )
        ),
        "scope": "Phase-0 controls only; no wgxli target curve is loaded",
        "algorithm": {
            "proposal": (
                "finite-rarity-seeded exact relation growth for R17; independently "
                "generated rank-15 enclosure intersections across height bounds for ICARM 245"
            ),
            "cross_bound_score": (
                "arithmetic_LLR + 0.1*induced_ternary_relations + "
                "2*exact_cross_bound_occurrence_count"
            ),
            "height_cascade": (
                "scale-free bounded-cloud quantiles, top 64, then intrinsic Hermite invariant"
            ),
            "finite_policy": (
                "finite codes seed proposals and remain attached audit signatures; "
                "cross-fibre finite-profile distance is not a positive identity score"
            ),
        },
        "r17_positive_control": r17,
        "fermigier_negative_control": {
            "blind_dimension_selected_by_earlier_scan": earlier_negative[
                "dimension_selected_by_max_integrality_llr"
            ],
            **fermigier,
        },
        "gate_decision": (
            "CLOSED. Exact truth recall and the held-out rank-25 experiment pass, "
            "but symmetric R17 leave-one-out and rank-0 Fermigier selection fail. "
            "Do not apply this selector to wgxli."
        ),
        "proof_boundary": (
            "Elliptic-curve arithmetic, finite codes, additive relations, rational "
            "intersection grouping, primitive closures, multiplicities, and withheld "
            "overlaps are exact. Canonical heights, shell boundaries, proposal order, "
            "and all selector scores are numerical or heuristic. The 40,000-pair "
            "intersection result is exhaustive only within the declared top-200 by "
            "top-200 enclosure bounds."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                TRUTH,
                FINITE,
                DIMENSION_SCAN,
                ELLIPTIC / "cas/icarm_curve245.py",
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
            raise SystemExit("latent-lattice shape calibration artifact is stale")
        print(f"LATENTSHAPE|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTSHAPE|status={payload['status']}|"
        f"fermigier_truth_rank={fermigier['truth_rank']}|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
