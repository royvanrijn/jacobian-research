#!/usr/bin/env python3
"""Calibrate the exact graph-walk selector on R17 and Fermigier controls."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from icarm_curve245 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as FERMIGIER_MODEL,
    POINTS as FERMIGIER_POINTS,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    exact_graph_walk_consensus,
    height_gram,
    primitive_hermite_signatures,
    primitive_span_basis,
    rational_nullspace,
    rescore_graph_walk_consensus,
    row_basis_coordinates,
    row_embedding_smith_invariant_factors,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
R17_CANDIDATES = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
R17_SHAPES = ARTIFACTS / "latent_lattice_joint_shape_ledger_v1.json.gz"
R17_DIAGNOSTIC = ARTIFACTS / "latent_lattice_joint_fingerprints_v1.json"
FERMIGIER_REPLAY = ARTIFACTS / "latent_lattice_fermigier_replay_v1.json.gz"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
DIMENSION_SCAN = ARTIFACTS / "latent_lattice_calibration_v2.json"
PINNED_R17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ARTIFACTS / "latent_lattice_graph_walk_calibration_v1.json"

DEFAULT_POOL = 64
DEFAULT_GAP = 0.005
DEFAULT_GRAPH_WEIGHT = 1.5
STABILITY_POOLS = (64, 80, 96, 112, 128)
STABILITY_GAPS = (0.004, 0.005, 0.006, 0.007)
STABILITY_WEIGHTS = (1.25, 1.5, 1.75)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def multiply_gram(change, gram):
    rank = len(change)
    return tuple(
        tuple(
            sum(
                change[left][a] * gram[a][b] * change[right][b]
                for a in range(rank)
                for b in range(rank)
            )
            for right in range(rank)
        )
        for left in range(rank)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    r17_candidates = json.loads(gzip.decompress(R17_CANDIDATES.read_bytes()))
    r17_shapes = json.loads(gzip.decompress(R17_SHAPES.read_bytes()))
    r17_diagnostic = json.loads(R17_DIAGNOSTIC.read_text())
    fermigier = json.loads(gzip.decompress(FERMIGIER_REPLAY.read_bytes()))
    truth_document = json.loads(TRUTH.read_text())
    dimension_document = json.loads(DIMENSION_SCAN.read_text())
    r17_truth_indices = {
        record["label"]: int(record["truth_source_index"])
        for record in r17_diagnostic["controls"]
    }
    r17_truth_bases = {
        record["label"]: primitive_span_basis(
            tuple(tuple(map(int, row)) for row in record["embedding_matrix_columns"])
        )
        for record in truth_document["positive_controls"]
    }
    fermigier_truth_record = next(
        record
        for record in truth_document["fermigier_family_controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    fermigier_truth_basis = primitive_span_basis(
        tuple(
            tuple(map(int, row))
            for row in fermigier_truth_record["embedding_matrix_columns"]
        )
    )
    fermigier_truth_key = rational_nullspace(fermigier_truth_basis)
    fermigier_truth_index = next(
        index
        for index, candidate in enumerate(fermigier["candidates"])
        if tuple(map(tuple, candidate["proposal"]["exact_annihilator_rows"]))
        == fermigier_truth_key
    )
    controls = []
    for fibre, shape_fibre in zip(r17_candidates["fibres"], r17_shapes["fibres"]):
        if fibre["label"] != shape_fibre["label"]:
            raise ArithmeticError("R17 candidate and shape ledgers differ")
        controls.append(
            {
                "label": fibre["label"],
                "family": "R17",
                "bases": [
                    tuple(map(tuple, candidate["basis_rows"]))
                    for candidate in fibre["candidates"]
                ],
                "shapes80": [
                    candidate["log_hermite_invariant"]
                    for candidate in shape_fibre["candidates"]
                ],
                "truth_index": r17_truth_indices[fibre["label"]],
            }
        )
    controls.append(
        {
            "label": "ICARM_245_Fermigier_negative_control",
            "family": "Fermigier_rank12",
            "bases": [
                tuple(map(tuple, candidate["proposal"]["primitive_basis_rows"]))
                for candidate in fermigier["candidates"]
            ],
            "shapes80": [
                candidate["primitive_hermite_signature"]["log_hermite_invariant"]
                for candidate in fermigier["candidates"]
            ],
            "truth_index": fermigier_truth_index,
        }
    )

    ledgers_by_pool = {}
    for pool in STABILITY_POOLS:
        ledgers_by_pool[pool] = []
        for control in controls:
            ledger = exact_graph_walk_consensus(
                control["bases"],
                control["shapes80"],
                pool_size=pool,
                shape_gap_threshold=DEFAULT_GAP,
                graph_weight=DEFAULT_GRAPH_WEIGHT,
                timeout=600,
            )
            ledgers_by_pool[pool].append(ledger)
            print(
                f"GRAPHWALKPROGRESS|pool={pool}|label={control['label']}|"
                f"selected={ledger.selected.source_index}|truth={control['truth_index']}",
                flush=True,
            )

    stability = []
    stability_failures = 0
    for gap in STABILITY_GAPS:
        for weight in STABILITY_WEIGHTS:
            selected_count = 0
            selections = []
            for pool in STABILITY_POOLS:
                for control, ledger in zip(controls, ledgers_by_pool[pool]):
                    selected = rescore_graph_walk_consensus(
                        ledger,
                        shape_gap_threshold=gap,
                        graph_weight=weight,
                    )
                    correct = selected == control["truth_index"]
                    selected_count += int(correct)
                    stability_failures += int(not correct)
                    selections.append(
                        {
                            "pool_size": pool,
                            "label": control["label"],
                            "selected_source_index": selected,
                            "is_withheld_truth": correct,
                        }
                    )
            stability.append(
                {
                    "shape_gap_threshold": f"{gap:.17g}",
                    "graph_weight": f"{weight:.17g}",
                    "exact_selection_count": selected_count,
                    "selection_count": len(selections),
                    "selections": selections,
                }
            )

    pinned_gram = tuple(
        tuple(map(int, line.split()))
        for line in PINNED_R17_GRAM.read_text().splitlines()
        if line.strip()
    )
    default_controls = []
    for control, ledger in zip(controls, ledgers_by_pool[DEFAULT_POOL]):
        selected = ledger.selected.source_index
        record = {
            "label": control["label"],
            "family": control["family"],
            "truth_source_index": control["truth_index"],
            "selected_source_index": selected,
            "selected_is_withheld_truth": selected == control["truth_index"],
            "selected_primitive_embedding_matrix_rows": [
                list(row) for row in control["bases"][selected]
            ],
            "default_ledger": ledger.to_record(),
        }
        if control["family"] == "R17":
            change = tuple(
                tuple(map(int, row))
                for row in row_basis_coordinates(
                    control["bases"][selected], r17_truth_bases[control["label"]]
                )
            )
            smith = row_embedding_smith_invariant_factors(change)
            if any(value != 1 for value in smith):
                raise ArithmeticError("selected R17 basis change is not unimodular")
            record.update(
                {
                    "selected_to_published_r17_basis_change": [
                        list(row) for row in change
                    ],
                    "basis_change_smith_invariant_factors": list(smith),
                    "pinned_r17_gram_in_selected_basis": [
                        list(row) for row in multiply_gram(change, pinned_gram)
                    ],
                }
            )
        else:
            record.update(
                {
                    "selected_exact_annihilator_rows": fermigier["candidates"][
                        selected
                    ]["proposal"]["exact_annihilator_rows"],
                    "postselection_generic_subgroup_index_in_primitive_closure": fermigier_truth_record[
                        "generic_subgroup_index_in_primitive_closure"
                    ],
                    "finite_index_blindly_inferred": False,
                }
            )
        default_controls.append(record)

    precision_controls = []
    for control in controls:
        if control["family"] == "R17":
            rank = int(control["label"].rsplit("_", 1)[1])
            module = importlib.import_module(f"elkies_rank{rank}")
            curve = EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS))
            points = tuple(module.POINTS)
        else:
            curve = EllipticCurve(tuple(FERMIGIER_MODEL))
            points = tuple(FERMIGIER_POINTS)
        gram120 = height_gram(curve, points, digits=120)
        signatures120 = primitive_hermite_signatures(
            gram120,
            control["bases"],
            digits=120,
            maximum_vectors=100_000,
            batch_size=64,
            timeout=600,
        )
        ledger120 = exact_graph_walk_consensus(
            control["bases"],
            [item.hermite.log_hermite_invariant for item in signatures120],
            pool_size=DEFAULT_POOL,
            shape_gap_threshold=DEFAULT_GAP,
            graph_weight=DEFAULT_GRAPH_WEIGHT,
            timeout=600,
        )
        precision_controls.append(
            {
                "label": control["label"],
                "selected_source_index_at_120_digits": ledger120.selected.source_index,
                "selected_is_withheld_truth": ledger120.selected.source_index
                == control["truth_index"],
                "shape_gap_at_80_digits": ledgers_by_pool[DEFAULT_POOL][
                    len(precision_controls)
                ].shape_gap,
                "shape_gap_at_120_digits": ledger120.shape_gap,
                "selector_mode_at_120_digits": ledger120.selector_mode,
            }
        )
        print(
            f"GRAPHWALKPRECISION|label={control['label']}|"
            f"selected={ledger120.selected.source_index}|truth={control['truth_index']}",
            flush=True,
        )

    dimension_selected = dimension_document["negative_control"][
        "dimension_selected_by_max_integrality_llr"
    ]
    default_exact = sum(
        record["selected_is_withheld_truth"] for record in default_controls
    )
    precision_exact = sum(
        record["selected_is_withheld_truth"] for record in precision_controls
    )
    dimension_window_pass = abs(dimension_selected - 12) <= 1
    passed = (
        dimension_window_pass
        and default_exact == len(controls)
        and precision_exact == len(controls)
        and stability_failures == 0
    )
    status = (
        "PASS_PHASE0_GRAPH_WALK_CONTROL_CALIBRATION_WITH_DIMENSION_WINDOW"
        if passed
        else "FAIL_PHASE0_GRAPH_WALK_CONTROL_CALIBRATION_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-graph-walk-calibration.v1",
        "status": status,
        "scope": "Phase-0 controls only; no wgxli target is loaded",
        "algorithm": {
            "candidate_dimension_policy": "blind 10..20 dimension scan before fixed-dimension replay",
            "shape_prefilter": f"top {DEFAULT_POOL} primitive Hermite candidates",
            "shape_gap_threshold": f"{DEFAULT_GAP:.17g}",
            "graph": "edge iff exact rational intersection dimension is at least k-1",
            "exact_graph_features": "triangle count and length-four walk count A^4*1",
            "graph_weight": f"{DEFAULT_GRAPH_WEIGHT:.17g}",
            "graph_score": (
                "graph_weight*(triangle rank percentile + length-four-walk rank "
                "percentile) + Hermite rank percentile"
            ),
        },
        "dimension_calibration": {
            "scanned_dimensions": list(range(10, 21)),
            "fermigier_estimated_dimension": dimension_selected,
            "accepted_calibration_window": [12, 13, 14],
            "truth_dimension_in_estimated_plus_or_minus_one_window": dimension_window_pass,
            "exact_graph_recovery_dimension": 12,
        },
        "default_controls": default_controls,
        "default_exact_selection_count": default_exact,
        "precision_120_digit_controls": precision_controls,
        "precision_120_digit_exact_selection_count": precision_exact,
        "stability_box": {
            "pool_sizes": list(STABILITY_POOLS),
            "shape_gap_thresholds": list(STABILITY_GAPS),
            "graph_weights": list(STABILITY_WEIGHTS),
            "parameter_configuration_count": len(stability),
            "control_selection_count": len(stability)
            * len(STABILITY_POOLS)
            * len(controls),
            "failure_count": stability_failures,
            "configurations": stability,
        },
        "recovered_abstract_positive_control_lattice": {
            "identification": "published_R17",
            "rank": 17,
            "pinned_gram_path": str(PINNED_R17_GRAM.relative_to(ROOT)),
            "pinned_gram_sha256": digest(PINNED_R17_GRAM),
            "identity_basis_gram": [list(row) for row in pinned_gram],
        },
        "gate_decision": (
            "PHASE_0 COMPONENT PASS. The corrected blind scan estimates Fermigier "
            "dimension 13, whose plus-or-minus-one window contains 12. Conditional "
            "on the scanned rank-12 ledger, the fixed selector recovers all four R17 subgroups and "
            "the primitive Fermigier rank-12 rational space throughout the declared "
            "stability box and at 120-digit height precision. The actual index-2^11 "
            "Fermigier sublattice and a unique end-to-end dimension are not inferred blindly. By user scope, do not load "
            "wgxli yet."
            if passed
            else "CLOSED. Modify the invariant before any target use."
        ),
        "proof_boundary": (
            "Candidate coordinates, Smith tests, rational intersection ranks, graph "
            "edges, triangle counts, length-four walk counts, selected embedding "
            "matrices, and R17 basis transports are exact within the pinned proposal "
            "bounds. Canonical heights, Hermite prefiltering, rank-percentile scores, "
            "thresholds, and the claim that calibration generalizes are numerical or "
            "heuristic. The graph selector is calibrated conditional on a candidate "
            "dimension; the corrected blind scan returns 13 rather than uniquely 12. "
            "Fermigier finite index 2^11 is a postselection truth audit, not a blind "
            "inference."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                R17_CANDIDATES,
                R17_SHAPES,
                R17_DIAGNOSTIC,
                FERMIGIER_REPLAY,
                TRUTH,
                DIMENSION_SCAN,
                PINNED_R17_GRAM,
                *(ELLIPTIC / "cas" / f"elkies_rank{rank}.py" for rank in range(25, 29)),
                ELLIPTIC / "cas/icarm_curve245.py",
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version()},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice graph-walk calibration artifact is stale")
        print(
            f"GRAPHWALK|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"GRAPHWALK|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
