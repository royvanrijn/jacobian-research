#!/usr/bin/env python3
"""Calibrate a local cross-dimension persistence test on ICARM 245."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from icarm_curve245 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    build_relation_complex,
    enumerate_short_vectors,
    exact_row_space_intersection,
    exact_span_mask,
    independent_row_basis,
    primitive_span_basis,
    rational_nullspace,
)
from latent_lattice.subspace import integrality_llr  # noqa: E402


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
GRAPH_CALIBRATION = ARTIFACTS / "latent_lattice_graph_walk_calibration_v1.json"
ACTIVE_SCAN = ARTIFACTS / "latent_lattice_calibration_v2.json"
FERMIGIER_REPLAY = ARTIFACTS / "latent_lattice_fermigier_replay_v1.json.gz"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_dimension_persistence_v1.json"
HEIGHT_BOUNDS = (28.0, 29.0)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    graph_document = json.loads(GRAPH_CALIBRATION.read_text())
    scan_document = json.loads(ACTIVE_SCAN.read_text())
    replay = json.loads(gzip.decompress(FERMIGIER_REPLAY.read_bytes()))
    selected_record = next(
        record
        for record in graph_document["default_controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    selected_core = tuple(
        tuple(map(int, row))
        for row in selected_record["selected_primitive_embedding_matrix_rows"]
    )
    if len(selected_core) != 12:
        raise ArithmeticError("graph-selected Fermigier core is no longer rank 12")

    active_bases = tuple(
        (
            int(record["dimension"]),
            tuple(tuple(map(int, row)) for row in record["primitive_basis_rows"]),
        )
        for record in scan_document["negative_control"]["dimension_scan"]
    )
    face_keys = {}
    for candidate in replay["candidates"]:
        other = tuple(
            tuple(map(int, row))
            for row in candidate["proposal"]["primitive_basis_rows"]
        )
        face = exact_row_space_intersection(selected_core, other)
        if len(face) != 11:
            continue
        key = rational_nullspace(face)
        face_keys.setdefault(key, face)
    face_bases = tuple(basis for _key, basis in sorted(face_keys.items()))

    curve = EllipticCurve(tuple(GENERAL_WEIERSTRASS_COEFFICIENTS))
    evaluations = []
    for height_bound in HEIGHT_BOUNDS:
        records = enumerate_short_vectors(
            curve,
            POINTS,
            height_bound=height_bound,
            digits=80,
            maximum_lines=100_000,
            materialize_points=True,
        )
        coordinates = tuple(record.coordinates for record in records)
        complex_ = build_relation_complex(coordinates)
        integral = np.asarray(
            [bool(record.arithmetic["integral"]) for record in records], dtype=bool
        )
        complex_index = {
            vertex: index for index, vertex in enumerate(complex_.vertices)
        }
        edges = np.asarray(complex_.ternary_relations, dtype=np.int64)

        def candidate_stats(basis_rows) -> dict[str, object]:
            basis = tuple(tuple(map(int, row)) for row in basis_rows)
            mask = exact_span_mask(coordinates, basis)
            relation_mask = np.zeros(len(complex_.vertices), dtype=bool)
            for input_index in np.flatnonzero(mask):
                relation_mask[complex_index[coordinates[input_index]]] = True
            relation_count = (
                int(np.sum(np.all(relation_mask[edges], axis=1)))
                if len(edges)
                else 0
            )
            return {
                "dimension": len(basis),
                "rational_basis_rows": [list(row) for row in basis],
                "support": int(np.sum(mask)),
                "integral_support": int(np.sum(mask & integral)),
                "integrality_llr": f"{integrality_llr(mask, integral):.17g}",
                "induced_ternary_relation_count": relation_count,
            }

        active_candidates = [
            {
                "source": "active_v2_dimension_winner",
                **candidate_stats(basis),
            }
            for _dimension, basis in active_bases
        ]
        core_stats = {
            "source": "exact_graph_walk_rank12_core",
            **candidate_stats(selected_core),
        }
        face_candidates = [
            {"source": "replay_graph_rank11_face", **candidate_stats(basis)}
            for basis in face_bases
        ]
        extension_keys = {}
        for vector in coordinates:
            basis = independent_row_basis(selected_core + (vector,))
            if len(basis) != 13:
                continue
            key = rational_nullspace(basis)
            extension_keys.setdefault(key, basis)
        extension_candidates = [
            {"source": "retained_short_ray_extension", **candidate_stats(basis)}
            for _key, basis in sorted(extension_keys.items())
        ]
        all_candidates = (
            active_candidates + [core_stats] + face_candidates + extension_candidates
        )
        ranked = sorted(
            all_candidates,
            key=lambda record: (
                -float(record["integrality_llr"]),
                -int(record["induced_ternary_relation_count"]),
                record["source"],
            ),
        )
        by_dimension = []
        best_by_dimension = {}
        for dimension in range(10, 21):
            candidates = [
                record
                for record in all_candidates
                if int(record["dimension"]) == dimension
            ]
            best = min(
                candidates,
                key=lambda record: (
                    -float(record["integrality_llr"]),
                    -int(record["induced_ternary_relation_count"]),
                    record["source"],
                ),
            )
            best_by_dimension[dimension] = best
            by_dimension.append(
                {
                    "dimension": dimension,
                    "candidate_count": len(candidates),
                    "best": best,
                }
            )
        persistence = {
            dimension: min(
                float(best_by_dimension[neighbour]["integrality_llr"])
                for neighbour in (dimension - 1, dimension, dimension + 1)
                if neighbour in best_by_dimension
            )
            for dimension in best_by_dimension
        }
        selected_dimension = max(
            persistence,
            key=lambda dimension: (
                persistence[dimension],
                float(best_by_dimension[dimension]["integrality_llr"]),
                -dimension,
            ),
        )
        selected = best_by_dimension[selected_dimension]
        evaluations.append(
            {
                "height_bound": height_bound,
                "short_vector_count": len(records),
                "rank11_face_count": len(face_candidates),
                "rank13_extension_count": len(extension_candidates),
                "selected_dimension": selected_dimension,
                "selected": selected,
                "persistence_scores": [
                    {
                        "dimension": dimension,
                        "three_level_bottleneck_integrality_llr": f"{persistence[dimension]:.17g}",
                    }
                    for dimension in sorted(persistence)
                ],
                "dimension_ledger": by_dimension,
                "top_twenty_by_raw_integrality_llr": ranked[:20],
            }
        )

    truth_document = json.loads(TRUTH.read_text())
    truth_record = next(
        record
        for record in truth_document["fermigier_family_controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    truth_basis = primitive_span_basis(
        tuple(tuple(map(int, row)) for row in truth_record["embedding_matrix_columns"])
    )
    truth_key = rational_nullspace(truth_basis)
    for evaluation in evaluations:
        evaluation["selected_is_withheld_truth_rational_space"] = (
            rational_nullspace(evaluation["selected"]["rational_basis_rows"])
            == truth_key
        )
    selected_is_truth_space = all(
        evaluation["selected_is_withheld_truth_rational_space"]
        for evaluation in evaluations
    )
    status = (
        "PASS_FERMIGIER_CROSS_DIMENSION_PERSISTENCE"
        if all(evaluation["selected_dimension"] == 12 for evaluation in evaluations)
        and selected_is_truth_space
        else "FAIL_FERMIGIER_CROSS_DIMENSION_PERSISTENCE_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-dimension-persistence.v1",
        "status": status,
        "scope": "Phase-0 ICARM-245 control only; no wgxli target is loaded",
        "bounds": {
            "height_bounds": list(HEIGHT_BOUNDS),
            "active_scan_dimensions": list(range(10, 21)),
            "rank11_faces": "all distinct intersections of the selected core with 128 replay candidates",
            "rank13_extensions": "all distinct spans of the selected core with one retained short ray",
        },
        "all_bounds_select_withheld_truth_rational_space": selected_is_truth_space,
        "active_scan_estimated_dimension_before_persistence": scan_document[
            "negative_control"
        ]["dimension_selected_by_max_integrality_llr"],
        "evaluations": evaluations,
        "gate_decision": (
            "PASS. Replacing the active rank-12 winner by the truth-free graph-walk "
            "core, auditing every retained-ray rank-13 extension, and maximizing the "
            "three-level bottleneck integrality enrichment changes the dimension "
            "estimate from 13 to the exact rank-12 Fermigier rational space at both "
            "height bounds."
            if status.startswith("PASS_")
            else "CLOSED. Cross-dimension persistence does not uniquely recover rank 12."
        ),
        "proof_boundary": (
            "Candidate coordinates, rational faces/extensions, "
            "cloud membership, integrality counts, and additive relation counts are "
            "exact within the declared height-28/29 clouds and replay ledger. Integrality "
            "LLR and its three-level bottleneck use as a dimension score are statistical. Rank-13 exhaustion "
            "covers extensions generated by one retained ray, not all rational "
            "rank-13 superlattices. The selected rank-12 embedding is primitive by the "
            "separate Smith audit; auxiliary face/extension bases are rational-space "
            "representatives rather than saturated embedding certificates. Withheld "
            "truth is used only after selection."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                GRAPH_CALIBRATION,
                ACTIVE_SCAN,
                FERMIGIER_REPLAY,
                TRUTH,
                ELLIPTIC / "cas/icarm_curve245.py",
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice dimension-persistence artifact is stale")
        print(
            f"DIMENSIONPERSISTENCE|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"DIMENSIONPERSISTENCE|status={status}|selected_dimensions="
        f"{','.join(str(item['selected_dimension']) for item in evaluations)}|"
        f"faces={len(face_bases)}|extensions="
        f"{','.join(str(item['rank13_extension_count']) for item in evaluations)}|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
