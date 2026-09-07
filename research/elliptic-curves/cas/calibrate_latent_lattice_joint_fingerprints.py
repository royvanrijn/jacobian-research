#!/usr/bin/env python3
"""Calibrate joint candidate-level fingerprints on all four R17 controls.

No wgxli target is loaded.  Each fibre independently builds the existing
finite-seeded rank-17 proposal ledger, converts every candidate to a compact
basis-independent relation/height/arithmetic/mod-2/mod-3 fingerprint, and
scores it by nearest candidates in the other three fibres.  Published R17
embeddings are revealed only after the joint scores are fixed.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    FiniteQuotientBlock,
    build_relation_complex,
    candidate_relation_fingerprint,
    enumerate_short_vectors,
    height_gram,
    independent_relation_growth_proposals,
    joint_nearest_candidate_scores,
    primitive_span_basis,
    rational_rank,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_joint_fingerprints_v1.json"
LEDGER_OUTPUT = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
HEIGHT_BOUNDS = {
    "rank_at_least_25": 40.0,
    "rank_at_least_26": 43.0,
    "rank_at_least_27": 52.0,
    "rank_at_least_28": 60.0,
}
SEED_EDGES = 3_000
FINGERPRINT_BOUNDS = {
    "quantiles": 16,
    "projective_multiplicities": 16,
    "finite_primes": [2, 3],
    "nearest_neighbour_chunk_size": 16,
}


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--ledger-output", type=Path, default=LEDGER_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    truth_document = json.loads(TRUTH.read_text())
    finite_document = json.loads(FINITE.read_text())
    truth_by_label = {
        record["label"]: record for record in truth_document["positive_controls"]
    }
    finite_by_label = {
        record["label"]: record
        for record in finite_document["controls"]
        if record["label"] in HEIGHT_BOUNDS
    }
    labels = tuple(HEIGHT_BOUNDS)
    states = []
    started = time.monotonic()
    for label in labels:
        stage_started = time.monotonic()
        rank = int(label.rsplit("_", 1)[1])
        module = importlib.import_module(f"elkies_rank{rank}")
        curve = EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS))
        points = tuple(module.POINTS)
        records = enumerate_short_vectors(
            curve,
            points,
            height_bound=HEIGHT_BOUNDS[label],
            digits=80,
            maximum_lines=100_000,
            materialize_points=False,
        )
        vectors = tuple(record.coordinates for record in records)
        heights = tuple(record.canonical_height for record in records)
        arithmetic = tuple(record.arithmetic for record in records)
        complex_ = build_relation_complex(vectors)
        ambient_gram = tuple(
            tuple(value for value in row)
            for row in height_gram(curve, points, digits=80)
        )
        development = tuple(
            finite_block(record)
            for record in finite_by_label[label]["development_blocks"]
        )
        proposals = independent_relation_growth_proposals(
            records,
            complex_,
            dimension=17,
            seed_edges=SEED_EDGES,
            priority_mode="finite",
            seed_strategy="stratified",
            finite_blocks=development,
        )
        fingerprints = tuple(
            candidate_relation_fingerprint(
                vectors,
                heights,
                arithmetic,
                proposal.inlier_indices,
                complex_,
                dimension=17,
                quantiles=FINGERPRINT_BOUNDS["quantiles"],
                projective_multiplicities=FINGERPRINT_BOUNDS[
                    "projective_multiplicities"
                ],
                finite_primes=FINGERPRINT_BOUNDS["finite_primes"],
            )
            for proposal in proposals
        )
        truth_basis = primitive_span_basis(
            tuple(
                tuple(map(int, row))
                for row in truth_by_label[label]["embedding_matrix_columns"]
            )
        )
        overlaps = tuple(
            overlap_dimension(proposal.basis_rows, truth_basis)
            for proposal in proposals
        )
        truth_indices = tuple(
            index for index, overlap in enumerate(overlaps) if overlap == 17
        )
        if len(truth_indices) != 1:
            raise ArithmeticError(
                f"{label} truth proposal multiplicity changed: {truth_indices}"
            )
        states.append(
            {
                "label": label,
                "records": records,
                "complex": complex_,
                "ambient_gram": ambient_gram,
                "proposals": proposals,
                "fingerprints": fingerprints,
                "overlaps": overlaps,
                "truth_index": truth_indices[0],
                "elapsed_seconds": time.monotonic() - stage_started,
            }
        )
        print(
            "LATENTJOINTPROGRESS|"
            f"label={label}|rays={len(records)}|proposals={len(proposals)}|"
            f"truth_index={truth_indices[0]}|"
            f"seconds={time.monotonic() - stage_started:.3f}",
            flush=True,
        )

    scores = joint_nearest_candidate_scores(
        tuple(state["fingerprints"] for state in states),
        chunk_size=FINGERPRINT_BOUNDS["nearest_neighbour_chunk_size"],
    )
    controls = []
    exact_selection_count = 0
    for fibre_index, (state, fibre_scores) in enumerate(zip(states, scores)):
        ranked_mean = sorted(
            fibre_scores,
            key=lambda score: (
                float(score.mean_distance),
                float(score.maximum_distance),
                -score.mutual_neighbour_count,
                score.candidate_index,
            ),
        )
        ranked_mutual = sorted(
            fibre_scores,
            key=lambda score: (
                -score.mutual_neighbour_count,
                float(score.mean_distance),
                float(score.maximum_distance),
                score.candidate_index,
            ),
        )
        truth_index = state["truth_index"]
        truth_mean_rank = next(
            rank
            for rank, score in enumerate(ranked_mean)
            if score.candidate_index == truth_index
        )
        truth_mutual_rank = next(
            rank
            for rank, score in enumerate(ranked_mutual)
            if score.candidate_index == truth_index
        )
        selected = ranked_mean[0]
        selected_overlap = state["overlaps"][selected.candidate_index]
        exact_selection_count += int(selected_overlap == 17)
        controls.append(
            {
                "label": state["label"],
                "ray_count": len(state["records"]),
                "ternary_relation_count": len(
                    state["complex"].ternary_relations
                ),
                "proposal_count": len(state["proposals"]),
                "truth_source_index": truth_index,
                "truth_mean_distance_rank": truth_mean_rank,
                "truth_mutual_first_rank": truth_mutual_rank,
                "truth_score": fibre_scores[truth_index].to_record(),
                "selected_source_index": selected.candidate_index,
                "selected_truth_overlap": selected_overlap,
                "selected_score": selected.to_record(),
                "top_ten_mean_distance": [
                    {
                        "source_index": score.candidate_index,
                        "truth_overlap": state["overlaps"][score.candidate_index],
                        "mean_distance": score.mean_distance,
                        "maximum_distance": score.maximum_distance,
                        "mutual_neighbour_count": score.mutual_neighbour_count,
                    }
                    for score in ranked_mean[:10]
                ],
                "ledger_seconds": f"{state['elapsed_seconds']:.17g}",
            }
        )
    passed = exact_selection_count == len(states)
    status = (
        "PASS_JOINT_R17_FINGERPRINT_SELECTION"
        if passed
        else "FAIL_JOINT_R17_FINGERPRINT_SELECTION_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    ledger_payload = {
        "schema": "elliptic-curves.latent-lattice-joint-fingerprint-ledger.v1",
        "scope": "Blind R17 control candidates; withheld truth overlaps excluded",
        "fingerprint_bounds": FINGERPRINT_BOUNDS,
        "feature_names": list(states[0]["fingerprints"][0].feature_names),
        "fibres": [
            {
                "label": state["label"],
                "ambient_rank": len(state["records"][0].coordinates),
                "ray_count": len(state["records"]),
                "ambient_height_gram": [
                    list(row) for row in state["ambient_gram"]
                ],
                "candidates": [
                    {
                        "source_index": index,
                        "basis_rows": [list(row) for row in proposal.basis_rows],
                        "inlier_indices": list(proposal.inlier_indices),
                        "ray_count": fingerprint.ray_count,
                        "ternary_relation_count": fingerprint.ternary_relation_count,
                        "scaled_relation_count": fingerprint.scaled_relation_count,
                        "integral_ray_count": fingerprint.integral_ray_count,
                        "feature_values": list(fingerprint.feature_values),
                    }
                    for index, (proposal, fingerprint) in enumerate(
                        zip(state["proposals"], state["fingerprints"])
                    )
                ],
            }
            for state in states
        ],
    }
    ledger_rendered = (
        json.dumps(ledger_payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    ledger_bytes = gzip.compress(ledger_rendered, compresslevel=9, mtime=0)
    ledger_sha256 = sha256(ledger_bytes).hexdigest()
    payload = {
        "schema": "elliptic-curves.latent-lattice-joint-fingerprints.v1",
        "status": status,
        "scope": "Phase-0 R17 rank-25--28 controls only; no wgxli target is loaded",
        "algorithm": {
            "proposal_dimension": 17,
            "finite_seed_edges": SEED_EDGES,
            "proposal_strategy": "independent finite-rarity stratified relation growth",
            "fingerprint": (
                "induced unit/scaled relation density and degree profiles; "
                "scale-free ray/relation height profiles; relative point complexity; "
                "mod-2/mod-3 projective-class multiplicity profiles"
            ),
            "standardization": "joint coordinatewise median/MAD with standard-deviation fallback",
            "selection": "minimum mean nearest-candidate distance to each other fibre",
            "fingerprint_bounds": FINGERPRINT_BOUNDS,
        },
        "controls": controls,
        "exact_truth_selection_count": exact_selection_count,
        "candidate_ledger": {
            "path": str(args.ledger_output.relative_to(ROOT)),
            "sha256": ledger_sha256,
            "compressed_bytes": len(ledger_bytes),
            "uncompressed_bytes": len(ledger_rendered),
        },
        "total_elapsed_seconds": f"{time.monotonic() - started:.17g}",
        "gate_decision": (
            "OPEN_FOR_FERMIGIER_CALIBRATION. All four R17 fibres select their exact "
            "withheld truth candidates."
            if passed
            else "CLOSED. Candidate-level joint selection does not choose the exact "
            f"R17 truth in all fibres ({exact_selection_count}/4); modify the "
            "fingerprint before any target use."
        ),
        "proof_boundary": (
            "Proposal coordinates, induced relation counts, projective-class "
            "multiplicities, rational truth overlaps, and candidate identities are "
            "exact within the declared finite clouds and proposal bounds. Canonical "
            "heights, arithmetic profiles, robust standardization, distances, and "
            "nearest-neighbour selection are numerical or heuristic. Withheld R17 "
            "embeddings are postselection diagnostics only."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                TRUTH,
                FINITE,
                *(ELLIPTIC / f"cas/elkies_rank{rank}.py" for rank in range(25, 29)),
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if (
            not args.output.exists()
            or args.output.read_text() != rendered
            or not args.ledger_output.exists()
            or args.ledger_output.read_bytes() != ledger_bytes
        ):
            raise SystemExit("latent-lattice joint-fingerprint artifact is stale")
        print(
            "LATENTJOINT|check=PASS|"
            f"sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.ledger_output.parent.mkdir(parents=True, exist_ok=True)
    args.ledger_output.write_bytes(ledger_bytes)
    args.output.write_text(rendered)
    print(
        f"LATENTJOINT|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
