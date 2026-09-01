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

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    build_relation_complex,
    enumerate_short_vectors,
    exact_span_mask,
    rational_rank,
    row_basis_coordinates,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
HEIGHT_BOUNDS = {25: 40.0, 26: 43.0, 27: 52.0, 28: 60.0}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


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
                "truth_subgroup_lines": len(coordinates),
                "truth_subgroup_line_rank": rational_rank(coordinates),
                "truth_subgroup_ternary_relations": len(
                    complex_.ternary_relations
                ),
                "truth_subgroup_relation_digest": complex_.canonical_digest,
            }
        )

    counts = Counter(
        vector for population in coefficient_sets.values() for vector in population
    )
    support_records = []
    for support in (4, 3, 2):
        vectors = tuple(sorted(vector for vector, count in counts.items() if count >= support))
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
            vector for vector, count in training_counts.items() if count >= 2
        }
        held_visible = training_core & coefficient_sets[held_label]
        leave_one_out.append(
            {
                "held_out_label": held_label,
                "training_two_of_three_ray_count": len(training_core),
                "training_two_of_three_rank": rational_rank(tuple(training_core)),
                "held_out_visible_ray_count": len(held_visible),
                "held_out_visible_rank": rational_rank(tuple(held_visible)),
                "held_out_visible_relation_digest": build_relation_complex(
                    tuple(held_visible)
                ).canonical_digest,
            }
        )

    passed = all(record["held_out_visible_rank"] == 17 for record in leave_one_out)
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-relation-consensus.v1",
        "status": (
            "PASS_R17_EXACT_RELATION_SIGNAL" if passed else "FAIL_R17_RELATION_SIGNAL"
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
        "gate_decision": (
            "CLOSED. Full-rank exact held-out signal exists in all four controls, "
            "but a basis-free matcher has not yet recovered the alignments blindly."
        ),
        "proof_boundary": (
            "Integer coordinate conversion, primitivity checks, coefficient-set "
            "intersections, ranks, and additive relations are exact. Canonical-height "
            "cutoffs are numerical. The corresponding generic bases are supplied by "
            "the withheld positive-control artifact, not discovered by this audit."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                TRUTH,
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
