#!/usr/bin/env python3
"""Calibrate exact proposal-intersection consensus on R17 controls only."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC))

from latent_lattice import exact_intersection_consensus  # noqa: E402
from latent_lattice import (  # noqa: E402
    primitive_span_basis,
    row_basis_coordinates,
    row_embedding_smith_invariant_factors,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
CANDIDATES = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
SHAPES = ARTIFACTS / "latent_lattice_joint_shape_ledger_v1.json.gz"
TRUTH_DIAGNOSTIC = ARTIFACTS / "latent_lattice_joint_fingerprints_v1.json"
CALIBRATION_TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
PINNED_R17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ARTIFACTS / "latent_lattice_intersection_consensus_v1.json"
POOL_SIZE = 64


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    candidates = json.loads(gzip.decompress(CANDIDATES.read_bytes()))
    shapes = json.loads(gzip.decompress(SHAPES.read_bytes()))
    truth = json.loads(TRUTH_DIAGNOSTIC.read_text())
    calibration_truth = json.loads(CALIBRATION_TRUTH.read_text())
    truth_bases = {
        record["label"]: primitive_span_basis(
            tuple(tuple(map(int, row)) for row in record["embedding_matrix_columns"])
        )
        for record in calibration_truth["positive_controls"]
    }
    pinned_gram = tuple(
        tuple(map(int, line.split()))
        for line in PINNED_R17_GRAM.read_text().splitlines()
        if line.strip()
    )
    truth_indices = {
        record["label"]: int(record["truth_source_index"])
        for record in truth["controls"]
    }
    controls = []
    exact_count = 0
    for fibre, shape_fibre in zip(candidates["fibres"], shapes["fibres"]):
        if fibre["label"] != shape_fibre["label"]:
            raise ArithmeticError("candidate and shape ledger fibre orders differ")
        ledger = exact_intersection_consensus(
            [candidate["basis_rows"] for candidate in fibre["candidates"]],
            [
                candidate["log_hermite_invariant"]
                for candidate in shape_fibre["candidates"]
            ],
            pool_size=POOL_SIZE,
            timeout=600,
        )
        ranked = sorted(
            ledger.candidates,
            key=lambda item: (
                -float(item.combined_score),
                -float(item.shape_value),
                item.source_index,
            ),
        )
        truth_index = truth_indices[fibre["label"]]
        truth_rank = next(
            (rank for rank, item in enumerate(ranked) if item.source_index == truth_index),
            None,
        )
        selected = ledger.selected
        selected_rows = tuple(
            map(tuple, fibre["candidates"][selected.source_index]["basis_rows"])
        )
        change = tuple(
            tuple(map(int, row))
            for row in row_basis_coordinates(
                selected_rows, truth_bases[fibre["label"]]
            )
        )
        change_smith = row_embedding_smith_invariant_factors(change)
        if any(value != 1 for value in change_smith):
            raise ArithmeticError("selected R17 basis change is not unimodular")
        selected_generic_gram = tuple(
            tuple(
                sum(
                    change[left][a] * pinned_gram[a][b] * change[right][b]
                    for a in range(17)
                    for b in range(17)
                )
                for right in range(17)
            )
            for left in range(17)
        )
        exact_count += int(selected.source_index == truth_index)
        controls.append(
            {
                "label": fibre["label"],
                "truth_source_index": truth_index,
                "truth_consensus_rank": truth_rank,
                "selected_source_index": selected.source_index,
                "selected_is_withheld_truth": selected.source_index == truth_index,
                "selected_primitive_embedding_matrix_rows": fibre["candidates"][
                    selected.source_index
                ]["basis_rows"],
                "selected_to_published_r17_basis_change": [
                    list(row) for row in change
                ],
                "basis_change_smith_invariant_factors": list(change_smith),
                "pinned_r17_gram_in_selected_basis": [
                    list(row) for row in selected_generic_gram
                ],
                "consensus": ledger.to_record(),
            }
        )
        print(
            f"LATENTCONSENSUSPROGRESS|label={fibre['label']}|"
            f"selected={selected.source_index}|truth={truth_index}|"
            f"truth_rank={truth_rank}",
            flush=True,
        )
    passed = exact_count == len(controls)
    status = (
        "PASS_R17_EXACT_INTERSECTION_CONSENSUS"
        if passed
        else "FAIL_R17_EXACT_INTERSECTION_CONSENSUS_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-intersection-consensus.v1",
        "status": status,
        "scope": "Phase-0 R17 rank-25--28 controls only; no wgxli target is loaded",
        "algorithm": {
            "shape_prefilter": f"top {POOL_SIZE} primitive Hermite invariants",
            "exact_relation": "rational intersection dimension at least k-1",
            "score": (
                "equal-weight sum of within-pool normalized Hermite extremality "
                "and exact codimension-one intersection fraction"
            ),
            "tie_break": "Hermite invariant, then source index",
        },
        "controls": controls,
        "exact_withheld_selection_count": exact_count,
        "recovered_abstract_lattice": {
            "identification": "published_R17 positive-control lattice",
            "rank": 17,
            "pinned_gram_path": str(PINNED_R17_GRAM.relative_to(ROOT)),
            "pinned_gram_sha256": digest(PINNED_R17_GRAM),
            "identity_basis_gram": [list(row) for row in pinned_gram],
        },
        "proof_boundary": (
            "Every pooled intersection rank and every reported embedding matrix is "
            "exact. Primitive status is inherited from the separately pinned Smith "
            "audit. Hermite values and the top-64 prefilter are numerical at 80-digit "
            "input precision. Withheld R17 indices are postselection diagnostics."
        ),
        "gate_decision": (
            "OPEN_FOR_FERMIGIER_CALIBRATION. All four R17 controls select the exact "
            "withheld primitive subgroup."
            if passed
            else "CLOSED. Modify the invariant before Fermigier or target use."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CANDIDATES,
                SHAPES,
                TRUTH_DIAGNOSTIC,
                CALIBRATION_TRUTH,
                PINNED_R17_GRAM,
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version()},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice intersection-consensus artifact is stale")
        print(
            f"LATENTCONSENSUS|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTCONSENSUS|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
