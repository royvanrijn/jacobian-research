#!/usr/bin/env python3
"""Score the pinned truth-free ICARM-245 replay and keep the target gate honest."""

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

from latent_lattice import (  # noqa: E402
    candidate_finite_signature_from_record,
    exact_intersection_consensus,
    finite_signature_distance,
    primitive_span_basis,
    rational_nullspace,
    rational_rank,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
REPLAY = ARTIFACTS / "latent_lattice_fermigier_replay_v1.json.gz"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
DIMENSION_SCAN = ARTIFACTS / "latent_lattice_calibration_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_fermigier_consensus_v1.json"
POOL_SIZE = 64


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def overlap_dimension(left, right) -> int:
    return len(left) + len(right) - rational_rank(tuple(left) + tuple(right))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    replay = json.loads(gzip.decompress(REPLAY.read_bytes()))
    truth_document = json.loads(TRUTH.read_text())
    finite_document = json.loads(FINITE.read_text())
    dimension_document = json.loads(DIMENSION_SCAN.read_text())
    bases = [
        tuple(map(tuple, candidate["proposal"]["primitive_basis_rows"]))
        for candidate in replay["candidates"]
    ]
    shapes = [
        candidate["primitive_hermite_signature"]["log_hermite_invariant"]
        for candidate in replay["candidates"]
    ]
    consensus = exact_intersection_consensus(
        bases, shapes, pool_size=POOL_SIZE, timeout=600
    )
    truth_record = next(
        record
        for record in truth_document["fermigier_family_controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    truth_basis = primitive_span_basis(
        tuple(tuple(map(int, row)) for row in truth_record["embedding_matrix_columns"])
    )
    truth_key = rational_nullspace(truth_basis)
    truth_index = next(
        index
        for index, candidate in enumerate(replay["candidates"])
        if tuple(map(tuple, candidate["proposal"]["exact_annihilator_rows"]))
        == truth_key
    )
    ranked = sorted(
        consensus.candidates,
        key=lambda item: (
            -float(item.combined_score),
            -float(item.shape_value),
            item.source_index,
        ),
    )
    truth_consensus_rank = next(
        rank for rank, item in enumerate(ranked) if item.source_index == truth_index
    )
    selected = consensus.selected.source_index
    references = [
        record
        for record in finite_document["controls"]
        if record["family"] == "Fermigier_rank12"
        and record["label"] != "ICARM_245_Fermigier_negative_control"
    ]
    finite_audits = {}
    for split in ("development", "held_out"):
        scores = []
        for candidate in replay["candidates"]:
            signature = candidate_finite_signature_from_record(
                candidate[f"{split}_finite_signature"]
            )
            distances = [
                finite_signature_distance(
                    signature,
                    candidate_finite_signature_from_record(
                        reference[f"{split}_truth_signature"]
                    ),
                    include_components=False,
                    active_prime_blocks_only=True,
                    allow_unmatched_blocks=True,
                )
                for reference in references
            ]
            scores.append(sum(distances) / len(distances))
        order = sorted(range(len(scores)), key=lambda index: (scores[index], index))
        finite_audits[split] = {
            "truth_zero_based_rank": order.index(truth_index),
            "truth_distance": f"{scores[truth_index]:.17g}",
            "selected_zero_based_rank": order.index(selected),
            "selected_distance": f"{scores[selected]:.17g}",
        }
    top = [
        {
            **item.to_record(),
            "withheld_truth_intersection_dimension": overlap_dimension(
                bases[item.source_index], truth_basis
            ),
        }
        for item in ranked[:16]
    ]
    earlier = dimension_document["negative_control"]
    dimension_pass = earlier["dimension_selected_by_max_integrality_llr"] == 12
    exact_selection = selected == truth_index
    status = (
        "PASS_FERMIGIER_EXACT_INTERSECTION_CONSENSUS"
        if dimension_pass and exact_selection
        else "PARTIAL_FERMIGIER_DIMENSION_PASS_SELECTOR_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-fermigier-consensus.v1",
        "status": status,
        "scope": "Phase-0 ICARM-245 control only; no wgxli target is loaded",
        "dimension_calibration": {
            "scan_dimensions": list(range(10, 21)),
            "selected_dimension": earlier["dimension_selected_by_max_integrality_llr"],
            "selected_is_true_fermigier_dimension": dimension_pass,
            "forced_rank17_baseline_truth_intersection_dimension": next(
                record["withheld_truth_intersection_dimension"]
                for record in earlier["dimension_scan"]
                if record["dimension"] == 17
            ),
        },
        "candidate_calibration": {
            "candidate_count": len(bases),
            "truth_source_index": truth_index,
            "truth_original_proposal_zero_based_rank": truth_index,
            "truth_consensus_zero_based_rank": truth_consensus_rank,
            "selected_source_index": selected,
            "selected_truth_intersection_dimension": overlap_dimension(
                bases[selected], truth_basis
            ),
            "selected_primitive_embedding_matrix_rows": [list(row) for row in bases[selected]],
            "top_sixteen": top,
            "consensus": consensus.to_record(),
        },
        "finite_code_audit": {
            "policy": (
                "development and held-out source-free sibling distances are reported "
                "but excluded from the selector because calibration is adverse"
            ),
            **finite_audits,
        },
        "gate_decision": (
            "OPEN. Rank 12 and the exact Fermigier subgroup are selected."
            if status.startswith("PASS_")
            else "CLOSED. The dimension scan selects 12 and exact consensus improves "
            f"the truth from zero-based rank {truth_index} to {truth_consensus_rank}, "
            "but it does not select the exact Fermigier subgroup. Do not apply this "
            "selector to wgxli."
        ),
        "proof_boundary": (
            "Dimension labels, candidate bases, rational intersections, finite-code "
            "restrictions, truth subspace intersections, and embedding matrices are "
            "exact within the pinned ledgers. Heights and Hermite ordering are numerical. "
            "Finite sibling distances and the combined selector are heuristic."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                REPLAY,
                TRUTH,
                FINITE,
                DIMENSION_SCAN,
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version()},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice Fermigier consensus artifact is stale")
        print(
            f"FERMIGIERCONSENSUS|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"FERMIGIERCONSENSUS|status={status}|truth_rank={truth_consensus_rank}|"
        f"selected={selected}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
