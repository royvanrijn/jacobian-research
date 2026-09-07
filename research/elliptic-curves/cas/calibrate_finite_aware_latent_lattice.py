#!/usr/bin/env python3
"""Calibrate source-free finite-code invariants on known lattice families.

No wgxli target is loaded by this program.  Finite blocks, short-vector
clouds, and blind proposal ledgers are constructed first.  Exact known
embeddings are revealed only afterward to measure proposal recall and family
cohesion.  Development and held-out reduction blocks are disjoint.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
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
    build_relation_complex,
    candidate_finite_signature,
    enumerate_short_vectors,
    finite_quotient_block,
    finite_signature_distance,
    independent_relation_growth_proposals,
    primitive_span_basis,
    rational_rank,
)


TRUTH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_truth_v1.json"
)
CURVE282_SOURCE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_7fff_zip_public_source_281_282_285_286.json"
)
FERMIGIER_RANK20_SOURCE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_finite_calibration_v1.json"
)
HEIGHT_BOUNDS = {
    "rank_at_least_25": 40.0,
    "rank_at_least_26": 43.0,
    "rank_at_least_27": 52.0,
    "rank_at_least_28": 60.0,
    "ICARM_245_Fermigier_negative_control": 28.0,
    "ICARM_282_Fermigier_sibling": 36.0,
    "Fermigier_u_28917_over_20_sibling": 140.0,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def primes_through(bound: int):
    for value in range(5, bound + 1):
        if all(value % divisor for divisor in range(2, int(value**0.5) + 1)):
            yield value


def curve282() -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...]]:
    document = json.loads(CURVE282_SOURCE.read_text())
    source = next(record for record in document["curves"] if int(record["id"]) == 282)
    model = tuple(Fraction(value) for value in source["ainvs"])
    points = tuple(
        (Fraction(x_value), Fraction(y_value)) for x_value, y_value in source["points"]
    )
    return model, points


def fermigier_rank20() -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...]]:
    document = json.loads(FERMIGIER_RANK20_SOURCE.read_text())
    model = tuple(
        Fraction(value)
        for value in document["models"]["global_minimal"]["coefficients"]
    )
    points = tuple(
        (
            Fraction(record["points"]["global_minimal"]["x"]),
            Fraction(record["points"]["global_minimal"]["y"]),
        )
        for record in document["imported_selected_twenty_basis"]["basis"]
    )
    return model, points


def public_controls() -> dict[str, tuple[tuple[object, ...], tuple[object, ...]]]:
    controls: dict[str, tuple[tuple[object, ...], tuple[object, ...]]] = {
        "ICARM_245_Fermigier_negative_control": (
            tuple(CURVE245_MODEL),
            tuple(CURVE245_POINTS),
        ),
        "ICARM_282_Fermigier_sibling": curve282(),
        "Fermigier_u_28917_over_20_sibling": fermigier_rank20(),
    }
    for rank in range(25, 29):
        module = importlib.import_module(f"elkies_rank{rank}")
        controls[f"rank_at_least_{rank}"] = (
            tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS),
            tuple(module.POINTS),
        )
    return controls


def truth_bases(truth: dict[str, object]) -> dict[str, tuple[tuple[int, ...], ...]]:
    records = list(truth["positive_controls"]) + list(truth["fermigier_family_controls"])
    answer = {}
    for record in records:
        raw = tuple(tuple(map(int, row)) for row in record["embedding_matrix_columns"])
        answer[record["label"]] = primitive_span_basis(raw)
    return answer


def quotient_ensemble(curve, points, *, per_prime: int, prime_bound: int):
    by_relation_prime = {2: [], 3: []}
    for reduction_prime in primes_through(prime_bound):
        for relation_prime in (2, 3):
            if reduction_prime == relation_prime:
                continue
            if len(by_relation_prime[relation_prime]) >= 2 * per_prime:
                continue
            try:
                block = finite_quotient_block(
                    curve, points, reduction_prime, relation_prime
                )
            except (ArithmeticError, ValueError):
                continue
            if block.quotient_dimension == 1:
                by_relation_prime[relation_prime].append(block)
        if all(len(blocks) >= 2 * per_prime for blocks in by_relation_prime.values()):
            break
    if any(len(blocks) < 2 * per_prime for blocks in by_relation_prime.values()):
        raise ArithmeticError("finite quotient scan did not fill both disjoint ensembles")
    development = tuple(
        block
        for relation_prime in (2, 3)
        for block in by_relation_prime[relation_prime][:per_prime]
    )
    held_out = tuple(
        block
        for relation_prime in (2, 3)
        for block in by_relation_prime[relation_prime][per_prime : 2 * per_prime]
    )
    return development, held_out


def overlap_dimension(left, right) -> int:
    return len(left) + len(right) - rational_rank(tuple(left) + tuple(right))


def pairwise_distances(records, key, *, active_only=False):
    answer = {}
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            name = f"{records[left]['label']}/{records[right]['label']}"
            answer[name] = f"{finite_signature_distance(records[left][key], records[right][key], include_components=False, active_prime_blocks_only=active_only, allow_unmatched_blocks=active_only):.17g}"
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--blocks-per-prime", type=int, default=3)
    parser.add_argument("--reduction-prime-bound", type=int, default=251)
    parser.add_argument("--seed-edges", type=int, default=3_000)
    args = parser.parse_args()
    if not 1 <= args.blocks_per_prime <= 4:
        raise SystemExit("--blocks-per-prime must lie in 1..4")

    truth_document = json.loads(TRUTH.read_text())
    withheld_bases = truth_bases(truth_document)
    controls = public_controls()
    computed = []
    rank25_state = None
    for label, (model, points) in controls.items():
        curve = EllipticCurve(model)
        development, held_out = quotient_ensemble(
            curve,
            points,
            per_prime=args.blocks_per_prime,
            prime_bound=args.reduction_prime_bound,
        )
        records = enumerate_short_vectors(
            curve,
            points,
            height_bound=HEIGHT_BOUNDS[label],
            digits=80,
            maximum_lines=100_000,
            materialize_points=False,
        )
        complex_ = build_relation_complex([record.coordinates for record in records])

        # Blind stage: neither the dimension nor the basis below is supplied
        # to the general selector.  This bounded truth-dimension proposal run
        # tests only whether finite seeding improves proposal recall.
        truth_basis = withheld_bases[label]
        proposals = independent_relation_growth_proposals(
            records,
            complex_,
            dimension=len(truth_basis),
            seed_edges=args.seed_edges,
            priority_mode="finite",
            seed_strategy="stratified",
            finite_blocks=development,
        )

        # Evaluation stage: reveal the withheld subspace only now.
        overlaps = [overlap_dimension(proposal.basis_rows, truth_basis) for proposal in proposals]
        development_signature = candidate_finite_signature(
            truth_basis, complex_, finite_blocks=development
        )
        held_out_signature = candidate_finite_signature(
            truth_basis, complex_, finite_blocks=held_out
        )
        family = "R17" if label.startswith("rank_at_least_") else "Fermigier_rank12"
        computed.append(
            {
                "label": label,
                "family": family,
                "ambient_rank": len(points),
                "truth_rank": len(truth_basis),
                "height_bound": HEIGHT_BOUNDS[label],
                "short_vector_lines": len(records),
                "complete_relation_edges": len(complex_.ternary_relations),
                "development_blocks": [block.to_record() for block in development],
                "held_out_blocks": [block.to_record() for block in held_out],
                "finite_seeded_proposal_count": len(proposals),
                "finite_seeded_selected_truth_overlap": overlaps[0] if overlaps else 0,
                "finite_seeded_maximum_truth_overlap": max(overlaps, default=0),
                "finite_seeded_exact_truth_blind_ranks": [
                    rank for rank, overlap in enumerate(overlaps) if overlap == len(truth_basis)
                ],
                "development_truth_signature_object": development_signature,
                "held_out_truth_signature_object": held_out_signature,
            }
        )
        if label == "rank_at_least_25":
            rank25_state = (proposals, complex_, development, held_out, truth_basis)

    r17 = [record for record in computed if record["family"] == "R17"]
    fermigier = [record for record in computed if record["family"] == "Fermigier_rank12"]
    cohesion = {
        "R17_development": pairwise_distances(r17, "development_truth_signature_object"),
        "R17_held_out": pairwise_distances(r17, "held_out_truth_signature_object"),
        "Fermigier_rank12_development": pairwise_distances(
            fermigier, "development_truth_signature_object"
        ),
        "Fermigier_rank12_held_out": pairwise_distances(
            fermigier, "held_out_truth_signature_object"
        ),
        "Fermigier_rank12_development_active_blocks_only": pairwise_distances(
            fermigier, "development_truth_signature_object", active_only=True
        ),
        "Fermigier_rank12_held_out_active_blocks_only": pairwise_distances(
            fermigier, "held_out_truth_signature_object", active_only=True
        ),
    }
    if rank25_state is None:
        raise ArithmeticError("rank-25 selector diagnostic state is missing")
    rank25_proposals, rank25_complex, rank25_development, rank25_held_out, rank25_truth = (
        rank25_state
    )
    r17_development_references = [
        record["development_truth_signature_object"]
        for record in r17
        if record["label"] != "rank_at_least_25"
    ]
    r17_held_out_references = [
        record["held_out_truth_signature_object"]
        for record in r17
        if record["label"] != "rank_at_least_25"
    ]
    rank25_scores = []
    for source_rank, proposal in enumerate(rank25_proposals):
        signature = candidate_finite_signature(
            proposal.basis_rows, rank25_complex, finite_blocks=rank25_development
        )
        score = sum(
            finite_signature_distance(signature, reference, include_components=False)
            for reference in r17_development_references
        ) / len(r17_development_references)
        rank25_scores.append((score, source_rank, proposal, signature))
    rank25_scores.sort(key=lambda item: (item[0], item[1]))
    truth_source_rank = next(
        rank
        for rank, proposal in enumerate(rank25_proposals)
        if overlap_dimension(proposal.basis_rows, rank25_truth) == len(rank25_truth)
    )
    truth_finite_rank = next(
        rank for rank, item in enumerate(rank25_scores) if item[1] == truth_source_rank
    )
    selected_score, selected_source_rank, selected_proposal, _selected_signature = (
        rank25_scores[0]
    )
    selected_held_out = candidate_finite_signature(
        selected_proposal.basis_rows,
        rank25_complex,
        finite_blocks=rank25_held_out,
    )
    truth_held_out = next(
        record["held_out_truth_signature_object"]
        for record in r17
        if record["label"] == "rank_at_least_25"
    )
    rank25_joint_selector_diagnostic = {
        "role": (
            "Leave-rank25-out calibration: rank every rank-25 proposal against the "
            "rank-26--28 development truth profiles, then evaluate on disjoint held-out blocks"
        ),
        "proposal_count": len(rank25_proposals),
        "selected_source_rank": selected_source_rank,
        "selected_development_distance": f"{selected_score:.17g}",
        "selected_withheld_truth_overlap": overlap_dimension(
            selected_proposal.basis_rows, rank25_truth
        ),
        "selected_held_out_distance": f"{sum(finite_signature_distance(selected_held_out, reference, include_components=False) for reference in r17_held_out_references) / len(r17_held_out_references):.17g}",
        "truth_source_rank": truth_source_rank,
        "truth_finite_profile_rank": truth_finite_rank,
        "truth_development_distance": f"{rank25_scores[truth_finite_rank][0]:.17g}",
        "truth_held_out_distance": f"{sum(finite_signature_distance(truth_held_out, reference, include_components=False) for reference in r17_held_out_references) / len(r17_held_out_references):.17g}",
        "top_ten_development_candidates": [
            {
                "finite_profile_rank": rank,
                "source_rank": item[1],
                "distance": f"{item[0]:.17g}",
                "support": item[2].support,
            }
            for rank, item in enumerate(rank25_scores[:10])
        ],
    }
    for record in computed:
        record["development_truth_signature"] = record.pop(
            "development_truth_signature_object"
        ).to_record()
        record["held_out_truth_signature"] = record.pop(
            "held_out_truth_signature_object"
        ).to_record()

    exact_r17 = sum(bool(record["finite_seeded_exact_truth_blind_ranks"]) for record in r17)
    fermigier245 = next(
        record for record in computed if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    passed = exact_r17 == len(r17) and bool(
        fermigier245["finite_seeded_exact_truth_blind_ranks"]
    )
    payload = {
        "schema": "elliptic-curves.latent-lattice-finite-calibration.v1",
        "status": "PASS_FINITE_PROPOSAL_RECALL" if passed else "FAIL_FINITE_PROPOSAL_RECALL",
        "scope": "Phase-0 controls only; no wgxli target curve is loaded",
        "algorithm": {
            "development_blocks": (
                f"first {args.blocks_per_prime} one-dimensional E(F_p)/ell E(F_p) "
                "blocks for each ell in {2,3}"
            ),
            "held_out_blocks": (
                f"next disjoint {args.blocks_per_prime} one-dimensional blocks for each ell"
            ),
            "reduction_prime_bound": args.reduction_prime_bound,
            "finite_seed_edges": args.seed_edges,
            "candidate_signature": (
                "candidate image ranks, unoriented class multiplicities, and induced "
                "unit/scaled relation-code types; raw primes and quotient bases excluded"
            ),
            "matching": "exact permutation matching among equal source-free block types",
        },
        "controls": computed,
        "truth_family_pairwise_profile_distances": cohesion,
        "rank25_joint_selector_diagnostic": rank25_joint_selector_diagnostic,
        "gate_decision": (
            "Finite awareness is not allowed onto unknown curves unless all four R17 "
            "truth spaces and the Fermigier rank-12 truth space occur in the bounded "
            "blind proposal ledgers."
        ),
        "proof_boundary": (
            "Finite quotient maps, candidate image ranks, class histograms, additive "
            "relations, rational overlaps, and the declared bounded enumerations are "
            "exact. Proposal ordering and canonical-height cutoffs are heuristic. "
            "Truth signatures and overlaps are post-selection calibration diagnostics, "
            "not inputs to the blind generator."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                TRUTH,
                CURVE282_SOURCE,
                FERMIGIER_RANK20_SOURCE,
                ELLIPTIC / "cas/icarm_curve245.py",
                *(ELLIPTIC / "cas" / f"elkies_rank{rank}.py" for rank in range(25, 29)),
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
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
        f"FINITECAL|status={payload['status']}|R17={exact_r17}/4|"
        f"F245max={fermigier245['finite_seeded_maximum_truth_overlap']}/12|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
