#!/usr/bin/env python3
"""Build the deterministic manifest for the structural elliptic-search lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from finite_quotient_escape import QuotientBlock, analyze_escape
from structural_search import (
    BranchDivisor,
    ExternalComputationTask,
    IntegralLattice,
    ProjectivePadicBall,
    ProjectiveRational,
    projective_congruence_lattice,
    ResidualSelmerBudget,
    TwistCharacter,
    V4CoverDecomposition,
    enumerate_isotropic_fibration_candidates,
    enumerate_k3_divisor_candidates,
    two_division_cubic,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic_structural_search_groundwork.json"
)
R20_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
FERMIGIER_GENERIC_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_fermigier_generic_rank_exact.json"
)
FERMIGIER_TRANSPORT_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_fermigier_exceptional_transport.json"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(record: dict[str, Any]) -> str:
    encoded = json.dumps(
        record, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _load_r20() -> dict[str, Any]:
    record = json.loads(R20_ARTIFACT.read_text())
    if record["global_arithmetic"]["rank_lower_bound"] != 20:
        raise RuntimeError("the canonical R20 rank lower bound changed")
    if record["identity"]["candidate_key"] != "fermigier-mestre-v1:u=28917/20":
        raise RuntimeError("the canonical R20 identity changed")
    return record


def build_manifest() -> dict[str, Any]:
    r20 = _load_r20()
    ainvs = tuple(
        r20["models"]["global_minimal"]["coefficients"]
    )
    relation_certificate = r20["imported_selected_twenty_basis"][
        "imported_ecsearch_cyclic_log_mod5_certificate"
    ]
    rows = tuple(tuple(row["logs"]) for row in relation_certificate["rows"])
    point_count = len(rows[0])
    if point_count != 20 or any(len(row) != point_count for row in rows):
        raise RuntimeError("the pinned R20 finite-reduction matrix changed shape")
    escape = analyze_escape(
        (
            QuotientBlock.build(
                modulus=5,
                rows=rows,
                column_count=point_count,
                source=(
                    "R20 pinned cyclic-log certificate at good reduction primes"
                ),
            ),
        ),
        known_column_count=12,
        candidate_labels=tuple(f"basis-column-{index}" for index in range(13, 21)),
    )
    if (
        escape.baseline_rank,
        escape.combined_rank,
        escape.marginal_dimension,
    ) != (12, 20, 8):
        raise RuntimeError("the R20 quotient-escape calibration changed")

    first = BranchDivisor.squarefree_polynomial(6, prefix="f1")
    second = BranchDivisor.squarefree_polynomial(6, prefix="f2")
    v4 = V4CoverDecomposition(first, second)
    if v4.quotient_genera != (2, 2, 5) or v4.cover_genus != 9:
        raise RuntimeError("the V4 genus calibration failed")

    # U is only a regression fixture for the lattice enumerator.  The actual
    # Fermigier NS Gram matrix remains an explicit first task below.
    hyperbolic_plane = IntegralLattice(
        gram=((0, 1), (1, 0)), labels=("F", "S")
    )
    lattice_examples = enumerate_k3_divisor_candidates(
        hyperbolic_plane,
        fiber_vector=(1, 0),
        coefficient_bound=2,
        fiber_degrees=(1, 2),
        self_intersections=(-2, 0, 2, 4),
    )
    isotropic_examples = enumerate_isotropic_fibration_candidates(
        hyperbolic_plane,
        ample_vector=(1, 1),
        coefficient_bound=2,
    )

    infinity_ball = ProjectivePadicBall(
        prime=7, exponent=2, chart="infinity", residue=0
    )
    infinity_example = ProjectiveRational.normalized(1, 49)
    if not infinity_ball.matches(infinity_example):
        raise RuntimeError("the projective p-adic infinity fixture failed")
    projective_lattice = projective_congruence_lattice(
        (
            ProjectivePadicBall(prime=5, exponent=1, chart="infinity", residue=0),
            ProjectivePadicBall(prime=7, exponent=1, chart="affine", residue=3),
        )
    )
    projective_example = ProjectiveRational.normalized(1, 5)
    if not projective_lattice.point_matches_charts(projective_example):
        raise RuntimeError("the mixed-chart projective lattice fixture failed")

    twist = TwistCharacter(
        "degree-six-test-character",
        BranchDivisor.squarefree_polynomial(6, prefix="twist"),
    )

    tasks = (
        ExternalComputationTask(
            task_id="R20-RESIDUAL-2SELMER",
            lane="residual Selmer and cubic class group",
            exact_inputs={
                "candidate": r20["identity"]["candidate_key"],
                "minimal_model": list(ainvs),
                "known_rank_lower_bound": 20,
                "two_division_cubic_coefficients_ascending": list(
                    two_division_cubic(ainvs)
                ),
                "known_basis_sha256": r20["imported_selected_twenty_basis"][
                    "canonical_point_sequence_sha256"
                ],
            },
            engine_options=("Magma", "SageMath", "PARI/GP plus class-group tooling"),
            required_outputs=(
                "2-Selmer dimension and rational 2-torsion dimension",
                "2-division cubic field and class-group 2-primary data",
                "Kummer images of all twenty pinned generators",
                "explicit residual locally soluble 2-cover representatives",
                "minimized and reduced cover equations",
            ),
            success_gate=(
                "either certify rank upper bound 20, or emit at least one "
                "residual cover class not explained by the pinned subgroup"
            ),
            claim_boundary=(
                "a residual Selmer class may lie in Sha; only a rational point "
                "on a residual cover creates a new Mordell--Weil direction"
            ),
        ),
        ExternalComputationTask(
            task_id="FERMIGIER-NS-LATTICE",
            lane="Neron--Severi and alternate K3 fibrations",
            exact_inputs={
                "family": "fermigier-mestre-v1",
                "known_arithmetic_NS_rank": 17,
                "trivial_lattice_rank": 5,
                "known_MW_rank": 12,
                "generic_rank_artifact": str(
                    FERMIGIER_GENERIC_ARTIFACT.relative_to(ROOT)
                ),
                "initial_fiber_degrees": [2, 3, 4],
                "initial_self_intersections": [-2, 0, 2, 4, 6],
            },
            engine_options=("SageMath", "Magma", "PARI/GP for auxiliary forms"),
            required_outputs=(
                "saturated Gram matrix for the known NS sublattice",
                "discriminant form and possible finite-index overlattices",
                "bounded primitive divisor classes by fiber degree and square",
                "effectiveness/nefness checks for retained multisections",
                "primitive isotropic classes defining alternate elliptic fibrations",
            ),
            success_gate=(
                "produce a geometrically realized low-genus multisection or an "
                "alternate fibration in which it becomes a section"
            ),
            claim_boundary=(
                "lattice enumeration alone does not prove effectiveness, "
                "irreducibility or existence of a rational curve"
            ),
        ),
        ExternalComputationTask(
            task_id="FERMIGIER-V4-PILOT",
            lane="V4 pair-cover Jacobian decomposition",
            exact_inputs={
                "transport_artifact": str(
                    FERMIGIER_TRANSPORT_ARTIFACT.relative_to(ROOT)
                ),
                "known_pair_cover_count": 3160,
                "calibration_quotient_genera": [2, 2, 5],
                "calibration_cover_genus": 9,
            },
            engine_options=("Magma", "SageMath"),
            required_outputs=(
                "one pinned pair cover with exact quotient equations",
                "Jacobian models for both genus-2 quotients and the genus-5 quotient",
                "Selmer/rank bounds for quotient Jacobians",
                "known-anchor divisor classes in each Jacobian",
                "Mordell--Weil sieve or Chabauty plan with explicit local data",
            ),
            success_gate=(
                "determine all rational points on one pair cover, or find a "
                "third rational point beyond the prescribed anchors"
            ),
            claim_boundary=(
                "the V4 isogeny decomposition reduces the problem but is not a "
                "rational-point theorem"
            ),
        ),
        ExternalComputationTask(
            task_id="TWIST-CHARACTER-ENGINE",
            lane="quadratic anti-invariant rank gain",
            exact_inputs={
                "starting_families": [
                    "fermigier-mestre-v1",
                    "mestre-(0,23,93,128,133,175)",
                    "mestre-(0,25,95,143,168,205)",
                ],
                "initial_branch_count_cap": 8,
                "initial_base_change_genus_cap": 3,
            },
            engine_options=("Singular", "SageMath", "Magma"),
            required_outputs=(
                "normalized squareclass d(T) with full branch divisor",
                "twisted surface model and singular-fiber inventory",
                "Frobenius/Picard upper filter",
                "exact section identity on E^(d) when found",
                "specialization independence test against the old subgroup",
            ),
            success_gate=(
                "certify a non-torsion twist section whose specialization escapes "
                "the old Mordell--Weil image"
            ),
            claim_boundary=(
                "making an old section divisible is saturation, not rank gain"
            ),
        ),
        ExternalComputationTask(
            task_id="MESTRE-MODULI-COMPONENTS",
            lane="algebraic family design in root moduli",
            exact_inputs={
                "constraints": [
                    "Mestre quartic obstruction",
                    "one or more affine companion-section equations",
                    "split-infinity squareclass condition",
                    "bounded discriminant-frontier complexity",
                ],
                "finite_field_primes": [101, 103, 107],
            },
            engine_options=("Singular", "Macaulay2", "SageMath"),
            required_outputs=(
                "primary decompositions modulo several good primes",
                "component dimensions/degrees and singular loci",
                "characteristic-zero lifted candidate components",
                "rational parametrizations or low-genus component models",
                "deduplication modulo affine root normalization and reflection",
            ),
            success_gate=(
                "produce a positive-dimensional component with provable extra "
                "sections and better conductor geometry than current templates"
            ),
            claim_boundary=(
                "a finite-field component is a search lead until lifted and "
                "verified over Q"
            ),
        ),
        ExternalComputationTask(
            task_id="ISOGENY-HOPPING",
            lane="isogeny-class point and descent transfer",
            exact_inputs={
                "candidate_ids": [
                    "fermigier-mestre-v1:u=28917/20",
                    "nagao-section7:T=5081/47",
                    "mestre-dsquare-four:u=197",
                ]
            },
            engine_options=("PARI/GP", "SageMath", "Magma"),
            required_outputs=(
                "complete Q-rational isogeny graph",
                "explicit maps in both directions",
                "minimal models and exact conductors for every vertex",
                "point-search and descent complexity metrics per vertex",
                "transported relations for every newly found point",
            ),
            success_gate=(
                "find a generator or complete a descent on an isogenous model "
                "that was inaccessible on the original model"
            ),
            claim_boundary=(
                "rank and conductor are isogeny invariants; easier coordinates "
                "do not themselves change either"
            ),
        ),
        ExternalComputationTask(
            task_id="PROJECTIVE-PADIC-INFINITY",
            lane="local conductor shaping on P1(Q_p)",
            exact_inputs={
                "families": ["fermigier-mestre-v1"],
                "initial_primes": [5, 7, 11, 13, 17, 19],
                "charts": ["affine T", "infinity S=1/T"],
                "known_infinity_fiber": "split I4 for Fermigier",
            },
            engine_options=("exact Python local tables", "PARI/GP verification"),
            required_outputs=(
                "compressed p-adic balls in both projective charts",
                "minimal-model scaling and Kodaira type on each ball",
                "homogeneous CRT/lattice reconstruction constraints",
                "exact conductor replay after specialization",
            ),
            success_gate=(
                "identify denominator-divisible residue balls with favorable "
                "discriminant valuation per conductor cost"
            ),
            claim_boundary=(
                "a raw valuation at infinity is not conductor data until the "
                "specialized global model is minimized"
            ),
        ),
    )

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "status": "groundwork implemented; structural computations pending",
        "claim_level": (
            "exact reusable arithmetic and deterministic task manifests only; "
            "no new rank, Selmer, Picard, class-group or rational-point theorem"
        ),
        "pinned_inputs": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in (
                R20_ARTIFACT,
                FERMIGIER_GENERIC_ARTIFACT,
                FERMIGIER_TRANSPORT_ARTIFACT,
            )
        },
        "exact_calibrations": {
            "R20_quotient_escape": escape.to_record(),
            "R20_two_division_cubic_coefficients_ascending": list(
                two_division_cubic(ainvs)
            ),
            "V4_pair_cover": v4.to_record(),
            "K3_lattice_enumerator_fixture": {
                "gram": [list(row) for row in hyperbolic_plane.gram],
                "multisection_candidates": [
                    candidate.to_record() for candidate in lattice_examples
                ],
                "isotropic_candidates": [
                    list(vector) for vector in isotropic_examples
                ],
                "fixture_only": True,
            },
            "projective_padic_infinity_fixture": {
                "ball": infinity_ball.to_record(),
                "point": infinity_example.to_record(),
                "matches": True,
                "mixed_chart_congruence_lattice": projective_lattice.to_record(),
                "mixed_chart_example_point": projective_example.to_record(),
            },
            "twist_character_fixture": twist.to_record(),
            "residual_selmer_schema_fixture": ResidualSelmerBudget(
                ell=2,
                selmer_dimension=22,
                known_free_rank=20,
                rational_ell_torsion_dimension=0,
            ).to_record(),
        },
        "priority_order": [
            "R20-RESIDUAL-2SELMER",
            "FERMIGIER-NS-LATTICE",
            "FERMIGIER-V4-PILOT",
            "TWIST-CHARACTER-ENGINE",
            "MESTRE-MODULI-COMPONENTS",
            "ISOGENY-HOPPING",
            "PROJECTIVE-PADIC-INFINITY",
        ],
        "tasks": [task.to_record() for task in tasks],
        "integration_rule": (
            "every point-producing lane must run finite-quotient escape before "
            "expensive height escalation, while every family-producing lane "
            "must run Frobenius/Picard and conductor-geometry filters before "
            "large specialization scans"
        ),
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas python3 "
            "elliptic-curves/cas/build_structural_search_groundwork.py"
        ),
    }
    manifest["result_sha256"] = canonical_digest(manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    manifest = build_manifest()
    serialized = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing pinned manifest: {args.output}")
        if args.output.read_text() != serialized:
            raise SystemExit("structural-search groundwork manifest is stale")
        print("STRUCTURAL_SEARCH_GROUNDWORK_CHECK_PASS")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
    print(f"STRUCTURAL_SEARCH_GROUNDWORK_WRITTEN={args.output}")
    print(f"RESULT_SHA256={manifest['result_sha256']}")


if __name__ == "__main__":
    main()
