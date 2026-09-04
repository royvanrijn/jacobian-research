#!/usr/bin/env python3
"""Freeze the sparse quotient-mask rank-32 protocol for ICARM curve 385."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
BLIND = ART / "curve385_iterated_half_lattice_blind_v1.json"
VERIFICATION = ART / "curve385_iterated_half_lattice_verification_v1.json"
PROFILE = ART / "curve385_quotient_weight_profile_v1.json"
POLICY = CAS / "curve385_sparse_quotient_policy.py"
RUNNER = CAS / "run_curve385_sparse_quotient_rank32_search.sage"
LEGACY_RUNNER = CAS / "run_curve385_iterated_half_lattice_search.sage"
OUTPUT = ART / "curve385_sparse_quotient_rank32_protocol_v1.json"

EXPECTED_BLIND_SHA256 = "356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502"
EXPECTED_VERIFICATION_SHA256 = "b281556f5d08250f67b69b2c62a640ac17ba4d03325e4402e85c7d60882c3ae5"
EXPECTED_PROFILE_SHA256 = "c321d1b40d9e5fc77ebff64e5d6584feeab5f503b13eadda4f6d524d0e38162a"

import sys

sys.path.insert(0, str(CAS))
from curve385_sparse_quotient_policy import (  # noqa: E402
    POLICY_DOMAIN,
    STAGE_SPECS,
    canonical_hash,
    stage_plan,
    validate_stage_plan,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def build_payload() -> dict[str, Any]:
    expected_hashes = {
        BLIND: EXPECTED_BLIND_SHA256,
        VERIFICATION: EXPECTED_VERIFICATION_SHA256,
        PROFILE: EXPECTED_PROFILE_SHA256,
    }
    for path, expected in expected_hashes.items():
        if digest(path) != expected:
            raise ArithmeticError(f"frozen input changed: {relative(path)}")
    blind = json.loads(BLIND.read_text())
    verification = json.loads(VERIFICATION.read_text())
    profile = json.loads(PROFILE.read_text())
    if blind.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("the blind M29 search is not frozen")
    if verification.get("status") != "PASS_BLIND_M29_EQUALS_DISPLAYED_PUBLIC_M29":
        raise ArithmeticError("the exact M29 public comparison is not closed")
    if profile.get("status") != "PASS_EXACT_POSTHOC_WEIGHT_PROFILE":
        raise ArithmeticError("the quotient-weight evidence gate is not closed")
    cumulative = profile["cylinder"]["cumulative_weight_profile"]
    if [row["quotient_rank_over_M20"] for row in cumulative] != [7, 9, 9]:
        raise ArithmeticError("the sparse quotient-weight evidence changed")

    plans = {str(bits): stage_plan(bits, 43) for bits in (12, 13, 14)}
    for bit_count, plan in plans.items():
        validate_stage_plan(plan, int(bit_count), 43)

    stage_rule = [
        {
            "index": index,
            "id": stage_id,
            "basis_id": basis_id,
            "basis_weight_shell": list(weights),
        }
        for index, (stage_id, basis_id, weights) in enumerate(STAGE_SPECS, start=1)
    ]
    definition = {
        "objective": (
            "find and exactly certify at least three directions beyond the frozen "
            "curve-385 M29 subgroup, yielding rank E(Q) at least 32"
        ),
        "starting_state": {
            "curve": "curve385-rank29",
            "generic_rank": 17,
            "discovered_rank": 29,
            "quotient_bit_count": 12,
            "source": relative(BLIND),
            "public_fixture_needed_for_selection": False,
        },
        "evidence_gate": {
            "source": relative(PROFILE),
            "curve385_M20_weight_one_quotient_rank": 7,
            "curve385_M20_weight_at_most_two_quotient_rank": 9,
            "curve385_M20_weight_three_marginal_quotient_rank": 0,
            "posthoc_not_success_probability": True,
            "basis_dependence_explicitly_audited": True,
        },
        "search_space": {
            "distinguished_generic_class_count": 43,
            "classes": "the frozen generic-deepest 43 masks",
            "full_12_bit_nonzero_word_chart_count": 176_085,
            "full_nonzero_enumeration_is_not_an_automatic_stage": True,
            "deduplication": (
                "an exact base-point chart already searched in the frozen source or an "
                "earlier sparse stage is recorded and not rerun"
            ),
        },
        "quotient_basis": {
            "natural_basis": (
                "greedy standard-unit complement to the exact generic M17 rows in "
                "the current discovered-group basis"
            ),
            "alternate_basis_domain": POLICY_DOMAIN,
            "alternate_bases": ["alternate-a", "alternate-b"],
            "alternate_basis_construction": (
                "SHA-256 candidates restricted to the quotient bit width, accepted "
                "greedily only on exact GF(2) rank gain"
            ),
        },
        "stage_rule": stage_rule,
        "restart_rule": (
            "after any exact rank gain or finite-index enlargement, finish classification, "
            "recompute the height lattice and quotient complement, and restart at stage 1; "
            "stop successfully as soon as certified rank is at least 32"
        ),
        "primary_campaign": {
            "maximum_stage_each_lattice_state": 2,
            "meaning": "natural quotient weight one, then natural quotient weight two",
            "alternate_bases_require_explicit_stage-limit escalation": True,
        },
        "point_search_budget": {
            "reduced_coordinate_height_bound_each_quartic": 100_000,
            "wall_timeout_seconds_each_quartic": 15,
            "gp_stack_bytes_each_quartic": 1_000_000_000,
            "retries": 0,
            "canonical_height_precision_decimal_digits": 110,
            "operative_CVP_rounding_scale": 1_000_000,
            "audit_CVP_rounding_scale": 100_000,
            "checkpoint_every_completed_searches": 10,
            "maximum_lattice_states": 4,
        },
        "acceptance": {
            "point_identity": "exact rational group-law transport to curve 385",
            "independence": "full finite-reduction certificate through prime 1000",
            "dependence": "exact integral relation in the current discovered-group basis",
            "finite_index": "exact primitive relation and Smith unimodular completion",
            "rank32": "a rank-at-least-32 finite-reduction certificate",
        },
        "fail_closed": {
            "timeout_or_PARl_failure": "stop the stage incomplete with no retry",
            "unclassified_point": "stop UNKNOWN",
            "completed_sparse_miss": "bounded negative result only",
            "stage_limit": "not a rank upper bound",
            "full_Mordell_Weil_saturation": "not claimed",
        },
    }
    definition_hash = canonical_hash(definition)
    input_paths = (
        BLIND,
        VERIFICATION,
        PROFILE,
        POLICY,
        RUNNER,
        LEGACY_RUNNER,
        Path(__file__),
    )
    claim_boundary = [
        "The completed 301-chart profile is an exact posthoc structural result and motivates this prospective ordering only.",
        "Quotient Hamming weight is basis-dependent; two deterministic alternate bases precede the natural weight-three stage.",
        "Every sparse stage is a bounded point search. Its failure to enlarge the discovered group is not a rank upper bound.",
        "Only exact rational point identities and full finite-reduction independence certificates can raise the rank lower bound.",
        "The protocol does not claim that the frozen M29 subgroup is saturated in E(Q).",
    ]
    return {
        "schema": "elliptic-curves.curve385-sparse-quotient-rank32-protocol.v1",
        "status": "FROZEN_BEFORE_SPARSE_RANK32_POINT_SEARCH",
        "protocol_definition": definition,
        "protocol_definition_hash": definition_hash,
        "starting_rank": 29,
        "target_rank": 32,
        "old_class_count": 43,
        "starting_basis_sha256": blind["iterations"][0]["basis_after_sha256"],
        "old_deep43_masks_sha256": blind["old_deep43"]["masks_sha256"],
        "stage_plans_by_quotient_bit_count": plans,
        "input_hashes": {relative(path): digest(path) for path in input_paths},
        "claim_boundary": claim_boundary,
        "generator": relative(Path(__file__)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    encoded = canonical_bytes(build_payload())
    if args.check:
        if not args.output.exists() or args.output.read_bytes() != encoded:
            raise SystemExit(f"stale or missing protocol artifact: {args.output}")
        payload = json.loads(encoded)
        print(
            f"C385SPARSEPROTOCOL|status=PASS|definition="
            f"{payload['protocol_definition_hash']}|sha256={sha256(encoded).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(encoded)
    payload = json.loads(encoded)
    print(
        f"C385SPARSEPROTOCOL|status=WROTE|definition="
        f"{payload['protocol_definition_hash']}|sha256={sha256(encoded).hexdigest()}"
    )


if __name__ == "__main__":
    main()
