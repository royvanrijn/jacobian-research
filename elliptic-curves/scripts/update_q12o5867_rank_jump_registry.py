#!/usr/bin/env python3
"""Admit an exact q12o5867 rank jump to its deduplicated registry."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_rank_jump_registry import (  # noqa: E402
    canonical_global_minimal_q_isomorphism_key,
    canonical_json_sha256,
    documented_q_isomorphism_keys,
    exact_admission_gate,
    merge_exact_entry,
    sha256_file,
    validate_registry,
)
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    verify_finite_quotient_certificate,
)


Q = Fraction
DEFAULT_REGISTRY = (
    REPOSITORY
    / "artifacts/generated-results/elliptic-curves/"
    "q12o5867_certified_rank_jump_registry_v1.json"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--specialization", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--certificate-output", type=Path, required=True)
    parser.add_argument("--reproducing-command", required=True)
    parser.add_argument("--overwrite-certificate", action="store_true")
    args = parser.parse_args()
    sys.set_int_max_str_digits(0)

    specialization = json.loads(args.specialization.read_text())
    probe = json.loads(args.probe.read_text())
    registry = json.loads(args.registry.read_text())
    validate_registry(registry)
    if specialization.get("status") != "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND":
        raise SystemExit("specialization lacks its exact rank-17 baseline certificate")
    probe_status = probe.get("status")
    if probe_status == "PASS_BOUNDED_MWRANK_PROBE" and probe.get("worker_result"):
        worker = probe["worker_result"]
        probe_engine = "mwrank"
    elif probe_status == "PASS_BOUNDED_RATPOINTS_PROBE":
        worker = {
            "finite_quotient_escape": probe["finite_quotient_escape"],
            "exact_nonbaseline_minimal_points": probe[
                "novel_points_up_to_sign_on_minimal_model"
            ],
        }
        probe_engine = "ratpoints"
    elif probe_status == "PASS_BOUNDED_SECTION_NORMALIZED_RATPOINTS_PROBE":
        worker = {
            "finite_quotient_escape": probe["finite_quotient_escape"],
            "exact_nonbaseline_minimal_points": probe[
                "novel_points_up_to_sign_on_minimal_model"
            ],
        }
        probe_engine = "section-normalized-ratpoints"
    else:
        raise SystemExit("probe did not complete with exact candidate data")
    escape = worker["finite_quotient_escape"]
    gain = int(escape["maximum_marginal_dimension"])
    if gain < 1:
        raise SystemExit("probe has no exact quotient gain beyond rank 17")
    labels = escape["independent_escape_basis_labels"]
    candidate_short = [
        tuple(Q(value) for value in point)
        for point in escape["candidate_points_on_certificate_short_model"]
    ]
    selected_indices = sorted(int(label.rsplit("-", 1)[1]) for label in labels)
    selected_candidates = [candidate_short[index] for index in selected_indices]
    candidate_minimal = worker["exact_nonbaseline_minimal_points"]
    selected_minimal_candidates = [
        candidate_minimal[index] for index in selected_indices
    ]
    specialization_certificate = specialization["finite_quotient_independence"]
    model = tuple(Q(value) for value in specialization_certificate["certificate_short_model"])
    baseline = [
        tuple(Q(value) for value in point)
        for point in specialization_certificate["points"]
    ]
    points = [*baseline, *selected_candidates]
    point_count = len(points)
    exact_admission_gate(
        certified_independent=True, point_count=point_count, quotient_gain=gain
    )
    successful_profile = next(
        profile
        for profile in escape["relation_prime_profiles"]
        if int(profile["marginal_dimension"]) == gain
    )
    relation_prime = int(successful_profile["modulus"])
    prime_bound = int(probe["bounds"]["reduction_prime_bound"])
    certificate = build_finite_quotient_certificate(
        model, points, relation_prime=relation_prime, prime_bound=prime_bound
    )
    verify_finite_quotient_certificate(model, points, certificate)
    exact_admission_gate(
        certified_independent=bool(certificate["certified_independent"]),
        point_count=point_count,
        quotient_gain=gain,
    )

    minimal_model = specialization["global_minimal_specialization"]["model"]
    key = canonical_global_minimal_q_isomorphism_key(minimal_model)
    generated = REPOSITORY / "artifacts/generated-results"
    elliptic_data = REPOSITORY / "elliptic-curves/data"
    archive = REPOSITORY / "archive"
    documented = documented_q_isomorphism_keys(
        (
            REPOSITORY / "MATH_STATUS.json",
            *generated.rglob("*.json"),
            *generated.rglob("*.json.gz"),
            *elliptic_data.rglob("*.json"),
            *elliptic_data.rglob("*.json.gz"),
            *archive.rglob("*.json"),
            *archive.rglob("*.json.gz"),
        ),
        excluded_paths=(
            args.registry,
            args.certificate_output,
            args.specialization,
            args.probe,
        ),
    )
    if key in documented:
        raise SystemExit(
            "curve is already documented in repository manifests/status: "
            + ", ".join(documented[key])
        )

    certificate_record = {
        "schema": "elliptic-curves.q12o5867-rank-jump-certificate.v1",
        "certificate_short_model": [str(value) for value in model],
        "points": [[str(x), str(y)] for x, y in points],
        "finite_quotient_certificate": certificate,
        "source_specialization_artifact": str(args.specialization.resolve()),
        "source_probe_artifact": str(args.probe.resolve()),
    }
    args.certificate_output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite_certificate else "x"
    with args.certificate_output.open(mode) as handle:
        json.dump(certificate_record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    entry = {
        "canonical_global_minimal_q_isomorphism_key": key,
        "global_minimal_model": minimal_model,
        "parameters": [specialization["parameter"]],
        "certified_independent": True,
        "exact_certified_rank_lower_bound": point_count,
        "exact_quotient_gain_beyond_generic_rank17": gain,
        "exact_minimal_model_points": [
            *specialization["global_minimal_specialization"]["points"],
            *selected_minimal_candidates,
        ],
        "artifacts": {
            "specialization": {
                "path": str(args.specialization.resolve()),
                "sha256": sha256_file(args.specialization),
            },
            "probe": {
                "path": str(args.probe.resolve()),
                "sha256": sha256_file(args.probe),
            },
            "exact_certificate": {
                "path": str(args.certificate_output.resolve()),
                "sha256": sha256_file(args.certificate_output),
                "canonical_json_sha256": canonical_json_sha256(certificate_record),
                "point_sequence_sha256": certificate["point_sequence_sha256"],
            },
        },
        "provenance": [
            {
                "reproducing_command": args.reproducing_command,
                "source": (
                    f"bounded q12o5867 {probe_engine} probe plus exact quotient replay"
                ),
            }
        ],
        "claim_boundary": "certified rank lower bound only; no upper bound claimed",
    }
    updated = merge_exact_entry(registry, entry)
    args.registry.write_text(json.dumps(updated, indent=2, sort_keys=True) + "\n")
    print(f"registered_key={key}")
    print(f"certified_rank_lower_bound={point_count}")
    print(f"registry={args.registry.resolve()}")


if __name__ == "__main__":
    main()
