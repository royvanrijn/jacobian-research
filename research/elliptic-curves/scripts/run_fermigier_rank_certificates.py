#!/usr/bin/env python3
"""Generate exact Fermigier section and E22 independence certificates."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
POINT_DATA = PROGRAM_ROOT / "data" / "fermigier_e22_points.json"
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_REPORTED_PARAMETER,
    thirteenth_visible_point,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    specialize_fermigier_rank_sections,
    write_json_exclusively,
)
from ecsearch.rank_certification import (  # noqa: E402
    build_independence_certificate,
)


def rational_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def build_manifest(maximum_reduction_prime: int = 1000) -> dict[str, object]:
    raw_point_data = POINT_DATA.read_bytes()
    point_data = json.loads(raw_point_data)
    e22_model = tuple(map(Fraction, point_data["weierstrass_coefficients"]))
    e22_points = tuple(
        (Fraction(point[0]), Fraction(point[1])) for point in point_data["points"]
    )

    specialization = specialize_fermigier_rank_sections(
        FERMIGIER_REPORTED_PARAMETER
    )
    section_certificate = build_independence_certificate(
        specialization.canonical_model,
        specialization.section_differences,
        relation_prime=5,
        maximum_reduction_prime=maximum_reduction_prime,
    )
    e22_certificate = build_independence_certificate(
        e22_model,
        e22_points,
        relation_prime=2,
        maximum_reduction_prime=maximum_reduction_prime,
    )
    thirteenth = thirteenth_visible_point(specialization.quartic_model)
    return {
        "schema": "elliptic-curves.fermigier-rank-certificates.v1",
        "claim_level": "exact_rank_lower_bound_replay",
        "generator": "elliptic-curves/scripts/run_fermigier_rank_certificates.py",
        "canonical_pinned_command": (
            "python3 elliptic-curves/scripts/run_fermigier_rank_certificates.py "
            "--output artifacts/generated-results/elliptic-curves/"
            "fermigier_rank_certificates_v1.json"
        ),
        "randomness": "none; reduction primes and generators are scanned deterministically",
        "input": {
            "published_points": "elliptic-curves/data/fermigier_e22_points.json",
            "published_points_sha256": hashlib.sha256(raw_point_data).hexdigest(),
            "source": point_data["source"],
        },
        "search_bound": {
            "maximum_reduction_prime": maximum_reduction_prime,
            "interpretation": (
                "certificate-construction bound only; the replay checks every stored row"
            ),
        },
        "generic_sections": {
            "adapter_parameter": rational_text(FERMIGIER_REPORTED_PARAMETER),
            "literal_shift": rational_text(FERMIGIER_E22_RECONSTRUCTION_SHIFT),
            "quartic_point_count": len(specialization.quartic_points),
            "section_difference_count": len(specialization.section_differences),
            "thirteenth_quartic_point": [
                rational_text(thirteenth[0]),
                rational_text(thirteenth[1]),
            ],
            "certificate": section_certificate.to_json_object(),
            "conclusion": (
                "the twelve specialized section differences are independent; "
                "because the thirteen points and covariant map are rational in the "
                "family parameter, the generic Mordell-Weil rank is at least twelve"
            ),
        },
        "published_e22_points": {
            "point_count": len(e22_points),
            "certificate": e22_certificate.to_json_object(),
            "conclusion": "the displayed E22 curve has Mordell-Weil rank at least 22",
        },
        "limitations": {
            "upper_bound": (
                "no unconditional rank upper bound is supplied; exact rank 22 remains "
                "conditional as recorded in the benchmark metadata"
            ),
            "saturation": (
                "the certificates prove independence of the displayed points, not "
                "that they form a saturated Mordell-Weil basis"
            ),
            "normalization": (
                "the unresolved factor-two discrepancy in the printed E22 family "
                "parameter is unchanged"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-reduction-prime", type=int, default=1000)
    args = parser.parse_args()
    result = build_manifest(args.maximum_reduction_prime)
    generic = result["generic_sections"]
    e22 = result["published_e22_points"]
    print(
        "FERMIGIER_RANK_CERTIFICATES "
        f"generic_differences={generic['section_difference_count']} "
        f"e22_points={e22['point_count']}"
    )
    if args.output is not None:
        text = json.dumps(result, indent=2, sort_keys=True) + "\n"
        write_json_exclusively(args.output, text)
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
