#!/usr/bin/env python3
"""Build a portable exact rank-29 certificate for the 2024 record curve.

The public announcement certified independence numerically through a canonical
height determinant.  This replay checks all coordinates exactly and supplies a
separate finite-reduction proof.  If a relation among the 29 points existed,
the full-rank images in ``prod E(F_p)/2E(F_p)`` would make every coefficient
even.  The irreducible reduced 2-division cubic at a separate prime proves
``E(Q)[2]=0``; division by two and infinite descent then kill the relation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

from elkies_klagsbrun_rank29 import (
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    PUBLISHED_POINTS,
    curve_discriminant,
    point_on_general_curve,
    point_on_short_curve,
    published_short_points,
    short_weierstrass_coefficients,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic-curves/elliptic_elkies_klagsbrun_rank29_certificate.json"
)
SOURCE_URLS = (
    "https://web.math.pmf.unizg.hr/~duje/tors/rk29.html",
    "https://mathoverflow.net/questions/477849/background-for-the-elkies-klagsbrun-curve-of-rank-29/478050",
)
DISCRIMINANT_FACTORIZATION = (
    (2, 19),
    (3, 7),
    (5, 7),
    (7, 4),
    (11, 5),
    (13, 3),
    (17, 4),
    (31, 3),
    (41, 2),
    (43, 2),
    (61, 2),
    (233, 1),
    (241, 2),
    (4139, 1),
    (678146849364709860535420504397393, 1),
    (159788990966780131363155786084695062643236502969, 1),
    (4402149008473369392540402625019227412319473055901, 1),
)


def rational_string(value: Any) -> str:
    return str(value)


def points_sha256(points: Any) -> str:
    payload = json.dumps(
        [[rational_string(x), rational_string(y)] for x, y in points],
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def build_artifact(*, certificate_prime_bound: int = 500) -> dict[str, Any]:
    if len(PUBLISHED_POINTS) != 29:
        raise AssertionError("the public data must contain exactly 29 points")
    if not all(point_on_general_curve(point) for point in PUBLISHED_POINTS):
        raise AssertionError("a published point is off the announced curve")

    short_coefficients = short_weierstrass_coefficients()
    short_points = published_short_points()
    if not all(point_on_short_curve(point) for point in short_points):
        raise AssertionError("the short-model transport failed")

    signatures = find_mod2_reduction_certificate(
        short_coefficients,
        short_points,
        prime_bound=certificate_prime_bound,
    )
    binary_rank = combined_mod2_rank(signatures, len(short_points))
    if binary_rank != 29:
        raise AssertionError(
            f"finite reductions reached rank {binary_rank}, not rank 29"
        )
    two_torsion_prime = find_two_torsion_certificate_prime(
        short_coefficients, prime_bound=certificate_prime_bound
    )

    discriminant = curve_discriminant()
    factor_product = -1
    for prime, exponent in DISCRIMINANT_FACTORIZATION:
        factor_product *= prime**exponent
    if factor_product != discriminant:
        raise AssertionError("the announced discriminant factorization changed")

    script_path = Path(__file__).resolve()
    data_path = script_path.with_name("elkies_klagsbrun_rank29.py")
    return {
        "schema_version": 1,
        "artifact_kind": "exact_elliptic_curve_rank_lower_bound_certificate",
        "status": "exact_unconditional_rank_at_least_29",
        "claim": {
            "curve": "y^2+x*y=x^3+A*x+B",
            "certified_algebraic_rank_lower_bound": 29,
            "target_rank_30_achieved": False,
            "conditional_exact_rank_statement_not_used": (
                "The public exact-rank-29 statement assumes GRH; this artifact "
                "uses only the unconditional lower-bound data."
            ),
        },
        "sources": list(SOURCE_URLS),
        "reproduction": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/verify_elkies_klagsbrun_rank29.py"
            ),
            "python": sys.version.split()[0],
            "certificate_prime_bound": certificate_prime_bound,
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
            "data_module_sha256": hashlib.sha256(data_path.read_bytes()).hexdigest(),
        },
        "general_weierstrass_coefficients": [
            rational_string(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS
        ],
        "integral_short_weierstrass_coefficients": [
            rational_string(value) for value in short_coefficients
        ],
        "short_model_transport": {
            "X": "36*x+3",
            "Y": "108*(2*y+x)",
            "isomorphism_over_Q": True,
        },
        "published_points": [
            {"x": rational_string(x), "y": rational_string(y)}
            for x, y in PUBLISHED_POINTS
        ],
        "published_points_sha256": points_sha256(PUBLISHED_POINTS),
        "exact_membership_checks_passed": 29,
        "discriminant": str(discriminant),
        "announced_discriminant_factorization": [
            {"factor": str(prime), "exponent": exponent}
            for prime, exponent in DISCRIMINANT_FACTORIZATION
        ],
        "factor_product_equals_discriminant": True,
        "finite_reduction_certificate": {
            "proof": (
                "Full F2 column rank forces every coefficient of a rational "
                "relation to be even.  The separate no-2-torsion certificate "
                "allows division by two; infinite descent proves independence."
            ),
            "two_torsion_certificate_prime": two_torsion_prime,
            "reduced_2_division_cubic_has_no_root": True,
            "certificate_primes": [signature.prime for signature in signatures],
            "signatures": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
            "combined_exact_rank_over_F2": binary_rank,
            "certified_algebraic_rank_lower_bound": 29,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    args = parser.parse_args()
    artifact = build_artifact(certificate_prime_bound=args.certificate_prime_bound)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    certificate = artifact["finite_reduction_certificate"]
    print(f"wrote {args.output}")
    print(
        "exact rank lower bound=29; certificate primes="
        f"{certificate['certificate_primes']}"
    )
    print("rank_30_target_hit=false")


if __name__ == "__main__":
    main()
