#!/usr/bin/env python3
"""Replay the four high-rank fibres disclosed in arXiv:2608.25406v1.

All four fibres are specialized from the compact published model, globally
minimalized with an exact PARI change, and supplied with an unconditional
finite-quotient independence certificate for the 17 generic sections.  At the
four parameters the minimal model is checked coefficient-for-coefficient
against the corresponding public Elkies curve in Dujella's rank-record table.
The generic 17 sections are then extended by a deterministic subset of those
public points.  One combined finite-quotient certificate proves quotient gains
8, 9, 10 and 11, rather than merely certifying the two point lists separately.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(ELLIPTIC_ROOT / "cas")]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
    weierstrass_invariants,
)
from elkies_rank25 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as RANK25_MODEL,
    POINTS as RANK25_POINTS,
)
from elkies_rank26 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as RANK26_MODEL,
    POINTS as RANK26_POINTS,
)
from elkies_rank27 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as RANK27_MODEL,
    POINTS as RANK27_POINTS,
)
from elkies_rank28 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as RANK28_MODEL,
    POINTS as RANK28_POINTS,
)


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
PARAMETERS = (
    ("rank_at_least_25", 25, -2, 377),
    ("rank_at_least_26", 26, -308, 251),
    ("rank_at_least_27", 27, 2456, 135),
    ("rank_at_least_28", 28, -9529, 5471),
)
PUBLIC_CONTROLS = {
    25: (RANK25_MODEL, RANK25_POINTS, "https://web.math.pmf.unizg.hr/~duje/tors/rk25.html"),
    26: (RANK26_MODEL, RANK26_POINTS, "https://web.math.pmf.unizg.hr/~duje/tors/rk26.html"),
    27: (RANK27_MODEL, RANK27_POINTS, "https://web.math.pmf.unizg.hr/~duje/tors/rk27.html"),
    28: (RANK28_MODEL, RANK28_POINTS, "https://web.math.pmf.unizg.hr/~duje/tors/rk28.html"),
}


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def model_record(model) -> list[str]:
    return [str(value) for value in model]


def certificate_for(model, points, *, prime_bound: int = 500):
    short_model, short_change = short_certificate_model(model)
    short_points = tuple(source_point_to_target(point, short_change) for point in points)
    certificate = build_finite_quotient_certificate(
        short_model, short_points, relation_prime=2, prime_bound=prime_bound
    )
    verify_finite_quotient_certificate(short_model, short_points, certificate)
    if not certificate["certified_independent"]:
        raise ArithmeticError("the expected mod-2 independence certificate failed")
    return short_model, certificate


def quotient_complement_certificate(model, generic_points, public_points, target_rank):
    """Select public points extending the fixed generic basis, then certify it.

    The all-points reduction matrix is scanned with the generic points first.
    Since its first seventeen pivot columns are therefore locked, the remaining
    pivot columns give a deterministic public complement.  The selected union
    is replayed through the ordinary independence verifier as a standalone
    certificate.
    """

    short_model, short_change = short_certificate_model(model)
    all_points = tuple(generic_points) + tuple(public_points)
    all_short_points = tuple(
        source_point_to_target(point, short_change) for point in all_points
    )
    discovery = build_finite_quotient_certificate(
        short_model, all_short_points, relation_prime=2, prime_bound=1000
    )
    pivots = tuple(discovery["pivot_columns_zero_based"])
    if discovery["combined_rank_over_relation_field"] != target_rank:
        raise ArithmeticError("the public control did not reach its advertised rank")
    if pivots[:17] != tuple(range(17)):
        raise ArithmeticError("the fixed generic basis was not preserved by pivoting")
    public_indices = tuple(index - 17 for index in pivots[17:])
    if len(public_indices) != target_rank - 17:
        raise ArithmeticError("the public complement has the wrong quotient dimension")
    combined_points = tuple(generic_points) + tuple(
        public_points[index] for index in public_indices
    )
    combined_short_points = tuple(
        source_point_to_target(point, short_change) for point in combined_points
    )
    certificate = build_finite_quotient_certificate(
        short_model, combined_short_points, relation_prime=2, prime_bound=1000
    )
    verify_finite_quotient_certificate(short_model, combined_short_points, certificate)
    if not certificate["certified_independent"]:
        raise ArithmeticError("the generic-plus-public complement was not certified")
    return short_model, public_indices, certificate


def main() -> None:
    sys.set_int_max_str_digits(0)
    data = load_q12o5867_data(MODEL, SECTIONS)
    rows = []
    for label, published_rank, numerator, denominator in PARAMETERS:
        specialization = evaluate_projective_specialization(
            data, numerator, denominator
        )
        minimal_model, minimal_change, minimal_metadata = global_minimal_model_with_change(
            specialization.model
        )
        minimal_points = tuple(
            source_point_to_target(point, minimal_change)
            for point in specialization.points
        )
        if any(
            not is_on_weierstrass_curve(minimal_model, point)
            for point in minimal_points
        ):
            raise AssertionError("a transported published section missed its minimal fibre")
        short_model, baseline_certificate = certificate_for(
            minimal_model, minimal_points
        )
        public_model, public_points, source_url = PUBLIC_CONTROLS[published_rank]
        if minimal_model != tuple(public_model):
            raise ArithmeticError(
                f"the rank-{published_rank} specialization is not its public minimal model"
            )
        if len(public_points) != published_rank or any(
            not is_on_weierstrass_curve(minimal_model, point)
            for point in public_points
        ):
            raise ArithmeticError(
                f"the public rank-{published_rank} point list failed exact replay"
            )
        public_short, public_certificate = certificate_for(
            minimal_model, public_points, prime_bound=1000
        )
        combined_short, complement_indices, combined_certificate = (
            quotient_complement_certificate(
                minimal_model, minimal_points, public_points, published_rank
            )
        )
        if public_short != combined_short or public_short != short_model:
            raise AssertionError("short certificate models unexpectedly differ")
        quotient_gain = len(complement_indices)
        row = {
            "label": label,
            "parameter": f"{numerator}/{denominator}",
            "projective_parameter": [numerator, denominator],
            "published_rank_lower_bound": published_rank,
            "minimal_model": model_record(minimal_model),
            "minimal_discriminant": str(
                weierstrass_invariants(minimal_model)["discriminant"]
            ),
            "minimalization": minimal_metadata,
            "generic_sections": {
                "count": 17,
                "short_certificate_model": model_record(short_model),
                "finite_quotient_independence": baseline_certificate,
            },
            "locally_certified_rank_lower_bound": published_rank,
            "public_positive_control": {
                "source_url": source_url,
                "exact_global_minimal_model_equality": True,
                "public_point_count": len(public_points),
                "public_basis_finite_quotient_independence": public_certificate,
                "selected_public_point_indices_one_based": [
                    index + 1 for index in complement_indices
                ],
                "quotient_gain_beyond_generic_rank_17": quotient_gain,
                "combined_generic_plus_complement_independence": combined_certificate,
            },
            "rank_validation_boundary": (
                f"An exact public rank-{published_rank} point list is imported from "
                "Dujella's table. A single finite-reduction matrix certifies the "
                f"generic 17 plus {quotient_gain} selected public directions. This is "
                "an unconditional rank lower bound, not an upper bound or exact rank."
            ),
        }
        rows.append(row)

    payload = {
        "schema": "elliptic-curves.elkies-2026-high-rank-positive-controls.v2",
        "status": "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS",
        "inputs": {
            str(MODEL): file_sha256(MODEL),
            str(SECTIONS): file_sha256(SECTIONS),
            **{
                f"elliptic-curves/cas/elkies_rank{rank}.py": file_sha256(
                    ELLIPTIC_ROOT / "cas" / f"elkies_rank{rank}.py"
                )
                for rank in range(25, 29)
            },
        },
        "coordinate": data.coordinate,
        "fibres": rows,
        "search_policy": {
            "use_as_calibration_anchors": True,
            "preferred_over_icarm_curve_inference": True,
            "rank28_is_calibration_not_open_four_point_search": True,
            "required_quotient_gains": [8, 9, 10, 11],
        },
        "proof_boundary": (
            "Unconditional finite-reduction rank lower bounds are 25,26,27,28. "
            "No exact-rank or Selmer upper-bound claim is made."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLS|local_bounds=25,26,27,28|quotient_gains=8,9,10,11|"
        f"models_exact=true|status={payload['status']}|output={OUTPUT}"
    )


if __name__ == "__main__":
    main()
