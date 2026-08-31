#!/usr/bin/env python3
"""Build eleven exact residual-Selmer positive controls at Elkies rank 28.

Each certified public complement point ``Q`` gives the cubic Kummer element
``alpha = X(Q)-theta`` on the monic descent model.  Its norm is exactly
``Z(Q)^2``, and the associated two-cover has the rational point
``[u:v:w:z]=[1:0:0:1]``.  The already replayed finite-reduction certificate
proves that these eleven point classes are independent modulo the generic
seventeen.

This is a lower-bound calibration for the descent machinery.  It supplies no
ambient ``K(S,2)`` enumeration, Selmer upper bound, or search authorization.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_public11_selmer_candidates_v1.json"
)
OUTPUT_SCHEMA = "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1"
PROTOCOL = "ELKIESR28PUBLICSELMER"

sys.path.insert(0, str(CAS))

from build_elkies_2026_rank28_relative_descent_magma import (  # noqa: E402
    KNOWN_QUOTIENT_GAIN,
    load_relative_input,
)
from build_q12o5867_bnf_free_signature import (  # noqa: E402
    evaluate_cubic,
    monic_cubic_coefficients,
    point_on_monic_cubic,
)


def build_controls() -> dict:
    source = load_relative_input()
    model = tuple(Fraction(value) for value in source.model)
    coefficients = monic_cubic_coefficients(model)
    candidates = []
    for index, point in enumerate(source.public_complement, start=1):
        x_coordinate, z_coordinate = point_on_monic_cubic(model, point)
        norm = evaluate_cubic(coefficients, x_coordinate)
        if norm != z_coordinate**2:
            raise ArithmeticError("a public complement point failed cubic transport")
        candidates.append(
            {
                "label": f"public-complement-Q{index}",
                "generator_coefficients": [str(x_coordinate), "-1", "0"],
                "norm": str(norm),
                "norm_square_root": str(z_coordinate),
                "normalization": "alpha=X(Q)-theta",
                "source_elliptic_point": [str(value) for value in point],
                "monic_cubic_point": [str(x_coordinate), str(z_coordinate)],
                "rational_cover_witness": ["1", "0", "0", "1"],
            }
        )
    if len(candidates) != KNOWN_QUOTIENT_GAIN:
        raise ArithmeticError("the public residual control count is no longer eleven")
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "EXACT_MW_KUMMER_POSITIVE_CONTROLS",
        "parameter": "-9529/5471",
        "minimal_model": [str(value) for value in model],
        "field_polynomial_ascending": [str(value) for value in coefficients],
        "candidate_count": len(candidates),
        "candidates": candidates,
        "independence_certificate": {
            "method": "replayed finite-reduction independence modulo the generic seventeen",
            "generic_rank": 17,
            "public_complement_quotient_gain": KNOWN_QUOTIENT_GAIN,
            "combined_rank": 28,
            "controls_sha256": source.controls_sha256,
            "generic_point_sequence_sha256": source.generic_point_sequence_sha256,
            "combined_point_sequence_sha256": source.combined_point_sequence_sha256,
        },
        "residual_two_selmer_dimension_lower_bound": KNOWN_QUOTIENT_GAIN,
        "expensive_search_authorized": False,
        "claim_boundary": [
            "Every candidate is the exact Kummer image X(Q)-theta of a certified rational point.",
            "Every norm is the displayed exact square Z(Q)^2 and every cover has the rational witness [1:0:0:1].",
            "The eleven classes are independent modulo the generic seventeen by the replayed finite-reduction certificate.",
            "This proves only a residual 2-Selmer lower bound of eleven at the positive-control fibre, not an upper bound or rank 32.",
            "No point search is authorized.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    output = build_controls()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"{PROTOCOL}|stage=complete|classes={output['candidate_count']}"
        f"|residual_lower_bound={output['residual_two_selmer_dimension_lower_bound']}"
        "|search_authorized=false"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
