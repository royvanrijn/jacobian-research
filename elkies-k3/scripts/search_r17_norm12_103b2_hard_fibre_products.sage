#!/usr/bin/env sage-python
"""Search the seven hard ``0x103b2`` fibre-product quotients exactly.

status: ACTIVE_SEARCH
claim: bounded exact genus-three quotient search for seven difficult partners
inputs: compiled R17 bisection quartics and the pointed-103b2 helper
outputs: artifacts/generated-results/elkies-k3-r17-norm12-103b2-hard-fibre-products-v1.json

If both quartic covers split at a rational parameter ``t``, their two square
coordinates multiply to a rational point on the genus-three quotient

    z^2 = f_103b2(t) f_partner(t).

Thus a bounded exact search on this quotient is a cheap necessary-condition
sieve before any expensive Mordell--Weil or Chabauty calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import PolynomialRing, QQ, pari


ROOT = Path(__file__).resolve().parents[2]
HELPER_PATH = ROOT / "elkies-k3/scripts/certify_r17_norm12_103b2_jacobian.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-hard-fibre-products-v1.json"
)
HARD_PARTNERS = (
    "norm12-orbit-135b7",
    "norm8-orbit-1723d",
    "norm8-orbit-12e28",
    "norm8-orbit-10092",
    "norm8-orbit-0a0e4",
    "norm12-orbit-06867",
    "norm12-orbit-08f72",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=300_000)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.height < 1:
        parser.error("height must be positive")

    helper = runpy.run_path(str(HELPER_PATH), run_name="r17_jacobian_helpers")
    splitting = json.loads(helper["SPLITTING"].read_text())
    records = {
        record["label"]: record for record in splitting["construction"]["records"]
    }
    missing = set((helper["TARGET_LABEL"], *HARD_PARTNERS)) - set(records)
    if missing:
        raise KeyError(f"missing cover records: {sorted(missing)}")

    ring = PolynomialRing(QQ, "t")
    source = helper["polynomial"](
        records[helper["TARGET_LABEL"]][
            "branch_polynomial_q_coefficients_low_to_high"
        ],
        ring,
    )
    results = []
    for label in HARD_PARTNERS:
        partner = helper["polynomial"](
            records[label]["branch_polynomial_q_coefficients_low_to_high"], ring
        )
        quotient = source * partner
        if quotient.degree() != 8 or not quotient.is_squarefree():
            raise ArithmeticError(f"{label}: quotient model is not squarefree of degree 8")

        raw_points = pari(quotient).hyperellratpoints(args.height)
        quotient_points = []
        simultaneous_splits = []
        for raw_parameter, raw_coordinate in raw_points:
            parameter = QQ(raw_parameter)
            coordinate = QQ(raw_coordinate)
            if coordinate**2 != quotient(parameter):
                raise ArithmeticError(f"{label}: PARI point failed exact replay")
            source_root = helper["rational_square_root"](source(parameter))
            partner_root = helper["rational_square_root"](partner(parameter))
            point_record = {
                "t": helper["rational_text"](parameter),
                "z": helper["rational_text"](coordinate),
                "source_value_is_square": source_root is not None,
                "partner_value_is_square": partner_root is not None,
            }
            quotient_points.append(point_record)
            if source_root is not None and partner_root is not None:
                simultaneous_splits.append({
                    **point_record,
                    "source_square_root": helper["rational_text"](source_root),
                    "partner_square_root": helper["rational_text"](partner_root),
                })
        results.append({
            "partner_cover": label,
            "quotient_genus": 3,
            "quotient_polynomial_coefficients_low_to_high": [
                helper["rational_text"](coefficient) for coefficient in quotient
            ],
            "affine_quotient_point_count": len(quotient_points),
            "affine_quotient_points": quotient_points,
            "simultaneous_split_count": len(simultaneous_splits),
            "simultaneous_splits": simultaneous_splits,
        })

    result = {
        "schema": "elkies-k3.r17-norm12-103b2-hard-fibre-products.v1",
        "status": "PASS_BOUNDED_GENUS_THREE_QUOTIENT_SEARCH",
        "inputs": {
            relative(helper["SPLITTING"]): digest(helper["SPLITTING"]),
            relative(HELPER_PATH): digest(HELPER_PATH),
        },
        "source_cover": helper["TARGET_LABEL"],
        "partner_covers": list(HARD_PARTNERS),
        "pari_hyperellratpoints_naive_height_bound": args.height,
        "results": results,
        "total_affine_quotient_point_count": sum(
            item["affine_quotient_point_count"] for item in results
        ),
        "total_simultaneous_split_count": sum(
            item["simultaneous_split_count"] for item in results
        ),
        "proof_boundary": (
            "For each partner, simultaneous splitting implies an affine rational point "
            "on z^2=f_103b2(t)f_partner(t). PARI/GP hyperellratpoints searched that "
            f"nonsingular genus-three model exactly to naive height {args.height}. "
            "This is an exhaustive bounded search, not a determination of all rational "
            "points on any quotient."
        ),
        "reproducing_command": (
            "sage -python "
            "elkies-k3/scripts/search_r17_norm12_103b2_hard_fibre_products.sage "
            f"--height {args.height}"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != serialized:
            raise ArithmeticError("stored fibre-product result differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17103B2FIBREPRODUCTS|"
        f"partners={len(results)}|"
        f"quotient_points={result['total_affine_quotient_point_count']}|"
        f"simultaneous_splits={result['total_simultaneous_split_count']}|"
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
