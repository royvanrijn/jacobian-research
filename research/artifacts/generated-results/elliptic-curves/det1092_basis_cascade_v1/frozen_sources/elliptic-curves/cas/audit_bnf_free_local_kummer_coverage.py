#!/usr/bin/env sage
"""Certify which stored local Kummer coordinates are covered by known points.

At an odd prime, the formal subgroup ``E_1(Q_p)`` is uniquely 2-divisible.
Thus ``dim_F2 E(Q_p)/2E(Q_p)`` is at most the 2-adic valuation of
``#E(Q_p)/E_1(Q_p)``.  The latter order is the Tamagawa number times the order
of the nonsingular special-fibre group.  If the stored Kummer images of known
rational points attain that upper bound, they span the *full* local Kummer
image at that prime.  This provides an exact BNF-free coverage certificate;
an inequality is reported as unresolved, never as a local condition.

The two-adic formal group is not 2-divisible, so this script deliberately
leaves two-adic coverage unresolved.  It also checks the real component bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ


PROTOCOL = "BNFFREELOCALCOVERAGE"
SIGNATURE_SCHEMA = "elliptic-curves.bnf-free-signature-map.v1"


def rational(value) -> QQ:
    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def f2_rank(masks: list[int], width: int) -> int:
    if not masks or width == 0:
        return 0
    return int(Matrix(GF(2), [[(mask >> index) & 1 for index in range(width)] for mask in masks]).rank())


def project(mask: int, indices: list[int]) -> int:
    return sum(((mask >> index) & 1) << output_index for output_index, index in enumerate(indices))


def local_nonsingular_order(curve, prime: int):
    """Return #E_ns(F_p) from the minimal reduction type at an odd prime."""
    local = curve.local_data(prime)
    if local.has_split_multiplicative_reduction():
        return ZZ(prime - 1), "split_multiplicative"
    if local.has_nonsplit_multiplicative_reduction():
        return ZZ(prime + 1), "nonsplit_multiplicative"
    if local.has_additive_reduction():
        return ZZ(prime), "additive"
    return ZZ(curve.change_ring(GF(prime)).cardinality()), "good"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    record = json.loads(args.signature_map.read_text())
    if not isinstance(record, dict) or record.get("schema") != SIGNATURE_SCHEMA:
        raise ValueError("expected a BNF-free signature-map object")
    coefficients = [rational(value) for value in record["defining_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("signature map must define a monic cubic")
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    cubic = sum(coefficient * x**index for index, coefficient in enumerate(coefficients))
    if not cubic.is_irreducible():
        raise ValueError("signature map cubic must be irreducible")
    curve = EllipticCurve(QQ, [0, coefficients[2], 0, coefficients[1], coefficients[0]])

    local_coordinates = record.get("local_coordinates")
    known = record.get("known_mw_images")
    if not isinstance(local_coordinates, list) or not isinstance(known, list):
        raise ValueError("signature map lacks local coordinates or known MW images")
    local_dimension = int(record.get("local_dimension", len(local_coordinates)))
    if local_dimension != len(local_coordinates):
        raise ValueError("signature local dimension does not match coordinate list")
    known_masks = []
    known_records = []
    for image in known:
        if not isinstance(image, dict):
            raise ValueError("known_mw_images must contain objects")
        mask = int(image["local"], 0) if isinstance(image["local"], str) else int(image["local"])
        if mask < 0 or mask >= 1 << local_dimension:
            raise ValueError("known local image lies outside local coordinate space")
        known_masks.append(mask)
        known_records.append({
            "label": str(image.get("label", len(known_records))),
            "local": hex(mask),
        })

    odd_by_prime: dict[int, list[int]] = {}
    two_indices = []
    real_indices = []
    for index, coordinate in enumerate(local_coordinates):
        if not isinstance(coordinate, dict):
            raise ValueError("local coordinate must be an object")
        kind = coordinate.get("kind")
        if kind in {"odd_valuation_parity", "odd_unit_squareclass"}:
            prime = int(coordinate["rational_prime"])
            odd_by_prime.setdefault(prime, []).append(index)
        elif kind in {"two_adic_product_basis", "two_adic_product_basis_extension"}:
            two_indices.append(index)
        elif kind == "real_sign":
            real_indices.append(index)
        else:
            raise ValueError(f"unsupported local coordinate kind {kind!r}")

    odd_records = []
    for prime, indices in sorted(odd_by_prime.items()):
        if prime == 2:
            raise ValueError("two-adic coordinates must use a two-adic kind")
        nonsingular_order, reduction_type = local_nonsingular_order(curve, prime)
        tamagawa = ZZ(curve.local_data(prime).tamagawa_number())
        upper_bound = int((tamagawa * nonsingular_order).valuation(2))
        rank = f2_rank([project(mask, indices) for mask in known_masks], len(indices))
        odd_records.append(
            {
                "rational_prime": prime,
                "coordinate_indices": indices,
                "reduction_type": reduction_type,
                "tamagawa_number": int(tamagawa),
                "nonsingular_special_fibre_order": int(nonsingular_order),
                "local_kummer_dimension_upper_bound": upper_bound,
                "known_kummer_projection_rank": rank,
                "classification": (
                    "CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE"
                    if rank == upper_bound
                    else "UNRESOLVED_LOCAL_KUMMER_IMAGE_COVERAGE"
                ),
            }
        )

    real_upper_bound = 1 if cubic.discriminant() > 0 else 0
    real_rank = f2_rank([project(mask, real_indices) for mask in known_masks], len(real_indices))
    if real_rank > real_upper_bound:
        raise ValueError("stored real Kummer images exceed the real component bound")
    real_record = {
        "coordinate_indices": real_indices,
        "local_kummer_dimension_upper_bound": real_upper_bound,
        "known_kummer_projection_rank": real_rank,
        "classification": (
            "CERTIFIED_FULL_REAL_KUMMER_IMAGE_COVERAGE"
            if real_rank == real_upper_bound
            else "UNRESOLVED_REAL_KUMMER_IMAGE_COVERAGE"
        ),
    }
    two_record = {
        "coordinate_indices": two_indices,
        "classification": "UNRESOLVED_TWO_ADIC_LOCAL_KUMMER_IMAGE_COVERAGE",
        "reason": "the p=2 formal subgroup is not uniquely 2-divisible",
    }
    output = {
        "protocol": "BNFFREELOCALCOVERAGE-v1",
        "status": "ODD_AND_REAL_LOCAL_COVERAGE_AUDIT_ONLY",
        "signature_map": {
            "path": str(args.signature_map.resolve()),
            "sha256": sha256(args.signature_map.read_bytes()).hexdigest(),
            "source": record.get("source"),
        },
        "signature_local_dimension": local_dimension,
        "known_mw_local_images": known_records,
        "known_mw_image_count": len(known_masks),
        "odd_places": odd_records,
        "real_place": real_record,
        "two_adic_place": two_record,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    certified_odd = sum(
        item["classification"] == "CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE"
        for item in odd_records
    )
    print(
        f"{PROTOCOL}|stage=complete|odd_places={len(odd_records)}"
        f"|certified_odd={certified_odd}"
        f"|real={real_record['classification']}"
        f"|two_adic={two_record['classification']}"
        "|status=PARTIAL_LOCAL_COVERAGE_ONLY",
        flush=True,
    )


if __name__ == "__main__":
    main()
