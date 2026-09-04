#!/usr/bin/env sage
"""Certify the class-group pressure already forced by known R17 points.

Status: ACTIVE_PROOF
Claim: unconditional lower bounds for auxiliary cubic class-group 2-ranks.
Inputs: pinned public-fibre and exact-lineage certificates.
Output: artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-v1.json.

For an integral elliptic curve with irreducible completed-square 2-division
cubic ``K = Q(zeta)``, a point has Kummer representative

    alpha(P) = 4*x(P) - zeta,

whose norm is a square.  Away from 2 and the discriminant, every valuation of
``alpha(P)`` is even.  Hence a combination of independent point classes whose
valuation parities also vanish at the bad primes determines a 2-torsion ideal
class.  Its norm is a positive square, so the kernel consists of norm-positive
unit squareclasses.  Since the cubic has odd degree and ``Norm(-1)=-1``, that
kernel has dimension ``r1+r2-1``.  This gives a rigorous lower bound for
``Cl(K)[2]`` without computing a BNF.

The replay also constructs the point half-ideals

    A_P = (d^2*alpha(P), d^3*4*(2*y+a1*x+a3)),

where ``d^2`` is the denominator of ``4*x``.  It verifies that
``A_P^2/(d^2*alpha(P))`` is supported only above the declared bad primes.
These ideals are exact seeds for a future S-class collector which quotients by
the known Mordell--Weil image before relation collection.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, pari, vector
from sage.version import version as sage_version


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
LINEAGE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
SOURCE = Path(__file__).resolve()
TARGET_IDS = (351, 356, 376, 377, 385)
GENERIC_RANK = 17
PUBLIC_STATUS = "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES"
LINEAGE_STATUS = "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS"
SCHEMA = "elkies-k3.r17-kummer-classgroup-pressure.v1"
STATUS = "PROVED_KUMMER_FORCED_CUBIC_CLASS_GROUP_2RANK_LOWER_BOUNDS"
PROTOCOL = "R17KUMMERCL2PRESSURE"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def divide_supported_rational(value: QQ, primes: list[ZZ]) -> list[dict[str, int]]:
    """Prove that a nonzero rational is supported on ``primes``."""

    numerator = abs(ZZ(value.numerator()))
    denominator = ZZ(value.denominator())
    support = []
    for prime in primes:
        exponent = 0
        while numerator % prime == 0:
            numerator //= prime
            exponent += 1
        while denominator % prime == 0:
            denominator //= prime
            exponent -= 1
        if exponent:
            support.append({"prime": int(prime), "exponent": exponent})
    if numerator != 1 or denominator != 1:
        raise ArithmeticError(
            f"rational ideal norm has support outside the bad set: "
            f"{numerator}/{denominator}"
        )
    return support


def certified_point_mod2_rank(independence: dict[str, Any], count: int) -> int:
    signature_rows = [
        [int(bit) for bit in row]
        for signature in independence["mod2_reduction_signatures"]
        for row in signature["rows"]
    ]
    columns = [
        [row[index] for row in signature_rows]
        for index in range(count)
    ]
    return int(Matrix(GF(2), columns).rank())


def half_ideal_record(nf, theta, ainvs, point, bad_primes):
    a1, _a2, a3, _a4, _a6 = ainvs
    x_coordinate, y_coordinate = map(QQ, point)
    scaled_x = 4 * x_coordinate
    denominator_root = ZZ(scaled_x.denominator()).isqrt()
    if denominator_root**2 != scaled_x.denominator():
        raise ArithmeticError("the denominator of 4*x is not a square")

    integral_alpha = (
        pari(str(denominator_root**2 * scaled_x))
        - denominator_root**2 * theta
    )
    integral_norm_root_q = QQ(
        denominator_root**3
        * 4
        * (2 * y_coordinate + a1 * x_coordinate + a3)
    )
    if integral_norm_root_q.denominator() != 1:
        raise ArithmeticError("the scaled Kummer norm root is not integral")
    integral_norm_root = pari(str(integral_norm_root_q))
    if QQ(str(pari.nfeltnorm(nf, integral_alpha))) != integral_norm_root_q**2:
        raise ArithmeticError("the integral Kummer norm identity failed")

    half_ideal = pari.idealadd(nf, integral_alpha, integral_norm_root)
    correction = pari.idealdiv(
        nf,
        pari.idealpow(nf, half_ideal, 2),
        integral_alpha,
    )
    correction_norm = QQ(str(pari.idealnorm(nf, correction)))
    correction_support = divide_supported_rational(correction_norm, bad_primes)
    return {
        "denominator_root_for_4x": str(denominator_root),
        "integral_alpha": str(integral_alpha),
        "integral_norm_root": str(integral_norm_root_q),
        "half_ideal_hnf": str(half_ideal),
        "localized_square_correction_ideal_hnf": str(correction),
        "localized_square_correction_norm": str(correction_norm),
        "localized_square_correction_rational_support": correction_support,
    }


def audit_curve(public_record, independence):
    curve_id = int(public_record["id"])
    ainvs = tuple(QQ(value) for value in public_record["ainvs"])
    curve = EllipticCurve(QQ, list(ainvs))
    points = public_record["points"]
    point_count = len(points)
    if point_count != int(public_record["snapshot_rank_lower_bound"]):
        raise ArithmeticError("the public point count and rank lower bound differ")
    if point_count != int(independence["proved_displayed_subgroup_rank"]):
        raise ArithmeticError("the lineage rank certificate changed")
    if certified_point_mod2_rank(independence, point_count) != point_count:
        raise ArithmeticError("the displayed points are not certified independent modulo 2")
    if any(curve(QQ(point[0]), QQ(point[1])) not in curve for point in points):
        raise ArithmeticError("a public point is not on its curve")

    a1, a2, a3, a4, a6 = ainvs
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    polynomial = z**3 + ZZ(b2) * z**2 + ZZ(8 * b4) * z + ZZ(16 * b6)
    if not polynomial.is_irreducible():
        raise ArithmeticError("the completed-square 2-division cubic is reducible")
    if abs(ZZ(polynomial.discriminant())) != 256 * abs(ZZ(curve.discriminant())):
        raise ArithmeticError("the completed-square discriminant identity failed")

    bad_primes = [ZZ(value) for value in public_record["bad_primes"]]
    remainder = abs(ZZ(curve.discriminant()))
    for prime in bad_primes:
        if not prime.is_prime(proof=True):
            raise ArithmeticError(f"the declared bad factor {prime} is not prime")
        while remainder % prime == 0:
            remainder //= prime
    if remainder != 1:
        raise ArithmeticError("the public bad-prime list is incomplete")

    pari.addprimes(bad_primes)
    nf = pari.nfinit([pari(polynomial), bad_primes])
    if list(pari.nfcertify(nf)):
        raise ArithmeticError("the cubic maximal order failed certification")
    signature = [int(value) for value in nf.nf_get_sign()]
    root_number = int(curve.root_number())
    total_two_selmer_parity = 0 if root_number == 1 else 1
    unit_squareclass_dimension = sum(signature)
    norm_square_unit_squareclass_dimension = unit_squareclass_dimension - 1
    if norm_square_unit_squareclass_dimension < 0:
        raise ArithmeticError("invalid cubic signature")
    theta = pari(f"Mod(z,{polynomial})")

    places = []
    place_records = []
    for prime in bad_primes:
        for index, place in enumerate(pari.idealprimedec(nf, prime), start=1):
            places.append(place)
            place_records.append(
                {
                    "column": len(places) - 1,
                    "rational_prime": str(prime),
                    "prime_ideal_index_one_based": index,
                    "ramification_index": int(place[2]),
                    "residue_degree": int(place[3]),
                    "prime_ideal": str(place),
                }
            )

    valuation_rows = []
    half_ideals = []
    for index, point in enumerate(points, start=1):
        x_coordinate, y_coordinate = map(QQ, point)
        alpha = pari(str(4 * x_coordinate)) - theta
        norm = QQ(str(pari.nfeltnorm(nf, alpha)))
        norm_root = 4 * (2 * y_coordinate + a1 * x_coordinate + a3)
        if norm != norm_root**2:
            raise ArithmeticError("a Kummer norm-square identity failed")
        valuation_rows.append(
            [int(pari.idealval(nf, alpha, place)) & 1 for place in places]
        )
        half_ideals.append(
            {
                "label": f"P{index}",
                **half_ideal_record(nf, theta, ainvs, point, bad_primes),
            }
        )

    valuation_matrix = Matrix(GF(2), valuation_rows)
    valuation_rank = int(valuation_matrix.rank())
    kernel = valuation_matrix.left_kernel()
    kernel_rows = [
        [int(value) for value in row]
        for row in kernel.basis()
    ]
    kernel_dimension = len(kernel_rows)
    class_group_lower = max(
        0, kernel_dimension - norm_square_unit_squareclass_dimension
    )
    generic_matrix = Matrix(GF(2), valuation_rows[:GENERIC_RANK])
    generic_rank = int(generic_matrix.rank())
    residual_rows = valuation_rows[GENERIC_RANK:]
    residual_rank = int(Matrix(GF(2), residual_rows).rank()) if residual_rows else 0
    residual_gain = point_count - GENERIC_RANK
    residual_valuation_rank_modulo_generic = valuation_rank - generic_rank
    everywhere_even_residual_dimension = (
        residual_gain - residual_valuation_rank_modulo_generic
    )
    residual_class_group_image_lower = max(
        0,
        everywhere_even_residual_dimension
        - norm_square_unit_squareclass_dimension,
    )
    residual_adjustments = []
    if residual_valuation_rank_modulo_generic == 0:
        transpose = generic_matrix.transpose()
        for offset, row in enumerate(residual_rows, start=GENERIC_RANK + 1):
            coefficients = transpose.solve_right(vector(GF(2), row))
            if vector(GF(2), coefficients) * generic_matrix != vector(GF(2), row):
                raise ArithmeticError("a residual valuation adjustment failed")
            residual_adjustments.append(
                {
                    "residual_label": f"P{offset}",
                    "mw17_labels_to_add": [
                        f"P{index + 1}"
                        for index, value in enumerate(coefficients)
                        if value
                    ],
                    "adjusted_bad_valuation_parity_row": [0] * len(places),
                }
            )

    return {
        "curve_id": curve_id,
        "role": "rank-29-record" if curve_id in (356, 385) else "exact-control",
        "global_minimal_model": [str(value) for value in ainvs],
        "point_count": point_count,
        "residual_gain_over_mw17": residual_gain,
        "two_division_cubic": str(polynomial),
        "field_signature": signature,
        "root_number": root_number,
        "rational_two_torsion_dimension": 0,
        "proved_total_two_selmer_dimension_mod_2": total_two_selmer_parity,
        "unit_squareclass_dimension": unit_squareclass_dimension,
        "norm_square_unit_squareclass_dimension": (
            norm_square_unit_squareclass_dimension
        ),
        "bad_rational_prime_count": len(bad_primes),
        "bad_prime_ideal_count": len(places),
        "bad_prime_ideal_columns": place_records,
        "known_point_bad_valuation_parity_rows": [
            {"label": f"P{index + 1}", "row": row}
            for index, row in enumerate(valuation_rows)
        ],
        "known_point_bad_valuation_rank": valuation_rank,
        "generic_mw17_bad_valuation_rank": generic_rank,
        "residual_point_bad_valuation_rank": residual_rank,
        "residual_bad_valuation_rank_modulo_mw17": (
            residual_valuation_rank_modulo_generic
        ),
        "everywhere_even_residual_quotient_dimension": (
            everywhere_even_residual_dimension
        ),
        "proved_adjusted_residual_class_group_image_dimension_lower_bound": (
            residual_class_group_image_lower
        ),
        "residual_adjustments_by_mw17": residual_adjustments,
        "everywhere_even_known_kummer_kernel": {
            "dimension": kernel_dimension,
            "basis_rows_in_P1_through_Pn_coordinates": kernel_rows,
        },
        "proved_class_group_2rank_lower_bound": class_group_lower,
        "point_half_ideals": half_ideals,
        "proof_identity": (
            "dim image(kernel(known Kummer classes -> bad-prime valuation "
            "parities) -> Cl(K)[2]) is at least "
            "kernel_dimension-(r1+r2-1); the extra subtraction is avoided "
            "because every Kummer norm is a positive square while Norm(-1)=-1"
        ),
    }


def build():
    public = json.loads(PUBLIC.read_text())
    lineage = json.loads(LINEAGE.read_text())
    if public.get("status") != PUBLIC_STATUS:
        raise ArithmeticError("the public-fibre certificate is not passing")
    if lineage.get("status") != LINEAGE_STATUS:
        raise ArithmeticError("the lineage certificate is not passing")
    public_by_id = {int(record["id"]): record for record in public["records"]}
    independence_by_id = {
        int(record["curve_id"]): record
        for record in lineage["displayed_point_independence"]
    }
    curves = [
        audit_curve(public_by_id[curve_id], independence_by_id[curve_id])
        for curve_id in TARGET_IDS
    ]
    return {
        "schema": SCHEMA,
        "status": STATUS,
        "summary": {
            "curve_ids": list(TARGET_IDS),
            "proved_class_group_2rank_lower_bounds": {
                str(record["curve_id"]): record[
                    "proved_class_group_2rank_lower_bound"
                ]
                for record in curves
            },
            "record_curve_bounds": {
                str(record["curve_id"]): record[
                    "proved_class_group_2rank_lower_bound"
                ]
                for record in curves
                if record["curve_id"] in (356, 385)
            },
            "operational_conclusion": (
                "A full BNF computes a large auxiliary cubic 2-class group which "
                "the known Kummer classes already force.  A residual descent should "
                "localize at S and quotient by the certified point half-ideals before "
                "collecting the remaining class relations."
            ),
        },
        "curves": curves,
        "method": {
            "kummer_representative": "4*x(P)-zeta",
            "unramified_outside": "2 and the elliptic discriminant",
            "valuation_matrix": "prime-ideal valuation parities at every bad rational prime",
            "class_group_map": (
                "an everywhere-even squareclass alpha maps to the class of the "
                "fractional ideal with square (alpha)"
            ),
            "kernel_upper_bound": (
                "the norm-positive unit squareclass dimension r1+r2-1; "
                "Norm(-1)=-1 in odd degree"
            ),
            "half_ideal": "(d^2*(4*x-zeta), d^3*4*(2*y+a1*x+a3))",
            "selmer_parity": (
                "the 2-parity theorem over Q identifies dim Sel_2(E/Q) "
                "modulo 2 with the exact root number when E(Q)[2]=0"
            ),
        },
        "claim_boundary": [
            "The displayed values are unconditional lower bounds for the full cubic ideal class group 2-ranks, not exact class groups.",
            "They are not lower bounds for the S-class group after the bad primes are inverted; bad-prime ideal classes may kill some or all of these directions.",
            "They do not compute a Selmer upper bound, a complete residual Selmer group, or an exact elliptic-curve rank.",
            "The point half-ideals are exact localized square roots and algorithmic seeds; relation collection and all local Selmer conditions remain to be completed.",
        ],
        "inputs": {
            str(PUBLIC.relative_to(ROOT)): digest(PUBLIC),
            str(LINEAGE.relative_to(ROOT)): digest(LINEAGE),
            str(SOURCE.relative_to(ROOT)): digest(SOURCE),
        },
        "software_assumptions": {
            "sage": str(sage_version),
            "pari": ".".join(str(part) for part in pari.version()),
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_kummer_classgroup_pressure.sage --check"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored Kummer/class-group pressure certificate differs")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"{PROTOCOL}|curves={len(document['curves'])}|status={document['status']}|"
        f"output={output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
