#!/usr/bin/env sage
"""Certify the complete rational CM locus on the det-1236 marked curve.

For

    C_1236 = X_0^6(103)/<w_618>,

Gonzalez--Rotger Corollary 5.14 bounds the class number of an order with a
rational CM image by two.  This replay checks the complete class-number-one
and class-number-two order lists, the local optimal-embedding factors, the
residue-field alternatives, and the CM multiplicities.  Exactly the orders
of discriminants -3, -43, and -67 survive, contributing respectively two,
four, and four rational points on C_1236.

The class-number lists and the Gonzalez--Rotger field formula are named
external theorem inputs.  All arithmetic substitutions and counts are exact.

Replay:
    sage elkies-k3/scripts/certify_det1236_rational_cm_locus.sage
    sage elkies-k3/scripts/certify_det1236_rational_cm_locus.sage check
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path.cwd().resolve()
if not (ROOT / "elkies-k3").is_dir():
    raise RuntimeError("run this certificate from the repository root")
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-det1236-rational-cm-locus-v1.json"
)

D = ZZ(6)
N = ZZ(103)
M = ZZ(618)

# Complete Watkins lists of imaginary quadratic orders of class number one
# and two, represented by (conductor, fundamental discriminant), as pinned in
# Padurariu--Saia GenusAtMost2 commit 6cc368fe37aa67187783118f18d149b2b1fd6230.
CLASS_NUMBER_ONE_ORDERS = [
    (1, -3), (2, -3), (3, -3), (1, -4), (2, -4),
    (1, -7), (2, -7), (1, -8), (1, -11), (1, -19),
    (1, -43), (1, -67), (1, -163),
]
CLASS_NUMBER_TWO_ORDERS = [
    (4, -3), (5, -3), (7, -3), (3, -4), (4, -4),
    (5, -4), (4, -7), (2, -8), (3, -8), (3, -11),
    (1, -15), (2, -15), (1, -20), (1, -24), (1, -35),
    (1, -40), (1, -51), (1, -52), (1, -88), (1, -91),
    (1, -115), (1, -123), (1, -148), (1, -187),
    (1, -232), (1, -235), (1, -267), (1, -403),
    (1, -427),
]


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order_record(class_number, conductor, field_discriminant):
    conductor = ZZ(conductor)
    field_discriminant = ZZ(field_discriminant)
    order_discriminant = conductor**2*field_discriminant
    assert len(BinaryQF_reduced_representatives(
        order_discriminant, primitive_only=True
    )) == class_number

    def eichler_symbol(prime):
        return ZZ(1) if conductor % prime == 0 else ZZ(
            kronecker(field_discriminant, prime)
        )

    local_factors = {
        "2": int(1-eichler_symbol(2)),
        "3": int(1-eichler_symbol(3)),
        "103": int(1+eichler_symbol(103)),
    }
    top_count = ZZ(class_number)*prod(local_factors.values())

    d_r = prod(
        prime for prime in (2, 3)
        if conductor % prime and kronecker(field_discriminant, prime) == -1
    )
    n_r = prod(
        prime for prime in (103,)
        if conductor % prime == 0
        or kronecker(field_discriminant, prime) == 1
    )
    nstar_r = prod(
        prime for prime in (103,)
        if conductor % prime and kronecker(field_discriminant, prime) == 1
    )
    m_r = gcd(M, abs(order_discriminant)//gcd(N, conductor))
    quotient = M//m_r
    d_r_nstar_r = ZZ(d_r*nstar_r)

    # Corollary 5.14: for h=1 the criterion below is exact.  For h=2 a
    # rational field can occur only in the two-involution case.  No nonempty
    # class-number-two row reaches that case, so no Artin-class ambiguity is
    # left in this instance.
    if class_number == 1:
        rational_image = bool(
            top_count
            and (d_r_nstar_r == 1 or quotient == d_r_nstar_r)
        )
    else:
        rational_image = bool(
            top_count and d_r_nstar_r == 1 and quotient == 1
        )

    marked_count = 0
    if rational_image:
        # The fixed locus of w_618 has order discriminant -2472.  The three
        # surviving orders are different, so w_618 acts freely on each locus.
        assert order_discriminant != -2472
        assert top_count % 2 == 0
        marked_count = top_count//2

    return {
        "class_number": int(class_number),
        "conductor": int(conductor),
        "field_discriminant": int(field_discriminant),
        "order_discriminant": int(order_discriminant),
        "local_embedding_factors_p2_p3_p103": local_factors,
        "top_curve_cm_points": int(top_count),
        "D_R": int(d_r),
        "N_R": int(n_r),
        "N_star_R": int(nstar_r),
        "m_R": int(m_r),
        "m_over_m_R": int(quotient),
        "rational_image_on_w618_quotient": rational_image,
        "marked_curve_rational_cm_points": int(marked_count),
    }


def build_payload():
    rows = []
    for class_number, orders in (
        (1, CLASS_NUMBER_ONE_ORDERS),
        (2, CLASS_NUMBER_TWO_ORDERS),
    ):
        rows.extend(
            order_record(class_number, conductor, field_discriminant)
            for conductor, field_discriminant in orders
        )

    rational_rows = [
        row for row in rows
        if row["rational_image_on_w618_quotient"]
    ]
    assert [row["order_discriminant"] for row in rational_rows] == [
        -3, -43, -67
    ]
    assert [row["marked_curve_rational_cm_points"] for row in rational_rows] == [
        2, 4, 4
    ]
    assert not any(
        row["class_number"] == 2
        and row["top_curve_cm_points"]
        and row["D_R"]*row["N_star_R"] == 1
        and row["m_over_m_R"] == 1
        for row in rows
    )
    total = sum(row["marked_curve_rational_cm_points"] for row in rational_rows)
    assert total == 10

    return {
        "schema": "elkies-k3.det1236-rational-cm-locus.v1",
        "status": "PASS_DET1236_COMPLETE_RATIONAL_CM_LOCUS",
        "marked_curve": "X_0^6(103)/<w_618>",
        "orders_checked": len(rows),
        "class_numbers_checked": [1, 2],
        "rational_cm_rows": rational_rows,
        "rational_cm_order_discriminants": [-3, -43, -67],
        "rational_cm_point_count": int(total),
        "higher_class_number_exclusion": (
            "Corollary 5.14 fixes the ring class field by at most one "
            "involution, except in one case by at most two involutions. "
            "The norm-m_R ideal class has order at most two. Hence a "
            "rational image forces h(R)<=2; the checked lists are complete."
        ),
        "cover_candidate_consequence": (
            "Any identified model of C_1236 with exactly two fixed and eight "
            "non-fixed rational points has no rational non-CM points: all ten "
            "are exhausted by the certified CM locus. A candidate twist with "
            "only the two fixed rational points cannot be C_1236."
        ),
        "external_inputs": [
            "Gonzalez--Rotger Corollary 5.14 and Proposition 5.6",
            "complete Watkins class-number-one/two quadratic-order lists",
        ],
        "reproduce": (
            "sage elkies-k3/scripts/"
            "certify_det1236_rational_cm_locus.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", choices=("write", "check"), default="write")
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n"
    if args.mode == "check":
        if not OUTPUT.is_file():
            raise FileNotFoundError(OUTPUT)
        if OUTPUT.read_text() != rendered:
            raise AssertionError("generated artifact changed: %s" % OUTPUT)
        print(json.dumps({
            "status": "PASS_DET1236_RATIONAL_CM_LOCUS_CHECK",
            "output": relative(OUTPUT),
            "sha256": digest(OUTPUT),
            "rational_cm_points": payload["rational_cm_point_count"],
        }, sort_keys=True))
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(rendered)
    print(json.dumps({
        "status": "WROTE_DET1236_RATIONAL_CM_LOCUS_CERTIFICATE",
        "output": relative(OUTPUT),
        "sha256": digest(OUTPUT),
        "rational_cm_points": payload["rational_cm_point_count"],
    }, sort_keys=True))


main()
