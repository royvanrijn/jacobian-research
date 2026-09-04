#!/usr/bin/env sage-python
"""Descend the four record-specific 074d9 bisections to twist sections.

For a lifted bisection point ``P`` over ``z^2=q(u)``, its trace
``T=P+sigma(P)`` lies in the native rank-17 group.  The anti-invariant point

    R = P-sigma(P) = 2*P-T

descends to the quadratic twist

    Y^2 = X^3 + q^2*A*X + q^3*B,

with ``X=q*x(R)`` and ``Y=q^2*coefficient_z(y(R))``.  This script reconstructs
the chord data from the canonical trace words, verifies all four identities
over QQ(u), and records the exact specialization class ``2*[P]`` in each
displayed record-fibre quotient.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
COVERS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
PRIORITY = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv"
CHORD_HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
CROSS_HELPER = ROOT / "elkies-k3/scripts/certify_r17_074d9_cross_fibre_bisection_transfer.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-record-twist-sections-v1.json"
)
RECORD_LABELS = (
    "074d9-orbit-04b07",
    "074d9-orbit-11a44",
    "074d9-orbit-11279",
    "074d9-orbit-080fa",
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def polynomial_text(polynomial) -> list[str]:
    return [str(value) for value in polynomial.list()]


def require_polynomial(value, ring, label):
    value = value.numerator() / value.denominator()
    if value.denominator().degree() != 0:
        raise ArithmeticError(f"{label} is not a polynomial")
    return ring(value)


def build_payload():
    lineage = json.loads(LINEAGE.read_text())
    covers = json.loads(COVERS.read_text())
    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("unexpected lineage status")
    if covers.get("status") != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER":
        raise ValueError("unexpected cover-certificate status")

    cross = runpy.run_path(str(CROSS_HELPER), run_name="r17_074d9_cross_helper")
    chord = runpy.run_path(str(CHORD_HELPER), run_name="r17_074d9_chord_helper")
    priority_rows = cross["parse_priority_rows"](PRIORITY)
    unused_grams = cross["representative_words"](lineage, priority_rows)
    words = unused_grams[3]
    ring, field, A, B, discriminant, curve, unused_basis, multiples = cross[
        "build_exact_context"
    ](lineage)

    selected = {
        row["label"]: (int(fibre["curve_id"]), QQ(fibre["parameter"]), row)
        for fibre in covers["fibres"]
        for row in fibre["records"]
        if row["label"] in RECORD_LABELS
    }
    if tuple(label for label in RECORD_LABELS if label not in selected):
        raise ArithmeticError("the cover certificate lost a record-specific twist")

    records = []
    for label in RECORD_LABELS:
        curve_id, parameter, cover = selected[label]
        index = int(cover["priority_index_zero_based"])
        trace = cross["trace_from_word"](words[index], curve, multiples)
        data = cross["exact_chord_data"](
            trace, A, B, discriminant, ring, field, chord
        )
        q = ring(data["q"])

        extension_ring = PolynomialRing(field, "z")
        z_polynomial = extension_ring.gen()
        extension = field.extension(z_polynomial**2 - q, "z")
        z = extension.gen()
        extended_curve = EllipticCurve(extension, [A, B])
        P = extended_curve(
            extension(data["x0"]) + extension(data["x1"]) * z,
            extension(data["y0"]) + extension(data["y1"]) * z,
        )
        sigma_P = extended_curve(
            extension(data["x0"]) - extension(data["x1"]) * z,
            extension(data["y0"]) - extension(data["y1"]) * z,
        )
        if P + sigma_P != extended_curve(extension(trace[0]), extension(trace[1])):
            raise ArithmeticError(f"{label}: chord trace identity failed")
        anti = P - sigma_P
        x_coefficients = anti[0].list()
        y_coefficients = anti[1].list()
        if len(x_coefficients) > 1 and any(x_coefficients[1:]):
            raise ArithmeticError(f"{label}: anti-invariant x is not invariant")
        if not y_coefficients or y_coefficients[0] != 0 or any(y_coefficients[2:]):
            raise ArithmeticError(f"{label}: anti-invariant y has the wrong character")
        X = require_polynomial(field(q) * field(x_coefficients[0]), ring, f"{label} X")
        Y = require_polynomial(field(q) ** 2 * field(y_coefficients[1]), ring, f"{label} Y")
        if Y**2 != X**3 + A * q**2 * X + B * q**3:
            raise ArithmeticError(f"{label}: twist equation failed")
        if X.degree() > 6 or Y.degree() > 9:
            raise ArithmeticError(f"{label}: section is not P.O=0")

        q_value = QQ(q(parameter))
        square_root = QQ(cover["canonical_square_root"])
        if square_root**2 != q_value:
            raise ArithmeticError(f"{label}: stored split square root changed")
        twist_curve = EllipticCurve(
            QQ,
            [QQ(A(parameter)) * q_value**2, QQ(B(parameter)) * q_value**3],
        )
        twist_point = twist_curve(QQ(X(parameter)), QQ(Y(parameter)))
        fibre_curve = EllipticCurve(QQ, [QQ(A(parameter)), QQ(B(parameter))])
        image = fibre_curve(
            twist_point[0] / q_value,
            twist_point[1] / (q_value * square_root),
        )
        positive = fibre_curve(
            QQ(cover["positive_chart_point"][0]),
            QQ(cover["positive_chart_point"][1]),
        )
        trace_at_fibre = fibre_curve(QQ(trace[0](parameter)), QQ(trace[1](parameter)))
        if image != 2 * positive - trace_at_fibre:
            raise ArithmeticError(f"{label}: specialized anti-invariant image changed")

        quotient_class = list(
            map(int, cover["exact_displayed_free_quotient_class"]["coordinates"])
        )
        image_class = [2 * value for value in quotient_class]
        records.append(
            {
                "label": label,
                "curve_id": curve_id,
                "lattice_orbit_mask": int(cover["lattice_orbit_mask"]),
                "priority_index_zero_based": index,
                "construction_chart": data["construction_chart"],
                "twist_model": "Y^2=X^3+A*q^2*X+B*q^3",
                "q_coefficients_low_to_high": polynomial_text(q),
                "X_coefficients_low_to_high": polynomial_text(X),
                "Y_coefficients_low_to_high": polynomial_text(Y),
                "degrees_X_Y": [int(X.degree()), int(Y.degree())],
                "P_dot_O": 0,
                "height_on_twist_surface": 6,
                "height_after_quadratic_base_change": 12,
                "record_parameter": str(parameter),
                "record_twist_point": [str(twist_point[0]), str(twist_point[1])],
                "record_fibre_image": [str(image[0]), str(image[1])],
                "displayed_quotient_basis": cover[
                    "exact_displayed_free_quotient_class"
                ]["basis"],
                "bisection_branch_quotient_class": quotient_class,
                "twist_section_image_quotient_class": image_class,
                "specialization_identity": "sp(R)=2*positive_branch-sp(trace)",
                "exact_identities_verified": True,
            }
        )

    return {
        "schema": "elkies-k3.r17-074d9-record-twist-sections.v1",
        "status": "PASS_EXACT_FOUR_RECORD_TWIST_SECTIONS_AND_SPECIALIZATIONS",
        "claim": (
            "The four rigid bisection points descend to exact P.O=0 sections on "
            "their singleton twists, and their record-fibre images are computed "
            "exactly in the displayed exceptional quotients."
        ),
        "proof_boundary": (
            "This proves one non-torsion direction and its specialization for each "
            "twist. It does not prove that the section is a full Mordell--Weil basis "
            "or give a rank upper bound."
        ),
        "records": records,
        "inputs": {
            relative(path): digest(path)
            for path in (LINEAGE, COVERS, PRIORITY, CHORD_HELPER, CROSS_HELPER)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale 074d9 record-twist section certificate")
        terminal = "PASS"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "R17074D9TWISTSECTIONS|records=4|degrees="
        + ",".join("/".join(map(str, row["degrees_X_Y"])) for row in payload["records"])
        + f"|status={terminal}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
