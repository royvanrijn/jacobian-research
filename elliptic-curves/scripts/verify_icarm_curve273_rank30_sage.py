#!/usr/bin/env sage-python
"""Independent Sage replay of the ICARM curve 273 rank-30 certificate."""

from __future__ import annotations

import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, QQ, matrix


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT / "cas"))

from icarm_curve273 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
    SHORT_POINTS,
    short_coefficients,
)


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve273_rank30_v1.json"
)


def sage_rational(value):
    return QQ(value.numerator) / QQ(value.denominator)


def reduce_rational(value, field):
    return field(value.numerator) / field(value.denominator)


def main() -> None:
    manifest = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    curve = EllipticCurve(
        QQ,
        [sage_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS],
    )
    rational_points = [
        curve(sage_rational(x_value), sage_rational(y_value))
        for x_value, y_value in POINTS
    ]
    assert len(rational_points) == 30
    assert list(curve.global_minimal_model().ainvs()) == list(curve.ainvs())
    assert str(curve.conductor()) == manifest["curve"]["conductor"]
    assert int(curve.root_number()) == manifest["curve"]["root_number"] == 1
    assert curve.torsion_subgroup().order() == manifest["curve"]["torsion_order"] == 1

    short_curve = EllipticCurve(
        QQ,
        [sage_rational(value) for value in short_coefficients()],
    )
    rows: list[list[int]] = []
    certificate_rows = manifest["independence_certificate"]["rows"]
    for expected in certificate_rows:
        prime = int(expected["prime"])
        field = GF(prime)
        reduction = EllipticCurve(field, short_curve.ainvs())
        group = reduction.abelian_group()
        invariants = tuple(int(value) for value in group.invariants())
        quotient_coordinates = [
            index for index, invariant in enumerate(invariants) if invariant % 2 == 0
        ]
        point_logs = []
        for x_value, y_value in SHORT_POINTS:
            reduced_point = reduction(
                reduce_rational(x_value, field),
                reduce_rational(y_value, field),
            )
            point_logs.append(group.discrete_log(reduced_point))
        for coordinate in quotient_coordinates:
            rows.append([int(log[coordinate]) % 2 for log in point_logs])
        assert reduction.cardinality() == expected["group_order"]
        assert len(quotient_coordinates) == expected["quotient_dimension"]

    binary_matrix = matrix(GF(2), rows)
    assert binary_matrix.ncols() == 30
    assert binary_matrix.rank() == 30

    witness_prime = manifest["independence_certificate"][
        "no_rational_2_torsion_witness_prime"
    ]
    witness_field = GF(witness_prime)
    division_polynomial = short_curve.division_polynomial(2).change_ring(
        witness_field
    )
    assert not division_polynomial.roots()
    print(
        "PASS independent Sage replay: 30 exact points, global minimal model, "
        "trivial torsion, and full-rank mod-2 reduction matrix"
    )


if __name__ == "__main__":
    main()
