#!/usr/bin/env sage -python
"""Freeze generic-subgroup-only inputs for three exact rank-29 controls.

This boundary builder may read the complete public fixtures.  Its output
contains only the curve, the seventeen exactly transported generic sections,
and their exact generic height form.  No exceptional point coordinate is
written to the blind-search input.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
WGXLI = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
CURVE12 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
DIRECT12 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
OUTPUT = ROOT / "elliptic-curves/data/half_lattice_rank29_control_inputs_v1.json"
sys.path.insert(0, str(CAS))

import elkies_klagsbrun_rank29
import icarm_curve356


GENERIC_WORDS = (
    ((2, 1),), ((3, 1),), ((4, 1),), ((5, 1),), ((8, 1),),
    ((11, 1),), ((13, 1),), ((15, 1),), ((16, 1),), ((17, 1),),
    ((1, 1), (2, -1)), ((1, 1), (6, -1)), ((1, 1), (7, -1)),
    ((1, 1), (9, -1)), ((1, 1), (10, -1)), ((1, 1), (12, -1)),
    ((1, 1), (14, -1)),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def short_data(ainvs, points):
    a1, a2, a3, a4, a6 = (QQ(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    model = (QQ(0), QQ(0), QQ(0), -27 * c4, -54 * c6)
    short_points = tuple(
        (36 * QQ(point[0]) + 3 * b2, 108 * (2 * QQ(point[1]) + a1 * QQ(point[0]) + a3))
        for point in points
    )
    return model, short_points


def point_records(points):
    return [[str(point[0]), str(point[1])] for point in points]


def main() -> None:
    wgxli = json.loads(WGXLI.read_text())
    public = json.loads(PUBLIC.read_text())
    curve12 = json.loads(CURVE12.read_text())
    direct12 = json.loads(DIRECT12.read_text())

    word_gram = Matrix(ZZ, wgxli["generic_basis"]["height_gram"])
    word_matrix = Matrix(
        ZZ,
        17,
        17,
        lambda row, column: next(
            (coefficient for index, coefficient in GENERIC_WORDS[row] if index - 1 == column),
            0,
        ),
    )
    if abs(word_matrix.det()) != 1:
        raise ArithmeticError("074d9 generic word change is no longer unimodular")
    inverse = word_matrix.inverse()
    first17_gram = inverse * word_gram * inverse.transpose()
    if any(value not in ZZ for value in first17_gram.list()) or first17_gram.det() != 948:
        raise ArithmeticError("transformed 074d9 height form failed")
    first17_gram = Matrix(ZZ, first17_gram)

    model356 = tuple(QQ(value) for value in icarm_curve356.short_coefficients())
    generic356 = tuple(tuple(QQ(value) for value in point) for point in icarm_curve356.SHORT_POINTS[:17])

    source385 = next(record for record in public["records"] if record["id"] == 385)
    model385, all385 = short_data(source385["ainvs"], source385["points"])
    generic385 = all385[:17]

    model12 = tuple(QQ(value) for value in elkies_klagsbrun_rank29.short_weierstrass_coefficients())
    all12 = tuple(
        tuple(QQ(value) for value in point)
        for point in elkies_klagsbrun_rank29.published_short_points()
    )
    curve12_matrix = Matrix(ZZ, curve12["specialized_generic_subgroup"]["coordinate_matrix_rows_in_ordered_29_public_points"])
    if curve12_matrix.dimensions() != (29, 17):
        raise ArithmeticError("curve12 generic coordinate matrix changed")
    E12 = EllipticCurve(QQ, list(model12))
    public12 = tuple(E12(point) for point in all12)
    generic12 = []
    for column in range(17):
        point = E12(0)
        for coefficient, basis_point in zip(curve12_matrix.column(column), public12):
            point += int(coefficient) * basis_point
        generic12.append((point[0], point[1]))
    gram12 = Matrix(ZZ, direct12["sections"]["height_gram"])
    if gram12.det() != 948:
        raise ArithmeticError("alternate-Q80 height form changed")

    cases = (
        ("curve12-2024-rank29", model12, tuple(generic12), gram12, "alternate-Q80/norm12-orbit-11952"),
        ("curve356-rank29", model356, generic356, first17_gram, "published-R17/norm12-orbit-074d9"),
        ("curve385-rank29", model385, generic385, first17_gram, "published-R17/norm12-orbit-074d9"),
    )
    for label, model, points, unused_gram, unused_lineage in cases:
        E = EllipticCurve(QQ, list(model))
        if len(points) != 17 or any(E(point) == E(0) for point in points):
            raise ArithmeticError(f"{label}: invalid generic point inventory")

    payload = {
        "schema": "elliptic-curves.half-lattice-rank29-control-inputs.v1",
        "status": "FROZEN_GENERIC_SUBGROUP_ONLY_NO_EXCEPTIONAL_COORDINATES",
        "boundary": {
            "builder_read_full_public_fixtures": True,
            "output_contains_exceptional_point_coordinates": False,
            "blind_search_must_not_import_public_fixture_modules": True,
        },
        "cases": [
            {
                "label": label,
                "lineage": lineage,
                "short_model": [str(value) for value in model],
                "generic_points": point_records(points),
                "generic_height_gram": [[int(value) for value in row] for row in gram.rows()],
                "generic_height_gram_determinant": int(gram.det()),
            }
            for label, model, points, gram, lineage in cases
        ],
        "source_hashes": {
            str(WGXLI.relative_to(ROOT)): digest(WGXLI),
            str(PUBLIC.relative_to(ROOT)): digest(PUBLIC),
            str(CURVE12.relative_to(ROOT)): digest(CURVE12),
            str(DIRECT12.relative_to(ROOT)): digest(DIRECT12),
            str((CAS / "elkies_klagsbrun_rank29.py").relative_to(ROOT)): digest(CAS / "elkies_klagsbrun_rank29.py"),
            str((CAS / "icarm_curve356.py").relative_to(ROOT)): digest(CAS / "icarm_curve356.py"),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "claim_boundary": "The three generic subgroup transports are exact and already certified by the cited source artifacts; the output discloses no displayed complement point.",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"RANK29HALFINPUT|status=PASS|cases=3|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
