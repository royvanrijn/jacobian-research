#!/usr/bin/env sage-python
"""Independently replay the exact curve-398 two-parent subgroup collision."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path

from sage.all import EllipticCurve, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_two_parent_collision_v1.json"
PUBLIC = ROOT / "elliptic-curves/cas/icarm_curve398.py"


def load(name: str, path: Path):
    return SourceFileLoader(name, str(path)).load_module()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()

    document = json.loads(args.artifact.read_text())
    public = load("curve398_two_parent_public", PUBLIC)
    if document["status"] != "PASS_EXACT_TWO_PARENT_COLLISION":
        raise ArithmeticError("two-parent artifact is not passing")
    for path, expected in document["inputs"].items():
        if digest(ROOT / path) != expected:
            raise ArithmeticError(f"two-parent input hash changed: {path}")
    generation = document["generation"]
    if digest(ROOT / generation["script"]) != generation["script_sha256"]:
        raise ArithmeticError("two-parent compiler hash changed")

    curve = EllipticCurve(QQ, list(public.GENERAL_WEIERSTRASS_COEFFICIENTS))
    public_points = tuple(curve(QQ(str(x)), QQ(str(y))) for x, y in public.POINTS)
    if len(public_points) != 30:
        raise ArithmeticError("public M30 basis length changed")

    first = matrix(ZZ, document["first_parent"]["matrix_16_by_30_rows"])
    second_record = document["second_parent"]
    second = matrix(ZZ, second_record["public_rank30_embedding"]["matrix_16_by_30_rows"])
    second_points = tuple(
        curve(QQ(record["specialized_public_point"]["x"]), QQ(record["specialized_public_point"]["y"]))
        for record in second_record["generic_mw16"]["records"]
    )
    for point, row in zip(second_points, second.rows()):
        replay = sum((int(coefficient) * basis for coefficient, basis in zip(row, public_points) if coefficient), curve(0))
        if replay != point:
            raise ArithmeticError("a second-parent point failed independent Sage group-law replay")

    combined = first.stack(second)
    collision = document["collision"]
    ranks = [int(first.rank()), int(second.rank()), int(combined.rank())]
    if ranks != [collision["rank_G1"], collision["rank_G2"], collision["rank_sum"]]:
        raise ArithmeticError("two-parent ranks changed")
    if collision["rank_intersection"] != ranks[0] + ranks[1] - ranks[2]:
        raise ArithmeticError("intersection dimension formula failed")

    second_in_first = matrix(ZZ, collision["G2_basis_rows_in_G1_basis"])
    first_in_second = matrix(ZZ, collision["G1_basis_rows_in_G2_basis"])
    if (
        second_in_first * first != second
        or first_in_second * second != first
        or second_in_first * first_in_second != 1
        or first_in_second * second_in_first != 1
        or abs(second_in_first.det()) != 1
        or abs(first_in_second.det()) != 1
    ):
        raise ArithmeticError("integral subgroup equality failed")

    smith = combined.smith_form()[0]
    diagonal = [abs(int(smith[index, index])) for index in range(ranks[2])]
    if diagonal != collision["smith_diagonal_nonzero"]:
        raise ArithmeticError("two-parent Smith diagonal changed")
    if collision["quotient_free_rank"] != 30 - ranks[2]:
        raise ArithmeticError("Smith quotient free rank changed")
    if collision["quotient_torsion_invariant_factors_nontrivial"] != [value for value in diagonal if value > 1]:
        raise ArithmeticError("Smith quotient torsion factors changed")
    expected_index = "infinite" if ranks[2] < 30 else str(abs(int(smith[:30, :30].det())))
    if collision["smith_index_in_public_M30"] != expected_index:
        raise ArithmeticError("Smith index changed")

    generic_gram = matrix(QQ, second_record["generic_mw16"]["height_gram"])
    if generic_gram.rank() != 16 or generic_gram.det() != 474:
        raise ArithmeticError("second generic MW16 certificate changed")
    print(
        "CURVE398COLLISIONVERIFY|G1=G2|intersection=16|sum=16|"
        "M30_quotient=Z^14|smith_index=infinite|sage_group_law=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
