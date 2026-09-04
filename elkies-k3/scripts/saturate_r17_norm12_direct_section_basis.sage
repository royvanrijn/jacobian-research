#!/usr/bin/env sage-python
"""Saturate an explicit norm-12 child section basis by exact halving."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[index]) for index in range(poly.degree() + 1)]


def rational_function_record(value):
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def parse_function(record, ring):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return ring.fraction_field()(numerator / denominator)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def binary_scalar_additions(coefficient: int) -> int:
    value = abs(int(coefficient))
    if value <= 1:
        return 0
    return value.bit_length() - 1 + value.bit_count() - 1


def combination_score(coefficients) -> tuple:
    values = tuple(map(int, coefficients))
    support = tuple(value for value in values if value)
    return (
        sum(binary_scalar_additions(value) for value in support)
        + max(0, len(support) - 1),
        len(support),
        max(map(abs, support)),
        sum(map(abs, support)),
        values,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    suffix = args.source_label.rsplit("-", 1)[1]
    input_path = args.input or ROOT / (
        f"artifacts/generated-results/elkies-k3-r17-norm12-orbit{suffix}-direct-fibration-v1.json"
    )
    output_path = args.output or ROOT / (
        f"artifacts/generated-results/elkies-k3-r17-norm12-orbit{suffix}-direct-fibration-saturated-v1.json"
    )
    source = json.loads(input_path.read_text())
    section_data = source["sections"]
    if (
        source["status"]
        != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_RANK17_SUBLATTICE"
        or section_data["status"] != "PASS_EXACT_RANK17_FINITE_INDEX_SUBLATTICE"
        or section_data["rank"] != 17
    ):
        raise ArithmeticError("input is not an exact finite-index rank-17 child marking")

    Ru = PolynomialRing(QQ, "u")
    Ku = Ru.fraction_field()
    A = Ru([QQ(value) for value in source["weierstrass_model"]["A_coefficients_low_to_high"]])
    B = Ru([QQ(value) for value in source["weierstrass_model"]["B_coefficients_low_to_high"]])
    curve = EllipticCurve(Ku, [A, B])
    records = section_data["records"]
    points = [
        curve(parse_function(record["X"], Ru), parse_function(record["Y"], Ru))
        for record in records
    ]
    coordinates = matrix(ZZ, section_data["coordinate_matrix_in_compiled_frame"])
    frame = matrix(ZZ, source["frame_certificate"]["frame_gram"])
    if len(points) != 17 or coordinates.nrows() != 17 or coordinates.ncols() != 17:
        raise ArithmeticError("input section basis has the wrong size")
    if abs(coordinates.det()) != int(section_data["index_in_saturated_mw_lattice"]):
        raise ArithmeticError("declared section index and coordinate determinant differ")

    inverse = coordinates.inverse()
    saturation_candidates = []
    short_vectors = matrix(ZZ, pari(frame).qfminim(4)[2])
    for column in short_vectors.columns():
        for candidate in (vector(ZZ, column), -vector(ZZ, column)):
            rational_coefficients = candidate * inverse
            doubled_coefficients = 2 * rational_coefficients
            if (
                any(value not in ZZ for value in rational_coefficients)
                and all(value in ZZ for value in doubled_coefficients)
            ):
                saturation_candidates.append(
                    (
                        int(candidate * frame * candidate),
                        combination_score(doubled_coefficients),
                        tuple(map(int, candidate)),
                        vector(ZZ, doubled_coefficients),
                    )
                )
    if not saturation_candidates:
        raise ArithmeticError("no order-two saturation coset found in frame coordinates")
    _, _, target_entries, coefficients = min(saturation_candidates)
    target = vector(ZZ, target_entries)

    doubled_point = sum(
        (coefficient * point for coefficient, point in zip(coefficients, points)),
        curve(0),
    )
    if doubled_point.is_zero() or doubled_point[1] == 0:
        raise ArithmeticError("selected double is unsuitable for rational halving")
    x_double, y_double = Ku(doubled_point[0]), Ku(doubled_point[1])
    Rz = PolynomialRing(Ku, "z")
    z = Rz.gen()
    duplication_quartic = (
        z**4
        - 4 * x_double * z**3
        - 2 * A * z**2
        - (4 * A * x_double + 8 * B) * z
        + A**2
        - 4 * B * x_double
    )
    linear_roots = []
    for factor, multiplicity in duplication_quartic.factor():
        if factor.degree() == 1:
            linear_roots.extend([-factor[0] / factor[1]] * multiplicity)
    if not linear_roots:
        raise ArithmeticError("duplication quartic has no rational linear factor")

    half_point = None
    for x_half in linear_roots:
        f_half = x_half**3 + A * x_half + B
        y_half = (
            (3 * x_half**2 + A) * (x_half - x_double) - 2 * f_half
        ) / (2 * y_double)
        candidate = curve(Ku(x_half), Ku(y_half))
        if 2 * candidate == doubled_point:
            half_point = candidate
            break
    if half_point is None:
        raise ArithmeticError("rational duplication root did not produce the selected half")

    saturation_record = {
        "basis_index": 17,
        "source": "exact rational halving of an index-two child-section combination",
        "doubled_section_coefficients_in_input_basis": list(map(int, coefficients)),
        "new_frame_coordinates": list(map(int, target)),
        "new_height": int(target * frame * target),
        "X": rational_function_record(Ku(half_point[0])),
        "Y": rational_function_record(Ku(half_point[1])),
        "equation_verified": True,
        "doubling_verified": True,
    }
    augmented_records = [copy.deepcopy(record) for record in records] + [saturation_record]
    augmented_rows = [vector(ZZ, row) for row in coordinates.rows()] + [target]
    selected_indices = None
    saturated_coordinates = None
    for removed in range(17):
        indices = [index for index in range(18) if index != removed]
        trial = matrix(ZZ, [augmented_rows[index] for index in indices])
        if abs(trial.det()) == 1:
            selected_indices = indices
            saturated_coordinates = trial
            break
    if selected_indices is None:
        raise ArithmeticError("halved section did not saturate the coordinate lattice")

    selected_records = [augmented_records[index] for index in selected_indices]
    for basis_index, record in enumerate(selected_records):
        record["basis_index"] = basis_index
    height_gram = saturated_coordinates * frame * saturated_coordinates.transpose()
    if height_gram.det() != 948:
        raise ArithmeticError("saturated section height determinant changed")

    result = copy.deepcopy(source)
    result["schema"] = f"elkies-k3.r17-norm12-orbit{suffix}-direct-fibration-saturated.v1"
    result["status"] = "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS"
    result["sections"] = {
        "status": "PASS_EXACT_SATURATED_RANK17_BASIS",
        "rank": 17,
        "basis_source_profile": {
            "retained_input_sections": 16,
            "exact_rational_halves": 1,
        },
        "saturation": {
            "input_index": int(abs(coordinates.det())),
            "target_frame_coordinates": list(map(int, target)),
            "doubled_section_coefficients_in_input_basis": list(map(int, coefficients)),
            "duplication_quartic_linear_factor_count": len(linear_roots),
            "removed_input_basis_index": next(
                index for index in range(17) if index not in selected_indices
            ),
        },
        "coordinate_matrix_in_compiled_frame": rows(saturated_coordinates),
        "coordinate_matrix_determinant": int(saturated_coordinates.det()),
        "height_gram": rows(height_gram),
        "height_gram_determinant": int(height_gram.det()),
        "roots_of_norm_two": 0,
        "records": selected_records,
    }
    result["inputs"] = {relative(input_path): digest(input_path)}
    result["reproducing_command"] = (
        "sage -python elkies-k3/scripts/saturate_r17_norm12_direct_section_basis.sage "
        f"--source-label {args.source_label} --input {relative(input_path)} "
        f"--output {relative(output_path)}"
    )
    result["proof_boundary"] = (
        "The displayed rational half doubles exactly to the declared combination of "
        "input sections. Its abstract frame coordinate enlarges the input lattice, and "
        "the resulting seventeen coordinate rows have determinant one. Hence the new "
        "equation section basis is saturated."
    )
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise ArithmeticError("stored saturated direct marking differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "R17NORM12SATURATE"
        f"|label={args.source_label}|input_index={abs(coordinates.det())}"
        f"|output_index={abs(saturated_coordinates.det())}"
        f"|half_height={target * frame * target}|output={relative(output_path)}"
    )


if __name__ == "__main__":
    main()
