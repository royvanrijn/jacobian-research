#!/usr/bin/env sage-python
"""Construct and specialize the complete norm-26 R17 trisection shell.

The complete degree-three coset census has 320 rational translation cosets
whose minimum norm is 26.  Up to inversion there are 160.  For a minimum
representative ``w`` let ``tau`` be the corresponding published R17 section.
The unique member of

    H^0(X, O_X(4 O + 20 F)) = <1, x, y, x^2>

through ``tau`` has coefficient-degree bounds ``20,16,14,12``.  Removing the
known section from its degree-four intersection with the generic elliptic
fibre produces an exact residual cubic.  This script constructs those cubics,
proves their generic irreducibility, and factors them exactly at the four
published rank-25--28 controls.

This is an equation and specialization certificate.  The separate quotient
analyser maps the split points into the known exceptional quotients.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    ZZ,
    QuadraticForm,
    lcm,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
DEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-degree3-deep-cosets-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-deep-trisections-v1.json"

CONTROLS = (
    ("rank_at_least_28", -9529, 5471, 28),
    ("rank_at_least_27", 2456, 135, 27),
    ("rank_at_least_26", -308, 251, 26),
    ("rank_at_least_25", -2, 377, 25),
)

BOUNDS = (20, 16, 14, 12)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def qtext(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def poly_record(polynomial) -> list[str]:
    polynomial = polynomial.parent()(polynomial)
    if not polynomial:
        return ["0"]
    return [qtext(polynomial[index]) for index in range(polynomial.degree() + 1)]


def rational_function_record(function) -> dict:
    return {
        "numerator_coefficients_low_to_high": poly_record(function.numerator()),
        "denominator_coefficients_low_to_high": poly_record(function.denominator()),
    }


def reconstruct_basis(R, A, B, section_data):
    points = []
    for expected_index, record in enumerate(section_data["sections"]):
        assert int(record["basis_index"]) == expected_index
        x_coordinate = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
        if expected_index == 0:
            y_coordinate = R(
                [QQ(value) for value in record["y_coefficients_low_to_high"]]
            )
        else:
            chord = record["chord"]
            reference_x, reference_y = points[int(chord["reference_basis_index"])]
            slope = R(
                [QQ(value) for value in chord["slope_coefficients_low_to_high"]]
            )
            y_coordinate = reference_y + slope * (x_coordinate - reference_x)
        assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
        points.append((x_coordinate, y_coordinate))
    return points


def quadratic_form_from_gram(G):
    coefficients = []
    for row in range(G.nrows()):
        for column in range(row, G.ncols()):
            coefficients.append(
                G[row, row] // 2 if row == column else G[row, column]
            )
    return QuadraticForm(ZZ, G.nrows(), coefficients)


def published_gram_and_foundry_map(foundry, target, pinned):
    basis_change = matrix(ZZ, target["pinned_identification"]["basis_change_matrix"])
    orientation = target["pinned_identification"]["gram_identity_orientation"]
    if orientation != "M^T*Gpub*M=Gpinned":
        raise ArithmeticError("the pinned/published Gram convention changed")
    published = matrix(
        ZZ, basis_change.transpose().inverse() * pinned * basis_change.inverse()
    )
    frame = next(
        item
        for ns_class in foundry["ns_classes"]
        for item in ns_class["frames"]
        if item["frame_id"] == "NS0001-F001"
    )
    frame_gram = matrix(ZZ, frame["gram"])
    reduced_gram = matrix(ZZ, frame["reduced_gram"])
    reduced_to_frame = matrix(ZZ, frame["reduced_basis_columns_in_frame_basis"])
    assert reduced_to_frame.transpose() * frame_gram * reduced_to_frame == reduced_gram
    isometry = quadratic_form_from_gram(published).is_globally_equivalent_to(
        quadratic_form_from_gram(frame_gram), return_matrix=True
    )
    if isometry is False:
        raise ArithmeticError("the published R17 control is no longer isometric")
    isometry = matrix(ZZ, isometry)
    if isometry.transpose() * published * isometry == frame_gram:
        frame_to_published = isometry
    elif isometry * published * isometry.transpose() == frame_gram:
        frame_to_published = isometry.transpose()
    else:
        raise ArithmeticError("qfisom returned an unrecognized orientation")
    assert frame_to_published.transpose() * published * frame_to_published == frame_gram
    return published, frame_gram, frame_to_published


def rational_bit_cost(value: str) -> int:
    numerator, separator, denominator = str(value).partition("/")
    cost = abs(int(numerator)).bit_length()
    if separator:
        cost += abs(int(denominator)).bit_length()
    return cost


def section_input_costs(section_data):
    direct_costs = []
    closures = []
    for index, record in enumerate(section_data["sections"]):
        fields = [record["x_coefficients_low_to_high"]]
        closure = {index}
        if index == 0:
            fields.append(record["y_coefficients_low_to_high"])
        else:
            chord = record["chord"]
            reference = int(chord["reference_basis_index"])
            closure.update(closures[reference])
            fields.append(chord["slope_coefficients_low_to_high"])
        direct_costs.append(
            sum(rational_bit_cost(value) for field in fields for value in field)
        )
        closures.append(closure)
    return direct_costs, closures


def scalar_additions(coefficient: int) -> int:
    value = abs(int(coefficient))
    if value <= 1:
        return 0
    return value.bit_length() - 1 + value.bit_count() - 1


def equation_score(published_vector, direct_costs, closures):
    support = [index for index, coefficient in enumerate(published_vector) if coefficient]
    dependency_closure = set()
    for index in support:
        dependency_closure.update(closures[index])
    additions = sum(scalar_additions(published_vector[index]) for index in support)
    additions += max(0, len(support) - 1)
    record = {
        "group_addition_upper_bound": additions,
        "support_count": len(support),
        "support_one_based": [index + 1 for index in support],
        "dependency_count": len(dependency_closure),
        "dependency_closure_one_based": [index + 1 for index in sorted(dependency_closure)],
        "coordinate_input_bits": sum(direct_costs[index] for index in dependency_closure),
        "maximum_absolute_coefficient": max(abs(value) for value in published_vector),
        "coefficient_l1": sum(abs(value) for value in published_vector),
    }
    key = (
        record["group_addition_upper_bound"],
        record["support_count"],
        record["dependency_count"],
        record["coordinate_input_bits"],
        record["maximum_absolute_coefficient"],
        record["coefficient_l1"],
        tuple(published_vector),
    )
    return key, record


def canonical_orientation(foundry_vector, frame_to_published, direct_costs, closures):
    published = frame_to_published * vector(ZZ, foundry_vector)
    choices = []
    for sign in (1, -1):
        candidate = tuple(sign * int(value) for value in published)
        key, score = equation_score(candidate, direct_costs, closures)
        choices.append((key, candidate, score, sign))
    return min(choices, key=lambda item: item[0])


def reciprocal_with_bound(polynomial, bound, R):
    polynomial = R(polynomial)
    if polynomial.degree() > bound:
        raise ArithmeticError("reciprocal degree exceeds its geometric bound")
    t = R.gen()
    return R(
        sum(
            polynomial[index] * t ** (bound - index)
            for index in range(polynomial.degree() + 1)
        )
    )


def invert_rational(function, weight, R, K):
    function = K(function)
    numerator = R(function.numerator())
    denominator = R(function.denominator())
    reversed_numerator = reciprocal_with_bound(numerator, numerator.degree(), R)
    reversed_denominator = reciprocal_with_bound(denominator, denominator.degree(), R)
    exponent = weight - numerator.degree() + denominator.degree()
    return K(t_power(R, exponent) * reversed_numerator / reversed_denominator)


def t_power(R, exponent):
    t = R.gen()
    if exponent >= 0:
        return t**exponent
    return R.fraction_field()(1) / t ** (-exponent)


def interpolation_relation(X, Y, A, B, R, K):
    denominator = R(X.denominator())
    if not denominator.is_square():
        raise ArithmeticError("trace x denominator is not a square")
    h = R(denominator.sqrt())
    h /= h.leading_coefficient()
    assert X.denominator() == h**2 and Y.denominator() == h**3
    if h.degree() != 11:
        raise ArithmeticError(f"trace denominator degree {h.degree()} is not 11")
    Nx = R(X * h**2)
    Ny = R(Y * h**3)
    terms = (h**4, Nx * h**2, Ny * h, Nx**2)
    assert all(term.degree() <= 64 for term in terms)
    columns = []
    for term, bound in zip(terms, BOUNDS):
        for power in range(bound + 1):
            polynomial = term * R.gen()**power
            columns.append([polynomial[index] for index in range(65)])
    interpolation = matrix(QQ, 65, len(columns), lambda row, column: columns[column][row])
    if interpolation.ncols() != 66 or interpolation.rank() != 65:
        raise ArithmeticError("the deep trisection interpolation is not uniquely one-dimensional")
    kernel = interpolation.right_kernel().basis()
    if len(kernel) != 1:
        raise ArithmeticError("unexpected interpolation kernel dimension")
    raw = kernel[0]
    common_denominator = lcm([value.denominator() for value in raw])
    primitive = [ZZ(common_denominator * value) for value in raw]
    common = 0
    for value in primitive:
        common = gcd(common, abs(int(value)))
    primitive = [value // common for value in primitive]
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = [-value for value in primitive]
    coefficients = []
    offset = 0
    for bound in BOUNDS:
        coefficients.append(R(primitive[offset : offset + bound + 1]))
        offset += bound + 1
    f0, f1, f2, f3 = coefficients
    assert f0 * h**4 + f1 * Nx * h**2 + f2 * Ny * h + f3 * Nx**2 == 0

    Px = PolynomialRing(K, "x")
    x = Px.gen()
    q = K(f0) + K(f1) * x + K(f3) * x**2
    total = q**2 - K(f2) ** 2 * (x**3 + K(A) * x + K(B))
    quotient, remainder = total.quo_rem(x - K(X))
    if remainder:
        raise ArithmeticError("the known trace section did not divide the quartic")
    if quotient.degree() != 3:
        raise ArithmeticError("the residual cover is not cubic")
    factor_degrees = sorted(int(factor.degree()) for factor, _ in quotient.factor())
    if factor_degrees != [3]:
        raise ArithmeticError(f"the deep residual cubic is reducible: {factor_degrees}")
    cover_denominator = R.one()
    for coefficient in quotient.list():
        cover_denominator = lcm(cover_denominator, R(coefficient.denominator()))
    cover = Px([K(cover_denominator * coefficient) for coefficient in quotient.list()])
    cover = Px([R(coefficient) for coefficient in cover.list()])
    # Scalar normalization is cosmetic; use a deterministic leading sign.
    if R(cover[cover.degree()]).leading_coefficient() < 0:
        cover = -cover
    discriminant = R(cover.discriminant())
    squarefree_discriminant = discriminant.squarefree_part()
    return {
        "h": h,
        "Nx": Nx,
        "Ny": Ny,
        "coefficients": coefficients,
        "cover": cover,
        "cover_denominator": cover_denominator,
        "discriminant": discriminant,
        "squarefree_discriminant": squarefree_discriminant,
        "interpolation_rank": interpolation.rank(),
    }


def specialize_record(data, controls, chart, R):
    f0, f1, f2, f3 = data["coefficients"]
    cover = data["cover"]
    output = []
    for label, numerator, denominator, known_rank in controls:
        parameter = QQ(numerator) / QQ(denominator)
        local_parameter = parameter if chart == "finite" else 1 / parameter
        specialized = PolynomialRing(QQ, "x")(
            [R(coefficient)(local_parameter) for coefficient in cover.list()]
        )
        factorization = specialized.factor()
        factor_degrees = sorted(
            [int(factor.degree()) for factor, exponent in factorization for _ in range(int(exponent))]
        )
        points = []
        for factor, exponent in factorization:
            if factor.degree() != 1:
                continue
            if exponent != 1:
                raise ArithmeticError("a control has a repeated rational cover root")
            x_value = -factor[0] / factor[1]
            denominator_y = f2(local_parameter)
            if denominator_y == 0:
                raise ArithmeticError("a split control lies in the f2=0 chart")
            y_value = -(
                f0(local_parameter)
                + f1(local_parameter) * x_value
                + f3(local_parameter) * x_value**2
            ) / denominator_y
            if chart == "finite":
                source_x = QQ(denominator**4) * x_value
                source_y = QQ(denominator**6) * y_value
            else:
                source_x = QQ(numerator**4) * x_value
                source_y = QQ(numerator**6) * y_value
            points.append([qtext(source_x), qtext(source_y)])
        output.append(
            {
                "label": label,
                "parameter": qtext(parameter),
                "known_rank_lower_bound": known_rank,
                "factor_degrees": factor_degrees,
                "rational_point_count": len(points),
                "projective_source_points": points,
            }
        )
    return output


def construct_one(row, context):
    foundry_vector = tuple(int(value) for value in row["minimum_representative"])
    key, published, score, sign = canonical_orientation(
        foundry_vector,
        context["frame_to_published"],
        context["direct_costs"],
        context["closures"],
    )
    published_vector = vector(ZZ, published)
    assert published_vector * context["published_gram"] * published_vector == 26
    tau = sum(
        (coefficient * point for coefficient, point in zip(published_vector, context["basis"])),
        context["curve"](0),
    )
    X, Y = context["K"](tau[0]), context["K"](tau[1])
    chart = "finite"
    x_denominator = context["R"](X.denominator())
    if not x_denominator.is_square():
        raise ArithmeticError("trace denominator is not a square")
    h = context["R"](x_denominator.sqrt())
    h /= h.leading_coefficient()
    if h.degree() != 11:
        chart = "inverted_at_infinity"
        X = invert_rational(X, 4, context["R"], context["K"])
        Y = invert_rational(Y, 6, context["R"], context["K"])
        A = reciprocal_with_bound(context["A"], 8, context["R"])
        B = reciprocal_with_bound(context["B"], 12, context["R"])
    else:
        A, B = context["A"], context["B"]
    data = interpolation_relation(X, Y, A, B, context["R"], context["K"])
    f0, f1, f2, f3 = data["coefficients"]
    cover = data["cover"]
    maximum_coefficient_bits = max(
        abs(QQ(value).numerator()).bit_length()
        + abs(QQ(value).denominator()).bit_length()
        for polynomial in (f0, f1, f2, f3)
        for value in polynomial.list()
    )
    cover_maximum_bits = max(
        abs(QQ(value).numerator()).bit_length()
        + abs(QQ(value).denominator()).bit_length()
        for polynomial in cover.list()
        for value in context["R"](polynomial).list()
    )
    return {
        "label": f"deep3-{int(row['residue_id']):08x}",
        "residue_id": int(row["residue_id"]),
        "residue_mod_3": list(map(int, row["residue"])),
        "foundry_frame_minimum_representative": [sign * value for value in foundry_vector],
        "published_basis_w": list(published),
        "minimum_norm": 26,
        "equation_complexity": {
            **score,
            "riemann_roch_coefficient_bounds": list(BOUNDS),
            "riemann_roch_unknown_count": 66,
            "interpolation_equation_count": 65,
            "interpolation_rank": data["interpolation_rank"],
            "maximum_interpolation_coefficient_bits": maximum_coefficient_bits,
            "maximum_cover_coefficient_bits": cover_maximum_bits,
        },
        "construction_chart": chart,
        "trace_section": {
            "h_coefficients_low_to_high": poly_record(data["h"]),
            "Nx_coefficients_low_to_high": poly_record(data["Nx"]),
            "Ny_coefficients_low_to_high": poly_record(data["Ny"]),
        },
        "riemann_roch_relation": {
            "identity": "f0+f1*x+f2*y+f3*x^2=0",
            "f_coefficients_low_to_high": [poly_record(value) for value in data["coefficients"]],
        },
        "cubic_cover": {
            "coefficients_in_x_low_to_high": [
                rational_function_record(context["K"](value)) for value in cover.list()
            ],
            "generic_factor_degrees_over_Qt": [3],
            "discriminant_degree": int(data["discriminant"].degree()),
            "squarefree_discriminant_degree": int(data["squarefree_discriminant"].degree()),
        },
        "specializations": specialize_record(data, CONTROLS, chart, context["R"]),
    }


def build_payload(start: int, limit: int | None):
    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    target = json.loads(TARGET.read_text())
    foundry = json.loads(FOUNDRY.read_text())
    deep = json.loads(DEEP.read_text())
    pinned = load_matrix(PINNED)
    published_gram, frame_gram, frame_to_published = published_gram_and_foundry_map(
        foundry, target, pinned
    )
    spectrum = next(item for item in deep["spectra"] if item["frame_id"] == "NS0001-F001")
    rows = spectrum["retained_inversion_representatives"]
    if len(rows) != 160 or spectrum["retained_full_coset_count_after_inversion"] != 320:
        raise ArithmeticError("the complete norm-26 shell is not pinned")
    selected = rows[start:] if limit is None else rows[start : start + limit]
    if not selected:
        raise ValueError("the selected deep-shell interval is empty")

    R = PolynomialRing(QQ, "t")
    K = R.fraction_field()
    A = R([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = R([QQ(value) for value in model["B_coefficients_low_to_high"]])
    basis_coordinates = reconstruct_basis(R, A, B, section_data)
    curve = EllipticCurve(K, [A, B])
    basis = [curve(K(x), K(y)) for x, y in basis_coordinates]
    direct_costs, closures = section_input_costs(section_data)
    context = {
        "R": R,
        "K": K,
        "A": A,
        "B": B,
        "curve": curve,
        "basis": basis,
        "published_gram": published_gram,
        "frame_gram": frame_gram,
        "frame_to_published": frame_to_published,
        "direct_costs": direct_costs,
        "closures": closures,
    }
    records = []
    for index, row in enumerate(selected, start=start):
        print(f"ELKIES2026DEEP3|stage=construct|index={index + 1}/160", flush=True)
        records.append(construct_one(row, context))
    return {
        "schema": "elkies-k3.elkies-2026-r17-deep-trisections.v1",
        "status": "PASS_EXACT_DEEP_TRISECTION_CONSTRUCTION_AND_SPECIALIZATION",
        "claim": (
            "Exact Riemann--Roch construction and rank-25--28 specialization of "
            "the selected inversion representatives in the complete norm-26 shell."
        ),
        "claim_boundary": (
            "Generic irreducibility and every displayed rational specialization are exact. "
            "The quotient directions and generated subgroup ranks are certified separately."
        ),
        "complete_shell": {
            "translation_cosets": 320,
            "inversion_representatives": 160,
            "minimum_norm": 26,
        },
        "selected_interval": {"start_zero_based": start, "record_count": len(records)},
        "records": records,
        "summary": {
            "split_cover_counts_rank28_to_rank25": [
                sum(record["specializations"][index]["rational_point_count"] > 0 for record in records)
                for index in range(4)
            ],
            "rational_point_counts_rank28_to_rank25": [
                sum(record["specializations"][index]["rational_point_count"] for record in records)
                for index in range(4)
            ],
        },
        "generation": {
            "command": (
                ".venv/bin/sage-python elkies-k3/scripts/"
                "construct_elkies_2026_deep_trisections.sage"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": {
                relative(path): digest(path)
                for path in (MODEL, SECTIONS, PINNED, TARGET, FOUNDRY, DEEP)
            },
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.start < 0 or (arguments.limit is not None and arguments.limit <= 0):
        parser.error("invalid selected interval")
    payload = build_payload(arguments.start, arguments.limit)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale deep-trisection certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "ELKIES2026DEEP3|records={}|splits={}|status={}".format(
            len(payload["records"]),
            ",".join(map(str, payload["summary"]["split_cover_counts_rank28_to_rank25"])),
            terminal,
        )
    )


if __name__ == "__main__":
    main()
