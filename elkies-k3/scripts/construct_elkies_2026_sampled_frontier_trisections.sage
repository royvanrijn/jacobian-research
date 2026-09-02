#!/usr/bin/env sage-python
"""Construct the deterministic norm-20 frontier R17 trisection sample.

This complements the complete 320-coset norm-26 deep shell with the pinned
1,025-coset inversion-closed degree-three sample.  The sample contains 138
rational frontier vertices of norm 20, or 69 inversion pairs.  Each trace is
removed from the unique section of O_X(4O+16F) through it, leaving a residual
cubic which is factored generically and at the rank-25--28 controls.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, lcm, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEEP_SCRIPT = ROOT / "elkies-k3/scripts/construct_elkies_2026_deep_trisections.sage"
DIVERSITY_SCRIPT = ROOT / "elkies-k3/scripts/analyze_r17_multisection_diversity.py"
DEEP = SourceFileLoader("r17_deep3_frontier_helpers", str(DEEP_SCRIPT)).load_module()
DIVERSITY = SourceFileLoader("r17_diversity_frontier_helpers", str(DIVERSITY_SCRIPT)).load_module()

MODEL = DEEP.MODEL
SECTIONS = DEEP.SECTIONS
PINNED = DEEP.PINNED
TARGET = DEEP.TARGET
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
DIVERSITY_CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-sampled-frontier-trisections-v1.json"

SAMPLE_SIZE = 1024
SAMPLE_SEED = 20260902
BOUNDS = (16, 12, 10, 8)


def interpolation_relation(X, Y, A, B, R, K):
    denominator = R(X.denominator())
    if not denominator.is_square():
        raise ArithmeticError("trace x denominator is not a square")
    h = R(denominator.sqrt())
    h /= h.leading_coefficient()
    assert X.denominator() == h**2 and Y.denominator() == h**3
    if h.degree() != 8:
        raise ArithmeticError(f"trace denominator degree {h.degree()} is not 8")
    Nx = R(X * h**2)
    Ny = R(Y * h**3)
    terms = (h**4, Nx * h**2, Ny * h, Nx**2)
    assert all(term.degree() <= 48 for term in terms)
    columns = []
    for term, bound in zip(terms, BOUNDS):
        for power in range(bound + 1):
            polynomial = term * R.gen()**power
            columns.append([polynomial[index] for index in range(49)])
    interpolation = matrix(
        QQ, 49, len(columns), lambda row, column: columns[column][row]
    )
    if interpolation.ncols() != 50 or interpolation.rank() != 49:
        raise ArithmeticError(
            f"frontier interpolation has rank {interpolation.rank()}, expected 49"
        )
    raw = interpolation.right_kernel().basis()[0]
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
    if remainder or quotient.degree() != 3:
        raise ArithmeticError("the trace did not leave a residual cubic")
    factor_degrees = sorted(
        int(factor.degree())
        for factor, exponent in quotient.factor()
        for _ in range(int(exponent))
    )
    cover_denominator = R.one()
    for coefficient in quotient.list():
        cover_denominator = lcm(cover_denominator, R(coefficient.denominator()))
    cover = Px([R(cover_denominator * coefficient) for coefficient in quotient.list()])
    if R(cover[cover.degree()]).leading_coefficient() < 0:
        cover = -cover
    discriminant = R(cover.discriminant())
    return {
        "h": h,
        "Nx": Nx,
        "Ny": Ny,
        "coefficients": coefficients,
        "cover": cover,
        "generic_factor_degrees": factor_degrees,
        "discriminant": discriminant,
        "squarefree_discriminant": discriminant.squarefree_part(),
        "interpolation_rank": interpolation.rank(),
    }


def specialize(data, chart, R):
    f0, f1, f2, f3 = data["coefficients"]
    output = []
    for label, numerator, denominator, known_rank in DEEP.CONTROLS:
        parameter = QQ(numerator) / QQ(denominator)
        local_parameter = parameter if chart == "finite" else 1 / parameter
        specialized = PolynomialRing(QQ, "x")(
            [R(coefficient)(local_parameter) for coefficient in data["cover"].list()]
        )
        factorization = specialized.factor()
        factor_degrees = sorted(
            int(factor.degree())
            for factor, exponent in factorization
            for _ in range(int(exponent))
        )
        points = []
        for factor, exponent in factorization:
            if factor.degree() != 1:
                continue
            if exponent != 1:
                raise ArithmeticError("a control has a repeated rational root")
            x_value = -factor[0] / factor[1]
            if f2(local_parameter) == 0:
                raise ArithmeticError("a split control lies outside the y chart")
            y_value = -(
                f0(local_parameter)
                + f1(local_parameter) * x_value
                + f3(local_parameter) * x_value**2
            ) / f2(local_parameter)
            if chart == "finite":
                source_x = QQ(denominator**4) * x_value
                source_y = QQ(denominator**6) * y_value
            else:
                source_x = QQ(numerator**4) * x_value
                source_y = QQ(numerator**6) * y_value
            points.append([DEEP.qtext(source_x), DEEP.qtext(source_y)])
        output.append(
            {
                "label": label,
                "parameter": DEEP.qtext(parameter),
                "known_rank_lower_bound": known_rank,
                "factor_degrees": factor_degrees,
                "rational_point_count": len(points),
                "projective_source_points": points,
            }
        )
    return output


def sampled_rows(short_gram):
    residues = DIVERSITY.aut_closed_sample(3, SAMPLE_SIZE, SAMPLE_SEED)
    if len(residues) != 1025:
        raise ArithmeticError("the pinned degree-three sample size changed")
    oracle = DIVERSITY.CosetOracle(
        [[int(value) for value in row] for row in short_gram.rows()], 3
    )
    rows = []
    histogram = Counter()
    for residue in residues:
        norm, representative, error = oracle.solve(residue)
        histogram[norm] += 1
        if norm != 20:
            continue
        negative = DIVERSITY.negate_residue(residue, 3)
        identifier = DIVERSITY.residue_id(residue, 3)
        negative_identifier = DIVERSITY.residue_id(negative, 3)
        if identifier > negative_identifier:
            continue
        rows.append(
            {
                "residue_id": identifier,
                "residue_mod_3": list(map(int, residue)),
                "minimum_representative": list(map(int, representative)),
                "minimum_norm": norm,
                "cvp_rounding_error": error,
            }
        )
    if histogram[20] != 138 or len(rows) != 69:
        raise ArithmeticError("the pinned norm-20 sample changed")
    return rows, histogram


def construct_one(row, context):
    short_vector = vector(ZZ, row["minimum_representative"])
    published_raw = short_vector * context["short_to_published"]
    choices = []
    for sign in (1, -1):
        published = tuple(sign * int(value) for value in published_raw)
        key, score = DEEP.equation_score(
            published, context["direct_costs"], context["closures"]
        )
        choices.append((key, published, score, sign))
    unused_key, published, score, sign = min(choices, key=lambda item: item[0])
    published_vector = vector(ZZ, published)
    assert published_vector * context["published_gram"] * published_vector == 20
    tau = sum(
        (coefficient * point for coefficient, point in zip(published_vector, context["basis"])),
        context["curve"](0),
    )
    X, Y = context["K"](tau[0]), context["K"](tau[1])
    denominator = context["R"](X.denominator())
    h = context["R"](denominator.sqrt())
    h /= h.leading_coefficient()
    chart = "finite"
    if h.degree() != 8:
        chart = "inverted_at_infinity"
        X = DEEP.invert_rational(X, 4, context["R"], context["K"])
        Y = DEEP.invert_rational(Y, 6, context["R"], context["K"])
        A = DEEP.reciprocal_with_bound(context["A"], 8, context["R"])
        B = DEEP.reciprocal_with_bound(context["B"], 12, context["R"])
    else:
        A, B = context["A"], context["B"]
    data = interpolation_relation(X, Y, A, B, context["R"], context["K"])
    maximum_relation_bits = max(
        abs(QQ(value).numerator()).bit_length()
        + abs(QQ(value).denominator()).bit_length()
        for polynomial in data["coefficients"]
        for value in polynomial.list()
    )
    maximum_cover_bits = max(
        abs(QQ(value).numerator()).bit_length()
        + abs(QQ(value).denominator()).bit_length()
        for coefficient in data["cover"].list()
        for value in context["R"](coefficient).list()
    )
    return {
        "label": f"sampled3f-{int(row['residue_id']):08x}",
        "residue_id": int(row["residue_id"]),
        "residue_mod_3": row["residue_mod_3"],
        "short_basis_minimum_representative": [
            sign * int(value) for value in short_vector
        ],
        "published_basis_w": list(published),
        "minimum_norm": 20,
        "equation_complexity": {
            **score,
            "riemann_roch_coefficient_bounds": list(BOUNDS),
            "riemann_roch_unknown_count": 50,
            "interpolation_equation_count": 49,
            "interpolation_rank": data["interpolation_rank"],
            "maximum_interpolation_coefficient_bits": maximum_relation_bits,
            "maximum_cover_coefficient_bits": maximum_cover_bits,
        },
        "construction_chart": chart,
        "trace_section": {
            "h_coefficients_low_to_high": DEEP.poly_record(data["h"]),
            "Nx_coefficients_low_to_high": DEEP.poly_record(data["Nx"]),
            "Ny_coefficients_low_to_high": DEEP.poly_record(data["Ny"]),
        },
        "riemann_roch_relation": {
            "identity": "f0+f1*x+f2*y+f3*x^2=0",
            "f_coefficients_low_to_high": [
                DEEP.poly_record(value) for value in data["coefficients"]
            ],
        },
        "cubic_cover": {
            "coefficients_in_x_low_to_high": [
                DEEP.rational_function_record(context["K"](value))
                for value in data["cover"].list()
            ],
            "generic_factor_degrees_over_Qt": data["generic_factor_degrees"],
            "discriminant_degree": int(data["discriminant"].degree()),
            "squarefree_discriminant_degree": int(
                data["squarefree_discriminant"].degree()
            ),
        },
        "specializations": specialize(data, chart, context["R"]),
    }


def build_payload(start, limit):
    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    target = json.loads(TARGET.read_text())
    pinned = DEEP.load_matrix(PINNED)
    short_coordinates = DEEP.load_matrix(SHORT_COORDS)
    short_gram = DEEP.load_matrix(SHORT_GRAM)
    assert short_gram == short_coordinates * pinned * short_coordinates.transpose()
    basis_change = matrix(ZZ, target["pinned_identification"]["basis_change_matrix"])
    published_gram = matrix(
        ZZ, basis_change.transpose().inverse() * pinned * basis_change.inverse()
    )
    short_to_published = short_coordinates * basis_change.transpose()
    assert short_to_published * published_gram * short_to_published.transpose() == short_gram
    rows, histogram = sampled_rows(short_gram)
    selected = rows[start:] if limit is None else rows[start : start + limit]
    if not selected:
        raise ValueError("the selected frontier interval is empty")
    R = PolynomialRing(QQ, "t")
    K = R.fraction_field()
    A = R([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = R([QQ(value) for value in model["B_coefficients_low_to_high"]])
    basis_coordinates = DEEP.reconstruct_basis(R, A, B, section_data)
    curve = EllipticCurve(K, [A, B])
    basis = [curve(K(x), K(y)) for x, y in basis_coordinates]
    direct_costs, closures = DEEP.section_input_costs(section_data)
    context = {
        "R": R,
        "K": K,
        "A": A,
        "B": B,
        "curve": curve,
        "basis": basis,
        "published_gram": published_gram,
        "short_to_published": short_to_published,
        "direct_costs": direct_costs,
        "closures": closures,
    }
    records = []
    for index, row in enumerate(selected, start=start):
        print(f"ELKIES2026SAMPLED3F|stage=construct|index={index + 1}/69", flush=True)
        records.append(construct_one(row, context))
    return {
        "schema": "elkies-k3.elkies-2026-r17-sampled-frontier-trisections.v1",
        "status": "PASS_EXACT_SAMPLED_FRONTIER_TRISECTION_CONSTRUCTION_AND_SPECIALIZATION",
        "claim": (
            "Exact construction and rank-25--28 specialization of one representative "
            "from every norm-20 inversion pair in the pinned degree-three sample."
        ),
        "claim_boundary": (
            "The norm-20 frontier has 18,024,296 translation cosets in the complete "
            "lattice census. The 1,025-coset equation universe here is deterministic "
            "and inversion-closed but is not a complete equation atlas."
        ),
        "tested_sample": {
            "degree": 3,
            "seed": SAMPLE_SEED,
            "requested_sample_size": SAMPLE_SIZE,
            "actual_inversion_closed_cosets": sum(histogram.values()),
            "minimum_norm_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "norm_20_rational_vertices": 138,
            "norm_20_inversion_representatives": 69,
        },
        "selected_interval": {"start_zero_based": start, "record_count": len(records)},
        "records": records,
        "summary": {
            "generic_factor_patterns": dict(
                sorted(
                    Counter(
                        "+".join(map(str, record["cubic_cover"]["generic_factor_degrees_over_Qt"]))
                        for record in records
                    ).items()
                )
            ),
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
                "construct_elkies_2026_sampled_frontier_trisections.sage"
            ),
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "inputs": {
                DEEP.relative(path): sha256(path.read_bytes()).hexdigest()
                for path in (
                    MODEL, SECTIONS, PINNED, TARGET, SHORT_COORDS, SHORT_GRAM,
                    DIVERSITY_CERTIFICATE,
                )
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
    payload = build_payload(arguments.start, arguments.limit)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale sampled-frontier-trisection certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "ELKIES2026SAMPLED3F|records={}|splits={}|patterns={}|status={}".format(
            len(payload["records"]),
            ",".join(map(str, payload["summary"]["split_cover_counts_rank28_to_rank25"])),
            payload["summary"]["generic_factor_patterns"],
            terminal,
        )
    )


if __name__ == "__main__":
    main()
