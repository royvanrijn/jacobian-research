#!/usr/bin/env sage-python
"""Compile alternate-Q80 norm-ten bisections from the complete priority table.

Each height-ten trace on the canonical alternate-Q80 24I1 model is converted
to its unique regular residual chord.  Proposition F1 then gives an exact
quadratic cover and a lifted section.  The script supports deterministic
intervals so the 39,147-record batch can be replayed in bounded chunks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, prime_range, vector,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
PRIORITY = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.tsv"
PRIORITY_CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.json"
ORBIT_TABLE = ROOT / "artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.tsv"
HISTORICAL_FRAME = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-v1.json"
HIDDEN_DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json"
HIDDEN_PRIORITY = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisection-priority-v1.tsv"
HIDDEN_PRIORITY_CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisection-priority-v1.json"
HIDDEN_DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisections-v1.json"
CONTENT_TRIAL_PRIMES = tuple(prime_range(2, 1001))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def polynomial_coefficients(polynomial) -> list[str]:
    if not polynomial:
        return ["0"]
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def integer_square_decomposition(integer):
    """Return exact ``square, residual`` with integer = square^2 * residual.

    This deliberately uses only bounded trial division followed by a perfect-
    square test.  It never starts a general factorization of a large residual.
    Leaving an opaque composite in ``residual`` is harmless: the downstream
    factorless collision partition compares rational constants by an exact
    square-ratio test.
    """

    residual = ZZ(integer)
    if residual <= 0:
        raise ArithmeticError("integer square decomposition requires a positive input")
    square = ZZ.one()
    retained = ZZ.one()
    for prime in CONTENT_TRIAL_PRIMES:
        exponent = 0
        while residual % prime == 0:
            residual //= prime
            exponent += 1
        square *= prime ** (exponent // 2)
        if exponent % 2:
            retained *= prime
    if residual.is_square():
        square *= residual.sqrt()
    else:
        retained *= residual
    if integer != square**2 * retained:
        raise ArithmeticError("bounded integer square decomposition identity failed")
    return square, retained


def primitive_square_normalization(polynomial, ring):
    """Remove detectable square content while retaining its exact squareclass.

    The trace-chord formulas produce ``q`` only up to a rational square.  In
    numerator and denominator can have roughly a thousand bits and usually
    contain enormous powers.  Bounded trial division plus perfect-square tests
    strip those powers without a potentially unbounded integer factorization.
    The retained integral content need not itself be squarefree, but it is an
    exact representative of the rational constant squareclass.
    """

    denominators = [coefficient.denominator() for coefficient in polynomial]
    common_denominator = lcm(denominators)
    integral_coefficients = [
        ZZ(coefficient * common_denominator) for coefficient in polynomial
    ]
    common_numerator = gcd(integral_coefficients)
    if common_numerator < 0:
        common_numerator = -common_numerator
    content = QQ(common_numerator) / common_denominator
    if content <= 0:
        raise ArithmeticError("trace-chord q has nonpositive rational content")
    primitive = ring(polynomial / content)
    if any(coefficient.denominator() != 1 for coefficient in primitive):
        raise ArithmeticError("normalized trace-chord q is not integral")
    if gcd([ZZ(coefficient) for coefficient in primitive]) not in (1, -1):
        raise ArithmeticError("normalized trace-chord q is not primitive")
    numerator_square, numerator_residual = integer_square_decomposition(common_numerator)
    denominator_square, denominator_residual = integer_square_decomposition(common_denominator)
    squareclass_content = numerator_residual * denominator_residual
    content_square_root = QQ(numerator_square) / (denominator_square * denominator_residual)
    normalized = ring(squareclass_content * primitive)
    if polynomial != content_square_root**2 * normalized:
        raise ArithmeticError("rational square normalization identity failed")
    return normalized, content_square_root, squareclass_content


def parse_vector(text: str):
    result = vector(ZZ, [ZZ(value) for value in text.split()])
    if len(result) != 17:
        raise ValueError("expected 17 integral coordinates")
    return result


def parse_rational_function(record, ring, field):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return field(numerator / denominator)


def compile_trace(X, Y, A, B, Delta, ring, field, helper):
    """Apply the exact finite chart, with the certified reciprocal fallback."""

    trace_chord_frame = helper["trace_chord_frame"]
    local_chord_data = helper["local_chord_data"]
    reciprocal_with_bound = helper["reciprocal_with_bound"]
    invert_rational = helper["invert_rational"]
    frame = trace_chord_frame(X, Y, ring)
    if frame["h"].degree() == 3:
        data = local_chord_data(X, Y, A, B, Delta, ring, field)
        data["construction_chart"] = "finite"
        return data
    A_inverse = reciprocal_with_bound(A, 8, ring)
    B_inverse = reciprocal_with_bound(B, 12, ring)
    Delta_inverse = reciprocal_with_bound(Delta, 24, ring)
    X_inverse = invert_rational(X, 4, ring, field)
    Y_inverse = invert_rational(Y, 6, ring, field)
    inverse = local_chord_data(
        X_inverse, Y_inverse, A_inverse, B_inverse, Delta_inverse, ring, field
    )
    data = {
        "h": reciprocal_with_bound(inverse["h"], 3, ring),
        "M": reciprocal_with_bound(inverse["M"], 5, ring),
        "q": reciprocal_with_bound(inverse["q"], 2, ring),
        "sum_x": reciprocal_with_bound(inverse["sum_x"], 4, ring),
        "product_x": reciprocal_with_bound(inverse["product_x"], 8, ring),
        "x0": reciprocal_with_bound(inverse["x0"], 4, ring),
        "x1": reciprocal_with_bound(inverse["x1"], 3, ring),
        "y0": reciprocal_with_bound(inverse["y0"], 6, ring),
        "y1": reciprocal_with_bound(inverse["y1"], 5, ring),
        "construction_chart": "inverted_at_infinity",
    }
    data["Nx"] = ring(X * data["h"] ** 2)
    data["Ny"] = ring(Y * data["h"] ** 3)
    data["branch_fibres_smooth"] = data["q"].gcd(Delta).degree() == 0
    if data["sum_x"] ** 2 - 4 * data["product_x"] != data["h"] ** 2 * data["q"]:
        raise ArithmeticError("reciprocal quadratic discriminant identity failed")
    if data["y0"] ** 2 + data["y1"] ** 2 * data["q"] != (
        data["x0"] ** 3
        + 3 * data["x0"] * data["x1"] ** 2 * data["q"]
        + A * data["x0"]
        + B
    ):
        raise ArithmeticError("reciprocal lifted-section constant identity failed")
    if 2 * data["y0"] * data["y1"] != (
        3 * data["x0"] ** 2 * data["x1"]
        + data["x1"] ** 3 * data["q"]
        + A * data["x1"]
    ):
        raise ArithmeticError("reciprocal lifted-section linear identity failed")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-label",
        choices=("norm12-orbit-11952", "norm12-orbit-103b2"),
        default="norm12-orbit-11952",
    )
    parser.add_argument("--priority-table", type=Path)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--verbose-records", action="store_true",
        help="retain the reconstructed trace, chord, and image quadratic in every record",
    )
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    is_alternate_target = arguments.source_label == "norm12-orbit-11952"
    direct_path = DIRECT if is_alternate_target else HIDDEN_DIRECT
    priority_certificate_path = (
        PRIORITY_CERTIFICATE if is_alternate_target else HIDDEN_PRIORITY_CERTIFICATE
    )
    default_priority = PRIORITY if is_alternate_target else HIDDEN_PRIORITY
    priority_table = arguments.priority_table or default_priority
    output = arguments.output or (
        DEFAULT_OUTPUT if is_alternate_target else HIDDEN_DEFAULT_OUTPUT
    )
    expected_count = 39147 if is_alternate_target else 39120
    if arguments.start < 0 or arguments.limit is not None and arguments.limit <= 0:
        parser.error("--start must be nonnegative and --limit must be positive")

    direct = json.loads(direct_path.read_text())
    priority_certificate = json.loads(priority_certificate_path.read_text())
    if direct["weierstrass_model"]["fibre_configuration"] != "24 I1":
        raise ArithmeticError("canonical alternate model is not 24I1")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("canonical alternate section basis is not saturated")
    expected_priority_status = (
        "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATION_PRIORITY"
        if is_alternate_target
        else "PASS_EXACT_COMPLETE_103B2_BISECTION_EQUATION_PRIORITY"
    )
    if priority_certificate["status"] != expected_priority_status:
        raise ArithmeticError("bisection priority certificate is not complete")
    with priority_table.open(newline="") as stream:
        all_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(all_rows) != expected_count:
        raise ArithmeticError(
            f"priority table does not contain all {expected_count} classes"
        )
    stop = len(all_rows) if arguments.limit is None else min(
        len(all_rows), arguments.start + arguments.limit
    )
    selected_rows = all_rows[arguments.start:stop]
    if not selected_rows:
        raise ArithmeticError("selected interval is empty")

    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    A = ring([QQ(value) for value in direct["weierstrass_model"]["A_coefficients_low_to_high"]])
    B = ring([QQ(value) for value in direct["weierstrass_model"]["B_coefficients_low_to_high"]])
    Delta = ring([QQ(value) for value in direct["weierstrass_model"]["discriminant_coefficients_low_to_high"]])
    if Delta != -16 * (4 * A**3 + 27 * B**2):
        raise ArithmeticError("stored alternate discriminant does not match A and B")
    E = EllipticCurve(field, [A, B])
    basis = []
    for expected_index, record in enumerate(direct["sections"]["records"]):
        if record["basis_index"] != expected_index or record["equation_verified"] is not True:
            raise ArithmeticError("alternate section record order or status changed")
        X = parse_rational_function(record["X"], ring, field)
        Y = parse_rational_function(record["Y"], ring, field)
        point = E(X, Y)
        basis.append(point)
    if len(basis) != 17:
        raise ArithmeticError("alternate equation basis does not have rank 17")
    height_gram = matrix(ZZ, direct["sections"]["height_gram"])

    maximum_coefficient = max(
        abs(int(value))
        for row in selected_rows
        for value in row["section_basis_w"].split()
    )
    multiples = []
    for point in basis:
        multiples.append(
            {coefficient: coefficient * point for coefficient in range(-maximum_coefficient, maximum_coefficient + 1)}
        )
    helper = runpy.run_path(str(HELPER))

    records = []
    chart_counts = {"finite": 0, "inverted_at_infinity": 0}
    for offset, row in enumerate(selected_rows):
        section_vector = parse_vector(row["section_basis_w"])
        direct_vector = parse_vector(
            row["direct_alternate_w"] if is_alternate_target else row["direct_hidden_w"]
        )
        source_vector = parse_vector(
            row["historical_alternate_w"] if is_alternate_target else row["short_basis_w"]
        )
        if section_vector * height_gram * section_vector != 10:
            raise ArithmeticError("priority trace word does not have height ten")
        trace = sum(
            (multiples[index][int(coefficient)] for index, coefficient in enumerate(section_vector)),
            E(0),
        )
        if trace.is_zero():
            raise ArithmeticError("height-ten trace unexpectedly vanished")
        X, Y = field(trace[0]), field(trace[1])
        data = compile_trace(X, Y, A, B, Delta, ring, field, helper)
        chart = data["construction_chart"]
        chart_counts[chart] += 1
        if not data["branch_fibres_smooth"]:
            raise ArithmeticError("alternate bisection branches over a singular fibre")
        h, Nx, Ny, M, raw_q = (data[key] for key in ("h", "Nx", "Ny", "M", "q"))
        q, q_content_square_root, q_squarefree_content = primitive_square_normalization(raw_q, ring)
        h = ring(h * q_content_square_root)
        Nx = ring(Nx * q_content_square_root**2)
        Ny = ring(Ny * q_content_square_root**3)
        M = ring(M * q_content_square_root)
        sum_x, product_x = data["sum_x"], data["product_x"]
        x0, x1, y0, y1 = (data[key] for key in ("x0", "x1", "y0", "y1"))
        x1 = ring(x1 * q_content_square_root)
        y1 = ring(y1 * q_content_square_root)
        if sum_x**2 - 4 * product_x != h**2 * q:
            raise ArithmeticError("normalized quadratic discriminant identity failed")
        if y0**2 + y1**2 * q != x0**3 + 3 * x0 * x1**2 * q + A * x0 + B:
            raise ArithmeticError("normalized lifted-section constant identity failed")
        if 2 * y0 * y1 != 3 * x0**2 * x1 + x1**3 * q + A * x1:
            raise ArithmeticError("normalized lifted-section linear identity failed")
        orbit = int(row["orbit_mask"], 0)
        label = (
            f"alternate-orbit-{orbit:05x}"
            if is_alternate_target else f"hidden-103b2-orbit-{orbit:05x}"
        )
        record = {
                "label": label,
                "lattice_orbit_mask": orbit,
                "section_basis_w": list(map(int, section_vector)),
                "priority_rank": int(row["priority_rank"]),
                "equation_complexity": {
                    key: int(row[key])
                    for key in (
                        "group_addition_upper_bound",
                        "support_count",
                        "maximum_absolute_coefficient",
                        "coefficient_l1",
                    )
                },
                "branch": {
                    "numerator_coefficients": polynomial_coefficients(q),
                    "denominator_coefficients": ["1"],
                },
                "residual_chord_certificate": {
                    "construction_chart": chart,
                    "raw_q_rational_square_multiplier": rational_text(q_content_square_root),
                    "normalized_q_integer_content_representative": rational_text(q_squarefree_content),
                    "squareclass_normalization_identity": "raw_q=(raw_q_rational_square_multiplier)^2*q",
                    "image_discriminant_identity_verified": True,
                    "branch_fibres_smooth": True,
                },
                "lifted_section": {
                    "field": "QQ(u,s), s^2=q(u)",
                    "x": "x0(u)+x1(u)*s",
                    "y": "y0(u)+y1(u)*s",
                    "x0_coefficients": polynomial_coefficients(x0),
                    "x1_coefficients": polynomial_coefficients(x1),
                    "y0_coefficients": polynomial_coefficients(y0),
                    "y1_coefficients": polynomial_coefficients(y1),
                    "constant_and_linear_identities_verified": True,
                    "two_branches_verified": True,
                    "anti_invariant_height": 12,
                },
            }
        if is_alternate_target:
            record["alternate_rank17_w"] = list(map(int, source_vector))
            record["direct_alternate_w"] = list(map(int, direct_vector))
        else:
            record["direct_hidden_w"] = list(map(int, direct_vector))
            record["short_basis_w"] = list(map(int, source_vector))
        if arguments.verbose_records:
            record["trace_section"] = {
                "h_coefficients": polynomial_coefficients(h),
                "Nx_coefficients": polynomial_coefficients(Nx),
                "Ny_coefficients": polynomial_coefficients(Ny),
            }
            record["residual_chord_reconstruction"] = {
                "slope": "M(u)/h(u)",
                "M_coefficients": polynomial_coefficients(M),
                "linear_congruence": "M*Nx+Ny == 0 mod h^2",
                "sum_x_coefficients": polynomial_coefficients(sum_x),
                "product_x_coefficients": polynomial_coefficients(product_x),
                "discriminant_identity": "sum_x^2-4*product_x=h^2*q",
            }
        records.append(record)
        completed = arguments.start + offset + 1
        if offset == 0 or (offset + 1) % 100 == 0 or completed == stop:
            print(
                f"ALTFULLBISECT|completed={completed}/{stop}|chunk={offset + 1}/{len(selected_rows)}|label={label}",
                flush=True,
            )

    input_paths = [HELPER, direct_path, priority_certificate_path, priority_table]
    if is_alternate_target:
        input_paths.extend([ORBIT_TABLE, HISTORICAL_FRAME])
    payload = {
        "schema": "elkies-k3.bisection-extension-input.v1",
        "artifact_schema": (
            "elkies-k3.r17-norm12-11952-alternate-bisections.v1"
            if is_alternate_target
            else "elkies-k3.r17-norm12-103b2-hidden-bisections.v1"
        ),
        "status": (
            "PASS_EXACT_ALTERNATE_BISECTION_EQUATION_CHUNK"
            if is_alternate_target
            else "PASS_EXACT_103B2_HIDDEN_BISECTION_EQUATION_CHUNK"
        ),
        "base_parameter": "u",
        "invariant_mw_rank": 17,
        "interval": {"start_zero_based": arguments.start, "stop_exclusive": stop},
        "bisections": records,
        "construction": {
            "method": "Proposition F1 exact regular residual chord",
            "record_count": len(records),
            "construction_chart_counts": chart_counts,
            "all_branch_fibres_smooth": True,
            "all_lifted_sections_verified": True,
            "verbose_records": bool(arguments.verbose_records),
        },
        "inputs": {
            relative(path): digest(path)
            for path in input_paths
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": ["Sage exact function fields", "exact elliptic group law"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/construct_r17_norm12_11952_alternate_bisections.sage "
            + ("" if is_alternate_target else "--source-label norm12-orbit-103b2 ")
            + f"--start {arguments.start} --limit {len(records)} --output {relative(output)}"
            + (" --verbose-records" if arguments.verbose_records else "")
        ),
        "proof_boundary": (
            "Every record in the declared interval has an exact quadratic relation, "
            "a squarefree degree-two branch polynomial coprime to the 24I1 discriminant, "
            f"and a lifted section verified coefficientwise. Complete {expected_count}-class coverage "
            "requires merging all disjoint intervals and running the squareclass checker."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored alternate bisection chunk differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "NORM12FULLBISECT|status={}|records={}|output={}".format(
            payload["status"], len(records), relative(output)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
