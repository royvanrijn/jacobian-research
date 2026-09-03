#!/usr/bin/env python3
"""Measure exact coefficient heights under two third-q12 field bases.

status: ACTIVE_COMPILER
claim: bounded 63-term comparison of theta^2, omega, and delta coordinates
inputs: exact connected pencil and certified biquadratic closure operands
outputs: elkies-k3-q80-third-q12-pencil-basis-heights-v1.json

Run this with the pinned Sage Python, which supplies gmpy2.  This is a
go/no-go compiler metric, not an equation reconstruction or a proof that one
basis is optimal under arbitrary field/model transformations.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from gmpy2 import gcd, isqrt, lcm, mpq, mpz


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OPERANDS = (
    ROOT
    / "artifacts/generated-results/"
    / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_PENCIL = (
    ROOT
    / "artifacts/generated-results/"
    / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    / "elkies-k3-q80-third-q12-pencil-basis-heights-v1.json"
)
COEFFICIENT = re.compile(
    r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$"
)


def rational(entry):
    return mpq(mpz(entry["numerator"]), mpz(entry["denominator"]))


def height_bits(value):
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def coordinate_height(values):
    return max(height_bits(value) for value in values)


def summarize(values):
    ordered = sorted(values)
    return {
        "minimum_bits": min(values),
        "median_bits": ordered[len(ordered) // 2],
        "maximum_bits": max(values),
        "sum_bits": sum(values),
    }


def rational_projective_height(coordinates):
    common_denominator = mpz(1)
    for value in coordinates:
        common_denominator = lcm(common_denominator, value.denominator)
    integer_coordinates = [
        value.numerator * (common_denominator // value.denominator)
        for value in coordinates
    ]
    content = mpz(0)
    for value in integer_coordinates:
        content = gcd(content, abs(value))
    primitive = [value // content for value in integer_coordinates]
    return {
        "coordinate_count": len(coordinates),
        "nonzero_coordinate_count": sum(value != 0 for value in coordinates),
        "common_denominator_bits": common_denominator.bit_length(),
        "integer_content_bits": content.bit_length(),
        "primitive_maximum_bits": max(abs(value).bit_length() for value in primitive),
    }


def source_record(path, encoded, data):
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "schema": data.get("schema"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
    parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write-artifact", action="store_true")
    args = parser.parse_args()

    operand_bytes = args.operands.read_bytes()
    pencil_bytes = args.pencil.read_bytes()
    operands = json.loads(operand_bytes)
    pencil = json.loads(pencil_bytes)
    field = operands["biquadratic_field"]
    q1 = rational(field["q1"])
    q2 = rational(field["q2"])
    product = q1 * q2
    numerator_root = isqrt(product.numerator)
    assert numerator_root * numerator_root == product.numerator
    product_denominator = product.denominator

    terms = pencil["moving_equation"]["terms_T_W_x_coefficient_1_r"]
    heights = {"theta_squared": [], "omega": [], "delta": []}
    coordinates = {"theta_squared": [], "omega": [], "delta": []}
    per_term = []
    for t_degree, w_degree, x_degree, encoded_coefficient in terms:
        assert len(encoded_coefficient) == 1
        match = COEFFICIENT.fullmatch(encoded_coefficient[0])
        assert match is not None
        theta2 = mpq(mpz(match[1]), mpz(match[2]))
        constant_sign = 1 if match[3] == "+" else -1
        constant = mpq(constant_sign * mpz(match[4]), mpz(match[5]))

        # theta^2 = (q1+q2) + omega/2 and
        # omega = 4*sqrt(N)*delta/D for delta^2=D.
        invariant = constant + theta2 * (q1 + q2)
        omega = theta2 / 2
        delta = 2 * theta2 * numerator_root / product_denominator
        current_coordinates = {
            "theta_squared": (constant, theta2),
            "omega": (invariant, omega),
            "delta": (invariant, delta),
        }
        current = {}
        for basis, values in current_coordinates.items():
            current[basis] = coordinate_height(values)
            heights[basis].append(current[basis])
            coordinates[basis].extend(values)
        per_term.append(
            {
                "degrees_T_W_x": [t_degree, w_degree, x_degree],
                "height_bits": current,
            }
        )

    summaries = {basis: summarize(values) for basis, values in heights.items()}
    rational_projective = {
        basis: rational_projective_height(values)
        for basis, values in coordinates.items()
    }
    omega_values = heights["omega"]
    delta_values = heights["delta"]
    delta_better = sum(d < o for d, o in zip(delta_values, omega_values))
    delta_equal = sum(d == o for d, o in zip(delta_values, omega_values))
    delta_worse = len(terms) - delta_better - delta_equal

    payload = {
        "schema": "elkies-k3-q80-third-q12-pencil-basis-heights-v1",
        "status": "PASS_EXACT_63_TERM_BOUNDED_METRIC",
        "inputs": {
            "operands": source_record(args.operands, operand_bytes, operands),
            "pencil": source_record(args.pencil, pencil_bytes, pencil),
            "checker": {
                "path": str(Path(__file__).resolve().relative_to(ROOT)),
                "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            },
        },
        "term_count": len(terms),
        "height_definition": (
            "maximum numerator/denominator bit length among the two rational "
            "coordinates of each coefficient, before any global model rescaling"
        ),
        "basis_relations": {
            "theta_squared": "theta^2",
            "omega": "theta^2=(q1+q2)+omega/2",
            "delta": (
                "delta^2=D, omega=4*sqrt(N)*delta/D for reduced q1*q2=N/D"
            ),
        },
        "summaries": summaries,
        "rational_projective_normalization": rational_projective,
        "delta_vs_omega_term_counts": {
            "lower": delta_better,
            "equal": delta_equal,
            "higher": delta_worse,
        },
        "maximum_delta_minus_omega_bits": (
            summaries["delta"]["maximum_bits"]
            - summaries["omega"]["maximum_bits"]
        ),
        "raw_conclusion_gate": (
            "DELTA_LOWERS_RAW_COORDINATE_HEIGHT"
            if summaries["delta"]["maximum_bits"]
            < summaries["omega"]["maximum_bits"]
            else "DELTA_DOES_NOT_LOWER_RAW_MAXIMUM_COORDINATE_HEIGHT"
        ),
        "rational_projective_delta_minus_omega_bits": (
            rational_projective["delta"]["primitive_maximum_bits"]
            - rational_projective["omega"]["primitive_maximum_bits"]
        ),
        "projective_conclusion_gate": (
            "DELTA_LOWERS_RATIONAL_PROJECTIVE_HEIGHT_ONLY_MARGINALLY"
            if rational_projective["delta"]["primitive_maximum_bits"]
            < rational_projective["omega"]["primitive_maximum_bits"]
            else "DELTA_DOES_NOT_LOWER_RATIONAL_PROJECTIVE_HEIGHT"
        ),
        "claim_boundary": (
            "The comparison covers exactly the 63 coefficients of the pinned "
            "moving equation. It does not optimize global rescaling, base PGL2, "
            "Weierstrass transformations, or integral ideals."
        ),
        "per_term": per_term,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/audit_q80_third_q12_pencil_basis_heights.py "
            "--write-artifact"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    if args.write_artifact:
        args.output.write_text(encoded)
        print(f"Q80Q12BASIS|artifact={args.output}|sha256={digest}|status=PASS_WRITE")
    elif args.check:
        assert args.output.read_text() == encoded
        print(f"Q80Q12BASIS|artifact={args.output}|sha256={digest}|status=PASS_CHECK")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
