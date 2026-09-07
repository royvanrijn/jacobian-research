#!/usr/bin/env python3
"""Audit an exact integral normalization of the third-q12 descent field.

status: ACTIVE_COMPILER
claim: identify an exact denominator-integral generator for QQ(sqrt(q1*q2))
inputs: certified characteristic-zero biquadratic closure operands
outputs: elkies-k3-q80-third-q12-descent-field-normalization-v1.json

The connected pencil is defined over QQ(sqrt(q1*q2)).  This dependency-free
checker reads the certified closure operands and tests whether the reduced
numerator or denominator of q1*q2 is a square.  No factorization is attempted.
"""

import argparse
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
import sys


sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/"
    / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    / "elkies-k3-q80-third-q12-descent-field-normalization-v1.json"
)


def rational(entry):
    return Fraction(int(entry["numerator"]), int(entry["denominator"]))


def is_square(value):
    root = isqrt(value)
    return root * root == value, root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write-artifact", action="store_true")
    args = parser.parse_args()

    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    field = source["biquadratic_field"]
    q1 = rational(field["q1"])
    q2 = rational(field["q2"])
    product = q1 * q2
    numerator_square, numerator_root = is_square(product.numerator)
    denominator_square, denominator_root = is_square(product.denominator)

    # This identity is the useful exact normalization:
    # sqrt(N/D) = sqrt(N)*sqrt(D)/D when N is a square.
    assert numerator_square
    assert not denominator_square
    assert Fraction(numerator_root, product.denominator) ** 2 * product.denominator == product

    payload = {
        "schema": "elkies-k3-q80-third-q12-descent-field-normalization-v1",
        "status": "PASS_EXACT_NO_FACTORIZATION",
        "source": {
            "path": str(args.input.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "schema": source.get("schema"),
        },
        "checker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "q1": {
            "numerator_bits": q1.numerator.bit_length(),
            "denominator_bits": q1.denominator.bit_length(),
        },
        "q2": {
            "numerator_bits": q2.numerator.bit_length(),
            "denominator_bits": q2.denominator.bit_length(),
        },
        "reduced_q1_q2": {
            "numerator": str(product.numerator),
            "denominator": str(product.denominator),
            "numerator_bits": product.numerator.bit_length(),
            "denominator_bits": product.denominator.bit_length(),
            "numerator_is_square": numerator_square,
            "denominator_is_square": denominator_square,
            "numerator_square_root": str(numerator_root),
            "denominator_floor_square_root": str(denominator_root),
        },
        "field_normalization": {
            "old_generator": "omega, omega^2 = 16*q1*q2",
            "new_generator": "delta, delta^2 = denominator(q1*q2)",
            "omega_in_delta_basis": (
                "omega = 4*sqrt(numerator(q1*q2))*delta/denominator(q1*q2)"
            ),
            "delta_in_omega_basis": (
                "delta = denominator(q1*q2)*omega/"
                "(4*sqrt(numerator(q1*q2)))"
            ),
            "same_quadratic_field": True,
            "warning": (
                "This removes the rational radicand denominator exactly but does "
                "not by itself bound or reduce transformed equation coefficients."
            ),
        },
        "reproduce": (
            "python3 elkies-k3/scripts/"
            "audit_q80_third_q12_descent_field_normalization.py --write-artifact"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(encoded.encode()).hexdigest()

    if args.write_artifact:
        args.output.write_text(encoded)
        print(f"Q80Q12FIELD|artifact={args.output}|sha256={digest}|status=PASS_WRITE")
    elif args.check:
        assert args.output.read_text() == encoded
        print(f"Q80Q12FIELD|artifact={args.output}|sha256={digest}|status=PASS_CHECK")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
