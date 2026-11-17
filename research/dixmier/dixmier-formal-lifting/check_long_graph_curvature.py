#!/usr/bin/env python3
"""Exact Weyl-ordered Long graph curvature; Fraction arithmetic only.

The canonical producer's long_images() is imported solely as polynomial
INPUT. No producer arithmetic, derivatives, Poisson brackets or PBW
multiplication are used to verify identities here.

Runs write a fresh temporary receipt by default. They do not overwrite
the canonical receipt or any proof/catalogue file.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
from math import factorial
from pathlib import Path
import platform
import tempfile

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent / "dixmier-audit" / "verify_long_pbw.py"
ZERO = (0, 0, 0, 0)
VARIABLES = ("X", "Q", "P", "Z")


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return {e: c for e, c in result.items() if c}


def scale(polynomial, scalar):
    scalar = Fraction(scalar)
    return {e: scalar * c for e, c in polynomial.items() if scalar * c}


def multiply(left, right):
    result = {}
    for e, a in left.items():
        for f, b in right.items():
            exponent = tuple(e[i] + f[i] for i in range(4))
            result[exponent] = result.get(exponent, Fraction(0)) + a * b
    return {e: c for e, c in result.items() if c}


def derivative(polynomial, counts):
    result = {}
    for exponent, coefficient in polynomial.items():
        if any(exponent[i] < counts[i] for i in range(4)):
            continue
        new_exponent = tuple(exponent[i] - counts[i] for i in range(4))
        multiplier = 1
        for i in range(4):
            for j in range(counts[i]):
                multiplier *= exponent[i] - j
        result[new_exponent] = coefficient * multiplier
    return {e: c for e, c in result.items() if c}


def bidifferential(left, right, order):
    """Pi^order, Pi=dP tensor dX-dX tensor dP+dZ tensor dQ-dQ tensor dZ."""
    result = {}
    for a in range(order + 1):
        for b in range(order - a + 1):
            for c in range(order - a - b + 1):
                d = order - a - b - c
                coefficient = Fraction(
                    factorial(order) * (-1) ** (b + d),
                    factorial(a) * factorial(b) * factorial(c) * factorial(d),
                )
                left_derivative = derivative(left, (b, d, a, c))
                right_derivative = derivative(right, (a, c, b, d))
                result = add(result, scale(
                    multiply(left_derivative, right_derivative), coefficient))
    return result


def evaluate(polynomial, point):
    result = Fraction(0)
    for exponent, coefficient in polynomial.items():
        term = coefficient
        for coordinate, power in zip(point, exponent):
            term *= Fraction(coordinate) ** power
        result += term
    return result


def rows(polynomial):
    return [[list(e), str(c)] for e, c in sorted(polynomial.items())]


def digest(polynomial):
    return hashlib.sha256(
        json.dumps(rows(polynomial), separators=(",", ":")).encode()
    ).hexdigest()


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    specification = importlib.util.spec_from_file_location("canonical_long_input", INPUT)
    require(specification is not None and specification.loader is not None,
            "canonical Long input can be loaded")
    producer = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(producer)
    raw_input = producer.long_images()  # INPUT ONLY; no other producer function is called.
    require(len(raw_input) == 4, "four full global input polynomials")
    images = []
    for polynomial in raw_input:
        converted = {tuple(e): Fraction(c) for e, c in polynomial.items() if c}
        require(all(len(e) == 4 and all(isinstance(n, int) and n >= 0 for n in e)
                    for e in converted), "actual polynomial source support")
        images.append(converted)
    labels = ("R", "T", "D", "S")
    full = dict(zip(labels, images))
    source = [
        {tuple(int(i == j) for i in range(4)): Fraction(1)}
        for j in range(4)
    ]

    # Positive canonical identity graph: all six relations are exact; all higher
    # bidifferentials vanish because each coordinate has ordinary degree one.
    identity_controls = []
    frame_controls = []
    for i, j in itertools.combinations(range(4), 2):
        expected = bidifferential(source[i], source[j], 1)
        identity_e3 = bidifferential(source[i], source[j], 3)
        require(not identity_e3, "identity graph has zero cubic Moyal curvature")
        identity_controls.append({
            "pair": [VARIABLES[i], VARIABLES[j]],
            "Poisson": rows(expected), "Pi3_zero": True,
            "all_higher_orders_zero_by_linear_degree": True,
        })
        actual = bidifferential(images[i], images[j], 1)
        require(actual == expected, "full global classical relation " + labels[i] + labels[j])
        frame_controls.append({
            "pair": [labels[i], labels[j]], "Poisson": rows(actual), "status": "PASS",
        })

    curvature = {}
    point = (1, 0, 0, 0)
    for left, right, expected_count, expected_value in (
        ("S", "T", 30, Fraction(199, 2)),
        ("D", "R", 30, Fraction(-27, 2)),
    ):
        # [f,g]_Moyal = h*Pi(f,g) + h^3*Pi^3(f,g)/24 + O(h^5).
        e3 = scale(bidifferential(full[left], full[right], 3), Fraction(1, 24))
        require(len(e3) == expected_count, "complete E3 support " + left + right)
        require(evaluate(e3, point) == expected_value, "nonzero E3 witness " + left + right)
        require(e3, "nonzero source curvature " + left + right)
        curvature[left + right] = {
            "definition": "E3=Pi^3(f,g)/24 for [f,g]_Moyal-h*{f,g}",
            "terms": len(e3), "polynomial": rows(e3),
            "polynomial_sha256": digest(e3),
            "point_source_order_X_Q_P_Z": list(point),
            "value": str(evaluate(e3, point)),
        }

    # R positional and S of momentum order two: the cubic Moyal term is zero,
    # although uncorrected normal PBW ordering has a different lower defect.
    require(not bidifferential(full["R"], full["S"], 3),
            "RS has zero cubic Weyl-order curvature")

    checker_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    input_hash = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    receipt = {
        "claim": "Canonical full Long graph has nonzero cubic Weyl curvature",
        "status": "PASS",
        "environment": {"python": platform.python_version(),
                        "arithmetic": "fractions.Fraction; complete finite coefficient support"},
        "canonical_input": {
            "path": "dixmier/dixmier-audit/verify_long_pbw.py",
            "callable": "long_images",
            "returned_order": list(labels), "source_order": list(VARIABLES),
            "sha256": input_hash,
            "use_boundary": "Input polynomial construction only; verifier algebra is independent",
        },
        "input_support": {
            label: {"terms": len(poly), "degree": max(map(sum, poly)),
                    "polynomial_sha256": digest(poly)}
            for label, poly in zip(labels, images)
        },
        "all_six_classical_full_frame_relations": frame_controls,
        "identity_graph_positive_controls": identity_controls,
        "weyl_curvature": curvature,
        "RS_Pi3_zero": True,
        "hashes": {
            "check_long_graph_curvature.py": checker_hash,
            "../dixmier-audit/verify_long_pbw.py": input_hash,
        },
        "coverage_boundary": (
            "Exact canonical global polynomial inputs, six Poisson identities, "
            "positive identity-graph control and complete Pi^3/24 curvature "
            "coefficients for ST/DR. This does not verify higher corrected "
            "quantum jets, finite graph quantization, source projectivity, "
            "BRST truncation or nonexistence of corrected finite lifts. "
            "The unchanged-constraint h-flatness and triangular arguments "
            "are separate written implications."
        ),
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    output = args.output_dir
    if output is None:
        output = Path(tempfile.mkdtemp(prefix="long-graph-curvature-"))
    output.mkdir(parents=True, exist_ok=True)
    destination = output / "long_graph_curvature_controls.json"
    destination.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: six global Long Poisson relations; identity graph; cubic Weyl curvature")
    print("ST: 30 terms; E3(1,0,0,0)=199/2")
    print("DR: 30 terms; E3(1,0,0,0)=-27/2")
    print("Receipt:", destination)


if __name__ == "__main__":
    main()
