#!/usr/bin/env sage
"""Materialize explicit 2-cover intersections of quadrics from cubic classes.

For ``E: y^2 = f(x)`` with monic cubic ``f`` and ``K = Q(theta)``, the
standard no-rational-2-torsion descent represents a class by ``alpha`` in
``K*/K*2`` and writes

    x - theta = alpha * (u + v theta + w theta^2)^2.

Equating the theta and theta^2 coefficients gives two affine quadrics.  Their
homogenisations in ``[u:v:w:z]`` are a genus-one 2-cover.  This script is
deliberately only a cover *builder*: square norm is necessary for a descent
class but does not certify local solubility, Selmer membership, or a rank
bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


PROTOCOL = "BNFFREE2COVER"
INPUT_SCHEMA = "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1"
OUTPUT_SCHEMA = "elliptic-curves.bnf-free-2cover-equations.v1"


def rational(value: str) -> QQ:
    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def multiply_mod_cubic(left, right, coefficients):
    """Multiply power-basis triples modulo a monic cubic, exactly."""
    a0, a1, a2, leading = coefficients
    if leading != 1:
        raise ValueError("the defining polynomial must be monic")
    product = [0] * 5
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product[left_index + right_index] += left_value * right_value

    # theta^3 = -a0 - a1 theta - a2 theta^2, and reduce theta^4 once more.
    return [
        product[0] - a0 * product[3] + a0 * a2 * product[4],
        product[1] - a1 * product[3] + (-a0 + a1 * a2) * product[4],
        product[2] - a2 * product[3] + (-a1 + a2**2) * product[4],
    ]


def cover_for(alpha, coefficients, ring):
    u, v, w, z = ring.gens()
    beta = [u, v, w]
    alpha_beta_square = multiply_mod_cubic(
        alpha,
        multiply_mod_cubic(beta, beta, coefficients),
        coefficients,
    )
    constant, theta_coefficient, theta_squared_coefficient = alpha_beta_square
    return {
        "theta_coefficient_plus_z_squared": str(theta_coefficient + z**2),
        "theta_squared_coefficient": str(theta_squared_coefficient),
        "constant_coefficient": str(constant),
        "affine_x_map": f"({constant})/z^2",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    record = json.loads(args.candidates.read_text())
    if not isinstance(record, dict) or record.get("schema") != INPUT_SCHEMA:
        raise ValueError("expected norm-filtered BNF-free squareclass candidates")
    coefficients = [rational(value) for value in record["field_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("the candidate field must have a monic cubic polynomial")

    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    covers = []
    for candidate in record.get("candidates", []):
        alpha = [rational(value) for value in candidate["generator_coefficients"]]
        if len(alpha) != 3:
            raise ValueError("a cubic candidate needs three power-basis coordinates")
        covers.append(
            {
                "label": str(candidate["label"]),
                "alpha_coefficients": [str(value) for value in alpha],
                "norm": str(candidate["norm"]),
                "quadrics": cover_for(alpha, coefficients, ring),
            }
        )

    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "explicit_two_covers_not_local_selmer_certificate",
        "field_polynomial_ascending": [str(value) for value in coefficients],
        "coordinate_names": ["u", "v", "w", "z"],
        "descent_equation": "x-theta = alpha*(u+v*theta+w*theta^2)^2",
        "covers": covers,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|covers={len(covers)}"
        "|status=EXPLICIT_COVERS_NOT_LOCAL_SELMER_CERTIFICATE",
        flush=True,
    )


if __name__ == "__main__":
    main()
