#!/usr/bin/env python3
"""Confluent-coordinate audit for the normalized tagged GVC(3) cubic.

The tagged lift has

    Lambda = d_t*d_z + B(d_t,d_y),  P = z(t-y) + C(t,y),

where B and C are binary cubics.  The factor-compatible locus is C=(t-y)Q.
For a normalized C, its transverse coordinate is epsilon=C(1,1), so write

    C = (t-y)(t^2 + q*t*y + r*y^2) + epsilon*y^3.

This is an invertible coordinate change on all normalized binary cubics, not
an ansatz restriction.  It is the direct analogue of a tuned leading-term
cancellation: epsilon=0 is the confluent face, while epsilon records the
first transverse correction.

The script forms the exact diagonal tagged moments over GF(101), eliminates
the four coefficients of B together with q, r, epsilon, and records whether
the global normalized cubic chart is empty through the selected cutoff.
This is a finite-characteristic bounded calculation only.
"""

from __future__ import annotations

import json
import subprocess
import sys
from math import comb, factorial
from pathlib import Path

from sympy.polys.domains import GF
from sympy.polys.rings import ring


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "three_variable_gvc_confluent_tag.json"
MODULUS = 101

BinaryForm = dict[int, object]


def multiply(left: BinaryForm, right: BinaryForm, zero) -> BinaryForm:
    result: BinaryForm = {}
    for left_t, left_coefficient in left.items():
        for right_t, right_coefficient in right.items():
            exponent = left_t + right_t
            result[exponent] = result.get(exponent, zero) + left_coefficient * right_coefficient
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def tagged_moment(order: int, symbol: BinaryForm, profile: BinaryForm, zero):
    """Evaluate the exact diagonal formula for L=t-y."""
    linear = {1: 1, 0: -1}
    linear_powers = [{0: 1}]
    for _ in range(order):
        linear_powers.append(multiply(linear_powers[-1], linear, zero))

    symbol_power: BinaryForm = {0: 1}
    profile_power: BinaryForm = {0: 1}
    result = zero
    for channel_count in range(order + 1):
        if channel_count:
            symbol_power = multiply(symbol_power, symbol, zero)
            profile_power = multiply(profile_power, profile, zero)
        polynomial = multiply(linear_powers[order - channel_count], profile_power, zero)
        pairing = zero
        for derivative_t, operator_coefficient in symbol_power.items():
            target_t = derivative_t + order - channel_count
            target_y = 3 * channel_count - derivative_t
            coefficient = polynomial.get(target_t, zero)
            if coefficient:
                pairing += operator_coefficient * coefficient * factorial(target_t) * factorial(target_y)
        result += comb(order, channel_count) ** 2 * factorial(order - channel_count) * pairing
    return result


def singular_expression(polynomial) -> str:
    return str(polynomial).replace("**", "^")


def singular_summary(equations: list[str], variables: str) -> tuple[int, int, str]:
    script = (
        f"ring R={MODULUS},({variables}),dp;\n"
        "option(redSB);\n"
        f"ideal I={','.join(equations)};\n"
        "ideal G=std(I);\n"
        "print(size(G));\n"
        "print(dim(G));\n"
        "reduce(1,G);\n"
        "quit;\n"
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        raise RuntimeError(completed.stdout + completed.stderr)
    return int(lines[-3]), int(lines[-2]), lines[-1]


def main() -> None:
    full_chart = "--full-chart" in sys.argv[1:]
    if full_chart:
        coefficient_ring, a0, a1, a2, a3, q, r, epsilon = ring(
            "a0,a1,a2,a3,q,r,epsilon", GF(MODULUS)
        )
        variables = "a0,a1,a2,a3,q,r,epsilon"
        transverse_coordinate = "epsilon (free)"
    else:
        coefficient_ring, a0, a1, a2, a3, q, r = ring(
            "a0,a1,a2,a3,q,r", GF(MODULUS)
        )
        epsilon = coefficient_ring.one
        variables = "a0,a1,a2,a3,q,r"
        transverse_coordinate = "epsilon=1 slice"
    symbol = {3: a0, 2: a1, 1: a2, 0: a3}
    # C=(t-y)(t^2+q*t*y+r*y^2)+epsilon*y^3.
    profile = {
        3: coefficient_ring.one,
        2: q - coefficient_ring.one,
        1: r - q,
        0: epsilon - r,
    }
    moments = [
        singular_expression(tagged_moment(order, symbol, profile, coefficient_ring.zero))
        for order in range(1, 8)
    ]
    summaries = {
        str(cutoff): singular_summary(moments[:cutoff], variables)
        for cutoff in (5, 6, 7)
    }
    artifact = {
        "format": "three-variable-gvc-confluent-tag-v1",
        "field": f"GF({MODULUS})",
        "normalization": "coefficient of t^3 in C is 1",
        "computed_chart": transverse_coordinate,
        "coordinates": {
            "C": "(t-y)(t^2+q*t*y+r*y^2)+epsilon*y^3",
            "epsilon": "C(1,1)",
            "epsilon_zero_locus": "factor-compatible cubic profiles",
            "coordinate_change": "c1=q-1, c2=r-q, c3=epsilon-r",
        },
        "orders": {
            cutoff: {
                "basis_size": summary[0],
                "dimension": summary[1],
                "normal_form_of_1": summary[2],
            }
            for cutoff, summary in summaries.items()
        },
        "status": (
            "exact modular calculation; the default epsilon=1 slice is a "
            "bounded transversal probe, while --full-chart requests the "
            "full normalized-cubic elimination; no characteristic-zero claim"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    for cutoff, summary in summaries.items():
        print(f"m<={cutoff}: basis_size={summary[0]}, dimension={summary[1]}, NF(1)={summary[2]}")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
