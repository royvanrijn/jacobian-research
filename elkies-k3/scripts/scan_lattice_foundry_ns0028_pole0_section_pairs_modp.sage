#!/usr/bin/env sage-python
"""Exhaust the two component-adapted pole-zero sections on NS0028 models.

For ``P``, node incidence at the finite I3 and I7 supports leaves three free
coefficients in ``X``.  For ``Q``, node incidence at the I8 support at infinity
leaves four.  The section equation determines ``Y`` up to sign, so exact
polynomial-square tests exhaust both charts without a tensor scan.  Candidate
pairs are retained only when their intersections occur on smooth fibres and
their exact intersection number gives the required off-diagonal Shioda height
``+/-1``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0028-source-ansatz-mod5.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod5.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def formal_center(A, B, point, precision):
    field = A.base_ring()
    base = A.parent()
    t = base.gen()
    shifted_A = base(A(t + point))
    shifted_B = base(B(t + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 2)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(6):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("formal center did not converge")
    return series_ring, center, node


def reversed_local(poly, weight, series_ring):
    """Return ``u^weight poly(1/u)`` in the local parameter at infinity."""
    return series_ring(
        sum(poly[index] * series_ring.gen() ** (weight - index) for index in range(poly.degree() + 1))
    )


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

input_path = arguments.input.resolve()
payload = json.loads(input_path.read_text())
if payload["schema"] != "elkies-k3.lattice-foundry-ns0028-source-ansatz-modp.v1":
    raise ValueError("unexpected NS0028 fibre-ansatz schema")
if payload["accounting"]["stored_examples"] != payload["accounting"]["squarefree_examples_with_signs"]:
    raise ValueError("section-pair scan requires every squarefree fibre model")
prime = int(payload["prime"])
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("--quadratic-twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

records = []
for example_index, example in enumerate(payload["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    series_zero, center_zero, node_zero = formal_center(A, B, field.zero(), 3)
    series_one, center_one, node_one = formal_center(A, B, field.one(), 7)

    infinity_ring = PowerSeriesRing(field, "u", default_prec=10)
    local_A_infinity = reversed_local(A, 8, infinity_ring)
    local_B_infinity = reversed_local(B, 12, infinity_ring)
    node_infinity = -field(3) * local_B_infinity[0] / (
        field(2) * local_A_infinity[0]
    )
    center_infinity = infinity_ring(node_infinity)
    for unused in range(6):
        center_infinity = (
            center_infinity + (-local_A_infinity / 3) / center_infinity
        ) / 2
    if (center_infinity**2 + local_A_infinity / 3).valuation() < 9:
        raise ArithmeticError("infinity formal center did not converge")

    # P passes through the finite nodes.  Interpolate the required two values
    # and add t(t-1) times an arbitrary quadratic.
    p_remainder = ring(node_zero + (node_one - node_zero) * t)
    p_sections = []
    for q_values in itertools.product(field, repeat=3):
        X = p_remainder + t * (t - 1) * ring(q_values)
        for Y in polynomial_roots(X**3 + A * X + B):
            shifted_X_zero = series_zero(ring(X(t)))
            shifted_Y_zero = series_zero(ring(Y(t)))
            shifted_X_one = series_one(ring(X(t + 1)))
            shifted_Y_one = series_one(ring(Y(t + 1)))
            depth_zero = min(
                int((shifted_X_zero - center_zero).valuation()),
                int(shifted_Y_zero.valuation()),
            )
            depth_one = min(
                int((shifted_X_one - center_one).valuation()),
                int(shifted_Y_one.valuation()),
            )
            local_X_infinity = reversed_local(X, 4, infinity_ring)
            local_Y_infinity = reversed_local(Y, 6, infinity_ring)
            smooth_infinity = not (
                local_X_infinity[0] == node_infinity
                and local_Y_infinity[0] == 0
            )
            if (depth_zero, depth_one) != (1, 1) or not smooth_infinity:
                continue
            p_sections.append(
                {
                    "X": X,
                    "Y": Y,
                    "component_depths_at_I3_I7": [depth_zero, depth_one],
                }
            )

    # Q passes through the node at infinity, so only its leading X
    # coefficient is fixed.
    q_sections = []
    for lower_values in itertools.product(field, repeat=4):
        X = ring(lower_values) + node_infinity * t**4
        for Y in polynomial_roots(X**3 + A * X + B):
            local_X_infinity = reversed_local(X, 4, infinity_ring)
            local_Y_infinity = reversed_local(Y, 6, infinity_ring)
            depth_infinity = min(
                int((local_X_infinity - center_infinity).valuation()),
                int(local_Y_infinity.valuation()),
            )
            smooth_zero = not (X(0) == node_zero and Y(0) == 0)
            smooth_one = not (X(1) == node_one and Y(1) == 0)
            if depth_infinity != 1 or not smooth_zero or not smooth_one:
                continue
            q_sections.append(
                {
                    "X": X,
                    "Y": Y,
                    "component_depth_at_I8": depth_infinity,
                }
            )

    pairs = []
    for p_index, P in enumerate(p_sections):
        for q_index, Q in enumerate(q_sections):
            common = (P["X"] - Q["X"]).gcd(P["Y"] - Q["Y"])
            if common.gcd(discriminant_core).degree() != 0:
                continue
            intersection = int(common.degree())
            height_pairing = 2 - intersection
            if height_pairing not in (-1, 1):
                continue
            pairs.append(
                {
                    "P_index": p_index,
                    "Q_index": q_index,
                    "intersection_on_smooth_fibres": intersection,
                    "shioda_height_pairing": height_pairing,
                    "height_gram": [
                        ["52/21", str(height_pairing)],
                        [str(height_pairing), "25/8"],
                    ],
                    "mw_regulator": "283/42",
                    "implied_NS_determinant": 1132,
                }
            )

    def serialize_section(section):
        return {
            key: value
            for key, value in section.items()
            if key not in ("X", "Y")
        } | {
            "X_coefficients_low_to_high": [int(value) for value in section["X"]],
            "Y_coefficients_low_to_high": [int(value) for value in section["Y"]],
        }

    records.append(
        {
            "example_index": example_index,
            "P_X_polynomials_scanned": prime**3,
            "Q_X_polynomials_scanned": prime**4,
            "P_marked_section_count": len(p_sections),
            "Q_marked_section_count": len(q_sections),
            "P_sections": [serialize_section(section) for section in p_sections],
            "Q_sections": [serialize_section(section) for section in q_sections],
            "marked_pair_count": len(pairs),
            "marked_pairs": pairs,
        }
    )

total_pairs = sum(record["marked_pair_count"] for record in records)
output = {
    "schema": "elkies-k3.lattice-foundry-ns0028-pole0-section-pairs-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_WITH_MARKED_SECTION_PAIRS"
        if total_pairs
        else "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_EMPTY_MARKED_PAIR_CHART"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "input": {
        "artifact": relative(input_path),
        "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
    },
    "scope": {
        "stored_fibre_models": len(records),
        "all_squarefree_fibre_models_stored": True,
        "P_X_polynomials_per_model": prime**3,
        "Q_X_polynomials_per_model": prime**4,
        "all_polynomial_Y_square_roots_retained": True,
        "fibre_ansatz_scan_exhausted": bool(payload["scan"]["exhausted"]),
        "pair_intersections_restricted_to_smooth_fibres": True,
    },
    "accounting": {
        "total_P_X_polynomials_scanned": len(records) * prime**3,
        "total_Q_X_polynomials_scanned": len(records) * prime**4,
        "models_with_P": sum(bool(record["P_marked_section_count"]) for record in records),
        "models_with_Q": sum(bool(record["Q_marked_section_count"]) for record in records),
        "models_with_marked_pairs": sum(bool(record["marked_pair_count"]) for record in records),
        "total_P_sections": sum(record["P_marked_section_count"] for record in records),
        "total_Q_sections": sum(record["Q_marked_section_count"] for record in records),
        "total_marked_pairs": total_pairs,
    },
    "models": records,
    "proof_boundary": {
        "proved": (
            "For every stored fibre model, both component-adapted pole-zero X "
            "charts and all polynomial Y square roots are exhausted. Retained "
            "pairs have exact component depths, smooth-fibre intersection number, "
            "MW height Gram, regulator, and NS determinant."
        ),
        "not_proved": (
            "A finite-field marked pair is not a rational source family, a "
            "characteristic-zero lift, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage"
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = arguments.output.resolve()
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0028 pole-zero section-pair scan is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0028POLE0PAIRS|"
    f"models={len(records)}|P={output['accounting']['total_P_sections']}|"
    f"Q={output['accounting']['total_Q_sections']}|pairs={total_pairs}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
