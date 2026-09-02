#!/usr/bin/env sage-python
"""Certify the rational G720-S0128 ``3I6+6I1`` marked source.

The input is the pinned 7-adic lift with free coordinate ``s6=10``.  Every
coordinate is rationally reconstructed, then checked over QQ against the full
55-equation marked system.  The script also certifies the fibre configuration,
split multiplicative supports, section component depths and Shioda height
matrix.  It proves a rank-19 Neron--Severi sublattice of determinant -720; an
independent Picard upper bound is deliberately outside this certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, PowerSeriesRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIFT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-s6-10-lift-v1.json"
)
DEFAULT_CLASSES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-ideal-source-isometries-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-source-qq-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize_polynomial(poly):
    return [str(value) for value in poly.list()]


def rational_reconstruct(coordinates, modulus):
    modulus = ZZ(modulus)
    return [ZZ(value).rational_reconstruction(modulus) for value in coordinates]


def formal_center(A, B, point, precision):
    base = A.parent()
    t = base.gen()
    shifted_A = base(A(t + point))
    shifted_B = base(B(t + point))
    node = -QQ(3) * shifted_B[0] / (QQ(2) * shifted_A[0])
    series_ring = PowerSeriesRing(QQ, "s", default_prec=precision + 3)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(8):
        center = (center + (-series_A / 3) / center) / 2
    assert (center**2 + series_A / 3).valuation() >= precision + 1
    return series_ring, center


def reversed_local(poly, weight, series_ring):
    u = series_ring.gen()
    return series_ring(
        sum(poly[index] * u ** (weight - index) for index in range(poly.degree() + 1))
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--classes", type=Path, default=DEFAULT_CLASSES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

lift_path = arguments.lift.resolve()
classes_path = arguments.classes.resolve()
output_path = arguments.output.resolve()
lift = json.loads(lift_path.read_text())
classes = json.loads(classes_path.read_text())

finite = lift["finite_precision_lift"]
assert lift["status"] == (
    "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
)
assert lift["prime"] == 7
assert finite["fixed_free_parameter_integer"] == 10
assert finite["achieved_precision_exponent"] >= 40
assert finite["all_55_residuals_zero_modulus"]
ideal_class = next(row for row in classes["classes"] if row["class_id"] == "G720-I001")
assert ideal_class["representative_source_id"] == "G720-S0128"
assert ideal_class["reduced_gram_row_count"] == 35

coordinates = rational_reconstruct(finite["coordinates_modulus"], finite["modulus"])
assert len(coordinates) == 46
a = coordinates[0:9]
b = coordinates[9:22]
p = coordinates[22:27]
q = coordinates[27:34]
r = coordinates[34:39]
s = coordinates[39:46]
assert s[-1] == 10

R = PolynomialRing(QQ, "t")
t = R.gen()
A, B = R(a), R(b)
X_P, Y_P = R(p), R(q)
X_Q, Y_Q = R(r), R(s)
D = 4 * A**3 + 27 * B**2

# Replay the full marked system over QQ, independently of rational reconstruction.
equations = []
equations.append(A[0] + 27)
equations.extend(D[index] for index in range(6))
equations.extend(D(t + 1)[index] for index in range(6))
equations.extend(D[index] for index in range(19, 25))
node_P = 2 * A * X_P + 3 * B
equations.extend(
    [
        node_P[0],
        Y_P[0],
        node_P(t + 1)[0],
        Y_P(1),
        node_P[12],
        node_P[11],
        node_P[10],
        Y_P[6],
        Y_P[5],
        Y_P[4],
    ]
)
equations.extend((Y_P**2 - X_P**3 - A * X_P - B)[index] for index in range(13))
equations.extend((Y_Q**2 - X_Q**3 - A * X_Q - B)[index] for index in range(13))
assert len(equations) == 55 and not any(equations)

# Exact semistable fibre configuration 3I6+6I1.
assert A.degree() == 8 and B.degree() == 12 and D.degree() == 18
fixed = t**6 * (t - 1) ** 6
assert D % fixed == 0
residual = (D // fixed).monic()
assert residual.degree() == 6 and residual.is_squarefree()
assert residual.gcd(t * (t - 1)) == 1
assert A.gcd(B) == 1

split_data = []
for support, av, bv in (
    ("0", A(0), B(0)),
    ("1", A(1), B(1)),
    ("infinity", A[8], B[12]),
):
    node = -3 * bv / (2 * av)
    tangent_square = QQ(3 * node)
    assert av and tangent_square.is_square()
    split_data.append(
        {
            "support": support,
            "kodaira": "I6",
            "node_x": str(node),
            "tangent_square": str(tangent_square),
            "tangent_slope": str(tangent_square.sqrt()),
            "split_over_Q": True,
        }
    )

# Recompute exact local component depths from the formal double-point centre.
precision = 6
zero_ring, zero_center = formal_center(A, B, QQ(0), precision)
one_ring, one_center = formal_center(A, B, QQ(1), precision)
infinity_ring = PowerSeriesRing(QQ, "u", default_prec=precision + 3)
infinity_A = reversed_local(A, 8, infinity_ring)
infinity_B = reversed_local(B, 12, infinity_ring)
infinity_center = infinity_ring(-QQ(3) * infinity_B[0] / (QQ(2) * infinity_A[0]))
for unused in range(8):
    infinity_center = (infinity_center + (-infinity_A / 3) / infinity_center) / 2
assert (infinity_center**2 + infinity_A / 3).valuation() >= precision + 1


def depths(X, Y):
    return [
        int(min((zero_ring(X(t)) - zero_center).valuation(), zero_ring(Y(t)).valuation())),
        int(min((one_ring(X(t + 1)) - one_center).valuation(), one_ring(Y(t + 1)).valuation())),
        int(
            min(
                (reversed_local(X, 4, infinity_ring) - infinity_center).valuation(),
                reversed_local(Y, 6, infinity_ring).valuation(),
            )
        ),
    ]


depth_P = depths(X_P, Y_P)
depth_Q = depths(X_Q, Y_Q)
assert depth_P == [1, 1, 3]
assert depth_Q == [0, 0, 0]

common = (X_P - X_Q).gcd(Y_P - Y_Q).monic()
assert common.degree() == 2 and common.is_squarefree() and common.gcd(D) == 1
component_correction_P = sum(QQ(depth * (6 - depth)) / 6 for depth in depth_P)
component_correction_Q = sum(QQ(depth * (6 - depth)) / 6 for depth in depth_Q)
height_P = QQ(4) - component_correction_P
height_Q = QQ(4) - component_correction_Q
height_cross = QQ(2) - common.degree()
height_gram = matrix(QQ, [[height_P, height_cross], [height_cross, height_Q]])
assert height_gram == matrix(QQ, [[QQ(5) / 6, 0], [0, 4]])
assert height_gram.det() == QQ(10) / 3

root_determinant = 6**3
ns_determinant = -QQ(root_determinant) * height_gram.det()
assert ns_determinant == -720

input_paths = (lift_path, classes_path)
payload = {
    "schema": "elkies-k3.golay-det720-3a5-source-qq.v1",
    "status": "PASS_EXACT_QQ_3I6_MW2_RANK19_SUBLATTICE_DET720",
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_golay_det720_3a5_source_qq.sage"
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "rational_reconstruction": {
        "prime": 7,
        "precision_exponent": finite["achieved_precision_exponent"],
        "fixed_free_coordinate": "s6=10",
        "all_46_coordinates_reconstructed": True,
        "all_55_equations_zero_over_Q": True,
        "maximum_numerator_or_denominator": int(
            max(max(abs(value.numerator()), value.denominator()) for value in coordinates)
        ),
    },
    "weierstrass_model": {
        "equation": "y^2=x^3+A(t)x+B(t)",
        "A_coefficients_low_to_high": serialize_polynomial(A),
        "B_coefficients_low_to_high": serialize_polynomial(B),
        "discriminant_core_factorization": str(D.factor()),
        "residual_I1_polynomial_monic": serialize_polynomial(residual),
        "fibre_profile": "3I6+6I1",
        "split_reducible_fibres": split_data,
    },
    "marked_sections": [
        {
            "name": "P",
            "X_coefficients_low_to_high": serialize_polynomial(X_P),
            "Y_coefficients_low_to_high": serialize_polynomial(Y_P),
            "component_depths_at_0_1_infinity": depth_P,
            "height": str(height_P),
        },
        {
            "name": "Q",
            "X_coefficients_low_to_high": serialize_polynomial(X_Q),
            "Y_coefficients_low_to_high": serialize_polynomial(Y_Q),
            "component_depths_at_0_1_infinity": depth_Q,
            "height": str(height_Q),
        },
    ],
    "section_pair": {
        "smooth_intersection_polynomial_monic": serialize_polynomial(common),
        "smooth_intersection_degree": int(common.degree()),
        "component_cross_correction": "0",
        "height_pairing": str(height_cross),
    },
    "lattice": {
        "trivial_root_type": "3A5",
        "trivial_root_rank": 15,
        "section_height_gram": [[str(value) for value in row] for row in height_gram.rows()],
        "section_height_determinant": str(height_gram.det()),
        "explicit_NS_sublattice_rank": 19,
        "explicit_NS_sublattice_determinant": int(ns_determinant),
        "matched_ideal_class": "G720-I001",
        "matched_source": "G720-S0128",
    },
    "proof_boundary": {
        "proved": (
            "The displayed rational Weierstrass model has split fibres 3I6+6I1 and "
            "the two displayed rational sections have height Gram diag(5/6,4), giving "
            "an explicit rank-19 Neron-Severi sublattice of determinant -720."
        ),
        "not_proved": (
            "This certificate alone does not prove geometric Picard rank exactly 19, "
            "MW saturation, identity with the abstract target fibration, a physical "
            "neighbour corridor, or a specialization rank jump."
        ),
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "GOLAY7203A5QQ|fibres=3I6+6I1|heights=5/6,4|"
    "NS_rank_at_least=19|NS_sublattice_det=-720|status=PASS"
)
