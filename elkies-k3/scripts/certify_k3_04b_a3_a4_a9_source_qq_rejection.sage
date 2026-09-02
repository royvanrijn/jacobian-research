#!/usr/bin/env sage-python
"""Reconstruct and saturation-reject the first rational determinant-500 seed.

The nonsquare GF(7) marked branch at free coordinate ``m4=-20`` reconstructs
to a small rational ``I4+I5+I10+5I1`` model.  Its displayed height-5/2
section is five times a rational pole-zero section of height 1/10.  Exact
discriminant-form enumeration shows that index five is the largest possible
even overlattice, so the rank-19 primitive closure has determinant 20 rather
than 500.  This is a useful rational K3, but not the intended foundry class.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, PowerSeriesRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIFT = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-marked-gf7-nonsquare-m0-free-minus20-p40-v1.json"
)
DEFAULT_SOURCES = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-source-qq-rejection-v1.json"
)
SOURCE_ID = "K3-04b86146cc6b284b-S0160"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(poly):
    return [str(value) for value in poly.list()]


def formal_center(A, B, point, precision):
    t = A.parent().gen()
    shifted_A = A.parent()(A(t + point))
    shifted_B = A.parent()(B(t + point))
    series_ring = PowerSeriesRing(QQ, "s", default_prec=precision + 3)
    center = series_ring(-3 * shifted_B[0] / (2 * shifted_A[0]))
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


def component_depths(A, B, X, Y):
    t = A.parent().gen()
    zero_ring, zero_center = formal_center(A, B, QQ(0), 8)
    one_ring, one_center = formal_center(A, B, QQ(1), 8)
    infinity_ring = PowerSeriesRing(QQ, "u", default_prec=11)
    infinity_A = reversed_local(A, 8, infinity_ring)
    infinity_B = reversed_local(B, 12, infinity_ring)
    infinity_center = infinity_ring(-3 * infinity_B[0] / (2 * infinity_A[0]))
    for unused in range(8):
        infinity_center = (
            infinity_center + (-infinity_A / 3) / infinity_center
        ) / 2
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


def isotropic_subgroups(gram):
    smith, left, right = gram.smith_form()
    assert smith == left * gram * right
    diagonal = [abs(ZZ(smith[index, index])) for index in range(gram.nrows())]
    active = [index for index, value in enumerate(diagonal) if value > 1]
    moduli = [diagonal[index] for index in active]
    left_inverse = left.inverse()
    gram_inverse = gram.inverse()
    zero = (0,) * len(active)

    def add(a, b):
        return tuple((a[index] + b[index]) % moduli[index] for index in range(len(active)))

    def multiply(integer, a):
        answer = zero
        for unused in range(integer):
            answer = add(answer, a)
        return answer

    def is_isotropic(coordinates):
        representative = vector(ZZ, gram.nrows())
        for index, value in zip(active, coordinates):
            representative[index] = value
        dual_numerator = left_inverse * representative
        norm = dual_numerator * gram_inverse * dual_numerator
        return norm.denominator() == 1 and norm.numerator() % 2 == 0

    isotropic = {
        coordinates
        for coordinates in itertools.product(*(range(modulus) for modulus in moduli))
        if is_isotropic(coordinates)
    }

    def closure(subgroup, element):
        order = next(
            integer
            for integer in range(1, 1 + max(moduli))
            if multiply(integer, element) == zero
        )
        return frozenset(
            add(member, multiply(integer, element))
            for member in subgroup
            for integer in range(order)
        )

    subgroups = {frozenset((zero,))}
    changed = True
    while changed:
        changed = False
        for subgroup in list(subgroups):
            for element in isotropic - subgroup:
                enlarged = closure(subgroup, element)
                if enlarged.issubset(isotropic) and enlarged not in subgroups:
                    subgroups.add(enlarged)
                    changed = True
    return [int(value) for value in diagonal if value > 1], subgroups


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
lift_path = arguments.lift.resolve()
sources_path = arguments.sources.resolve()
output_path = arguments.output.resolve()
lift = json.loads(lift_path.read_text())
sources = json.loads(sources_path.read_text())
if lift["finite_precision_lift"].get("fixed_free_parameter_integer") != -20:
    raise ValueError("QQ reconstruction requires the m4=-20 lift")
modulus = ZZ(lift["finite_precision_lift"]["modulus"])
coordinates = [
    ZZ(value).rational_reconstruction(modulus)
    for value in lift["finite_precision_lift"]["coordinates_modulus"]
]

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R(coordinates[:9])
B = R(coordinates[9:22])
C = t + coordinates[22]
X_numerator = R(coordinates[23:30])
Y_numerator = R(coordinates[30:40])
section_residual = (
    Y_numerator**2
    - X_numerator**3
    - A * X_numerator * C**4
    - B * C**6
)
assert section_residual == 0

D = 4 * A**3 + 27 * B**2
residual, remainder = D.quo_rem(t**4 * (t - 1) ** 5)
assert remainder == 0 and residual.degree() == 5
assert residual(0) and residual(1) and residual.gcd(residual.derivative()).degree() == 0
split_nodes = []
for label, local_A, local_B in (
    ("0", A(0), B(0)),
    ("1", A(1), B(1)),
    ("infinity", A[8], B[12]),
):
    node = -3 * local_B / (2 * local_A)
    assert QQ(3 * node).is_square()
    split_nodes.append({"support": label, "node_x": str(node), "tangent_square": str(3 * node)})

E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
P = E(K(X_numerator / C**2), K(Y_numerator / C**3))
S = PolynomialRing(K, "x")
x = S.gen()
torsion_factor_degrees = {}
for order, polynomial in (
    (2, S(E.two_division_polynomial(x))),
    (3, S(E.division_polynomial(3, x))),
    (5, S(E.division_polynomial(5, x))),
):
    torsion_factor_degrees[str(order)] = [
        int(factor.degree()) for factor, multiplicity in polynomial.factor()
    ]
assert torsion_factor_degrees == {"2": [3], "3": [4], "5": [12]}

preimage_factor_degrees = {}
linear_fifth = None
for order in (2, 5):
    preimage = S(E._multiple_x_numerator(order, x)) - P[0] * S(
        E._multiple_x_denominator(order, x)
    )
    factors = list(preimage.factor())
    preimage_factor_degrees[str(order)] = [int(factor.degree()) for factor, exponent in factors]
    if order == 5:
        linear_fifth = next(factor for factor, exponent in factors if factor.degree() == 1)
assert preimage_factor_degrees == {"2": [4], "5": [1, 24]}
X_R = K(-linear_fifth[0] / linear_fifth[1])
Y_R = (X_R**3 + K(A) * X_R + K(B)).sqrt()
R_point = E(X_R, Y_R)
if 5 * R_point != P:
    R_point = -R_point
assert 5 * R_point == P
X_R_poly = R(X_R)
Y_R_poly = R(R_point[1])
depths_R = component_depths(A, B, X_R_poly, Y_R_poly)
assert depths_R == [2, 1, 3]
height_R = QQ(4) - (QQ(1) + QQ(4) / 5 + QQ(21) / 10)
assert height_R == QQ(1) / 10 and 25 * height_R == QQ(5) / 2

source = next(row for row in sources["sources"] if row["source_id"] == SOURCE_ID)["source"]
displayed_frame = matrix(ZZ, source["gram"])
smith_invariants, subgroups = isotropic_subgroups(displayed_frame)
index_histogram = Counter(len(subgroup) for subgroup in subgroups)
assert smith_invariants == [5, 5, 20]
assert index_histogram == Counter({5: 6, 1: 1})
assert max(index_histogram) == 5
assert source["root_determinant"] * height_R == 20

payload = {
    "schema": "elkies-k3.k3-04b-a3-a4-a9-source-qq-rejection.v1",
    "status": "PASS_EXACT_QQ_POINT_REJECTED_PRIMITIVE_CLOSURE_DET20_FIFTH_ROOT",
    "inputs": {
        relative(lift_path): digest(lift_path),
        relative(sources_path): digest(sources_path),
    },
    "weierstrass_model": {
        "equation": "y^2=x^3+A(t)*x+B(t)",
        "A_coefficients_low_to_high": serialize(A),
        "B_coefficients_low_to_high": serialize(B),
        "discriminant_factorization": (
            "65536*t^4*(t-1)^5*(t^5-t^4/2+3*t^3/4+73*t^2/8+9*t/4-81/8)"
        ),
        "split_fibres": ["I4", "I5", "I10"],
        "residual_squarefree_I1_fibres": 5,
        "split_node_data": split_nodes,
    },
    "displayed_section": {
        "C_coefficients_low_to_high": serialize(C),
        "X_numerator_coefficients_low_to_high": serialize(X_numerator),
        "Y_numerator_coefficients_low_to_high": serialize(Y_numerator),
        "component_depths_at_I4_I5_I10": [2, 0, 5],
        "height": "5/2",
        "pole_order": 1,
    },
    "saturation": {
        "torsion_division_polynomial_factor_degrees": torsion_factor_degrees,
        "displayed_section_preimage_factor_degrees": preimage_factor_degrees,
        "rational_fifth_root": {
            "X_coefficients_low_to_high": serialize(X_R_poly),
            "Y_coefficients_low_to_high": serialize(Y_R_poly),
            "component_depths_at_I4_I5_I10": depths_R,
            "height": "1/10",
            "pole_order": 0,
            "exact_relation": "5*R=P",
        },
        "displayed_frame_smith_invariants": smith_invariants,
        "isotropic_subgroup_order_histogram": {
            str(key): value for key, value in sorted(index_histogram.items())
        },
        "maximum_even_overlattice_index": 5,
        "rank19_primitive_closure_determinant": 20,
    },
    "proof_boundary": {
        "proved": (
            "All forty coordinates reconstruct over Q and satisfy the section and "
            "fibre identities exactly. The three reducible fibres are split and "
            "the residual quintic is squarefree. The displayed section is exactly "
            "five times the rational height-1/10 section. Index five is maximal "
            "among even overlattices of the displayed determinant-500 frame, so "
            "its rank-19 primitive closure has determinant 20."
        ),
        "not_proved": (
            "The full geometric Picard rank of this rational K3 is not computed. "
            "The rejection needs only the explicit index-five enlargement. No "
            "rational point in the intended primitive determinant-500 locus or "
            "neighbour corridor is constructed."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_04b_a3_a4_a9_source_qq_rejection.sage"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("determinant-500 QQ saturation rejection is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print(
    "K304BEQQREJECT|fibres=I4+I5+I10+5I1|P_divisibility=5|"
    "primitive_closure_det=20|status=PASS",
    flush=True,
)
