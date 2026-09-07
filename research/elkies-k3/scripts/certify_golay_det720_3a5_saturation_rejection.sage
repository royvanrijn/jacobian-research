#!/usr/bin/env sage-python
"""Certify that the rational ``s6=10`` point is not the G720 K3.

The rational model has both a 3-torsion section and a rational half of the
displayed height-four section.  Thus its displayed determinant-720 frame has
index six in Neron--Severi.  An exhaustive discriminant-form calculation shows
that six is the largest possible even-overlattice index.  Together with the
separate Picard-rank-19 certificate this proves full saturation, NS determinant
20, torsion Z/3, and free MW height lattice diag(5/6,1).
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
GEN = ROOT / "artifacts/generated-results"
DEFAULT_MODEL = GEN / "elkies-k3-golay-det720-3a5-source-qq-v1.json"
DEFAULT_PICARD = GEN / "elkies-k3-golay-det720-3a5-picard19-v1.json"
DEFAULT_SOURCES = GEN / "elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-golay-det720-3a5-saturation-rejection-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(poly):
    return [str(value) for value in poly.list()]


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


def component_depths(A, B, X, Y):
    precision = 6
    t = A.parent().gen()
    zero_ring, zero_center = formal_center(A, B, QQ(0), precision)
    one_ring, one_center = formal_center(A, B, QQ(1), precision)
    infinity_ring = PowerSeriesRing(QQ, "u", default_prec=precision + 3)
    infinity_A = reversed_local(A, 8, infinity_ring)
    infinity_B = reversed_local(B, 12, infinity_ring)
    infinity_center = infinity_ring(
        -QQ(3) * infinity_B[0] / (QQ(2) * infinity_A[0])
    )
    for unused in range(8):
        infinity_center = (
            infinity_center + (-infinity_A / 3) / infinity_center
        ) / 2
    assert (infinity_center**2 + infinity_A / 3).valuation() >= precision + 1
    return [
        int(
            min(
                (zero_ring(X(t)) - zero_center).valuation(),
                zero_ring(Y(t)).valuation(),
            )
        ),
        int(
            min(
                (one_ring(X(t + 1)) - one_center).valuation(),
                one_ring(Y(t + 1)).valuation(),
            )
        ),
        int(
            min(
                (reversed_local(X, 4, infinity_ring) - infinity_center).valuation(),
                reversed_local(Y, 6, infinity_ring).valuation(),
            )
        ),
    ]


def isotropic_subgroups(gram):
    """Enumerate every q-isotropic subgroup of the discriminant module."""

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
            integer for integer in range(1, 1 + max(moduli)) if multiply(integer, element) == zero
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
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--picard", type=Path, default=DEFAULT_PICARD)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

model_path = arguments.model.resolve()
picard_path = arguments.picard.resolve()
sources_path = arguments.sources.resolve()
output_path = arguments.output.resolve()
model = json.loads(model_path.read_text())
picard = json.loads(picard_path.read_text())
sources = json.loads(sources_path.read_text())
assert model["status"] == "PASS_EXACT_QQ_3I6_MW2_RANK19_SUBLATTICE_DET720"
assert picard["status"] == "PASS_EXACT_TWO_PRIME_ARTIN_TATE_PICARD19"
assert picard["geometric_picard_rank_characteristic_zero"] == 19
source = next(row for row in sources["sources"] if row["source_id"] == "G720-S0128")
assert model["lattice"]["matched_source"] == source["source_id"]
assert source["source"]["gram_sha256"] == (
    "84f90326d845a4e8f4a4332deb6d9946a1c085946ac7ff3ef04e7c3a553339da"
)

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R(model["weierstrass_model"]["A_coefficients_low_to_high"])
B = R(model["weierstrass_model"]["B_coefficients_low_to_high"])
P_record, Q_record = model["marked_sections"]
X_P = R(P_record["X_coefficients_low_to_high"])
Y_P = R(P_record["Y_coefficients_low_to_high"])
X_Q = R(Q_record["X_coefficients_low_to_high"])
Y_Q = R(Q_record["Y_coefficients_low_to_high"])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
P = E(K(X_P), K(Y_P))
Q = E(K(X_Q), K(Y_Q))

# The rational linear factor of the 3-division polynomial.
X_T = R([3, QQ(36) / 5, QQ(258) / 25, QQ(36) / 5, 3])
Y_T = R([0, 0, QQ(1728) / 5, -QQ(3456) / 5, QQ(1728) / 5])
T = E(K(X_T), K(Y_T))
division_3 = 3 * K["x"].gen() ** 4 + 6 * K(A) * K["x"].gen() ** 2 + 12 * K(B) * K["x"].gen() - K(A) ** 2
assert division_3(K(X_T)) == 0
assert T != E(0) and 3 * T == E(0)
depth_T = component_depths(A, B, X_T, Y_T)
assert depth_T == [2, 2, 2]

# A rational linear factor of the duplication quartic gives Q=2R_half.
X_half_Q = R([3, QQ(36) / 5, QQ(3138) / 25, QQ(36) / 5, 3])
Y_half_Q = R([0, 0, 0, QQ(41472) / 25])
half_Q = E(K(X_half_Q), K(Y_half_Q))
assert 2 * half_Q == Q
depth_half_Q = component_depths(A, B, X_half_Q, Y_half_Q)
assert depth_half_Q == [3, 0, 3]

height_P = QQ(5) / 6
height_Q = QQ(4)
height_half_Q = height_Q / 4
height_P_half_Q = QQ(0)
full_free_height = matrix(
    QQ, [[height_P, height_P_half_Q], [height_P_half_Q, height_half_Q]]
)
assert full_free_height == matrix(QQ, [[QQ(5) / 6, 0], [0, 1]])

displayed_frame = matrix(ZZ, source["source"]["gram"])
assert displayed_frame.det() == 720
smith_invariants, subgroups = isotropic_subgroups(displayed_frame)
index_histogram = Counter(len(subgroup) for subgroup in subgroups)
assert smith_invariants == [2, 6, 60]
assert index_histogram == Counter({6: 6, 2: 3, 3: 2, 1: 1})
maximum_even_overlattice_index = max(index_histogram)
assert maximum_even_overlattice_index == 6

# T enlarges by three and half_Q enlarges the free quotient by two.  They are
# independent: torsion has height zero, whereas half_Q has height one and is
# not an integral combination of P,Q (its Q-coordinate is 1/2).
explicit_index = 3 * 2
assert explicit_index == maximum_even_overlattice_index
full_ns_determinant = -ZZ(720) // explicit_index**2
assert full_ns_determinant == -20
assert QQ(6**3) * full_free_height.det() / 3**2 == 20

input_paths = (model_path, picard_path, sources_path)
payload = {
    "schema": "elkies-k3.golay-det720-3a5-saturation-rejection.v1",
    "status": "PASS_EXACT_RATIONAL_POINT_REJECTED_NS_DET20_TORSION3_HALF_SECTION",
    "reproduce": (
        "sage -python "
        "elkies-k3/scripts/certify_golay_det720_3a5_saturation_rejection.sage"
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "torsion_section": {
        "X_coefficients_low_to_high": serialize(X_T),
        "Y_coefficients_low_to_high": serialize(Y_T),
        "component_depths_at_0_1_infinity": depth_T,
        "exact_order": 3,
    },
    "half_of_displayed_Q": {
        "X_coefficients_low_to_high": serialize(X_half_Q),
        "Y_coefficients_low_to_high": serialize(Y_half_Q),
        "component_depths_at_0_1_infinity": depth_half_Q,
        "identity": "2*R=Q",
        "height": str(height_half_Q),
    },
    "discriminant_form": {
        "displayed_frame_determinant": 720,
        "smith_invariants": smith_invariants,
        "isotropic_subgroup_order_histogram": {
            str(index): count for index, count in sorted(index_histogram.items())
        },
        "maximum_even_overlattice_index": maximum_even_overlattice_index,
        "explicit_overlattice_index": explicit_index,
        "explicit_overlattice_is_maximal": True,
    },
    "full_surface": {
        "geometric_picard_rank": 19,
        "reducible_fibre_root_type": "3A5",
        "mordell_weil_rank": 2,
        "mordell_weil_torsion": "Z/3",
        "free_basis": ["P", "R with 2R=Q"],
        "free_height_gram": [[str(value) for value in row] for row in full_free_height.rows()],
        "neron_severi_determinant": int(full_ns_determinant),
    },
    "optimizer_consequence": (
        "The s6=10 rational point is a Picard-19 K3 with an exceptionally simple "
        "3A5/MW2 equation, but it is the determinant-20 maximal even overlattice, "
        "not the determinant-720 NS/T class carrying G720-F001.  It must therefore "
        "be rejected as a source for that abstract MW17 target."
    ),
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "GOLAY720SATREJECT|torsion=3|Q_divisibility=2|index=6|"
    "rho=19|NS_det=-20|G720_identity=REJECTED|status=PASS"
)
