#!/usr/bin/env sage
"""Construct the fixed-corridor terminal reverse-neighbour section over QQ.

The preferred q12/orbit5867 endpoint is an exact rootless Jacobian with a
seventeen-section basis integrally identified with pinned R17.  The historical
fixed corridor has a certified final transition

    A1/MW16 --q6 orbit2247--> rootless/MW17.

Invert that transition, compose the two pinned-R17 markings, and construct the
section P for the reverse fibre

    F_A1 = O + P - 2 F,       P.O = 6.

Only exact elliptic group law and integral matrix arithmetic are used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
BRIDGE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
HEIGHTS = LOCAL / "q12o5867-rootless-height-basis-qq.json"
SECTIONS = LOCAL / "q12o5867-rootless-selected-basis-qq.json"
MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = LOCAL / "fixed-final-a1-horizontal-from-q12-endpoint-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(M):
    return [[int(value) for value in row] for row in M.rows()]


def coeffs(poly):
    return [str(value) for value in poly.list()]


def maximum_rational_bits(polys):
    answer = 0
    for poly in polys:
        for value in poly:
            value = QQ(value)
            answer = max(
                answer,
                abs(value.numerator()).nbits(),
                value.denominator().nbits(),
            )
    return int(answer)


started = time.monotonic()
bridge = read_json(BRIDGE)
heights = read_json(HEIGHTS)
lifted = read_json(SECTIONS)
model = read_json(MODEL)
assert bridge["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert heights["status"] == "PASS_EXACT_QQ_Q12O5867_ROOTLESS_RANK17_HEIGHT_BASIS_PINNED"
assert lifted["status"] == "PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS"
assert model["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"

# Rows of A are the pinned-rootless basis expressed in the final A1 basis.
A = matrix(ZZ, bridge["final_a1_to_pinned_r17_transition"])
assert abs(A.det()) == 1
Ainv = A.inverse().change_ring(ZZ)
e = identity_matrix(ZZ, 19)
parent_fibre_pinned = vector(ZZ, e.row(0) * Ainv)
parent_zero_pinned = vector(ZZ, (-e.row(0) + e.row(1)) * Ainv)

pinned = load_matrix(PINNED)
w = vector(ZZ, parent_fibre_pinned[2:])
assert list(parent_fibre_pinned[:2]) == [4, 2]
assert w * pinned * w == 16

# A rootless section with MW vector w is [7,1,w], hence
# [4,2,w] = O + [7,1,w] - 2F and P.O=(16-4)/2=6.
target_section_class = vector(ZZ, [7, 1] + list(w))
assert parent_fibre_pinned == vector(ZZ, [-1, 1] + [0] * 17) + target_section_class - 2 * vector(ZZ, [1, 0] + [0] * 17)

# Rows of C are the exact q12 endpoint basis sections in pinned MW coords.
C = matrix(ZZ, heights["basis_to_pinned_rank17"])
assert C.det() == -1
word = vector(ZZ, w * C.inverse().change_ring(ZZ))
assert word * C == w

R = PolynomialRing(QQ, "u")
K = R.fraction_field()
Acurve = R([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
Bcurve = R([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])


def point_from_record(record):
    section = record["section"]
    x = K(R([QQ(value) for value in section["x_coefficients_low_to_high"]]))
    y = K(R([QQ(value) for value in section["y_coefficients_low_to_high"]]))
    assert y ** 2 == x ** 3 + K(Acurve) * x + K(Bcurve)
    return x, y


def neg(P):
    return None if P is None else (P[0], -P[1])


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1 ** 2 + K(Acurve)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope ** 2 - x1 - x2
    y3 = slope * (x1 - x3) - y1
    assert y3 ** 2 == x3 ** 3 + K(Acurve) * x3 + K(Bcurve)
    return x3, y3


def multiply(P, coefficient):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return neg(multiply(P, -coefficient))
    answer = None
    addend = P
    while coefficient:
        if coefficient & 1:
            answer = add(answer, addend)
        coefficient >>= 1
        if coefficient:
            addend = add(addend, addend)
    return answer


points = [point_from_record(record) for record in lifted["sections"]]


def compose(section_word):
    summands = [
        multiply(P, coefficient)
        for P, coefficient in zip(points, section_word)
        if coefficient
    ]
    # Pairwise reduction keeps intermediate rational functions smaller than a
    # fixed left-associated sum.
    while len(summands) > 1:
        next_level = []
        for index in range(0, len(summands), 2):
            if index + 1 == len(summands):
                next_level.append(summands[index])
            else:
                next_level.append(add(summands[index], summands[index + 1]))
        summands = next_level
    return summands[0]


target = compose(word)
assert target is not None
x, y = target
assert y ** 2 == x ** 3 + K(Acurve) * x + K(Bcurve)

xn, xd = R(x.numerator()), R(x.denominator())
yn, yd = R(y.numerator()), R(y.denominator())
assert xd.is_square() and yd == xd.sqrt() ** 3
Z = xd.sqrt()
assert Z.degree() == 6
assert [xn.degree(), xd.degree()] == [16, 12]
assert [yn.degree(), yd.degree()] == [24, 18]

# The historical A1 zero is another short exact endpoint word.  It has class
# [2,1,w0], hence P.O=1 and height six on the rootless model.
zero_w = vector(ZZ, parent_zero_pinned[2:])
assert zero_w * pinned * zero_w == 6
zero_word = vector(ZZ, zero_w * C.inverse().change_ring(ZZ))
assert zero_word * C == zero_w
fixed_zero = compose(zero_word)
assert fixed_zero is not None
zero_x, zero_y = fixed_zero
assert zero_y ** 2 == zero_x ** 3 + K(Acurve) * zero_x + K(Bcurve)
zero_xn, zero_xd = R(zero_x.numerator()), R(zero_x.denominator())
zero_yn, zero_yd = R(zero_y.numerator()), R(zero_y.denominator())
assert zero_xd.is_square() and zero_yd == zero_xd.sqrt() ** 3
zero_Z = zero_xd.sqrt()
assert zero_Z.degree() == 1

# The effective nonidentity A1 component is the negative of the stored simple
# root row after transport to the pinned rootless frame.  It is a polynomial
# P.O=0 section on the rootless source and becomes vertical for F_A1.
effective_component_class = -vector(ZZ, e.row(2) * Ainv)
assert list(effective_component_class[:2]) == [1, 1]
component_w = vector(ZZ, effective_component_class[2:])
assert component_w * pinned * component_w == 4
component_word = vector(ZZ, component_w * C.inverse().change_ring(ZZ))
component_point = compose(component_word)
assert component_point is not None
component_x, component_y = component_point
assert component_y ** 2 == component_x ** 3 + K(Acurve) * component_x + K(Bcurve)
component_xn, component_xd = R(component_x.numerator()), R(component_x.denominator())
component_yn, component_yd = R(component_y.numerator()), R(component_y.denominator())
assert component_xd.degree() == component_yd.degree() == 0

payload = {
    "schema": "elkies-k3.fixed-final-a1-horizontal-from-q12-endpoint-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_HORIZONTAL",
    "fixed_edge": {
        "forward": "A1/MW16 --q6 orbit2247--> rootless/MW17",
        "reverse_source": "exact q12/orbit5867 rootless Jacobian",
        "parent_fibre_in_pinned_rootless_coordinates": list(map(int, parent_fibre_pinned)),
        "parent_zero_in_pinned_rootless_coordinates": list(map(int, parent_zero_pinned)),
        "target_section_class_in_pinned_coordinates": list(map(int, target_section_class)),
        "target_section_word_in_exact_q12_basis": list(map(int, word)),
        "fixed_A1_zero_word_in_exact_q12_basis": list(map(int, zero_word)),
        "effective_A1_component_word_in_exact_q12_basis": list(map(int, component_word)),
        "identity": "F_A1 = O + P - 2F",
        "P_dot_O": 6,
        "height": 16,
    },
    "section": {
        "x_numerator_coefficients_low_to_high": coeffs(xn),
        "x_denominator_coefficients_low_to_high": coeffs(xd),
        "y_numerator_coefficients_low_to_high": coeffs(yn),
        "y_denominator_coefficients_low_to_high": coeffs(yd),
        "Z_coefficients_low_to_high": coeffs(Z),
        "degrees_X_Y_Z": [int(xn.degree()), int(yn.degree()), int(Z.degree())],
        "degrees_x_numerator_denominator": [int(xn.degree()), int(xd.degree())],
        "degrees_y_numerator_denominator": [int(yn.degree()), int(yd.degree())],
        "maximum_rational_bits": maximum_rational_bits((xn, xd, yn, yd)),
        "literal_weierstrass_identity": True,
    },
    "fixed_A1_zero_on_rootless_source": {
        "x_numerator_coefficients_low_to_high": coeffs(zero_xn),
        "x_denominator_coefficients_low_to_high": coeffs(zero_xd),
        "y_numerator_coefficients_low_to_high": coeffs(zero_yn),
        "y_denominator_coefficients_low_to_high": coeffs(zero_yd),
        "Z_coefficients_low_to_high": coeffs(zero_Z),
        "P_dot_O": 1,
        "height": 6,
        "maximum_rational_bits": maximum_rational_bits((zero_xn, zero_xd, zero_yn, zero_yd)),
        "literal_weierstrass_identity": True,
    },
    "effective_A1_component_on_rootless_source": {
        "class_in_pinned_coordinates": list(map(int, effective_component_class)),
        "x_coefficients_low_to_high": coeffs(component_xn / component_xd[0]),
        "y_coefficients_low_to_high": coeffs(component_yn / component_yd[0]),
        "P_dot_O": 0,
        "height": 4,
        "maximum_rational_bits": maximum_rational_bits((component_xn, component_yn)),
        "literal_weierstrass_identity": True,
    },
    "method": {
        "integral_marking_composition": True,
        "exact_group_law": True,
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (BRIDGE, HEIGHTS, SECTIONS, MODEL, PINNED)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (BRIDGE, HEIGHTS, SECTIONS, MODEL, PINNED)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDFINALA1HORIZONTAL|word={}|degrees=({},{},{})|PdotO=6|bits={}|"
    "seconds={:.3f}|status={}|output={}".format(
        ",".join(map(str, word)), xn.degree(), yn.degree(), Z.degree(),
        payload["section"]["maximum_rational_bits"],
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
