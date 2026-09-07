#!/usr/bin/env sage -python
"""Apply the compact q4/orbit164 scaling to its exact binary quartic.

If U=c*t and the child Weierstrass coordinates scale by x_old=s^2*x,
y_old=s^3*y, then dividing the quartic by s^2 gives invariants scaled by
s^-4 and s^-6.  The resulting quartic has exactly the compact child as its
invariant Jacobian.  No elimination or factorization is used.
"""

import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
RAW = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
COMPACT = LOCAL / "q4o164-compact-weierstrass-qq.json"
PARENT_COMPACT = LOCAL / "q4o1584-compact-weierstrass-qq.json"
OUTPUT = LOCAL / "q4o164-compact-binary-quartic-qq.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


raw = json.loads(RAW.read_text())
compact = json.loads(COMPACT.read_text())
parent_compact = json.loads(PARENT_COMPACT.read_text())
assert raw["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert compact["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert parent_compact["status"] == "PASS_EXACT_QQ_Q4O1584_COMPACT_WEIERSTRASS_NORMALIZATION"

RU = PolynomialRing(QQ, "U")
U = RU.gen()
RT = PolynomialRing(QQ, "t")
t = RT.gen()
c_scale = QQ(compact["exact_coordinate_change"]["c"])
s_scale = QQ(compact["exact_coordinate_change"]["s"])
raw_coefficients = [
    RU([QQ(value) for value in row])
    for row in raw["quartic"]["coefficients_in_T_low_to_high"]
]
parent_a = QQ(parent_compact["exact_coordinate_change"]["a"])
parent_c = QQ(parent_compact["exact_coordinate_change"]["c"])
base_scaled = [RT(value(c_scale * t)) for value in raw_coefficients]
# Substitute T_old=parent_a+parent_c*V.  A binary quartic's I,J invariants
# acquire determinant powers parent_c^4,parent_c^6, so division by
# parent_c^2 (as well as the child ordinate scale s) preserves the compact
# invariant Jacobian.
coefficients = [RT.zero()] * 5
for old_degree, polynomial in enumerate(base_scaled):
    for new_degree in range(old_degree + 1):
        coefficients[new_degree] += (
            polynomial
            * ZZ(old_degree).binomial(new_degree)
            * parent_a**(old_degree-new_degree)
            * parent_c**new_degree
            / (s_scale**2 * parent_c**2)
        )
e, d, cc, b, a = coefficients
I = 12 * a * e - 3 * b * d + cc**2
J = 72 * a * cc * e + 9 * b * cc * d - 27 * a * d**2 - 27 * b**2 * e - 2 * cc**3
A = RT([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = RT([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
assert -27 * I == A and -27 * J == B

payload = {
    "schema": "elkies-k3.q4o164-compact-binary-quartic-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_COMPACT_BINARY_QUARTIC",
    "equation": "W^2=sum(q_i(t)*V^i,i=0..4)",
    "coefficients_in_T_low_to_high": [
        [str(value) for value in polynomial.list()] for polynomial in coefficients
    ],
    "coefficient_degrees_in_t": [int(value.degree()) for value in coefficients],
    "maximum_rational_bits": max(bits(value) for polynomial in coefficients for value in polynomial),
    "exact_coordinate_change": {
        "old_base": "U=c*t",
        "old_quartic_variable": "T_old=parent_a+parent_c*V",
        "quartic_ordinate": "W_old=s*parent_c*W",
        "c": str(c_scale),
        "s": str(s_scale),
        "parent_a": str(parent_a),
        "parent_c": str(parent_c),
    },
    "invariant_jacobian": {
        "A_equals_compact_child_A": True,
        "B_equals_compact_child_B": True,
    },
    "method": {"large_Groebner_required": False, "factorization_required": False},
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (RAW, COMPACT, PARENT_COMPACT)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (RAW, COMPACT, PARENT_COMPACT)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164QUARTICCOMPACT|degrees={}|bits={}|status={}|output={}".format(
        payload["coefficient_degrees_in_t"], payload["maximum_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
