#!/usr/bin/env sage -python
"""Compress the exact q4/orbit208 3A3 Weierstrass equation over QQ.

Send the three finite I4 supports to 1, 0, and infinity, with the second
old-I6 I4 (the q4/orbit323 middle-component seed fibre) at t=0.  Apply the
corresponding weighted Weierstrass change and normalize its repeated root to
x=3.  This exact PGL2 base change reduces million-bit source coefficients to
a few hundred bits.  No elimination or Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
OUTPUT = LOCAL / "q4o208-compact-weierstrass-qq.json"

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def poly_bits(polynomial):
    return max(bits(value) for value in polynomial)


def polynomial_record(polynomial):
    return [str(value) for value in polynomial.list()]


model = json.loads(MODEL.read_text())
marking = json.loads(MARKING.read_text())
assert model["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"

RU = PolynomialRing(QQ, "U")
U = RU.gen()
RT = PolynomialRing(QQ, "t")
t = RT.gen()
A_old = RU([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B_old = RU([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
fibre_items = list(marking["physical_fibres"].items())
supports_old = [QQ(record["support"]) for unused, record in fibre_items]
nodes_old = [-3 * B_old(support) / (2 * A_old(support)) for support in supports_old]
assert all(
    A_old(support) == -3 * node**2 and B_old(support) == 2 * node**3
    for support, node in zip(supports_old, nodes_old)
)

# The permutation (1,0,2) puts the p=59-good selected fibre at zero, the
# first physical I4 at one, and the special I4 at infinity.
i, j, k = 1, 0, 2
s0, s1, s2 = [supports_old[index] for index in (i, j, k)]
a = (s1 - s0) * s2
b = -(s1 - s2) * s0
c = s1 - s0
d = -(s1 - s2)
N = a * t + b
D = c * t + d
assert a * d - b * c

# U=N/D.  Multiplication by D^8 and D^12 is the O(4),O(6) base change.
# At t=0 the repeated root is d^4*nodes_old[i]; normalize it to 3.
assert (nodes_old[i] / 3).is_square()
m = d**2 * (nodes_old[i] / 3).sqrt()
A = RT(sum(A_old[degree] * N**degree * D**(8-degree) for degree in range(9)) / m**4)
B = RT(sum(B_old[degree] * N**degree * D**(12-degree) for degree in range(13)) / m**6)
assert A.degree() == 8 and B.degree() == 12
assert A[0] == -27 and B[0] == 54

Delta_old = -16 * (4 * A_old**3 + 27 * B_old**2)
Delta = -16 * (4 * A**3 + 27 * B**2)
Delta_transformed = RT(
    sum(Delta_old[degree] * N**degree * D**(24-degree) for degree in range(25))
    / m**12
)
assert Delta == Delta_transformed
assert Delta.degree() == 20
assert Delta.valuation() == 4
assert Delta(t + 1).valuation() == 4
finite_nodal = Delta // (t**4 * (t - 1)**4)
assert finite_nodal.degree() == 12
assert finite_nodal.gcd(finite_nodal.derivative()).degree() == 0
assert finite_nodal(0) and finite_nodal(1)


def load_rational(record):
    numerator = RU([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RU([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return numerator, denominator


def transform_coordinate(record, weight):
    numerator, denominator = load_rational(record)
    numerator_homogeneous = sum(
        numerator[degree] * N**degree * D**(weight-degree)
        for degree in range(len(numerator.list()))
    )
    denominator_substituted = sum(
        denominator[degree] * N**degree * D**(denominator.degree()-degree)
        for degree in range(len(denominator.list()))
    )
    value = RT(numerator_homogeneous) / (m**(weight//2) * RT(denominator_substituted))
    return value


section_old = marking["old_A11_component_7_on_C5_pointed_child"]
x_section = transform_coordinate(section_old["x"], 4)
y_section = transform_coordinate(section_old["y"], 6)
assert y_section**2 == x_section**3 + A * x_section + B

payload = {
    "schema": "elkies-k3.h3-q4o208-compact-weierstrass-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION",
    "compact_model": {
        "equation": "y^2 = x^3 + A(t)*x + B(t)",
        "A_coefficients_low_to_high": polynomial_record(A),
        "B_coefficients_low_to_high": polynomial_record(B),
        "degrees_A_B_Delta": [8, 12, 20],
        "maximum_A_rational_bits": poly_bits(A),
        "maximum_B_rational_bits": poly_bits(B),
        "reducible_fibres": [
            {"label": fibre_items[i][0], "kodaira": "I4", "support": "0"},
            {"label": fibre_items[j][0], "kodaira": "I4", "support": "1"},
            {"label": fibre_items[k][0], "kodaira": "I4", "support": "infinity"},
        ],
        "finite_nodal_factor_degree": 12,
        "finite_nodal_factor_squarefree": True,
        "ADE": "3A3",
        "MW_rank_if_rho19": 8,
    },
    "exact_coordinate_change": {
        "base": "U=(a*t+b)/(c*t+d)",
        "x": "x_old=m^2*x/(c*t+d)^4",
        "y": "y_old=m^3*y/(c*t+d)^6",
        "a": str(a), "b": str(b), "c": str(c), "d": str(d), "m": str(m),
        "old_support_indices_to_new_1_0_infinity": [j, i, k],
        "selected_old_support_index": i,
        "selected_new_support": "0",
        "selected_new_repeated_root": "3",
    },
    "transported_exact_section": {
        "source_label": "old_A11_component_7_on_C5_pointed_child",
        "x_coefficients_low_to_high": polynomial_record(RT(x_section)),
        "y_coefficients_low_to_high": polynomial_record(RT(y_section)),
        "exact_compact_model_identity": True,
    },
    "coefficient_growth": {
        "old_maximum_A_B_rational_bits": max(poly_bits(A_old), poly_bits(B_old)),
        "new_maximum_A_rational_bits": poly_bits(A),
        "new_maximum_B_rational_bits": poly_bits(B),
    },
    "method": {
        "large_Groebner_required": False,
        "polynomial_factorization_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This is an exact QQ PGL2/weighted-Weierstrass isomorphism of the certified "
        "q4/orbit208 model. It verifies the discriminant transform, the two finite I4 "
        "orders, the I4 fibre at infinity, residual squarefreeness, and one transported "
        "exact section. The q4/orbit323 horizontal remains a separate equation gate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, MARKING)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, MARKING)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O208COMPACT|A_bits={}|B_bits={}|old_bits={}|status={}|output={}".format(
        poly_bits(A), poly_bits(B), payload["coefficient_growth"]["old_maximum_A_B_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
