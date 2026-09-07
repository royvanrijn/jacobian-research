#!/usr/bin/env sage -python
"""Compress the exact q4/orbit1584 Weierstrass equation over QQ.

Translate its finite I4 support to zero, put the preceding finite I2 support
at one, and scale x,y so the I4 node is x=3.  The huge displayed source
coefficients are almost entirely coordinate growth: this exact change lowers
them to a few hundred bits and prepares small resolved-Hensel section lifts.
No elimination or Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
OUTPUT = LOCAL / "q4o1584-compact-weierstrass-qq.json"

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def poly_bits(polynomial):
    return max(bits(value) for value in polynomial)


model = json.loads(MODEL.read_text())
assert model["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"

RU = PolynomialRing(QQ, "U")
U = RU.gen()
RT = PolynomialRing(QQ, "t")
t = RT.gen()
A_old = RU([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B_old = RU([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
records = model["child"]["finite_reducible_fibres"]
assert [record["kodaira"] for record in records] == ["I2", "I2", "I2", "I4"]
supports_old = [QQ(record["support"]) for record in records]

a = supports_old[3]
c = supports_old[2] - a
node_old = -3 * B_old(a) / (2 * A_old(a))
assert A_old(a) == -3 * node_old**2 and B_old(a) == 2 * node_old**3
assert (node_old / 3).is_square()
m = (node_old / 3).sqrt()

A = RT(A_old(a + c * t)) / m**4
B = RT(B_old(a + c * t)) / m**6
assert A.degree() == 6 and B.degree() == 9
supports = [(value - a) / c for value in supports_old]
assert supports[2:] == [1, 0]

Delta_old = -16 * (4 * A_old**3 + 27 * B_old**2)
Delta = -16 * (4 * A**3 + 27 * B**2)
assert Delta == RT(Delta_old(a + c * t)) / m**12
assert Delta.degree() == 18
orders = [2, 2, 2, 4]
for support, order in zip(supports, orders):
    assert Delta(t + support).valuation() == order
nodal = Delta
for support, order in zip(supports, orders):
    nodal //= (t - support)**order
assert nodal.degree() == 8
assert nodal.gcd(nodal.derivative()).degree() == 0
assert all(nodal(support) for support in supports)

payload = {
    "schema": "elkies-k3.q4o1584-compact-weierstrass-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O1584_COMPACT_WEIERSTRASS_NORMALIZATION",
    "compact_model": {
        "equation": "y^2 = x^3 + A(t)*x + B(t)",
        "A_coefficients_low_to_high": [str(value) for value in A.list()],
        "B_coefficients_low_to_high": [str(value) for value in B.list()],
        "degrees_A_B_Delta": [6, 9, 18],
        "maximum_A_rational_bits": poly_bits(A),
        "maximum_B_rational_bits": poly_bits(B),
        "finite_reducible_fibres": [
            {"kodaira": record["kodaira"], "delta_order": order, "support": str(support)}
            for record, order, support in zip(records, orders, supports)
        ],
        "infinity": {"kodaira": "I0*", "orders_A_B_Delta": [2, 3, 6]},
        "finite_nodal_factor_degree": 8,
        "finite_nodal_factor_squarefree": True,
        "ADE": "D4+A3+3A1",
        "MW_rank_if_rho19": 7,
    },
    "exact_coordinate_change": {
        "base": "U = a+c*t",
        "x": "x_old = m^2*x",
        "y": "y_old = m^3*y",
        "a": str(a),
        "c": str(c),
        "m": str(m),
        "old_I4_node": str(node_old),
        "new_I4_node": "3",
    },
    "coefficient_growth": {
        "old_maximum_A_B_rational_bits": int(model["child"]["maximum_A_B_rational_bits"]),
        "new_maximum_A_rational_bits": poly_bits(A),
        "new_maximum_B_rational_bits": poly_bits(B),
    },
    "method": {
        "large_Groebner_required": False,
        "polynomial_factorization_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This is an exact QQ Weierstrass isomorphism of the certified q4/orbit1584 model. "
        "The discriminant transform, all finite reducible orders, residual squarefreeness, "
        "and the I0* order profile at infinity are checked exactly."
    ),
    "next_required": (
        "Recover one of the two parent-degree-one sections that meets every reducible-fibre "
        "node, then transport it through the q4/orbit164 quartic map."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): sha256(MODEL)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O1584COMPACT|A_bits={}|B_bits={}|status={}|output={}".format(
        poly_bits(A), poly_bits(B), payload["status"], OUTPUT,
    ),
    flush=True,
)
