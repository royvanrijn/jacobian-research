#!/usr/bin/env sage -python
"""Point the physical q4/orbit164 quartic at old_A11_component_8 over QQ.

C8 is the identity component of the unique finite I4 fibre on the q4/o1584
parent.  Specializing the stored q4/o164 quartic at that support gives an
exact polynomial square.  Of its two signs, one reconstructs a variable
nonsingular point on the nodal cubic and the other reconstructs the node;
the variable branch is therefore exactly C8.  Translating the quartic at
that branch supplies the pointed q4/o164 Jacobian.  No Groebner basis is
used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, factorial


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PARENT = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
PARENT_MARKING = LOCAL / "q4o1584-second-affine-equation-marking-qq.json"
Q164 = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
C8_MARKING = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
OUTPUT = LOCAL / "q4o164-c8-equation-marking-qq.json"
INPUTS = (PARENT, PARENT_MARKING, Q164, C8_MARKING)

started = time.monotonic()


def log(stage, **fields):
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O164C8MARKQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


def rational_bits(values):
    entries = [QQ(value) for value in values]
    return max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in entries
    )


parent = json.loads(PARENT.read_text())
parent_marking = json.loads(PARENT_MARKING.read_text())
q164 = json.loads(Q164.read_text())
c8_marking = json.loads(C8_MARKING.read_text())
assert parent["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"
assert parent_marking["status"] == "PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING"
assert q164["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert c8_marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert c8_marking["zero"] == "old_A11_component_8"

RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
RT = PolynomialRing(KU, "T")
T = RT.gen()


def read_polynomial(record):
    numerator = RU([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RU([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    quotient = numerator / denominator
    if quotient.denominator() != 1:
        raise ArithmeticError("stored q4/o1584 section is not polynomial")
    return RU(quotient)


def record(value):
    value = KU(value)
    numerator = RU(value.numerator())
    denominator = RU(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    if denominator.leading_coefficient() < 0:
        numerator, denominator = -numerator, -denominator
    return {
        "numerator_coefficients_low_to_high": [str(entry) for entry in numerator.list()],
        "denominator_coefficients_low_to_high": [str(entry) for entry in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


A_parent = RU([QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]])
B_parent = RU([QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]])
c0_point = parent_marking["old_A11_component_0_on_second_affine_pointed_child"]
PX = read_polynomial(c0_point["x"])
PY = read_polynomial(c0_point["y"])
if PY**2 != PX**3 + A_parent * PX + B_parent:
    raise ArithmeticError("C0 misses the q4/o1584 parent model")

i4_records = [
    item for item in parent["child"]["finite_reducible_fibres"]
    if item["kodaira"] == "I4"
]
if len(i4_records) != 1:
    raise ArithmeticError("q4/o1584 parent does not have a unique finite I4")
t0 = QQ(i4_records[0]["support"])

quartic_coefficients = [
    RU([QQ(value) for value in coefficients])
    for coefficients in q164["quartic"]["coefficients_in_T_low_to_high"]
]
quartic = RT([KU(value) for value in quartic_coefficients])
special = RU(quartic(KU(t0)))


def polynomial_square_root(value):
    value = RU(value)
    shift = next(QQ(integer) for integer in range(8) if value(QQ(integer)) != 0)
    shifted = RU(value(U + shift))
    constant = QQ(shifted[0])
    if not constant.is_square():
        raise ArithmeticError("specialized quartic has nonsquare constant after shift")
    coefficients = [QQ(constant.sqrt())]
    root_degree = value.degree() // 2
    for degree in range(1, root_degree + 1):
        known = sum(
            (coefficients[left] * coefficients[degree-left]
             for left in range(1, degree)),
            QQ.zero(),
        )
        coefficients.append(QQ((shifted[degree] - known) / (2 * coefficients[0])))
    shifted_root = RU(coefficients)
    root = RU(shifted_root(U - shift))
    if root**2 != value:
        raise ArithmeticError("coefficient-recursion square root failed")
    return root


W_root = polynomial_square_root(special)
linear = RU(q164["resolved_RR"]["interpolating_linear_polynomial"])
denominator = RU(q164["resolved_RR"]["common_denominator"])
m_special = RU(linear(t0) + U * denominator(t0))
node = QQ(-3 * B_parent(t0) / (2 * A_parent(t0)))
if node**3 + A_parent(t0) * node + B_parent(t0):
    raise ArithmeticError("finite I4 support has the wrong cubic node")

recovered = []
for ordinate in (W_root, -W_root):
    x_value = RU((denominator(t0) * ordinate - PX(t0) + m_special**2) / 2)
    y_value = RU(m_special * (x_value - PX(t0)) - PY(t0))
    if y_value**2 != x_value**3 + A_parent(t0) * x_value + B_parent(t0):
        raise ArithmeticError("specialized branch misses the q4/o1584 nodal cubic")
    recovered.append((ordinate, x_value, y_value))

variable_indices = [
    index for index, (unused, x_value, unused_y) in enumerate(recovered)
    if x_value.degree() > 0
]
node_indices = [
    index for index, (unused, x_value, y_value) in enumerate(recovered)
    if x_value == node and y_value == 0
]
if len(variable_indices) != 1 or len(node_indices) != 1 or variable_indices == node_indices:
    raise ArithmeticError("C8 identity branch is not uniquely separated from the I4 node")
c8_index = variable_indices[0]
W_c8, x_c8_parent, y_c8_parent = recovered[c8_index]
log(
    "C8_BRANCH", selected_index=c8_index,
    parent_x_degree=x_c8_parent.degree(), node_index=node_indices[0],
)

# Point the q4/o164 quartic at (T,W)=(t0,W_c8(U)).
translated = [
    KU(quartic.derivative(order)(KU(t0)) / factorial(order))
    for order in range(5)
]
e, d, c, b, a = translated
W0 = KU(W_c8)
if e != W0**2:
    raise ArithmeticError("C8 is not an exact point on the q4/o164 quartic")
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
a4 = -4 * W0**2 * a
a6 = a2 * a4
b2 = a1**2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3**2 + 4 * a6
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
pointed_A = -c4 / 48
pointed_B = -c6 / 864
A_short = KU(RU([QQ(value) for value in q164["child"]["minimal_A_coefficients_low_to_high"]]))
B_short = KU(RU([QQ(value) for value in q164["child"]["minimal_B_coefficients_low_to_high"]]))
if 81 * pointed_A != A_short or 729 * pointed_B != B_short:
    raise ArithmeticError("C8-pointed quartic misses the stored q4/o164 short Jacobian")
E = EllipticCurve(KU, [0, 0, 0, A_short, B_short])

# The opposite sign above the same old-base coordinate has the exact limit
# x=-a2, y=a1*a2-a3 on the pointed generalized model.
xg = -a2
yg = a1 * a2 - a3
x_opposite = KU(9 * (xg + b2 / 12))
y_opposite = KU(27 * (yg + (a1 * xg + a3) / 2))
opposite_point = E(x_opposite, y_opposite)
if opposite_point.is_zero():
    raise ArithmeticError("opposite I4-node branch collapsed to the selected C8 zero")
log("POINTED_JACOBIAN", A=True, B=True, opposite=True)

payload = {
    "schema": "elkies-k3.q4o164-c8-equation-marking-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING",
    "selected_zero": {
        "label": "old_A11_component_8",
        "parent_fibre": "unique finite I4 identity component",
        "parent_base_support": str(t0),
        "selected_sign_index": c8_index,
        "quartic_ordinate": record(W0),
        "recovered_parent_x": record(KU(x_c8_parent)),
        "recovered_parent_y": record(KU(y_c8_parent)),
        "parent_point_is_variable_nonsingular": True,
        "opposite_sign_recovers_parent_node": True,
        "maps_to": "point at infinity on the C8-pointed generalized Weierstrass model",
    },
    "pointed_jacobian": {
        "standard_short_scaling": "x_short=9*x_pointed_short, y_short=27*y_pointed_short",
        "exact_A_identity": True,
        "exact_B_identity": True,
    },
    "opposite_constant_support_section": {
        "x": record(x_opposite),
        "y": record(y_opposite),
        "exact_child_identity": True,
    },
    "coefficient_growth": {
        "C8_quartic_ordinate_maximum_rational_bits": rational_bits(W_root.list()),
        "opposite_section_maximum_rational_bits": rational_bits(
            list(x_opposite.numerator()) + list(x_opposite.denominator())
            + list(y_opposite.numerator()) + list(y_opposite.denominator())
        ),
    },
    "method": {
        "large_Groebner_required": False,
        "polynomial_factorization_required_for_square_root": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "next_required": (
        "Recover the promoted q8/o376 P.O=4 horizontal section, then compile its "
        "expected 12-dimensional RR ambient space."
    ),
    "proof_boundary": (
        "The unique finite-I4 identity branch is exactly C8, the opposite sign "
        "recovers the contracted node, and pointing at C8 gives exactly the stored "
        "q4/o164 short Jacobian. This attaches the q8/o376 source zero; its P.O=4 "
        "horizontal section and q8 equation remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("DONE", zero="C8", status=payload["status"], output=OUTPUT)
