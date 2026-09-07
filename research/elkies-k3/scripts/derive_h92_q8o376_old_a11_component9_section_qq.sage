#!/usr/bin/env sage -python
"""Recover old_A11_component_9 on the P1229-pointed q8/o376 child.

The curve is the affine component F+e1 of the second compact I2 of the
q4/o164 model, at T=1.  Normalize that nodal cubic, restrict the stored
resolved q8 pencil, and use the restriction to choose the sign of the square
quartic value at T=1.  The stored generalized P1229-pointed quartic map then
gives the child section by direct substitution.  No elimination is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
CLASSIFIER = LOCAL / "q12o5867-p0-shell-lattice-classification-mod89.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q8o376-old-a11-component9-section-qq.json",
)
args = parser.parse_args()
OUTPUT = args.output.resolve()
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
model = json.loads(MODEL.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
marking = json.loads(MARKING.read_text())
classifier = json.loads(CLASSIFIER.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"

source_class = marking["equation_explicit_curves_in_child"]["old_A11_component_9"]
assert source_class == marking["equation_explicit_curves_in_child"]["old_zero"]
assert source_class == [1, 0, 0, 1, 0] + [0]*14
class460 = classifier["lattice_shell"]["classes"][460]
assert class460["class_index"] == 460
assert class460["q4o164_parent_curve"] == source_class
assert class460["q4o164_parent_degree"] == 0
assert class460["q4o164_parent_a_minus_b"] == 1

RT = PolynomialRing(QQ, "T")
T = RT.gen()
RU = PolynomialRing(QQ, "u")
u = RU.gen()
KU = RU.fraction_field()
TUP = PolynomialRing(RU, "T")


def polynomial(values):
    return RT([QQ(value) for value in values])


def rational(record):
    return polynomial(record["numerator_coefficients_low_to_high"])/polynomial(
        record["denominator_coefficients_low_to_high"]
    )


def u_rational(record):
    numerator = RU([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RU([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return KU(numerator)/KU(denominator)


def rational_record(value, ring):
    value = ring.fraction_field()(value)
    return {
        "numerator_coefficients_low_to_high": [str(x) for x in value.numerator().list()],
        "denominator_coefficients_low_to_high": [str(x) for x in value.denominator().list()],
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree())
        ],
    }


def rational_bits(values):
    values = [QQ(value) for value in values]
    return int(max(
        max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
        for value in values
    ))


compact = model["compact_model"]
A = polynomial(compact["A_coefficients_low_to_high"])
B = polynomial(compact["B_coefficients_low_to_high"])
assert compact["reducible_fibres"][1] == {"kodaira": "I2", "support": "1"}

Hx = rational(horizontal["section"]["x"])
Hy = rational(horizontal["section"]["y"])
X, denominator_x = Hx.numerator(), Hx.denominator()
Y, denominator_y = Hy.numerator(), Hy.denominator()
Z = RT.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic()**(int(exponent)//2)
assert Z**2 == denominator_x and Z**3 == denominator_y
assert Y**2 == X**3+A*X*Z**4+B*Z**6

resolved_pairs = []
for record in q8["resolved_RR"]["resolved_basis_pairs"]:
    resolved_pairs.append((
        polynomial(record["AA_coefficients_low_to_high"]),
        polynomial(record["BB_coefficients_low_to_high"]),
    ))
(AA0, BB0), (AA1, BB1) = resolved_pairs

quartic = TUP([
    RU([QQ(value) for value in coefficient])
    for coefficient in q8["quartic"]["coefficients_in_old_T_low_to_high"]
])

# Recover the stored quartic square factor by the same bounded factorization
# used in the q8 compiler.  This only fixes the sign convention for the
# resolved restriction; it is not an elimination calculation.
TK = PolynomialRing(KU, "old_T")


def lift_TK(poly):
    return TK([KU(value) for value in RT(poly).list()])


aa = lift_TK(AA0)+KU(u)*lift_TK(AA1)
bb = lift_TK(BB0)+KU(u)*lift_TK(BB1)
XX, YY, ZZpoly, AApoly = map(lift_TK, (X, Y, Z, A))
raw = (
    aa**4-6*XX*aa**2*bb**2+8*YY*aa*bb**3
    -3*XX**2*bb**4-4*AApoly*bb**4*ZZpoly**4
)
after_collision, remainder = raw.quo_rem(ZZpoly**4)
assert not remainder
square_quotient, remainder = TUP(after_collision).quo_rem(quartic)
assert not remainder
factorization = square_quotient.factor()
unit = QQ(factorization.unit())
assert unit.is_square()
square_factor = TUP(unit.sqrt())
for factor, exponent in factorization:
    assert int(exponent) % 2 == 0
    square_factor *= factor**(int(exponent)//2)
assert quartic*square_factor**2 == TUP(after_collision)


def evaluate_u(poly, value, field):
    answer = field.zero()
    for coefficient in reversed(RU(poly).list()):
        answer = answer*value+field(QQ(coefficient))
    return answer


def evaluate_bivariate(poly, old_t_value, u_value, field):
    answer = field.zero()
    for coefficient in reversed(TUP(poly).list()):
        answer = answer*old_t_value+evaluate_u(coefficient, u_value, field)
    return answer


# Normalize the affine component of the second I2 nodal cubic.  At T=1 the
# cubic is (x-node)^2*(x+2*node), with normalization
# x=r^2-2*node, y=r*(r^2-3*node).
Kr = FunctionField(QQ, "r")
r = Kr.gen()
A1, B1 = QQ(A(1)), QQ(B(1))
node = QQ(-3*B1/(2*A1))
affine_x = r**2-2*node
affine_y = r*(r**2-3*node)
assert affine_y**2 == affine_x**3+A1*affine_x+B1
horizontal_x_at_1 = QQ(X(1))/QQ(Z(1))**2
horizontal_y_at_1 = QQ(Y(1))/QQ(Z(1))**3
slope = (affine_y+horizontal_y_at_1)/(affine_x-horizontal_x_at_1)
restrictions = [
    Kr(AA(1))+Kr(BB(1))*Kr(Z(1))*slope for AA, BB in resolved_pairs
]
u_of_r = -restrictions[0]/restrictions[1]
assert max(u_of_r.numerator().degree(), u_of_r.denominator().degree()) == 1
aa_on_component = Kr(AA0(1))+u_of_r*Kr(AA1(1))
bb_on_component = Kr(BB0(1))+u_of_r*Kr(BB1(1))
square_on_component = evaluate_bivariate(square_factor, Kr(1), u_of_r, Kr)
W_on_component = (
    bb_on_component**2
    *(2*affine_x+horizontal_x_at_1-slope**2)
    / square_on_component
)
assert W_on_component**2 == evaluate_bivariate(quartic, Kr(1), u_of_r, Kr)

quartic_at_second_I2 = RU(sum(quartic.list()))
assert KU(quartic_at_second_I2).is_square()
selected_W = KU(quartic_at_second_I2).sqrt()
assert W_on_component == evaluate_u(selected_W, u_of_r, Kr)

# Apply the generalized P1229-pointed map at the constant old coordinate T=1.
p1229_support = QQ(q8["preferred_pointed_zero"]["old_base_coordinate"])
q_origin = u_rational(q8["preferred_pointed_zero"]["quartic_ordinate"])
a1, a2, a3, unused_a4, unused_a6 = [
    u_rational(record)
    for record in q8["preferred_pointed_zero"]["generalized_weierstrass_a_invariants"]
]
d_coefficient = a1*q_origin
c_coefficient = a2+d_coefficient**2/(4*q_origin**2)
b2 = a1**2+4*a2
local_T = KU(QQ(1)-p1229_support)
x_generalized = (
    2*q_origin*(selected_W+q_origin)+d_coefficient*local_T
)/local_T**2
y_generalized = (
    4*q_origin**2*(selected_W+q_origin)
    +2*q_origin*(d_coefficient*local_T+c_coefficient*local_T**2)
    -d_coefficient**2*local_T**2/(2*q_origin)
)/local_T**3
child_x = KU(9*(x_generalized+b2/12))
child_y = KU(27*(y_generalized+(a1*x_generalized+a3)/2))

A_child = RU([QQ(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]])
B_child = RU([QQ(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]])
assert child_y**2 == child_x**3+KU(A_child)*child_x+KU(B_child)
assert child_x.denominator() == 1 and child_y.denominator() == 1
assert (child_x.numerator().degree(), child_y.numerator().degree()) == (4, 6)

# Literal inverse-map regression pins the constant parent coordinate and the
# chosen second-I2 quartic sign independently of the child equation check.
inverse_x_generalized = child_x/9-b2/12
inverse_y_generalized = child_y/27-(a1*inverse_x_generalized+a3)/2
inverse_local_T = 2*q_origin*(inverse_x_generalized+a2)/inverse_y_generalized
inverse_old_T = p1229_support+inverse_local_T
inverse_W = (
    inverse_x_generalized*inverse_local_T**2-d_coefficient*inverse_local_T
)/(2*q_origin)-q_origin
assert inverse_old_T == 1
assert inverse_W == selected_W

prime = ZZ(137)
F = GF(prime)


def reduce_qq(value):
    value = QQ(value)
    assert value.denominator() % prime
    return F(value.numerator())/F(value.denominator())


mod137_x = [int(reduce_qq(value)) for value in child_x.numerator().list()]
mod137_y = [int(reduce_qq(value)) for value in child_y.numerator().list()]
FU = PolynomialRing(F, "u")
mod_A = FU([reduce_qq(value) for value in A_child.list()])
mod_B = FU([reduce_qq(value) for value in B_child.list()])
mod_x = FU(mod137_x)
mod_y = FU(mod137_y)
assert mod_y**2 == mod_x**3+mod_A*mod_x+mod_B

inputs = (Q8, MODEL, HORIZONTAL, MARKING, CLASSIFIER)
payload = {
    "schema": "elkies-k3.h92-q8o376-old-a11-component9-section-qq.v1",
    "status": "PASS_EXACT_QQ_Q8O376_OLD_A11_COMPONENT9_SECTION",
    "curve_identification": {
        "label": "old_A11_component_9",
        "lattice_class_index": 460,
        "q4o164_parent_curve": source_class,
        "q4o164_parent_degree": 0,
        "q4o164_parent_a_minus_b": 1,
        "q4_fibre_component": "F+e1, affine component of the second compact I2",
        "old_base_support": "1",
    },
    "second_I2_resolved_restriction": {
        "nodal_cubic_node": str(node),
        "normalization": {
            "x_of_r": str(affine_x),
            "y_of_r": str(affine_y),
        },
        "new_q8_base_u_of_r": str(u_of_r),
        "new_q8_base_degree": 1,
        "quartic_W_of_r": str(W_on_component),
        "selected_quartic_sign": 1,
        "sign_identity": "W_on_affine_component = +sqrt(quartic(T=1,u)) after the stored square-factor convention",
        "exact_resolved_pencil_restriction": True,
        "exact_quartic_identity": True,
    },
    "pointed_quartic": {
        "origin_label": "P1229",
        "origin_old_base_support": str(p1229_support),
        "source_old_base_coordinate": "1",
        "source_quartic_ordinate": rational_record(selected_W, RU),
        "map": (
            "L=T-T0; xg=(2*q0*(W+q0)+d*L)/L^2; "
            "yg=(4*q0^2*(W+q0)+2*q0*(d*L+c*L^2)-d^2*L^2/(2*q0))/L^3; "
            "x=9*(xg+b2/12), y=27*(yg+(a1*xg+a3)/2)"
        ),
        "literal_inverse_old_base_coordinate": "1",
        "literal_inverse_quartic_ordinate": True,
    },
    "section": {
        "x": rational_record(child_x, RU),
        "y": rational_record(child_y, RU),
        "degrees_x_y": [4, 6],
        "P_dot_O": 0,
        "exact_child_weierstrass_identity": True,
        "maximum_rational_bits": rational_bits(
            child_x.numerator().list()+child_y.numerator().list()
        ),
    },
    "mod137_regression": {
        "x_coefficients_low_to_high": mod137_x,
        "y_coefficients_low_to_high": mod137_y,
        "exact_child_weierstrass_identity": True,
    },
    "method": {
        "construction": "second-I2 nodal normalization, resolved-pencil restriction, generalized P1229-pointed quartic map",
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
    "proof_boundary": (
        "The exact class-460 curve F+e1 is restricted through the stored q8/o376 resolved pencil, "
        "its quartic sign is fixed on the normalized affine component of the second q4 I2, and the "
        "resulting P1229-pointed child section satisfies the exact QQ and mod-137 child equations. "
        "This is a coordinate/marking certificate, not a q12 resolved-RR construction."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q8O376C9QQ|degrees=4,6|P.O=0|sign=+|mod137=True|runtime={:.3f}|status={}|output={}".format(
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
