#!/usr/bin/env sage -python
"""Construct the degree-one parent branch of the q12/o5867 compiler mod 131.

The promoted four-section word contains one curve of degree one over the
q4/o164 base.  In the selected marked embedding its Mordell--Weil class is an
integral word in the already exact B0,...,B7 sections and the saturated
C8-opposite point.  Compute that word by exact lattice linear algebra, evaluate
it by elliptic group law modulo 131, restrict the resolved q8/o376 pencil to
the resulting section, and apply the pointed-quartic map based at P1229.

This gives a literal polynomial P.O=0 section on the q8 child without a
Groebner basis.  It is a good-reduction construction seed, not a QQ lift.
"""

import hashlib
import json
import time
from math import comb
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
HEIGHT = LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
FRONTIER = GEN / "elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json"
OUTPUT = LOCAL / "q12o5867-degree1-compiler-branch-mod131.json"
INPUTS = (MODEL, BASIS, HEIGHT, C8, HORIZONTAL, Q8, FRONTIER)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
prime = ZZ(131)
F = GF(prime)
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()
U = PolynomialRing(F, "u")
u = U.gen()
KU = U.fraction_field()
TU = PolynomialRing(KU, "T")
T_u = TU.gen()
TUP = PolynomialRing(U, "T")

model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
height = json.loads(HEIGHT.read_text())
c8 = json.loads(C8.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
q8 = json.loads(Q8.read_text())
frontier = json.loads(FRONTIER.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert height["status"] == "PASS_EXACT_QQ_Q4O164_FOURFOLD_HEIGHT_GRAM_AND_C8_MARKED_EMBEDDING_CENSUS"
assert c8["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert frontier["status"] == "PASS_EXACT_ROOTLESS_P0_SECTION_WORD_FRONTIER"


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


def polynomial(values):
    return R([reduce_qq(value) for value in values])


def rational_coordinate(record):
    return K(polynomial(record["numerator_coefficients_low_to_high"])) / K(
        polynomial(record["denominator_coefficients_low_to_high"])
    )


def substitute_scaled(function, scalar):
    return K(R(function.numerator()(scalar*T))) / K(R(function.denominator()(scalar*T)))


compact = model["compact_model"]
A = polynomial(compact["A_coefficients_low_to_high"])
B = polynomial(compact["B_coefficients_low_to_high"])


def checked_point(x_coordinate, y_coordinate):
    point = (K(x_coordinate), K(y_coordinate))
    assert point[1]**2 == point[0]**3+K(A)*point[0]+K(B)
    return point


def point_neg(point):
    return None if point is None else (point[0], -point[1])


def point_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+K(A))/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return checked_point(x3, slope*(x1-x3)-y1)


def point_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return point_mul(-coefficient, point_neg(point))
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = point_add(answer, addend)
        addend = point_add(addend, addend)
        coefficient >>= 1
    return answer


# Select q12/o5867 and its unique degree-one parent branch.
o5867_fibre = [6, 2, -16, 18, 2, 5, -11, 2, 13, 11, -6, -41, 6, 5, -14, 0, -3, -1, -2]
target = next(item for item in frontier["targets"] if item.get("fibre") == o5867_fibre)
compiler = target["best_four_P0_word_by_parent_a_minus_b"]
degree_one = [item for item in compiler["new_sections"] if item["q4o164_parent_degree"] == 1]
assert len(degree_one) == 1
branch = degree_one[0]
assert branch["current_4A1_P_dot_O"] == 0
assert branch["current_4A1_component_pairings"] == [0, 0, 0, 0]
parent_curve = vector(ZZ, branch["q4o164_parent_curve"])
target_tail = vector(QQ, parent_curve[-9:])

# Solve in the selected marked equation basis.  The B-lattice has index three;
# the exact 3*C8opp relation clears the fractional word.
marked = height["marked_embedding_enumeration"]
embedding = next(item for item in marked["embeddings"] if item["embedding_index"] == 15)
B_rows = matrix(QQ, embedding["rows_B0_through_B7_in_marked_MW9"])
H_tail = vector(QQ, marked["q8_horizontal_marked_MW9_tail"])
equation_rows = B_rows.stack(matrix(QQ, [H_tail]))
rational_word = equation_rows.transpose().solve_right(target_tail)
assert rational_word[-1] == 0
c8_relation = vector(ZZ, [-2, -3, -4, 3, -2, 2, -1, -2])
c8_basis_coordinates = vector(QQ, [QQ(value)/3 for value in c8_relation])
c8_coefficient = ZZ(2)
integral_B_word = vector(QQ, rational_word[:8])-c8_coefficient*c8_basis_coordinates
assert all(value in ZZ for value in integral_B_word)
integral_B_word = vector(ZZ, integral_B_word)
assert list(integral_B_word) == [0, 1, 2, -2, 1, -2, 1, 1]
assert target_tail == integral_B_word*B_rows+c8_coefficient*c8_basis_coordinates*B_rows

basis_points = [
    checked_point(
        polynomial(item["x_coefficients_low_to_high"]),
        polynomial(item["y_coefficients_low_to_high"]),
    )
    for item in basis["resolved_hensel"]["sections"]
]
c8_record = c8["opposite_constant_support_section"]
c8_old = (rational_coordinate(c8_record["x"]), rational_coordinate(c8_record["y"]))
base_scale = reduce_qq(model["exact_coordinate_change"]["c"])
xy_scale = reduce_qq(model["exact_coordinate_change"]["s"])
c8_opposite = checked_point(
    substitute_scaled(c8_old[0], base_scale)/xy_scale**2,
    substitute_scaled(c8_old[1], base_scale)/xy_scale**3,
)
relation_point = None
for coefficient, point in zip(c8_relation, basis_points):
    relation_point = point_add(relation_point, point_mul(coefficient, point))
assert point_mul(3, c8_opposite) == relation_point
parent_point = point_mul(c8_coefficient, c8_opposite)
for coefficient, point in zip(integral_B_word, basis_points):
    parent_point = point_add(parent_point, point_mul(coefficient, point))
parent_x, parent_y = parent_point
assert (parent_x.numerator().degree(), parent_x.denominator().degree()) == (6, 2)
assert (parent_y.numerator().degree(), parent_y.denominator().degree()) == (8, 3)

# Recover the q8 horizontal projective triple X,Y,Z.
X = polynomial(horizontal["section"]["x"]["numerator_coefficients_low_to_high"])
denominator_x = polynomial(horizontal["section"]["x"]["denominator_coefficients_low_to_high"])
Y = polynomial(horizontal["section"]["y"]["numerator_coefficients_low_to_high"])
denominator_y = polynomial(horizontal["section"]["y"]["denominator_coefficients_low_to_high"])
Z = R.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic()**(int(exponent)//2)
assert Z**2 == denominator_x and Z**3 == denominator_y

resolved_pairs = []
for item in q8["resolved_RR"]["resolved_basis_pairs"]:
    resolved_pairs.append((
        polynomial(item["AA_coefficients_low_to_high"]),
        polynomial(item["BB_coefficients_low_to_high"]),
    ))
(AA0, BB0), (AA1, BB1) = resolved_pairs

# Reconstruct only the small univariate square factor of the chord radicand.
def lift_TU(poly):
    return TU([KU(value) for value in R(poly).list()])


aa = lift_TU(AA0)+KU(u)*lift_TU(AA1)
bb = lift_TU(BB0)+KU(u)*lift_TU(BB1)
X_u, Y_u, Z_u, A_u = map(lift_TU, (X, Y, Z, A))
raw = (
    aa**4-6*X_u*aa**2*bb**2+8*Y_u*aa*bb**3
    -3*X_u**2*bb**4-4*A_u*bb**4*Z_u**4
)
after_collision, remainder = raw.quo_rem(Z_u**4)
assert not remainder
# Pin the exact QQ normalization reduced modulo 131, because the stored P1229
# ordinate uses that square-root sign and scalar.  Only factor the square
# quotient, whose factors have even multiplicity.
quartic = TUP([
    U([reduce_qq(value) for value in coefficient])
    for coefficient in q8["quartic"]["coefficients_in_old_T_low_to_high"]
])
square_quotient, square_remainder = TUP(after_collision).quo_rem(quartic)
assert not square_remainder
square_factorization = square_quotient.factor()
assert all(int(exponent) % 2 == 0 for unused, exponent in square_factorization)
square_unit = F(square_factorization.unit())
assert square_unit.is_square()
square_factor = TUP(square_unit.sqrt())
for factor, exponent in square_factorization:
    square_factor *= factor**(int(exponent)//2)
assert TUP(after_collision) == quartic*square_factor**2


def evaluate_bivariate(poly, old_t_value, u_value):
    answer = K.zero()
    for coefficient in reversed(TUP(poly).list()):
        coefficient_value = K.zero()
        for scalar in reversed(U(coefficient).list()):
            coefficient_value = coefficient_value*u_value+K(scalar)
        answer = answer*old_t_value+coefficient_value
    return answer


horizontal_x = K(X)/K(Z**2)
horizontal_y = K(Y)/K(Z**3)
slope = (parent_y+horizontal_y)/(parent_x-horizontal_x)
restrictions = [K(AA)+K(BB)*K(Z)*slope for AA, BB in resolved_pairs]
new_base = -restrictions[0]/restrictions[1]
assert max(new_base.numerator().degree(), new_base.denominator().degree()) == 1
aa_on_curve = K(AA0)+new_base*K(AA1)
bb_on_curve = K(BB0)+new_base*K(BB1)
square_on_curve = evaluate_bivariate(square_factor, K(T), new_base)
W_on_curve = (
    bb_on_curve**2*(2*parent_x+horizontal_x-slope**2)/square_on_curve
)
assert W_on_curve**2 == evaluate_bivariate(quartic, K(T), new_base)

# Invert the Mobius new-base map and apply the standard pointed-quartic map.
V = PolynomialRing(F, "v")
v = V.gen()
L = V.fraction_field()
numerator = new_base.numerator()
denominator = new_base.denominator()
n0 = numerator[0]
n1 = numerator[1] if numerator.degree() else F.zero()
d0 = denominator[0]
d1 = denominator[1] if denominator.degree() else F.zero()
old_base_of_v = L(n0-v*d0)/L(v*d1-n1)


def substitute_old_base(function):
    function = K(function)

    def evaluate(poly):
        answer = L.zero()
        for coefficient in reversed(R(poly).list()):
            answer = answer*old_base_of_v+L(coefficient)
        return answer

    return evaluate(function.numerator())/evaluate(function.denominator())


def evaluate_u_record(record):
    numerator_v = V([reduce_qq(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator_v = V([reduce_qq(value) for value in record["denominator_coefficients_low_to_high"]])
    return L(numerator_v)/L(denominator_v)


old_t_point = old_base_of_v
W_point = substitute_old_base(W_on_curve)
quartic_coefficients_v = []
for coefficient in TUP(quartic).list():
    quartic_coefficients_v.append(L(V([F(value) for value in U(coefficient).list()])))
p1229_support = F(QQ(q8["preferred_pointed_zero"]["old_base_coordinate"]))
translated = []
for new_degree in range(5):
    translated.append(sum(
        quartic_coefficients_v[old_degree]*F(comb(old_degree, new_degree))
        * L(p1229_support)**(old_degree-new_degree)
        for old_degree in range(new_degree, 5)
    ))
e, d, c, b, a = translated
q_origin = evaluate_u_record(q8["preferred_pointed_zero"]["quartic_ordinate"])
assert e == q_origin**2
a1 = d/q_origin
a2 = c-d**2/(4*q_origin**2)
a3 = 2*q_origin*b
b2 = a1**2+4*a2

child_A = V([reduce_qq(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]])
child_B = V([reduce_qq(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]])


def pointed_image(ordinate):
    local_x = old_t_point-L(p1229_support)
    x_generalized = (2*q_origin*(ordinate+q_origin)+d*local_x)/local_x**2
    y_generalized = (
        4*q_origin**2*(ordinate+q_origin)
        +2*q_origin*(d*local_x+c*local_x**2)
        -d**2*local_x**2/(2*q_origin)
    )/local_x**3
    answer_x = 9*(x_generalized+b2/12)
    answer_y = 27*(y_generalized+(a1*x_generalized+a3)/2)
    assert answer_y**2 == answer_x**3+L(child_A)*answer_x+L(child_B)
    return answer_x, answer_y


# Factoring the square quotient determines its root only up to sign.  The
# marked P.O=0 condition selects the sign corresponding to the named curve.
pointed_trials = [pointed_image(W_point), pointed_image(-W_point)]
polynomial_trials = [
    (index, trial) for index, trial in enumerate(pointed_trials)
    if trial[0].denominator().degree() == 0
    and trial[1].denominator().degree() == 0
]
assert len(polynomial_trials) == 1
square_root_sign_index, (child_x, child_y) = polynomial_trials[0]
child_x_poly = V(child_x)
child_y_poly = V(child_y)
assert child_x_poly.degree() <= 4 and child_y_poly.degree() <= 6

# The lattice profile says identity at all four I2 fibres; verify the three
# finite nodes and the scaled node at infinity literally.  Since every entry
# is zero, no ordering convention is needed here.
component_profile = []
RX = PolynomialRing(F, "x")
xvar = RX.gen()
child_delta = -F(16)*(4*child_A**3+27*child_B**2)
finite_i2_supports = sorted(
    (-factor[0]/factor[1] for factor, exponent in child_delta.factor()
     if int(exponent) == 2 and factor.degree() == 1),
    key=int,
)
assert len(finite_i2_supports) == 3
for support in finite_i2_supports:
    cubic = xvar**3+child_A(support)*xvar+child_B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    node = -repeated[0]/repeated[1]
    component_profile.append(int(child_x_poly(support) == node and child_y_poly(support) == 0))
infinity_A = child_A[8]
infinity_B = child_B[12]
infinity_cubic = xvar**3+infinity_A*xvar+infinity_B
infinity_repeated = infinity_cubic.gcd(infinity_cubic.derivative())
assert infinity_repeated.degree() == 1
infinity_node = -infinity_repeated[0]/infinity_repeated[1]
infinity_x = child_x_poly[4] if child_x_poly.degree() == 4 else F.zero()
infinity_y = child_y_poly[6] if child_y_poly.degree() == 6 else F.zero()
component_profile.append(int(infinity_x == infinity_node and infinity_y == 0))
assert component_profile == branch["current_4A1_component_pairings"]

payload = {
    "schema": "elkies-k3.h92-q12o5867-degree1-compiler-branch-mod131.v1",
    "status": "PASS_MOD131_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD",
    "prime": int(prime),
    "compiler_branch": {
        "index_in_promoted_four_section_word": compiler["new_sections"].index(branch),
        "q4o164_parent_curve": list(map(int, parent_curve)),
        "q4o164_parent_degree": 1,
        "q4o164_parent_a_minus_b": 1,
        "current_4A1_P_dot_O": 0,
        "current_4A1_component_pairings": component_profile,
    },
    "exact_lattice_to_equation_word": {
        "selected_marked_embedding_index": 15,
        "C8opposite_coefficient": int(c8_coefficient),
        "B0_through_B7_coefficients": list(map(int, integral_B_word)),
        "identity": "2*C8opp+B1+2*B2-2*B3+B4-2*B5+B6+B7",
        "exact_marked_MW9_tail": [str(value) for value in target_tail],
    },
    "parent_section_mod131": {
        "x_numerator_coefficients_low_to_high": list(map(int, parent_x.numerator())),
        "x_denominator_coefficients_low_to_high": list(map(int, parent_x.denominator())),
        "y_numerator_coefficients_low_to_high": list(map(int, parent_y.numerator())),
        "y_denominator_coefficients_low_to_high": list(map(int, parent_y.denominator())),
        "degrees_x_num_x_den_y_num_y_den": [
            int(parent_x.numerator().degree()), int(parent_x.denominator().degree()),
            int(parent_y.numerator().degree()), int(parent_y.denominator().degree()),
        ],
        "exact_weierstrass_identity_mod131": True,
    },
    "q8_restriction_mod131": {
        "new_base_numerator_coefficients_low_to_high": list(map(int, new_base.numerator())),
        "new_base_denominator_coefficients_low_to_high": list(map(int, new_base.denominator())),
        "new_base_degree": 1,
        "square_root_sign_index": int(square_root_sign_index),
        "quartic_ordinate_numerator_degree": int(W_on_curve.numerator().degree()),
        "quartic_ordinate_denominator_degree": int(W_on_curve.denominator().degree()),
        "exact_quartic_identity_mod131": True,
    },
    "q8_child_section_mod131": {
        "x_coefficients_low_to_high": list(map(int, child_x_poly.list())),
        "y_coefficients_low_to_high": list(map(int, child_y_poly.list())),
        "degrees_x_y": [int(child_x_poly.degree()), int(child_y_poly.degree())],
        "P_dot_O": 0,
        "exact_weierstrass_identity_mod131": True,
    },
    "method": {
        "large_Groebner_required": False,
        "nonlinear_system_solve_required": False,
        "construction": "exact marked group word, mod-131 elliptic group law, resolved-q8 restriction, pointed-quartic map",
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "This exactly identifies the promoted degree-one parent class with stored equation points and constructs its "
        "P.O=0 q8-child section at the pinned good prime. It is a modular seed only; characteristic-zero reconstruction, "
        "the other three compiler branches, and q12 resolved RR remain open."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867DEG1MOD131|parent_degrees=6,2,8,3|child_degrees={},{}|profile={}|"
    "status={}|output={}".format(
        child_x_poly.degree(), child_y_poly.degree(), component_profile,
        payload["status"], OUTPUT,
    ),
    flush=True,
)
