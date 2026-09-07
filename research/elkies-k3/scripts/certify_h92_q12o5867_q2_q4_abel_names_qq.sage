#!/usr/bin/env sage -python
"""Certify the exact QQ Abel names of the lifted q12/o5867 Q2 and Q4 seeds.

Invert each exact degree-(4,6) q8-child section through the pointed binary
quartic and the resolved q4 chord.  Its old-base map T(u) has degree two.
Over QQ(s), reduce the recovered parent point modulo

    numerator(T(u)) - s*denominator(T(u)),

use the explicit quadratic conjugation, and add the two conjugate parent
points.  The resulting point is compared literally with the requested exact
word on the q4/o164 parent.  This is symmetric quadratic arithmetic only; no
Groebner basis or elimination is used.

The script also evaluates both lifted sections at every exact q8 I2 node.
"""

import argparse
import hashlib
import json
import time
from math import comb
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
LIFTS = LOCAL / "q12o5867-abel-trace-named-seeds-qq.json"
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
SHELL = LOCAL / "q12o5867-p0-shell-word-fingerprints-mod89.json"
INPUTS = (Q8, MODEL, BASIS, C8, HORIZONTAL, LIFTS, Q3, SHELL)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q12o5867-q2-q4-abel-names-qq.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def function_record(function):
    function = function.parent()(function)
    numerator = function.numerator()
    denominator = function.denominator()
    values = list(numerator.list())+list(denominator.list())
    return {
        "numerator_coefficients_low_to_high": [str(value) for value in numerator.list()],
        "denominator_coefficients_low_to_high": [str(value) for value in denominator.list()],
        "numerator_denominator_degrees": [int(numerator.degree()), int(denominator.degree())],
        "maximum_rational_bits": max(map(coefficient_bits, values)),
    }


started = time.monotonic()
q8 = json.loads(Q8.read_text())
model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
c8 = json.loads(C8.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
lifts = json.loads(LIFTS.read_text())
q3_exact = json.loads(Q3.read_text())
shell = json.loads(SHELL.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert c8["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert lifts["status"] == "PASS_EXACT_QQ_Q12O5867_Q2_AND_Q4_MODULAR_ABEL_SEED_SECTIONS"
assert q3_exact["status"] == "PASS_EXACT_QQ_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD"

PT = PolynomialRing(QQ, "T")
T = PT.gen()
PK = PT.fraction_field()
U = PolynomialRing(QQ, "u")
u = U.gen()
KU = U.fraction_field()
TU = PolynomialRing(KU, "T")
Tvar = TU.gen()
TUP = PolynomialRing(U, "T")


def parent_poly(values):
    return PT([QQ(value) for value in values])


def child_poly(values):
    return U([QQ(value) for value in values])


def parent_function(record):
    return PK(parent_poly(record["numerator_coefficients_low_to_high"])) / PK(
        parent_poly(record["denominator_coefficients_low_to_high"])
    )


def child_function(record):
    return KU(child_poly(record["numerator_coefficients_low_to_high"])) / KU(
        child_poly(record["denominator_coefficients_low_to_high"])
    )


parent_A = parent_poly(model["compact_model"]["A_coefficients_low_to_high"])
parent_B = parent_poly(model["compact_model"]["B_coefficients_low_to_high"])
child_A = child_poly(q8["child"]["minimal_A_coefficients_low_to_high"])
child_B = child_poly(q8["child"]["minimal_B_coefficients_low_to_high"])


def checked_parent_point(x_coordinate, y_coordinate):
    point = PK(x_coordinate), PK(y_coordinate)
    assert point[1]**2 == point[0]**3+PK(parent_A)*point[0]+PK(parent_B)
    return point


def parent_neg(point):
    return None if point is None else (point[0], -point[1])


def parent_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+PK(parent_A))/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return checked_parent_point(x3, slope*(x1-x3)-y1)


def parent_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return parent_mul(-coefficient, parent_neg(point))
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = parent_add(answer, addend)
        addend = parent_add(addend, addend)
        coefficient >>= 1
    return answer


parent_basis = [
    checked_parent_point(
        parent_poly(section["x_coefficients_low_to_high"]),
        parent_poly(section["y_coefficients_low_to_high"]),
    )
    for section in basis["resolved_hensel"]["sections"]
]
H = checked_parent_point(
    parent_function(horizontal["section"]["x"]),
    parent_function(horizontal["section"]["y"]),
)
c8_record = c8["opposite_constant_support_section"]
c8_old = parent_function(c8_record["x"]), parent_function(c8_record["y"])
base_scale = QQ(model["exact_coordinate_change"]["c"])
xy_scale = QQ(model["exact_coordinate_change"]["s"])


def substitute_scaled(function):
    return PK(PT(function.numerator()(base_scale*T))) / PK(
        PT(function.denominator()(base_scale*T))
    )


C8opposite = checked_parent_point(
    substitute_scaled(c8_old[0])/xy_scale**2,
    substitute_scaled(c8_old[1])/xy_scale**3,
)


def word_point(terms):
    answer = None
    for coefficient, point in terms:
        answer = parent_add(answer, parent_mul(coefficient, point))
    return answer


word_formulas = {
    "Q2": "H+B0+B2-B3+B4",
    "Q4": "H-C8opp-B1-B2+B3+B5-B6-B7",
    "Q3": "2C8opp+B1+2B2-2B3+B4-2B5+B6+B7",
}
target_words = {
    "Q2": word_point((
        (1, H), (1, parent_basis[0]), (1, parent_basis[2]),
        (-1, parent_basis[3]), (1, parent_basis[4]),
    )),
    "Q4": word_point((
        (1, H), (-1, C8opposite), (-1, parent_basis[1]),
        (-1, parent_basis[2]), (1, parent_basis[3]),
        (1, parent_basis[5]), (-1, parent_basis[6]), (-1, parent_basis[7]),
    )),
    "Q3": word_point((
        (2, C8opposite), (1, parent_basis[1]), (2, parent_basis[2]),
        (-2, parent_basis[3]), (1, parent_basis[4]), (-2, parent_basis[5]),
        (1, parent_basis[6]), (1, parent_basis[7]),
    )),
}

# Reconstruct the exact normalized quartic square factor from the stored RR
# basis.  This is the same checked identity used by the modular Abel script.
XH = parent_poly(horizontal["section"]["x"]["numerator_coefficients_low_to_high"])
denominator_x = parent_poly(horizontal["section"]["x"]["denominator_coefficients_low_to_high"])
YH = parent_poly(horizontal["section"]["y"]["numerator_coefficients_low_to_high"])
denominator_y = parent_poly(horizontal["section"]["y"]["denominator_coefficients_low_to_high"])
Z = PT.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic()**(int(exponent)//2)
assert Z**2 == denominator_x and Z**3 == denominator_y
resolved_pairs = [
    (
        parent_poly(item["AA_coefficients_low_to_high"]),
        parent_poly(item["BB_coefficients_low_to_high"]),
    )
    for item in q8["resolved_RR"]["resolved_basis_pairs"]
]
(AA0, BB0), (AA1, BB1) = resolved_pairs


def lift_TU(poly):
    return TU([KU(value) for value in PT(poly).list()])


aa = lift_TU(AA0)+KU(u)*lift_TU(AA1)
bb = lift_TU(BB0)+KU(u)*lift_TU(BB1)
X_u, Y_u, Z_u, A_u = map(lift_TU, (XH, YH, Z, parent_A))
raw = (
    aa**4-6*X_u*aa**2*bb**2+8*Y_u*aa*bb**3
    -3*X_u**2*bb**4-4*A_u*bb**4*Z_u**4
)
after_collision, remainder = raw.quo_rem(Z_u**4)
assert not remainder
quartic = TUP([
    U([QQ(value) for value in coefficient])
    for coefficient in q8["quartic"]["coefficients_in_old_T_low_to_high"]
])
square_quotient, remainder = TUP(after_collision).quo_rem(quartic)
assert not remainder
factorization = square_quotient.factor()
assert all(int(exponent) % 2 == 0 for unused, exponent in factorization)
unit = QQ(factorization.unit())
assert unit.is_square()
square_factor = TUP(unit.sqrt())
for factor, exponent in factorization:
    square_factor *= factor**(int(exponent)//2)
assert TUP(after_collision) == quartic*square_factor**2


def evaluate_bivariate(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(TUP(poly).list()):
        answer = answer*old_t+KU(U(coefficient))
    return answer


def evaluate_parent_poly(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(PT(poly).list()):
        answer = answer*old_t+KU(coefficient)
    return answer


def compose_parent(function, old_t):
    function = PK(function)
    return evaluate_parent_poly(function.numerator(), old_t) / evaluate_parent_poly(
        function.denominator(), old_t
    )


pointing = q8["preferred_pointed_zero"]
q_origin = child_function(pointing["quartic_ordinate"])
p1229_support = QQ(pointing["old_base_coordinate"])
quartic_coefficients = [KU(U(coefficient)) for coefficient in TUP(quartic).list()]
translated = []
for new_degree in range(5):
    translated.append(sum(
        quartic_coefficients[old_degree]*QQ(comb(old_degree, new_degree))
        * KU(p1229_support)**(old_degree-new_degree)
        for old_degree in range(new_degree, 5)
    ))
e, d, c, b, a = translated
assert e == q_origin**2
a1 = d/q_origin
a2 = c-d**2/(4*q_origin**2)
a3 = 2*q_origin*b
b2 = a1**2+4*a2


def invert_child_section(X, Y):
    X = KU(X)
    Y = KU(Y)
    x_general = X/9-b2/12
    y_general = Y/27-(a1*x_general+a3)/2
    assert y_general
    local_t = 2*q_origin*(x_general+a2)/y_general
    old_t = KU(p1229_support)+local_t
    ordinate = (x_general*local_t**2-d*local_t)/(2*q_origin)-q_origin
    assert ordinate**2 == evaluate_bivariate(quartic, old_t)
    AA = evaluate_parent_poly(AA0, old_t)+KU(u)*evaluate_parent_poly(AA1, old_t)
    BB = evaluate_parent_poly(BB0, old_t)+KU(u)*evaluate_parent_poly(BB1, old_t)
    ZZ = evaluate_parent_poly(Z, old_t)
    slope = -AA/(BB*ZZ)
    square = evaluate_bivariate(square_factor, old_t)
    hx = compose_parent(H[0], old_t)
    hy = compose_parent(H[1], old_t)
    old_x = (ordinate*square/BB**2-hx+slope**2)/2
    old_y = slope*(old_x-hx)-hy
    assert old_y**2 == old_x**3+compose_parent(parent_A, old_t)*old_x+compose_parent(parent_B, old_t)
    restrictions = [
        evaluate_parent_poly(AAi, old_t)+evaluate_parent_poly(BBi, old_t)*ZZ*slope
        for AAi, BBi in resolved_pairs
    ]
    assert KU(u) == -restrictions[0]/restrictions[1]
    return old_t, ordinate, old_x, old_y


def rational_degree(value):
    value = KU(value)
    return max(value.numerator().degree(), value.denominator().degree())


# The factorization determines the square factor only up to sign.  Pin that
# sign over QQ using the already-certified degree-one Q3 inverse, rather than
# inheriting an implementation-dependent square-root choice.
q3_section = q3_exact["section"]
Q3_child = (
    child_poly(q3_section["x_coefficients_low_to_high"]),
    child_poly(q3_section["y_coefficients_low_to_high"]),
)


def q3_inverse_matches():
    q3_old_t, unused_q3_ordinate, q3_old_x, q3_old_y = invert_child_section(*Q3_child)
    return (
        rational_degree(q3_old_t) == 1
        and q3_old_x == compose_parent(target_words["Q3"][0], q3_old_t)
        and q3_old_y == compose_parent(target_words["Q3"][1], q3_old_t)
    )


square_factor_sign_from_initial_sqrt = 1
if not q3_inverse_matches():
    square_factor = -square_factor
    square_factor_sign_from_initial_sqrt = -1
    assert q3_inverse_matches()


# Symmetric quadratic trace over QQ(s), without constructing a number-field
# extension.  Every rational function is reduced to x0+x1*u modulo the monic
# quadratic defining the two conjugate inverse points.
S = PolynomialRing(QQ, "s")
s = S.gen()
KS = S.fraction_field()
US = PolynomialRing(KS, "u")
us = US.gen()


def lift_U_to_US(poly):
    return US([KS(value) for value in U(poly).list()])


def reduce_KU_mod_quadratic(function, quadratic):
    function = KU(function)
    numerator = lift_U_to_US(function.numerator())
    denominator = lift_U_to_US(function.denominator())
    inverse = denominator.inverse_mod(quadratic)
    answer = (numerator*inverse) % quadratic
    assert answer.degree() <= 1
    return answer


def parent_function_at_s(function):
    function = PK(function)

    def evaluate(poly):
        answer = KS.zero()
        for coefficient in reversed(PT(poly).list()):
            answer = answer*KS(s)+KS(coefficient)
        return answer

    return evaluate(function.numerator())/evaluate(function.denominator())


def quadratic_abel_trace(old_t, old_x, old_y):
    old_t = KU(old_t)
    equation = lift_U_to_US(old_t.numerator())-KS(s)*lift_U_to_US(old_t.denominator())
    assert equation.degree() == 2
    equation = equation/equation[2]
    discriminant = equation[1]**2-4*equation[0]
    assert discriminant
    x_reduced = reduce_KU_mod_quadratic(old_x, equation)
    y_reduced = reduce_KU_mod_quadratic(old_y, equation)
    assert x_reduced.degree() == 1 and x_reduced[1]
    conjugate_root_sum = -equation[1]
    chord_slope = y_reduced[1]/x_reduced[1]
    x_trace = chord_slope**2-(2*x_reduced[0]+x_reduced[1]*conjugate_root_sum)
    y_trace_reduced = (
        chord_slope*(x_reduced-US(x_trace))-y_reduced
    ) % equation
    assert y_trace_reduced.degree() <= 0
    y_trace = KS(y_trace_reduced[0])
    assert y_trace**2 == x_trace**3+parent_function_at_s(parent_A)*x_trace+parent_function_at_s(parent_B)
    return {
        "quadratic": equation,
        "discriminant": discriminant,
        "x_reduced": x_reduced,
        "y_reduced": y_reduced,
        "chord_slope": chord_slope,
        "trace": (KS(x_trace), y_trace),
    }


def i2_profile(X, Y):
    X, Y = U(X), U(Y)
    entries = []
    profile = []
    for fibre in q8["child"]["finite_reducible_fibres"]:
        factor = U(fibre["factor"])
        assert factor.degree() == 1
        support = -factor[0]/factor[1]
        node = -3*child_B(support)/(2*child_A(support))
        x_value, y_value = X(support), Y(support)
        meets_node = x_value == node and y_value == 0
        profile.append(int(meets_node))
        entries.append({
            "support": str(support),
            "node_x": str(node),
            "section_x": str(x_value),
            "section_y": str(y_value),
            "meets_singular_node": bool(meets_node),
        })
    infinity_node = -3*child_B[12]/(2*child_A[8])
    infinity_x = X[4] if X.degree() == 4 else QQ.zero()
    infinity_y = Y[6] if Y.degree() == 6 else QQ.zero()
    infinity_meets = infinity_x == infinity_node and infinity_y == 0
    profile.append(int(infinity_meets))
    entries.append({
        "support": "infinity",
        "node_x": str(infinity_node),
        "section_x": str(infinity_x),
        "section_y": str(infinity_y),
        "meets_singular_node": bool(infinity_meets),
    })
    return profile, entries


specs = (
    (
        "Q2", "Q2_unique_rank12_seed", [0, 0, 0, 0], [0, 0, 0, 0],
    ),
    (
        "Q4", "Q4_candidate1_shell220_rank12_seed", [1, 0, 0, 0], [0, 1, 0, 0],
    ),
)
certificates = {}
for name, lift_key, expected_equation_profile, expected_marked_profile in specs:
    section_record = lifts["sections"][lift_key]["section"]
    X = child_poly(section_record["x_coefficients_low_to_high"])
    Y = child_poly(section_record["y_coefficients_low_to_high"])
    assert Y**2 == X**3+child_A*X+child_B
    old_t, ordinate, old_x, old_y = invert_child_section(X, Y)
    assert rational_degree(old_t) == 2
    trace_data = quadratic_abel_trace(old_t, old_x, old_y)
    target = tuple(parent_function_at_s(value) for value in target_words[name])
    if trace_data["trace"] != target:
        trace_x, trace_y = trace_data["trace"]
        print(
            "Q12O5867ABELQQDIAG|name={}|x_equal={}|y_equal={}|negative_equal={}|"
            "trace_degrees={}/{},{}//{}|target_degrees={}/{},{}//{}".format(
                name, int(trace_x == target[0]), int(trace_y == target[1]),
                int(trace_x == target[0] and trace_y == -target[1]),
                trace_x.numerator().degree(), trace_x.denominator().degree(),
                trace_y.numerator().degree(), trace_y.denominator().degree(),
                target[0].numerator().degree(), target[0].denominator().degree(),
                target[1].numerator().degree(), target[1].denominator().degree(),
            ), flush=True,
        )
        raise ArithmeticError("{} exact quadratic Abel trace missed target word".format(name))
    equation_profile, node_evaluations = i2_profile(X, Y)
    assert equation_profile == expected_equation_profile
    shell_branch = shell["branches"][name]
    assert shell_branch["expected_equation_component_profile"] == expected_equation_profile
    assert shell_branch["expected_marked_component_profile"] == expected_marked_profile

    certificates[name] = {
        "lift_key": lift_key,
        "target_parent_word": word_formulas[name],
        "inverse": {
            "old_base_map": function_record(old_t),
            "old_base_map_degree": 2,
            "pointed_quartic_identity": True,
            "recovered_parent_weierstrass_identity": True,
            "resolved_RR_parameter_identity": True,
        },
        "quadratic_Abel_trace": {
            "defining_polynomial_coefficients_low_to_high_in_u": [
                str(trace_data["quadratic"][index]) for index in range(3)
            ],
            "defining_polynomial_discriminant_nonzero": True,
            "conjugate_addition_descends_to_QQ_s": True,
            "literal_target_parent_word_identity": True,
            "trace_x": function_record(trace_data["trace"][0]),
            "trace_y": function_record(trace_data["trace"][1]),
        },
        "q8_I2_profile": {
            "support_order": [
                str(-U(fibre["factor"])[0]/U(fibre["factor"])[1])
                for fibre in q8["child"]["finite_reducible_fibres"]
            ]+["infinity"],
            "exact_equation_profile": equation_profile,
            "expected_marked_profile": expected_marked_profile,
            "equation_to_marked_profile_transport_from_mod89_artifact": True,
            "node_evaluations": node_evaluations,
        },
        "name_status": (
            "exact_QQ_Abel_name_and_q8_I2_profile"
            if name == "Q2" else
            "exact_QQ_Abel_name_and_q8_I2_profile_for_candidate1_shell220"
        ),
        "class_scope": {
            "exact_q8_section_curve_certified": True,
            "generic_parent_MW_Abel_class_certified": True,
            "full_parent_multisection_divisor_NS_class_certified": False,
            "resolved_q4_vertical_component_data_computed": False,
            "vertical_component_ambiguity": (
                "No second p=89 Abel-named Q2 seed survives in the recorded shell, but the "
                "parent degree-two multisection's resolved q4 vertical correction was not "
                "computed, so the Abel trace alone does not exclude an NS-vertical translate."
                if name == "Q2" else
                "Unresolved: p=89 also has the Abel-equivalent Q4 candidate7/shell269, and "
                "the present exact trace/profile data do not compare their resolved q4 "
                "vertical corrections."
            ),
        },
    }
    print(
        "Q12O5867ABELQQ|name={}|degree=2|profile={}|target_word=1|status=PASS".format(
            name, ",".join(map(str, equation_profile))
        ),
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h92-q12o5867-q2-q4-abel-names-qq.v1",
    "status": "PASS_EXACT_QQ_Q12O5867_Q2_Q4_ABEL_NAMES_AND_Q8_I2_PROFILES",
    "certificates": certificates,
    "method": {
        "construction": "exact pointed-quartic inverse, resolved chord inverse, and symmetric quadratic conjugate addition over QQ(s)",
        "quartic_square_factor_sign_pinned_by_exact_Q3_inverse": True,
        "square_factor_sign_relative_to_initial_sqrt": square_factor_sign_from_initial_sqrt,
        "large_Groebner_required": False,
        "elimination_required": False,
        "finite_fibre_sampling_used_for_names": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "The Q2 and requested Q4 candidate1/shell220 lifts now have literal exact QQ Abel-trace "
        "identities and exact q8 I2 node profiles. The Q4 result certifies that lifted branch's "
        "generic parent Abel word, but this artifact does not yet compute resolved q4 multisection "
        "component/vertical intersection data; consequently it does not assert that Abel-equivalent "
        "multisections with the same q8 profile are globally identical in NS. The q12 resolved RR "
        "pencil and endpoint equation remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ABELQQ|names=2|status={}|output={}".format(
        payload["status"], OUTPUT
    ),
    flush=True,
)
