#!/usr/bin/env sage -python
"""Lift q4/orbit323 via resolved polynomial summands and exact group law.

The compact source has bad auxiliary reduction at p=59, making the direct
simple-pole chart singular.  Its one-node and two-node polynomial charts are
resolved and mostly regular.  Lift those smaller charts separately, add the
prescribed modular pairs over Q_59(t), and rationally reconstruct only their
simple-pole sums.  All retained sections are verified by literal substitution
over QQ.  No elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
ONE_NODE = LOCAL / "q4o208-q4o323-one-node-f1-mod59.json"
SUMS = LOCAL / "q4o208-q4o323-horizontal-sums-mod59.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=500)
parser.add_argument("--output", type=Path, default=LOCAL / "q4o208-q4o323-horizontal-via-summands-qq.json")
args = parser.parse_args()
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
PRECISION = int(args.precision)
PRIME = ZZ(59)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


compact = json.loads(COMPACT.read_text())
one_node = json.loads(ONE_NODE.read_text())
sums = json.loads(SUMS.read_text())
marking = json.loads(MARKING.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert one_node["status"] == "PASS_MODP_Q4O323_ONE_NODE_ADJACENT_COMPONENT_SECTION_SCAN"
assert sums["status"] == "PASS_MOD59_Q4O323_SIMPLE_POLE_SUM_SEEDS"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"

RQ = PolynomialRing(QQ, "t")
tq = RQ.gen()
A_QQ = RQ([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
node_zero_QQ = QQ(3)
node_one_QQ = -3*B_QQ(1)/(2*A_QQ(1))
node_infinity_QQ = -3*B_QQ[12]/(2*A_QQ[8])
assert A_QQ(1) == -3*node_one_QQ**2 and B_QQ(1) == 2*node_one_QQ**3
assert A_QQ[8] == -3*node_infinity_QQ**2 and B_QQ[12] == 2*node_infinity_QQ**3


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


F = GF(PRIME)
RF = PolynomialRing(F, "t")
tf = RF.gen()
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
node_zero_F = reduce_qq(node_zero_QQ, F)
node_one_F = reduce_qq(node_one_QQ, F)
node_infinity_F = reduce_qq(node_infinity_QQ, F)
assert node_one_F == 0

K = Qp(PRIME, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
node_zero = K(node_zero_QQ)
node_one = K(node_one_QQ)
node_infinity = K(node_infinity_QQ)


def one_unpack(values, ring, node):
    variable = ring.gen()
    X = ring([node] + list(values[:4]))
    Q = ring(list(values[4:10]))
    return X, variable*Q


def two_unpack(values, ring, node_one_value, node_infinity_value):
    variable = ring.gen()
    x0, x1, x2 = values[:3]
    x3 = node_one_value-node_infinity_value-x0-x1-x2
    X = ring([x0, x1, x2, x3, node_infinity_value])
    Q = ring(list(values[3:8]))
    return X, (variable-1)*Q


def equation_residual(X, Y, surface_A, surface_B, field):
    equation = Y**2-X**3-surface_A*X-surface_B
    return vector(field, [equation[degree] for degree in range(13)])


def one_jacobian(values, ring, surface_A, node):
    X, Y = one_unpack(values, ring, node)
    variable = ring.gen()
    dx, dy = -3*X**2-surface_A, 2*Y
    derivatives = [dx*variable**degree for degree in range(1, 5)]
    derivatives.extend(dy*variable**degree for degree in range(1, 7))
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), 13, 10, lambda row, column: (
        derivatives[column][row] if row <= derivatives[column].degree() else zero
    ))


def two_jacobian(values, ring, surface_A, node_one_value, node_infinity_value):
    X, Y = two_unpack(values, ring, node_one_value, node_infinity_value)
    variable = ring.gen()
    dx, dy = -3*X**2-surface_A, 2*Y
    derivatives = [dx*(variable**degree-variable**3) for degree in range(3)]
    derivatives.extend(dy*(variable-1)*variable**degree for degree in range(5))
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), 13, 8, lambda row, column: (
        derivatives[column][row] if row <= derivatives[column].degree() else zero
    ))


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def newton_lift(seed, unpack_padic, jacobian_finite, jacobian_padic):
    J_F = jacobian_finite(seed)
    variable_count = len(seed)
    if J_F.rank() != variable_count:
        return None
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    determinant = int(matrix(F, [J_F.row(row) for row in pivot_rows]).det())
    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < PRECISION:
        working_precision = min(2*known_precision, PRECISION)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        X, Y = unpack_padic(values)
        full = equation_residual(X, Y, A, B, K)
        square = matrix(K, [jacobian_padic(values).row(row) for row in pivot_rows])
        correction = square.solve_right(-vector(K, [full[row] for row in pivot_rows]))
        values += correction
        X_after, Y_after = unpack_padic(values)
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_residual_valuation_after": int(minimum_valuation(
                equation_residual(X_after, Y_after, A, B, K)
            )),
        })
        known_precision = working_precision
    return values, pivot_rows, determinant, iterations


def transform_one_record(record):
    change = compact["exact_coordinate_change"]
    a, b, c, d, m = [reduce_qq(change[key], F) for key in ("a", "b", "c", "d", "m")]
    N, D = a*tf+b, c*tf+d
    X = RF(sum(F(value)*N**degree*D**(4-degree)
               for degree, value in enumerate(record["x_coefficients_low_to_high"]))/m**2)
    Y = RF(sum(F(value)*N**degree*D**(6-degree)
               for degree, value in enumerate(record["y_coefficients_low_to_high"]))/m**3)
    return X, Y


one_finite = [transform_one_record(record) for record in one_node["sections"]]
one_representatives = {}
one_lifts = {}
one_diagnostics = []
for X, Y in one_finite:
    key = tuple(X.list())
    if key in one_representatives:
        continue
    Q, remainder = Y.quo_rem(tf)
    assert remainder == 0
    seed = vector(F, X.list()[1:5] + [F.zero()]*(5-len(X.list()))
                  + Q.list() + [F.zero()]*(6-len(Q.list())))
    result = newton_lift(
        seed,
        lambda values: one_unpack(values, RT, node_zero),
        lambda values: one_jacobian(values, RF, A_F, node_zero_F),
        lambda values: one_jacobian(values, RT, A, node_zero),
    )
    assert result is not None
    values, rows, determinant, iterations = result
    X_padic, Y_padic = one_unpack(values, RT, node_zero)
    one_representatives[key] = Y
    one_lifts[key] = (X_padic, Y_padic)
    one_diagnostics.append({"pivot_rows": rows, "determinant_mod59": determinant, "iterations": iterations})

KT = RT.fraction_field()
E_padic = EllipticCurve(KT, [0, 0, 0, KT(A), KT(B)])

# The exact transported C7 section has relative physical component labels
# (-1,+1,2) at compact (t=0,t=1,infinity).  Hence the marked target's
# relative labels (-1,-1,-1) have multipliers (1,3,*) against C7 at 0,1.
c7_record = compact["transported_exact_section"]
C7_X = RT([K(QQ(value)) for value in c7_record["x_coefficients_low_to_high"]])
C7_Y = RT([K(QQ(value)) for value in c7_record["y_coefficients_low_to_high"]])
C7 = E_padic.point([KT(C7_X), KT(C7_Y), KT(1)], check=False)


def transform_exact_old_coordinate(record, weight):
    change = compact["exact_coordinate_change"]
    a_QQ, b_QQ, c_QQ, d_QQ, m_QQ = [QQ(change[key]) for key in ("a", "b", "c", "d", "m")]
    N = K(a_QQ)*t+K(b_QQ)
    D = K(c_QQ)*t+K(d_QQ)
    numerator = [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    denominator = [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    numerator_value = sum(K(value)*N**degree*D**(weight-degree)
                          for degree, value in enumerate(numerator))
    denominator_degree = len(denominator)-1
    denominator_value = sum(K(value)*N**degree*D**(denominator_degree-degree)
                            for degree, value in enumerate(denominator))
    return KT(numerator_value)/(K(m_QQ)**(weight//2)*KT(denominator_value))


inherited_labels = (
    "first_I6_affine_component_on_C5_pointed_child",
    "old_A11_component_7_on_C5_pointed_child",
    "second_I6_affine_component_on_C5_pointed_child",
)
inherited_points = []
for label in inherited_labels:
    record = marking[label]
    x_value = transform_exact_old_coordinate(record["x"], 4)
    y_value = transform_exact_old_coordinate(record["y"], 6)
    inherited_points.append(E_padic.point([x_value, y_value, KT(1)], check=False))
expected_inherited_intersections = (3, 2, 2)


def section_intersection(left, right):
    difference = left-right
    if difference.is_zero():
        return -2
    x_value = KT(difference[0])
    numerator_degree = int(x_value.numerator().degree())
    denominator_degree = int(x_value.denominator().degree())
    if denominator_degree % 2:
        return None
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    if infinity_excess % 2:
        return None
    return denominator_degree//2 + infinity_excess//2


def finite_node_hit_padic(point, support, node):
    if point.is_zero():
        return False
    x_value, y_value = map(KT, point[:2])
    x_denominator = x_value.denominator()(support)
    y_denominator = y_value.denominator()(support)
    if not x_denominator or not y_denominator:
        return False
    x_residual = x_value(support)-node
    y_residual = y_value(support)
    threshold = PRECISION-20
    x_valuation = x_residual.valuation() if x_residual else PRECISION
    y_valuation = y_residual.valuation() if y_residual else PRECISION
    return x_valuation >= threshold and y_valuation >= threshold


one_points = []
for X, Y in one_finite:
    key = tuple(X.list())
    X_padic, Y_padic = one_lifts[key]
    sign = 1 if Y == one_representatives[key] else -1
    one_points.append(E_padic.point([KT(X_padic), KT(sign*Y_padic), KT(1)], check=False))


def finite_polynomial(record):
    denominator = RF(record["denominator_coefficients_low_to_high"])
    assert denominator.degree() == 0
    return RF(record["numerator_coefficients_low_to_high"])/denominator[0]


two_finite = [
    (finite_polynomial(record["x"]), finite_polynomial(record["y"]))
    for record in sums["two_node_sections"]
]
two_representatives = {}
two_lifts = {}
two_diagnostics = []
singular_two_x = 0
for X, Y in two_finite:
    key = tuple(X.list())
    if key in two_representatives:
        continue
    Q, remainder = Y.quo_rem(tf-1)
    assert remainder == 0
    seed = vector(F, [X[index] for index in range(3)]
                  + Q.list() + [F.zero()]*(5-len(Q.list())))
    result = newton_lift(
        seed,
        lambda values: two_unpack(values, RT, node_one, node_infinity),
        lambda values: two_jacobian(values, RF, A_F, node_one_F, node_infinity_F),
        lambda values: two_jacobian(values, RT, A, node_one, node_infinity),
    )
    if result is None:
        singular_two_x += 1
        continue
    values, rows, determinant, iterations = result
    X_padic, Y_padic = two_unpack(values, RT, node_one, node_infinity)
    two_representatives[key] = Y
    two_lifts[key] = (X_padic, Y_padic)
    two_diagnostics.append({"pivot_rows": rows, "determinant_mod59": determinant, "iterations": iterations})

two_points = []
for X, Y in two_finite:
    key = tuple(X.list())
    if key not in two_lifts:
        two_points.append(None)
        continue
    X_padic, Y_padic = two_lifts[key]
    sign = 1 if Y == two_representatives[key] else -1
    two_points.append(E_padic.point([KT(X_padic), KT(sign*Y_padic), KT(1)], check=False))

# Resolve only the component tests actually needed.  At t=0 reduction is
# nodal and finite-field group law is cheap: the desired one-node summand has
# the same component as C7.  At t=1 the desired two-node summand is the
# inverse component of C7, so Q+C7 is on the identity component.  The other
# summand must be identity at the opposite fibre.
KF = RF.fraction_field()
E_F = EllipticCurve(KF, [0, 0, 0, KF(A_F), KF(B_F)])
C7_F = E_F(KF(RF([reduce_qq(value, F) for value in c7_record["x_coefficients_low_to_high"]])),
           KF(RF([reduce_qq(value, F) for value in c7_record["y_coefficients_low_to_high"]])))


def finite_node_hit_finite(point, support, node):
    if point.is_zero():
        return False
    x_value, y_value = map(KF, point[:2])
    return (
        x_value.denominator()(support) != 0
        and y_value.denominator()(support) != 0
        and x_value(support) == node and y_value(support) == 0
    )


one_eligible = set()
for index, ((X_F, Y_F), point_padic) in enumerate(zip(one_finite, one_points)):
    point_F = E_F(KF(X_F), KF(Y_F))
    same_as_c7_at_zero = not finite_node_hit_finite(point_F-C7_F, F(0), node_zero_F)
    identity_at_one = not finite_node_hit_padic(point_padic, K(1), node_one)
    if same_as_c7_at_zero and identity_at_one:
        one_eligible.add(index)

two_eligible = set()
for index, ((X_F, Y_F), point_padic) in enumerate(zip(two_finite, two_points)):
    if point_padic is None:
        continue
    point_F = E_F(KF(X_F), KF(Y_F))
    identity_at_zero = not finite_node_hit_finite(point_F, F(0), node_zero_F)
    inverse_c7_at_one = not finite_node_hit_padic(point_padic+C7, K(1), node_one)
    if identity_at_zero and inverse_c7_at_one:
        two_eligible.add(index)


def reconstruct(value):
    digits = PRECISION-12
    modulus = PRIME**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


exact_lifts = []
attempted_pairs = 0
component_matched_pairs = 0
intersection_matched_pairs = 0
intersection_fingerprints = Counter()
for pair_index, source in enumerate(sums["horizontal_sum_seeds"]):
    one_index = int(source["one_node_index"])
    two_index = int(source["two_node_index"])
    one_point = one_points[one_index]
    two_point = two_points[two_index]
    if two_point is None:
        continue
    attempted_pairs += 1
    if one_index not in one_eligible or two_index not in two_eligible:
        continue
    point = one_point+two_point
    component_matched_pairs += 1
    inherited_intersections = tuple(
        section_intersection(point, inherited) for inherited in inherited_points
    )
    intersection_fingerprints[str(inherited_intersections)] += 1
    if inherited_intersections != expected_inherited_intersections:
        continue
    intersection_matched_pairs += 1
    x_value, y_value = map(KT, point[:2])
    x_denominator = x_value.denominator().monic()
    if x_denominator.degree() != 2:
        continue
    z = -x_denominator[1]/2
    Z = t-z
    if x_denominator != Z**2:
        continue
    X = RT(x_value*Z**2)
    Y = RT(y_value*Z**3)
    if X.degree() > 6 or Y.degree() > 9:
        continue
    padic_values = [z] + X.list() + [K.zero()]*(7-len(X.list()))
    padic_values += Y.list() + [K.zero()]*(10-len(Y.list()))
    try:
        exact = [reconstruct(value) for value in padic_values]
    except ArithmeticError:
        continue
    z_QQ = exact[0]
    Z_QQ = tq-z_QQ
    X_QQ = RQ(exact[1:8])
    Y_QQ = RQ(exact[8:18])
    if Y_QQ**2 != X_QQ**3+A_QQ*X_QQ*Z_QQ**4+B_QQ*Z_QQ**6:
        continue
    exact_lifts.append({
        "pair_index": pair_index,
        "source_one_node_index": int(source["one_node_index"]),
        "source_two_node_index": int(source["two_node_index"]),
        "z": str(z_QQ),
        "Z_coefficients_low_to_high": [str(value) for value in Z_QQ.list()],
        "X_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "maximum_rational_bits": max(map(coefficient_bits, exact)),
        "exact_compact_weierstrass_identity": True,
        "negative_section_also_certified": True,
    })

unique = {}
for record in exact_lifts:
    key = (record["z"], tuple(record["X_coefficients_low_to_high"]))
    unique.setdefault(key, record)
exact_lifts = list(unique.values())

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-via-summands-qq.v1",
    "status": (
        "PASS_EXACT_QQ_Q4O323_SIMPLE_POLE_HORIZONTAL_VIA_SUMMANDS"
        if exact_lifts else "PASS_NO_QQ_Q4O323_SUM_AFTER_RESOLVED_HENSEL"
    ),
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "resolved_hensel": {
        "one_node_signed_sections": len(one_points),
        "one_node_x_branches": len(one_lifts),
        "two_node_signed_sections": len(two_points),
        "two_node_regular_x_branches": len(two_lifts),
        "two_node_singular_x_branches": singular_two_x,
        "attempted_modular_sum_pairs": attempted_pairs,
        "eligible_one_node_signed_indices": sorted(one_eligible),
        "eligible_two_node_signed_indices": sorted(two_eligible),
        "c7_component_matched_pairs": component_matched_pairs,
        "inherited_intersection_matched_pairs": intersection_matched_pairs,
        "component_matched_inherited_intersection_fingerprints": dict(intersection_fingerprints),
        "exact_QQ_simple_pole_x_branches": len(exact_lifts),
        "one_node_diagnostics": one_diagnostics,
        "two_node_diagnostics": two_diagnostics,
    },
    "exact_QQ_horizontal_sections": exact_lifts,
    "method": {"large_Groebner_required": False, "runtime_seconds": time.monotonic()-started},
    "proof_boundary": (
        "Regular resolved one-node and two-node charts were Newton lifted separately and "
        "added over Q_59(t). Listed sums rationally reconstruct and satisfy the compact "
        "simple-pole Weierstrass identity exactly over QQ. Marked lattice matching and the "
        "q4/orbit323 genus-one pencil/Jacobian remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (COMPACT, ONE_NODE, SUMS, MARKING)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (COMPACT, ONE_NODE, SUMS, MARKING)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323SUMMANDLIFT|one={}|two_regular={}|pairs={}|component_matched={}|intersection_matched={}|qq={}|status={}|output={}".format(
        len(one_lifts), len(two_lifts), attempted_pairs, component_matched_pairs,
        intersection_matched_pairs, len(exact_lifts), payload["status"], OUTPUT,
    ), flush=True,
)
