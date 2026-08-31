#!/usr/bin/env sage -python
"""Lift q4/orbit323 via resolved polynomial summands and exact group law.

The compact source has bad auxiliary reduction at p=59, making the direct
simple-pole chart singular.  Its one-node and two-node polynomial charts are
resolved and mostly regular.  Lift those smaller charts separately, require
convergence of the full equation residual, reconstruct and verify each
polynomial summand over QQ, and only then add the prescribed pairs.  Compare
the resulting simple-pole sums with both signs of the resolved eleven-section
catalog.  No elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
ONE_NODE = LOCAL / "q4o208-q4o323-one-node-f1-mod59.json"
SUMS = LOCAL / "q4o208-q4o323-horizontal-sums-mod59.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
RESOLVED_HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
HORIZONTAL_MARKING = LOCAL / "q4o208-q4o323-horizontal-marking-qq.json"

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
resolved_horizontal = json.loads(RESOLVED_HORIZONTAL.read_text())
horizontal_marking = json.loads(HORIZONTAL_MARKING.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert one_node["status"] == "PASS_MODP_Q4O323_ONE_NODE_ADJACENT_COMPONENT_SECTION_SCAN"
assert sums["status"] == "PASS_MOD59_Q4O323_SIMPLE_POLE_SUM_SEEDS"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert resolved_horizontal["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert horizontal_marking["status"] == "PASS_EXACT_QQ_Q4O323_LIFTED_SHELL_EXCLUDES_TARGET"

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
    X_final, Y_final = unpack_padic(values)
    final_residual_valuation = int(minimum_valuation(
        equation_residual(X_final, Y_final, A, B, K)
    ))
    return (
        values, pivot_rows, determinant, iterations,
        final_residual_valuation >= PRECISION-12,
    )


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
one_full_residual_rejected = 0
one_seen = set()
for X, Y in one_finite:
    key = tuple(X.list())
    if key in one_seen:
        continue
    one_seen.add(key)
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
    values, rows, determinant, iterations, full_residual_converged = result
    one_diagnostics.append({
        "pivot_rows": rows,
        "determinant_mod59": determinant,
        "iterations": iterations,
        "full_residual_converged": full_residual_converged,
    })
    if not full_residual_converged:
        one_full_residual_rejected += 1
        continue
    X_padic, Y_padic = one_unpack(values, RT, node_zero)
    one_representatives[key] = Y
    one_lifts[key] = (X_padic, Y_padic)

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
    if key not in one_lifts:
        one_points.append(None)
        continue
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
two_full_residual_rejected = 0
two_seen = set()
for X, Y in two_finite:
    key = tuple(X.list())
    if key in two_seen:
        continue
    two_seen.add(key)
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
    values, rows, determinant, iterations, full_residual_converged = result
    two_diagnostics.append({
        "pivot_rows": rows,
        "determinant_mod59": determinant,
        "iterations": iterations,
        "full_residual_converged": full_residual_converged,
    })
    if not full_residual_converged:
        two_full_residual_rejected += 1
        continue
    X_padic, Y_padic = two_unpack(values, RT, node_one, node_infinity)
    two_representatives[key] = Y
    two_lifts[key] = (X_padic, Y_padic)

two_points = []
for X, Y in two_finite:
    key = tuple(X.list())
    if key not in two_lifts:
        two_points.append(None)
        continue
    X_padic, Y_padic = two_lifts[key]
    sign = 1 if Y == two_representatives[key] else -1
    two_points.append(E_padic.point([KT(X_padic), KT(sign*Y_padic), KT(1)], check=False))

def reconstruct(value):
    digits = PRECISION-12
    modulus = PRIME**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


KQ = RQ.fraction_field()
E_QQ = EllipticCurve(KQ, [0, 0, 0, KQ(A_QQ), KQ(B_QQ)])


def reconstruct_polynomial(value, coefficient_count):
    coefficients = list(value.list())
    coefficients += [K.zero()]*(coefficient_count-len(coefficients))
    return RQ([reconstruct(coefficient) for coefficient in coefficients])


def reconstruct_polynomial_branches(lifts, x_count, y_count):
    exact = {}
    for key, (X_padic, Y_padic) in lifts.items():
        try:
            X_QQ = reconstruct_polynomial(X_padic, x_count)
            Y_QQ = reconstruct_polynomial(Y_padic, y_count)
        except ArithmeticError:
            continue
        if Y_QQ**2 != X_QQ**3+A_QQ*X_QQ+B_QQ:
            continue
        exact[key] = (X_QQ, Y_QQ)
    return exact


one_exact_branches = reconstruct_polynomial_branches(one_lifts, 5, 7)
two_exact_branches = reconstruct_polynomial_branches(two_lifts, 5, 6)


def signed_exact_points(finite_points, representatives, exact_branches):
    points = []
    polynomial_coordinates = []
    for X_finite, Y_finite in finite_points:
        key = tuple(X_finite.list())
        if key not in exact_branches:
            points.append(None)
            polynomial_coordinates.append(None)
            continue
        X_QQ, representative_Y_QQ = exact_branches[key]
        sign = 1 if Y_finite == representatives[key] else -1
        Y_QQ = sign*representative_Y_QQ
        points.append(E_QQ.point([KQ(X_QQ), KQ(Y_QQ), KQ(1)], check=False))
        polynomial_coordinates.append((X_QQ, Y_QQ))
    return points, polynomial_coordinates


one_points_QQ, one_polynomials_QQ = signed_exact_points(
    one_finite, one_representatives, one_exact_branches
)
two_points_QQ, two_polynomials_QQ = signed_exact_points(
    two_finite, two_representatives, two_exact_branches
)


def transform_exact_old_coordinate_qq(record, weight):
    change = compact["exact_coordinate_change"]
    a, b, c, d, m = [QQ(change[key]) for key in ("a", "b", "c", "d", "m")]
    numerator_linear = a*tq+b
    denominator_linear = c*tq+d
    numerator = [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    denominator = [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    numerator_value = sum(
        value*numerator_linear**degree*denominator_linear**(weight-degree)
        for degree, value in enumerate(numerator)
    )
    denominator_degree = len(denominator)-1
    denominator_value = sum(
        value*numerator_linear**degree*denominator_linear**(denominator_degree-degree)
        for degree, value in enumerate(denominator)
    )
    return KQ(numerator_value)/(m**(weight//2)*KQ(denominator_value))


inherited_points_QQ = []
for label in inherited_labels:
    record = marking[label]
    inherited_points_QQ.append(E_QQ(
        transform_exact_old_coordinate_qq(record["x"], 4),
        transform_exact_old_coordinate_qq(record["y"], 6),
    ))


def exact_section_intersection(left, right):
    difference = left-right
    if difference.is_zero():
        return -2
    x_value = KQ(difference[0])
    numerator_degree = int(x_value.numerator().degree())
    denominator_degree = int(x_value.denominator().degree())
    assert denominator_degree % 2 == 0
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    assert infinity_excess % 2 == 0
    return denominator_degree//2 + infinity_excess//2


C7_QQ = E_QQ.point([
    KQ(RQ([QQ(value) for value in c7_record["x_coefficients_low_to_high"]])),
    KQ(RQ([QQ(value) for value in c7_record["y_coefficients_low_to_high"]])),
    KQ(1),
], check=False)


def exact_finite_node_hit(point, support, node):
    if point.is_zero():
        return False
    x_value, y_value = map(KQ, point[:2])
    return (
        x_value.denominator()(support) != 0
        and y_value.denominator()(support) != 0
        and x_value(support) == node and y_value(support) == 0
    )


# These component filters are now evaluated only on literal QQ sections.
one_eligible = set()
for index, point in enumerate(one_points_QQ):
    if point is None:
        continue
    same_as_c7_at_zero = not exact_finite_node_hit(point-C7_QQ, QQ(0), node_zero_QQ)
    identity_at_one = not exact_finite_node_hit(point, QQ(1), node_one_QQ)
    if same_as_c7_at_zero and identity_at_one:
        one_eligible.add(index)

two_eligible = set()
for index, point in enumerate(two_points_QQ):
    if point is None:
        continue
    identity_at_zero = not exact_finite_node_hit(point, QQ(0), node_zero_QQ)
    inverse_c7_at_one = not exact_finite_node_hit(point+C7_QQ, QQ(1), node_one_QQ)
    if identity_at_zero and inverse_c7_at_one:
        two_eligible.add(index)


def polynomial_record(pair):
    X_QQ, Y_QQ = pair
    return {
        "X_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "exact_compact_weierstrass_identity": True,
    }


def resolved_catalog_matches(X_QQ, Y_QQ, z_QQ):
    matches = []
    for record in resolved_horizontal["exact_QQ_horizontal_sections"]:
        if QQ(record["z"]) != z_QQ:
            continue
        catalog_X = RQ([QQ(value) for value in record["X_coefficients_low_to_high"]])
        catalog_Y = RQ([QQ(value) for value in record["Y_coefficients_low_to_high"]])
        if catalog_X != X_QQ:
            continue
        if catalog_Y == Y_QQ:
            matches.append({"branch_index": int(record["branch_index"]), "Y_sign": 1})
        if catalog_Y == -Y_QQ:
            matches.append({"branch_index": int(record["branch_index"]), "Y_sign": -1})
    return matches


exact_lifts = []
all_simple_pole_sums = []
attempted_pairs = 0
component_matched_pairs = 0
intersection_matched_pairs = 0
intersection_fingerprints = Counter()
for pair_index, source in enumerate(sums["horizontal_sum_seeds"]):
    one_index = int(source["one_node_index"])
    two_index = int(source["two_node_index"])
    one_point = one_points_QQ[one_index]
    two_point = two_points_QQ[two_index]
    if one_point is None or two_point is None:
        continue
    attempted_pairs += 1
    if one_index not in one_eligible or two_index not in two_eligible:
        continue
    point = one_point+two_point
    component_matched_pairs += 1
    inherited_intersections = tuple(
        exact_section_intersection(point, inherited)
        for inherited in inherited_points_QQ
    )
    intersection_fingerprints[str(inherited_intersections)] += 1
    x_value, y_value = map(KQ, point[:2])
    x_denominator = x_value.denominator().monic()
    if x_denominator.degree() != 2:
        continue
    z = -x_denominator[1]/2
    Z = tq-z
    if x_denominator != Z**2:
        continue
    X = RQ(x_value*Z**2)
    Y = RQ(y_value*Z**3)
    if X.degree() > 6 or Y.degree() > 9:
        continue
    z_QQ = QQ(z)
    Z_QQ = RQ(Z)
    X_QQ = RQ(X)
    Y_QQ = RQ(Y)
    if Y_QQ**2 != X_QQ**3+A_QQ*X_QQ*Z_QQ**4+B_QQ*Z_QQ**6:
        continue
    catalog_matches = resolved_catalog_matches(X_QQ, Y_QQ, z_QQ)
    record = {
        "pair_index": pair_index,
        "source_one_node_index": one_index,
        "source_two_node_index": two_index,
        "z": str(z_QQ),
        "Z_coefficients_low_to_high": [str(value) for value in Z_QQ.list()],
        "X_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "summands": {
            "one_node": polynomial_record(one_polynomials_QQ[one_index]),
            "two_node": polynomial_record(two_polynomials_QQ[two_index]),
        },
        "inherited_intersections": list(inherited_intersections),
        "resolved_11_section_catalog_matches_both_Y_signs": catalog_matches,
        "novel_relative_to_resolved_11_section_catalog": not catalog_matches,
        "maximum_rational_bits": max(
            [coefficient_bits(z_QQ)]
            + [coefficient_bits(value) for value in X_QQ]
            + [coefficient_bits(value) for value in Y_QQ]
        ),
        "exact_compact_weierstrass_identity": True,
        "negative_section_also_certified": True,
    }
    all_simple_pole_sums.append(record)
    if inherited_intersections != expected_inherited_intersections:
        continue
    intersection_matched_pairs += 1
    exact_lifts.append(record)

unique = {}
for record in exact_lifts:
    key = (record["z"], tuple(record["X_coefficients_low_to_high"]))
    unique.setdefault(key, record)
exact_lifts = list(unique.values())


def exact_polynomial_records(polynomials):
    return [
        {
            "signed_index": index,
            "X_coefficients_low_to_high": [str(value) for value in X.list()],
            "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        }
        for index, value in enumerate(polynomials) if value is not None
        for X, Y in (value,)
    ]


payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-via-summands-qq.v2",
    "status": (
        "PASS_EXACT_QQ_Q4O323_SUM_REPRODUCES_BRANCH33_INVERSE_TARGET_EXCLUDED"
        if exact_lifts
        else "PASS_NO_QQ_Q4O323_SUM_AFTER_EXACT_SUMMAND_RECONSTRUCTION"
    ),
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "resolved_hensel": {
        "one_node_mod59_signed_sections": len(one_points),
        "one_node_full_residual_x_branches": len(one_lifts),
        "one_node_full_residual_rejected_x_branches": one_full_residual_rejected,
        "one_node_exact_QQ_x_branches": len(one_exact_branches),
        "one_node_exact_QQ_signed_indices": [
            index for index, point in enumerate(one_points_QQ) if point is not None
        ],
        "two_node_mod59_signed_sections": len(two_points),
        "two_node_full_residual_x_branches": len(two_lifts),
        "two_node_full_residual_rejected_x_branches": two_full_residual_rejected,
        "two_node_exact_QQ_x_branches": len(two_exact_branches),
        "two_node_exact_QQ_signed_indices": [
            index for index, point in enumerate(two_points_QQ) if point is not None
        ],
        "two_node_singular_x_branches": singular_two_x,
        "stored_seed_pairs_with_two_exact_QQ_summands": attempted_pairs,
        "eligible_exact_QQ_one_node_signed_indices": sorted(one_eligible),
        "eligible_exact_QQ_two_node_signed_indices": sorted(two_eligible),
        "exact_QQ_component_matched_pairs": component_matched_pairs,
        "inherited_intersection_matched_pairs": intersection_matched_pairs,
        "exact_QQ_component_matched_inherited_intersection_fingerprints": dict(intersection_fingerprints),
        "all_exact_QQ_simple_pole_sum_count": len(all_simple_pole_sums),
        "profile_matched_exact_QQ_simple_pole_sum_count": len(exact_lifts),
        "one_node_diagnostics": one_diagnostics,
        "two_node_diagnostics": two_diagnostics,
    },
    "exact_QQ_polynomial_sections": {
        "one_node_signed": exact_polynomial_records(one_polynomials_QQ),
        "two_node_signed": exact_polynomial_records(two_polynomials_QQ),
    },
    "all_exact_QQ_simple_pole_sums": all_simple_pole_sums,
    "exact_QQ_horizontal_sections": exact_lifts,
    "resolved_11_section_catalog_comparison": {
        "both_Y_signs_tested": True,
        "catalog_section_count": len(resolved_horizontal["exact_QQ_horizontal_sections"]),
        "profile_matched_sum_count": len(exact_lifts),
        "profile_matched_sums_novel_to_catalog": all(
            record["novel_relative_to_resolved_11_section_catalog"] for record in exact_lifts
        ),
        "branch33_target_hits_in_complete_graph_solutions": int(
            horizontal_marking["lattice_match"][
                "branch33_target_hits_in_complete_graph_solutions"
            ]
        ),
        "full_graph_excludes_q4o323_target": (
            horizontal_marking["lattice_match"][
                "branch33_target_hits_in_complete_graph_solutions"
            ] == 0
        ),
    },
    "method": {"large_Groebner_required": False, "runtime_seconds": time.monotonic()-started},
    "proof_boundary": (
        "Only branches whose full thirteen-equation residual converges are retained. Each "
        "polynomial summand is rationally reconstructed and verified over QQ before exact "
        "group addition. Listed sums satisfy the compact simple-pole Weierstrass identity "
        "exactly. The unique (3,2,2) sum is exactly branch 33 with the opposite Y sign, not a "
        "new section. The complete signed eleven-section graph has zero q4/orbit323 target "
        "hits, so this reconstruction does not recover the marked target. A different "
        "presentation or the missing Mordell--Weil direction remains necessary; no q4/orbit323 "
        "pencil or child Jacobian is certified."
    ),
    "inputs": {
        "paths": [
            str(path.relative_to(ROOT))
            for path in (
                COMPACT, ONE_NODE, SUMS, MARKING, RESOLVED_HORIZONTAL,
                HORIZONTAL_MARKING,
            )
        ],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                COMPACT, ONE_NODE, SUMS, MARKING, RESOLVED_HORIZONTAL,
                HORIZONTAL_MARKING,
            )
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323SUMMANDLIFT|one={}|two_regular={}|pairs={}|component_matched={}|intersection_matched={}|qq={}|status={}|output={}".format(
        len(one_lifts), len(two_lifts), attempted_pairs, component_matched_pairs,
        intersection_matched_pairs, len(exact_lifts), payload["status"], OUTPUT,
    ), flush=True,
)
