#!/usr/bin/env sage -python
"""Lift and identify the q4/orbit164 parent-degree-three branch over QQ.

The mod-31 scan leaves two x-coordinate branches (with both y signs) after
imposing the infinity-I4 node and apparent parent degree three.  Fix the exact
infinity-node equations x_4=node and y_6=0, Newton lift each ten-variable
branch, rationally reconstruct, and retain exact QQ sections.  Exact group law
then identifies the rational lift inside the known rank-eight subgroup.  This
corrects the earlier, invalid identification by coarse node/parent-degree
fingerprints alone.

Only linear Newton systems are solved; no Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
QUARTIC = LOCAL / "q4o164-compact-binary-quartic-qq.json"
POINTING = LOCAL / "q4o164-c8-equation-marking-qq.json"
AUDIT = LOCAL / "q4o164-zero-one-node-parent-degree-audit.json"
DEGREE_ONE = LOCAL / "q4o1584-degree1-sections-to-q4o164-qq.json"
BASIS8 = LOCAL / "q4o164-integral-basis-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--modular", type=Path,
    default=LOCAL / "q4o164-missing-direction-parent-degree3-mod41.json",
)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o164-missing-direction-parent-degree3-qq.json",
)
args = parser.parse_args()
MODULAR = args.modular if args.modular.is_absolute() else ROOT / args.modular
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
INPUTS = (MODEL, QUARTIC, POINTING, MODULAR, AUDIT, DEGREE_ONE, BASIS8)
PRECISION = 180

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


model = json.loads(MODEL.read_text())
quartic_data = json.loads(QUARTIC.read_text())
pointing = json.loads(POINTING.read_text())
modular = json.loads(MODULAR.read_text())
audit = json.loads(AUDIT.read_text())
degree_one = json.loads(DEGREE_ONE.read_text())
basis8 = json.loads(BASIS8.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert quartic_data["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_BINARY_QUARTIC"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
assert modular["status"] == "PASS_MODP_Q4O164_MISSING_DIRECTION_PARENT_DEGREE_THREE_CANDIDATES"
assert audit["status"] == "PASS_EXACT_Q4O164_ZERO_ONE_NODE_PARENT_DEGREE_AUDIT"
assert degree_one["status"] == "PASS_EXACT_QQ_Q4O164_TWO_PRIMITIVE_ONE_NODE_SECTIONS"
assert basis8["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
PRIME = ZZ(modular["prime"])

RQ = PolynomialRing(QQ, "t")
tq = RQ.gen()
KQ = RQ.fraction_field()
A_QQ = RQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])

RX = PolynomialRing(QQ, "x")
xvar = RX.gen()
cubic_infinity = xvar**3 + A_QQ[8] * xvar + B_QQ[12]
repeated_infinity = cubic_infinity.gcd(cubic_infinity.derivative())
assert repeated_infinity.degree() == 1
node_infinity_QQ = QQ(-repeated_infinity[0] / repeated_infinity[1])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(PRIME)
RF = PolynomialRing(F, "t")
tf = RF.gen()
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
assert reduce_qq(node_infinity_QQ, F) == modular["search"]["infinity_node_x_leading_coefficient"]

K = Qp(PRIME, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
node_infinity = K(node_infinity_QQ)


def section_polynomials(values, ring, node):
    return ring(list(values[:4]) + [node]), ring(list(values[4:]) + [0])


def residual(values):
    X, Y = section_polynomials(values, RT, node_infinity)
    equation = Y**2 - X**3 - A * X - B
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A, variable):
    X, Y = section_polynomials(values, ring, variable)
    dx = -3 * X**2 - surface_A
    dy = 2 * Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), [[
        dx[degree-shift] if 0 <= degree-shift <= dx.degree() else zero
        for shift in range(4)
    ] + [
        dy[degree-shift] if 0 <= degree-shift <= dy.degree() else zero
        for shift in range(6)
    ] for degree in range(13)])


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def rational_reconstruct(value):
    modulus = PRIME ** (PRECISION - 12)
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


lifted = []
seen_x = set()
for modular_index, record in enumerate(modular["sections"]):
    x_seed = list(record["x_coefficients_low_to_high"])
    if tuple(x_seed) in seen_x:
        continue
    seen_x.add(tuple(x_seed))
    y_seed = list(record["y_coefficients_low_to_high"])
    assert len(x_seed) == 5 and len(y_seed) <= 6
    seed = vector(F, x_seed[:4] + y_seed + [0] * (6 - len(y_seed)))
    XF, YF = section_polynomials(seed, RF, reduce_qq(node_infinity_QQ, F))
    assert YF**2 == XF**3 + A_F * XF + B_F
    JF = jacobian(seed, RF, A_F, reduce_qq(node_infinity_QQ, F))
    rank = int(JF.rank())
    assert rank == 10
    pivot_rows = list(map(int, JF.transpose().pivots()))
    assert len(pivot_rows) == 10
    determinant = int(matrix(F, [JF.row(row) for row in pivot_rows]).det())
    assert determinant

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < PRECISION:
        working_precision = min(2 * known_precision, PRECISION)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        full = residual(values)
        chosen = vector(K, [full[row] for row in pivot_rows])
        J = jacobian(values, RT, A, node_infinity)
        square = matrix(K, [J.row(row) for row in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
            "minimum_correction_valuation": int(minimum_valuation(correction)),
        })
        known_precision = working_precision

    try:
        reconstructed = [rational_reconstruct(value) for value in values]
    except ArithmeticError:
        continue
    X, Y = section_polynomials(reconstructed, RQ, node_infinity_QQ)
    if Y**2 != X**3 + A_QQ * X + B_QQ:
        continue
    lifted.append({
        "modular_candidate_index": modular_index,
        "ordinary_resolved_jacobian_rank": rank,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_modp": determinant,
        "iterations": iterations,
        "x_coefficients_low_to_high": [str(value) for value in X.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "maximum_x_rational_bits": max(map(coefficient_bits, X)),
        "maximum_y_rational_bits": max(map(coefficient_bits, Y)),
        "exact_compact_weierstrass_identity": True,
    })

assert lifted, "no modular candidate rationally reconstructed to an exact QQ section"

# Exact inverse C8-pointed map and lattice fingerprints.
quartic = [RQ([QQ(value) for value in row]) for row in quartic_data["coefficients_in_T_low_to_high"]]
raw_W0_record = pointing["selected_zero"]["quartic_ordinate"]
RU = PolynomialRing(QQ, "U")
raw_W0_numerator = RU([QQ(value) for value in raw_W0_record["numerator_coefficients_low_to_high"]])
raw_W0_denominator = RU([QQ(value) for value in raw_W0_record["denominator_coefficients_low_to_high"]])
child_c = QQ(model["exact_coordinate_change"]["c"])
child_s = QQ(model["exact_coordinate_change"]["s"])
parent_c = QQ(quartic_data["exact_coordinate_change"]["parent_c"])
W0 = KQ(RQ(raw_W0_numerator(child_c * tq))) / (
    KQ(RQ(raw_W0_denominator(child_c * tq))) * child_s * parent_c
)
e, d, c, b, a = map(KQ, quartic)
assert e == W0**2
a1 = d / W0
a2 = c - d**2 / (4 * W0**2)
a3 = 2 * W0 * b
b2 = a1**2 + 4 * a2

supports = []
nodes = []
for fibre in model["compact_model"]["reducible_fibres"]:
    if fibre["support"] == "infinity":
        supports.append(None)
        cubic = cubic_infinity
    else:
        support = QQ(fibre["support"])
        supports.append(support)
        cubic = xvar**3 + A_QQ(support) * xvar + B_QQ(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    nodes.append(QQ(-repeated[0] / repeated[1]))

for record in lifted:
    X = RQ([QQ(value) for value in record["x_coefficients_low_to_high"]])
    Y = RQ([QQ(value) for value in record["y_coefficients_low_to_high"]])
    hits = []
    for support, node in zip(supports, nodes):
        hits.append(
            X[4] == node and Y[6] == 0
            if support is None else X(support) == node and Y(support) == 0
        )
    x_general = KQ(X) / 9 - b2 / 12
    y_general = KQ(Y) / 27 - (a1 * x_general + a3) / 2
    V = KQ(2 * W0 * (x_general + a2) / y_general)
    W = KQ((x_general * V**2 - d * V) / (2 * W0) - W0)
    assert W**2 == sum(KQ(value) * V**index for index, value in enumerate(quartic))
    record["exact_node_hits"] = hits
    record["parent_base_numerator_coefficients_low_to_high"] = [str(value) for value in V.numerator().list()]
    record["parent_base_denominator_coefficients_low_to_high"] = [str(value) for value in V.denominator().list()]
    record["q4o1584_parent_degree"] = int(max(V.numerator().degree(), V.denominator().degree()))

profile_matches = [
    record for record in lifted
    if record["exact_node_hits"] == [False, False, False, True]
    and record["q4o1584_parent_degree"] == 3
]

# Two classes with this coarse fingerprint were already transported exactly
# from q4/orbit1584.  Remove those literal points (up to sign); the remaining
# exact branch is the primitive ninth-coordinate class.
known_points = []
for row in degree_one["degree_one_sections"]:
    known_points.append((
        RQ([QQ(value) for value in row["compact_child_x"]["numerator_coefficients_low_to_high"]]),
        RQ([QQ(value) for value in row["compact_child_y"]["numerator_coefficients_low_to_high"]]),
    ))
for record in profile_matches:
    X = RQ([QQ(value) for value in record["x_coefficients_low_to_high"]])
    Y = RQ([QQ(value) for value in record["y_coefficients_low_to_high"]])
    matches = [index for index, (known_x, known_y) in enumerate(known_points) if X == known_x and Y in (known_y, -known_y)]
    record["known_transported_degree_one_matches_up_to_sign"] = matches
selected = [record for record in profile_matches if not record["known_transported_degree_one_matches_up_to_sign"]]
assert len(selected) == 1, (
    f"expected one new exact branch after removing transported sections; "
    f"lifted={len(lifted)}, profile_matches={len(profile_matches)}, selected={len(selected)}, "
    f"profiles={[(row['exact_node_hits'], row['q4o1584_parent_degree']) for row in lifted]}"
)

# Coarse node and inverse-parent-degree data do not attach an equation section
# to a marked component orientation.  Identify this branch by literal group
# law instead.  The exact relation also proves that it is not a ninth MW
# direction.
E_QQ = EllipticCurve(KQ, [0, 0, 0, KQ(A_QQ), KQ(B_QQ)])


def exact_point(record):
    return E_QQ(
        KQ(RQ([QQ(value) for value in record["x_coefficients_low_to_high"]])),
        KQ(RQ([QQ(value) for value in record["y_coefficients_low_to_high"]])),
    )


basis_points = [exact_point(record) for record in basis8["resolved_hensel"]["sections"]]
selected_point = exact_point(selected[0])
assert selected_point == 2 * basis_points[0] + basis_points[5] + basis_points[7]


def fourfold_height(point):
    """Clear every q4/o164 component group and read height from pole growth."""
    fourfold = 4 * point
    x_coordinate, y_coordinate = fourfold[0], fourfold[1]
    x_numerator, x_denominator = x_coordinate.numerator(), x_coordinate.denominator()
    y_numerator, y_denominator = y_coordinate.numerator(), y_coordinate.denominator()
    pole_degree = max(x_denominator.degree(), x_numerator.degree() - 4)
    assert pole_degree >= 0 and pole_degree % 2 == 0
    assert x_denominator**3 == y_denominator**2
    degrees = (
        x_numerator.degree(), x_denominator.degree(),
        y_numerator.degree(), y_denominator.degree(),
    )
    assert degrees == (
        4 + pole_degree, pole_degree,
        6 + 3 * pole_degree // 2, 3 * pole_degree // 2,
    )
    return {
        "fourfold_compact_degrees_x_num_x_den_y_num_y_den": list(map(int, degrees)),
        "cleared_pole_degree": int(pole_degree),
        "canonical_height": str(QQ(4 + pole_degree) / 16),
    }


selected_height_audit = fourfold_height(selected_point)
assert selected_height_audit["canonical_height"] == "13/4"

payload = {
    "schema": "elkies-k3.q4o164-parent-degree3-rank8-relation-qq.v3",
    "status": "PASS_EXACT_QQ_Q4O164_PARENT_DEGREE3_SECTION_RANK8_RELATION",
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "exact_rational_lifts": lifted,
    "selected_section": selected[0],
    "exact_group_law_identification": {
        "relation": "N=2*B0+B5+B7",
        "basis_order": "resolved_hensel.sections in q4o164-integral-basis-qq.json",
        "component_profile": {
            "finite_I2_labels": [0, 0],
            "finite_I4_label": 0,
            "infinity_I4_label": "1_or_3",
            "orientation_status": "conditional on the unresolved marked embedding",
        },
        "fourfold_pole_height_audit": selected_height_audit,
        "shioda_height": "13/4",
        "withdrawn_coarse_claim": (
            "The former [0,0,0,2] node-incidence profile and height 3 did not "
            "resolve the infinity-I4 tangent branch and are invalid."
        ),
        "outside_known_rank_eight_subgroup": False,
        "literal_QQ_function_field_identity": True,
    },
    "method": {
        "large_Groebner_required": False,
        "fixed_exact_node_equations": ["x_4=node_infinity", "y_6=0"],
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "The selected polynomial section is exact over QQ, but literal group law puts it in "
        "the known rank-eight subgroup. Fourfold pole growth gives height 13/4 and withdraws "
        "the former coarse infinity-I4 component-2 profile and height 3. A resolved local "
        "component anchor is still needed to choose its odd I4 orientation. Recovering the "
        "actual ninth direction and q8/orbit376 horizontal remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164MISSINGLIFT|exact_lifts={}|selected_mod={}|bits={}/{}|parent_degree={}|nodes={}|status={}|output={}".format(
        len(lifted), selected[0]["modular_candidate_index"],
        selected[0]["maximum_x_rational_bits"], selected[0]["maximum_y_rational_bits"],
        selected[0]["q4o1584_parent_degree"], selected[0]["exact_node_hits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
