#!/usr/bin/env sage -python
"""Compile and mark the exact resolved RR model for q8/orbit376.

status: PROVED_EXACT
claim: exact QQ 12-to-4-to-2 RR plane, 4A1 Jacobian, and P1229 zero
outputs: artifacts/local/elkies-k3/q4o164-q8o376-smooth-rr-qq.json

The exact inherited-P1/Abel reconstruction supplies a section

    x_H=X/Z^2, y_H=Y/Z^3,  deg(X,Y,Z)=(12,18,4).

For the q8 marked chord, the global smooth ambient is

    a=AA/Z^2, deg AA<=8;  b=BB/Z, deg BB<=2.

Regularity at the four smooth collisions is the polynomial congruence
``AA*X=BB*Y mod Z^2``.  Inverting X modulo Z^2 solves it directly and leaves
the exact four-plane.  Split-toric quotient rows on the two old I4 chains
leave h0=2.  The selected marked embedding orients the T=0 components, and an
exact rational arc on the first finite-I2 exceptional conic points the child
at P1229.  No Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import (
    Conic, FunctionField, GF, LaurentSeriesRing, PolynomialRing, QQ, ZZ,
    factorial, inverse_mod, matrix,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
HEIGHT_AUDIT = LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json"
SELECTOR = LOCAL / "q4o164-q8o376-horizontal-from-abel-trace-mod131.json"
SOURCE_MARKING = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
EDGE = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-certificate.json"
COST = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-old_a11_component_8-q8d2-cap10000-mw13-equation-cost.json"
SOURCE_FRAME = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-frame.txt"
OUTPUT = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
INPUTS = (
    MODEL, HORIZONTAL, BASIS, HEIGHT_AUDIT, SELECTOR, SOURCE_MARKING,
    EDGE, COST, SOURCE_FRAME,
)

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join("{}={}".format(key, value) for key, value in fields.items())
    print(
        "Q8O376SMOOTHQQ|stage={}|elapsed={:.3f}{}".format(
            stage, time.monotonic()-started, "|"+suffix if suffix else ""
        ),
        flush=True,
    )


def rational_bits(values):
    values = [QQ(value) for value in values]
    return int(max(
        max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
        for value in values
    ))


def coefficients(poly):
    return [str(value) for value in poly.list()]


model = json.loads(MODEL.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
basis = json.loads(BASIS.read_text())
height_audit = json.loads(HEIGHT_AUDIT.read_text())
selector = json.loads(SELECTOR.read_text())
source_marking = json.loads(SOURCE_MARKING.read_text())
edge = json.loads(EDGE.read_text())
cost = json.loads(COST.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert horizontal["selected_identity"] == "H=T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7"
assert horizontal["checks"]["exact_QQ_weierstrass_identity"]
assert horizontal["checks"]["canonical_height"] == "11"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert height_audit["status"] == "PASS_EXACT_QQ_Q4O164_FOURFOLD_HEIGHT_GRAM_AND_C8_MARKED_EMBEDDING_CENSUS"
assert selector["status"] == "PASS_EXACT_MODP_Q4O164_Q8O376_HORIZONTAL_FROM_ABEL_TRACE"
assert selector["prime"] == 131 and selector["selected"]["embedding_index"] == 15
assert source_marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert source_marking["zero"] == "old_A11_component_8"
assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert edge["candidate_id"] == {"label": "q8o376", "q": 8, "old_fibre_degree": 2}
candidate = next(
    item for item in cost["retained_candidates"]
    if item["candidate_id"] == {"q": 8, "old_fibre_degree": 2, "orbit_index": 376}
)
assert candidate["fibre"] == edge["source_to_child_basis"][0]
assert candidate["horizontal"]["P_dot_O"] == 4
assert candidate["horizontal"]["vertical"] == [1, 0, 0, 1, 1, 1, 1, 1]
source_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in SOURCE_FRAME.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
source_cartan = source_frame[:8, :8]
first_a3 = (0, 3, 4)
second_a3 = (5, 6, 7)
A3 = matrix(ZZ, ((2, -1, 0), (-1, 2, -1), (0, -1, 2)))
assert source_cartan.matrix_from_rows_and_columns(first_a3, first_a3) == A3
assert source_cartan.matrix_from_rows_and_columns(second_a3, second_a3) == A3
assert not any(source_cartan.matrix_from_rows_and_columns(first_a3, second_a3).list())

# The exact mod-131 degree fingerprint selects marked embedding 15.  Its
# fibre symmetry swaps the two I4 blocks into equation order but does not
# reverse the marked second block.  Since equation profiles order the finite
# I4 before the I4 at infinity, the compact T=0 chain is therefore the
# ordered marked block (5,6,7).
embedding15 = next(
    item for item in height_audit["marked_embedding_enumeration"]["embeddings"]
    if item["embedding_index"] == 15
)
assert embedding15["compatible_fibre_symmetries_swapI2_swapI4_reverseI4a_reverseI4b"] == [
    [True, True, False, False],
    [True, True, True, False],
]
finite_i4_root_indices = second_a3
finite_i4_component_labels = (
    "old_A11_component_6",
    "second_old_I6_I4_missing_component",
    "old_A11_component_5",
)

# The two marked A1 roots are indices 1 and 2.  Their affine components are
# old_zero=F+e_1 and P1229=F+e_2.  Embedding 15 swaps the two I2 factors into
# equation order, so P1229 is the nonidentity component of the first compact
# finite I2.
explicit_source_curves = source_marking["equation_explicit_curves_in_child"]
p1229_source = explicit_source_curves["P1229"]
old_zero_source = explicit_source_curves["old_zero"]
assert p1229_source[:5] == [1, 0, 0, 0, 1] and not any(p1229_source[5:])
assert old_zero_source[:5] == [1, 0, 0, 1, 0] and not any(old_zero_source[5:])
compact_reducible_fibres = model["compact_model"]["reducible_fibres"]
assert [item["kodaira"] for item in compact_reducible_fibres] == ["I2", "I2", "I4", "I4"]
assert compact_reducible_fibres[2]["support"] == "0"
assert compact_reducible_fibres[3]["support"] == "infinity"
p1229_old_base_support = QQ(compact_reducible_fibres[0]["support"])
assert p1229_old_base_support == QQ(25281)/QQ(168246841)

R = PolynomialRing(QQ, "T")
T = R.gen()
compact = model["compact_model"]
A = R([QQ(value) for value in compact["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in compact["B_coefficients_low_to_high"]])


def rational_function(record):
    numerator = R([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = R([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return numerator, denominator


X, denominator_x = rational_function(horizontal["section"]["x"])
Y, denominator_y = rational_function(horizontal["section"]["y"])
assert (X.degree(), denominator_x.degree(), Y.degree(), denominator_y.degree()) == (12, 8, 18, 12)
assert denominator_x.is_monic() and denominator_y.is_monic()

# Recover the common monic denominator Z without a general polynomial solve.
factorization = list(denominator_x.factor())
assert all(int(exponent) % 2 == 0 for unused, exponent in factorization)
Z = R.one()
for factor, exponent in factorization:
    Z *= factor.monic() ** (int(exponent)//2)
assert Z.degree() == 4 and Z**2 == denominator_x and Z**3 == denominator_y
assert Y**2 == X**3 + A*X*Z**4 + B*Z**6
log("LOAD", degrees="12,18,4", coefficient_bits=rational_bits(X.list()+Y.list()+Z.list()))

# Exact 12 -> 4 smooth collision module.
aa_degree = 8
bb_degree = 2
collision_modulus = Z**2
assert collision_modulus.degree() == 8
X_inverse = R(inverse_mod(X, collision_modulus))
assert (X*X_inverse) % collision_modulus == 1

pairs = [(Z**2, R.zero())]
for degree in range(bb_degree+1):
    BB = T**degree
    AA = R((BB*Y*X_inverse) % collision_modulus)
    assert AA.degree() <= 7
    assert (AA*X-BB*Y) % collision_modulus == 0
    pairs.append((AA, BB))

ambient_rows = matrix(QQ, [
    [AA[index] for index in range(aa_degree+1)]
    + [BB[index] for index in range(bb_degree+1)]
    for AA, BB in pairs
])
assert ambient_rows.rank() == 4

# Independently form all eight scalar rows and compare the complete kernel.
ambient = [
    (T**degree, R.zero()) for degree in range(aa_degree+1)
] + [
    (R.zero(), T**degree) for degree in range(bb_degree+1)
]
remainders = [R((AA*X-BB*Y) % collision_modulus) for AA, BB in ambient]
condition_matrix = matrix(QQ, [
    [remainder[degree] for remainder in remainders]
    for degree in range(collision_modulus.degree())
])
assert condition_matrix.nrows() == 8 and condition_matrix.ncols() == 12
assert condition_matrix.rank() == 8
kernel = condition_matrix.right_kernel().basis_matrix()
assert kernel.nrows() == 4
assert ambient_rows.row_space() == kernel.row_space()
log("SMOOTH_RR", ambient=12, rank=8, kernel=4, status="PASS_EXACT")


# Resolve the two old I4 fibres.  The lattice vertical divisor contains the
# complete nonidentity A3 chain at each one, so each connected chain supplies
# one quotient line.  We compute the traces in the standard split-I_n toric
# chart, entirely over QQ.
Kr = FunctionField(QQ, "r")
r = Kr.gen()
PREC = 12
LS = LaurentSeriesRing(Kr, "s", default_prec=PREC)
s = LS.gen()


def weighted_local(poly, weight, infinity):
    poly = R(poly)
    if not infinity:
        return sum((LS(Kr(value))*s**degree for degree, value in enumerate(poly.list())), LS.zero())
    return sum(
        (LS(Kr(value))*s**(weight-degree) for degree, value in enumerate(poly.list())),
        LS.zero(),
    )


def finite_support_local(poly, support):
    translated = R(poly)(T+QQ(support))
    return sum(
        (LS(Kr(value))*s**degree for degree, value in enumerate(translated.list())),
        LS.zero(),
    )


def newton_sqrt(value, constant_root):
    answer = LS(Kr(constant_root))
    for unused in range(5):
        answer = (answer+value/answer)/2
    assert (answer**2-value).valuation() >= PREC-4
    return answer


# B0 has exact equation component profile (0,1,1,0), so it meets component 1
# of the finite I4.  Its local branch fixes which of the two split-toric rho
# signs agrees with the marked chain orientation.
b0_record = basis["resolved_hensel"]["sections"][0]
assert b0_record["component_profile"] == [0, 1, 1, 0]
B0X = R([QQ(value) for value in b0_record["x_coefficients_low_to_high"]])
B0Y = R([QQ(value) for value in b0_record["y_coefficients_low_to_high"]])
assert B0Y**2 == B0X**3 + A*B0X + B


def finite_i4_rho_orientation():
    Aloc = weighted_local(A, 8, False)
    Bloc = weighted_local(B, 12, False)
    Xsec = weighted_local(B0X, 4, False)
    Ysec = weighted_local(B0Y, 6, False)
    A0, B0 = QQ(Aloc[0]), QQ(Bloc[0])
    node = QQ(-3*B0/(2*A0))
    center = newton_sqrt(-Aloc/3, node)
    nodal_error = center**3+Aloc*center+Bloc
    rho0 = QQ((3*node).sqrt())
    trials = []
    for rho_start in (rho0, -rho0):
        rho = newton_sqrt(Xsec+2*center, rho_start)
        ww = rho*(Xsec-center)
        assert (rho**3-3*center*rho-ww).valuation() >= PREC-4
        assert (Ysec**2-ww**2-nodal_error).valuation() >= PREC-4
        left_valuation = int((Ysec-ww).valuation())
        right_valuation = int((Ysec+ww).valuation())
        trials.append((rho_start, left_valuation, right_valuation))
    selected = [
        rho_start for rho_start, left_valuation, right_valuation in trials
        if (left_valuation, right_valuation) == (1, 3)
    ]
    assert len(selected) == 1
    return selected[0], trials


finite_rho_component_one, finite_rho_trials = finite_i4_rho_orientation()


def functional_rows(values):
    nonzero = [value for value in values if value]
    common = Kr.one()
    for value in nonzero:
        common = common.lcm(value.denominator())
    numerators = [(value*common).numerator() for value in values]
    degree = max((value.degree() for value in numerators if value), default=-1)
    rows = []
    for coefficient_degree in range(degree+1):
        row = tuple(
            QQ(value[coefficient_degree]) if value and coefficient_degree <= value.degree() else QQ(0)
            for value in numerators
        )
        if any(row):
            pivot = next(value for value in row if value)
            normalized = tuple(value/pivot for value in row)
            if normalized not in rows:
                rows.append(normalized)
    return rows


def i4_trace_rows(infinity):
    Aloc = weighted_local(A, 8, infinity)
    Bloc = weighted_local(B, 12, infinity)
    Zloc = weighted_local(Z, 4, infinity)
    Xloc = weighted_local(X, 12, infinity)
    Yloc = weighted_local(Y, 18, infinity)
    A0, B0 = QQ(Aloc[0]), QQ(Bloc[0])
    node = QQ(-3*B0/(2*A0))
    assert node**3+A0*node+B0 == 0 and 3*node**2+A0 == 0
    center = newton_sqrt(-Aloc/3, node)
    nodal_error = center**3+Aloc*center+Bloc
    assert int(nodal_error.valuation()) == 4
    unit = nodal_error/s**4
    rho_square = QQ(3*node)
    assert rho_square.is_square()
    rho0 = QQ(rho_square.sqrt())
    Hxloc = Xloc/Zloc**2
    Hyloc = Yloc/Zloc**3
    all_trials = []
    for rho_start in (rho0, -rho0):
        trace_rows = []
        diagnostics = []
        for component in (1, 2, 3):
            left = LS(r)*s**component
            right = unit*s**(4-component)/LS(r)
            yy = (left+right)/2
            ww = (right-left)/2
            rho = LS(Kr(rho_start))
            for unused in range(5):
                rho -= (rho**3-3*center*rho-ww)/(3*rho**2-3*center)
            assert (rho**3-3*center*rho-ww).valuation() >= PREC-5
            xx = center+ww/rho
            mloc = (yy+Hyloc)/(xx-Hxloc)
            valuation = int(mloc.valuation())
            if valuation < 0:
                diagnostics.append((component, valuation, "pole"))
                continue
            mres = Kr(mloc[0]) if valuation == 0 else Kr.zero()
            values = []
            for AA, BB in pairs:
                AAloc = weighted_local(AA, 8, infinity)
                BBloc = weighted_local(BB, 2, infinity)
                values.append(Kr(AAloc[0]) + Kr(BBloc[0])*Kr(Zloc[0])*mres)
            local_rows = functional_rows(values)
            trace_rows.extend(local_rows)
            diagnostics.append((component, valuation, str(mres), local_rows))
        unique_rows = []
        for row in trace_rows:
            if row not in unique_rows:
                unique_rows.append(row)
        all_trials.append((rho_start, unique_rows, diagnostics))
    return all_trials


def i4_component_local_series(infinity, rho_start, component):
    """Return the exact split-toric generic point and chord slope."""
    Aloc = weighted_local(A, 8, infinity)
    Bloc = weighted_local(B, 12, infinity)
    Zloc = weighted_local(Z, 4, infinity)
    Xloc = weighted_local(X, 12, infinity)
    Yloc = weighted_local(Y, 18, infinity)
    A0, B0 = QQ(Aloc[0]), QQ(Bloc[0])
    node = QQ(-3*B0/(2*A0))
    center = newton_sqrt(-Aloc/3, node)
    nodal_error = center**3+Aloc*center+Bloc
    unit = nodal_error/s**4
    left = LS(r)*s**component
    right = unit*s**(4-component)/LS(r)
    yy = (left+right)/2
    ww = (right-left)/2
    rho = LS(Kr(rho_start))
    for unused in range(5):
        rho -= (rho**3-3*center*rho-ww)/(3*rho**2-3*center)
    assert (rho**3-3*center*rho-ww).valuation() >= PREC-5
    xx = center+ww/rho
    mloc = (yy+Yloc/Zloc**3)/(xx-Xloc/Zloc**2)
    return {
        "Aloc": Aloc, "Bloc": Bloc, "Zloc": Zloc,
        "Hxloc": Xloc/Zloc**2, "Hyloc": Yloc/Zloc**3,
        "x": xx, "y": yy, "m": mloc,
    }


vertical_rows = []
vertical_diagnostics = []
for fibre_label, infinity in (("zero", False), ("infinity", True)):
    trials = i4_trace_rows(infinity)
    trial_row_sets = {tuple(rows) for unused, rows, unused_diag in trials}
    assert len(trial_row_sets) == 1
    selected_rows = list(trial_row_sets.pop())
    assert len(selected_rows) == 1
    vertical_rows.append(selected_rows[0])
    for rho_start, trace_rows, diagnostics in trials:
        print(
            "Q8O376I4TRACE|fibre={}|rho={}|rows={}|diag={}".format(
                fibre_label, rho_start, trace_rows, diagnostics,
            ),
            flush=True,
        )
    vertical_diagnostics.append({
        "fibre": fibre_label,
        "orientation_independent": True,
        "row_on_smooth_basis": [str(value) for value in selected_rows[0]],
    })

vertical_matrix = matrix(QQ, vertical_rows)
assert vertical_matrix.rank() == 2
resolved_kernel = vertical_matrix.right_kernel().basis_matrix()
assert resolved_kernel.nrows() == resolved_kernel.rank() == 2
assert vertical_matrix*resolved_kernel.transpose() == matrix(QQ, 2, 2)
resolved_pairs = []
for row in resolved_kernel.rows():
    AA = sum((row[index]*pairs[index][0] for index in range(4)), R.zero())
    BB = sum((row[index]*pairs[index][1] for index in range(4)), R.zero())
    resolved_pairs.append((AA, BB))
AA0, BB0 = resolved_pairs[0]
AA1, BB1 = resolved_pairs[1]
log("RESOLVED_RR", smooth=4, vertical_rank=2, h0=2, status="PASS_EXACT")
print(
    "Q8O376RESOLVEDB|BB0={}|BB1={}".format(BB0, BB1),
    flush=True,
)

# Compile the chord radicand over QQ(u), then inspect its univariate square
# factors.  This is a degree-at-most-16 gcd/factor calculation, not a surface
# elimination.
UQ = PolynomialRing(QQ, "u")
u = UQ.gen()
KU = UQ.fraction_field()
TU = PolynomialRing(KU, "T")
tt = TU.gen()


def lift_to_TU(poly):
    return TU([KU(value) for value in R(poly).list()])


aa = lift_to_TU(AA0) + KU(u)*lift_to_TU(AA1)
bb = lift_to_TU(BB0) + KU(u)*lift_to_TU(BB1)
X_u, Y_u, Z_u, A_u = map(lift_to_TU, (X, Y, Z, A))
raw = (
    aa**4 - 6*X_u*aa**2*bb**2 + 8*Y_u*aa*bb**3
    - 3*X_u**2*bb**4 - 4*A_u*bb**4*Z_u**4
)
after_collision, collision_remainder = raw.quo_rem(Z_u**4)
assert not collision_remainder
factorization_after_collision = list(after_collision.factor())
print(
    "Q8O376RADICAND|degree={}|factors={}".format(
        after_collision.degree(),
        [(factor.degree(), int(exponent)) for factor, exponent in factorization_after_collision],
    ),
    flush=True,
)

# Refactor in QQ[u][T] so the squareclass keeps its global scalar and gives a
# polynomial binary quartic directly.
TUP = PolynomialRing(UQ, "T")
after_polynomial = TUP(after_collision)
factorization_polynomial = after_polynomial.factor()
factor_list_polynomial = list(factorization_polynomial)
assert sorted(
    (factor.degree(), int(exponent))
    for factor, exponent in factor_list_polynomial
) == [(1, 2), (4, 1), (4, 2)]
odd_factors = [
    factor for factor, exponent in factor_list_polynomial if int(exponent) % 2
]
assert len(odd_factors) == 1 and odd_factors[0].degree() == 4
quartic = TUP(factorization_polynomial.unit()*odd_factors[0])
square_factor = TUP.one()
for factor, exponent in factor_list_polynomial:
    square_factor *= factor**(int(exponent)//2)
assert after_polynomial == quartic*square_factor**2
e, d, c, b, a = quartic.list()
I = UQ(12*a*e-3*b*d+c**2)
J = UQ(72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3)
A_child = UQ(-27*I)
B_child = UQ(-27*J)
Delta_child = UQ(-16*(4*A_child**3+27*B_child**2))
assert (A_child.degree(), B_child.degree(), Delta_child.degree()) == (8, 12, 22)
delta_factorization = list(Delta_child.factor())
double_factors = [
    factor for factor, exponent in delta_factorization
    if int(exponent) == 2 and factor.degree() == 1
]
nodal_factors = [
    factor for factor, exponent in delta_factorization
    if int(exponent) == 1
]
assert len(double_factors) == 3
assert len(nodal_factors) == 1 and nodal_factors[0].degree() == 16
assert nodal_factors[0].is_squarefree()
assert all(A_child.gcd(factor) == B_child.gcd(factor) == 1 for factor, unused in delta_factorization)
assert A_child.leading_coefficient() and B_child.leading_coefficient()
# Degree 22 gives the fourth I2 at infinity.
assert (24-Delta_child.degree(), 8-A_child.degree(), 12-B_child.degree()) == (2, 0, 0)

# Both old I4 chains give rational points on the quartic.  Point at T=0 and
# verify that the standard degree-one pointed-quartic construction has exactly
# the invariant Jacobian above, rather than the multiplication-by-two covariant.
assert KU(e).is_square() and KU(a).is_square()
zero_coordinate = KU(0)
zero_ordinate = KU(e).sqrt()
assert zero_ordinate**2 == KU(e)


def evaluate_u_polynomial_series(poly, value):
    answer = LS.zero()
    for coefficient in reversed(UQ(poly).list()):
        answer = answer*value + LS(Kr(QQ(coefficient)))
    return answer


def evaluate_bivariate_series(poly, value, old_t_support=QQ(0)):
    answer = LS.zero()
    old_t_series = LS(Kr(QQ(old_t_support)))+s
    for old_t_degree, coefficient in enumerate(TUP(poly).list()):
        answer += evaluate_u_polynomial_series(coefficient, value)*old_t_series**old_t_degree
    return answer


def evaluate_u_function_at_r(value, u_of_r):
    value = KU(value)

    def evaluate_polynomial(poly):
        answer = Kr.zero()
        for coefficient in reversed(UQ(poly).list()):
            answer = answer*u_of_r + Kr(QQ(coefficient))
        return answer

    return evaluate_polynomial(value.numerator())/evaluate_polynomial(value.denominator())


def rational_map_degree(value):
    value = Kr(value)
    return max(value.numerator().degree(), value.denominator().degree())


# Restrict the resolved pencil and quartic square root to each component of
# the ordered finite I4.  This gives the exact new-base map u(r) and fixes the
# quartic sign on every degree-one component.
finite_component_maps = []
for component, label in enumerate(finite_i4_component_labels, start=1):
    local = i4_component_local_series(False, finite_rho_component_one, component)
    restrictions = []
    for AA, BB in resolved_pairs:
        AAloc = weighted_local(AA, 8, False)
        BBloc = weighted_local(BB, 2, False)
        restrictions.append(AAloc+BBloc*local["Zloc"]*local["m"])
    u_series = -restrictions[0]/restrictions[1]
    assert int(u_series.valuation()) == 0
    u_of_r = Kr(u_series[0])
    u_degree = rational_map_degree(u_of_r)
    aa_series = (
        weighted_local(AA0, 8, False)
        + u_series*weighted_local(AA1, 8, False)
    )
    bb_series = (
        weighted_local(BB0, 2, False)
        + u_series*weighted_local(BB1, 2, False)
    )
    assert (aa_series+bb_series*local["Zloc"]*local["m"]).valuation() >= PREC-5
    square_series = evaluate_bivariate_series(square_factor, u_series)
    quartic_series = evaluate_bivariate_series(quartic, u_series)
    W_series = (
        bb_series**2
        * (2*local["x"]+local["Hxloc"]-local["m"]**2)
        / square_series
    )
    print(
        "Q8O376FINITEI4|component={}|label={}|restriction_vals={}|u_val={}|"
        "square_val={}|W_val={}".format(
            component, label, [int(value.valuation()) for value in restrictions],
            int(u_series.valuation()), int(square_series.valuation()),
            int(W_series.valuation()),
        ),
        flush=True,
    )
    assert (W_series**2-quartic_series).valuation() >= PREC-6
    W_valuation = int(W_series.valuation())
    W_of_r = Kr(W_series[W_valuation])
    if u_degree == 1:
        assert W_valuation == 0
        positive_of_r = evaluate_u_function_at_r(zero_ordinate, u_of_r)
        if W_of_r == positive_of_r:
            quartic_sign = 1
        elif W_of_r == -positive_of_r:
            quartic_sign = -1
        else:
            quartic_sign = 0
    else:
        quartic_sign = 0
    finite_component_maps.append({
        "component_number": component,
        "label": label,
        "new_base_degree": int(u_degree),
        "u_of_r": str(u_of_r),
        "W_of_r": str(W_of_r),
        "quartic_sign_against_selected_sqrt": quartic_sign,
    })

degree_one_finite_components = [
    item for item in finite_component_maps if item["new_base_degree"] == 1
]
assert len(degree_one_finite_components) == 2
assert sorted(item["quartic_sign_against_selected_sqrt"] for item in degree_one_finite_components) == [-1, 1]
selected_origin_component = next(
    item for item in degree_one_finite_components
    if item["quartic_sign_against_selected_sqrt"] == 1
)
opposite_origin_component = next(
    item for item in degree_one_finite_components
    if item["quartic_sign_against_selected_sqrt"] == -1
)
print(
    "Q8O376ORIGIN|selected={}|opposite={}|rho_component1={}".format(
        selected_origin_component["label"], opposite_origin_component["label"],
        finite_rho_component_one,
    ),
    flush=True,
)


# P1229 is the nonidentity component of the first compact I2.  Its toric
# restriction is a degree-one point of the same quartic.  Select the exact
# square-root sign and point the quartic there, avoiding any large section
# translation on the child Jacobian.
p1229_Aloc = finite_support_local(A, p1229_old_base_support)
p1229_Bloc = finite_support_local(B, p1229_old_base_support)
p1229_node = QQ(-3*QQ(p1229_Bloc[0])/(2*QQ(p1229_Aloc[0])))
p1229_center = newton_sqrt(-p1229_Aloc/3, p1229_node)
p1229_nodal_error = p1229_center**3+p1229_Aloc*p1229_center+p1229_Bloc
assert int(p1229_nodal_error.valuation()) == 2
p1229_conic_constant = QQ(p1229_nodal_error[2])
p1229_exceptional_conic = Conic(
    QQ, [-3*p1229_node, QQ(1), -p1229_conic_constant],
)
p1229_has_point, p1229_conic_point = p1229_exceptional_conic.has_rational_point(point=True)
assert p1229_has_point
p1229_conic_coordinates = [QQ(value) for value in p1229_conic_point]
p1229_X1 = p1229_conic_coordinates[0]/p1229_conic_coordinates[2]
p1229_Y1 = p1229_conic_coordinates[1]/p1229_conic_coordinates[2]
assert p1229_Y1**2 == 3*p1229_node*p1229_X1**2+p1229_conic_constant
p1229_x_series = p1229_center+s*LS(Kr(p1229_X1))
p1229_y_square_over_s2 = (
    p1229_x_series**3+p1229_Aloc*p1229_x_series+p1229_Bloc
)/s**2
p1229_y_series = s*newton_sqrt(p1229_y_square_over_s2, Kr(p1229_Y1))
p1229_Zloc = finite_support_local(Z, p1229_old_base_support)
p1229_Hxloc = finite_support_local(X, p1229_old_base_support)/p1229_Zloc**2
p1229_Hyloc = finite_support_local(Y, p1229_old_base_support)/p1229_Zloc**3
p1229_m_series = (
    p1229_y_series+p1229_Hyloc
)/(p1229_x_series-p1229_Hxloc)
p1229_local = {
    "Zloc": p1229_Zloc,
    "Hxloc": p1229_Hxloc,
    "Hyloc": p1229_Hyloc,
    "x": p1229_x_series,
    "y": p1229_y_series,
    "m": p1229_m_series,
}
p1229_restrictions = []
for AA, BB in resolved_pairs:
    AAloc = finite_support_local(AA, p1229_old_base_support)
    BBloc = finite_support_local(BB, p1229_old_base_support)
    p1229_restrictions.append(
        AAloc+BBloc*p1229_local["Zloc"]*p1229_local["m"]
    )
p1229_u_series = -p1229_restrictions[0]/p1229_restrictions[1]
assert int(p1229_u_series.valuation()) == 0
p1229_u_of_r = Kr(p1229_u_series[0])
assert candidate["exact_curve_degrees"]["P1229"] == 1
p1229_aa_series = (
    finite_support_local(AA0, p1229_old_base_support)
    + p1229_u_series*finite_support_local(AA1, p1229_old_base_support)
)
p1229_bb_series = (
    finite_support_local(BB0, p1229_old_base_support)
    + p1229_u_series*finite_support_local(BB1, p1229_old_base_support)
)
assert (
    p1229_aa_series
    + p1229_bb_series*p1229_local["Zloc"]*p1229_local["m"]
).valuation() >= PREC-5
p1229_square_series = evaluate_bivariate_series(
    square_factor, p1229_u_series, p1229_old_base_support,
)
p1229_quartic_series = evaluate_bivariate_series(
    quartic, p1229_u_series, p1229_old_base_support,
)
p1229_W_series = (
    p1229_bb_series**2
    * (2*p1229_local["x"]+p1229_local["Hxloc"]-p1229_local["m"]**2)
    / p1229_square_series
)
assert int(p1229_W_series.valuation()) == 0
assert (p1229_W_series**2-p1229_quartic_series).valuation() >= PREC-6
p1229_W_of_r = Kr(p1229_W_series[0])
p1229_special = KU(TUP(quartic)(KU(p1229_old_base_support)))
assert p1229_special.is_square()
p1229_positive_ordinate = KU(p1229_special.sqrt())
if p1229_W_of_r == evaluate_u_function_at_r(p1229_positive_ordinate, p1229_u_of_r):
    p1229_ordinate = p1229_positive_ordinate
    p1229_quartic_sign = 1
elif p1229_W_of_r == -evaluate_u_function_at_r(p1229_positive_ordinate, p1229_u_of_r):
    p1229_ordinate = -p1229_positive_ordinate
    p1229_quartic_sign = -1
else:
    raise ArithmeticError("P1229 toric component does not select a quartic square-root sign")


def pointed_model_at(old_t_coordinate, ordinate):
    quartic_over_KU = TU(quartic)
    translated = [
        KU(quartic_over_KU.derivative(order)(KU(old_t_coordinate))/factorial(order))
        for order in range(5)
    ]
    ee, dd, cc, bb_value, aa_value = translated
    ordinate = KU(ordinate)
    assert ee == ordinate**2
    aa1 = dd/ordinate
    aa2 = cc-dd**2/(4*ordinate**2)
    aa3 = 2*ordinate*bb_value
    aa4 = -4*ordinate**2*aa_value
    aa6 = aa2*aa4
    bb2 = aa1**2+4*aa2
    bb4 = 2*aa4+aa1*aa3
    bb6 = aa3**2+4*aa6
    cc4 = bb2**2-24*bb4
    cc6 = -bb2**3+36*bb2*bb4-216*bb6
    return {
        "a_invariants": (aa1, aa2, aa3, aa4, aa6),
        "pointed_A": -cc4/48,
        "pointed_B": -cc6/864,
    }


p1229_pointed_model = pointed_model_at(p1229_old_base_support, p1229_ordinate)
assert 81*p1229_pointed_model["pointed_A"] == KU(A_child)
assert 729*p1229_pointed_model["pointed_B"] == KU(B_child)
print(
    "Q8O376P1229ZERO|support={}|conic_point={}|sign={}|u_degree=1|A=True|B=True".format(
        p1229_old_base_support, p1229_conic_point, p1229_quartic_sign,
    ),
    flush=True,
)
a1 = KU(d)/zero_ordinate
a2 = KU(c)-KU(d)**2/(4*zero_ordinate**2)
a3 = 2*zero_ordinate*KU(b)
a4 = -4*zero_ordinate**2*KU(a)
a6 = a2*a4
b2_pointed = a1**2+4*a2
b4_pointed = 2*a4+a1*a3
b6_pointed = a3**2+4*a6
c4_pointed = b2_pointed**2-24*b4_pointed
c6_pointed = -b2_pointed**3+36*b2_pointed*b4_pointed-216*b6_pointed
pointed_A = -c4_pointed/48
pointed_B = -c6_pointed/864
assert 81*pointed_A == KU(A_child)
assert 729*pointed_B == KU(B_child)
print(
    "Q8O376CHILDPREFLIGHT|quartic_degree={}|coeff_degrees={}|Adeg={}|Bdeg={}|Ddeg={}|"
    "delta_factors={}|constant_square={}|leading_square={}".format(
        quartic.degree(), [value.degree() for value in quartic.list()],
        A_child.degree(), B_child.degree(), Delta_child.degree(),
        [(factor.degree(), int(exponent)) for factor, exponent in delta_factorization],
        bool(KU(e).is_square()), bool(KU(a).is_square()),
    ),
    flush=True,
)

# Every CRT prime is a good-reduction regression for the same four-plane.
prime_regressions = []
for prime in horizontal["primes"]:
    p = ZZ(prime)
    F = GF(p)

    def red(value):
        value = QQ(value)
        assert value.denominator() % p
        return F(value.numerator())/F(value.denominator())

    reduced_conditions = matrix(F, [[red(value) for value in row] for row in condition_matrix.rows()])
    reduced_basis = matrix(F, [[red(value) for value in row] for row in ambient_rows.rows()])
    assert reduced_conditions.rank() == 8
    assert reduced_basis.rank() == 4
    assert reduced_conditions*reduced_basis.transpose() == matrix(F, 8, 4)
    prime_regressions.append(int(p))


def u_rational_record(value):
    value = KU(value)
    return {
        "numerator_coefficients_low_to_high": coefficients(UQ(value.numerator())),
        "denominator_coefficients_low_to_high": coefficients(UQ(value.denominator())),
    }

payload = {
    "schema": "elkies-k3.h3-q4o164-q8o376-resolved-rr-qq.v1",
    "status": "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO",
    "selected_horizontal": {
        "identity": horizontal["selected_identity"],
        "P_dot_O": 4,
        "height": "11",
        "projective_degrees_X_Y_Z": [12, 18, 4],
        "exact_weierstrass_identity": True,
    },
    "smooth_RR": {
        "ambient": "a=AA/Z^2, deg AA<=8; b=BB/Z, deg BB<=2",
        "ambient_dimension": 12,
        "collision_modulus": "Z^2",
        "collision_modulus_degree": 8,
        "condition": "AA*X-BB*Y == 0 mod Z^2",
        "condition_rank": 8,
        "kernel_dimension": 4,
        "basis_pairs": [{
            "AA_coefficients_low_to_high": [str(value) for value in AA.list()],
            "BB_coefficients_low_to_high": [str(value) for value in BB.list()],
        } for AA, BB in pairs],
        "maximum_basis_rational_bits": rational_bits(ambient_rows.list()),
        "regression_primes": prime_regressions,
    },
    "resolved_RR": {
        "lattice_divisor": {
            "fibre_in_source": candidate["fibre"],
            "horizontal_section_in_source": candidate["horizontal"]["section"],
            "vertical_coefficients": candidate["horizontal"]["vertical"],
            "connected_full_A3_root_indices": [list(first_a3), list(second_a3)],
        },
        "vertical_quotients": vertical_diagnostics,
        "vertical_condition_rank_on_smooth_four_plane": 2,
        "resolved_kernel_on_smooth_basis": [
            [str(value) for value in row] for row in resolved_kernel.rows()
        ],
        "resolved_basis_pairs": [{
            "AA_coefficients_low_to_high": coefficients(AA),
            "BB_coefficients_low_to_high": coefficients(BB),
        } for AA, BB in resolved_pairs],
        "dimensions": {"ambient": 12, "smooth": 4, "h0": 2},
    },
    "quartic": {
        "equation": "W^2 = sum_i quartic_coefficients[i](u)*T^i",
        "degree_in_old_T": 4,
        "coefficients_in_old_T_low_to_high": [coefficients(value) for value in quartic.list()],
        "coefficient_degrees_in_u": [int(value.degree()) for value in quartic.list()],
        "raw_after_collision_degree": int(after_collision.degree()),
        "raw_factor_degrees_and_exponents": [
            [int(factor.degree()), int(exponent)]
            for factor, exponent in factor_list_polynomial
        ],
        "maximum_rational_bits": rational_bits(
            coefficient for value in quartic.list() for coefficient in value.list()
        ),
    },
    "pointed_origin": {
        "old_base_coordinate": "T=0",
        "quartic_ordinate": u_rational_record(zero_ordinate),
        "exact_quartic_identity": True,
        "degree_one_pointed_quartic_map": True,
        "short_scaling": "x_short=9*x_pointed, y_short=27*y_pointed",
        "exact_A_identity": "81*A_pointed = A_child",
        "exact_B_identity": "729*B_pointed = B_child",
        "lattice_label": selected_origin_component["label"],
        "opposite_quartic_sign_lattice_label": opposite_origin_component["label"],
        "finite_I4_marking": {
            "equation_support": "T=0",
            "selected_marked_embedding_index": 15,
            "marked_root_indices_in_component_order": list(finite_i4_root_indices),
            "component_labels_in_order": list(finite_i4_component_labels),
            "B0_component_profile": b0_record["component_profile"],
            "rho_component_one": str(finite_rho_component_one),
            "rho_orientation_trials": [
                [str(rho_start), left_valuation, right_valuation]
                for rho_start, left_valuation, right_valuation in finite_rho_trials
            ],
            "component_maps": finite_component_maps,
        },
    },
    "preferred_pointed_zero": {
        "lattice_label": "P1229",
        "old_base_coordinate": str(p1229_old_base_support),
        "source_component_identification": (
            "nonidentity component of the first compact finite I2; exact marked "
            "A1 swap in selected embedding 15"
        ),
        "source_exact_q8_degree": int(candidate["exact_curve_degrees"]["P1229"]),
        "exceptional_conic": {
            "equation": "Y1^2 = 3*node*X1^2 + c0",
            "node": str(p1229_node),
            "c0": str(p1229_conic_constant),
            "rational_point_X1_Y1": [str(p1229_X1), str(p1229_Y1)],
        },
        "local_new_base_value": str(p1229_u_of_r),
        "quartic_sign": p1229_quartic_sign,
        "quartic_ordinate": u_rational_record(p1229_ordinate),
        "degree_one_pointed_quartic_map": True,
        "generalized_weierstrass_a_invariants": [
            u_rational_record(value)
            for value in p1229_pointed_model["a_invariants"]
        ],
        "short_scaling": "x_short=9*x_pointed, y_short=27*y_pointed",
        "exact_A_identity": "81*A_P1229_pointed = A_child",
        "exact_B_identity": "729*B_P1229_pointed = B_child",
        "maximum_pointing_rational_bits": rational_bits(
            coefficient
            for value in (
                (p1229_ordinate,) + p1229_pointed_model["a_invariants"]
            )
            for polynomial in (value.numerator(), value.denominator())
            for coefficient in polynomial.list()
        ),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "degrees_A_B_Delta": [8, 12, 22],
        "finite_reducible_fibres": [{
            "factor": str(factor.monic()),
            "kodaira": "I2",
            "delta_order": 2,
        } for factor in double_factors],
        "finite_nodal_factor": str(nodal_factors[0].monic()),
        "finite_nodal_factor_degree": 16,
        "finite_nodal_factor_squarefree": True,
        "infinity": {"kodaira": "I2", "orders_A_B_Delta": [0, 0, 2]},
        "fibre_profile": "4I2+16I1",
        "ADE": "4A1",
        "root_data": [4, 8, 16],
        "MW_rank_if_rho19": 13,
        "euler_number": 24,
        "maximum_A_B_rational_bits": rational_bits(A_child.list()+B_child.list()),
    },
    "method": {
        "large_Groebner_required": False,
        "nonlinear_characteristic_zero_solve_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "proof_boundary": (
        "The exact QQ q8/o376 horizontal, 12-to-4 smooth saturation, two orientation-independent "
        "connected-A3 quotient rows, h0=2 plane, quartic, degree-one pointed map, Jacobian, "
        "4I2+16I1 fibres, Euler 24 and 4A1/MW13 in the pinned rank-19 NS lattice are certified. "
        "The T=0 rational origin is assigned exactly to its oriented inherited degree-one "
        "component by the selected marked embedding, the exact B0 component-one tangent, and "
        "split-toric component maps. The lattice-selected P1229 curve is identified with the "
        "first compact I2 component, its exact exceptional-conic arc selects the quartic sign, "
        "and pointing there gives the same short child invariants. Thus the q8 equation is now "
        "marked with the preferred P1229 zero. The four q12/o5867 compiler-section equations and "
        "its resolved RR pencil remain open."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q8O376RRQQ|ambient=12|smooth=4|vertical_rank=2|h0=2|quartic=4|"
    "fibres=4I2+16I1|ADE=4A1|status={}|output={}".format(payload["status"], OUTPUT),
    flush=True,
)
