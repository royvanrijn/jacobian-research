#!/usr/bin/env sage
"""Point the corrected q323 quartic at old_A11_component_2.

status: ACTIVE_PROOF
claim: exact component-2 split-I4 sign and pointed A3+2A2 Weierstrass model
inputs: q4/o208 compact model, exact q323 horizontal, exact q323 RR quartic
outputs: artifacts/local/elkies-k3/q4o323-component2-pointing-qq.json

The old t=0 fibre is split I4.  On its orientation-invariant middle
component use the toric arc left=r*s^2, right=kappa*r^-1*s^2.  Evaluating the
exact q323 pencil and normalized chord ordinate on that arc distinguishes the
two square roots of Q(0,u).  The resulting rational quartic point is converted
to the global short Jacobian by the standard pointed-binary-quartic formulas.
No elimination or Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import (
    LaurentSeriesRing, PolynomialRing, QQ, ZZ,
)


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
if not (ROOT / "MATH_STATUS.json").is_file():
    ROOT = Path.cwd().resolve()
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
RR = LOCAL / "q4o208-q4o323-a3-2a2-rr-qq.json"
OUTPUT = LOCAL / "q4o323-component2-pointing-qq.json"
INPUTS = (COMPACT, HORIZONTAL, RR)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_function_record(value, ring):
    value = ring.fraction_field()(value)
    return {
        "numerator_coefficients_low_to_high": [str(item) for item in value.numerator().list()],
        "denominator_coefficients_low_to_high": [str(item) for item in value.denominator().list()],
    }


compact = json.loads(COMPACT.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
rr = json.loads(RR.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING"
assert rr["status"] == "PASS_EXACT_QQ_Q4O323_A3_2A2_RR_AND_JACOBIAN"

R = PolynomialRing(QQ, "t")
t = R.gen()
RU = PolynomialRing(QQ, "u")
u = RU.gen()
KU = RU.fraction_field()
TK = PolynomialRing(KU, "t")
tt = TK.gen()

A = R([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
x_record = horizontal["exact_QQ_horizontal"]["x"]
y_record = horizontal["exact_QQ_horizontal"]["y"]
X = R([QQ(value) for value in x_record["numerator_coefficients_low_to_high"]])
Z2 = R([QQ(value) for value in x_record["denominator_coefficients_low_to_high"]])
Y = R([QQ(value) for value in y_record["numerator_coefficients_low_to_high"]])
Z3 = R([QQ(value) for value in y_record["denominator_coefficients_low_to_high"]])
assert Z2.is_square()
Z = Z2.sqrt()
assert Z**3 == Z3
assert Y**2 == X**3 + A*X*Z**4 + B*Z**6

kernel = [[QQ(value) for value in row] for row in rr["resolved_RR"]["kernel_basis"]]
assert len(kernel) == 2 and all(len(row) == 5 for row in kernel)
pairs = []
for row in kernel:
    pairs.append((R(row[:4]), QQ(row[4])))
(a0, b0), (a1, b1) = pairs
a_u = TK(a0) + KU(u)*TK(a1)
b_u = KU(b0) + KU(u)*KU(b1)
raw = (
    a_u**4 - 6*TK(X)*a_u**2*b_u**2 + 8*TK(Y)*a_u*b_u**3
    - 3*TK(X)**2*b_u**4 - 4*TK(A)*b_u**4*TK(Z)**4
)
after_pole, remainder = raw.quo_rem(TK(Z)**4)
assert not remainder
factorization = after_pole.factor()
linear_squares = [factor for factor, exponent in factorization if factor.degree() == 1 and exponent == 2]
quartics = [factor for factor, exponent in factorization if factor.degree() == 4 and exponent == 1]
assert len(linear_squares) == 2 and len(quartics) == 1
square_factor = linear_squares[0]*linear_squares[1]

quartic = TK([
    KU(RU([QQ(value) for value in coefficient]))
    for coefficient in rr["quartic"]["coefficients_in_t_low_to_high"]
])
quartic_factor_scale = KU(quartics[0].leading_coefficient()/quartic.leading_coefficient())
assert quartics[0] == quartic_factor_scale*quartic
normalizing_unit = KU(after_pole/(square_factor**2*quartic))
assert normalizing_unit.is_square()
unit_root = KU(normalizing_unit.sqrt())
assert after_pole == (unit_root*square_factor)**2*quartic

e, d, c, b, a = quartic.list()
e_poly = RU(e)
assert e_poly.is_square()
L0 = e_poly.sqrt()
assert L0.degree() == 1 and L0**2 == e_poly

# Exact split-I4 component-2 arc at t=0.
r_ring = PolynomialRing(QQ, "r")
r = r_ring.gen()
Kr = r_ring.fraction_field()
PRECISION = 12
LS = LaurentSeriesRing(Kr, "s", default_prec=PRECISION)
s = LS.gen()


def shifted(poly):
    return LS([Kr(value) for value in poly.list()])


def newton_sqrt(value, constant_root):
    answer = LS(Kr(constant_root))
    for unused in range(6):
        answer = (answer + value/answer)/2
    assert (answer**2-value).valuation() >= PRECISION-3
    return answer


Aloc, Bloc = shifted(A), shifted(B)
Xloc, Yloc, Zloc = shifted(X), shifted(Y), shifted(Z)
A0, B0 = QQ(A(0)), QQ(B(0))
node = QQ(-3*B0/(2*A0))
assert node**3+A0*node+B0 == 0 and 3*node**2+A0 == 0
center = newton_sqrt(-Aloc/3, node)
nodal_error = center**3+Aloc*center+Bloc
assert int(nodal_error.valuation()) == 4
kappa = nodal_error/s**4
rho_square = QQ(3*node)
assert rho_square.is_square()
rho0 = QQ(rho_square.sqrt())


def eval_ku_series(value, u_series):
    value = KU(value)
    numerator = RU(value.numerator())(u_series)
    denominator = RU(value.denominator())(u_series)
    return numerator/denominator


def eval_tk_series(poly, u_series):
    return sum(
        (eval_ku_series(poly[index], u_series)*s**index for index in range(poly.degree()+1)),
        LS.zero(),
    )


orientation_records = []
selected_signs = []
u_degrees = None
for rho_start in (rho0, -rho0):
    left = LS(Kr(r))*s**2
    right = kappa*s**2/LS(Kr(r))
    yy = (left+right)/2
    ww = (right-left)/2
    rho = LS(Kr(rho_start))
    for unused in range(6):
        rho -= (rho**3-3*center*rho-ww)/(3*rho**2-3*center)
    xx = center+ww/rho
    assert (yy**2-(xx**3+Aloc*xx+Bloc)).valuation() >= PRECISION-3
    Hx, Hy = Xloc/Zloc**2, Yloc/Zloc**3
    # The quartic radicand used by the RR compiler is the residual quadratic
    # for the line through H: slope=(y-y_H)/(x-x_H).  Its linear term is
    # +8*y_H*slope, matching the stored raw polynomial.
    slope = (yy-Hy)/(xx-Hx)
    scaled_slope = Zloc*slope
    a0loc, a1loc = shifted(a0), shifted(a1)
    u_loc = (Kr(b0)*scaled_slope-a0loc)/(a1loc-Kr(b1)*scaled_slope)
    print(
        "Q4O323C2DIAG|rho={}|slope_val={}|scaled_val={}|num_val={}|den_val={}|u_val={}".format(
            rho_start, slope.valuation(), scaled_slope.valuation(),
            (Kr(b0)*scaled_slope-a0loc).valuation(),
            (a1loc-Kr(b1)*scaled_slope).valuation(), u_loc.valuation(),
        ),
        flush=True,
    )
    assert int(u_loc.valuation()) == 0
    u_of_r = Kr(u_loc[0])
    current_u_degrees = (
        int(u_of_r.numerator().degree()), int(u_of_r.denominator().degree())
    )
    if u_degrees is None:
        u_degrees = current_u_degrees

    branch = 2*xx+Hx-slope**2
    b_loc = Kr(b0)+u_loc*Kr(b1)
    square_loc = eval_tk_series(square_factor, u_loc)
    unit_loc = eval_ku_series(unit_root, u_loc)
    w_loc = b_loc**2*branch/(square_loc*unit_loc)
    assert int(w_loc.valuation()) == 0
    w_of_r = Kr(w_loc[0])
    l_of_r = Kr(L0(u_of_r))
    ratio = Kr(w_of_r/l_of_r)
    assert ratio in (Kr.one(), -Kr.one())
    selected_signs.append(int(ratio))
    orientation_records.append({
        "rho_start": str(rho_start),
        "u_of_r": rational_function_record(u_of_r, r_ring),
        "u_degrees": list(current_u_degrees),
        "w_of_r": rational_function_record(w_of_r, r_ring),
        "sign_relative_to_stored_L0": int(ratio),
    })

# Component 2 is fixed by reversing the I4 orientation, so both tangent signs
# must give the same quartic branch.
assert len(set(selected_signs)) == 1
component2_sign = selected_signs[0]
w0 = KU(component2_sign)*KU(L0)
assert w0**2 == e

# Pointed quartic -> generalized Weierstrass -> global short model.
a_1 = d/w0
a_2 = c-d**2/(4*w0**2)
a_3 = 2*w0*b
a_4 = -4*w0**2*a
a_6 = a_2*a_4
b_2 = a_1**2+4*a_2
b_4 = 2*a_4+a_1*a_3
b_6 = a_3**2+4*a_6
c_4 = b_2**2-24*b_4
c_6 = -b_2**3+36*b_2*b_4-216*b_6
pointed_A = -c_4/48
pointed_B = -c_6/864
A_child = KU(RU([QQ(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]]))
B_child = KU(RU([QQ(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]]))
assert 81*pointed_A == A_child
assert 729*pointed_B == B_child

# The opposite t=0 branch becomes an exact section after pointing at C2.
x_general = -a_2
y_general = a_1*a_2-a_3
x_opposite = KU(9*(x_general+b_2/12))
y_opposite = KU(27*(y_general+(a_1*x_general+a_3)/2))
assert y_opposite**2 == x_opposite**3+A_child*x_opposite+B_child

payload = {
    "schema": "elkies-k3.h3-q4o323-component2-pointing-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING",
    "zero": "old_A11_component_2",
    "old_fibre": "second_old_I6_I4",
    "old_component_index": int(2),
    "quartic_square_at_t0": {
        "L0_coefficients_low_to_high": [str(value) for value in L0.list()],
        "selected_sign": int(component2_sign),
    },
    "split_I4_component2_arc": {
        "orientation_invariant": True,
        "orientation_records": orientation_records,
    },
    "pointed_generalized_weierstrass": {
        "a1": rational_function_record(a_1, RU),
        "a2": rational_function_record(a_2, RU),
        "a3": rational_function_record(a_3, RU),
        "a4": rational_function_record(a_4, RU),
        "a6": rational_function_record(a_6, RU),
    },
    "global_short_model": {
        "A_coefficients_low_to_high": [str(value) for value in RU(A_child).list()],
        "B_coefficients_low_to_high": [str(value) for value in RU(B_child).list()],
        "pointing_scales": {"A": int(81), "B": int(729)},
    },
    "opposite_t0_branch_section": {
        "x": rational_function_record(x_opposite, RU),
        "y": rational_function_record(y_opposite, RU),
    },
    "quartic_map": {
        "rr_pairs": [
            {"a_coefficients_low_to_high": [str(value) for value in apoly.list()], "b": str(bvalue)}
            for apoly, bvalue in pairs
        ],
        "square_factor_coefficients_in_t_low_to_high": [
            rational_function_record(value, RU) for value in square_factor.list()
        ],
        "normalizing_unit_root": rational_function_record(unit_root, RU),
    },
    "checks": {
        "component2_sign_orientation_independent": True,
        "quartic_point_exact": True,
        "global_A_identity": True,
        "global_B_identity": True,
        "opposite_branch_section_exact": True,
        "large_Groebner_required": False,
    },
    "runtime_seconds": float(time.monotonic()-started),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
    "proof_boundary": (
        "Exact component-2 zero, quartic sign, pointed global short equation, and opposite "
        "t=0 section. Recovery of the outgoing physical q207 horizontal remains separate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323C2POINTQQ|sign={}|u_degrees={}/{}|opposite_x={}/{}|opposite_y={}/{}|"
    "status={}|output={}".format(
        component2_sign,
        u_degrees[0], u_degrees[1],
        x_opposite.numerator().degree(), x_opposite.denominator().degree(),
        y_opposite.numerator().degree(), y_opposite.denominator().degree(),
        payload["status"], OUTPUT,
    ),
    flush=True,
)
