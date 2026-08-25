#!/usr/bin/env sage -python
"""Smoke-test the inherited-P1 degree-seven Abel trace on q4/orbit164.

The exact polynomial P1 section on the 2A5 model is carried through the
q4/orbit208, q4/orbit1584, and q4/orbit164 binary quartics.  On a specialized
q4/orbit164 fibre the resulting degree-seven divisor is Abel-reduced by the
unique relation in

    L(8 O) = <1,x,x^2,x^3,x^4,y,xy,x^2y>.

This is a finite-field construction smoke test.  It uses only rational
function arithmetic and a 7-by-8 kernel; no Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
P146 = LOCAL / "q24-2a5-p146-p1307-scaled-x-qq.json"
P1229 = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
Q208 = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
Q208_MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
Q1584 = LOCAL / "q4o208-physical-q4o1584-rr-qq.json"
Q1584_MARKING = LOCAL / "q4o1584-second-affine-equation-marking-qq.json"
Q164 = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
Q164_MARKING = LOCAL / "q4o164-c8-equation-marking-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=131)
parser.add_argument("--tau", type=int, default=2)
parser.add_argument(
    "--interpolate", action="store_true",
    help="interpolate the expected (32,28)/(48,42) modular trace section",
)
parser.add_argument(
    "--output",
    type=Path,
)
parser.add_argument(
    "--good-fibre-limit", type=int,
    help=(
        "stop after this many good fibres; interpolation requires at least 92. "
        "The default scans the whole prime field"
    ),
)
args = parser.parse_args()
if args.good_fibre_limit is not None and args.good_fibre_limit < 92:
    parser.error("--good-fibre-limit must be at least 92")

started = time.monotonic()
load(str(ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"))


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O164P1TRACE|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


data_p146 = json.loads(P146.read_text())
data_p1229 = json.loads(P1229.read_text())
q208 = json.loads(Q208.read_text())
mark208 = json.loads(Q208_MARKING.read_text())
q1584 = json.loads(Q1584.read_text())
mark1584 = json.loads(Q1584_MARKING.read_text())
q164 = json.loads(Q164.read_text())
mark164 = json.loads(Q164_MARKING.read_text())
assert data_p146["status"] == "PASS_EXACT_QQ_P146_AND_P1307_SHORT_MW_WORDS"
assert data_p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert q208["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert mark208["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert q1584["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN"
assert mark1584["status"] == "PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING"
assert q164["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert mark164["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"
log("INPUTS")

p = ZZ(args.prime)
F = GF(p)
R = PolynomialRing(F, "r")
r = R.gen()
K = R.fraction_field()


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad reduction at p={p}")
    return F(value.numerator()) / F(value.denominator())


def poly(values):
    return R([red(value) for value in values])


def rf_record(record):
    numerator = poly(record["numerator_coefficients_low_to_high"])
    denominator = poly(record["denominator_coefficients_low_to_high"])
    return K(numerator) / K(denominator)


def evaluate(polynomial, value):
    return sum(K(coefficient)*value**index for index, coefficient in enumerate(polynomial))


def substitute(value, argument):
    value = K(value)
    return evaluate(R(value.numerator()), argument) / evaluate(R(value.denominator()), argument)


def qq_function_list(texts, variable_name):
    QR = PolynomialRing(QQ, variable_name)
    variable = QR.gen()
    QK = QR.fraction_field()
    FR = PolynomialRing(F, variable_name)
    answer = []
    for text in texts:
        value = QK(sage_eval(str(text), locals={variable_name: variable}))
        numerator = FR([red(coefficient) for coefficient in QR(value.numerator()).list()])
        denominator = FR([red(coefficient) for coefficient in QR(value.denominator()).list()])
        answer.append((numerator, denominator))
    return answer


def evaluate_function(record, value):
    numerator, denominator = record
    nv = sum(K(coefficient)*value**index for index, coefficient in enumerate(numerator))
    dv = sum(K(coefficient)*value**index for index, coefficient in enumerate(denominator))
    if not dv:
        raise ZeroDivisionError("transition coefficient denominator vanished")
    return nv/dv


def pointed_quartic_point(
    coefficients, coordinate, ordinate, zero_coordinate, zero_ordinate,
    global_scale=1,
):
    """Use the certified degree-one pointed-quartic map.

    The invariant covariant map is a 2-cover map and therefore multiplies the
    divisor class.  Translating by a rational quartic zero instead preserves
    the primitive point represented by the multisection.
    """

    ZR = PolynomialRing(K, "z")
    quartic = ZR([K(value) for value in coefficients])
    assert ordinate**2 == quartic(coordinate)
    zero_coordinate = K(zero_coordinate)
    zero_ordinate = K(zero_ordinate)
    assert zero_ordinate**2 == quartic(zero_coordinate)
    translated = [
        K(quartic.derivative(order)(zero_coordinate) / factorial(order))
        for order in range(5)
    ]
    e, d, c, b, a = translated
    a1 = d/zero_ordinate
    a2 = c-d**2/(4*zero_ordinate**2)
    a3 = 2*zero_ordinate*b
    a4 = -4*zero_ordinate**2*a
    a6 = a2*a4
    b2 = a1**2+4*a2
    local = coordinate-zero_coordinate
    if not local:
        raise ArithmeticError("test multisection coincides with the pointed zero")
    x_general = (2*zero_ordinate*(ordinate+zero_ordinate)+d*local)/local**2
    y_general = (
        4*zero_ordinate**2*(ordinate+zero_ordinate)
        + 2*zero_ordinate*d*local
        + (2*zero_ordinate*c-d**2/(2*zero_ordinate))*local**2
    )/local**3
    assert (
        y_general**2+a1*x_general*y_general+a3*y_general
        == x_general**3+a2*x_general**2+a4*x_general+a6
    )
    global_scale = K(global_scale)
    x_short = global_scale**2*9*(x_general+b2/12)
    y_short = global_scale**3*27*(y_general+(a1*x_general+a3)/2)
    return K(x_short), K(y_short)


# Exact P1 and P1229 on the 2A5 equation, reduced directly from QQ.
p1 = data_p146["polynomial_inputs"]["P1"]
x0 = K(poly(p1["X_coefficients_low_to_high"]))
y0 = K(poly(p1["Y_coefficients_low_to_high"]))
p1229 = data_p1229["P1229"]
px0 = K(poly(p1229["X_coefficients_low_to_high"]))
py0 = K(poly(p1229["Y_coefficients_low_to_high"]))

# q4/o208: U=f1/f0 and the two old-I6 factors are the removed square.
kernel208 = matrix(F, [[red(value) for value in row] for row in q208["resolved_RR"]["kernel_basis"]])
m0 = (y0+py0)/(x0-px0)
functions208 = []
for row in kernel208.rows():
    functions208.append(K(row[0]+row[1]*r+row[2]*r**2)+K(row[3])*m0)
u1 = functions208[1]/functions208[0]
old_square = R.one()
for root in q208["resolved_RR"]["old_I6_roots"]:
    old_square *= r-red(root)
w1 = (2*x0+px0-m0**2)/K(old_square)
quartic208 = qq_function_list(q208["quartic"]["coefficients_in_old_T_low_to_high"], "U")
coeff208 = [evaluate_function(value, u1) for value in quartic208]
b0, b1 = K(kernel208[0, 3]), K(kernel208[1, 3])
zero208 = red(mark208["resolved_C5_slice"]["parent_I6_base_value"])
wzero208 = substitute(
    rf_record(mark208["resolved_C5_slice"]["quartic_ordinate_on_C5"]), u1
)
x1, y1 = pointed_quartic_point(
    coeff208, K(r), w1, zero208, wzero208,
    global_scale=(u1*b0-b1)**2,
)
log("Q4O208", degree=max(u1.numerator().degree(), u1.denominator().degree()))

# q4/o1584: resolved double branch at the second inherited I4.
point1 = mark208["first_I6_affine_component_on_C5_pointed_child"]
px1 = rf_record(point1["x"])
py1 = rf_record(point1["y"])
m1 = (y1+substitute(py1, u1))/(x1-substitute(px1, u1))
s1584 = red(q1584["resolved_RR"]["support"])
c0 = red(q1584["resolved_RR"]["double_branch_c0"])
c1 = red(q1584["resolved_RR"]["first_jet_c1"])
local1 = u1-s1584
u2 = (m1-c0-c1*local1)/local1**2
w2 = (2*x1+substitute(px1, u1)-m1**2)/local1**2
coeff1584 = [
    evaluate(poly(values), u2)
    for values in q1584["quartic"]["coefficients_in_L_low_to_high"]
]
lzero1584 = substitute(rf_record(mark1584["selected_zero"]["old_base_coordinate"]), u2)
wzero1584 = substitute(rf_record(mark1584["selected_zero"]["quartic_ordinate"]), u2)
x2, y2 = pointed_quartic_point(
    coeff1584, local1, w2, lzero1584, wzero1584
)
log("Q4O1584", degree=max(u2.numerator().degree(), u2.denominator().degree()))

# q4/o164: interpolate the two resolved I2 branch values.
point2 = mark1584["old_A11_component_0_on_second_affine_pointed_child"]
px2_raw = rf_record(point2["x"])
py2_raw = rf_record(point2["y"])
px2 = substitute(px2_raw, u2)
py2 = substitute(py2_raw, u2)
m2 = (y2+py2)/(x2-px2)
linear = poly(q164["resolved_RR"]["interpolating_linear_polynomial"])
denominator = poly(q164["resolved_RR"]["common_denominator"])
linear_value = evaluate(linear, u2)
denominator_value = evaluate(denominator, u2)
u3 = (m2-linear_value)/denominator_value
w3 = (2*x2+px2-m2**2)/denominator_value
coeff164 = [
    evaluate(poly(values), u3)
    for values in q164["quartic"]["coefficients_in_T_low_to_high"]
]
zero164 = red(mark164["selected_zero"]["parent_base_support"])
wzero164 = substitute(rf_record(mark164["selected_zero"]["quartic_ordinate"]), u3)
x3, y3 = pointed_quartic_point(coeff164, u2, w3, zero164, wzero164)
A3 = poly(q164["child"]["minimal_A_coefficients_low_to_high"])
B3 = poly(q164["child"]["minimal_B_coefficients_low_to_high"])
assert y3**2 == x3**3 + evaluate(A3, u3)*x3 + evaluate(B3, u3)
degree_u3 = max(u3.numerator().degree(), u3.denominator().degree())
if degree_u3 != 7:
    raise ArithmeticError(f"inherited P1 has final degree {degree_u3}, expected 7")
log("Q4O164", degree=degree_u3)


def reduce_mod_H(value, H):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    if denominator.gcd(H).degree() != 0:
        raise ZeroDivisionError("denominator not invertible on specialized degree-seven fibre")
    return (numerator*denominator.inverse_mod(H)) % H


def newton_power_sums(polynomial):
    degree = polynomial.degree()
    assert polynomial[degree] == 1
    sums = [F(degree)]
    for order in range(1, degree):
        total = F(order)*polynomial[degree-order]
        for index in range(1, order):
            total += polynomial[degree-index]*sums[order-index]
        sums.append(-total)
    return sums


def abel_trace_at(tau):
    tau = F(tau)
    Atau = A3(tau)
    Btau = B3(tau)
    E_tau = EllipticCurve(F, [0, 0, 0, Atau, Btau])
    if not E_tau.discriminant():
        raise ArithmeticError("selected q4/o164 fibre is singular")
    H = R(u3.numerator()-tau*u3.denominator())
    if H.degree() != 7:
        raise ArithmeticError("specialized inherited divisor dropped degree")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ArithmeticError("specialized inherited divisor is not etale")
    xA = reduce_mod_H(x3, H)
    yA = reduce_mod_H(y3, H)
    if (yA**2-xA**3-Atau*xA-Btau) % H:
        raise ArithmeticError("transported divisor misses specialized q4/o164 fibre")

    # L(8O): 1,x,x^2,x^3,x^4,y,xy,x^2y.
    xp = [R.one()]
    for unused in range(4):
        xp.append((xp[-1]*xA) % H)
    columns = xp + [(yA*xp[index]) % H for index in range(3)]
    evaluation = matrix(F, 7, 8, lambda row, column: columns[column][row])
    kernel = evaluation.right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError(f"L(8O) trace kernel dimension {kernel.nrows()}")
    relation = kernel[0]
    XR = PolynomialRing(F, "X")
    XX = XR.gen()
    Afun = sum(relation[index]*XX**index for index in range(5))
    Bfun = sum(relation[5+index]*XX**index for index in range(3))
    intersection = Afun**2-(XX**3+Atau*XX+Btau)*Bfun**2
    if intersection.degree() != 8:
        raise ArithmeticError("residual L(8O) intersection does not have degree eight")
    root_sum = -intersection[7]/intersection[8]
    power_sums = newton_power_sums(H)
    divisor_x_sum = sum(xA[index]*power_sums[index] for index in range(7))
    residual_x = root_sum-divisor_x_sum
    if not Bfun(residual_x):
        raise ArithmeticError("residual point has vanishing B coordinate")
    residual_y = -Afun(residual_x)/Bfun(residual_x)
    abel_trace = -E_tau(residual_x, residual_y)
    if abel_trace.is_zero():
        raise ArithmeticError("inherited P1 Abel trace is zero")
    trace_x, trace_y = map(F, abel_trace.xy())
    return trace_x, trace_y, H, relation


def interpolate_rational(points, coordinate_index, numerator_degree, denominator_degree):
    needed = numerator_degree+denominator_degree+1
    if len(points) <= needed:
        raise ArithmeticError(
            f"need more than {needed} good fibres for interpolation and holdout"
        )
    training = points[:needed]
    rows = []
    for row in training:
        base, value = row[0], row[coordinate_index]
        rows.append(
            [base**index for index in range(numerator_degree+1)]
            + [-value*base**index for index in range(denominator_degree+1)]
        )
    kernel = matrix(F, rows).right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError(
            f"rational interpolation kernel dimension {kernel.nrows()}"
        )
    relation = kernel[0]
    coefficients = list(relation)
    numerator = R(coefficients[:numerator_degree+1])
    denominator = R(coefficients[numerator_degree+1:])
    if denominator.degree() != denominator_degree:
        raise ArithmeticError("interpolated denominator dropped degree")
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    for base, xvalue, yvalue in points:
        value = (xvalue, yvalue)[coordinate_index-1]
        if not denominator(base) or numerator(base) != value*denominator(base):
            raise ArithmeticError("rational interpolation failed a holdout fibre")
    return numerator, denominator, needed, len(points)-needed


def ff_record(numerator, denominator):
    return {
        "numerator_coefficients_low_to_high": [int(value) for value in numerator.list()],
        "denominator_coefficients_low_to_high": [int(value) for value in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


if args.interpolate:
    samples = []
    rejected = []
    for value in range(int(p)):
        try:
            trace_x, trace_y, unused_H, unused_relation = abel_trace_at(F(value))
            samples.append((F(value), trace_x, trace_y))
            if args.good_fibre_limit is not None and len(samples) >= args.good_fibre_limit:
                break
        except (ArithmeticError, ZeroDivisionError, ValueError) as error:
            rejected.append({"tau": value, "reason": str(error)})
    xnum, xden, xtrain, xholdout = interpolate_rational(samples, 1, 32, 28)
    ynum, yden, ytrain, yholdout = interpolate_rational(samples, 2, 48, 42)
    trace_x_function = K(xnum)/K(xden)
    trace_y_function = K(ynum)/K(yden)
    if trace_y_function**2 != (
        trace_x_function**3+K(A3)*trace_x_function+K(B3)
    ):
        raise ArithmeticError("interpolated Abel trace misses the q4/o164 equation")
    log(
        "INTERPOLATE", good=len(samples), rejected=len(rejected),
        x_degrees="32/28", y_degrees="48/42",
        holdouts=f"{xholdout}/{yholdout}",
    )
    payload = {
        "schema": "elkies-k3.h3-q4o164-inherited-p1-abel-trace-section-modp.v1",
        "status": "PASS_MODP_Q4O164_INHERITED_P1_DEGREE7_ABEL_TRACE_SECTION",
        "prime": int(p),
        "degree_after_three_q4_maps": int(degree_u3),
        "good_fibre_count": len(samples),
        "rejected_fibres": rejected,
        "trace_section": {
            "x": ff_record(xnum, xden),
            "y": ff_record(ynum, yden),
            "exact_modp_weierstrass_identity": True,
            "training_fibres": {"x": xtrain, "y": ytrain},
            "holdout_fibres": {"x": xholdout, "y": yholdout},
        },
        "method": {
            "large_Groebner_required": False,
            "fibrewise_kernel_shape": [7, 8],
            "interpolation_degrees_from_lattice_profile": {
                "x": [32, 28], "y": [48, 42], "P_dot_O": 14,
            },
            "runtime_seconds": time.monotonic()-started,
            "good_fibre_limit": args.good_fibre_limit,
        },
        "proof_boundary": (
            "This is an exact section over GF(p)(t), reconstructed from independent "
            "fibrewise degree-seven Abel reductions and verified on the q4/orbit164 "
            "equation. It is not a characteristic-zero lift and does not by itself "
            "identify the q8/orbit376 horizontal class."
        ),
    }
else:
    tau = F(args.tau)
    trace_x, trace_y, H, relation = abel_trace_at(tau)
    log("ABEL_TRACE", tau=int(tau), x=int(trace_x), y=int(trace_y), kernel="7x8")
    payload = {
        "schema": "elkies-k3.h3-q4o164-inherited-p1-abel-trace-smoke-modp.v1",
        "status": "PASS_MODP_Q4O164_INHERITED_P1_DEGREE7_ABEL_TRACE_SMOKE",
        "prime": int(p),
        "tau": int(tau),
        "degree_after_three_q4_maps": int(degree_u3),
        "specialized_fibre_polynomial_coefficients_low_to_high": [int(value) for value in H.list()],
        "L8O_kernel": [int(value) for value in relation],
        "abel_trace": {"x": int(trace_x), "y": int(trace_y)},
        "method": {
            "large_Groebner_required": False,
            "kernel_shape": [7, 8],
            "runtime_seconds": time.monotonic()-started,
        },
        "proof_boundary": (
            "This is one exact finite-field fibre of the inherited-P1 degree-seven "
            "Abel reduction. It does not interpolate the modular trace section, identify "
            "the q8/orbit376 horizontal, or lift either object to characteristic zero."
        ),
    }
input_paths = (
    P146, P1229, Q208, Q208_MARKING, Q1584, Q1584_MARKING, Q164, Q164_MARKING,
)
payload["inputs"] = {
    "paths": [str(path.relative_to(ROOT)) for path in input_paths],
    "sha256": {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in input_paths
    },
}
default_name = (
    f"q4o164-inherited-p1-abel-trace-section-mod{int(p)}.json"
    if args.interpolate else
    f"q4o164-inherited-p1-abel-trace-smoke-mod{int(p)}.json"
)
output_arg = args.output or (LOCAL / default_name)
output = output_arg if output_arg.is_absolute() else ROOT / output_arg
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("DONE", status=payload["status"], output=output)
