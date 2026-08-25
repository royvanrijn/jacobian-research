#!/usr/bin/env sage -python
"""Sample selected low-degree D12-to-A11 traces over one good finite field.

The default ``bridge-m`` curve set restricts the certified orbit64 RR pencil
to E5 and identity-shell curves S7, S14, S17.  The corrected bridge set uses
E5, S0 and S1 after the pinned-zero audit reversed the E5 trace marking.  The
``q8-zero-target`` set uses
S2, S5, S6, S7, S8, S13, S17, which are precisely the traces needed for the C10
zero correction and the q8 target word.  Both quartic branches are mapped
through the global pointed-quartic isomorphism.  Both signs of the exact
pointed-opposite section are retained; q8 mode also records both signs of the
new exact bridge M.

Only univariate finite-field arithmetic and small linear kernels are used.
This is construction data, not a characteristic-zero arrow certificate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, is_prime, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--start-tau", type=int, default=2)
parser.add_argument("--samples", type=int, default=80)
parser.add_argument(
    "--curve-set",
    choices=("bridge-m", "bridge-m-corrected", "q8-zero-target"),
    default="bridge-m",
)
parser.add_argument(
    "--output",
    type=Path,
    default=None,
)
args = parser.parse_args()

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")

A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
CANDIDATES = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
E5 = LOCAL / "q24-d12-missing-e5-section-qq.json"
POINTED = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
BRIDGE_M = LOCAL / "q24-a11-bridge-m-section-marked-qq.json"
if args.curve_set in ("bridge-m", "bridge-m-corrected"):
    INPUTS = (A11, PARENT, CANDIDATES, ZERO, MODEL, E5, POINTED)
    output = args.output or LOCAL / f"q24-a11-e5-bridge-traces-mod{p}.json"
else:
    INPUTS = (A11, PARENT, CANDIDATES, ZERO, MODEL, POINTED, BRIDGE_M)
    output = args.output or LOCAL / f"q24-a11-q8-zero-target-traces-mod{p}.json"
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

a11 = json.loads(A11.read_text())
parent = json.loads(PARENT.read_text())
candidates = json.loads(CANDIDATES.read_text())
zero = json.loads(ZERO.read_text())
model = json.loads(MODEL.read_text())["exact_model"]
e5 = json.loads(E5.read_text())
pointed = json.loads(POINTED.read_text())
bridge_m = json.loads(BRIDGE_M.read_text()) if args.curve_set == "q8-zero-target" else None
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert candidates["status"] == "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
if args.curve_set in ("bridge-m", "bridge-m-corrected"):
    assert e5["status"] == "PASS_EXACT_Q24_D12_MISSING_E5_SECTION_QQ"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"
if bridge_m is not None:
    assert bridge_m["status"] == "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ"

F = GF(p)
RV = PolynomialRing(F, "V")
V = RV.gen()
KV = RV.fraction_field()
RT = PolynomialRing(F, "T")
T = RT.gen()
KT = RT.fraction_field()
UQ = PolynomialRing(QQ, "u")
KUQ = UQ.fraction_field()
VQ = PolynomialRing(QQ, "V")
KVQ = VQ.fraction_field()
TQ = PolynomialRing(QQ, "T")


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


def red_poly(values, ring):
    return ring([red(value) for value in values])


def red_rational(text, variable):
    source_ring = PolynomialRing(QQ, variable)
    source = source_ring.fraction_field()(str(text))
    target_ring = RV if variable == "V" else RT
    numerator = target_ring([red(value) for value in source.numerator().list()])
    denominator = target_ring([red(value) for value in source.denominator().list()])
    return target_ring.fraction_field()(numerator) / target_ring.fraction_field()(denominator)


def evaluate_polynomial(poly, argument):
    return sum(poly[index] * argument**index for index in range(poly.degree() + 1))


def evaluate_u_rational(value, argument):
    value = KUQ(value)
    numerator = RV([red(item) for item in value.numerator().list()])
    denominator = RV([red(item) for item in value.denominator().list()])
    return KV(evaluate_polynomial(numerator, argument)) / KV(evaluate_polynomial(denominator, argument))


def evaluate_T(value, argument):
    value = KT(value)
    return KV(evaluate_polynomial(RT(value.numerator()), argument)) / KV(
        evaluate_polynomial(RT(value.denominator()), argument)
    )


def rational_square_root(value, ring):
    value = ring.fraction_field()(value)
    if not value.is_square():
        raise ArithmeticError("restricted quartic value is not a rational-function square")
    answer = value.sqrt()
    assert answer**2 == value
    return answer


def map_degree(value, ring):
    value = ring.fraction_field()(value)
    return max(int(value.numerator().degree()), int(value.denominator().degree()))


def denominator_power_root(value, exponent):
    denominator = RV(KV(value).denominator()).monic()
    answer = RV.one()
    for factor, multiplicity in denominator.factor():
        if int(multiplicity) % exponent:
            return None
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


# -------------------------------------------------------------------------
# Rebuild the selected orbit64 pencil on the minimal D12 model.
# -------------------------------------------------------------------------
u_of_V = red_rational(a11["coordinate_change"]["u_of_V"], "V")
x_scale = red_rational(a11["coordinate_change"]["x_scale"], "V")
y_scale = red_rational(a11["coordinate_change"]["y_scale"], "V")
A12 = red_poly(parent["child"]["minimal_A_coefficients_low_to_high"], RV)
B12 = red_poly(parent["child"]["minimal_B_coefficients_low_to_high"], RV)

selected = candidates["candidates"][a11["selected_section"]["candidate_index"]]
Xu = red_poly(selected["X_coefficients_low_to_high"], RV)
Yu = red_poly(selected["Y_coefficients_low_to_high"], RV)
Zu = red_poly(selected["Z_coefficients_low_to_high"], RV)
x_parent = x_scale * evaluate_u_rational(KUQ(UQ([QQ(v) for v in selected["X_coefficients_low_to_high"]])) / KUQ(UQ([QQ(v) for v in selected["Z_coefficients_low_to_high"]])**2), u_of_V)
y_parent = y_scale * evaluate_u_rational(KUQ(UQ([QQ(v) for v in selected["Y_coefficients_low_to_high"]])) / KUQ(UQ([QQ(v) for v in selected["Z_coefficients_low_to_high"]])**3), u_of_V)
assert y_parent**2 == x_parent**3 + KV(A12) * x_parent + KV(B12)

denominator = RV(x_parent.denominator()).monic()
Z_parent = RV.one()
for factor, multiplicity in denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("selected parent x denominator is not a square")
    Z_parent *= factor.monic() ** (int(multiplicity) // 2)
X_parent = RV(x_parent * Z_parent**2)
Y_parent = RV(y_parent * Z_parent**3)
alpha_collision = red(model["I8star_root"])
modulus_collision = Z_parent**2
X_inverse = X_parent.inverse_mod(modulus_collision)
pairs = []
for BB in (RV.one(), V):
    AA = RV((BB * Y_parent * X_inverse) % modulus_collision)
    AA -= AA(alpha_collision) / Z_parent(alpha_collision) ** 2 * Z_parent**2
    pairs.append((KV(AA) / KV(Z_parent**2), KV(BB) / KV(Z_parent)))
(a0, b0), (a1, b1) = pairs

# A11 binary quartic in old-base V with coefficients in F_p(T).
quartic_coefficients = []
for text_value in a11["quartic"]["coefficients_in_T_low_to_high"]:
    source = TQ(str(text_value))
    quartic_coefficients.append(RT([red(value) for value in source.list()]))

A11poly = red_poly(a11["child"]["minimal_A_coefficients_low_to_high"], RT)
B11poly = red_poly(a11["child"]["minimal_B_coefficients_low_to_high"], RT)

i8 = next(row for row in parent["child"]["finite_fibres"] if row["kodaira"] == "I8*")
i8_factor = VQ(str(i8["factor"]))
alpha = red(-i8_factor[0] / i8_factor[1])
q_squared = sum(KT(poly) * alpha**index for index, poly in enumerate(quartic_coefficients))
q_anchor = rational_square_root(q_squared, RT)
e, d0, c0, b0q, a0q = [KT(poly) for poly in quartic_coefficients]
anchor_a = a0q
anchor_b = b0q + 4 * alpha * a0q
anchor_c = c0 + 3 * alpha * b0q + 6 * alpha**2 * a0q
anchor_d = d0 + 2 * alpha * c0 + 3 * alpha**2 * b0q + 4 * alpha**3 * a0q


def quartic_value(Tmap):
    return sum(evaluate_T(poly, Tmap) * V**index for index, poly in enumerate(quartic_coefficients))


def point_to_A11(Tmap, Wmap):
    q = evaluate_T(q_anchor, Tmap)
    aa = evaluate_T(anchor_a, Tmap)
    bb = evaluate_T(anchor_b, Tmap)
    cc = evaluate_T(anchor_c, Tmap)
    dd = evaluate_T(anchor_d, Tmap)
    uu = KV(V - alpha)
    xg = (2 * q * (Wmap + q) + dd * uu) / uu**2
    yg = (
        4 * q**2 * (Wmap + q)
        + 2 * q * (dd * uu + cc * uu**2)
        - dd**2 * uu**2 / (2 * q)
    ) / uu**3
    aa1 = dd / q
    aa2 = cc - dd**2 / (4 * q**2)
    aa3 = 2 * q * bb
    b2 = aa1**2 + 4 * aa2
    x = KV(9 * (xg + b2 / 12))
    y = KV(27 * (yg + (aa1 * xg + aa3) / 2))
    assert y**2 == x**3 + evaluate_T(KT(A11poly), Tmap) * x + evaluate_T(KT(B11poly), Tmap)
    return x, y


def d12_exact_section(section):
    X = red_poly(section["X_coefficients_low_to_high"], RV)
    Y = red_poly(section["Y_coefficients_low_to_high"], RV)
    Z = red_poly(section["Z_coefficients_low_to_high"], RV)
    x = KV(X) / KV(Z**2)
    y = KV(Y) / KV(Z**3)
    assert y**2 == x**3 + KV(A12) * x + KV(B12)
    return x, y


def d12_shell_section(index):
    row = zero["sections"][index]
    x_u = UQ([QQ(value) for value in row["x_coefficients_low_to_high"]])
    y_u = UQ([QQ(value) for value in row["y_coefficients_low_to_high"]])
    x = x_scale * evaluate_u_rational(KUQ(x_u), u_of_V)
    y = y_scale * evaluate_u_rational(KUQ(y_u), u_of_V)
    assert y**2 == x**3 + KV(A12) * x + KV(B12)
    return x, y


def prepare_curve(name, x_curve, y_curve, expected_degree):
    chord = (y_curve + x_parent.parent()(y_parent)) / (x_curve - x_parent)
    Tmap = KV((a1 + b1 * chord) / (a0 + b0 * chord))
    if map_degree(Tmap, RV) != expected_degree:
        raise ArithmeticError(f"{name}: new-base degree {map_degree(Tmap, RV)}, expected {expected_degree}")
    Wplus = rational_square_root(quartic_value(Tmap), RV)
    branches = []
    for sign, Wmap in ((1, Wplus), (-1, -Wplus)):
        x11, y11 = point_to_A11(Tmap, Wmap)
        pole_order = None
        if expected_degree == 1:
            Zx = denominator_power_root(x11, 2)
            Zy = denominator_power_root(y11, 3)
            if Zx is not None and Zy is not None and Zx == Zy:
                pole_order = int(Zx.degree())
        branches.append((sign, x11, y11, pole_order))
    return {"name": name, "degree": expected_degree, "Tmap": Tmap, "branches": branches}


if args.curve_set == "bridge-m":
    curves = [
        prepare_curve("E5", *d12_exact_section(e5["section"]), 13),
        prepare_curve("S7", *d12_shell_section(7), 1),
        prepare_curve("S14", *d12_shell_section(14), 3),
        prepare_curve("S17", *d12_shell_section(17), 1),
    ]
elif args.curve_set == "bridge-m-corrected":
    curves = [
        prepare_curve("E5", *d12_exact_section(e5["section"]), 13),
        prepare_curve("S0", *d12_shell_section(0), 6),
        prepare_curve("S1", *d12_shell_section(1), 4),
        prepare_curve("S7", *d12_shell_section(7), 1),
    ]
else:
    curves = [
        prepare_curve(f"S{index}", *d12_shell_section(index), degree)
        for index, degree in (
            (2, 2),
            (5, 3),
            (6, 9),
            (7, 1),
            (8, 2),
            (13, 3),
            (17, 1),
        )
    ]


def reduce_modulus(value, modulus):
    value = KV(value)
    numerator = RV(value.numerator())
    denominator = RV(value.denominator())
    if denominator.gcd(modulus).degree() != 0:
        raise ZeroDivisionError("curve denominator meets trace fibre")
    return RV((numerator * denominator.inverse_mod(modulus)) % modulus)


def newton_power_sums(monic):
    degree = monic.degree()
    sums = [F(degree)]
    for exponent in range(1, degree):
        total = F(exponent) * monic[degree - exponent]
        for index in range(1, exponent):
            total += monic[degree - index] * sums[exponent - index]
        sums.append(-total)
    return sums


def trace_at(curve_record, branch_record, tau_integer):
    degree = curve_record["degree"]
    sign, xmap, ymap, unused_pole_order = branch_record
    tau = F(tau_integer)
    Tmap = curve_record["Tmap"]
    equation = RV(Tmap.numerator() - tau * Tmap.denominator())
    if equation.gcd(equation.derivative()).degree() != 0:
        raise ArithmeticError("trace fibre is not etale")
    modulus = equation.monic()
    if modulus.degree() != degree:
        raise ArithmeticError("trace fibre has wrong degree")
    x_residue = reduce_modulus(xmap, modulus)
    y_residue = reduce_modulus(ymap, modulus)
    Avalue, Bvalue = A11poly(tau), B11poly(tau)
    if (y_residue**2 - x_residue**3 - Avalue * x_residue - Bvalue) % modulus:
        raise ArithmeticError("A11 multisection identity failed modulo trace fibre")
    curve = EllipticCurve(F, [0, 0, 0, Avalue, Bvalue])
    if degree == 1:
        point = curve(x_residue[0], y_residue[0])
    else:
        trace_order = degree + 1
        x_count = trace_order // 2 + 1
        y_count = (trace_order - 3) // 2 + 1
        powers = [RV.one()]
        for unused in range(max(x_count, y_count) - 1):
            powers.append((powers[-1] * x_residue) % modulus)
        columns = powers[:x_count] + [(y_residue * powers[index]) % modulus for index in range(y_count)]
        assert len(columns) == trace_order
        evaluation = matrix(F, degree, trace_order, lambda row, column: columns[column][row])
        kernel = evaluation.right_kernel().basis_matrix()
        if kernel.nrows() != 1:
            raise ArithmeticError("trace kernel is not one-dimensional")
        relation = kernel[0]
        RX = PolynomialRing(F, "Xtrace")
        XX = RX.gen()
        afun = sum(relation[index] * XX**index for index in range(x_count))
        bfun = sum(relation[x_count + index] * XX**index for index in range(y_count))
        intersection = afun**2 - (XX**3 + Avalue * XX + Bvalue) * bfun**2
        if intersection.degree() != trace_order:
            raise ArithmeticError("trace residual intersection has wrong degree")
        residual_root_sum = -intersection[trace_order - 1] / intersection[trace_order]
        power_sums = newton_power_sums(modulus)
        trace_x = sum(x_residue[index] * power_sums[index] for index in range(degree))
        residual_x = residual_root_sum - trace_x
        residual_b = bfun(residual_x)
        if not residual_b:
            raise ArithmeticError("trace residual B value vanished")
        residual_y = -afun(residual_x) / residual_b
        point = -curve(residual_x, residual_y)
    return {
        "curve": curve_record["name"],
        "quartic_branch_sign": sign,
        "tau": int(tau),
        "degree": degree,
        "trace_space": f"L({degree + 1}O)" if degree > 1 else "direct degree-one point",
        "AJ_x": int(point[0]),
        "AJ_y": int(point[1]),
    }


exact_sections = {"Q": pointed["section"]}
if bridge_m is not None:
    exact_sections["M"] = bridge_m["section"]
exact_section_polynomials = {
    name: (
        red_poly(section["X_coefficients_low_to_high"], RT),
        red_poly(section["Y_coefficients_low_to_high"], RT),
        red_poly(section["Z_coefficients_low_to_high"], RT),
    )
    for name, section in exact_sections.items()
}

samples = []
tau_integer = args.start_tau
while len(samples) < args.samples:
    if tau_integer >= p:
        raise ArithmeticError("ran out of finite-field specializations")
    tau = F(tau_integer)
    try:
        row = {"tau": int(tau), "points": []}
        for curve_record in curves:
            for branch_record in curve_record["branches"]:
                row["points"].append(trace_at(curve_record, branch_record, tau_integer))
        EC = EllipticCurve(F, [0, 0, 0, A11poly(tau), B11poly(tau)])
        for name, (XX, YY, ZZpoly) in exact_section_polynomials.items():
            if not ZZpoly(tau):
                raise ZeroDivisionError(f"exact section {name} has a pole at sample")
            point = EC(XX(tau) / ZZpoly(tau) ** 2, YY(tau) / ZZpoly(tau) ** 3)
            for sign, signed_point in ((1, point), (-1, -point)):
                row["points"].append(
                    {
                        "curve": name,
                        "quartic_branch_sign": sign,
                        "tau": int(tau),
                        "degree": 0,
                        "trace_space": "exact section",
                        "AJ_x": int(signed_point[0]),
                        "AJ_y": int(signed_point[1]),
                    }
                )
    except (ArithmeticError, ZeroDivisionError, ValueError):
        tau_integer += 1
        continue
    samples.append(row)
    tau_integer += 1

if args.curve_set == "bridge-m":
    schema = "elkies-k3.h3-q24-a11-e5-bridge-traces-modp.v1"
    status = "PASS_Q24_A11_E5_BRIDGE_TRACE_SAMPLES_MODP"
elif args.curve_set == "bridge-m-corrected":
    schema = "elkies-k3.h3-q24-a11-corrected-bridge-traces-modp.v1"
    status = "PASS_Q24_A11_CORRECTED_BRIDGE_TRACE_SAMPLES_MODP"
else:
    schema = "elkies-k3.h3-q24-a11-q8-zero-target-traces-modp.v1"
    status = "PASS_Q24_A11_Q8_ZERO_TARGET_TRACE_SAMPLES_MODP"

payload = {
    "schema": schema,
    "status": status,
    "curve_set": args.curve_set,
    "prime": int(p),
    "sample_count": len(samples),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "curve_degrees": {row["name"]: row["degree"] for row in curves},
    "degree_one_branch_pole_orders": {
        row["name"]: {str(branch[0]): branch[3] for branch in row["branches"]}
        for row in curves if row["degree"] == 1
    },
    "trace_spaces": {
        row["name"]: f"L({row['degree'] + 1}O)"
        for row in curves
        if row["degree"] > 1
    },
    "exact_sections": sorted(exact_sections),
    "samples": samples,
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact over one good finite field. Both quartic branches and both Q signs are retained; "
        "the subsequent univariate reconstruction must select the pinned equation orientation."
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11TRACES|set={}|prime={}|samples={}|degrees={}|branches=2|status={}".format(
        args.curve_set,
        p,
        len(samples),
        ",".join(str(row["degree"]) for row in curves),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{output.resolve()}", flush=True)
