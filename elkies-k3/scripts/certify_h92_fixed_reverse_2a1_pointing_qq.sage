#!/usr/bin/env sage
"""Point and mark the fixed reverse-compiled 2A1 equation over QQ.

The prescribed 2A1 zero is the nonidentity component of the old A1 fibre.
Use its split-nodal exceptional conic to select the exact sign of the quartic
point at the old I2 support.  The two transported effective 2A1 roots are
exact sections on the A1 source; restricting them to the pencil contracts
them to the two roots of the child I2-support polynomial.
"""

import argparse
import hashlib
import json
import time
from math import factorial
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--surface", type=Path, default=LOCAL / "fixed-final-a1-reverse-rr-qq.json")
parser.add_argument("--curves", type=Path, default=LOCAL / "fixed-reverse-2a1-horizontal-from-a1-qq.json")
parser.add_argument("--rr", type=Path, default=LOCAL / "fixed-reverse-2a1-rr-qq.json")
parser.add_argument("--output", type=Path, default=LOCAL / "fixed-reverse-2a1-pointing-qq.json")
parser.add_argument("--surface-status", default="PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_RR_JACOBIAN")
parser.add_argument("--curves-status", default="PASS_EXACT_QQ_FIXED_REVERSE_2A1_HORIZONTAL_ON_A1")
parser.add_argument("--rr-status", default="PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN")
parser.add_argument("--result-status", default="PASS_EXACT_QQ_FIXED_REVERSE_2A1_POINTING")
parser.add_argument("--edge", default="2A1/MW15 --q4 orbit981--> A1/MW16")
parser.add_argument("--old-i2-support")
parser.add_argument("--horizontal-roots-key", default="effective_2A1_roots_on_A1_source")
parser.add_argument("--horizontal-root-class-key", default="class_in_A1_coordinates")
parser.add_argument("--expected-child-i2-count", type=int, default=2)
parser.add_argument("--remaining-vertical-root-count", type=int, default=0)
parser.add_argument(
    "--fixed-zero-source",
    default="effective nonidentity component of the old A1 I2 fibre",
)
args = parser.parse_args()
SURFACE = args.surface.resolve()
CURVES = args.curves.resolve()
RR = args.rr.resolve()
OUTPUT = args.output.resolve()


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coeffs(poly):
    return [str(value) for value in poly.list()]


def rational_record(value):
    return {
        "numerator_coefficients_low_to_high": coeffs(value.numerator()),
        "denominator_coefficients_low_to_high": coeffs(value.denominator()),
        "degrees_numerator_denominator": [int(value.numerator().degree()), int(value.denominator().degree())],
    }


started = time.monotonic()
surface = read_json(SURFACE)
curves = read_json(CURVES)
rr = read_json(RR)
assert surface["status"] == args.surface_status
assert curves["status"] == args.curves_status
assert rr["status"] == args.rr_status

R = PolynomialRing(QQ, "s")
s = R.gen()
K = R.fraction_field()
T = PolynomialRing(QQ, "t")
t = T.gen()
KT = T.fraction_field()
U = PolynomialRing(T, "s")
ss = U.gen()
A = R(surface["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface["child"]["minimal_B_coefficients_low_to_high"])
target = curves["section"]
X = R(target["x_numerator_coefficients_low_to_high"])
Y = R(target["y_numerator_coefficients_low_to_high"])
Z = R(target["Z_coefficients_low_to_high"])
Hx = K(X) / K(Z ** 2)
Hy = K(Y) / K(Z ** 3)
basis = rr["smooth_RR"]["basis_pairs"]
AA0 = R(basis[0]["AA_coefficients_low_to_high"])
BB0 = R(basis[0]["BB_coefficients_low_to_high"])
AA1 = R(basis[1]["AA_coefficients_low_to_high"])
BB1 = R(basis[1]["BB_coefficients_low_to_high"])
quartic = sum(
    (T(values) * ss ** degree for degree, values in enumerate(rr["binary_quartic"]["coefficients_in_old_u_low_to_high"])),
    U.zero(),
)
square_factor = sum(
    (T(values) * ss ** degree for degree, values in enumerate(rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"])),
    U.zero(),
)


def newton_sqrt(value, initial, precision=14):
    PS = value.parent()
    answer = PS(initial)
    for unused in range(6):
        answer = (answer + value / answer) / 2
        answer = PS(answer.add_bigoh(precision))
    assert (answer ** 2 - value).valuation() >= precision - 2
    return answer


if args.old_i2_support is None:
    old_i2 = R(surface["child"]["finite_fibres"][0]["factor_coefficients_low_to_high"])
    assert old_i2.degree() == 1
    support = -old_i2[0] / old_i2[1]
else:
    support = QQ(args.old_i2_support)
# At the old I2 support the quartic is a polynomial square over QQ(t).
# Recover the square root coefficient-by-coefficient, avoiding factorization.
quartic_special = U(quartic)(U(support))
quartic_special = T(quartic_special)


def polynomial_square_root(poly):
    poly = T(poly)
    assert poly.degree() % 2 == 0
    degree = poly.degree() // 2
    leading = QQ(poly.leading_coefficient())
    assert leading.is_square()
    root_coefficients = [QQ(0)] * (degree + 1)
    root_coefficients[degree] = leading.sqrt()
    for index in range(degree - 1, -1, -1):
        target_degree = degree + index
        known = sum(
            root_coefficients[left] * root_coefficients[target_degree - left]
            for left in range(index + 1, degree + 1)
            if 0 <= target_degree - left <= degree
        )
        root_coefficients[index] = (poly[target_degree] - known) / (2 * root_coefficients[degree])
    root = T(root_coefficients)
    assert root ** 2 == poly
    return root


positive_ordinate = KT(polynomial_square_root(quartic_special))


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


def modular_sign(prime):
    field = GF(prime)
    RF = PolynomialRing(field, "s")
    TF = PolynomialRing(field, "t")
    UF = PolynomialRing(TF, "s")
    PS = PowerSeriesRing(field, "h", default_prec=14)
    h = PS.gen()

    def reduce_r(poly):
        return RF([reduce_qq(value, field) for value in R(poly).list()])

    def reduce_t(poly):
        return TF([reduce_qq(value, field) for value in T(poly).list()])

    try:
        support_f = reduce_qq(support, field)
        Af, Bf, Xf, Yf, Zf = map(reduce_r, (A, B, X, Y, Z))
        AA0f, BB0f, AA1f, BB1f = map(reduce_r, (AA0, BB0, AA1, BB1))
        quartic_f = UF([reduce_t(value) for value in U(quartic).list()])
        square_f = UF([reduce_t(value) for value in U(square_factor).list()])
        positive_f = reduce_t(T(positive_ordinate))
    except ZeroDivisionError:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=denominator", flush=True)
        return None

    def local(poly):
        answer = PS.zero()
        for coefficient in reversed(RF(poly).list()):
            answer = answer * (PS(support_f) + h) + PS(coefficient)
        return answer

    Aloc, Bloc = local(Af), local(Bf)
    if not Aloc[0]:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=A0", flush=True)
        return None
    node = -3 * Bloc[0] / (2 * Aloc[0])
    if node ** 2 != -Aloc[0] / 3:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=node", flush=True)
        return None
    center = newton_sqrt(-Aloc / 3, node)
    nodal_error = center ** 3 + Aloc * center + Bloc
    if nodal_error.valuation() != 2:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=nodal_val_{nodal_error.valuation()}", flush=True)
        return None
    conic_constant = nodal_error[2]
    conic_point = None
    for x1 in field:
        square = 3 * node * x1 ** 2 + conic_constant
        if square.is_square():
            conic_point = (x1, square.sqrt())
            break
    if conic_point is None:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=conic", flush=True)
        return None
    x1, y1 = conic_point
    component_x = center + h * PS(x1)
    component_y = h * newton_sqrt((component_x ** 3 + Aloc * component_x + Bloc) / h ** 2, y1)
    Zloc = local(Zf)
    Hxloc = local(Xf) / Zloc ** 2
    Hyloc = local(Yf) / Zloc ** 3
    slope = (component_y + Hyloc) / (component_x - Hxloc)
    r0 = local(AA0f) + local(BB0f) * Zloc * slope
    r1 = local(AA1f) + local(BB1f) * Zloc * slope
    if not r1:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=r1", flush=True)
        return None
    base = -r0 / r1
    if base.valuation() != 0:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=base_val_{base.valuation()}", flush=True)
        return None

    def eval_t(poly):
        answer = PS.zero()
        for coefficient in reversed(TF(poly).list()):
            answer = answer * base + PS(coefficient)
        return answer

    def eval_bivariate(poly):
        answer = PS.zero()
        for coefficient in reversed(UF(poly).list()):
            answer = answer * (PS(support_f) + h) + eval_t(coefficient)
        return answer

    bb_value = local(BB0f) + base * local(BB1f)
    W = bb_value ** 2 * (2 * component_x + Hxloc - slope ** 2) / eval_bivariate(square_f)
    identity_valuation = (W ** 2 - eval_bivariate(quartic_f)).valuation()
    if identity_valuation < 10:
        print(f"FIXEDREVERSE2A1SIGN|p={prime}|reject=quartic_val_{identity_valuation}", flush=True)
        return None
    positive = eval_t(positive_f)
    plus_valuation = (W - positive).valuation()
    minus_valuation = (W + positive).valuation()
    print(f"FIXEDREVERSE2A1SIGN|p={prime}|quartic_val={identity_valuation}|plus={plus_valuation}|minus={minus_valuation}", flush=True)
    # The arc is transverse to the exceptional component, so only its h=0
    # value lies on the component; valuation one is the exact sign gate.
    if plus_valuation >= 1:
        return 1, [int(x1), int(y1)]
    if minus_valuation >= 1:
        return -1, [int(x1), int(y1)]
    return None


sign_certificate = None
for sign_prime in (131, 137, 139, 149, 151, 157, 163, 167, 173, 179):
    sign_certificate = modular_sign(sign_prime)
    if sign_certificate is not None:
        break
if sign_certificate is None:
    raise ArithmeticError("no good small prime selected the exceptional-component sign")
quartic_sign, modular_conic_point = sign_certificate
ordinate = positive_ordinate if quartic_sign == 1 else -positive_ordinate


def eval_u_at_support(poly):
    return KT(U(poly)(U(support)))


translated = [eval_u_at_support(quartic.derivative(order)) / factorial(order) for order in range(5)]
ee, dd, cc, bb, aa = translated
assert ee == ordinate ** 2
a1 = dd / ordinate
a2 = cc - dd ** 2 / (4 * ordinate ** 2)
a3 = 2 * ordinate * bb
a4 = -4 * ordinate ** 2 * aa
a6 = a2 * a4
b2 = a1 ** 2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3 ** 2 + 4 * a6
c4 = b2 ** 2 - 24 * b4
c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
A_pointed = -c4 / 48
B_pointed = -c6 / 864
A_child = T(rr["child"]["minimal_A_coefficients_low_to_high"])
B_child = T(rr["child"]["minimal_B_coefficients_low_to_high"])
assert 81 * A_pointed == KT(A_child)
assert 729 * B_pointed == KT(B_child)

# Contract the two exact effective root sections to the two child I2 supports.
child_i2 = T(rr["child"]["finite_fibres"][0]["factor_coefficients_low_to_high"])
assert child_i2.degree() == args.expected_child_i2_count and child_i2.is_squarefree()


def section_from_record(record):
    x = K(R(record["x_numerator_coefficients_low_to_high"])) / K(R(record["x_denominator_coefficients_low_to_high"]))
    y = K(R(record["y_numerator_coefficients_low_to_high"])) / K(R(record["y_denominator_coefficients_low_to_high"]))
    assert y ** 2 == x ** 3 + K(A) * x + K(B)
    return x, y


def restrict_section(P):
    px, py = P
    m = (py + Hy) / (px - Hx)
    base = -(K(AA0) + K(BB0 * Z) * m) / (K(AA1) + K(BB1 * Z) * m)
    assert base.numerator().degree() == base.denominator().degree() == 0
    return QQ(base)


root_supports = [
    restrict_section(section_from_record(row["section"]))
    for row in curves[args.horizontal_roots_key]
]
assert len(set(root_supports)) == len(root_supports)
assert all(child_i2(value) == 0 for value in root_supports)
assert len(root_supports) + args.remaining_vertical_root_count == args.expected_child_i2_count

remaining_factor = child_i2
for root_support in root_supports:
    quotient, remainder = remaining_factor.quo_rem(t - root_support)
    assert not remainder
    remaining_factor = quotient
assert remaining_factor.degree() == args.remaining_vertical_root_count

payload = {
    "schema": "elkies-k3.fixed-reverse-semistable-pointing-qq.v1",
    "status": args.result_status,
    "fixed_zero": {
        "source_curve": args.fixed_zero_source,
        "old_I2_support": str(support),
        "sign_selection_prime": int(sign_prime),
        "exceptional_conic_point_mod_p": modular_conic_point,
        "quartic_sign": quartic_sign,
        "quartic_ordinate": rational_record(ordinate),
        "degree_one_pointed_quartic_map": True,
        "short_scaling": "x_short=9*x_pointed, y_short=27*y_pointed",
        "exact_A_identity": "81*A_pointed=A_child",
        "exact_B_identity": "729*B_pointed=B_child",
    },
    "effective_horizontal_components": [
        {
            "source_class": row[args.horizontal_root_class_key],
            "child_I2_support": str(root_support),
            "new_base_degree": 0,
            "role": "nonidentity component relative to the prescribed fixed zero",
        }
        for row, root_support in zip(curves[args.horizontal_roots_key], root_supports)
    ],
    "remaining_vertical_components": {
        "count": args.remaining_vertical_root_count,
        "child_I2_support_factor_coefficients_low_to_high": coeffs(remaining_factor),
        "identification": "remaining child I2 supports after all exact horizontal root curves are contracted",
    },
    "marking": {
        "fixed_transition": args.edge,
        "prescribed_zero_pointed": True,
        "horizontal_roots_bound_to_distinct_equation_I2_fibres": True,
        "pinned_R17_transport_inherited": True,
    },
    "method": {
        "split_nodal_exceptional_conic": True,
        "exact_curve_restrictions": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SURFACE, CURVES, RR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SURFACE, CURVES, RR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSEPOINT|edge={}|old_I2_component_zero=1|sign={}|A=1|B=1|"
    "horizontal_child_I2_supports={}|remaining={}|seconds={:.3f}|status={}|output={}".format(
        args.edge, quartic_sign, ",".join(map(str, root_supports)), args.remaining_vertical_root_count,
        payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
