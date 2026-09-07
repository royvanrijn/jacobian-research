#!/usr/bin/env sage -python
"""Construct the pinned A11 zero from the exact low-trace word.

On the A11 equation (whose stored zero is C10), the pinned zero is

    O_pinned = S5 - 2*S7 - S17 - Qminus.

All four sections are already exact and marked.  This script performs only
univariate QQ(T) elliptic-curve arithmetic, verifies the resulting section
identity and I12 profile, and compares its reduction with the independently
reconstructed pinned modular zero.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, PowerSeriesRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--marked-reduction",
    type=Path,
    default=LOCAL / "q24-a11-pinned-zero-section-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-pinned-zero-section-qq.json",
)
args = parser.parse_args()

A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
S5 = LOCAL / "q24-a11-s5-trace-section-qq.json"
DEGREE_ONE = LOCAL / "q24-a11-degree-one-shell-sections-qq.json"
COVARIANTS = LOCAL / "q24-a11-degree-one-shell-covariants-qq.json"
QPOINT = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
ROUTE = GENERATED / "elkies-k3-h3-q24-a11-q8-zero-translation-route.json"
MARKED = args.marked_reduction.resolve()
INPUTS = (A11, S5, DEGREE_ONE, COVARIANTS, QPOINT, ROUTE, MARKED)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

a11 = json.loads(A11.read_text())
s5 = json.loads(S5.read_text())
degree_one = json.loads(DEGREE_ONE.read_text())
covariants = json.loads(COVARIANTS.read_text())
qpoint = json.loads(QPOINT.read_text())
route = json.loads(ROUTE.read_text())
marked = json.loads(MARKED.read_text())
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert s5["status"] == "PASS_EXACT_MARKED_A11_S5_TRACE_SECTION_QQ"
assert degree_one["status"] == "PASS_EXACT_MARKED_A11_DEGREE_ONE_SHELL_SECTIONS_QQ"
assert covariants["status"] == "PASS_EXACT_A11_DEGREE_ONE_CURVE_COVARIANTS_QQ"
assert qpoint["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"
assert route["status"] == "PASS_EXACT_A11_Q8_ZERO_TRANSLATION_ROUTE"
assert marked["status"] == "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in a11["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in a11["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])


def section_point(record):
    X = R([QQ(value) for value in record["X_coefficients_low_to_high"]])
    Y = R([QQ(value) for value in record["Y_coefficients_low_to_high"]])
    Z = R([QQ(value) for value in record["Z_coefficients_low_to_high"]])
    return E(K(X) / K(Z**2), K(Y) / K(Z**3))


def rational_record(record):
    numerator = R([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = R(
        [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    )
    return K(numerator) / K(denominator)


def power_root(poly, exponent):
    poly = R(poly)
    leading = QQ(poly.leading_coefficient())
    if exponent % 2 == 0 and not leading.is_square():
        raise ArithmeticError("denominator leading coefficient is not a square")
    answer = R.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


def normalized_section(point):
    x, y = map(K, point.xy())
    Zx = power_root(x.denominator(), 2)
    Zy = power_root(y.denominator(), 3)
    if Zx != Zy:
        raise ArithmeticError("x and y denominators give different section Z")
    Z = Zx
    X = R(x * Z**2)
    Y = R(y * Z**3)
    if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
        raise ArithmeticError("normalized section missed A11")
    return X, Y, Z


degree_one_by_index = {
    int(row["equation_shell_index"]): row["section"]
    for row in degree_one["sections"]
}
S5point = section_point(s5["section"])
S17point = section_point(degree_one_by_index[17])
C7row = next(
    row for row in covariants["points"]
    if row["shell_kind"] == "identity" and int(row["equation_shell_index"]) == 7
)
C7 = E(
    rational_record(C7row["covariant_point"]["x"]),
    rational_record(C7row["covariant_point"]["y"]),
)

# The fixed-pointed degree-one construction gives C7=-2*S7-Qminus exactly.
# Using the already constructed covariant avoids a costly redundant doubling.
partial = S5point + C7
print("A11PINNEDZEROQQ|stage=S5_plus_C7|status=PASS", flush=True)
Opinned = partial - S17point
print("A11PINNEDZEROQQ|stage=minus_S17|status=PASS", flush=True)
section = normalized_section(Opinned)
X, Y, Z = section
if [int(X.degree()), int(Y.degree()), int(Z.degree())] != [10, 15, 3]:
    raise ArithmeticError("pinned zero has the wrong exact degrees")

# Exact formal I12 profile.
i12_factor = R(
    next(
        row["factor"]
        for row in a11["child"]["discriminant_factorization"]
        if int(row["multiplicity"]) == 12
    )
)
if i12_factor.degree() != 1:
    raise ArithmeticError("A11 I12 factor is not linear")
beta = -QQ(i12_factor[0]) / QQ(i12_factor[1])
PS = PowerSeriesRing(QQ, "s", default_prec=15)
s = PS.gen()
A_series = PS(A(T + beta))
B_series = PS(B(T + beta))
center = PS(-3 * B_series[0] / (2 * A_series[0]))
for unused in range(7):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 14:
    raise ArithmeticError("I12 formal center did not converge")


def component_depth(point):
    x, y = map(K, point.xy())
    xs = PS(x(T + beta))
    ys = PS(y(T + beta))
    return min(int((xs - center).valuation()), int(ys.valuation()), 6)


zero_depth = component_depth(Opinned)
if zero_depth != 2:
    raise ArithmeticError("exact pinned-zero component marking failed")

# Full exact reduction comparison at the pinned good prime.
p = ZZ(marked["prime"])
F = GF(p)
RF = PolynomialRing(F, "T")
KF = RF.fraction_field()
AF = RF([F(value.numerator()) / F(value.denominator()) for value in A.list()])
BF = RF([F(value.numerator()) / F(value.denominator()) for value in B.list()])
EF = EllipticCurve(KF, [0, 0, 0, KF(AF), KF(BF)])


def reduce_exact(poly):
    return RF([F(value.numerator()) / F(value.denominator()) for value in poly.list()])


Xp, Yp, Zp = map(reduce_exact, section)
reduced = EF(KF(Xp) / KF(Zp**2), KF(Yp) / KF(Zp**3))
target = marked["section"]
Xm = RF(target["X_coefficients_low_to_high"])
Ym = RF(target["Y_coefficients_low_to_high"])
Zm = RF(target["Z_coefficients_low_to_high"])
if reduced != EF(KF(Xm) / KF(Zm**2), KF(Ym) / KF(Zm**3)):
    raise ArithmeticError("exact pinned zero missed the marked modular section")


def record(value):
    XX, YY, ZZpoly = value
    return {
        "X_coefficients_low_to_high": [str(item) for item in XX.list()],
        "Y_coefficients_low_to_high": [str(item) for item in YY.list()],
        "Z_coefficients_low_to_high": [str(item) for item in ZZpoly.list()],
        "degrees_X_Y_Z": [int(XX.degree()), int(YY.degree()), int(ZZpoly.degree())],
        "exact_weierstrass_identity": True,
    }


payload = {
    "schema": "elkies-k3.h3-q24-a11-pinned-zero-section-qq.v1",
    "status": "PASS_EXACT_MARKED_A11_PINNED_ZERO_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "lattice_word": "S5 - 2*S7 - S17 - Qminus",
    "section": record(section),
    "marking": {
        "P_dot_equation_zero": int(Z.degree()),
        "I12_component_depth_up_to_negation": int(zero_depth),
        "marked_reduction_prime": int(p),
        "exact_reduction_matches_marked_zero": True,
    },
    "method": (
        "exact A11 group law on S5+C7-S17, using the fixed-pointed covariant "
        "identity C7=-2*S7-Qminus"
    ),
    "large_Groebner_required": False,
    "proof_boundary": (
        "This constructs the exact pinned zero on A11. The second q8 horizontal "
        "section and the two-dimensional H0 basis remain to be constructed exactly."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11PINNEDZEROQQ|degrees={}|PO={}|I12depth={}|prime={}|status={}".format(
        tuple(payload["section"]["degrees_X_Y_Z"]),
        Z.degree(),
        zero_depth,
        p,
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
