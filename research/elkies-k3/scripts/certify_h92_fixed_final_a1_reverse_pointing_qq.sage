#!/usr/bin/env sage
"""Point the reverse-compiled fixed A1 equation at its prescribed zero.

Restrict the exact smooth-chord pencil and quartic square root to two exact
curves on the q12 rootless source:

* the historical fixed-corridor A1 zero (degree one over the new base);
* the effective A1 root component (contracted to the unique I2 support).

The degree-one quartic point gives a generalized Weierstrass model whose short
invariants agree exactly with the compiled A1 Jacobian.  No elimination or
Groebner basis is used.
"""

import hashlib
import json
import time
from math import factorial
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SOURCE = LOCAL / "q12o5867-smooth-rr-qq.json"
CURVES = LOCAL / "fixed-final-a1-horizontal-from-q12-endpoint-qq.json"
RR = LOCAL / "fixed-final-a1-reverse-rr-qq.json"
OUTPUT = LOCAL / "fixed-final-a1-reverse-pointing-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coeffs(poly):
    return [str(value) for value in poly.list()]


def rational_record(value):
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": coeffs(value.numerator()),
        "denominator_coefficients_low_to_high": coeffs(value.denominator()),
        "degrees_numerator_denominator": [int(value.numerator().degree()), int(value.denominator().degree())],
    }


def rational_bits(values):
    answer = 0
    for value in values:
        for coefficient in value.numerator().list() + value.denominator().list():
            coefficient = QQ(coefficient)
            answer = max(answer, abs(ZZ(coefficient.numerator())).nbits(), ZZ(coefficient.denominator()).nbits())
    return int(answer)


started = time.monotonic()
source = read_json(SOURCE)
curves = read_json(CURVES)
rr = read_json(RR)
assert source["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
assert curves["status"] == "PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_HORIZONTAL"
assert rr["status"] == "PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_RR_JACOBIAN"

R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()
S = PolynomialRing(QQ, "s")
s = S.gen()
KS = S.fraction_field()
U = PolynomialRing(S, "u")
uu = U.gen()

A = R(source["child"]["minimal_A_coefficients_low_to_high"])
B = R(source["child"]["minimal_B_coefficients_low_to_high"])
horizontal = curves["section"]
X = R(horizontal["x_numerator_coefficients_low_to_high"])
Y = R(horizontal["y_numerator_coefficients_low_to_high"])
Z = R(horizontal["Z_coefficients_low_to_high"])
Hx = K(X) / K(Z ** 2)
Hy = K(Y) / K(Z ** 3)

basis = rr["smooth_RR"]["basis_pairs"]
AA0 = R(basis[0]["AA_coefficients_low_to_high"])
BB0 = R(basis[0]["BB_coefficients_low_to_high"])
AA1 = R(basis[1]["AA_coefficients_low_to_high"])
BB1 = R(basis[1]["BB_coefficients_low_to_high"])

quartic = U.zero()
for degree, values in enumerate(rr["binary_quartic"]["coefficients_in_old_u_low_to_high"]):
    quartic += S(values) * uu ** degree
square_factor = U.zero()
for degree, values in enumerate(rr["binary_quartic"]["square_factor_coefficients_in_old_u_low_to_high"]):
    square_factor += S(values) * uu ** degree
assert quartic.degree() == 4


def source_point(record, polynomial=False):
    if polynomial:
        x = K(R(record["x_coefficients_low_to_high"]))
        y = K(R(record["y_coefficients_low_to_high"]))
    else:
        x = K(R(record["x_numerator_coefficients_low_to_high"])) / K(R(record["x_denominator_coefficients_low_to_high"]))
        y = K(R(record["y_numerator_coefficients_low_to_high"])) / K(R(record["y_denominator_coefficients_low_to_high"]))
    assert y ** 2 == x ** 3 + K(A) * x + K(B)
    return x, y


def new_base_and_ordinate(Q):
    qx, qy = Q
    slope = (qy + Hy) / (qx - Hx)
    restriction0 = K(AA0) + K(BB0 * Z) * slope
    restriction1 = K(AA1) + K(BB1 * Z) * slope
    new_base = -restriction0 / restriction1

    def evaluate_s_polynomial_at_k(poly):
        answer = K.zero()
        for coefficient in reversed(S(poly).list()):
            answer = answer * new_base + K(coefficient)
        return answer

    def evaluate_bivariate_at_k(poly):
        answer = K.zero()
        for coefficient in reversed(U(poly).list()):
            answer = answer * K(u) + evaluate_s_polynomial_at_k(coefficient)
        return answer

    bb = K(BB0) + new_base * K(BB1)
    square_value = evaluate_bivariate_at_k(square_factor)
    ordinate = bb ** 2 * (2 * qx + Hx - slope ** 2) / square_value
    assert ordinate ** 2 == evaluate_bivariate_at_k(quartic)
    return new_base, ordinate


fixed_zero = source_point(curves["fixed_A1_zero_on_rootless_source"])
new_base, ordinate_u = new_base_and_ordinate(fixed_zero)
new_base_num = R(new_base.numerator())
new_base_den = R(new_base.denominator())
assert max(new_base_num.degree(), new_base_den.degree()) == 1

# Invert s=(n0+n1*u)/(d0+d1*u) exactly.
n0, n1 = new_base_num[0], new_base_num[1]
d0, d1 = new_base_den[0], new_base_den[1]
u_of_s = KS(n0 - s * d0) / KS(s * d1 - n1)
assert u_of_s.numerator().degree() <= 1 and u_of_s.denominator().degree() <= 1


def evaluate_k_at_ks(value):
    value = K(value)

    def evaluate(poly):
        answer = KS.zero()
        for coefficient in reversed(R(poly).list()):
            answer = answer * u_of_s + KS(coefficient)
        return answer

    return evaluate(value.numerator()) / evaluate(value.denominator())


ordinate = evaluate_k_at_ks(ordinate_u)


def evaluate_u_polynomial_at_ks(poly):
    answer = KS.zero()
    for coefficient in reversed(U(poly).list()):
        answer = answer * u_of_s + KS(S(coefficient))
    return answer


assert ordinate ** 2 == evaluate_u_polynomial_at_ks(quartic)

# Standard degree-one pointed-quartic model at (u(s), W(s)).
translated = [
    evaluate_u_polynomial_at_ks(quartic.derivative(order)) / factorial(order)
    for order in range(5)
]
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
A_child = S(rr["child"]["minimal_A_coefficients_low_to_high"])
B_child = S(rr["child"]["minimal_B_coefficients_low_to_high"])
assert 81 * A_pointed == KS(A_child)
assert 729 * B_pointed == KS(B_child)

# The transported effective A1 root is contracted by the pencil to precisely
# the unique I2 support.  Its quartic ordinate identity certifies that this is
# a genuine component of the singular genus-one fibre.
component = source_point(curves["effective_A1_component_on_rootless_source"], polynomial=True)
component_base, component_ordinate = new_base_and_ordinate(component)
assert component_base.denominator().degree() == 0 and component_base.numerator().degree() == 0
component_support = QQ(component_base)
repeated = S(rr["child"]["finite_fibres"][0]["factor_coefficients_low_to_high"])
assert repeated.degree() == 1 and repeated(component_support) == 0

payload = {
    "schema": "elkies-k3.fixed-final-a1-reverse-pointing-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_POINTING",
    "fixed_zero": {
        "new_base_map": rational_record(new_base),
        "new_base_degree": 1,
        "old_base_as_function_of_new_base": rational_record(u_of_s),
        "quartic_ordinate": rational_record(ordinate),
        "exact_quartic_identity": True,
        "degree_one_pointed_quartic_map": True,
        "a_invariants": [rational_record(value) for value in (a1, a2, a3, a4, a6)],
        "short_scaling": "x_short=9*x_pointed, y_short=27*y_pointed",
        "exact_A_identity": "81*A_pointed=A_child",
        "exact_B_identity": "729*B_pointed=B_child",
        "maximum_rational_bits": rational_bits((new_base, u_of_s, ordinate, a1, a2, a3, a4, a6)),
    },
    "effective_A1_component": {
        "new_base_support": str(component_support),
        "new_base_degree": 0,
        "equals_unique_I2_support": True,
        "exact_quartic_identity": True,
        "source_class": curves["effective_A1_component_on_rootless_source"]["class_in_pinned_coordinates"],
        "role": "unique nonidentity A1 component relative to the prescribed pointed zero",
    },
    "marking": {
        "fixed_transition": "A1/MW16 --q6 orbit2247--> rootless/MW17",
        "prescribed_zero_pointed": True,
        "unique_A1_root_bound_to_equation": True,
        "pinned_R17_transport_inherited": True,
    },
    "method": {
        "exact_curve_restriction": True,
        "exact_degree_one_pointing": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SOURCE, CURVES, RR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SOURCE, CURVES, RR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDFINALA1POINT|zero_degree=1|I2_component_degree=0|I2_support={}|"
    "A=1|B=1|bits={}|seconds={:.3f}|status={}|output={}".format(
        component_support, payload["fixed_zero"]["maximum_rational_bits"],
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
