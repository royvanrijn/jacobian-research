#!/usr/bin/env sage
"""Compact the exact fixed-reverse 4A1 model by its four rational I2 fibres.

The current base coordinate has about 20,000-bit I2 supports.  Send three to
0, 1 and infinity.  The fourth becomes 923/3815.  Removing the induced
constant Weierstrass scale turns the million-bit A/B normalization into an
exact QQ model with small primitive coefficients.  This is an explicit base
change and Weierstrass isomorphism, not a twist.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "fixed-reverse-4a1-rr-qq.json"
POINTING = LOCAL / "fixed-reverse-4a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-4a1-compact-crossratio-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return list(map(str, poly.list()))


def rational_bits(values):
    return max(
        max(abs(QQ(value).numerator()).nbits(), QQ(value).denominator().nbits())
        for value in values
    )


def rational_content(poly):
    denominator = ZZ.one()
    for value in poly:
        denominator = denominator.lcm(QQ(value).denominator())
    numerator = ZZ.zero()
    for value in poly:
        numerator = numerator.gcd(ZZ(QQ(value) * denominator))
    return QQ(abs(numerator)) / denominator


started = time.monotonic()
model = read_json(MODEL)
pointing = read_json(POINTING)
assert model["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN"
assert pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING"

Old = PolynomialRing(QQ, "t")
New = PolynomialRing(QQ, "T")
T = New.gen()
A_old = Old(model["child"]["minimal_A_coefficients_low_to_high"])
B_old = Old(model["child"]["minimal_B_coefficients_low_to_high"])

supports = [
    QQ(record["child_I2_support"])
    for record in pointing["effective_horizontal_components"]
]
remaining_factor = Old([
    QQ(value) for value in
    pointing["remaining_vertical_components"][
        "child_I2_support_factor_coefficients_low_to_high"
    ]
])
supports.extend(remaining_factor.roots(QQ, multiplicities=False))
assert len(supports) == 4
r0, r1, rinf, r4 = supports

# T=(t-r0)(r1-rinf)/((t-rinf)(r1-r0)); normalize the inverse denominator
# to have constant coefficient one.
left = r1 - r0
right = r1 - rinf
base_numerator = New([r0, -left * rinf / right])
base_denominator = New([1, -left / right])
assert base_numerator / base_denominator == (
    (left * rinf * T - right * r0) / (left * T - right)
)


def cross_ratio(value):
    if value == rinf:
        return None
    return QQ((value - r0) * (r1 - rinf) / ((value - rinf) * (r1 - r0)))


assert cross_ratio(r0) == 0 and cross_ratio(r1) == 1
fourth_support = cross_ratio(r4)
assert fourth_support == QQ(923) / 3815


def transform_binary(poly, weight):
    return sum((
        QQ(coefficient) * base_numerator ** degree
        * base_denominator ** (weight - degree)
        for degree, coefficient in enumerate(Old(poly).list())
    ), New.zero())


A_raw = transform_binary(A_old, 8)
B_raw = transform_binary(B_old, 12)
content_A = rational_content(A_raw)
content_B = rational_content(B_raw)
assert content_B ** 2 / content_A ** 3 == QQ(4) / 4563

# content_B/content_A = 2*u^2.  This choice yields compact contents
# 4563 and 9126 and, crucially, an isomorphism over QQ rather than a twist.
u_squared = content_B / (2 * content_A)
root_numerator, exact_numerator = u_squared.numerator().nth_root(
    2, truncate_mode=True
)
root_denominator, exact_denominator = u_squared.denominator().nth_root(
    2, truncate_mode=True
)
assert exact_numerator and exact_denominator
u = QQ(root_numerator) / root_denominator
assert u ** 2 == u_squared

A_compact = New(A_raw / u ** 4)
B_compact = New(B_raw / u ** 6)
assert rational_content(A_compact) == 4563
assert rational_content(B_compact) == 9126
assert A_compact.degree() == 8 and B_compact.degree() == 12

Delta = New(-16 * (4 * A_compact ** 3 + 27 * B_compact ** 2))
finite_factors = [
    (factor.monic(), int(exponent)) for factor, exponent in Delta.factor()
]
finite_factor_pattern = sorted(
    (factor.degree(), exponent) for factor, exponent in finite_factors
)
assert finite_factor_pattern == [
    (1, 2), (1, 2), (1, 2), (16, 1),
], finite_factor_pattern
assert Delta.degree() == 22  # the fourth I2 is at infinity
finite_i2 = [factor for factor, exponent in finite_factors if exponent == 2]
assert New.prod(finite_i2).monic() == New(
    T * (T - 1) * (T - fourth_support)
).monic()

payload = {
    "schema": "elkies-k3.fixed-reverse-4a1-compact-crossratio-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_COMPACT_CROSSRATIO",
    "base_change": {
        "meaning": "old_t=base_numerator(T)/base_denominator(T)",
        "base_numerator_coefficients_low_to_high": coefficients(base_numerator),
        "base_denominator_coefficients_low_to_high": coefficients(base_denominator),
        "old_I2_supports_sent_to_0_1_infinity": list(map(str, (r0, r1, rinf))),
        "fourth_compact_I2_support": str(fourth_support),
    },
    "weierstrass_isomorphism": {
        "meaning": (
            "X_raw=base_denominator^4*x_old, Y_raw=base_denominator^6*y_old; "
            "X_raw=u^2*x_compact, Y_raw=u^3*y_compact"
        ),
        "u": str(u),
        "u_squared": str(u_squared),
    },
    "compact_model": {
        "A_coefficients_low_to_high": coefficients(A_compact),
        "B_coefficients_low_to_high": coefficients(B_compact),
        "Delta_coefficients_low_to_high": coefficients(Delta),
        "degrees_A_B_Delta": [8, 12, 22],
        "finite_I2_supports": ["0", "1", str(fourth_support)],
        "infinity": "I2",
        "fibres": "4I2+16I1",
        "ADE": "4A1",
        "maximum_A_B_rational_bits": rational_bits(
            list(A_compact) + list(B_compact)
        ),
    },
    "checks": {
        "exact_base_cross_ratio": True,
        "exact_QQ_weierstrass_isomorphism": True,
        "constant_twist_avoided": True,
        "fibre_euler_number": 24,
        "large_Groebner_required": False,
    },
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, POINTING)
    },
    "runtime_seconds": time.monotonic() - started,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE4A1COMPACT|fourth={}|bits={}|fibres=4I2+16I1|seconds={:.3f}|status={}|output={}".format(
        fourth_support,
        payload["compact_model"]["maximum_A_B_rational_bits"],
        payload["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
