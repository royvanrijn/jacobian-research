#!/usr/bin/env sage
"""Verify the two CM Inose equations anchoring the H2 Kumar family."""

from sage.all import Infinity, PolynomialRing, QQ, factor, hilbert_class_polynomial


R = PolynomialRing(QQ, "T")
T = R.gen()

# Calibrate the inverse Kumar map against Elkies's published discriminant-6
# family.  Here z is the Shimura parameter called b in the source, while the
# middle Kumar coefficient is named b_kumar to avoid a collision.
S = PolynomialRing(QQ, "z")
z = S.gen()
a, a_prime = -3 * z, S(1)
b_kumar, b_prime, b_double_prime = -2 * z**2, -(z + 1), -z**3
I2 = -24 * b_prime / a_prime
I4 = -12 * a
I6 = 96 * a * b_prime / a_prime - 36 * b_kumar
I10 = -4 * a_prime * b_double_prime
assert (I2, I4, I6, I10) == (
    24 * z + 24,
    36 * z,
    72 * z * (5 * z + 4),
    4 * z**3,
)
print(
    f"KUMARCMINOSE|calibration=D6|I2={I2}|I4={I4}|I6={I6}|I10={I10}",
    flush=True,
)


def order_at_zero(poly):
    poly = R(poly)
    if not poly:
        return +Infinity
    return min(index for index, coefficient in enumerate(poly.list()) if coefficient)


def inose_anchor(discriminant, A, B, expected_residual):
    hilbert = R(hilbert_class_polynomial(discriminant))
    j_sum = -hilbert[1] if hilbert.degree() == 2 else QQ(0)
    j_product = hilbert[0] if hilbert.degree() == 2 else QQ(0)
    assert A**3 == j_product / 12**6
    assert B**2 == 1 - j_sum / 12**3 + j_product / 12**6

    f = -3 * A * T**4
    g = T**5 * (T**2 - 2 * B * T + 1)
    delta = -16 * (4 * f**3 + 27 * g**2)
    residual = R(delta / (-432 * T**10))
    assert order_at_zero(delta) == 10
    assert 24 - delta.degree() == 10
    assert residual == expected_residual
    print(
        f"KUMARCMINOSE|Delta={discriminant}|H={hilbert}|A={A}|B={B}|"
        f"f={f}|g={g}|residualDelta={factor(residual)}|fibers_at_0_inf=II*,II*",
        flush=True,
    )
    return f, g, delta


# H_-3(X)=X gives j1=j2=0, hence A=0 and B=1.  The residual square
# produces the extra IV fiber and exactly Utsumi No.1.
f3, g3, delta3 = inose_anchor(-3, QQ(0), QQ(1), (T - 1) ** 4)
assert f3 == 0
assert g3 == T**5 * (T - 1) ** 2

# H_-24 has conjugate roots 2417472 +/- 1707264*sqrt(2).  Their symmetric
# functions give A^3=17^3 and B^2=46^2, so this rational Inose model is the
# marked E8+E8/MW diag(4,6) CM anchor, up to quadratic twist.
residual24 = (T**2 - 92 * T + 1) ** 2 - 4 * 17**3 * T**2
f24, g24, delta24 = inose_anchor(-24, QQ(17), QQ(46), residual24)
assert residual24.is_squarefree()

print("KUMARCMINOSE|status=PASS", flush=True)
