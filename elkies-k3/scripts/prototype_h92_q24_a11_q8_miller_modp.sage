#!/usr/bin/env sage -python
"""Prototype the A11 q8 pencil by an exact modular Miller function.

Let ``M`` be the already exact target-coset bridge, let ``R=P12-M`` be the
marked degree-one residual point, and put ``R'=R+O_pinned``.  On every old
A11 Weierstrass fibre the quotient of Miller addition functions

    u = g(M,R) / g(R,O_pinned)

has horizontal divisor ``M+R'-O_pinned-P12``.  Substitution into the
Weierstrass equation has a U-independent base-point factor; after exact
saturation the quotient is quadratic in x.  Its discriminant is therefore
the degree-two cover of the old base that should give the q8 neighbour.

This script performs that construction over the pinned good-reduction field,
extracts the squarefree quartic, and classifies its short Jacobian.  It is a
construction regression, not a characteristic-zero H0 or equation proof.
No Groebner basis is used.
"""

import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
BRIDGE = LOCAL / "q24-a11-bridge-m-section-marked-qq.json"
ZERO = LOCAL / "q24-a11-pinned-zero-section-qq.json"
MARKED = LOCAL / "q24-a11-q8-horizontal-points-mod100003.json"
OUTPUT = LOCAL / "q24-a11-q8-miller-pencil-mod100003.json"

for path in (MODEL, BRIDGE, ZERO, MARKED):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
bridge = json.loads(BRIDGE.read_text())
zero = json.loads(ZERO.read_text())
marked = json.loads(MARKED.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert bridge["status"] == "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ"
assert zero["status"] == "PASS_EXACT_MARKED_A11_PINNED_ZERO_SECTION_QQ"
assert marked["status"] == "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP"

p = ZZ(marked["prime"])
assert p == 100003
F = GF(p)
FT = PolynomialRing(F, "T")
T0 = FT.gen()
K = FT.fraction_field()
T = K(T0)


def red_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError("bad denominator at pinned prime")
    return F(value.numerator()) / F(value.denominator())


def exact_poly(values):
    return FT([red_q(value) for value in values])


def modular_poly(values):
    return FT([F(value) for value in values])


def load_point(section, exact):
    loader = exact_poly if exact else modular_poly
    X = loader(section["X_coefficients_low_to_high"])
    Y = loader(section["Y_coefficients_low_to_high"])
    Z = loader(section["Z_coefficients_low_to_high"])
    return K(X) / K(Z) ** 2, K(Y) / K(Z) ** 3


A = K(exact_poly(model["child"]["minimal_A_coefficients_low_to_high"]))
B = K(exact_poly(model["child"]["minimal_B_coefficients_low_to_high"]))
xM, yM = load_point(bridge["section"], exact=True)
xO, yO = load_point(zero["section"], exact=True)
residual_record = marked["residual_P12_minus_M"]["section"]
xR, yR = load_point(residual_record, exact=False)
xP, yP = load_point(marked["q8_target"]["section"], exact=False)

for label, xx, yy in (("M", xM, yM), ("O_pinned", xO, yO), ("R", xR, yR), ("P12", xP, yP)):
    if yy**2 != xx**3 + A * xx + B:
        raise ArithmeticError(f"{label} misses the reduced A11 equation")

# Check the marked group-law identity before using its Miller function.
slope = (yR - yM) / (xR - xM)
x_sum = slope**2 - xM - xR
y_sum = slope * (xM - x_sum) - yM
if (x_sum, y_sum) != (xP, yP):
    raise ArithmeticError("the marked identity M+R=P12 failed")

# The desired horizontal poles are O_pinned and P12, not the Weierstrass
# origin.  Put R'=R+O_pinned.  The quotient of Miller addition functions
#
#     g(M,R) / g(R,O_pinned)
#
# has divisor M+R'-O_pinned-P12.  In line/vertical factors it is
#
#     line(M,R)*(x-x(R')) / ((x-x(P12))*line(R,O_pinned)).
def group_add(x1, y1, x2, y2):
    addition_slope = (y2 - y1) / (x2 - x1)
    x3 = addition_slope**2 - x1 - x2
    y3 = addition_slope * (x1 - x3) - y1
    return x3, y3


xRprime, yRprime = group_add(xR, yR, xO, yO)
x_total_left, y_total_left = group_add(xM, yM, xRprime, yRprime)
x_total_right, y_total_right = group_add(xO, yO, xP, yP)
if (x_total_left, y_total_left) != (x_total_right, y_total_right):
    raise ArithmeticError("four-point Miller divisor does not have group sum zero")

RU = PolynomialRing(K, "U")
Uring = RU.gen()
RX = PolynomialRing(RU, "x")
x = RX.gen()
cubic = x**3 + RX(A) * x + RX(B)


def line_constant(x1, y1, x2, y2):
    line_slope = (y2 - y1) / (x2 - x1)
    # line(A,B) = y + returned polynomial.
    return RX(-y1) - RX(line_slope) * (x - RX(x1))


b_mr = line_constant(xM, yM, xR, yR)
b_ro = line_constant(xR, yR, xO, yO)
a_numerator = x - RX(xRprime)
b_numerator = b_mr * (x - RX(xRprime))
a_denominator = x - RX(xP)
b_denominator = b_ro * (x - RX(xP))

# U*denominator-numerator = a*y+b.  Eliminating y gives b^2-a^2*cubic.
a_linear = RX(Uring) * a_denominator - a_numerator
b_linear = RX(Uring) * b_denominator - b_numerator
raw_cover = b_linear**2 - a_linear**2 * cubic

# Saturate the U-independent base points coefficientwise.  The remaining
# polynomial must be quadratic in x because the q8 fibre has old-base degree 2.
KX = PolynomialRing(K, "xbase")
xbase = KX.gen()
by_u_degree = {}
for x_degree, coefficient_u in enumerate(raw_cover.list()):
    for u_degree, coefficient in enumerate(coefficient_u.list()):
        by_u_degree[u_degree] = by_u_degree.get(u_degree, KX.zero()) + K(coefficient) * xbase**x_degree
base_factor = None
for coefficient_x in by_u_degree.values():
    if coefficient_x:
        base_factor = coefficient_x if base_factor is None else base_factor.gcd(coefficient_x)
if base_factor is None:
    raise ArithmeticError("four-point Miller elimination vanished")
base_factor = base_factor.monic()
quadratic, remainder = raw_cover.quo_rem(RX(base_factor))
if remainder or quadratic.degree() != 2:
    raise ArithmeticError(
        f"four-point Miller saturation left x-degree {quadratic.degree()} "
        f"after base-factor degree {base_factor.degree()}"
    )
qa, qb, qc = quadratic[2], quadratic[1], quadratic[0]
cover_discriminant = qb**2 - 4 * qa * qc

# Clear only the common old-base denominator.  Multiplying numerator and
# denominator gives the exact square class in F[T,U].
common_denominator = FT.one()
for coefficient in cover_discriminant.list():
    common_denominator = common_denominator.lcm(coefficient.denominator())
BTU = PolynomialRing(F, names=("Tb", "U"))
Tb, U = BTU.gens()
cover_numerator = BTU.zero()
for u_degree, coefficient in enumerate(cover_discriminant.list()):
    cleared = coefficient.numerator() * (common_denominator // coefficient.denominator())
    cover_numerator += BTU(cleared(Tb)) * U**u_degree
square_class_polynomial = cover_numerator * BTU(common_denominator(Tb))
factorization = square_class_polynomial.factor()
odd_part = BTU(factorization.unit())
factor_profile = []
for factor, exponent in factorization:
    exponent = int(exponent)
    factor_profile.append((str(factor), int(factor.degree(Tb)), int(factor.degree(U)), exponent))
    if exponent % 2:
        odd_part *= factor

# Regard the remaining square class as a polynomial in the old base Tb over F(U).
FU = PolynomialRing(F, "U")
U0 = FU.gen()
KU = FU.fraction_field()
RT = PolynomialRing(KU, "T")
Tv = RT.gen()
quartic = RT.zero()
for (degree_t, degree_u), coefficient in odd_part.dict().items():
    quartic += KU(FU(coefficient) * U0**degree_u) * Tv**degree_t
if quartic.degree() != 4:
    raise ArithmeticError(
        f"Miller square class has old-base degree {quartic.degree()}, not 4; "
        f"factor profile={factor_profile}"
    )

coefficients = list(quartic.list()) + [KU.zero()] * 5
e, d, c, b, a = coefficients[:5]
Iinv = 12 * a * e - 3 * b * d + c**2
Jinv = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
jacA = -27 * Iinv
jacB = -27 * Jinv
jacDelta = -16 * (4 * jacA**3 + 27 * jacB**2)


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return n - 1, n, n, f"I{n}"
    if ord_delta == 2:
        return 0, 1, 2, "II"
    if ord_delta == 3:
        return 1, 2, 3, "III"
    if ord_delta == 4:
        return 2, 3, 4, "IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, 4, 6, "I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta - 6)
        return n + 4, 4, n + 6, f"I{n}*"
    if ord_delta == 8:
        return 6, 3, 8, "IV*"
    if ord_delta == 9:
        return 7, 2, 9, "III*"
    if ord_delta == 10:
        return 8, 1, 10, "II*"
    raise ArithmeticError((ord_a, ord_b, ord_delta))


def valuation(value, factor):
    return int(value.numerator().valuation(factor) - value.denominator().valuation(factor))


factors = set()
for value in (jacA, jacB, jacDelta):
    for polynomial in (value.numerator(), value.denominator()):
        factors.update(factor for factor, unused in polynomial.factor())

finite_fibres = []
root_rank = 0
root_determinant = 1
euler = 0
for factor in sorted(factors, key=str):
    raw = (valuation(jacA, factor), valuation(jacB, factor), valuation(jacDelta, factor))
    scale = min(raw[0] // 4, raw[1] // 6, raw[2] // 12)
    orders = (raw[0] - 4 * scale, raw[1] - 6 * scale, raw[2] - 12 * scale)
    if orders[2] <= 0:
        continue
    rank, determinant, local_euler, kind = kodaira_data(*orders)
    degree = int(factor.degree())
    root_rank += degree * rank
    root_determinant *= determinant**degree
    euler += degree * local_euler
    finite_fibres.append(
        {"factor": str(factor), "degree": degree, "orders_A_B_Delta": list(orders), "kind": kind}
    )

raw_infinity = tuple(
    int(value.denominator().degree() - value.numerator().degree())
    for value in (jacA, jacB, jacDelta)
)
scale_infinity = min(raw_infinity[0] // 4, raw_infinity[1] // 6, raw_infinity[2] // 12)
infinity_orders = tuple(
    raw_infinity[index] - (4, 6, 12)[index] * scale_infinity for index in range(3)
)
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    rank, determinant, local_euler, infinity_kind = kodaira_data(*infinity_orders)
    root_rank += rank
    root_determinant *= determinant
    euler += local_euler

status = (
    "PASS_MODP_A11_Q8_MILLER_2A5_REGRESSION"
    if (root_rank, root_determinant, euler) == (10, 36, 24)
    else "FAIL_MODP_A11_Q8_MILLER_CHILD_SIGNATURE"
)
payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-miller-pencil-modp.v1",
    "status": status,
    "prime": int(p),
    "inputs": {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (MODEL, BRIDGE, ZERO, MARKED)
    },
    "marked_group_law": "M+(P12-M)=P12",
    "miller_function": "u=g(M,R)/g(R,O_pinned)",
    "base_point_factor_degree_in_x": int(base_factor.degree()),
    "quadratic_after_base_point_saturation": True,
    "cover_square_class_factor_profile": factor_profile,
    "quartic": {
        "old_base_degree": int(quartic.degree()),
        "coefficients_low_to_high": [str(value) for value in quartic.list()],
    },
    "jacobian": {
        "A": str(jacA),
        "B": str(jacB),
        "finite_fibres": finite_fibres,
        "infinity_orders_A_B_Delta": list(infinity_orders),
        "infinity_kind": infinity_kind,
        "root_rank": int(root_rank),
        "root_determinant": int(root_determinant),
        "euler_number": int(euler),
        "MW_rank_if_rho19": int(17 - root_rank),
    },
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact arithmetic over the pinned good-reduction field verifies the Miller divisor identity, "
        "base-point saturation, binary quartic, and modular Jacobian fibre signature. Characteristic-zero "
        "resolved H0 normalization, exact quartic/Jacobian, and marked transport remain open."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"A11Q8MILLER|prime={p}|quartic_degree={quartic.degree()}|root_rank={root_rank}|"
    f"root_det={root_determinant}|euler={euler}|MW={17-root_rank}|status={status}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
