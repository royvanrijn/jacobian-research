#!/usr/bin/env sage -python
"""Enumerate zero-pole sections on the exact q24-derived D12 model mod p.

This is the small-prime seed stage for the corrected orbit42 lift.  It
reconstructs the exact characteristic-zero D12 model directly from the
certified q24 child, reduces it at a requested good prime, and solves the
degree-(4,6) polynomial-section system.  Every returned section is checked on
the curve and its full 13-by-12 coefficient Jacobian rank is recorded.

The output is modular seed data only.  It does not identify MW classes or
construct the orbit42 pencil.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, matrix, sage_eval


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
EXACT_D12 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

p = ZZ(args.prime)
if not p.is_prime() or p in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
F = GF(p)

exact = json.loads(EXACT_D12.read_text())
assert exact["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"

VR = PolynomialRing(QQ, "V")
V = VR.gen()
VF = VR.fraction_field()
RQ = PolynomialRing(QQ, "u")
u = RQ.gen()
KQ = RQ.fraction_field()


def parse_vf(text):
    return VF(sage_eval(str(text), locals={"V": V}))


def monic_power_root(poly, exponent):
    poly = RQ(poly)
    if not poly:
        raise ArithmeticError("zero has no power root")
    lc = QQ(poly.leading_coefficient())
    target = RQ(poly / lc)
    if target.degree() == 0:
        if target != 1:
            raise ArithmeticError("constant monic normalization failed")
        return RQ.one()
    if target.degree() % exponent:
        raise ArithmeticError("power-root degree is not divisible by exponent")
    degree = target.degree() // exponent
    coefficients = [QQ.zero()] * degree + [QQ.one()]
    for j in range(degree - 1, -1, -1):
        current = RQ(coefficients)
        k = (exponent - 1) * degree + j
        coefficients[j] = QQ(target[k] - (current**exponent)[k]) / QQ(exponent)
    root = RQ(coefficients)
    if root**exponent != target:
        raise ArithmeticError("exact monic power-root reconstruction failed")
    return root


raw = exact["jacobian_raw"]
Aorig = parse_vf(raw["A"])
Borig = parse_vf(raw["B"])
Jmap = VF(6912) * Aorig**3 / (VF(4) * Aorig**3 + VF(27) * Borig**2)
N = VR(Jmap.numerator())
Den = VR(Jmap.denominator())
common = N.gcd(Den)
if common.degree() > 0:
    N //= common
    Den //= common
if (N.degree(), Den.degree()) != (18, 18):
    raise ArithmeticError("unexpected exact D12 j-map degrees")

repeated = VR(Den)
for unused in range(7):
    repeated = repeated.gcd(repeated.derivative())
if repeated.degree() != 1:
    raise ArithmeticError("could not isolate the exact I8* factor")
repeated = repeated.monic()
r = -QQ(repeated[0]) / QQ(repeated[1])

RT = PolynomialRing(QQ, "T")
T = RT.gen()


def invpoly(poly):
    poly = VR(poly)
    return RT(sum(
        QQ(poly[i]) * (QQ(r) * T + 1)**i * T**(18 - i)
        for i in range(poly.degree() + 1)
    ))


Pj = invpoly(N)
Qj = invpoly(Den)
common = Pj.gcd(Qj)
if common.degree() > 0:
    Pj //= common
    Qj //= common
if (Pj.degree(), Qj.degree()) != (18, 10):
    raise ArithmeticError("unexpected inverted D12 j-map degrees")
lc = Qj.leading_coefficient()
Pj /= lc
Qj /= lc

center = -QQ(Qj[9]) / QQ(10)
RS = PolynomialRing(QQ, "S")
S = RS.gen()
P1 = RS(Pj(S + center))
Q1 = RS(Qj(S + center))
if not (Q1.is_monic() and Q1[9] == 0 and Q1[8] and Q1[7]):
    raise ArithmeticError("exact centered D12 normalization failed")

base_scale = QQ(Q1[7]) / QQ(Q1[8])
P2 = RQ(P1(base_scale * u))
Q2 = RQ(Q1(base_scale * u))
lc2 = Q2.leading_coefficient()
P2 /= lc2
Q2 /= lc2
if not (
    Q2.is_monic() and Q2.degree() == 10 and Q2[9] == 0
    and Q2[8] == Q2[7]
):
    raise ArithmeticError("exact scaled D12 normalization failed")

a = monic_power_root(P2, 3)
b = monic_power_root(P2 - QQ(1728) * Q2, 2)
if (a.degree(), b.degree()) != (6, 9):
    raise ArithmeticError("unexpected canonical D12 coefficient degrees")

vmap = KQ(r) + KQ.one() / (KQ(base_scale) * KQ(u) + KQ(center))


def eval_v_rational(value, argument):
    value = VF(value)
    return KQ(value.numerator()(argument)) / KQ(value.denominator()(argument))


Aeval = eval_v_rational(Aorig, vmap)
Beval = eval_v_rational(Borig, vmap)
Acan = KQ(-QQ(3) * a)
Bcan = KQ(QQ(2) * b)
cA = Aeval / Acan
cB = Beval / Bcan
wfun = cB / cA
if cA != wfun**2 or cB != wfun**3:
    raise ArithmeticError("exact D12 Weierstrass scaling failed")

wn = RQ(wfun.numerator())
wd = RQ(wfun.denominator())
sn = monic_power_root(RQ(wn / wn.leading_coefficient()), 2)
sd = monic_power_root(RQ(wd / wd.leading_coefficient()), 2)
square_part = KQ(sn) / KQ(sd)
Dfun = wfun / square_part**2
Dnum = RQ(Dfun.numerator())
Dden = RQ(Dfun.denominator())
if Dnum.degree() > 0 or Dden.degree() > 0:
    raise ArithmeticError("exact D12 twist is not constant")
twist = QQ(Dnum[0]) / QQ(Dden[0])
At = RQ(-QQ(3) * twist**2 * a)
Bt = RQ(QQ(2) * twist**3 * b)


def red_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(f"bad reduction denominator at p={p}")
    return F(ZZ(value.numerator())) / F(denominator)


R = PolynomialRing(F, "u")
uf = R.gen()
A = R([red_q(v) for v in At.list()])
B = R([red_q(v) for v in Bt.list()])

print(
    "Q42ZPSMALL|stage=MODEL|"
    f"prime={p}|Adeg={A.degree()}|Bdeg={B.degree()}|"
    f"I8star={int(red_q(r))}|twist={int(red_q(twist))}|status=PASS",
    flush=True,
)


def solve_dx4():
    answer = []
    print(
        "Q42ZPSMALL|stage=DX4|method=FIXED_LEADING_PARAMETER|"
        f"branches={p-1}|status=BEGIN",
        flush=True,
    )
    for s_integer in range(1, int(p)):
        s_value = F(s_integer)
        names = ("x3", "x2", "x1", "x0")
        SR = PolynomialRing(F, names=names, order="degrevlex")
        x3, x2, x1, x0 = SR.gens()
        K = SR.fraction_field()
        U = PolynomialRing(K, "z")
        z = U.gen()
        AA = U([K(v) for v in A.list()])
        BB = U([K(v) for v in B.list()])
        x = (
            K(s_value**2) * z**4 + K(x3) * z**3 + K(x2) * z**2
            + K(x1) * z + K(x0)
        )
        rhs = x**3 + AA * x + BB
        if rhs.degree() > 12:
            raise ArithmeticError("degree-(4,6) RHS exceeds degree 12")
        ys = {6: K(s_value**3)}
        for degree in range(11, 5, -1):
            j = degree - 6
            known = sum(
                ys[i] * ys[degree - i]
                for i in ys
                if (degree - i) in ys and i != 6 and (degree - i) != 6
            )
            ys[j] = (K(rhs[degree]) - known) / (K(2) * ys[6])
        y = sum(ys[i] * z**i for i in range(7))
        residual = y**2 - rhs
        equations = [SR(K(residual[k]).numerator()) for k in range(6)]
        solutions = SR.ideal(equations).variety()
        if solutions:
            print(
                "Q42ZPSMALL_BRANCH|"
                f"s={s_integer}|solutions={len(solutions)}|status=PASS",
                flush=True,
            )
        for solution in solutions:
            values = {g: F(solution[g]) for g in SR.gens()}
            xx = R([
                values[x0], values[x1], values[x2], values[x3], s_value**2
            ])
            rhs_value = xx**3 + A * xx + B
            yy = [F.zero()] * 7
            yy[6] = s_value**3
            for degree in range(11, 5, -1):
                j = degree - 6
                known = sum(
                    yy[i] * yy[degree - i]
                    for i in range(7)
                    if 0 <= degree - i < 7 and i != 6 and (degree - i) != 6
                )
                yy[j] = (rhs_value[degree] - known) / (F(2) * yy[6])
            yy = R(yy)
            if yy**2 != rhs_value:
                raise ArithmeticError("modular dx4 section identity failed")
            answer.append((xx, yy))
    print(
        "Q42ZPSMALL|stage=DX4|method=FIXED_LEADING_PARAMETER|"
        f"raw_solutions={len(answer)}|status=PASS",
        flush=True,
    )
    return answer


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def coefficient_jacobian(x_value, y_value):
    derivative_x = -3 * x_value**2 - A
    derivative_y = 2 * y_value
    columns = []
    for power in range(5):
        columns.append(padded(uf**power * derivative_x, 13))
    for power in range(7):
        columns.append(padded(uf**power * derivative_y, 13))
    return matrix(F, columns).transpose()


sections = solve_dx4()
unique = {}
for x_value, y_value in sections:
    unique[(tuple(x_value.list()), tuple(y_value.list()))] = (x_value, y_value)
sections = list(unique.values())

records = []
for index, (x_value, y_value) in enumerate(sections):
    rank = int(coefficient_jacobian(x_value, y_value).rank())
    records.append({
        "index": index,
        "x_coefficients_low_to_high": [int(v) for v in x_value.list()],
        "y_coefficients_low_to_high": [int(v) for v in y_value.list()],
        "coefficient_jacobian_rank": rank,
    })
    print(
        "Q42ZPSMALL_SECTION|"
        f"index={index}|xdeg={x_value.degree()}|ydeg={y_value.degree()}|"
        f"jacobian_rank={rank}|status=PASS_SECTION",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-zero-pole-smallprime.v1",
    "status": "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS",
    "prime": int(p),
    "input": str(EXACT_D12.relative_to(ROOT)),
    "exact_model": {
        "A_coefficients_low_to_high": [str(v) for v in At.list()],
        "B_coefficients_low_to_high": [str(v) for v in Bt.list()],
        "I8star_root": str(r),
        "center": str(center),
        "base_scale": str(base_scale),
        "twist": str(twist),
    },
    "section_count": len(records),
    "isolated_count": sum(r["coefficient_jacobian_rank"] == 12 for r in records),
    "sections": records,
    "proof_boundary": (
        "Exact reduction of the certified D12 parent and modular polynomial-"
        "section enumeration only. No characteristic-zero section, MW matching, "
        "orbit42 pencil, or A11 child is claimed."
    ),
}

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-orbit42-zero-pole-seeds-mod-{p}.json"
)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42ZPSMALL_RESULT|"
    f"prime={p}|sections={len(records)}|"
    f"isolated={payload['isolated_count']}|"
    "status=PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS",
    flush=True,
)
