#!/usr/bin/env sage
"""
Exact characteristic-zero q6_7774 experiment.

Stage 1:
  * load the certified normalized orbit1222 parent + exact P3;
  * identify exact I7/I2 factors lifting V=6,5 mod 73;
  * derive the selected resolved chord residues:
      - I7: branch tangent at the Weierstrass node (P3 specializes to the node);
      - I2: node-to-P3 chord value (P3 is smooth there);
  * reconstruct the linear correction a(V);
  * verify a(V) mod 73 == 4V+64;
  * form T=(m+a(V))/((V-alpha7)(V-alpha2));
  * obtain the genus-one quartic and Jacobian by binary-quartic invariants;
  * classify finite fibres and check the infinity valuation;
  * audit reduction against the pinned GF(73) q6_7774 model up to constant
    Weierstrass scaling.

This deliberately does NOT certify h0(D)=2 yet.  The next gate is the full
resolved I7 quotient/jet matrix; this script first pins the exact parameter and
binary-quartic child independently.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, GF, PolynomialRing, QuadraticField, prod
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
OUT = ROOT / "artifacts" / "generated-results"
OUT.mkdir(parents=True, exist_ok=True)

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(DATA / "q80-orbit1222-char0" / "q80_char0_orbit1222_P1_P3_normalized.sage"))

# The normalized exact parent uses K=QQ(j), j^2=-3 and R=K[V].
assert j**2 == -3

# Exact discriminant support from the independently certified orbit1222 model.
P7 = V^2 + (-QQ(65223)/11438*j - QQ(97929)/11438)*V \
     - QQ(397791)/45752*j - QQ(930933)/45752
P2 = V^3 + (-QQ(2354225337)/232790428*j + QQ(6190374699)/232790428)*V^2 \
     + (-QQ(23116822551)/465580856*j - QQ(87776231337)/465580856)*V \
     + QQ(20457191439)/465580856*j + QQ(120333555774)/58197607

# The two exact roots selected by the pinned j -> 17 (mod 73) place.
alpha7 = (-QQ(135) + QQ(9)*j) / 76
alpha2 = (QQ(18765) - QQ(15471)*j) / 5668
assert P7(alpha7) == 0
assert P2(alpha2) == 0

Fp = GF(73)
JMOD = Fp(17)
assert JMOD**2 == Fp(-3)

def red_q(q):
    q = QQ(q)
    return Fp(q.numerator()) / Fp(q.denominator())

def red_k(z):
    z = K(z)
    cc = list(z)
    cc += [QQ(0)] * (2-len(cc))
    return red_q(cc[0]) + red_q(cc[1])*JMOD

assert red_k(alpha7) == 6
assert red_k(alpha2) == 5

def exact_node_x(r):
    ar = K(A(r))
    br = K(B(r))
    if ar == 0:
        raise ArithmeticError("unexpected c4=0 at multiplicative fibre")
    x0 = -3*br/(2*ar)
    assert x0^3 + ar*x0 + br == 0
    assert 3*x0^2 + ar == 0
    return K(x0)

x07 = exact_node_x(alpha7)
x02 = exact_node_x(alpha2)

# IMPORTANT: the modular c7=58 is NOT the nodal tangent slope.
# It is the free coordinate of the connected-A6 quotient line
#     (1, alpha7, alpha7^2, c7).
# The old modular scan selected it globally from the square-factor/quartic
# condition. First identify the constant Weierstrass gauge of this normalized
# exact parent at the pinned p=73 place, then recover c7 exactly in the same way.

FpV0 = PolynomialRing(Fp, "V0")
v0 = FpV0.gen()

def red_poly0(f):
    f = R(f)
    return FpV0([red_k(c) for c in f.list()])

A_pin = FpV0([23,17,63,2,58,33,47,16,6])
B_pin = FpV0([43,47,57,50,8,54,20,14,45,61,64,0,33])
PX_pin = FpV0([47,69,7,4,51])
PY_pin = FpV0([22,2,48,23,8,47,35])

Ared, Bred = red_poly0(A), red_poly0(B)
PXred, PYred = red_poly0(P3x), red_poly0(P3y)

gauge_hits = []
for u in Fp:
    if u == 0:
        continue
    if Ared != u^4*A_pin or Bred != u^6*B_pin or PXred != u^2*PX_pin:
        continue
    for eps in (1, -1):
        if PYred == Fp(eps)*u^3*PY_pin:
            gauge_hits.append((u, eps))

chord_scale_values = sorted(set(int(Fp(eps)*u) for u, eps in gauge_hits))
if len(chord_scale_values) != 1:
    raise ArithmeticError("parent gauges do not give a unique chord scaling: {}".format(gauge_hits))

chord_scale73 = Fp(chord_scale_values[0])
u73, section_sign = next(
    (u, eps) for u, eps in gauge_hits
    if Fp(eps)*u == chord_scale73
)
scale_int = int(chord_scale73)
if scale_int > 36:
    scale_int -= 73
gauge = K(scale_int)

print(
    "Q80Q67774GAUGE|u73={}|section_sign={}|chord_scale73={}|gauge={}|"
    "status=PASS_PARENT_GAUGE".format(
        int(u73), section_sign, int(chord_scale73), gauge
    )
)

# The selected modular branch is the chord through -P3:
#     m=(y+P3y)/(x-P3x).
# At the selected I2, P3 is smooth, so this has an intrinsic nodal value.
x3_2 = K(P3x(alpha2))
y3_2 = K(P3y(alpha2))
if x3_2 == x02 and y3_2 == 0:
    raise ArithmeticError("P3 unexpectedly specializes to the selected I2 node")
c2 = K(y3_2/(x02-x3_2))
target_c2_73 = chord_scale73*Fp(62)
assert red_k(c2) == target_c2_73

print(
    "Q80Q67774I2ROW|c2={}|mod73={}|normalized_mod73=62|"
    "status=PASS_EXACT_I2_ROW".format(c2, int(red_k(c2)))
)

# Cheap exact c7 candidate from the resolved nodal branch.
# The normalized raw quotient coordinate must reduce to 38*58 = 14 mod 73.
target_c7_73 = chord_scale73*Fp(58)
sq7 = K(3*x07)
if not sq7.is_square():
    raise ArithmeticError("selected I7 nodal tangent is not defined over K")
r7 = K(sq7.sqrt())
c7_candidates = (r7,-r7)
c7 = next((c for c in c7_candidates if red_k(c) == target_c7_73), None)
if c7 is None:
    raise ArithmeticError("no exact I7 branch reduces to raw quotient residue {}: {}".format(int(target_c7_73),[(str(c),int(red_k(c))) for c in c7_candidates]))
print("Q80Q67774C7CANDIDATE|c7={}|mod73={}|normalized_mod73=58|status=PASS_EXACT_C7_CANDIDATE".format(c7,int(red_k(c7))))

# m = T*d - a, so m(alpha_i) = -a(alpha_i) = c_i.
# Interpolate the unique exact linear a(V).
slope = K(((-c7)-(-c2))/(alpha7-alpha2))
intercept = K(-c7 - slope*alpha7)
a_corr = R(slope*V + intercept)
d_old = R((V-alpha7)*(V-alpha2))
assert a_corr(alpha7) == -c7
assert a_corr(alpha2) == -c2

FpV = PolynomialRing(Fp, "V")
v73 = FpV.gen()
def red_poly_v(f):
    f = R(f)
    return FpV([red_k(c) for c in f.list()])

a73 = red_poly_v(a_corr)
d73 = red_poly_v(d_old)
assert a73 == chord_scale73*(4*v73 + 64)
assert d73 == (v73-6)*(v73-5)

print("Q80Q67774ROOTS|alpha7={}|alpha2={}|status=PASS_EXACT_FACTORS".format(alpha7, alpha2))
print("Q80Q67774RESIDUES|c7={}|c2={}|c7mod73={}|c2mod73={}|status=PASS_CHILD_SELECTED_RESIDUES".format(
    c7, c2, int(red_k(c7)), int(red_k(c2))
))
print("Q80Q67774PARAMETER|d={}|a={}|normalized_a={}|mod73_a={}|status=PASS_EXACT_PARAMETER".format(
    gauge*d_old, a_corr, a_corr/gauge, a73
))

# Binary-quartic compiler.  New T is a coefficient-field variable and V
# remains the quartic coordinate.
KTpoly = PolynomialRing(K, "T")
T0 = KTpoly.gen()
KT = KTpoly.fraction_field()
RV = PolynomialRing(KT, "V")
vv = RV.gen()
T = KT(T0)

def lift_v(f):
    f = R(f)
    return RV([KT(c) for c in f.list()])

d_l = lift_v(gauge*d_old)
a_l = lift_v(a_corr)
A_l = lift_v(A)
B_l = lift_v(B)
P3x_l = lift_v(P3x)
P3y_l = -lift_v(P3y)

hop = compile_degree_two_chord_hop(
    RV,
    T,
    d_l, 0,          # s0=d
    a_l, 1,          # s1=a+m
    P3x_l,
    P3y_l,
    A_l,
    B_l,
)

# This is exactly m=T*d-a.
assert hop["chord"] == T*d_l-a_l
quartic = hop["binary_quartic"]
assert quartic.degree() == 4

Aj = KT(hop["jacobian_a"])
Bj = KT(hop["jacobian_b"])
Dj = KT(-16*(4*Aj^3+27*Bj^2))

print("Q80Q67774QUARTIC|degree={}|quartic={}|status=PASS_BINARY_QUARTIC".format(
    quartic.degree(), quartic
))
print("Q80Q67774JACOBIAN_RAW|A={}|B={}|status=PASS_BINARY_QUARTIC_INVARIANTS".format(Aj, Bj))

# Finite minimization / Kodaira classification.  This may be the expensive
# part, but it is far cheaper than the old characteristic-zero normalization.
classification = classify_finite_short_weierstrass_fibres(KTpoly, Aj, Bj)
fibres = [
    {
        "factor": str(item["factor"]),
        "degree": item["degree"],
        "kodaira": item["kodaira"],
        "orders": list(item["minimal_orders"]),
    }
    for item in classification["finite_fibres"]
]
infty = classification["infinity_boundary"]

print("Q80Q67774FIBRES|finite={}|infinity={}".format(fibres, infty))

# Expected global pattern: finite 2 I5 + 2 I2 + 4 I1, infinity I6.
finite_degree_totals = {}
for item in classification["finite_fibres"]:
    symbol = item["kodaira"]
    finite_degree_totals[symbol] = finite_degree_totals.get(symbol, 0) + int(item["degree"])
finite_ok = finite_degree_totals == {"I5": 2, "I2": 2, "I1": 4}
infinity_ok = tuple(infty["normalized_orders"]) == (0, 0, 6)

# Reduction K(T) -> F73(T).
FpTpoly = PolynomialRing(Fp, "T")
t73 = FpTpoly.gen()
FpT = FpTpoly.fraction_field()

def red_kt_poly(poly):
    poly = KTpoly(poly)
    return FpTpoly([red_k(c) for c in poly.list()])

def red_kt(frac):
    frac = KT(frac)
    num = red_kt_poly(frac.numerator())
    den = red_kt_poly(frac.denominator())
    if den == 0:
        raise ZeroDivisionError("denominator vanished modulo 73")
    return FpT(num)/FpT(den)

Aj73 = red_kt(Aj)
Bj73 = red_kt(Bj)

Atarget = FpTpoly(
    46*t73^8 + 5*t73^7 + 16*t73^6 + 44*t73^5 + 6*t73^4
    + 13*t73^3 + t73^2 + t73
)
Btarget = FpTpoly(
    54*t73^12 + 58*t73^11 + 48*t73^10 + 16*t73^9 + 42*t73^8
    + 67*t73^7 + 25*t73^6 + 19*t73^5 + 27*t73^4 + 45*t73^3
    + 61*t73^2 + 44*t73 + 49
)

# Binary-quartic invariants may differ from the pinned equation by a constant
# x=u^2 X, y=u^3 Y scaling.  Search the 72 possible constants exactly.
scale = None
for u in Fp:
    if u == 0:
        continue
    if Aj73 == FpT(u^4*Atarget) and Bj73 == FpT(u^6*Btarget):
        scale = int(u)
        break
    if FpT(u^4)*Aj73 == FpT(Atarget) and FpT(u^6)*Bj73 == FpT(Btarget):
        scale = -int(u)  # negative marker means inverse orientation
        break

reduction_ok = scale is not None
print("Q80Q67774REDUCTION|scale_marker={}|status={}".format(
    scale, "PASS_GF73_7774" if reduction_ok else "FAIL_GF73_7774"
))

status = (
    "PASS_EXACT_Q6_7774_STAGE1"
    if finite_ok and infinity_ok and reduction_ok
    else "PARTIAL_Q6_7774_STAGE1"
)

result = {
    "status": status,
    "field": "QQ(sqrt(-3))",
    "p73_embedding": "j->17",
    "alpha7": str(alpha7),
    "alpha2": str(alpha2),
    "alpha7_mod73": int(red_k(alpha7)),
    "alpha2_mod73": int(red_k(alpha2)),
    "c7": str(c7),
    "c2": str(c2),
    "c7_mod73": int(red_k(c7)),
    "c2_mod73": int(red_k(c2)),
    "denominator_d": str(gauge*d_old),
    "correction_a": str(a_corr),
    "correction_a_mod73": str(a73),
    "quartic": str(quartic),
    "jacobian_a_raw": str(Aj),
    "jacobian_b_raw": str(Bj),
    "finite_fibres": fibres,
    "infinity": {
        "raw_orders": [int(x) for x in infty["raw_orders"]],
        "scaling": int(infty["scaling"]),
        "normalized_orders": [int(x) for x in infty["normalized_orders"]],
    },
    "finite_pattern_ok": finite_ok,
    "infinity_I6_ok": infinity_ok,
    "gf73_reduction_ok": reduction_ok,
    "gf73_scale_marker": scale,
    "next_gate": "resolved I7 connected quotient/jet RR matrix and h0(D)=2",
}

out = OUT / "q80-q6-7774-char0-stage1.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print("Q80Q67774|artifact={}|status={}".format(out, status))

if status != "PASS_EXACT_Q6_7774_STAGE1":
    raise SystemExit(2)
