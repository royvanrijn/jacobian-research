#!/usr/bin/env sage
"""
Exact characteristic-zero probe for the final Q80 q6 RR pencil.

Inputs:
  * selected exact q4-candidate1 parent;
  * the pinned GF(73) final-q6 certificate.

Strategy:
  1. identify the exact parent's modular base coordinate with the historical
     pinned coordinate by matching the five uniquely typed fibres
     I2,I3,I4,I5,I6 (a unique PGL2 map);
  2. transport the pinned H=P2-P3 section through that base map and
     Weierstrass gauge to obtain the correct modular seed in our exact gauge;
  3. Hensel-lift the section subject to its three node hits I3,I4,I6,
     recognize only the two x-parameters in QQ(sqrt(-3)), and recover y
     afterwards as an exact polynomial square;
  4. derive the whole-A4 row exactly at I5;
  5. resolve the I6 torically and select the connected-A5 quotient line
     reducing to the transported pinned residue -4;
  6. certify the resulting 4 -> 2 exact RR kernel.

This intentionally stops before child compilation.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, GF, Qp, PolynomialRing, QuadraticField, FunctionField,
    LaurentSeriesRing, matrix, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ART = ROOT / "artifacts" / "generated-results" / "q80-q4-candidate1-char0-family.json"
OUT = ROOT / "artifacts" / "generated-results"

if not ART.exists():
    raise SystemExit(f"missing candidate1 family artifact: {ART}")

data = json.loads(ART.read_text())
assert data["status"] == "PASS_EXACT_Q4_CANDIDATE1_FAMILY"
selected = data["selected"]
assert selected["target_candidate1"]

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()
R = PolynomialRing(K, "W")
W = R.gen()
KF = R.fraction_field()

A = R(sage_eval(selected["minimal_A"], locals={"j":j, "W":W}))
B = R(sage_eval(selected["minimal_B"], locals={"j":j, "W":W}))
assert A.degree() == 8 and B.degree() == 12

Fp = GF(73)
JMOD = Fp(17)
assert JMOD^2 == Fp(-3)
Rp = PolynomialRing(Fp, "W")
wp = Rp.gen()
FU = Rp.fraction_field()

def red_q(c):
    c = QQ(c)
    return Fp(c.numerator())/Fp(c.denominator())

def red_k(c):
    c = K(c)
    cc = list(c)+[QQ(0),QQ(0)]
    return red_q(cc[0])+JMOD*red_q(cc[1])

def red_poly(f):
    f = R(f)
    return Rp([red_k(c) for c in f.list()])

def eval_poly(poly, value):
    result = value.parent()(0)
    for coefficient in reversed(list(poly)):
        result = result*value + value.parent()(coefficient)
    return result

# ---------------------------------------------------------------------------
# 1. Exact finite fibres and the unique modular PGL2 identification.
# ---------------------------------------------------------------------------
fibres = {}
for item in selected["finite_fibres"]:
    symbol = str(item["kodaira"])
    if symbol not in ("I2","I3","I4","I5","I6") or int(item["degree"]) != 1:
        continue
    factor = R(sage_eval(item["factor"], locals={"j":j, "W":W}))
    assert factor.degree() == 1
    root = K(-factor[0]/factor[1])
    if symbol in fibres:
        raise ArithmeticError(f"duplicate {symbol} fibre in selected candidate1 parent")
    fibres[symbol] = (root,factor)

assert set(fibres) == {"I2","I3","I4","I5","I6"}

exact_mod_root = {symbol:Fp(red_k(root)) for symbol,(root,unused) in fibres.items()}
pinned_root = {
    "I2":Fp(60),
    "I3":Fp(23),
    "I4":Fp(24),
    "I5":Fp(25),
    "I6":Fp(47),
}

# Solve phi(w)=(a*w+b)/(c*w+d) from three typed fibres.
rows = []
for symbol in ("I2","I3","I4"):
    x = exact_mod_root[symbol]
    y = pinned_root[symbol]
    rows.append([x,1,-y*x,-y])
Mmob = matrix(Fp,rows)
ker = Mmob.right_kernel()
assert ker.dimension() == 1
a,b,c,d = ker.basis()[0]
assert a*d-b*c

def phi_value(x):
    den = c*x+d
    if den == 0:
        return None
    return (a*x+b)/den

for symbol in pinned_root:
    got = phi_value(exact_mod_root[symbol])
    assert got == pinned_root[symbol], (symbol,got,pinned_root[symbol])

den = Rp(c*wp+d)
phi = FU(a*wp+b)/FU(den)

print(
    "Q80FINALBASE|mobius=({}*W+{})/({}*W+{})|roots={}|"
    "status=PASS_TYPED_FIBRE_PGL2".format(
        int(a),int(b),int(c),int(d),
        {k:(int(exact_mod_root[k]),int(pinned_root[k])) for k in sorted(pinned_root)}
    ),
    flush=True,
)

# Historical pinned candidate1 parent.
A_pin = Rp([29,63,55,53,9,49,66,72,32])
B_pin = Rp([35,6,56,66,55,35,37,34,22,62,50,65,18])

A_base = FU(den^8) * eval_poly(A_pin,phi)
B_base = FU(den^12) * eval_poly(B_pin,phi)
A_red = FU(red_poly(A))
B_red = FU(red_poly(B))

rA = A_red/A_base
rB = B_red/B_base
u2 = rB/rA
assert u2^2 == rA
assert u2^3 == rB
assert rB.is_square()
u3root = rB.sqrt()

# Pinned final horizontal H=P2-P3.
Hx_pin = Rp([15,29,9,59,4])
Hy_pin = Rp([7,61,7,3,14,57,8])

X_base = FU(den^4) * eval_poly(Hx_pin,phi)
Y_base = FU(den^6) * eval_poly(Hy_pin,phi)
X_seed_f = FU(u2*X_base)

def as_poly(value):
    value = FU(value)
    q = value.denominator()
    if q.degree() != 0:
        return None
    return Rp(value.numerator()/q[0])

X_seed = as_poly(X_seed_f)
if X_seed is None or X_seed.degree() > 4:
    raise ArithmeticError(
        "transported pinned final Hx is not polynomial of degree <=4 in exact gauge"
    )

seed_options = []
for u3 in (u3root,-u3root):
    Y_seed = as_poly(FU(u3*Y_base))
    if Y_seed is not None and Y_seed.degree() <= 6:
        seed_options.append((u3,Y_seed))
if not seed_options:
    raise ArithmeticError(
        "neither transported pinned Hy sign is polynomial of degree <=6 in exact gauge"
    )

# Use the first sign as the pinned continuation; the other is its -H mate.
u3,Y_seed = seed_options[0]

print(
    "Q80FINALMODSEED|X={}|Y={}|sign_options={}|"
    "status=PASS_TRANSPORTED_PINNED_HORIZONTAL".format(
        X_seed,Y_seed,len(seed_options)
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Hensel-lift the exact H from the three node hits I3,I4,I6.
# ---------------------------------------------------------------------------
def node_x(symbol):
    root = fibres[symbol][0]
    ar,br = K(A(root)),K(B(root))
    assert ar != 0
    x0 = K(-3*br/(2*ar))
    assert x0^3+ar*x0+br == 0
    assert 3*x0^2+ar == 0
    return x0

hit_symbols = ("I3","I4","I6")
points = [(fibres[s][0],node_x(s)) for s in hit_symbols]
N0 = R.lagrange_polynomial(points)
assert N0.degree() <= 2

Fhit = R.one()
for s in hit_symbols:
    Fhit *= W-fibres[s][0]
assert Fhit.degree() == 3

N0p = red_poly(N0)
Fhitp = red_poly(Fhit)

qseed,rem = (X_seed-N0p).quo_rem(Fhitp)
assert rem == 0 and qseed.degree() <= 1

yseed_q,rem = Y_seed.quo_rem(Fhitp)
assert rem == 0 and yseed_q.degree() <= 3

qseed = Rp(qseed)
yseed_q = Rp(yseed_q)
seed = [
    qseed[0] if qseed.degree() >= 0 else Fp(0),
    qseed[1] if qseed.degree() >= 1 else Fp(0),
] + [
    yseed_q[i] if i <= yseed_q.degree() else Fp(0)
    for i in range(4)
]


# ---------------------------------------------------------------------------
# Recover the final horizontal indirectly from easier exact MW sections.
# ---------------------------------------------------------------------------

from itertools import combinations
from sage.all import EllipticCurve

def poly_square_roots(poly):
    poly = R(poly)
    if poly == 0:
        return (R.zero(),)
    fac = poly.factor()
    unit = K(fac.unit())
    if not unit.is_square():
        return ()
    root = R(unit.sqrt())
    for factor, exponent in fac:
        exponent = int(exponent)
        if exponent % 2:
            return ()
        root *= R(factor) ** (exponent // 2)
    assert root^2 == poly
    return (root, -root)

def exact_hitset(Hx,Hy):
    out = []
    for symbol in fibres:
        root = fibres[symbol][0]
        if Hx(root) == node_x(symbol) and Hy(root) == 0:
            out.append(symbol)
    return tuple(sorted(out))

def candidate_from_fixed_x(Hx, Fhit, label, source_parameter=None):
    Hx = R(Hx)
    Fhit = R(Fhit)
    rhs = R(Hx^3 + A*Hx + B)
    quotient, remainder = rhs.quo_rem(Fhit^2)
    if remainder:
        return []
    roots = poly_square_roots(quotient)
    answer = []
    for G in roots:
        Hy = R(Fhit*G)
        if Hx.degree() > 4 or Hy.degree() > 6:
            continue
        if Hy^2 != rhs:
            continue
        answer.append({
            "Hx":Hx,
            "Hy":Hy,
            "hitset":exact_hitset(Hx,Hy),
            "source_subset":tuple(label),
            "parameter":source_parameter,
        })
    return answer

symbols = ("I2","I3","I4","I5","I6")
if set(fibres.keys()) != set(symbols):
    raise ArithmeticError(f"unexpected candidate1 reducible fibre keys: {tuple(fibres.keys())}")

pool = []

# Five-node sections: X is uniquely interpolated.
for subset in combinations(symbols,5):
    points = [(fibres[s][0],node_x(s)) for s in subset]
    N0s = R.lagrange_polynomial(points)
    F = R.one()
    for symbol in subset:
        F *= W-fibres[symbol][0]
    found = candidate_from_fixed_x(N0s,F,subset)
    pool.extend(found)
    print(
        "Q80MWBASISSEARCH|hits={}|mode=five|candidates={}|"
        "status=PASS_EXACT_HIGH_INCIDENCE_SUBSET".format(
            subset,len(found)
        ),
        flush=True,
    )

# Four-node sections: one X parameter c.
for subset in combinations(symbols,4):
    points = [(fibres[s][0],node_x(s)) for s in subset]
    N0s = R.lagrange_polynomial(points)
    assert N0s.degree() <= 3
    F = R.one()
    for symbol in subset:
        F *= W-fibres[symbol][0]
    F = R(F)
    assert F.degree() == 4

    CR = PolynomialRing(K,"c")
    cc = CR.gen()
    RUC = PolynomialRing(CR,"W")

    def lift_c(poly):
        poly = R(poly)
        return RUC([CR(value) for value in poly.list()])

    Xc = lift_c(N0s) + cc*lift_c(F)
    rhsc = RUC(Xc^3 + lift_c(A)*Xc + lift_c(B))
    qc, remc = rhsc.quo_rem(lift_c(F)^2)

    remainder_constraints = [
        CR(coefficient)
        for coefficient in remc.list()
        if CR(coefficient)
    ]

    constraint = None
    for equation in remainder_constraints:
        constraint = equation if constraint is None else constraint.gcd(equation)

    if constraint is None:
        qcoeff = [
            CR(qc[i]) if i <= qc.degree() else CR(0)
            for i in range(5)
        ]
        q0,q1,q2,q3,q4 = qcoeff
        square_eqs = []
        if q4:
            FC = CR.fraction_field()
            e1 = FC(q3)/(2*FC(q4))
            e0 = (FC(q2)/FC(q4)-e1^2)/2
            square_eqs = [
                CR((FC(q1)/FC(q4)-2*e1*e0).numerator()),
                CR((FC(q0)/FC(q4)-e0^2).numerator()),
            ]
        square_eqs = [f for f in square_eqs if f]
        for equation in square_eqs:
            constraint = equation if constraint is None else constraint.gcd(equation)

    if constraint is None or constraint.degree() <= 0:
        print(
            "Q80MWBASISSEARCH|hits={}|mode=four|constraint=none|"
            "status=SKIP_NO_FINITE_C_CONDITION".format(subset),
            flush=True,
        )
        continue

    constraint = CR(constraint.monic())
    c_values = []
    nonlinear = []
    for factor, unused_exp in constraint.factor():
        factor = CR(factor)
        if factor.degree() == 1:
            value = K(-factor[0]/factor[1])
            if value not in c_values:
                c_values.append(value)
        else:
            nonlinear.append((factor.degree(),str(factor)))

    found = []
    for c_value in c_values:
        Hx0 = R(N0s+c_value*F)
        found.extend(
            candidate_from_fixed_x(
                Hx0,F,subset,source_parameter=c_value
            )
        )
    pool.extend(found)

    print(
        "Q80MWBASISSEARCH|hits={}|mode=four|constraint_degree={}|"
        "linear_roots={}|nonlinear={}|candidates={}|"
        "status=PASS_EXACT_HIGH_INCIDENCE_SUBSET".format(
            subset,constraint.degree(),len(c_values),
            tuple(nonlinear),len(found)
        ),
        flush=True,
    )

unique = {}
for row in pool:
    key = (str(row["Hx"]),str(row["Hy"]))
    unique[key] = row
pool = sorted(
    unique.values(),
    key=lambda row:(str(row["Hx"]),str(row["Hy"]))
)

print(
    "Q80MWBASISPOOL|section_count={}|hitsets={}|"
    "status=PASS_EXACT_HIGH_INCIDENCE_SECTION_POOL".format(
        len(pool),
        tuple(row["hitset"] for row in pool),
    ),
    flush=True,
)

if len(pool) < 2:
    raise ArithmeticError(
        "high-incidence exact section pool has fewer than two sections"
    )

KF = R.fraction_field()
E = EllipticCurve(KF,[0,0,0,KF(A),KF(B)])

A73 = red_poly(A)
B73 = red_poly(B)
F73W = Rp.fraction_field()
E73 = EllipticCurve(F73W,[0,0,0,F73W(A73),F73W(B73)])

H73 = E73(F73W(X_seed),F73W(Y_seed))
assert not H73.is_zero()

exact_points = []
mod_points = []

for index,row in enumerate(pool):
    Hx0,Hy0 = row["Hx"],row["Hy"]
    Pexact = E(KF(Hx0),KF(Hy0))
    Pmod = E73(F73W(red_poly(Hx0)),F73W(red_poly(Hy0)))
    exact_points.append(Pexact)
    mod_points.append(Pmod)

    print(
        "Q80MWBASISCAND|index={}|hits={}|Hx={}|Hy={}|"
        "status=PASS_EXACT_HIGH_INCIDENCE_SECTION".format(
            index,row["hitset"],Hx0,Hy0
        ),
        flush=True,
    )

matches = []
for i in range(len(pool)):
    for jj in range(len(pool)):
        if i == jj:
            continue
        diff73 = mod_points[i]-mod_points[jj]
        sign = 0
        if diff73 == H73:
            sign = +1
        elif diff73 == -H73:
            sign = -1
        if not sign:
            continue

        diff = exact_points[i]-exact_points[jj]
        if diff.is_zero():
            continue
        hx,hy = diff.xy()
        hx = KF(hx)
        hy = KF(hy)

        polynomial = (
            hx.denominator() in K
            and hy.denominator() in K
        )
        hx_poly = None
        hy_poly = None
        if polynomial:
            hx_poly = R(hx)
            hy_poly = R(hy)
            polynomial = (
                hx_poly.degree() <= 4
                and hy_poly.degree() <= 6
            )

        matches.append({
            "i":i,
            "j":jj,
            "sign":sign,
            "exact":diff,
            "hx":hx,
            "hy":hy,
            "polynomial":bool(polynomial),
            "hx_poly":hx_poly,
            "hy_poly":hy_poly,
        })

        print(
            "Q80MWBASISDIFF|i={}|j={}|sign={}|polynomial={}|"
            "Hx={}|Hy={}|status=PASS_MOD73_P2_MINUS_P3_MATCH".format(
                i,jj,sign,int(polynomial),hx,hy
            ),
            flush=True,
        )

poly_matches = [row for row in matches if row["polynomial"]]
if not poly_matches:
    raise ArithmeticError(
        "found no exact high-incidence section difference reducing to the "
        "historical P2-P3 seed with P.O=0 polynomial coordinates"
    )

positive = [row for row in poly_matches if row["sign"] == +1]
selected = positive[0] if positive else poly_matches[0]

Hx = selected["hx_poly"]
Hy = selected["hy_poly"]
assert Hy^2 == Hx^3+A*Hx+B
assert red_poly(Hx) == X_seed
assert red_poly(Hy) == Y_seed

print(
    "Q80FINALHBYBASIS|i={}|j={}|pool={}|Hx={}|Hy={}|hits={}|"
    "status=PASS_EXACT_FINAL_HORIZONTAL_BY_MW_BASIS_DIFFERENCE".format(
        selected["i"],selected["j"],len(pool),
        Hx,Hy,exact_hitset(Hx,Hy),
    ),
    flush=True,
)

artifact = OUT / "q80-final-q6-char0-horizontal-via-basis-sections.json"
artifact.write_text(json.dumps({
    "status":"PASS_EXACT_FINAL_Q6_HORIZONTAL_VIA_BASIS_SECTIONS",
    "field":"QQ(sqrt(-3))",
    "method":(
        "exact >=4-node P.O=0 section recovery; identify a section "
        "difference by the transported historical GF73 P2-P3 marking"
    ),
    "pool_count":len(pool),
    "pool":[{
        "index":index,
        "Hx":str(row["Hx"]),
        "Hy":str(row["Hy"]),
        "hitset":list(row["hitset"]),
        "source_subset":list(row["source_subset"]),
        "parameter":None if row["parameter"] is None else str(row["parameter"]),
    } for index,row in enumerate(pool)],
    "selected":{
        "i":selected["i"],
        "j":selected["j"],
        "Hx":str(Hx),
        "Hy":str(Hy),
        "hitset":list(exact_hitset(Hx,Hy)),
        "reduction_Hx":str(red_poly(Hx)),
        "reduction_Hy":str(red_poly(Hy)),
    },
    "historical_regression":{
        "H_seed_X":str(X_seed),
        "H_seed_Y":str(Y_seed),
        "relation":"section_i - section_j = transported historical P2-P3 mod 73",
    },
    "next":"feed this exact H directly into the existing A4/A5 resolved RR compiler",
},indent=2,default=int)+"\n")

print(
    "Q80MWBASISFINAL|artifact={}|status="
    "PASS_EXACT_FINAL_Q6_HORIZONTAL_VIA_BASIS_SECTIONS".format(artifact),
    flush=True,
)
