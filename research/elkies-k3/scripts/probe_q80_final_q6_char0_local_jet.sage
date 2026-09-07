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
# 2. Exact final-horizontal FAMILY over K = QQ(sqrt(-3)).
#
# Historical construction order matters here. The three node hits do NOT
# isolate the final horizontal over GF(73): the original CM24 calculation
# found six polynomial sections (three +/- pairs) before the A4+A5 quotient
# conditions selected the target pair. Accordingly we solve the complete
# characteristic-zero family first and use p=73 only afterwards as the
# marking/regression prime.
#
# P.O=0 gives deg X<=4, deg Y<=6. The horizontal hits exactly I3,I4,I6.
# Interpolating the three node x-values leaves precisely two x-parameters:
#
#     X = N0 + F*(c0+c1*W),       F=(W-r3)(W-r4)(W-r6).
#
# Since Y vanishes at all three nodes,
#
#     Y = F*G,  deg G<=3.
#
# Hence
#
#     (X^3+A*X+B)/F^2 = G^2.
#
# We eliminate G structurally. For a degree-six polynomial
#
#     q = q6 W^6 + ... + q0 = q6*(W^3+e2 W^2+e1 W+e0)^2
#
# the top four coefficients determine e2,e1,e0 rationally in c0,c1;
# the remaining three coefficients are exact polynomial conditions.
# Together with any F^2-remainder equations this gives a zero-dimensional
# two-variable problem.
# ---------------------------------------------------------------------------

# Keep the rank-5 historical diagnostic: it is expected and is NOT used for
# lifting.
Xmod = N0p + Fhitp*(seed[0]+seed[1]*wp)
Ymod = Fhitp*sum(seed[2+i]*wp^i for i in range(4))
Ebase = -3*Xmod^2-red_poly(A)
columns73 = [
    Ebase*Fhitp,
    Ebase*Fhitp*wp,
] + [
    2*Ymod*Fhitp*wp^i for i in range(4)
]
J73 = matrix(Fp,13,6,lambda r,col:
    columns73[col][r] if r <= columns73[col].degree() else 0
)
rank73 = int(J73.rank())
print(
    "Q80FINALFAMILY73|rank={}|nullity={}|"
    "status=PASS_HISTORICAL_NONISOLATED_SEED".format(rank73,6-rank73),
    flush=True,
)

Cpar = PolynomialRing(K, names=("c0","c1"), order="lex")
c0v,c1v = Cpar.gens()

# ---------------------------------------------------------------------------
# Structural quotient: avoid expanding (N0+F*L)^3 in K[c0,c1,W].
#
# P(N0) is divisible by F^2 and P_x(N0)=3*N0^2+A is divisible by F at
# the three I3/I4/I6 double roots.  Therefore
#
#   Q = (X^3+A*X+B)/F^2
#     = Q0 + Q1*L + 3*N0*L^2 + F*L^3,
#
# with Q0,Q1 univariate over K and L=c0+c1*W.
# ---------------------------------------------------------------------------
P0 = R(N0^3 + A*N0 + B)
Q0,rem0 = P0.quo_rem(Fhit^2)
if rem0:
    raise ArithmeticError(
        "structural final-horizontal identity failed: F^2 does not divide P(N0)"
    )

Px0 = R(3*N0^2 + A)
Q1,rem1 = Px0.quo_rem(Fhit)
if rem1:
    raise ArithmeticError(
        "structural final-horizontal identity failed: F does not divide P_x(N0)"
    )

Q2 = R(3*N0)
Q3 = R(Fhit)

print(
    "Q80FINALFAMILYSTRUCT|degQ0={}|degQ1={}|degQ2={}|degQ3={}|"
    "status=PASS_STRUCTURAL_PREQUOTIENT".format(
        Q0.degree(),Q1.degree(),Q2.degree(),Q3.degree()
    ),
    flush=True,
)

def rcoeff(poly,index):
    poly = R(poly)
    if index < 0 or index > poly.degree():
        return K(0)
    return K(poly[index])

# Coefficients of L, L^2, L^3 in W, kept directly in K[c0,c1].
L1 = (c0v,c1v)
L2 = (c0v^2, 2*c0v*c1v, c1v^2)
L3 = (
    c0v^3,
    3*c0v^2*c1v,
    3*c0v*c1v^2,
    c1v^3,
)

qcoeff = []
for degree in range(7):
    value = Cpar(rcoeff(Q0,degree))
    for shift,monomial in enumerate(L1):
        value += Cpar(rcoeff(Q1,degree-shift))*monomial
    for shift,monomial in enumerate(L2):
        value += Cpar(rcoeff(Q2,degree-shift))*monomial
    for shift,monomial in enumerate(L3):
        value += Cpar(rcoeff(Q3,degree-shift))*monomial
    qcoeff.append(Cpar(value))

q0c,q1c,q2c,q3c,q4c,q5c,q6c = qcoeff
if not q6c:
    raise ArithmeticError("final horizontal square quotient has q6=0 identically")

# ---------------------------------------------------------------------------
# Direct polynomial square conditions.
#
# If
#   q/q6 = (W^3 + e2 W^2 + e1 W + e0)^2,
#
# then with
#   A1 = 4*q6*q4 - q5^2,
#   B1 = 8*q6^2*q3 - A1*q5,
#
# clearing denominators in the remaining W^2,W^1,W^0 equations gives:
#
#   64 q6^3 q2 - A1^2 - 4 q5 B1 = 0
#   64 q6^4 q1 - A1 B1          = 0
#  256 q6^5 q0 - B1^2           = 0.
#
# This is exactly the old fraction-field recurrence, but creates no rational
# functions and requires no multivariate numerator normalization.
# ---------------------------------------------------------------------------
A1sq = Cpar(4*q6c*q4c - q5c^2)
B1sq = Cpar(8*q6c^2*q3c - A1sq*q5c)

equations = [
    Cpar(64*q6c^3*q2c - A1sq^2 - 4*q5c*B1sq),
    Cpar(64*q6c^4*q1c - A1sq*B1sq),
    Cpar(256*q6c^5*q0c - B1sq^2),
]

# Primitive scalar normalization only; do not invoke fraction fields.
normalized = []
for equation in equations:
    if not equation:
        raise ArithmeticError("one structural square equation vanished identically")
    normalized.append(
        Cpar(equation / equation.leading_coefficient())
    )
equations = normalized
equations.sort(key=lambda f:(int(f.total_degree()),len(f.dict()),str(f)))

print(
    "Q80FINALFAMILYEQ|equations={}|degrees={}|terms={}|"
    "status=PASS_EXACT_TWO_PARAMETER_SYSTEM_STRUCTURAL".format(
        len(equations),
        tuple(int(f.total_degree()) for f in equations),
        tuple(len(f.dict()) for f in equations),
    ),
    flush=True,
)


# ---------------------------------------------------------------------------
# Local primary-component diagnostic at the historical p=73 seed.
#
# The reduced (c0,c1) Jacobian has rank one.  Put t along its tangent line
# and n transverse to it.  One equation has unit derivative in n, so solve
# it uniquely for n=n(t) in F_73[[t]], then inspect the two remaining
# obstruction series.  Their first nonzero orders give the local
# intersection multiplicity seen at the historical seed.
# ---------------------------------------------------------------------------

from sage.all import PowerSeriesRing

seed_c0 = Fp(seed[0])
seed_c1 = Fp(seed[1])

def eval_mod73(poly,x,y):
    total = Fp(0)
    for exponent,coefficient in Cpar(poly).dict().items():
        e0,e1 = map(int,tuple(exponent))
        total += red_k(coefficient)*(x^e0)*(y^e1)
    return Fp(total)

J2 = matrix(Fp,[
    [
        eval_mod73(eq.derivative(c0v),seed_c0,seed_c1),
        eval_mod73(eq.derivative(c1v),seed_c0,seed_c1),
    ]
    for eq in equations
])
assert J2.rank() == 1

ker = J2.right_kernel().basis()[0]
if ker[1] == 0:
    raise ArithmeticError(
        "historical tangent is vertical in c1; local diagnostic expects c1 parameter"
    )
ker = vector(Fp,[ker[0]/ker[1],Fp(1)])
tangent_c0 = Fp(ker[0])

# Choose an equation with nonzero normal derivative d/dc0.
normal_equation = next(
    i for i in range(len(equations))
    if J2[i,0] != 0
)
other_equations = tuple(
    i for i in range(len(equations))
    if i != normal_equation
)

print(
    "Q80FINALLOCALJETSETUP|seed=({}, {})|jacobian={}|"
    "tangent=({},1)|normal_equation={}|others={}|"
    "status=PASS_LOCAL_COORDINATES".format(
        int(seed_c0),int(seed_c1),
        tuple(tuple(int(v) for v in row) for row in J2.rows()),
        int(tangent_c0),normal_equation,other_equations,
    ),
    flush=True,
)

PREC = 80
PS = PowerSeriesRing(Fp,"t",default_prec=PREC)
tt = PS.gen()

def eval_series(poly,c0s,c1s):
    total = PS(0)
    for exponent,coefficient in Cpar(poly).dict().items():
        e0,e1 = map(int,tuple(exponent))
        total += PS(red_k(coefficient))*(c0s^e0)*(c1s^e1)
    return total.add_bigoh(PREC)

def derivative_series(poly,variable,c0s,c1s):
    return eval_series(poly.derivative(variable),c0s,c1s)

# Coordinates:
#   c1 = seed_c1 + t
#   c0 = seed_c0 + tangent_c0*t + n(t)
#
# By construction n has no constant or linear term.
c1s = PS(seed_c1)+tt
nseries = PS(0)

for iteration in range(8):
    c0s = PS(seed_c0)+PS(tangent_c0)*tt+nseries
    f = eval_series(equations[normal_equation],c0s,c1s)
    df = derivative_series(
        equations[normal_equation],c0v,c0s,c1s
    )
    if df[0] == 0:
        raise ArithmeticError(
            "normal derivative lost unit constant during local Newton solve"
        )
    correction = (f/df).add_bigoh(PREC)
    nseries = (nseries-correction).add_bigoh(PREC)
    residual = eval_series(
        equations[normal_equation],
        PS(seed_c0)+PS(tangent_c0)*tt+nseries,
        c1s,
    )
    valuation = (
        PREC if residual == 0 else int(residual.valuation())
    )
    print(
        "Q80FINALLOCALJETNEWTON|iteration={}|normal_residual_order={}|"
        "n_order={}".format(
            iteration,valuation,
            PREC if nseries == 0 else int(nseries.valuation())
        ),
        flush=True,
    )
    if valuation >= PREC-5:
        break

c0s = PS(seed_c0)+PS(tangent_c0)*tt+nseries
normal_residual = eval_series(
    equations[normal_equation],c0s,c1s
)
assert normal_residual == 0 or normal_residual.valuation() >= PREC-5

obstructions = []
for index in other_equations:
    series = eval_series(equations[index],c0s,c1s)
    order = PREC if series == 0 else int(series.valuation())
    lead = None if series == 0 else int(series[order])
    obstructions.append((index,order,lead,series))
    print(
        "Q80FINALLOCALJETOBS|equation={}|order={}|leading={}|"
        "series={}|status=PASS_LOCAL_OBSTRUCTION".format(
            index,order,lead,
            series.add_bigoh(min(PREC,order+12)) if series != 0 else 0,
        ),
        flush=True,
    )

finite_orders = [
    order for unused_i,order,unused_l,unused_s in obstructions
    if order < PREC
]
if not finite_orders:
    raise ArithmeticError(
        "both residual equations vanish to precision 80; local component "
        "is higher-dimensional or multiplicity exceeds diagnostic precision"
    )

local_order = min(finite_orders)

# Compare leading residual lines after dividing by the minimal t-power.
leading_rows = []
for index,order,lead,series in obstructions:
    if order == local_order:
        leading_rows.append((index,lead))

print(
    "Q80FINALLOCALJETFINAL|tangent=({},1)|normal_equation={}|"
    "obstruction_orders={}|local_order={}|leading_rows={}|"
    "n_series={}|status=PASS_LOCAL_PRIMARY_JET".format(
        int(tangent_c0),normal_equation,
        tuple((i,o,l) for i,o,l,unused_s in obstructions),
        local_order,tuple(leading_rows),
        nseries.add_bigoh(min(PREC,local_order+12)),
    ),
    flush=True,
)
