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
# Exact cubic-orbit recovery by modular Groebner + CRT.
#
# Historical GF(73) data has three x-solutions. Do NOT assume they lift as
# three K-rational points. Recover the whole K-scheme
#
#     g(c1)=0,    c0=h(c1),    deg(g)=3, deg(h)<3
#
# from small split primes. Exact substitution into the original three
# characteristic-zero equations is the final acceptance gate.
# ---------------------------------------------------------------------------

from sage.all import prime_range

def coeff_parts(value):
    value = K(value)
    parts = list(value)+[QQ(0),QQ(0)]
    return QQ(parts[0]),QQ(parts[1])

def qmod(value,p):
    value = QQ(value)
    den = int(value.denominator()) % p
    if den == 0:
        raise ZeroDivisionError
    return (
        (int(value.numerator()) % p)
        * pow(den,-1,p)
    ) % p

def kmod(value,p,jres):
    aa,bb = coeff_parts(value)
    return (qmod(aa,p)+jres*qmod(bb,p)) % p

def reduce_equation(eq,p,jres,Pq,xq,zq):
    ans = Pq.zero()
    for exponent,coefficient in Cpar(eq).dict().items():
        e0,e1 = map(int,tuple(exponent))
        ans += Pq.base_ring()(kmod(coefficient,p,jres))*xq^e0*zq^e1
    return ans

def groebner_graph(p,jres):
    Fq = GF(p)
    Pq = PolynomialRing(Fq,names=("x","z"),order="lex")
    xq,zq = Pq.gens()

    try:
        eqs = [
            reduce_equation(eq,p,jres,Pq,xq,zq)
            for eq in equations
        ]
    except ZeroDivisionError:
        return None

    if any(not eq for eq in eqs):
        return None

    I = Pq.ideal(eqs)
    G = tuple(Pq(g) for g in I.groebner_basis())

    Rz = PolynomialRing(Fq,"z")
    zz = Rz.gen()

    univariate = []
    linear_x = []
    for poly in G:
        dx = int(poly.degree(xq))
        dz = int(poly.degree(zq))
        if dx == 0 and dz > 0:
            val = Rz.zero()
            for (ex,ez),coefficient in poly.dict().items():
                if int(ex) != 0:
                    raise AssertionError
                val += Rz.base_ring()(coefficient)*zz^int(ez)
            univariate.append(val)
        elif dx == 1:
            linear_x.append(poly)

    if not univariate or not linear_x:
        return None

    g = min(
        univariate,
        key=lambda f:(int(f.degree()),len(f.dict()),str(f))
    ).monic()

    if g.degree() != 3:
        return {
            "status":"wrong_degree",
            "degree":int(g.degree()),
            "gb_size":len(G),
        }

    relation = min(
        linear_x,
        key=lambda f:(len(f.dict()),int(f.total_degree()),str(f))
    )

    a = Rz.zero()
    b = Rz.zero()
    bad = False
    for (ex,ez),coefficient in relation.dict().items():
        ex = int(ex); ez = int(ez)
        if ex == 1:
            a += Rz.base_ring()(coefficient)*zz^ez
        elif ex == 0:
            b += Rz.base_ring()(coefficient)*zz^ez
        else:
            bad = True
            break
    if bad or not a:
        return None

    d,scoef,unused = a.xgcd(g)
    if d.degree() != 0:
        return None
    ainv = (scoef/d[0]).quo_rem(g)[1]
    h = (-b*ainv).quo_rem(g)[1]

    for eq in eqs:
        image = Rz.zero()
        for (ex,ez),coefficient in eq.dict().items():
            image += (
                Rz.base_ring()(coefficient)
                * h^int(ex)
                * zz^int(ez)
            )
        if image.quo_rem(g)[1] != 0:
            return None

    return {
        "status":"ok",
        "g":g,
        "h":h,
        "gb_size":len(G),
        "g_factor_degrees":tuple(
            sorted(int(f.degree()) for f,unused_e in g.factor())
        ),
    }

def crt_update(old_residue,old_modulus,new_residue,p):
    if old_modulus == 1:
        return int(new_residue)%p, p
    correction = (
        (int(new_residue)-old_residue)
        * pow(old_modulus % p,-1,p)
    ) % p
    value = old_residue + old_modulus*correction
    modulus = old_modulus*p
    return value % modulus, modulus

names = ("g0","g1","g2","h0","h1","h2")
crt_state = {
    (name,part):(0,1)
    for name in names
    for part in ("a","b")
}

def coeff_list3(poly):
    return [
        poly[i] if i <= poly.degree() else poly.base_ring()(0)
        for i in range(3)
    ]

def combine_embeddings(v_plus,v_minus,p,j_plus,j_minus):
    denom = (j_plus-j_minus) % p
    b = ((int(v_plus)-int(v_minus))*pow(denom,-1,p)) % p
    a = (int(v_plus)-b*j_plus) % p
    return a,b

def try_reconstruct():
    values = {}
    common_modulus = None
    for name in names:
        parts = []
        for part in ("a","b"):
            residue,modulus = crt_state[(name,part)]
            if common_modulus is None:
                common_modulus = modulus
            if modulus != common_modulus or modulus <= 1:
                return None
            try:
                rr = QQ(ZZ(residue).rational_reconstruction(ZZ(modulus)))
            except Exception:
                return None
            parts.append(rr)
        values[name] = K(parts[0]+parts[1]*j)

    ZK = PolynomialRing(K,"z")
    zK = ZK.gen()
    gK = ZK([
        values["g0"],values["g1"],values["g2"],K(1)
    ])
    hK = ZK([
        values["h0"],values["h1"],values["h2"]
    ])

    hom = Cpar.hom([hK,zK],ZK)
    for eq in equations:
        image = ZK(hom(eq))
        if image.quo_rem(gK)[1] != 0:
            return None

    return gK,hK,common_modulus

good_primes = []
bad_primes = []
recovered = None

for p in prime_range(7,600):
    p = int(p)
    if p == 73 or p % 3 != 1:
        continue

    roots_j = [
        value for value in range(p)
        if (value*value+3) % p == 0
    ]
    if len(roots_j) != 2:
        continue
    j0,j1 = roots_j

    plus = groebner_graph(p,j0)
    minus = groebner_graph(p,j1)

    if (
        plus is None or minus is None
        or plus.get("status") != "ok"
        or minus.get("status") != "ok"
    ):
        bad_primes.append(p)
        if plus is not None or minus is not None:
            print(
                "Q80FINALCUBICMOD|p={}|plus={}|minus={}|"
                "status=SKIP_BAD_MODULAR_GRAPH".format(
                    p,plus,minus
                ),
                flush=True,
            )
        continue

    gp = coeff_list3(plus["g"])
    gm = coeff_list3(minus["g"])
    hp = coeff_list3(plus["h"])
    hm = coeff_list3(minus["h"])

    residues = {}
    for name,vp,vm in zip(("g0","g1","g2"),gp,gm):
        residues[name] = combine_embeddings(vp,vm,p,j0,j1)
    for name,vp,vm in zip(("h0","h1","h2"),hp,hm):
        residues[name] = combine_embeddings(vp,vm,p,j0,j1)

    for name,(ares,bres) in residues.items():
        old,mod = crt_state[(name,"a")]
        crt_state[(name,"a")] = crt_update(old,mod,ares,p)
        old,mod = crt_state[(name,"b")]
        crt_state[(name,"b")] = crt_update(old,mod,bres,p)

    good_primes.append(p)

    current_modulus = crt_state[("g0","a")][1]
    print(
        "Q80FINALCUBICMOD|p={}|j={}|g_factors_plus={}|g_factors_minus={}|"
        "gb_sizes={},{}|good_primes={}|crt_bits={}|"
        "status=PASS_CUBIC_GRAPH_MODP".format(
            p,(j0,j1),
            plus["g_factor_degrees"],minus["g_factor_degrees"],
            plus["gb_size"],minus["gb_size"],
            len(good_primes),ZZ(current_modulus).nbits(),
        ),
        flush=True,
    )

    recovered = try_reconstruct()
    if recovered is not None:
        break

if recovered is None:
    raise ArithmeticError(
        "modular cubic graph was found, but CRT modulus through p<600 "
        "was insufficient for exact K reconstruction"
    )

gK,hK,crt_modulus = recovered

print(
    "Q80FINALCUBICCRT|good_primes={}|bad_primes={}|crt_bits={}|"
    "g={}|h={}|factorization={}|"
    "status=PASS_EXACT_CUBIC_HORIZONTAL_X_ORBIT".format(
        tuple(good_primes),tuple(bad_primes),
        ZZ(crt_modulus).nbits(),gK,hK,gK.factor()
    ),
    flush=True,
)

def red_K_to_73(value):
    return Fp(red_k(value))

g73 = Rp([red_K_to_73(c) for c in gK.list()])
h73 = Rp([red_K_to_73(c) for c in hK.list()])
seed_c0 = Fp(seed[0])
seed_c1 = Fp(seed[1])

assert g73(seed_c1) == 0
assert h73(seed_c1) == seed_c0

print(
    "Q80FINALCUBIC73|seed_c0={}|seed_c1={}|g73={}|h73={}|"
    "status=PASS_HISTORICAL_SEED_IN_EXACT_CUBIC_ORBIT".format(
        int(seed_c0),int(seed_c1),g73,h73
    ),
    flush=True,
)

L = K.extension(gK,"alpha")
alpha = L.gen()
c1L = alpha
c0L = L(hK(alpha))

RL = PolynomialRing(L,"W")
WL = RL.gen()

def lift_to_L(poly):
    poly = R(poly)
    return RL([L(c) for c in poly.list()])

AL = lift_to_L(A)
BL = lift_to_L(B)
N0L = lift_to_L(N0)
FhitL = lift_to_L(Fhit)

HxL = RL(N0L+FhitL*(c0L+c1L*WL))
rhsL = RL(HxL^3+AL*HxL+BL)
Y2L,remL = rhsL.quo_rem(FhitL^2)
assert remL == 0
assert Y2L.degree() <= 6

print(
    "Q80FINALCUBICFIELD|relative_degree={}|absolute_degree={}|"
    "Hx_degree={}|Y2_degree={}|"
    "status=PASS_EXACT_HORIZONTAL_X_OVER_CUBIC_FIELD".format(
        L.relative_degree(),L.absolute_degree(),
        HxL.degree(),Y2L.degree()
    ),
    flush=True,
)

artifact = OUT / "q80-final-q6-char0-cubic-horizontal-orbit.json"
artifact.write_text(json.dumps({
    "status":"PASS_EXACT_FINAL_Q6_CUBIC_HORIZONTAL_X_ORBIT",
    "base_field":"QQ(sqrt(-3))",
    "g_c1":str(gK),
    "c0_of_c1":str(hK),
    "g_factorization":str(gK.factor()),
    "relative_degree":int(L.relative_degree()),
    "absolute_degree":int(L.absolute_degree()),
    "good_primes":good_primes,
    "bad_primes":bad_primes,
    "crt_modulus_bits":int(ZZ(crt_modulus).nbits()),
    "historical_mod73_seed":[int(seed_c0),int(seed_c1)],
    "Hx_over_relative_field":str(HxL),
    "Y2_over_relative_field":str(Y2L),
    "next":(
        "recover signed y over this relative field (or quadratic extension), "
        "then apply exact A4/A5 quotient conditions"
    ),
},indent=2,default=int)+"\n")

print(
    "Q80FINALCUBICFINAL|json={}|"
    "status=PASS_EXACT_FINAL_Q6_CUBIC_HORIZONTAL_ORBIT".format(
        artifact
    ),
    flush=True,
)
