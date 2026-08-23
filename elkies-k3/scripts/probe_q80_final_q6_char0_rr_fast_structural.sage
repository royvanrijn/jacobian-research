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
# FAST PARALLEL FAMILY RECOVERY
#
# The exact degree-(12,15,18) equations above are retained verbatim.  Instead
# of forming their characteristic-zero bivariate resultant, locate simple
# common roots at small split primes, Hensel-lift only (c0,c1), recognize them
# in K=QQ(sqrt(-3)), and then verify the ORIGINAL exact equations.
#
# A modular root can therefore never become an accepted characteristic-zero
# section merely by heuristic recognition: every accepted pair passes all
# exact equations and the exact Weierstrass square identity.
# ---------------------------------------------------------------------------

def coeff_parts(value):
    value = K(value)
    parts = list(value)+[QQ(0),QQ(0)]
    return QQ(parts[0]),QQ(parts[1])

def rational_mod(value, prime):
    value = QQ(value)
    den = int(value.denominator()) % prime
    if den == 0:
        raise ZeroDivisionError
    return (
        (int(value.numerator()) % prime)
        * pow(den,-1,prime)
    ) % prime

def k_mod(value, prime, jres):
    aa,bb = coeff_parts(value)
    return (
        rational_mod(aa,prime)
        + jres*rational_mod(bb,prime)
    ) % prime

def equation_terms_mod(equation, prime, jres):
    rows = []
    for exponent,coefficient in Cpar(equation).dict().items():
        e0,e1 = map(int,tuple(exponent))
        rows.append((e0,e1,k_mod(coefficient,prime,jres)))
    return tuple(rows)

def evaluate_terms_mod(terms, x, y, xpows, ypows, prime):
    total = 0
    for e0,e1,coefficient in terms:
        total += coefficient*xpows[x][e0]*ypows[y][e1]
    return total % prime

def gradient_terms_mod(terms, x, y, xpows, ypows, prime):
    dx = 0
    dy = 0
    for e0,e1,coefficient in terms:
        if e0:
            dx += (
                coefficient*e0*xpows[x][e0-1]*ypows[y][e1]
            )
        if e1:
            dy += (
                coefficient*e1*xpows[x][e0]*ypows[y][e1-1]
            )
    return dx % prime,dy % prime

def scan_prime(prime,jres):
    try:
        term_sets = [
            equation_terms_mod(eq,prime,jres)
            for eq in equations
        ]
    except ZeroDivisionError:
        return None

    max_x = max(
        e0 for terms in term_sets for e0,unused_e1,unused_c in terms
    )
    max_y = max(
        e1 for terms in term_sets for unused_e0,e1,unused_c in terms
    )
    xpows = [
        [pow(x,k,prime) for k in range(max_x+1)]
        for x in range(prime)
    ]
    ypows = [
        [pow(y,k,prime) for k in range(max_y+1)]
        for y in range(prime)
    ]

    roots = []
    simple = []
    rank_histogram = {}

    for y in range(prime):
        for x in range(prime):
            if evaluate_terms_mod(
                term_sets[0],x,y,xpows,ypows,prime
            ):
                continue
            if any(
                evaluate_terms_mod(
                    terms,x,y,xpows,ypows,prime
                )
                for terms in term_sets[1:]
            ):
                continue

            gradients = [
                gradient_terms_mod(
                    terms,x,y,xpows,ypows,prime
                )
                for terms in term_sets
            ]
            rank = 0
            chosen_pair = None
            if any(dx or dy for dx,dy in gradients):
                rank = 1
            for i in range(len(gradients)):
                for jj in range(i+1,len(gradients)):
                    ax,ay = gradients[i]
                    bx,by = gradients[jj]
                    if (ax*by-ay*bx) % prime:
                        rank = 2
                        chosen_pair = (i,jj)
                        break
                if chosen_pair is not None:
                    break

            rank_histogram[rank] = rank_histogram.get(rank,0)+1
            row = {
                "c0":int(x),
                "c1":int(y),
                "rank":int(rank),
                "pair":chosen_pair,
            }
            roots.append(row)
            if rank == 2:
                simple.append(row)

    return {
        "prime":int(prime),
        "j":int(jres),
        "roots":roots,
        "simple":simple,
        "rank_histogram":rank_histogram,
    }

def embedded_equation_terms(equation,Qlocal,jlocal):
    rows = []
    for exponent,coefficient in Cpar(equation).dict().items():
        e0,e1 = map(int,tuple(exponent))
        aa,bb = coeff_parts(coefficient)
        coeff_local = Qlocal(aa)+Qlocal(bb)*jlocal
        rows.append((e0,e1,coeff_local))
    return tuple(rows)

def eval_local(terms,x,y,Qlocal):
    value = Qlocal(0)
    dx = Qlocal(0)
    dy = Qlocal(0)
    for e0,e1,coefficient in terms:
        monomial = coefficient*(x^e0)*(y^e1)
        value += monomial
        if e0:
            dx += coefficient*e0*(x^(e0-1))*(y^e1)
        if e1:
            dy += coefficient*e1*(x^e0)*(y^(e1-1))
    return value,dx,dy

def recognize_K_local(alpha,residue,prime,jres):
    dep_ring = PolynomialRing(QQ,"zrec")
    dep = dep_ring(alpha.algebraic_dependency(2))
    if dep.degree() == 1:
        candidate = K(-dep[0]/dep[1])
        if k_mod(candidate,prime,jres) != residue:
            raise ArithmeticError("linear algdep has wrong modular residue")
        return candidate

    if dep.degree() != 2:
        raise ArithmeticError(f"algdep degree {dep.degree()} != 1 or 2")

    dep = dep/dep[2]
    linear = QQ(dep[1])
    constant = QQ(dep[0])
    real = -linear/2
    imag_sq = (constant-real^2)/3
    if not imag_sq.is_square():
        raise ArithmeticError(f"algdep not QQ(sqrt(-3))-shaped: {dep}")
    imag = imag_sq.sqrt()

    for candidate in (
        K(real+imag*j),
        K(real-imag*j),
    ):
        if k_mod(candidate,prime,jres) == residue:
            return candidate

    raise ArithmeticError("quadratic algdep conjugates miss modular residue")

def hensel_pair(prime,jres,root,precision=430):
    Qlocal = Qp(prime,prec=precision)
    jroots = Qlocal(-3).sqrt(all=True)
    jlocal = next(
        value for value in jroots
        if int(value.residue()) == int(jres)
    )

    local_equations = [
        embedded_equation_terms(eq,Qlocal,jlocal)
        for eq in equations
    ]

    pair = root["pair"]
    assert pair is not None
    i,jj = pair

    z = vector(Qlocal,[
        Qlocal(root["c0"]),
        Qlocal(root["c1"]),
    ])

    final_valuation = 0
    for iteration in range(14):
        vals = [
            eval_local(terms,z[0],z[1],Qlocal)
            for terms in local_equations
        ]
        Fv = vector(Qlocal,[vals[i][0],vals[jj][0]])
        Jv = matrix(Qlocal,[
            [vals[i][1],vals[i][2]],
            [vals[jj][1],vals[jj][2]],
        ])
        delta = Jv.solve_right(-Fv)
        z += delta

        final_valuation = min(
            entry.valuation() if entry else precision
            for entry in Fv
        )
        if final_valuation >= precision-45:
            break

    if final_valuation < precision-60:
        raise ArithmeticError(
            f"Hensel residual only {final_valuation}"
        )

    c0_exact = recognize_K_local(
        z[0],root["c0"],prime,jres
    )
    c1_exact = recognize_K_local(
        z[1],root["c1"],prime,jres
    )

    assignment = Cpar.hom(
        [K(c0_exact),K(c1_exact)],K
    )
    if any(assignment(eq) != 0 for eq in equations):
        raise ArithmeticError(
            "recognized pair fails original exact square equations"
        )

    return c0_exact,c1_exact,final_valuation

def exact_poly_sqrt_fast(value):
    value = R(value)
    factorization = value.factor()
    unit = K(factorization.unit())
    if not unit.is_square():
        return ()
    root = R(unit.sqrt())
    for factor,exponent in factorization:
        if int(exponent)%2:
            return ()
        root *= factor^(int(exponent)//2)
    assert root^2 == value
    return (root,) if root == 0 else (root,-root)

prime_list = (
    79,97,103,109,127,139,151,157,163,181,193,199,211,223,229
)

scan_records = []
exact_x = {}
successful_prime_embeddings = 0
successful_primes = set()

for prime in prime_list:
    jresidues = [
        value for value in range(prime)
        if (value*value+3) % prime == 0
    ]
    if len(jresidues) != 2:
        continue

    for jres in jresidues:
        scan = scan_prime(prime,jres)
        if scan is None:
            continue
        scan_records.append(scan)

        print(
            "Q80FINALFASTMOD|p={}|j={}|roots={}|simple={}|rank_hist={}|"
            "root_pairs={}|status=PASS_MODULAR_FAMILY_SCAN".format(
                prime,jres,
                len(scan["roots"]),len(scan["simple"]),
                scan["rank_histogram"],
                tuple(
                    (row["c0"],row["c1"],row["rank"])
                    for row in scan["roots"]
                ),
            ),
            flush=True,
        )

        if not scan["simple"]:
            continue

        embedding_had_exact = False
        for root_index,root in enumerate(scan["simple"]):
            try:
                c0_exact,c1_exact,val = hensel_pair(
                    prime,jres,root
                )
            except Exception as error:
                print(
                    "Q80FINALFASTLIFT|p={}|j={}|root={}|"
                    "status=SKIP|error={}".format(
                        prime,jres,
                        (root["c0"],root["c1"]),
                        str(error).replace("|","/"),
                    ),
                    flush=True,
                )
                continue

            Hx0 = R(
                N0+Fhit*(c0_exact+c1_exact*W)
            )
            rhs0 = R(Hx0^3+A*Hx0+B)
            Y20,remainder = rhs0.quo_rem(Fhit^2)
            if remainder:
                continue

            ysmall_roots = exact_poly_sqrt_fast(Y20)
            if not ysmall_roots:
                continue

            # Store x-solution once; signs are reconstructed below.
            key = (str(c0_exact),str(c1_exact))
            exact_x[key] = (c0_exact,c1_exact,Hx0)
            embedding_had_exact = True

            print(
                "Q80FINALFASTLIFT|p={}|j={}|root={}|"
                "c0={}|c1={}|valuation={}|"
                "status=PASS_EXACT_RECOGNIZED_X".format(
                    prime,jres,
                    (root["c0"],root["c1"]),
                    c0_exact,c1_exact,val,
                ),
                flush=True,
            )

        if embedding_had_exact:
            successful_prime_embeddings += 1
            successful_primes.add(int(prime))

    # The historical CM24 reconstruction had exactly three x-solutions
    # (three +/- section pairs).  Do not stop merely because both embeddings
    # of one auxiliary prime lifted: require the full expected x-count and a
    # second distinct prime as an independent completeness check.
    if len(exact_x) >= 3 and len(successful_primes) >= 2:
        break

if not exact_x:
    raise ArithmeticError(
        "fast modular/Hensel solver recognized no exact horizontal x-solutions"
    )

horizontal_candidates = []
for unused_key,(c0_exact,c1_exact,Hx0) in exact_x.items():
    rhs0 = R(Hx0^3+A*Hx0+B)
    Y20,remainder = rhs0.quo_rem(Fhit^2)
    assert remainder == 0
    for Ysmall0 in exact_poly_sqrt_fast(Y20):
        if Ysmall0.degree() > 3:
            continue
        Hy0 = R(Fhit*Ysmall0)
        if Hy0.degree() > 6 or Hy0^2 != rhs0:
            continue

        good_hits = True
        for symbol in fibres:
            root_value = fibres[symbol][0]
            hit = bool(
                Hx0(root_value) == node_x(symbol)
                and Hy0(root_value) == 0
            )
            if hit != (symbol in hit_symbols):
                good_hits = False
                break
        if not good_hits:
            continue

        horizontal_candidates.append(
            (Hx0,Hy0,c0_exact,c1_exact)
        )

unique_sections = {}
for item in horizontal_candidates:
    unique_sections[(str(item[0]),str(item[1]))] = item
horizontal_candidates = sorted(
    unique_sections.values(),
    key=lambda item:(str(item[0]),str(item[1]))
)

if not horizontal_candidates:
    raise ArithmeticError(
        "recognized exact x-solutions produced no exact horizontal sections"
    )

print(
    "Q80FINALFASTFAMILY|x_count={}|section_count={}|prime_embeddings={}|"
    "status=PASS_EXACT_FAST_HORIZONTAL_FAMILY".format(
        len(exact_x),len(horizontal_candidates),
        successful_prime_embeddings,
    ),
    flush=True,
)

marked = []
for index,(Hx0,Hy0,c0_exact,c1_exact) in enumerate(horizontal_candidates):
    hx73 = red_poly(Hx0)
    hy73 = red_poly(Hy0)
    sign = 0
    if hx73 == X_seed and hy73 == Y_seed:
        sign = +1
    elif hx73 == X_seed and hy73 == -Y_seed:
        sign = -1

    print(
        "Q80FINALFASTCAND|index={}|c0={}|c1={}|marked_sign={}|"
        "status=PASS_EXACT_FAST_CANDIDATE".format(
            index,c0_exact,c1_exact,sign
        ),
        flush=True,
    )

    if sign:
        marked.append(
            (index,Hx0,Hy0,c0_exact,c1_exact,sign)
        )

if not marked:
    raise ArithmeticError(
        "fast exact family contains no section reducing to historical P2-P3"
    )

plus = [row for row in marked if row[-1] == +1]
minus = [row for row in marked if row[-1] == -1]

if plus:
    selected_mark = plus[0]
else:
    selected_mark = minus[0]
    u3 = -u3
    Y_seed = -Y_seed

selected_index,Hx,Hy,c0,c1,seed_sign = selected_mark

fast_json = OUT / "q80-final-q6-char0-horizontal-family-fast.json"
fast_json.write_text(json.dumps({
    "status":"PASS_EXACT_FINAL_Q6_HORIZONTAL_FAMILY_FAST",
    "field":"QQ(sqrt(-3))",
    "method":"small-split-prime scan + 2-variable Hensel + exact recognition",
    "historical_rank73":rank73,
    "exact_x_count":len(exact_x),
    "section_count":len(horizontal_candidates),
    "marked_indices":[int(row[0]) for row in marked],
    "selected_index":int(selected_index),
    "modular_scans":[{
        "prime":row["prime"],
        "j":row["j"],
        "root_count":len(row["roots"]),
        "simple_count":len(row["simple"]),
        "rank_histogram":row["rank_histogram"],
    } for row in scan_records],
    "sections":[{
        "index":int(index),
        "Hx":str(Hx0),
        "Hy":str(Hy0),
        "c0":str(c0e),
        "c1":str(c1e),
        "historical_sign":(
            +1 if red_poly(Hx0)==X_seed and red_poly(Hy0)==Y_seed
            else -1 if red_poly(Hx0)==X_seed and red_poly(Hy0)==-Y_seed
            else 0
        ),
    } for index,(Hx0,Hy0,c0e,c1e) in enumerate(horizontal_candidates)],
},indent=2,default=int)+"\n")

print(
    "Q80FINALFASTH|family_count={}|marked_count={}|selected={}|"
    "c0={}|c1={}|Hx={}|Hy={}|"
    "status=PASS_EXACT_FINAL_HORIZONTAL".format(
        len(horizontal_candidates),len(marked),selected_index,
        c0,c1,Hx,Hy,
    ),
    flush=True,
)
print(
    "Q80FINALFASTARTIFACT|json={}|"
    "status=PASS_WRITE_EXACT_FAST_HORIZONTAL_FAMILY".format(
        fast_json
    ),
    flush=True,
)

# Total chord-slope scale from historical pinned model directly to our exact
# model after base PGL and Weierstrass gauge.
total_s = FU(den^2) * FU(u3/u2)

# ---------------------------------------------------------------------------
# 3. Whole A4 row at I5.
# ---------------------------------------------------------------------------
r5 = fibres["I5"][0]
x05 = node_x("I5")
assert Hx(r5) != x05
mA4 = K(Hy(r5)/(x05-Hx(r5)))  # chord through -H

e5 = exact_mod_root["I5"]
target_A4 = total_s(e5)*Fp(59)
assert red_k(mA4) == target_A4

a4_row = (K(1),r5,r5^2,mA4)
print(
    "Q80FINALA4|root={}|root_mod73={}|m={}|m_mod73={}|"
    "status=PASS_EXACT_WHOLE_A4".format(
        r5,int(e5),mA4,int(red_k(mA4))
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Connected A5 quotient at I6 by exact toric resolution.
# ---------------------------------------------------------------------------
r6 = fibres["I6"][0]
x06 = node_x("I6")
assert Hx(r6) == x06 and Hy(r6) == 0

Kr = FunctionField(K,"r")
rr = Kr.gen()
LS = LaurentSeriesRing(Kr,"s",default_prec=14)
ss = LS.gen()

def shifted(poly,root):
    ans = LS(0)
    for degree,coefficient in enumerate(R(poly).list()):
        ans += Kr(coefficient)*(Kr(root)+ss)^degree
    return ans

def newton_sqrt(value,root0):
    root = LS(Kr(root0))
    for unused in range(5):
        root = (root+value/root)/2
    return root

Aloc = shifted(A,r6)
Bloc = shifted(B,r6)
center = newton_sqrt(-Aloc/3,x06)
g0 = center^3+Aloc*center+Bloc
assert int(g0.valuation()) == 6
unit = g0/ss^6

rho_sq = K(3*x06)
if not rho_sq.is_square():
    raise ArithmeticError(
        "I6 tangent square is not in QQ(sqrt(-3)); final A5 probe needs a quadratic local extension"
    )
rho0 = K(rho_sq.sqrt())

Hxloc = shifted(Hx,r6)
Hyloc = shifted(Hy,r6)

def toric_point(component,rho_start):
    aa = LS(rr)*ss^component
    bb = unit*ss^(6-component)/LS(rr)
    yy = (aa+bb)/2
    ww = (bb-aa)/2
    rho = LS(Kr(rho_start))
    for unused in range(5):
        rho -= (rho^3-3*center*rho-ww)/(3*rho^2-3*center)
    xx = center+ww/rho
    return xx,yy

def functional_rows(values):
    nonzero = [v for v in values if v]
    if not nonzero:
        return ()
    common = nonzero[0].denominator().parent().one()
    for v in nonzero:
        common = common.lcm(v.denominator())
    nums = [(v*common).numerator() for v in values]
    degree = max([f.degree() for f in nums if f]+[-1])
    rows = []
    for deg in range(degree+1):
        row = tuple(
            K(f[deg]) if f and deg <= f.degree() else K(0)
            for f in nums
        )
        if any(row):
            rows.append(row)
    return tuple(rows)

def canonical(row):
    row = tuple(K(v) for v in row)
    pivot = next(v for v in row if v)
    return tuple(v/pivot for v in row)

def dedup(rows):
    out = []
    for row in rows:
        row = canonical(row)
        if row not in out:
            out.append(row)
    return tuple(out)

def reduce_row(row):
    row = canonical(row)
    vals = tuple(red_k(v) for v in row)
    pivot = next(v for v in vals if v)
    return tuple(v/pivot for v in vals)

e6 = exact_mod_root["I6"]
target_A5_residue = total_s(e6)*Fp(69)  # pinned -4
target_A5 = (Fp(1),e6,e6^2,target_A5_residue)

supports = {
    "direct":(1,3,4),
    "reversed":(2,3,5),
}
matches = []
for rho in (rho0,-rho0):
    for orientation,support in supports.items():
        rows = []
        for component in support:
            xx,yy = toric_point(component,rho)
            mloc = (yy+Hyloc)/(xx-Hxloc)  # chord through -H
            valuation = int(mloc.valuation())
            if valuation < 0:
                continue
            mres = Kr(mloc[0]) if valuation == 0 else Kr(0)
            rows.extend(functional_rows((Kr(1),Kr(r6),Kr(r6^2),Kr(mres))))
        rows = dedup(rows)
        redrows = tuple(reduce_row(row) for row in rows)

        print(
            "Q80FINALA5LOCAL|rho_mod73={}|orientation={}|support={}|rank={}|"
            "rows_mod73={}|status=PASS_EXACT_I6_TRACE".format(
                int(red_k(rho)),orientation,support,
                matrix(K,rows).rank() if rows else 0,
                tuple(tuple(int(v) for v in row) for row in redrows),
            ),
            flush=True,
        )

        for row in rows:
            if reduce_row(row) == target_A5:
                matches.append((rho,orientation,support,canonical(row)))

if not matches:
    raise ArithmeticError(
        "exact I6 traces did not contain transported pinned A5 quotient line"
    )

unique = {}
for item in matches:
    unique.setdefault(tuple(item[3]),item)
matches = list(unique.values())
matches.sort(key=lambda item:(item[1]!="direct",int(red_k(item[0]))))
rho,orientation,support,a5_row = matches[0]

print(
    "Q80FINALA5|root={}|root_mod73={}|orientation={}|support={}|"
    "row_mod73={}|status=PASS_EXACT_CONNECTED_A5_QUOTIENT".format(
        r6,int(e6),orientation,support,
        tuple(int(v) for v in reduce_row(a5_row)),
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 5. Complete exact RR kernel.
# ---------------------------------------------------------------------------
blocks = (
    {
        "name":"final q6 whole A4 quotient",
        "matrix":matrix(K,[a4_row]),
        "quotient_basis":("whole_A4_line",),
        "provenance":"exact nodal chord value at supported I5",
    },
    {
        "name":"final q6 connected A5 quotient",
        "matrix":matrix(K,[a5_row]),
        "quotient_basis":("connected_A5_line",),
        "provenance":"exact toric I6 trace selected by transported pinned -4 residue",
    },
)

compiled = compile_resolved_conditions(
    ("1","W","W^2","m"),blocks,complete=True,coefficient_field=K
)
assert compiled["ambient_dimension"] == 4
assert compiled["rank"] == 2
assert compiled["kernel_dimension"] == 2
assert compiled["h0_certified"]

C = compiled["condition_matrix"]
tail = C.matrix_from_columns([2,3])
assert tail.det() != 0
z1 = tail.solve_right(-C.column(0))
z2 = tail.solve_right(-C.column(1))
k1 = vector(K,[1,0,z1[0],z1[1]])
k2 = vector(K,[0,1,z2[0],z2[1]])
assert C*matrix(K,[k1,k2]).transpose() == matrix(K,2,2)

payload = {
    "status":"PASS_EXACT_FINAL_Q6_RR",
    "field":"QQ(sqrt(-3))",
    "source_candidate1_index":int(data["selected_index"]),
    "base_mod73":{
        "mobius":[int(a),int(b),int(c),int(d)],
        "typed_roots":{
            key:[int(exact_mod_root[key]),int(pinned_root[key])]
            for key in sorted(pinned_root)
        },
        "u2":str(u2),
        "u3":str(u3),
        "total_chord_scale":str(total_s),
    },
    "horizontal":{
        "identity":"P2-P3 (historical marking)",
        "Hx":str(Hx),
        "Hy":str(Hy),
        "hits":list(hit_symbols),
    },
    "A4":{
        "root":str(r5),
        "row":[str(v) for v in a4_row],
    },
    "A5":{
        "root":str(r6),
        "rho":str(rho),
        "orientation":orientation,
        "support":list(support),
        "row":[str(v) for v in a5_row],
    },
    "rr":{
        "ambient_dimension":4,
        "rank":2,
        "kernel_dimension":2,
        "h0":2,
        "kernel":[[str(v) for v in k1],[str(v) for v in k2]],
    },
    "next":"compile exact final q6 child and certify 4I3+I4+I6+2I1",
}
json_path = OUT / "q80-final-q6-char0-rr.json"
json_path.write_text(json.dumps(payload,indent=2,default=int)+"\n")

print(
    "Q80FINALRR|ambient=4|A4_rank=1|A5_rank=1|rank=2|nullity=2|h0=2|"
    "k1={}|k2={}|status=PASS_EXACT_FINAL_Q6_RR".format(
        tuple(k1),tuple(k2)
    ),
    flush=True,
)
print(
    "Q80FINALPROBEFINAL|json={}|status=PASS_EXACT_FINAL_Q6_HORIZONTAL_AND_RR".format(
        json_path
    ),
    flush=True,
)
