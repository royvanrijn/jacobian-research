#!/usr/bin/env sage
"""
Probe the exact characteristic-zero q4 candidate1 RR pencil on the certified
q4_6855 child.

This stage intentionally stops before child construction.  It checks:

  * the exact finite-minimal 6855 model and its GF73 rational gauge bridge;
  * deterministic recovery of the chosen H=-P3 section from its four
    singular-point hits;
  * exact identification of H with the pinned modular representative;
  * the supported whole-A3 chord quotient at the I4 reducing to U=66;
  * the D4 outer-complement c=0 condition at the I0* reducing to U=48;
  * the resulting 4 -> 2 exact RR kernel.

No p-adics, algdep, or toric scan.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, GF, PolynomialRing, QuadraticField, matrix, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ART = ROOT / "artifacts" / "generated-results" / "q80-q4-6855-char0-resolved-rr.json"

if not ART.exists():
    raise SystemExit(f"missing q4_6855 artifact: {ART}")

data = json.loads(ART.read_text())
assert data["status"] == "PASS_EXACT_Q4_6855_RESOLVED_RR"

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()
R = PolynomialRing(K, "U")
U = R.gen()
KF = R.fraction_field()
Uf = KF(U)

child = data["child"]
A = R(sage_eval(child["minimal_A"], locals={"j":j, "U":U}))
B = R(sage_eval(child["minimal_B"], locals={"j":j, "U":U}))
s_min = KF(sage_eval(child["finite_scaling_unit"], locals={"j":j, "U":Uf}))

Fp = GF(73)
JMOD = Fp(17)
assert JMOD^2 == Fp(-3)

def red_q(c):
    c = QQ(c)
    return Fp(c.numerator()) / Fp(c.denominator())

def red_k(c):
    c = K(c)
    cc = list(c) + [QQ(0), QQ(0)]
    return red_q(cc[0]) + JMOD*red_q(cc[1])

Rp = PolynomialRing(Fp, "U")
up = Rp.gen()
FU = Rp.fraction_field()

def red_poly(f):
    f = R(f)
    return Rp([red_k(c) for c in f.list()])

def red_frac(f):
    f = KF(f)
    num = red_poly(f.numerator())
    den = red_poly(f.denominator())
    if den == 0:
        raise ZeroDivisionError("denominator dies modulo 73")
    return FU(num)/FU(den)

# Pinned minimal q4_6855 parent from the candidate1 certificate.
A_pin = Rp([35,65,22,17,8,8,34,31,33,10,71,52,29,7,35,17])
B_pin = Rp([3,39,8,24,18,29,17,35,8,38,58,18,59,60,20,25,19,17,53,69,65,46,5,21,49])

g73 = Fp(int(data["parent"]["gf73_gauge"]))
raw_bridge = FU(g73^2) / FU((up+25)^4)
bridge = FU(raw_bridge * red_frac(s_min))

assert red_frac(KF(A)) == bridge^4 * FU(A_pin)
assert red_frac(KF(B)) == bridge^6 * FU(B_pin)

print(
    "Q80CAND1PARENT|degA={}|degB={}|bridge={}|"
    "status=PASS_EXACT_6855_PARENT_BRIDGE".format(
        A.degree(), B.degree(), bridge
    ),
    flush=True,
)

# Recover the exact reducible fibres and match them by their GF73 roots.
fibres = {}
for item in child["finite_fibres"]:
    symbol = str(item["kodaira"])
    if symbol not in ("I0*","I4","I2") or int(item["degree"]) != 1:
        continue
    f = R(sage_eval(item["factor"], locals={"j":j, "U":U}))
    assert f.degree() == 1
    root = K(-f[0]/f[1])
    fibres[int(red_k(root))] = (symbol, root, f)

expected = {48:"I0*",49:"I0*",2:"I4",66:"I4",32:"I2"}
assert {res:symbol for res,(symbol,unused_root,unused_f) in fibres.items()} == expected

def node_x(residue):
    symbol, root, unused_factor = fibres[residue]
    ar, br = K(A(root)), K(B(root))
    if symbol == "I0*":
        # In the minimal short model the D4 Weierstrass singular point is
        # (x,y)=(0,0); A and B vanish to orders >=2,>=3.
        assert ar == 0 and br == 0
        return K(0)
    assert ar != 0
    x0 = K(-3*br/(2*ar))
    assert x0^3 + ar*x0 + br == 0
    assert 3*x0^2 + ar == 0
    return x0

# The candidate1 horizontal is P.O=0, hence deg X<=4 and deg Y<=6.
# It hits I4@2, I2@32 and both D4 singular points @48,49.
#
# The previous probe incorrectly forced deg X<=3 because the *mod-73 pinned*
# representative has no U^4 term.  In characteristic zero there is one extra
# coefficient:
#
#   X = N0 + c*Fhit,  deg Fhit=4.
#
# At an I0* singular point a nonidentity section has ord(x)>=1, ord(y)>=2.
# Therefore the degree-six y-coordinate is forced up to scalar:
#
#   Y = lambda*(U-r2)*(U-r32)*(U-r48)^2*(U-r49)^2.
#
# This reduces the exact reconstruction to one variable c.  Comparing
# X^3+A*X+B with lambda^2*Ybase^2 gives a univariate condition in c.
hit_residues = (2,32,48,49)
points = [(fibres[res][1], node_x(res)) for res in hit_residues]
N0 = R.lagrange_polynomial(points)
assert N0.degree() <= 3

Fhit = R.one()
for res in hit_residues:
    Fhit *= U-fibres[res][1]
assert Fhit.degree() == 4 and Fhit.leading_coefficient() == 1

Ybase = (
    (U-fibres[2][1])
    * (U-fibres[32][1])
    * (U-fibres[48][1])^2
    * (U-fibres[49][1])^2
)
Ybase = R(Ybase)
assert Ybase.degree() == 6 and Ybase.leading_coefficient() == 1

CR = PolynomialRing(K, "c")
cc = CR.gen()
RUC = PolynomialRing(CR, "U")
uc = RUC.gen()

def lift_c(poly):
    poly = R(poly)
    return RUC([CR(value) for value in poly.list()])

Xc = lift_c(N0) + cc*lift_c(Fhit)
Ac = lift_c(A)
Bc = lift_c(B)
Ybc = lift_c(Ybase)

rhs_c = RUC(Xc^3 + Ac*Xc + Bc)
assert rhs_c.degree() <= 12

# Ybase is monic degree 6, so lambda^2 is the U^12 coefficient.
mu_c = CR(rhs_c[12])
difference = RUC(rhs_c - mu_c*Ybc^2)

constraint = None
for coeff in difference.list():
    coeff = CR(coeff)
    if not coeff:
        continue
    constraint = coeff if constraint is None else constraint.gcd(coeff)

if constraint is None or constraint.degree() <= 0:
    raise ArithmeticError(
        "candidate1 one-variable node reconstruction produced no finite c condition"
    )

constraint = CR(constraint.monic())
print(
    "Q80CAND1HSEARCH|constraint_degree={}|factorization={}|"
    "status=PASS_EXACT_ONE_VARIABLE_REDUCTION".format(
        constraint.degree(), constraint.factor()
    ),
    flush=True,
)

c_values = []
for factor,unused_exp in constraint.factor():
    factor = CR(factor)
    if factor.degree() != 1:
        continue
    value = K(-factor[0]/factor[1])
    if value not in c_values:
        c_values.append(value)

if not c_values:
    raise ArithmeticError(
        "candidate1 horizontal condition has no QQ(sqrt(-3))-rational roots"
    )

horizontal_candidates = []
for c_value in c_values:
    Hx0 = R(N0 + c_value*Fhit)
    rhs0 = R(Hx0^3 + A*Hx0 + B)
    mu = K(rhs0[12])  # Ybase is monic.
    if not mu or not mu.is_square():
        continue
    lam0 = K(mu.sqrt())

    for lam in (lam0,-lam0):
        Hy0 = R(lam*Ybase)
        if Hy0^2 != rhs0:
            continue

        if Hx0.degree() > 4 or Hy0.degree() > 6:
            continue

        good = True
        for residue in (2,32,48,49,66):
            root = fibres[residue][1]
            hit = bool(Hx0(root) == node_x(residue) and Hy0(root) == 0)
            if hit != (residue in hit_residues):
                good = False
                break
        if not good:
            continue

        d4vals = {}
        for residue in (48,49):
            factor = fibres[residue][2]
            vx = int(Hx0.valuation(factor))
            vy = int(Hy0.valuation(factor))
            if vx < 1 or vy < 2:
                good = False
                break
            d4vals[residue] = (vx,vy)
        if not good:
            continue

        horizontal_candidates.append((Hx0,Hy0,c_value,lam,d4vals))

unique = {}
for item in horizontal_candidates:
    unique[(str(item[0]),str(item[1]))] = item
horizontal_candidates = sorted(
    unique.values(), key=lambda item: (str(item[0]),str(item[1]))
)
if not horizontal_candidates:
    raise ArithmeticError("no exact candidate1 P.O=0 horizontal survived")

print(
    "Q80CAND1HSEARCH|c_roots={}|horizontal_count={}|"
    "status=PASS_EXACT_HORIZONTAL_FAMILY".format(
        len(c_values),len(horizontal_candidates)
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# Exact local quotient blocks for every surviving horizontal.
# ---------------------------------------------------------------------------
r66 = fibres[66][1]
x066 = node_x(66)
r48 = fibres[48][1]

records = []
for index,(Hx,Hy,c_value,lam,d4vals) in enumerate(horizontal_candidates):
    assert Hx(r66) != x066
    cA3 = K(Hy(r66)/(x066-Hx(r66)))

    a3_row = (K(1),r66,r66^2,cA3)
    d4_row = (K(1),r48,r48^2,K(0))
    C = matrix(K,[a3_row,d4_row])
    assert C.rank() == 2

    blocks = (
        {
            "name":"candidate1 whole A3 connected quotient",
            "matrix":matrix(K,[a3_row]),
            "quotient_basis":("connected_A3_line",),
            "provenance":(
                "exact nodal chord value at supported I4; "
                "whole A3 coefficients (-1,-1,-1)"
            ),
        },
        {
            "name":"candidate1 D4 outer-complement quotient",
            "matrix":matrix(K,[d4_row]),
            "quotient_basis":("D4_outer_complement",),
            "provenance":(
                "exact I0* ramified-chart valuation input; structural "
                "outer-complement residue c=0 for (-1,0,-1,-1)"
            ),
        },
    )

    compiled = compile_resolved_conditions(
        ("1","U","U^2","m"), blocks, complete=True, coefficient_field=K
    )
    assert compiled["ambient_dimension"] == 4
    assert compiled["rank"] == 2
    assert compiled["kernel_dimension"] == 2
    assert compiled["h0_certified"]

    tail = C.matrix_from_columns([2,3])
    assert tail.det() != 0
    z1 = tail.solve_right(-C.column(0))
    z2 = tail.solve_right(-C.column(1))
    k1 = vector(K,[1,0,z1[0],z1[1]])
    k2 = vector(K,[0,1,z2[0],z2[1]])
    assert C*matrix(K,[k1,k2]).transpose() == matrix(K,2,2)

    cA3_mod73 = None
    try:
        cA3_mod73 = int(red_k(cA3))
    except Exception:
        pass

    print(
        "Q80CAND1RRCAND|index={}|c={}|lambda={}|A3={}|A3_mod73={}|"
        "D4vals={}|k1={}|k2={}|status=PASS_EXACT_CANDIDATE1_RR".format(
            index,c_value,lam,cA3,cA3_mod73,d4vals,tuple(k1),tuple(k2)
        ),
        flush=True,
    )

    records.append({
        "index":int(index),
        "Hx":str(Hx),
        "Hy":str(Hy),
        "c":str(c_value),
        "lambda":str(lam),
        "D4_valuations":{str(key):list(value) for key,value in d4vals.items()},
        "A3_root":str(r66),
        "A3_chord_residue":str(cA3),
        "A3_mod73":cA3_mod73,
        "D4_root":str(r48),
        "D4_residue":"0",
        "kernel":[[str(x) for x in k1],[str(x) for x in k2]],
        "h0":2,
    })

out = ROOT / "artifacts" / "generated-results" / "q80-q4-candidate1-char0-rr-family.json"
out.write_text(json.dumps({
    "status":"PASS_EXACT_CANDIDATE1_RR_FAMILY",
    "parent_bridge_mod73":str(bridge),
    "horizontal_constraint":str(constraint),
    "horizontal_count":len(records),
    "candidates":records,
    "next":(
        "compile the finitely many exact quartics and select/certify the "
        "candidate1 child by exact fibre type; do not use the bad U=48 "
        "reduction as a section-definition condition"
    ),
},indent=2,default=int)+"\n")

print(
    "Q80CAND1PROBEFINAL|horizontal_count={}|json={}|"
    "status=PASS_EXACT_CANDIDATE1_HORIZONTAL_AND_RR_FAMILY".format(
        len(records),out
    ),
    flush=True,
)

