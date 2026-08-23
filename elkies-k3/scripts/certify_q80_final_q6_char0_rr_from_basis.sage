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
# 2. Exact final horizontal from the recovered MW-basis difference artifact.
#
# The direct three-node section scheme is deliberately bypassed.  The exact
# section was reconstructed as a difference of easier high-incidence sections
# and identified only by reduction to the transported historical P2-P3 seed.
# ---------------------------------------------------------------------------
HART = OUT / "q80-final-q6-char0-horizontal-via-basis-sections.json"
if not HART.exists():
    raise SystemExit(f"missing recovered final-horizontal artifact: {HART}")

hdata = json.loads(HART.read_text())
assert hdata["status"] == "PASS_EXACT_FINAL_Q6_HORIZONTAL_VIA_BASIS_SECTIONS"
hsel = hdata["selected"]

Hx = R(sage_eval(hsel["Hx"], locals={"j":j, "W":W}))
Hy = R(sage_eval(hsel["Hy"], locals={"j":j, "W":W}))

assert Hx.degree() <= 4
assert Hy.degree() <= 6
assert Hy^2 == Hx^3 + A*Hx + B

def node_x(symbol):
    root = fibres[symbol][0]
    ar,br = K(A(root)),K(B(root))
    assert ar != 0
    x0 = K(-3*br/(2*ar))
    assert x0^3+ar*x0+br == 0
    assert 3*x0^2+ar == 0
    return x0

hit_symbols = ("I3","I4","I6")
for symbol in fibres:
    root = fibres[symbol][0]
    hit = bool(Hx(root) == node_x(symbol) and Hy(root) == 0)
    assert hit == (symbol in hit_symbols), (symbol, hit)

hx73 = red_poly(Hx)
hy73 = red_poly(Hy)
if hx73 != X_seed:
    raise ArithmeticError(
        "basis-recovered exact H does not reduce to transported historical Hx"
    )
if hy73 == -Y_seed:
    # Keep the downstream historical chord-slope scale aligned with the
    # selected exact sign.
    u3 = -u3
    Y_seed = -Y_seed
elif hy73 != Y_seed:
    raise ArithmeticError(
        "basis-recovered exact H does not reduce to either historical Hy sign"
    )

print(
    "Q80FINALHART|artifact={}|Hx={}|Hy={}|hits={}|"
    "status=PASS_LOAD_EXACT_FINAL_HORIZONTAL_BY_MW_BASIS_DIFFERENCE".format(
        HART,Hx,Hy,hit_symbols
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

# Only a leading local jet is needed for the connected-A5 quotient residue.
# Precision 9 is enough: the I6 smoothing term starts at s^6, while every
# supported component has local order <= 5. Newton doubles the s-adic
# accuracy, so three iterations from the correct constant root give >=8
# correct orders.
LS = LaurentSeriesRing(Kr,"s",default_prec=9)
ss = LS.gen()

def shifted(poly,root):
    ans = LS(0)
    power = LS(1)
    base = LS(Kr(root)) + ss
    for coefficient in R(poly).list():
        ans += Kr(coefficient)*power
        power *= base
    return ans

def newton_sqrt(value,root0):
    root = LS(Kr(root0))
    for unused in range(3):
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

def toric_leading_residue(component,rho_start):
    """
    Return (valuation(m), residue(m)) for m=(y+H_y)/(x-H_x) on one
    resolved I6 component.

    Crucially, never divide Laurent series: when val(num)=val(den), the
    residue is just num[v]/den[v].  This avoids the symbolic K(r)[[s]]
    quotient that dominated the previous runtime.
    """
    aa = LS(rr)*ss^component
    bb = unit*ss^(6-component)/LS(rr)
    yy = (aa+bb)/2
    ww = (bb-aa)/2

    rho_series = LS(Kr(rho_start))
    for unused in range(3):
        rho_series -= (
            rho_series^3-3*center*rho_series-ww
        )/(3*rho_series^2-3*center)

    xx = center+ww/rho_series
    numerator = yy+Hyloc
    denominator = xx-Hxloc

    vn = int(numerator.valuation())
    vd = int(denominator.valuation())
    valuation = vn-vd

    if valuation < 0:
        return valuation,None
    if valuation > 0:
        return valuation,Kr(0)

    # Exact residue in K(r); no Laurent-series division.
    residue = Kr(numerator[vn])/Kr(denominator[vd])
    return 0,residue

def reduce_a5_row(row):
    vals = tuple(red_k(v) for v in row)
    pivot = next(v for v in vals if v)
    return tuple(v/pivot for v in vals)

e6 = exact_mod_root["I6"]
target_A5_residue = total_s(e6)*Fp(69)  # historical -4, transported gauge
target_A5 = (Fp(1),e6,e6^2,target_A5_residue)

def matching_rows_from_residue(mres):
    """
    If mres=N(r)/D(r), clearing denominators in
        a0+a1*r6+a2*r6^2+a3*mres
    gives coefficient rows
        (D_k, r6*D_k, r6^2*D_k, N_k).
    We need only rows whose canonical reduction is the pinned target, so
    inspect N_k/D_k directly instead of constructing/lcm'ing every row.
    """
    mres = Kr(mres)
    if not mres:
        return ()

    num = mres.numerator()
    den = mres.denominator()
    answer = []

    for degree in range(int(den.degree())+1):
        dk = K(den[degree])
        if not dk:
            continue
        nk = K(num[degree]) if degree <= num.degree() else K(0)
        c_exact = K(nk/dk)
        if red_k(c_exact) != target_A5_residue:
            continue

        row = (K(1),K(r6),K(r6^2),c_exact)
        assert reduce_a5_row(row) == target_A5
        answer.append((degree,row))

    return tuple(answer)

supports = {
    "direct":(1,3,4),
    "reversed":(2,3,5),
}

selected_A5 = None

# Stop as soon as the exact coefficient row reducing to the pinned quotient
# is found.  The old code evaluated all 12 component/sign/orientation cases,
# expanded all K(r) numerator/denominator coefficients, and deduplicated them.
for rho_candidate in (rho0,-rho0):
    if selected_A5 is not None:
        break
    for orientation,support_candidate in supports.items():
        if selected_A5 is not None:
            break

        for component in support_candidate:
            valuation,mres = toric_leading_residue(component,rho_candidate)

            if valuation < 0:
                print(
                    "Q80FINALA5FAST|rho_mod73={}|orientation={}|component={}|"
                    "mvaluation={}|status=SKIP_POLE".format(
                        int(red_k(rho_candidate)),orientation,component,valuation
                    ),
                    flush=True,
                )
                continue

            matches_here = matching_rows_from_residue(mres)
            print(
                "Q80FINALA5FAST|rho_mod73={}|orientation={}|component={}|"
                "mvaluation={}|mres_num_degree={}|mres_den_degree={}|"
                "matching_coefficients={}|status=PASS_FAST_I6_COMPONENT".format(
                    int(red_k(rho_candidate)),orientation,component,valuation,
                    -1 if mres is None or not mres else int(mres.numerator().degree()),
                    -1 if mres is None or not mres else int(mres.denominator().degree()),
                    tuple(degree for degree,unused in matches_here),
                ),
                flush=True,
            )

            if matches_here:
                coefficient_degree,a5_row = matches_here[0]
                selected_A5 = (
                    rho_candidate,orientation,support_candidate,
                    component,coefficient_degree,a5_row
                )
                break

if selected_A5 is None:
    raise ArithmeticError(
        "fast exact I6 leading-jet traces did not contain the transported "
        "pinned A5 quotient line"
    )

rho,orientation,support,matched_component,coefficient_degree,a5_row = selected_A5

print(
    "Q80FINALA5|root={}|root_mod73={}|orientation={}|support={}|"
    "matched_component={}|coefficient_degree={}|row_mod73={}|"
    "status=PASS_EXACT_CONNECTED_A5_QUOTIENT".format(
        r6,int(e6),orientation,support,matched_component,coefficient_degree,
        tuple(int(v) for v in reduce_a5_row(a5_row)),
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
