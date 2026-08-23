#!/usr/bin/env sage
"""
Exact characteristic-zero q4_1938 resolved-RR compiler.

Historical construction replay:
  raw 5D chord ambient:
      f = A(T)/q^2 + B*m/q, deg A <= 3, B constant,
      m = chord through -H;
  smooth H.O=1 saturation:
      A*N_H - B*M_H == 0 mod q^2                 (rank 2);
  vertical A4=(-1,-1,-1,0):
      one connected quotient line on the selected I5 fibre (rank 1).

The resulting 3x5 exact matrix must have nullity 2 and reduce to the pinned
GF(73) q4_1938 kernel.  The same exact pencil is then compiled by binary
quartic invariants and checked against the pinned q4_1938 child.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, GF, PolynomialRing, FunctionField, LaurentSeriesRing,
    matrix, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
QDIR = DATA / "q80-q6-7774-char0"
ART = ROOT / "artifacts" / "generated-results" / "q80-q6-7774-char0-resolved-rr.json"
HFILE = QDIR / "q80_char0_q6_7774_H1938.sage"
OUT = ROOT / "artifacts" / "generated-results"
OUT.mkdir(parents=True, exist_ok=True)

if not ART.exists():
    raise SystemExit(f"missing certified 7774 artifact: {ART}")
if not HFILE.exists():
    raise SystemExit(f"missing exact q4_1938 horizontal: {HFILE}")

data = json.loads(ART.read_text())
assert data["status"] == "PASS_EXACT_Q6_7774_RESOLVED_RR"

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))
load(str(HFILE))

# HFILE defines K,j,R,T,A,B,q,N_H,M_H,H_x,H_y.
Delta = -16*(4*A^3+27*B^2)
assert (A.degree(), B.degree(), Delta.degree()) == (8,12,18)

Fp = GF(73)
JMOD = Fp(17)
assert JMOD^2 == Fp(-3)

def red_q(c):
    c = QQ(c)
    return Fp(c.numerator())/Fp(c.denominator())

def red_k(c):
    c = K(c)
    cc = list(c)+[QQ(0),QQ(0)]
    return red_q(cc[0]) + JMOD*red_q(cc[1])

Rp = PolynomialRing(Fp, "T")
tp = Rp.gen()

def red_poly(f):
    f = R(f)
    return Rp([red_k(c) for c in f.list()])

scale_marker = int(data["child"]["gf73_scale_marker"])
if scale_marker <= 0:
    raise ArithmeticError(f"unexpected inverse 7774 scale marker {scale_marker}")
u73 = Fp(scale_marker)

# Exact selected I5 roots from the already-certified 7774 fibre factorization.
i5 = {}
for item in data["child"]["finite_fibres"]:
    if item["kodaira"] != "I5" or int(item["degree"]) != 1:
        continue
    fac = R(sage_eval(item["factor"], locals={"j":j, "T":T}))
    root = K(-fac[0]/fac[1])
    i5[int(red_k(root))] = root
assert set(i5) == {14,42}
r14, r42 = i5[14], i5[42]

# ---------------------------------------------------------------------------
# 1. Exact smooth P.O=1 saturated module.
#
# Repo formula for a=A/q^2, b=B/q in the saturated frame
# <1,(m-y_P/x_P)/q> is
#
#   A*(den(p)/q) + B*num(p) == 0 mod q^2.
#
# Here the marked point is -H, hence p=-H_y/H_x and this becomes
#
#   A*N_H - B*M_H == 0 mod q^2.
# ---------------------------------------------------------------------------
pole = K(-q[0]/q[1])
assert red_k(pole) == 47
assert q.degree() == 1 and q.leading_coefficient() == 1
assert N_H(pole) != 0 and M_H(pole) != 0

ambient = ("A0","A1","A2","A3","B")
smooth_exprs = tuple(R(T^i*N_H) for i in range(4)) + (R(-M_H),)

smooth_value = tuple(K(f(pole)) for f in smooth_exprs)
smooth_deriv = tuple(K(f.derivative()(pole)) for f in smooth_exprs)
smooth_matrix = matrix(K, [smooth_value, smooth_deriv])
assert smooth_matrix.rank() == 2

print(
    "Q801938RRSMOOTH|pole={}|pole_mod73=47|rank=2|"
    "status=PASS_EXACT_SMOOTH_PO1_SATURATION".format(pole),
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Exact resolved I5 / connected-A4 trace.
# ---------------------------------------------------------------------------
def node_x(root):
    ar, br = K(A(root)), K(B(root))
    x0 = K(-3*br/(2*ar))
    assert x0^3 + ar*x0 + br == 0
    assert 3*x0^2 + ar == 0
    return x0

x014 = node_x(r14)
assert K(H_x(r14)) == x014
assert K(H_y(r14)) == 0

Kr = FunctionField(K, "r")
rr = Kr.gen()
PREC = 13
LS = LaurentSeriesRing(Kr, "s", default_prec=PREC)
ss = LS.gen()

def shifted(poly, root):
    ans = LS(0)
    for degree, coefficient in enumerate(R(poly).list()):
        ans += Kr(coefficient)*(Kr(root)+ss)^degree
    return ans

def newton_sqrt(value, root0):
    root = LS(Kr(root0))
    for _ in range(5):
        root = (root+value/root)/2
    if (root^2-value).valuation() < PREC-4:
        raise ArithmeticError("I5 center square-root precision failure")
    return root

Aloc = shifted(A,r14)
Bloc = shifted(B,r14)
center = newton_sqrt(-Aloc/3, x014)
g0 = center^3+Aloc*center+Bloc
assert int(g0.valuation()) == 5
unit = g0/ss^5

rho_sq = K(3*x014)
if not rho_sq.is_square():
    raise ArithmeticError("selected I5 is not split over QQ(sqrt(-3))")
rho0 = K(rho_sq.sqrt())
rho_choices = (rho0,-rho0)

qloc = shifted(q,r14)
Hxloc = shifted(N_H,r14)/qloc^2
Hyloc = shifted(M_H,r14)/qloc^3

def toric_point(component, rho_start):
    aa = LS(rr)*ss^component
    bb = unit*ss^(5-component)/LS(rr)
    yy = (aa+bb)/2
    ww = (bb-aa)/2

    rho = LS(Kr(rho_start))
    for _ in range(5):
        residual = rho^3-3*center*rho-ww
        derivative = 3*rho^2-3*center
        rho -= residual/derivative
    if (rho^3-3*center*rho-ww).valuation() < PREC-5:
        raise ArithmeticError("I5 toric rho precision failure")
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
    for d in range(degree+1):
        row = tuple(
            K(f[d]) if f and d <= f.degree() else K(0)
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
    vals = tuple(red_k(v) for v in canonical(row))
    pivot = next(v for v in vals if v)
    return tuple(v/pivot for v in vals)

# Pinned row is in the unscaled GF73 model.  Our exact 7774 model has
# x=u^2*x_pin, y=u^3*y_pin, hence m_exact=u*m_pin.
target73 = (
    Fp(1),Fp(14),Fp(50),Fp(43),u73*Fp(46)
)

supports = {
    "direct": (1,2,3),
    "reversed": (2,3,4),
}
trials = []
selected = []
for rho in rho_choices:
    for orientation,support in supports.items():
        rows = []
        component_dump = []
        for component in support:
            xx,yy = toric_point(component,rho)
            # Historical sign=+1 uses the chord through -H.
            mloc = (yy+Hyloc)/(xx-Hxloc)
            mv = int(mloc.valuation())
            if mv < 0:
                component_dump.append((component,mv,"POLE"))
                continue
            mres = Kr(mloc[0]) if mv == 0 else Kr(0)
            vals = (
                Kr(1),Kr(r14),Kr(r14^2),Kr(r14^3),
                Kr(q(r14))*mres,
            )
            local_rows = functional_rows(vals)
            rows.extend(local_rows)
            component_dump.append(
                (component,mv,str(mres),
                 tuple(tuple(int(red_k(v)) for v in reduce_row(row))
                       for row in local_rows))
            )

        rows = dedup(rows)
        rank = matrix(K,rows).rank() if rows else 0
        redrows = tuple(reduce_row(row) for row in rows)
        trials.append((rho,orientation,support,rank,redrows,component_dump))

        print(
            "Q801938RRA4LOCAL|rho_mod73={}|orientation={}|support={}|rank={}|"
            "rows_mod73={}|status=PASS_EXACT_I5_TRACE".format(
                int(red_k(rho)),orientation,support,rank,
                tuple(tuple(int(x) for x in row) for row in redrows)
            ),
            flush=True,
        )

        for row in rows:
            if reduce_row(row) == target73:
                selected.append((rho,orientation,support,canonical(row),component_dump))

if not selected:
    raise ArithmeticError(
        "resolved I5 traces did not contain the pinned connected-A4 quotient line; "
        "inspect Q801938RRA4LOCAL diagnostics"
    )

# Deduplicate exact quotient rows, preferring direct orientation.
unique = {}
for item in selected:
    unique.setdefault(tuple(item[3]),item)
selected = list(unique.values())
selected.sort(key=lambda z: (z[1] != "direct", int(red_k(z[0]))))
rho,orientation,support,a4_row,a4_diag = selected[0]

assert reduce_row(a4_row) == target73
print(
    "Q801938RRA4|root={}|root_mod73=14|orientation={}|support={}|"
    "row_mod73={}|status=PASS_EXACT_CONNECTED_A4_QUOTIENT".format(
        r14,orientation,support,
        tuple(int(x) for x in reduce_row(a4_row))
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Complete exact RR matrix and historically normalized kernel.
# ---------------------------------------------------------------------------
smooth_block = {
    "name":"q4_1938 smooth P.O=1 saturation",
    "matrix":smooth_matrix,
    "quotient_basis":("value","first_jet"),
    "provenance":(
        "saturated frame <1,(m-y_P/x_P)/q> for P=-H; "
        "A*N_H-B*M_H == 0 mod q^2"
    ),
}
a4_block = {
    "name":"q4_1938 connected A4 quotient",
    "matrix":matrix(K,[a4_row]),
    "quotient_basis":("connected_A4_line",),
    "provenance":(
        "exact toric I5 resolution; selected connected quotient line is the "
        "exact resolved trace reducing to the pinned modular row"
    ),
}

compiled = compile_resolved_conditions(
    ambient,(smooth_block,a4_block),complete=True,coefficient_field=K
)
assert compiled["ambient_dimension"] == 5
assert compiled["rank"] == 3
assert compiled["kernel_dimension"] == 2
assert compiled["h0_certified"]

C = compiled["condition_matrix"]
tail = C.matrix_from_columns([2,3,4])
assert tail.det() != 0
z1 = tail.solve_right(-C.column(0))
z2 = tail.solve_right(-C.column(1))
k1 = vector(K,[1,0,z1[0],z1[1],z1[2]])
k2 = vector(K,[0,1,z2[0],z2[1],z2[2]])
kernel = matrix(K,[k1,k2])
assert C*kernel.transpose() == matrix(K,C.nrows(),2)
assert kernel.rank() == 2

def reduce_vector(v):
    return tuple(red_k(x) for x in v)

expected_k1 = (Fp(1),Fp(0),Fp(51),Fp(40),Fp(23)/u73)
expected_k2 = (Fp(0),Fp(1),Fp(36),Fp(27),Fp(48)/u73)
assert reduce_vector(k1) == expected_k1
assert reduce_vector(k2) == expected_k2

print(
    "Q801938RR|ambient=5|smooth_rank=2|A4_rank=1|rank=3|nullity=2|h0=2|"
    "status=PASS_EXACT_Q4_1938_RESOLVED_RR",
    flush=True,
)
print(
    "Q801938RRKERNEL|k1={}|k2={}|status=PASS_PINNED_KERNEL_LIFT".format(
        tuple(k1),tuple(k2)
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Compile the child from THIS exact resolved pencil.
# ---------------------------------------------------------------------------
Sring = PolynomialRing(K,"S")
S0 = Sring.gen()
KS = Sring.fraction_field()
RT = PolynomialRing(KS,"T")
tt = RT.gen()
S = KS(S0)

def lift_poly(f):
    return RT([KS(c) for c in R(f).list()])

qL = lift_poly(q)
NL = lift_poly(N_H)
ML = lift_poly(M_H)
AL = lift_poly(A)
BL = lift_poly(B)

def chord_coefficients(row):
    Apoly = R(sum(row[i]*T^i for i in range(4)))
    return (
        RT(lift_poly(Apoly))/qL^2,
        KS(row[4])/qL,
    )

(a0,b0) = chord_coefficients(k1)
(a1,b1) = chord_coefficients(k2)

hop = compile_degree_two_chord_hop(
    RT,S,a0,b0,a1,b1,
    NL/qL^2,-ML/qL^3,
    AL,BL,
)
quartic = hop["binary_quartic"]
assert quartic.degree() == 4

Aj = KS(hop["jacobian_a"])
Bj = KS(hop["jacobian_b"])
classification = classify_finite_short_weierstrass_fibres(Sring,Aj,Bj)

finite_totals = {}
finite_dump = []
for item in classification["finite_fibres"]:
    symbol = str(item["kodaira"])
    finite_totals[symbol] = finite_totals.get(symbol,0)+int(item["degree"])
    finite_dump.append({
        "factor":str(item["factor"]),
        "degree":int(item["degree"]),
        "kodaira":symbol,
        "orders":[int(x) for x in item["minimal_orders"]],
    })
infty = classification["infinity_boundary"]
infty_orders = tuple(int(x) for x in infty["normalized_orders"])

assert finite_totals == {"I5":2,"I4":2,"I2":1,"I1":4}
assert infty_orders[2] == 0

# ---------------------------------------------------------------------------
# 5. Exact reduction to pinned GF73 q4_1938 child.
# ---------------------------------------------------------------------------
FSp = PolynomialRing(Fp,"S")
sp = FSp.gen()
FS = FSp.fraction_field()

def red_s_poly(poly):
    poly = Sring(poly)
    return FSp([red_k(c) for c in poly.list()])

def red_s(frac):
    frac = KS(frac)
    num = red_s_poly(frac.numerator())
    den = red_s_poly(frac.denominator())
    if den == 0:
        raise ZeroDivisionError("q4_1938 child denominator vanished modulo 73")
    return FS(num)/FS(den)

Aj73,Bj73 = red_s(Aj),red_s(Bj)

Apin_num = (
    34*sp^8+32*sp^7+15*sp^6+58*sp^5+17*sp^4+52*sp^3
    +48*sp^2+22*sp+38
)
Apin_den = (
    sp^8+69*sp^7+7*sp^6+66*sp^5+50*sp^4+53*sp^3
    +5*sp^2+41*sp+2
)
Bpin_num = (
    11*sp^12+6*sp^11+5*sp^10+3*sp^9+64*sp^8+53*sp^7
    +23*sp^6+42*sp^5+32*sp^4+37*sp^3+27*sp^2+43*sp+38
)
Bpin_den = (
    sp^12+67*sp^11+53*sp^10+9*sp^9+72*sp^8+30*sp^7
    +19*sp^6+44*sp^5+41*sp^4+72*sp^3+33*sp^2+70*sp+64
)
Apin = FS(Apin_num)/FS(Apin_den)
Bpin = FS(Bpin_num)/FS(Bpin_den)

child_scale = None
for v in Fp:
    if not v:
        continue
    if Aj73 == v^4*Apin and Bj73 == v^6*Bpin:
        child_scale = int(v)
        break
    if v^4*Aj73 == Apin and v^6*Bj73 == Bpin:
        child_scale = -int(v)
        break

if child_scale is None:
    raise ArithmeticError(
        "exact q4_1938 child has correct fibres but did not reduce to pinned "
        "GF73 equation in the pinned S coordinate"
    )

print(
    "Q801938RRCHILD|quartic_degree=4|finite={}|infinity=smooth|"
    "gf73_scale_marker={}|status=PASS_EXACT_Q4_1938_CHILD".format(
        finite_totals,child_scale
    ),
    flush=True,
)

payload = {
    "status":"PASS_EXACT_Q4_1938_RESOLVED_RR",
    "field":"QQ(sqrt(-3))",
    "horizontal_file":str(HFILE.relative_to(ROOT)),
    "smooth_collision":{
        "pole":str(pole),
        "rank":2,
        "congruence":"A*N_H-B*M_H == 0 mod q^2",
        "marked_point":"-H",
    },
    "connected_A4":{
        "I5_root":str(r14),
        "I5_mod73":14,
        "orientation":orientation,
        "support":[int(x) for x in support],
        "rho":str(rho),
        "row":[str(x) for x in a4_row],
        "row_mod73":[int(x) for x in reduce_row(a4_row)],
        "diagnostics":a4_diag,
    },
    "rr":{
        "ambient_dimension":5,
        "smooth_rank":2,
        "A4_rank":1,
        "condition_rank":3,
        "kernel_dimension":2,
        "h0":2,
        "kernel":[[str(x) for x in k1],[str(x) for x in k2]],
    },
    "child":{
        "quartic":str(quartic),
        "jacobian_A":str(Aj),
        "jacobian_B":str(Bj),
        "finite_fibres":finite_dump,
        "infinity_orders":[int(x) for x in infty_orders],
        "global_pattern":"2 I5 + 2 I4 + I2 + 4 I1",
        "gf73_scale_marker":int(child_scale),
    },
    "next":"q4_6855 characteristic-zero propagation",
}

json_path = OUT / "q80-q4-1938-char0-resolved-rr.json"
json_path.write_text(json.dumps(payload,indent=2,default=int)+"\n")

note_path = QDIR / "Q80_CHAR0_Q4_1938_RESOLVED_RR_CERTIFICATE.md"
note_path.write_text(
    "# Q80 q4_1938 — exact characteristic-zero resolved RR certificate\n\n"
    "Status: **PASS_EXACT_Q4_1938_RESOLVED_RR**\n\n"
    "The historical modular construction was replayed exactly over "
    "`QQ(sqrt(-3))`:\n\n"
    "- raw 5D marked-chord ambient;\n"
    "- smooth `P.O=1` saturated collision block of rank 2 using the chord "
    "through `-H`;\n"
    "- one connected A4 quotient on the exact I5 fibre reducing to `T=14`;\n"
    "- total rank 3 and exact `h0=2`;\n"
    "- binary-quartic Jacobian from that same certified pencil;\n"
    "- fibres `2 I5 + 2 I4 + I2 + 4 I1`;\n"
    f"- pinned GF(73) regression, constant scale marker `{child_scale}`.\n\n"
    "Next: q4_6855 characteristic-zero propagation.\n"
)

print(
    "Q801938RRFNAL|json={}|certificate={}|"
    "status=PASS_EXACT_Q4_1938_RESOLVED_RR".format(json_path,note_path),
    flush=True,
)
