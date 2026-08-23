#!/usr/bin/env sage -python
"""
H3-03 q24 -> D12 modular degree-two RR preflight.

Purpose
-------
Do NOT reconstruct q24 in characteristic zero.  Use the already-certified
q24 point over GF(100003)(U) and the exact D13 parent to test the actual
Riemann--Roch architecture forced by the construction files.

Certified inputs:
  D24 = O + P24 - 7 F + V_D13,
  P24.O = 24,
  x(P24)=X/Z^2 with deg (X,Z)=(52,24),
  y(P24)=Y/Z^3 with deg Y=78,
  old-fibre degree = 2.

For the marked chord
    m = (y + y(P24))/(x - x(P24)),
the smooth P24.O collision frame gives
    a=A/Z^2,  b=B/Z,
with congruence
    A*X - B*Y == 0 mod Z^2.

Assigning the -7F twist to infinity, and using ord_inf(m)=-2, gives
    deg A <= 41,
    deg B <= 15.
Thus the predicted ambient has 42+16 = 58 columns.  The collision quotient
has dimension 48.  This script checks its exact modular rank and reports the
remaining dimension before the single I9* resolved block is imposed.

It also profiles the actual I9* local model and q24 chord jets.  Those jet
ranks are DIAGNOSTIC ONLY: they are not substituted for a resolved I9*
component cover.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
        h/"projects"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

MOD=LOCAL/f"q24-degree46-direct-global-mod-{args.prime}.json"
TRANS=LOCAL/"q8-q24-physical-to-equation-translation.json"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((p for p in q8_candidates if p.exists() and json.loads(p.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"),None)
if Q8 is None:
    raise SystemExit("No passing exact D13 q8 child artifact")

for path in (MOD,TRANS,Q8):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
trans=json.loads(TRANS.read_text())
q8=json.loads(Q8.read_text())

assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert trans["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

qe=trans["q24_equation"]
assert qe["child_root_lattice"]=="D12"
assert qe["child_root_data"]==[12,264,4]
assert qe["MW_rank_if_rho19"]==5
assert qe["height"]=="52"
assert qe["D13_local_correction"]=="0"
assert qe["P_dot_O"]==24
assert qe["vertical_fibre_coefficient"]==-7
assert len(qe["vertical_root_coefficients"])==13

p=ZZ(args.prime)
F=GF(p)
R=PolynomialRing(F,"U")
U=R.gen()
K=R.fraction_field()


def red_q(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError(f"denominator divisible by {p}")
    return F(ZZ(q.numerator()))/F(d)


def red_poly(values):
    return R([red_q(QQ(v)) for v in values])

child=q8["child"]
A=red_poly(child["minimal_A_coefficients_low_to_high"])
B=red_poly(child["minimal_B_coefficients_low_to_high"])
Delta=-16*(4*A**3+27*B**2)

sec=mod["section_mod_p"]
Z=R([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
X=R([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
Y=R([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
assert Z.degree()==24 and X.degree()==52 and Y.degree()==78
assert Z.is_monic()
assert X.gcd(Z).degree()==0 and Y.gcd(Z).degree()==0
assert Y**2 == X**3 + A*X*Z**4 + B*Z**6

xP=K(X)/K(Z**2)
yP=K(Y)/K(Z**3)
assert yP**2==xP**3+K(A)*xP+K(B)

print(
    "Q24D12PREFLIGHT_INPUT|"
    "old_degree=2|PdotO=24|height=52|corr=0|vertical_F=-7|"
    f"vertical_roots={','.join(map(str,qe['vertical_root_coefficients']))}|"
    "status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Infinity envelope from the certified fibre twist.
# ---------------------------------------------------------------------------
# On a short elliptic K3 x,y have base weights 4,6.  The marked q24 point
# has exactly the same leading weights (52-48=4, 78-72=6), hence m has
# base-infinity order -2.  A coefficient section a+b*m must vanish to order
# 7 because the chosen divisor representative has -7F at infinity.
required=7
m_inf=-2
Amax=2*Z.degree()-required
Bmax=Z.degree()+m_inf-required
assert Amax==41 and Bmax==15
ambient=[("A",i) for i in range(Amax+1)] + [("B",i) for i in range(Bmax+1)]
assert len(ambient)==58

print(
    "Q24D12PREFLIGHT_INFINITY|"
    f"required={required}|m_order={m_inf}|Amax={Amax}|Bmax={Bmax}|"
    f"ambient={len(ambient)}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Smooth P.O collision quotient: A*X-B*Y == 0 mod Z^2.
# ---------------------------------------------------------------------------
modulus=Z**2
assert modulus.degree()==48

cols=[]
for kind,i in ambient:
    if kind=="A":
        cols.append((U**i*X)%modulus)
    else:
        cols.append((-U**i*Y)%modulus)

C=matrix(F,modulus.degree(),len(cols),lambda row,col: cols[col][row])
rank=C.rank()
nullity=C.ncols()-rank
print(
    "Q24D12PREFLIGHT_COLLISION|"
    f"rows={C.nrows()}|cols={C.ncols()}|rank={rank}|nullity={nullity}|"
    f"expected_rank=48|expected_nullity=10|"
    f"status={'PASS' if rank==48 and nullity==10 else 'UNEXPECTED'}",
    flush=True,
)
if rank!=48 or nullity!=10:
    raise ArithmeticError("smooth-collision architecture differs from 58 -> 10 prediction")

K10=C.right_kernel().basis_matrix()
assert K10.dimensions()==(10,58)

# ---------------------------------------------------------------------------
# Actual I9* local fibre and q24 marked point.
# ---------------------------------------------------------------------------
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
# The stored factor is over QQ. Parse it in QQ[U] first, then reduce its
# coefficients modulo p. Parsing the string directly in GF(p)[U] makes
# Sage mix a finite-field generator with Rational Field constants.
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=R([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]

print(
    "Q24D12PREFLIGHT_I9FACTOR|"
    f"QQ_factor={fQ}|modp_factor={f}|base={int(alpha)}|status=PASS",
    flush=True,
)

S=PolynomialRing(F,"u")
u=S.gen()
KS=S.fraction_field()

def shift(poly):
    return S(poly(alpha+u))

Al=shift(A); Bl=shift(B); Dl=shift(Delta)
assert Al.valuation()==2 and Bl.valuation()==3 and Dl.valuation()==15

Zl=shift(Z); Xl=shift(X); Yl=shift(Y)
assert Zl[0] != 0  # P24 does not collide with O at the I9* place.
xloc=KS(Xl)/KS(Zl**2)
yloc=KS(Yl)/KS(Zl**3)
assert xloc.valuation()>=0 and yloc.valuation()>=0
x0=F(xloc(0)); y0=F(yloc(0))
assert y0**2==x0**3
smooth_reduction = bool(x0 or y0)
assert smooth_reduction

# P24 has trivial D13 correction, hence meets the identity component.  On the
# singular affine cubic at u=0 this appears as a smooth finite reduction in
# this chart (the zero section itself is at infinity).
node_m = -yloc/xloc
assert node_m.valuation()>=0

print(
    "Q24D12PREFLIGHT_I9STAR|"
    f"base={int(alpha)}|orders=2,3,15|Z_unit=1|"
    f"P_reduction={int(x0)},{int(y0)}|smooth={int(smooth_reduction)}|"
    f"m_node={int(F(node_m(0)))}|status=PASS_LOCAL_INPUT",
    flush=True,
)

# ---------------------------------------------------------------------------
# Diagnostic only: image of the 10-dimensional collision kernel in ordinary
# u-jets at the singular affine chart.  The true next certificate must replace
# these jets by restrictions to the RESOLVED I9* component charts.
# ---------------------------------------------------------------------------
def local_rational_from_row(row):
    AA=R.zero(); BB=R.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A": AA += row[j]*U**i
        else: BB += row[j]*U**i
    a=K(AA)/K(Z**2)
    b=K(BB)/K(Z)
    return a,b

local_pairs=[local_rational_from_row(row) for row in K10.rows()]
jet_ranks={}
for n in range(1,13):
    values=[]
    for a,b in local_pairs:
        an=KS(shift(R(a.numerator())))/KS(shift(R(a.denominator())))
        bn=KS(shift(R(b.numerator())))/KS(shift(R(b.denominator())))
        fn=an+bn*node_m
        # denominator is a unit here; extract first n Taylor coefficients.
        vals=[]
        den=S(fn.denominator()); num=S(fn.numerator())
        invden=den.inverse_mod(u**n)
        rem=(num*invden)%(u**n)
        vals=[rem[i] for i in range(n)]
        values.append(vals)
    M=matrix(F,n,10,lambda i,j:values[j][i])
    jet_ranks[n]=int(M.rank())

print(
    "Q24D12PREFLIGHT_JETS|"
    + "|".join(f"j{n}={r}" for n,r in jet_ranks.items())
    + "|scope=DIAGNOSTIC_NOT_RESOLVED|status=PASS",
    flush=True,
)

needed=8
print(
    "Q24D12PREFLIGHT_TARGET|"
    f"after_collision={nullity}|required_final=2|"
    f"required_I9star_codimension={needed}|"
    "next=RESOLVE_ACTUAL_I9STAR_AND_DERIVE_CONNECTED_CHORD_QUOTIENT|"
    "status=PASS_ACTIONABLE_GATE",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-d12-modp-rr-preflight.v1",
    "status":"PASS_H3_Q24_D12_MODP_RR_PREFLIGHT",
    "prime":int(p),
    "divisor":{
        "old_fibre_degree":2,
        "P_dot_O":24,
        "height":"52",
        "D13_local_correction":"0",
        "vertical_fibre_coefficient":-7,
        "vertical_root_coefficients":qe["vertical_root_coefficients"],
    },
    "infinity":{
        "required_order":required,
        "marked_chord_order":m_inf,
        "A_max_degree":Amax,
        "B_max_degree":Bmax,
        "ambient_dimension":58,
    },
    "smooth_collision":{
        "modulus":"Z^2",
        "degree":48,
        "rank":int(rank),
        "nullity_after_collision":int(nullity),
    },
    "I9star":{
        "base":int(alpha),
        "orders_A_B_Delta":[2,3,15],
        "Z_unit":True,
        "P_reduction":[int(x0),int(y0)],
        "marked_chord_node_residue":int(F(node_m(0))),
        "ordinary_jet_ranks_diagnostic":jet_ranks,
        "required_resolved_codimension":8,
    },
    "next":(
        "Resolve the actual I9* fibre of the canonical D13 equation, identify "
        "the transported q24 vertical class in the geometric component basis, "
        "restrict the 10-dimensional post-collision space to the required "
        "resolved charts, and certify rank 8 / final h0=2."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-d12-rr-preflight-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D12PREFLIGHT_RESULT|ambient=58|collision=48|post_collision=10|"
    "I9star_needed=8|target_h0=2|status=PASS_H3_Q24_D12_MODP_RR_PREFLIGHT",
    flush=True,
)
