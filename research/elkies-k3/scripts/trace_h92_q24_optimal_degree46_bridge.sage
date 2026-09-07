#!/usr/bin/env sage -python
"""
Construct and directly trace the optimal equation-level H92 q24 bridge.

Equation search result:
    q6 old-zero MW word  W = (-2,-1,1)
    q8 degree                46
    AJ(W)                    (0,1,-1,1)
    q24                      AJ(W) + 2*G1
up to simultaneous global inversion.

Since standard_zero has old-group coordinate (-2,1,0),
the standard q6 MW coordinate of W is
    (0,-2,1) = Qmap + S3.

Thus W is already explicit on the exact q6 Weierstrass equation.

This script:
  1. constructs W = Qmap + S3 exactly over QQ(T);
  2. replays the repaired q8 RR pencil and proves deg(U|W)=46;
  3. modulo p at U=tau, maps all 46 points directly through the anchored
     II*_E8_1 quartic -> canonical D13 birational map;
  4. sums them with L(47 O), no covariant 2-cover and no halving;
  5. forms both orientation-compatible q24 candidates
         AJ_direct + 2 G1
         AJ_direct - 2 G1
     because the global quartic square-root sign may negate AJ_direct;
  6. when the earlier independent q24 modular artifact exists, selects the
     correct candidate by matching its x-coordinate.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, sage_eval
)


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
        try: c=c.resolve()
        except Exception: continue
        if c in seen: continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate repo")


def qpoly(R, values):
    return R([QQ(v) for v in values])


def qrat(K,R,data,nk,dk):
    return K(qpoly(R,data[nk]))/K(qpoly(R,data[dk]))


def rf_payload(v,R):
    return {
        "numerator_coefficients_low_to_high":[str(x) for x in R(v.numerator()).list()],
        "denominator_coefficients_low_to_high":[str(x) for x in R(v.denominator()).list()],
    }


def monic_power_root(value, exponent):
    out=value.parent().one()
    for f,m in value.factor():
        assert m % exponent == 0
        out *= f.monic()**(m//exponent)
    return out.monic()


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--tau",type=int,default=2)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"
CORE=ROOT/"elkies-k3/scripts/elliptic_neighbor_compiler.sage"

Q6=GEN/"elkies-k3-h92-q6-child-jacobian.json"
ZERO=GEN/"elkies-k3-h92-q6-child-zero-section.json"
COMP=GEN/"elkies-k3-h92-q6-child-e7-infinity-sections.json"
S3BR=LOCAL/"q6-third-to-q8-bridge.json"
ANCHOR=LOCAL/"q8-d13-branch-anchor.json"
EQNS=LOCAL/"q8-equation-ns-divisor.json"
BACK=LOCAL/"q8-q24-canonical-backtrack.json"
G3ART=LOCAL/"q8-d13-g3-from-e77-bisection.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    p for p in q8_candidates
    if p.exists()
    and "rr" in json.loads(p.read_text())
    and "kernel_polynomials" in json.loads(p.read_text()).get("rr",{})
    and "child" in json.loads(p.read_text())
),None)
if Q8 is None:
    raise SystemExit("No complete corrected q8 child artifact")

for pth in (CORE,Q6,ZERO,COMP,S3BR,ANCHOR,EQNS,Q8,G3ART):
    if not pth.exists():
        raise SystemExit(f"missing {pth}")

OUT=args.output.resolve() if args.output else LOCAL/f"q24-degree46-direct-trace-mod-{args.prime}-tau-{args.tau}.json"

scope={}
exec(compile(CORE.read_text(),str(CORE),"exec"),scope)
squarefree_binary_quartic=scope["squarefree_binary_quartic"]

q6=json.loads(Q6.read_text())
zero=json.loads(ZERO.read_text())
comp=json.loads(COMP.read_text())
s3br=json.loads(S3BR.read_text())
anchor=json.loads(ANCHOR.read_text())
eqns=json.loads(EQNS.read_text())
q8=json.loads(Q8.read_text())
g3art=json.loads(G3ART.read_text())

assert q6["status"]=="PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"]=="PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert comp["status"]=="PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert s3br["status"]=="PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert anchor["status"]=="PASS_EXACT_D13_BRANCH_ANCHOR"
assert eqns["status"]=="PASS_EXACT_Q8_EQUATION_NS_DIVISOR"
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert g3art["status"]=="PASS_EXACT_D13_G3_FROM_E77_BISECTION"

# ===========================================================================
# 1. Exact QQ(T) bridge W = Qmap + S3.
# ===========================================================================

RQ=PolynomialRing(QQ,"T"); TQ=RQ.gen(); KQ=RQ.fraction_field()
model=q6["minimal_short_weierstrass"]
A6q=qpoly(RQ,model["A_coefficients_low_to_high"])
B6q=qpoly(RQ,model["B_coefficients_low_to_high"])
EQ6=EllipticCurve(KQ,[0,0,0,KQ(A6q),KQ(B6q)])

zd=zero["section"]
Pold=EQ6(
    qrat(KQ,RQ,zd,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
    qrat(KQ,RQ,zd,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
)

entries={e["sign"]:e for e in comp["sections"]}
points={
    sign:EQ6(
        qrat(KQ,RQ,e,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
        qrat(KQ,RQ,e,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
    )
    for sign,e in entries.items()
}
affine=points[comp["source"]["affine_E7_sign"]]
e77=points[comp["source"]["E7_7_sign"]]

Pmap=e77-Pold
Qmap=e77-affine

s3d=s3br["third_section_canonical_q6"]
S3=EQ6(
    qrat(KQ,RQ,s3d["x"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
    qrat(KQ,RQ,s3d["y"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
)

W=Qmap-S3
assert W in EQ6 and not W.is_zero()
wxq,wyq=W.xy()
assert wyq**2==wxq**3+KQ(A6q)*wxq+KQ(B6q)

# Repaired q8 marked section and RR pencil.
md=q8["marking"]["section"]
sxq=qrat(KQ,RQ,md,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high")
syq=qrat(KQ,RQ,md,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high")
Smark=EQ6(sxq,syq)
assert Smark==Pmap+Qmap

nxq,dxq=RQ(sxq.numerator()),RQ(sxq.denominator())
nyq,dyq=RQ(syq.numerator()),RQ(syq.denominator())
hq=monic_power_root(dxq,2)
assert hq==monic_power_root(dyq,3) and hq.degree()==10

iiq=RQ(next(x for x in q6["finite_fibres"] if x["kodaira"]=="II*")["factor"]).monic()
ivq=RQ(next(x for x in q6["finite_fibres"] if x["kodaira"]=="IV*")["factor"]).monic()
Mq=(iiq**2*ivq**2).monic()
normalizerq=(nyq*dxq*(hq*dyq).inverse_mod(nxq)).mod(nxq)
pfunq=-syq/sxq
rhoq=(normalizerq*nxq.inverse_mod(Mq)).mod(Mq)

pairsq=[]
for entry in q8["rr"]["kernel_polynomials"]:
    sp=RQ(entry["s"]); tp=RQ(entry["t"])
    Bc=KQ(sp)/KQ(hq)
    Ac=(
        -KQ(sp)*pfunq/KQ(hq)
        -KQ(sp)*KQ(normalizerq)/KQ(nxq)
        +KQ(sp*rhoq)+KQ(tp*Mq)
    )
    pairsq.append((Ac,Bc))
(A0q,B0q),(A1q,B1q)=pairsq

mWq=(wyq+syq)/(wxq-sxq)
UWq=KQ((A1q+B1q*mWq)/(A0q+B0q*mWq))
Unq,Udq=RQ(UWq.numerator()),RQ(UWq.denominator())
assert Unq.gcd(Udq) in QQ
degree=max(Unq.degree(),Udq.degree())
assert degree==46,degree

print(
    "Q24D46_BRIDGE|"
    "old_mw=-2,-1,-1|standard_mw=0,-2,-1|formula=Qmap-S3|"
    f"x={RQ(wxq.numerator()).degree()}/{RQ(wxq.denominator()).degree()}|"
    f"y={RQ(wyq.numerator()).degree()}/{RQ(wyq.denominator()).degree()}|"
    f"q8={Unq.degree()}/{Udq.degree()}|degree={degree}|"
    "status=PASS_EXACT_DEGREE46",
    flush=True,
)

# ===========================================================================
# 2. Reduce to GF(p), specialize U=tau, and directly trace 46 points.
# ===========================================================================

p=ZZ(args.prime)
if not p.is_prime() or p in (2,3):
    raise ValueError("prime must be odd and != 3")
F=GF(p)
R=PolynomialRing(F,"T"); T=R.gen(); K=R.fraction_field()

def modq(x):
    x=QQ(x); d=ZZ(x.denominator())
    if d%p==0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {x}")
    return F(ZZ(x.numerator()))/F(d)

def red_poly(poly):
    poly=RQ(poly)
    return R([modq(v) for v in poly.list()])

def red_rf(v):
    v=KQ(v)
    return K(red_poly(v.numerator()))/K(red_poly(v.denominator()))

A6=red_poly(A6q); B6=red_poly(B6q)
sx,sy=red_rf(sxq),red_rf(syq)
wx,wy=red_rf(wxq),red_rf(wyq)
A0,B0=red_rf(A0q),red_rf(B0q)
A1,B1=red_rf(A1q),red_rf(B1q)
Un,Ud=red_poly(Unq),red_poly(Udq)

tau=F(args.tau)

A13p=[modq(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]]
B13p=[modq(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]]
def eval_coeff(vals,t):
    return sum(v*t**i for i,v in enumerate(vals))

A13=eval_coeff(A13p,tau)
B13=eval_coeff(B13p,tau)
E13=EllipticCurve(F,[0,0,0,A13,B13])
if not E13.discriminant():
    raise ArithmeticError("singular D13 specialization")

H=R(Un-tau*Ud)
if H.degree()!=46:
    raise ArithmeticError(f"degree-46 fibre dropped to {H.degree()}")
H=H.monic()
if H.gcd(H.derivative()).degree()!=0:
    raise ArithmeticError("degree-46 bridge fibre is not etale")

print(
    f"Q24D46_TRACE|prime={p}|tau={int(tau)}|degree=46|"
    "stage=setup|status=PASS",
    flush=True,
)

q8_m=-(A1-tau*A0)/(B1-tau*B0)
radicand=q8_m**4-6*sx*q8_m**2-8*sy*q8_m-3*sx**2-4*K(A6)
quartic,square_factor=squarefree_binary_quartic(radicand,R)
if quartic.degree()!=4:
    raise ArithmeticError("q8 quartic degree dropped")

def mod_H(v):
    v=K(v)
    num,den=R(v.numerator()),R(v.denominator())
    if den.gcd(H).degree()!=0:
        raise ZeroDivisionError("denominator not invertible modulo bridge fibre")
    return (num*den.inverse_mod(H))%H

mW=(wy+sy)/(wx-sx)
if mod_H(q8_m-mW):
    raise ArithmeticError("q8 chord mismatch on degree-46 bridge")
wW=(2*wx+sx-q8_m**2)/square_factor
wA=mod_H(wW)
if (wA*wA-quartic)%H:
    raise ArithmeticError("bridge quartic square-root mismatch")

# Exact branch anchor specialized at U=tau.
RUQ=PolynomialRing(QQ,"U"); UQ=RUQ.gen(); KUQ=RUQ.fraction_field()
def parse_u(text):
    return KUQ(sage_eval(str(text),locals={"U":UQ}))
def spec_u(text):
    v=parse_u(text)
    num,den=RUQ(v.numerator()),RUQ(v.denominator())
    n=sum(modq(c)*tau**i for i,c in enumerate(num.list()))
    d=sum(modq(c)*tau**i for i,c in enumerate(den.list()))
    if not d:
        raise ZeroDivisionError("anchor denominator vanished")
    return n/d

tii=modq(QQ(anchor["zero"]["old_base_T"]))
coef=anchor["quartic_to_anchor"]["shifted_coefficients"]
aa=spec_u(coef["a_r4"]); bb=spec_u(coef["b_r3"])
cc=spec_u(coef["c_r2"]); dd=spec_u(coef["d_r1"])

rpoly=T-tii
branch_poly=dd*rpoly+cc*rpoly**2+bb*rpoly**3+aa*rpoly**4
scale=quartic[4]/branch_poly[4]
if quartic != scale*branch_poly:
    raise ArithmeticError("direct and anchored quartics are not scalar-equivalent")
if not scale.is_square():
    raise ArithmeticError("quartic scale nonsquare")
wbA=wA/scale.sqrt()
if (wbA*wbA-branch_poly)%H:
    raise ArithmeticError("branch W conversion failed")

rA=rpoly%H
if rA.gcd(H).degree()!=0:
    raise ZeroDivisionError("bridge meets branch zero")
rinv=rA.inverse_mod(H)
Xa=(dd*rinv)%H
Ya=(dd*wbA*rinv**2)%H

urst=anchor["anchor_to_canonical"]["urst"]
u=spec_u(urst[0]); rr=spec_u(urst[1]); ss=spec_u(urst[2]); tt=spec_u(urst[3])
if not u:
    raise ZeroDivisionError("anchor u vanished")

xA=((Xa-rr)/(u**2))%H
yA=((Ya-ss*(Xa-rr)-tt)/(u**3))%H
if (yA*yA-xA*xA*xA-A13*xA-B13)%H:
    Ya=(-dd*wbA*rinv**2)%H
    yA=((Ya-ss*(Xa-rr)-tt)/(u**3))%H
    if (yA*yA-xA*xA*xA-A13*xA-B13)%H:
        raise ArithmeticError("direct bridge images miss D13 child")
    w_sign=-1
else:
    w_sign=1

# L(47 O): x^0..x^23 and y*x^0..x^22.
one=R.one()
xp=[one]
for unused in range(23):
    xp.append((xp[-1]*xA)%H)
assert len(xp)==24
columns=list(xp)+[(yA*xp[e])%H for e in range(23)]
assert len(columns)==47

Eval=matrix(F,46,47,lambda row,col:columns[col][row])
ker=Eval.right_kernel().basis_matrix()
if ker.nrows()!=1:
    raise ArithmeticError(f"L(47O) trace kernel dimension {ker.nrows()}")
rel=ker[0]

XR=PolynomialRing(F,"X"); X=XR.gen()
Afun=sum(rel[i]*X**i for i in range(24))
Bfun=sum(rel[24+i]*X**i for i in range(23))
Rint=Afun**2-(X**3+A13*X+B13)*Bfun**2
if Rint.degree()!=47:
    raise ArithmeticError(f"residual intersection degree {Rint.degree()}, expected 47")
root_sum=-Rint[46]/Rint[47]

def newton_power_sums(poly):
    n=poly.degree()
    assert poly[n]==1
    sums=[F(n)]
    for k in range(1,n):
        total=F(k)*poly[n-k]
        for j in range(1,k):
            total += poly[n-j]*sums[k-j]
        sums.append(-total)
    return sums

ps=newton_power_sums(H)
trace_x=sum(xA[i]*ps[i] for i in range(46))
xQ=root_sum-trace_x
bQ=Bfun(xQ)
if not bQ:
    raise ArithmeticError("residual B(x_Q)=0")
yQ=-Afun(xQ)/bQ
AJ=-E13(xQ,yQ)

if AJ.is_zero():
    raise ArithmeticError("AJ(W) unexpectedly zero")
ajx,ajy=AJ.xy()

print(
    "Q24D46_AJ|"
    f"AJ_x={int(ajx)}|AJ_y={int(ajy)}|W_sign={w_sign}|"
    "method=direct_L47|status=PASS_DIRECT_AJ",
    flush=True,
)

# ===========================================================================
# 3. Solve the q24 correction directly on the canonical D13 elliptic curve.
# ===========================================================================

if not BACK.exists():
    raise SystemExit(
        "independent q8-q24-canonical-backtrack.json is required for the "
        "equation-level correction solve"
    )

back=json.loads(BACK.read_text())
if back.get("status")!="PASS_EXPLICIT_Q24_MODP_FROM_AJ_G1_G3":
    raise ArithmeticError(
        f"unexpected independent q24 status: {back.get('status')}"
    )

def canonical_point_from_exact_payload(rec):
    if rec.get("zero"):
        return E13(0)
    return E13(spec_u(rec["x"]), spec_u(rec["y"]))

canon=g3art["canonical_D13"]
G1=canonical_point_from_exact_payload(canon["G1"])
G3=canonical_point_from_exact_payload(canon["G3"])
assert not G1.is_zero() and not G3.is_zero()

qref=back["q24_modp"]

def eval_mod_coeffs(vals):
    return sum(F(int(v))*tau**i for i,v in enumerate(vals))

rx=(
    eval_mod_coeffs(qref["x_numerator_coefficients_low_to_high"])
    / eval_mod_coeffs(qref["x_denominator_coefficients_low_to_high"])
)
ry=(
    eval_mod_coeffs(qref["y_numerator_coefficients_low_to_high"])
    / eval_mod_coeffs(qref["y_denominator_coefficients_low_to_high"])
)
Qref=E13(rx,ry)
assert not Qref.is_zero()

BOUND=64
solutions=[]
g3_multiples={b:b*G3 for b in range(-BOUND,BOUND+1)}
for asign in (+1,-1):
    AJo=asign*AJ
    for qsign in (+1,-1):
        target=qsign*Qref-AJo
        for acoef in range(-BOUND,BOUND+1):
            rem=target-acoef*G1
            for bcoef,Pb in g3_multiples.items():
                if Pb==rem:
                    solutions.append({
                        "AJ_sign":asign,
                        "q24_sign":qsign,
                        "G1":acoef,
                        "G3":bcoef,
                    })

uniq={}
for sol in solutions:
    key=(sol["AJ_sign"],sol["q24_sign"],sol["G1"],sol["G3"])
    uniq[key]=sol
solutions=list(uniq.values())
solutions.sort(
    key=lambda r:(
        abs(r["G1"])+abs(r["G3"]),
        abs(r["G1"]),abs(r["G3"]),
        -r["q24_sign"],-r["AJ_sign"],
    )
)

print(
    "Q24D46_CORRECTION_SEARCH|"
    f"solutions={len(solutions)}|bound={BOUND}|"
    f"reference={int(rx)},{int(ry)}|"
    f"G1={int(G1[0])},{int(G1[1])}|"
    f"G3={int(G3[0])},{int(G3[1])}|"
    "status=PASS",
    flush=True,
)

for rank,sol in enumerate(solutions[:12],1):
    print(
        "Q24D46_CORRECTION|"
        f"rank={rank}|AJ_sign={sol['AJ_sign']:+d}|"
        f"q24_sign={sol['q24_sign']:+d}|"
        f"add={sol['G1']}*G1+{sol['G3']}*G3|"
        f"l1={abs(sol['G1'])+abs(sol['G3'])}|"
        "status=EXACT_GROUP_IDENTITY",
        flush=True,
    )

if not solutions:
    raise ArithmeticError(
        "degree-46 AJ differs from independent q24 by no small G1/G3 combination"
    )

best=solutions[0]
Qselected=best["AJ_sign"]*AJ + best["G1"]*G1 + best["G3"]*G3
assert Qselected == best["q24_sign"]*Qref
qx,qy=Qselected.xy()

payload={
    "schema":"elkies-k3.h92-q24-degree46-direct-trace-modp.v2",
    "status":"PASS_DEGREE46_DIRECT_Q24_GROUPLAW_CLOSURE",
    "prime":int(p),
    "tau":int(tau),
    "bridge":{
        "q6_old_zero_mw":[-2,-1,-1],
        "q6_standard_mw":[0,-2,-1],
        "formula":"Qmap-S3",
        "q8_degree":46,
        "q6_point":{
            "x":rf_payload(wxq,RQ),
            "y":rf_payload(wyq,RQ),
        },
        "q8_parameter":{
            "numerator_coefficients_low_to_high":[str(v) for v in Unq.list()],
            "denominator_coefficients_low_to_high":[str(v) for v in Udq.list()],
        },
    },
    "direct_AJ":{
        "x":int(ajx),
        "y":int(ajy),
        "W_sign":int(w_sign),
        "method":"II*_E8_1 branch-zero direct L(47O) trace",
    },
    "independent_q24":{
        "x":int(rx),
        "y":int(ry),
        "source":"q8-q24-canonical-backtrack.json",
    },
    "correction_solutions":solutions,
    "selected":{
        **best,
        "formula":(
            f"{best['AJ_sign']}*AJ(Qmap-S3)"
            f"{best['G1']:+d}*G1{best['G3']:+d}*G3"
            f" = {best['q24_sign']}*q24"
        ),
        "q24_x":int(qx),
        "q24_y":int(qy),
    },
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D46_RESULT|"
    f"degree=46|AJ={int(ajx)},{int(ajy)}|"
    f"AJ_sign={best['AJ_sign']:+d}|q24_sign={best['q24_sign']:+d}|"
    f"add={best['G1']}*G1+{best['G3']}*G3|"
    "status=PASS_DEGREE46_DIRECT_Q24_GROUPLAW_CLOSURE",
    flush=True,
)
