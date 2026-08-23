#!/usr/bin/env sage -python
"""
Recover the H92 q24 section EXACTLY over QQ(U) by rational specialization and
interpolation of the proven degree-46 direct trace.

No p-adics. No Hensel. No coefficient reconstruction from modular residues.

Certified construction:
    W = Qmap - S3                  (exact over QQ(T))
    deg_q8(W) = 46
    q24 = AJ_II*(W) + 2*G1

The direct II*_E8_1 branch-zero map is exact over QQ(U).  For many small
rational U=tau, compute q24(tau) exactly over QQ by:

    degree-46 fibre -> direct canonical D13 points -> L(47O) trace -> +2G1.

The known q24 profile
    x(U): 52/48
requires 102 coefficients up to scale, hence 101 exact samples determine x.
We collect a few extra samples as certification.

After rational interpolation:
    D(U) must be Z(U)^2 with deg Z=24,
    RHS = X^3 + A X Z^4 + B Z^6 must be Y(U)^2 with deg Y=78.

All retained exact samples are replayed, then the result is reduced modulo
100003 and required to equal the independently certified degree-46 modular
q24 section coefficient-for-coefficient.

This is a direct exactification of the geometric construction.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, matrix, sage_eval
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
parser.add_argument("--start",type=int,default=2)
parser.add_argument("--samples",type=int,default=108)
parser.add_argument("--scan-limit",type=int,default=300)
parser.add_argument("--sample-output",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
if args.samples < 102:
    raise ValueError("need at least 102 exact samples; 108 is recommended")

ROOT=locate_repo(args.repo)
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"

Q6=GEN/"elkies-k3-h92-q6-child-jacobian.json"
ZERO=GEN/"elkies-k3-h92-q6-child-zero-section.json"
COMP=GEN/"elkies-k3-h92-q6-child-e7-infinity-sections.json"
S3BR=LOCAL/"q6-third-to-q8-bridge.json"
ORIENT=LOCAL/"q8-global-orientation.json"
ANCHOR=LOCAL/"q8-d13-branch-anchor.json"
G3ART=LOCAL/"q8-d13-g3-from-e77-bisection.json"
MOD=LOCAL/"q24-degree46-direct-global-mod-100003.json"

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

for path in (Q6,ZERO,COMP,S3BR,ORIENT,ANCHOR,G3ART,MOD,Q8):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

SAMPLE_OUTPUT=(
    args.sample_output.resolve()
    if args.sample_output
    else LOCAL/"q24-degree46-exact-qq-samples.json"
)
OUTPUT=(
    args.output.resolve()
    if args.output
    else LOCAL/"q8-q24-horizontal-section-qq.json"
)

q6=json.loads(Q6.read_text())
zero=json.loads(ZERO.read_text())
comp=json.loads(COMP.read_text())
s3br=json.loads(S3BR.read_text())
orient=json.loads(ORIENT.read_text())
anchor=json.loads(ANCHOR.read_text())
g3art=json.loads(G3ART.read_text())
modart=json.loads(MOD.read_text())
q8=json.loads(Q8.read_text())

assert q6["status"]=="PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"]=="PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert comp["status"]=="PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert s3br["status"]=="PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert orient["status"]=="PASS_EXACT_Q8_GLOBAL_ORIENTATION"
assert orient["schema"]=="elkies-k3.h92-q8-global-orientation.v2"
assert anchor["status"]=="PASS_EXACT_D13_BRANCH_ANCHOR"
assert g3art["status"]=="PASS_EXACT_D13_G3_FROM_E77_BISECTION"
assert modart["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

global_sign=ZZ(modart["direct_global_map"]["global_scale_root_sign"])
assert global_sign in (-1,1)

# ===========================================================================
# Exact q6 bridge W = Qmap-S3 and exact q8 parameter U(T).
# ===========================================================================

RT=PolynomialRing(QQ,"T")
T=RT.gen()
KT=RT.fraction_field()

def poly(vals):
    return RT([QQ(v) for v in vals])

def rat(data,nk,dk):
    return KT(poly(data[nk]))/KT(poly(data[dk]))

A6=poly(q6["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B6=poly(q6["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
E6=EllipticCurve(KT,[0,0,0,KT(A6),KT(B6)])

zd=zero["section"]
Oold=E6(
    rat(zd,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
    rat(zd,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
)

entries={e["sign"]:e for e in comp["sections"]}
points={
    sign:E6(
        rat(e,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
        rat(e,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
    )
    for sign,e in entries.items()
}
affine=points[comp["source"]["affine_E7_sign"]]
e77=points[comp["source"]["E7_7_sign"]]
Pmap=e77-Oold
Qmap=e77-affine

sd=s3br["third_section_canonical_q6"]
S3=E6(
    rat(sd["x"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
    rat(sd["y"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
)

W=Qmap-S3
assert W in E6 and not W.is_zero()
wx,wy=W.xy()

md=q8["marking"]["section"]
sx=rat(md,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high")
sy=rat(md,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high")
assert E6(sx,sy)==Pmap+Qmap

def monic_power_root(v,e):
    out=v.parent().one()
    for f,m in v.factor():
        assert m%e==0
        out*=f.monic()**(m//e)
    return out.monic()

nx,dx=RT(sx.numerator()),RT(sx.denominator())
ny,dy=RT(sy.numerator()),RT(sy.denominator())
h=monic_power_root(dx,2)
assert h==monic_power_root(dy,3)

ii=RT(next(x for x in q6["finite_fibres"] if x["kodaira"]=="II*")["factor"]).monic()
iv=RT(next(x for x in q6["finite_fibres"] if x["kodaira"]=="IV*")["factor"]).monic()
M=(ii**2*iv**2).monic()

normalizer=(ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
pfun=-sy/sx
rho=(normalizer*nx.inverse_mod(M)).mod(M)

pairs=[]
for entry in q8["rr"]["kernel_polynomials"]:
    sp=RT(entry["s"])
    tp=RT(entry["t"])
    Bc=KT(sp)/KT(h)
    Ac=(
        -KT(sp)*pfun/KT(h)
        -KT(sp)*KT(normalizer)/KT(nx)
        +KT(sp*rho)+KT(tp*M)
    )
    pairs.append((Ac,Bc))
(A0,B0),(A1,B1)=pairs

mW=(wy+sy)/(wx-sx)
UW=KT((A1+B1*mW)/(A0+B0*mW))
Un,Ud=RT(UW.numerator()),RT(UW.denominator())
assert Un.gcd(Ud) in QQ
assert max(Un.degree(),Ud.degree())==46

print(
    "Q24QQ_BRIDGE|formula=Qmap-S3|degree=46|"
    f"x={RT(wx.numerator()).degree()}/{RT(wx.denominator()).degree()}|"
    f"y={RT(wy.numerator()).degree()}/{RT(wy.denominator()).degree()}|"
    "status=PASS",
    flush=True,
)

# ===========================================================================
# Exact global branch-zero map over QQ(U).
# ===========================================================================

RU=PolynomialRing(QQ,"U")
U=RU.gen()
KU=RU.fraction_field()
TU=PolynomialRing(KU,"T")
TT=TU.gen()
KTU=TU.fraction_field()

def ku_from_record(rec):
    return KU(
        RU([QQ(v) for v in rec["numerator_coefficients_low_to_high"]])
    )/KU(
        RU([QQ(v) for v in rec["denominator_coefficients_low_to_high"]])
    )

def tu_from_record(rec):
    return TU([ku_from_record(c) for c in rec["coefficients_low_to_high"]])

def ktu_from_record(rec):
    return KTU(tu_from_record(rec["numerator"]))/KTU(
        tu_from_record(rec["denominator"])
    )

transport=orient["generic_transport"]
quartic_generic=tu_from_record(transport["quartic_in_T_over_QQ_U"])
square_factor_generic=ktu_from_record(transport["square_factor_in_QQ_U_T"])
assert quartic_generic.degree()==4

def parse_u(text):
    return KU(sage_eval(str(text),locals={"U":U}))

tii=QQ(anchor["zero"]["old_base_T"])
coef=anchor["quartic_to_anchor"]["shifted_coefficients"]
aa=parse_u(coef["a_r4"])
bb=parse_u(coef["b_r3"])
cc=parse_u(coef["c_r2"])
dd=parse_u(coef["d_r1"])

rgen=TT-KU(tii)
branch_generic=TU(dd*rgen+cc*rgen**2+bb*rgen**3+aa*rgen**4)
scale_generic=KU(quartic_generic[4]/branch_generic[4])
assert quartic_generic==TU(scale_generic)*branch_generic

def qq_sqrt(q):
    q=QQ(q)
    if q<0:
        return None
    n=ZZ(q.numerator())
    d=ZZ(q.denominator())
    if not n.is_square() or not d.is_square():
        return None
    return QQ(n.sqrt())/QQ(d.sqrt())

def polynomial_sqrt(poly):
    poly=RU(poly)
    if not poly:
        return RU.zero()
    fac=poly.factor()
    ur=qq_sqrt(QQ(fac.unit()))
    if ur is None:
        return None
    out=RU(ur)
    for f,m in fac:
        if m%2:
            return None
        out*=f**(m//2)
    assert out**2==poly
    return out

def rational_sqrt(v):
    v=KU(v)
    nr=polynomial_sqrt(v.numerator())
    dr=polynomial_sqrt(v.denominator())
    if nr is None or dr is None or not dr:
        return None
    out=KU(nr)/KU(dr)
    assert out**2==v
    if QQ(out.numerator().leading_coefficient())<0:
        out=-out
    return out

scale_root=global_sign*rational_sqrt(scale_generic)
if scale_root is None:
    raise ArithmeticError("global quartic scale has no rational square root")

u_iso,rr_iso,ss_iso,tt_iso=[
    parse_u(x) for x in anchor["anchor_to_canonical"]["urst"]
]

A13=RU([QQ(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]])
B13=RU([QQ(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]])

G1rec=g3art["canonical_D13"]["G1"]
G1x=parse_u(G1rec["x"])
G1y=parse_u(G1rec["y"])

print(
    "Q24QQ_MAP|direct_IIstar=PASS|"
    f"scale_root_sign={global_sign:+d}|status=PASS",
    flush=True,
)

# ===========================================================================
# Exact specialization / L(47O) trace.
# ===========================================================================

def spec_KU(v,tau):
    v=KU(v)
    n=RU(v.numerator())(tau)
    d0=RU(v.denominator())(tau)
    if not d0:
        raise ZeroDivisionError("KU pole")
    return QQ(n)/QQ(d0)

def spec_TU(poly,tau):
    poly=TU(poly)
    return RT([spec_KU(c,tau) for c in poly.list()])

def spec_KTU(v,tau):
    v=KTU(v)
    n=spec_TU(v.numerator(),tau)
    d0=spec_TU(v.denominator(),tau)
    if not d0:
        raise ZeroDivisionError("KTU pole")
    return KT(n)/KT(d0)

def newton_power_sums(poly):
    n=poly.degree()
    assert poly[n]==1
    sums=[QQ(n)]
    for k in range(1,n):
        total=QQ(k)*poly[n-k]
        for j in range(1,k):
            total += poly[n-j]*sums[k-j]
        sums.append(-total)
    return sums

def exact_sample(tau):
    tau=QQ(tau)
    a13=QQ(A13(tau))
    b13=QQ(B13(tau))
    E13=EllipticCurve(QQ,[0,0,0,a13,b13])
    if not E13.discriminant():
        raise ArithmeticError("singular D13")

    H=RT(Un-tau*Ud)
    if H.degree()!=46:
        raise ArithmeticError(f"degree drop {H.degree()}")
    H=H.monic()
    if H.gcd(H.derivative()).degree()!=0:
        raise ArithmeticError("non-etale")

    quartic=spec_TU(quartic_generic,tau)
    sqfactor=spec_KTU(square_factor_generic,tau)
    scaleroot=spec_KU(scale_root,tau)
    d=spec_KU(dd,tau)
    uu=spec_KU(u_iso,tau)
    rr=spec_KU(rr_iso,tau)
    ss=spec_KU(ss_iso,tau)
    tt=spec_KU(tt_iso,tau)
    gx=spec_KU(G1x,tau)
    gy=spec_KU(G1y,tau)

    q8m=-(A1-tau*A0)/(B1-tau*B0)

    def modH(v):
        v=KT(v)
        n=RT(v.numerator())
        d0=RT(v.denominator())
        if d0.gcd(H).degree()!=0:
            raise ZeroDivisionError("noninvertible mod H")
        return (n*d0.inverse_mod(H))%H

    msec=(wy+sy)/(wx-sx)
    if modH(q8m-msec):
        raise ArithmeticError("chord mismatch")

    wsec=(2*wx+sx-q8m**2)/sqfactor
    wA=modH(wsec)
    if (wA*wA-quartic)%H:
        raise ArithmeticError("quartic sqrt mismatch")

    wb=wA/scaleroot
    rA=(T-tii)%H
    if rA.gcd(H).degree()!=0:
        raise ArithmeticError("meets branch zero")
    rinv=rA.inverse_mod(H)

    Xa=(d*rinv)%H
    Ya=(d*wb*rinv**2)%H
    xA=((Xa-rr)/(uu**2))%H
    yA=((Ya-ss*(Xa-rr)-tt)/(uu**3))%H

    if (yA*yA-xA*xA*xA-a13*xA-b13)%H:
        raise ArithmeticError("canonical D13 miss")

    xp=[RT.one()]
    for unused in range(23):
        xp.append((xp[-1]*xA)%H)
    cols=list(xp)+[(yA*xp[e])%H for e in range(23)]
    assert len(cols)==47

    Eval=matrix(QQ,46,47,lambda row,col: cols[col][row])
    ker=Eval.right_kernel().basis_matrix()
    if ker.nrows()!=1:
        raise ArithmeticError(f"L47 kernel {ker.nrows()}")
    rel=ker[0]

    XR=PolynomialRing(QQ,"X")
    Xv=XR.gen()
    Afun=sum(rel[i]*Xv**i for i in range(24))
    Bfun=sum(rel[24+i]*Xv**i for i in range(23))
    Rint=Afun**2-(Xv**3+a13*Xv+b13)*Bfun**2
    if Rint.degree()!=47:
        raise ArithmeticError(f"residual degree {Rint.degree()}")

    root_sum=-Rint[46]/Rint[47]
    ps=newton_power_sums(H)
    trace_x=sum(xA[i]*ps[i] for i in range(46))
    xQ=root_sum-trace_x
    if not Bfun(xQ):
        raise ArithmeticError("trace B zero")
    yQ=-Afun(xQ)/Bfun(xQ)

    AJ=-E13(xQ,yQ)
    G1=E13(gx,gy)
    Q24=AJ+2*G1
    if Q24.is_zero():
        raise ArithmeticError("q24 zero")
    qx,qy=Q24.xy()
    return QQ(qx),QQ(qy)

# ===========================================================================
# Resumable exact samples.
# ===========================================================================

samples=[]
attempted=0
skips={}

if SAMPLE_OUTPUT.exists():
    try:
        saved=json.loads(SAMPLE_OUTPUT.read_text())
        if (
            saved.get("schema")=="elkies-k3.h92-q24-exact-qq-samples.v1"
            and int(saved.get("start",-1))==int(args.start)
        ):
            samples=[
                (QQ(row["U"]),QQ(row["x"]),QQ(row["y"]))
                for row in saved.get("samples",[])
            ]
            skips=dict(saved.get("skip_counts",{}))
            attempted=int(saved.get("attempted",0))
            print(
                f"Q24QQ_SAMPLE|resume={len(samples)}|attempted={attempted}|status=PASS",
                flush=True,
            )
    except Exception:
        samples=[]

used={int(u0) for u0,unused_x,unused_y in samples}
candidate=args.start
if used:
    candidate=max(used)+1

def save_samples():
    SAMPLE_OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    SAMPLE_OUTPUT.write_text(json.dumps({
        "schema":"elkies-k3.h92-q24-exact-qq-samples.v1",
        "status":"PARTIAL" if len(samples)<args.samples else "PASS_EXACT_SAMPLES",
        "start":int(args.start),
        "requested":int(args.samples),
        "attempted":int(attempted),
        "skip_counts":skips,
        "samples":[
            {"U":str(u0),"x":str(x0),"y":str(y0)}
            for u0,x0,y0 in samples
        ],
    },indent=2,sort_keys=True)+"\n")

while len(samples)<args.samples and attempted<args.scan_limit:
    tau=ZZ(candidate)
    candidate+=1
    attempted+=1
    if int(tau) in used:
        continue

    try:
        x0,y0=exact_sample(tau)
    except Exception as exc:
        reason=type(exc).__name__+":"+str(exc)
        skips[reason]=skips.get(reason,0)+1
        continue

    samples.append((QQ(tau),x0,y0))
    used.add(int(tau))
    save_samples()

    n=len(samples)
    if n<=5 or n%5==0 or n==args.samples:
        print(
            "Q24QQ_SAMPLE|"
            f"count={n}|U={tau}|"
            f"x_bits={max(abs(ZZ(x0.numerator())).nbits(),abs(ZZ(x0.denominator())).nbits())}|"
            f"y_bits={max(abs(ZZ(y0.numerator())).nbits(),abs(ZZ(y0.denominator())).nbits())}|"
            "status=PASS",
            flush=True,
        )

if len(samples)<args.samples:
    save_samples()
    raise RuntimeError(
        f"only {len(samples)} exact samples after {attempted}; skips={skips}"
    )

print(
    f"Q24QQ_SAMPLE|good={len(samples)}|attempted={attempted}|"
    "stage=collection|status=PASS",
    flush=True,
)

# ===========================================================================
# Exact rational interpolation x=N/D, degrees 52/48.
# ===========================================================================

NUM=52
DEN=48
IM=matrix(
    QQ,
    len(samples),
    (NUM+1)+(DEN+1),
    lambda row,col: (
        samples[row][0]**col
        if col<=NUM
        else -samples[row][1]*samples[row][0]**(col-(NUM+1))
    ),
)
IK=IM.right_kernel().basis_matrix()
if IK.nrows()!=1:
    raise ArithmeticError(f"exact 52/48 interpolation kernel {IK.nrows()}")

rv=IK[0]
X=RU(list(rv[:NUM+1]))
D=RU(list(rv[NUM+1:]))
if not D:
    raise ArithmeticError("x denominator vanished")
scale=D.leading_coefficient()
X/=scale
D/=scale
if X.gcd(D).degree()!=0:
    raise ArithmeticError("x interpolation not reduced")

assert X.degree()==52 and D.degree()==48 and D.is_monic()

for u0,x0,unused_y in samples:
    assert D(u0) and X(u0)/D(u0)==x0

Z=polynomial_sqrt(D)
if Z is None:
    raise ArithmeticError("exact denominator is not a square")
if Z.leading_coefficient()<0:
    Z=-Z
assert Z.degree()==24 and Z.leading_coefficient()==1
assert Z**2==D

RHS=X**3+A13*X*Z**4+B13*Z**6
Y=polynomial_sqrt(RHS)
if Y is None:
    raise ArithmeticError("exact Weierstrass RHS is not a square")
assert Y.degree()==78

direct=opposite=True
for u0,unused_x,y0 in samples:
    pred=Y(u0)/Z(u0)**3
    direct &= pred==y0
    opposite &= -pred==y0
if direct==opposite:
    raise ArithmeticError("exact Y orientation unresolved")
if opposite:
    Y=-Y

assert Y**2==X**3+A13*X*Z**4+B13*Z**6
for u0,x0,y0 in samples:
    assert X(u0)/Z(u0)**2==x0
    assert Y(u0)/Z(u0)**3==y0

print(
    "Q24QQ_INTERP|x=52/48|Z=24|y=78/72|"
    "identity=PASS|samples=PASS|status=PASS_EXACT_SECTION",
    flush=True,
)

# ===========================================================================
# Independent mod-p coefficient check.
# ===========================================================================

p=ZZ(modart["prime"])
F=__import__("sage.all",fromlist=["GF"]).GF(p)
RF=PolynomialRing(F,"U")

def red(poly):
    vals=[]
    for q in poly.list():
        q=QQ(q)
        den=ZZ(q.denominator())
        if den%p==0:
            raise ArithmeticError("exact coefficient denominator divisible by check prime")
        vals.append(F(ZZ(q.numerator()))/F(den))
    return RF(vals)

sec=modart["section_mod_p"]
Zm=RF([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
Xm=RF([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
Ym=RF([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])

assert red(Z)==Zm
assert red(X)==Xm
assert red(Y)==Ym

print(
    "Q24QQ_CROSSCHECK|prime=100003|Z=IDENTICAL|X=IDENTICAL|Y=IDENTICAL|"
    "status=PASS_MODULAR_REDUCTION",
    flush=True,
)

max_num_bits=max(
    abs(ZZ(q.numerator())).nbits()
    for q in list(Z)+list(X)+list(Y)
)
max_den_bits=max(
    abs(ZZ(q.denominator())).nbits()
    for q in list(Z)+list(X)+list(Y)
)

payload={
    "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.v4",
    "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
    "zero":"II*_E8_1_branch_anchor",
    "formula":"Q24 = AJ(Qmap-S3) + 2*G1",
    "bridge":{
        "formula":"Qmap-S3",
        "q6_standard_mw":[0,-2,-1],
        "q8_degree":46,
    },
    "profile":{
        "P_dot_O":24,
        "height":"52",
        "D13_local_correction":"0",
        "Z_degree":24,
        "X_degree":52,
        "Y_degree":78,
        "x_degrees":[52,48],
        "y_degrees":[78,72],
    },
    "section":{
        "Z_coefficients_low_to_high":[str(v) for v in Z.list()],
        "X_coefficients_low_to_high":[str(v) for v in X.list()],
        "Y_coefficients_low_to_high":[str(v) for v in Y.list()],
        "x_numerator_coefficients_low_to_high":[str(v) for v in X.list()],
        "x_denominator_coefficients_low_to_high":[str(v) for v in (Z**2).list()],
        "y_numerator_coefficients_low_to_high":[str(v) for v in Y.list()],
        "y_denominator_coefficients_low_to_high":[str(v) for v in (Z**3).list()],
        "x":"X/Z^2",
        "y":"Y/Z^3",
    },
    "verification":{
        "exact_weierstrass_identity":True,
        "exact_QQ_samples":len(samples),
        "replays_all_exact_samples":True,
        "reduction_matches_degree46_modular_section":True,
    },
    "exactification":{
        "method":"exact QQ specialization of degree-46 direct L(47O) trace + rational interpolation",
        "p_adics_used":False,
        "hensel_used":False,
        "sample_artifact":str(SAMPLE_OUTPUT.relative_to(ROOT)),
        "max_numerator_bits":int(max_num_bits),
        "max_denominator_bits":int(max_den_bits),
    },
    "inputs":{
        "orientation":str(ORIENT.relative_to(ROOT)),
        "branch_anchor":str(ANCHOR.relative_to(ROOT)),
        "modular_crosscheck":str(MOD.relative_to(ROOT)),
    },
    "next":(
        "Construct the exact q24 isotropic divisor on D13, solve its resolved "
        "2-dimensional Riemann-Roch space, eliminate to the genus-one quartic, "
        "and certify the D12/MW5 Jacobian."
    ),
}

OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(f"OUTPUT|{OUTPUT}",flush=True)
print(
    "Q24QQ_RESULT|degree46=1|exact=1|x=52/48|y=78/72|"
    f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
    "status=PASS_EXACT_Q24_FROM_QQ_TRACE_INTERPOLATION",
    flush=True,
)
