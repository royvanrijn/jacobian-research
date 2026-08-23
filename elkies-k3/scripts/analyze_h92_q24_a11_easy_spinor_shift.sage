#!/usr/bin/env sage -python
"""
Normalize the pointed q24 D12 model at its I8* fibre and look for easy
zero-pole translations of the four A11 target sections.

For each target P and n in [-4,4], profile P+n*Q where
Q=A0-R3 is the already-explicit P.O=0 pointed-spinor section.
"""
import argparse, json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, binomial, lcm, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)
F=GF(p)

PROF=OUTDIR/f"pointed-d12-a11-profile-p{p}.json"
if not PROF.exists():
    raise SystemExit(f"missing {PROF}")
d=json.loads(PROF.read_text())
assert d["status"]=="PASS_Q24_POINTED_D12_A11_PROFILE"

RV=PolynomialRing(F,"V")
V=RV.gen()
K=RV.fraction_field()

def rf(rec):
    return K(RV([F(v) for v in rec["num"]]))/K(
        RV([F(v) for v in rec["den"]])
    )

pt=d["pointed_quartic"]
A=rf(pt["compiler_jacobian_A"])
B=rf(pt["compiler_jacobian_B"])
Qx=rf(pt["opposite_section"]["x"])
Qy=rf(pt["opposite_section"]["y"])
assert Qy**2==Qx**3+A*Qx+B

def exact_root(P,e):
    P=RV(P)
    out=RV.one()
    for fac,m in P.factor():
        m=int(m)
        if m%e:
            return None
        out*=fac.monic()**(m//e)
    return out.monic()

# Normalize common line-bundle denominator.
An=RV(A.numerator()); Ad=RV(A.denominator())
Bn=RV(B.numerator()); Bd=RV(B.denominator())
Xn=RV(Qx.numerator()); Xd=RV(Qx.denominator())
Yn=RV(Qy.numerator()); Yd=RV(Qy.denominator())

# Strip denominator units first.
for n0,d0,name in ((An,Ad,"A"),(Bn,Bd,"B"),(Xn,Xd,"X"),(Yn,Yd,"Y")):
    if not d0:
        raise ArithmeticError(f"{name} denominator zero")

def monic_den(n0,d0):
    lc=d0.leading_coefficient()
    return RV(n0/lc),RV(d0/lc)

An,Ad=monic_den(An,Ad)
Bn,Bd=monic_den(Bn,Bd)
Xn,Xd=monic_den(Xn,Xd)
Yn,Yd=monic_den(Yn,Yd)

h=exact_root(Ad,4)
assert h is not None
assert h==exact_root(Bd,6)==exact_root(Xd,2)==exact_root(Yd,3)
assert Ad==h**4 and Bd==h**6 and Xd==h**2 and Yd==h**3

# Polynomial numerator model after clearing line-bundle h.
A0=An; B0=Bn; X0=Xn; Y0=Yn
assert Y0**2==X0**3+A0*X0+B0

Delta=-16*(4*A0**3+27*B0**2)
stars=[(fac,int(e)) for fac,e in Delta.factor() if int(e)==14 and fac.degree()==1]
if len(stars)!=1:
    raise ArithmeticError(f"expected unique linear I8*^14 factor, got {Delta.factor()}")
fstar=stars[0][0].monic()
beta=-fstar[0]

T=RV.gen()

def infinity_transform(P,w):
    P=RV(P)
    out=RV.zero()
    for i,c in enumerate(P.list()):
        if not c:
            continue
        for j in range(i+1):
            exponent=w-j
            if exponent<0:
                raise ArithmeticError((P.degree(),w,i,j))
            out += c*F(binomial(i,j))*beta**(i-j)*T**exponent
    return RV(out)

Ai=infinity_transform(A0,8)
Bi=infinity_transform(B0,12)
Xi=infinity_transform(X0,4)
Yi=infinity_transform(Y0,6)

assert Ai.degree()<=6 and Bi.degree()<=9
assert Xi.degree()<=4 and Yi.degree()<=6
assert Yi**2==Xi**3+Ai*Xi+Bi

print(
    "Q24POINTED_I8INF|"
    f"prime={p}|beta={int(beta)}|"
    f"Adeg={Ai.degree()}|Bdeg={Bi.degree()}|"
    f"Qxdeg={Xi.degree()}|Qydeg={Yi.degree()}|"
    "status=PASS_POLYNOMIAL_Q",
    flush=True,
)

# Lattice profile.
frame_path=ROOT/d["R3_zero_lattice_marking"].get(
    "frame","artifacts/local/elkies-k3/q24-downstream-lift/d12-c10a-zero-frame.txt"
)
# profile artifact currently does not embed frame; use original scan.
SCAN=OUTDIR/f"d12-to-a11-equation-friendly-p{p}.json"
scan=json.loads(SCAN.read_text())
frame_path=ROOT/scan["frame"]

G=matrix(ZZ,[
    [ZZ(v) for v in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
R=G[:12,:12]
C=G[:12,12:]
Tail=G[12:,12:]
H=Tail-C.transpose()*R.inverse()*C
q=vector(ZZ,d["R3_zero_lattice_marking"]["explicit_A0_minus_R3_mw"])

def class_order(dual):
    o=ZZ(1)
    for x in dual:
        o=lcm(o,ZZ(QQ(x).denominator()))
    return o

def mod2(x):
    x=QQ(x)
    return x-2*(x/2).floor()

def profile(z):
    z=vector(ZZ,z)
    hgt=QQ(z*H*z)
    base=vector(ZZ,[0]*12+list(z))
    dual=vector(QQ,base*G[:,:12])*R.inverse()
    order=class_order(dual)
    raw=QQ(dual*R*dual)
    poss=[]
    for corr in (QQ(0),QQ(1),QQ(3)):
        if mod2(corr)!=mod2(raw):
            continue
        po=(hgt+corr-4)/2
        if po in ZZ and po>=0:
            poss.append((ZZ(po),corr))
    if len(poss)>1 and order==2:
        v=[x for x in poss if x[1]==1]
        if v: poss=v
    if len(poss)!=1:
        return {"height":hgt,"order":order,"possible":poss}
    po,corr=poss[0]
    return {
        "height":hgt,"order":order,"correction":corr,"P_dot_O":po,
    }

rows=[]
easy=[]
for target in d["A11_targets"]:
    t=vector(ZZ,target["mw_projection"])
    for n in range(-4,5):
        z=t+n*q
        pr=profile(z)
        row={
            "orbit_index":target["orbit_index"],
            "n":n,
            "mw":list(map(int,z)),
            "height":str(pr["height"]),
            "class_order":int(pr["order"]),
            "P_dot_O":None if "P_dot_O" not in pr else int(pr["P_dot_O"]),
            "correction":None if "correction" not in pr else str(pr["correction"]),
        }
        rows.append(row)
        if row["P_dot_O"]==0:
            easy.append(row)
            print(
                "Q24A11_EASYSHIFT|"
                f"orbit={target['orbit_index']}|n={n}|"
                f"mw={','.join(map(str,z))}|height={pr['height']}|"
                f"corr={pr['correction']}|status=PASS_ZERO_POLE",
                flush=True,
            )

if not easy:
    # Still report best shifts by pole order.
    good=[r for r in rows if r["P_dot_O"] is not None]
    good.sort(key=lambda r:(r["P_dot_O"],abs(r["n"]),r["orbit_index"]))
    for r in good[:8]:
        print(
            "Q24A11_SHIFTBEST|"
            f"orbit={r['orbit_index']}|n={r['n']}|"
            f"PdotO={r['P_dot_O']}|height={r['height']}|corr={r['correction']}|"
            f"mw={','.join(map(str,r['mw']))}|status=PROFILE",
            flush=True,
        )

def polyrec(P):
    P=RV(P)
    return {
        "degree":int(P.degree()),
        "coefficients_low_to_high":[int(x) for x in P.list()],
    }

payload={
    "schema":"elkies-k3.h3-q24-a11-easy-spinor-shift.v1",
    "status":(
        "PASS_Q24_A11_ZERO_POLE_TRANSLATION"
        if easy else
        "Q24_A11_NO_ZERO_POLE_TRANSLATION_BY_KNOWN_Q"
    ),
    "prime":int(p),
    "I8star_infinity":{
        "beta":int(beta),
        "A":polyrec(Ai),
        "B":polyrec(Bi),
        "known_Q_x":polyrec(Xi),
        "known_Q_y":polyrec(Yi),
    },
    "known_Q_mw":list(map(int,q)),
    "profiles":rows,
    "zero_pole_translations":easy,
    "route_end":"R17",
}
OUT=args.output.resolve() if args.output else OUTDIR/f"a11-easy-spinor-shift-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24A11_EASYSHIFT_RESULT|"
    f"zero_pole={len(easy)}|status={payload['status']}",
    flush=True,
)
