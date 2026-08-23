#!/usr/bin/env sage -python
"""
Point the q24 binary quartic at one of its two q24-degree-one old-I9* points
and profile the four native D12 -> A11 targets in that correct group law.

Route pinned by this artifact:
  H3 -> E8+E6 -> D13 --q24--> D24eq/D12
     -> A11 -> 2A5 -> 3A3 -> A3+2A2 -> 5A1 -> 4A1
     -> 3A1 -> 2A1 -> A1 --q6--> rootless/MW17 -> R17.

This script is intentionally only the D12 marking/profile gate.  It does not
change route selection.
"""
import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, lcm, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

p=ZZ(args.prime)
F=GF(p)

SIG=LOCAL/f"q24-orbit85-d12-signature-mod-{p}.json"
SCAN=OUTDIR/f"d12-to-a11-equation-friendly-p{p}.json"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((x for x in q8_candidates if x.exists()),None)

for path in (SIG,SCAN):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")
if Q8 is None:
    raise SystemExit("missing exact q8/D13 child")

sig=json.loads(SIG.read_text())
scan=json.loads(SCAN.read_text())
q8=json.loads(Q8.read_text())

assert sig["status"] in (
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
)
assert scan["status"] in (
    "Q24_D12_A11_NEEDS_SECOND_EXPLICIT_MW_DIRECTION",
    "PASS_Q24_D12_A11_EXPLICIT_MARKED_SECTION",
)

# -------------------------------------------------------------------------
# 1. Old I9* point and q24 quartic over GF(p)(V).
# -------------------------------------------------------------------------
i9=next(x for x in q8["child"]["finite_fibres"] if x["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
UQ=RQ.gen()
fQ=RQ(str(i9["factor"]))
assert fQ.degree()==1
alphaQ=-fQ[0]/fQ[1]
den=ZZ(alphaQ.denominator())%p
assert den
alpha=F(ZZ(alphaQ.numerator())%p)/F(den)

RV=PolynomialRing(F,"V")
V=RV.gen()
K=RV.fraction_field()
RT=PolynomialRing(K,"T")
T=RT.gen()

def rf(rec):
    return K(RV([F(v) for v in rec["num"]]))/K(
        RV([F(v) for v in rec["den"]])
    )

quartic=RT([rf(rec) for rec in sig["quartic_coefficients"]])
assert quartic.degree()==4
jacA=rf(sig["jacobian_A"])
jacB=rf(sig["jacobian_B"])

shifted=RT(quartic(T+K(alpha)))
r2=K(shifted[0])

def poly_sqrt(P):
    P=RV(P)
    if not P:
        return RV.zero()
    fac=P.factor()
    unit=F(fac.unit())
    if not unit.is_square():
        return None
    out=RV(unit.sqrt())
    for f,e in fac:
        if int(e)%2:
            return None
        out*=f**(int(e)//2)
    assert out**2==P
    return out

def rat_sqrt(value):
    value=K(value)
    if not value:
        return K.zero()
    nr=poly_sqrt(value.numerator())
    dr=poly_sqrt(value.denominator())
    if nr is None or dr is None or not dr:
        return None
    out=K(nr)/K(dr)
    assert out**2==value
    return out

r=rat_sqrt(r2)
if r is None or not r:
    raise SystemExit("q24 degree-one quartic value is not a nonzero square")

# -------------------------------------------------------------------------
# 2. Pointed quartic -> Weierstrass, exactly as q32 pointed-spinor method.
#
# Choose (T,w)=(0,+r) as zero.  The opposite point maps to Q.
# -------------------------------------------------------------------------
a=K(shifted[4])
b=K(shifted[3])
c=K(shifted[2])
d=K(shifted[1])

a1=d/r
a2=c-d**2/(K(4)*r**2)
a3=K(2)*r*b
a4=-K(4)*r**2*a
a6=a2*a4

b2=a1**2+4*a2
b4=2*a4+a1*a3
b6=a3**2+4*a6
b8=a1**2*a6+4*a2*a6-a1*a3*a4+a2*a3**2-a4**2
c4p=b2**2-24*b4
c6p=-b2**3+36*b2*b4-216*b6
deltap=-b2**2*b8-8*b4**3-27*b6**2+9*b2*b4*b6
assert deltap

# Opposite q24 degree-one point relative to chosen zero.
xg=-a2
yg=a1*a2-a3
assert yg**2+a1*xg*yg+a3*yg == xg**3+a2*xg**2+a4*xg+a6

xs=xg+b2/K(12)
ys=yg+(a1*xg+a3)/K(2)
shortA=-c4p/K(48)
shortB=-c6p/K(864)
assert ys**2==xs**3+shortA*xs+shortB

# q24 compiler convention has the same pointed elliptic curve after u=3.
xj=9*xs
yj=27*ys
assert jacA==81*shortA
assert jacB==729*shortB
E=EllipticCurve(K,[0,0,0,jacA,jacB])
Q=E(xj,yj)

# Switching the chosen quartic point must negate Q.
rn=-r
na1=d/rn
na2=c-d**2/(K(4)*rn**2)
na3=K(2)*rn*b
nxg=-na2
nyg=na1*na2-na3
nb2=na1**2+4*na2
nxs=nxg+nb2/K(12)
nys=nyg+(na1*nxg+na3)/K(2)
assert 9*nxs==xj and 27*nys==-yj

def norm_rf(value):
    value=K(value)
    n=RV(value.numerator())
    d0=RV(value.denominator())
    lc=d0.leading_coefficient()
    n/=lc
    d0/=lc
    return {
        "num_degree":int(n.degree()),
        "den_degree":int(d0.degree()),
        "num":[int(x) for x in n.list()],
        "den":[int(x) for x in d0.list()],
    }

print(
    "Q24POINTED_D12|"
    f"prime={p}|alpha={int(alpha)}|"
    f"Qxdeg={RV(xj.numerator()).degree()}/{RV(xj.denominator()).degree()}|"
    f"Qydeg={RV(yj.numerator()).degree()}/{RV(yj.denominator()).degree()}|"
    "sign_swap=1|status=PASS_POINTED_Q24_D12",
    flush=True,
)

# -------------------------------------------------------------------------
# 3. R3-zero native D12 frame from the previous scanner.
# -------------------------------------------------------------------------
frame_path=ROOT/scan["frame"]
search_path=ROOT/scan["search_artifact"]
if not frame_path.exists() or not search_path.exists():
    raise SystemExit("missing R3-zero D12 frame/search artifact")

G=matrix(ZZ,[
    [ZZ(v) for v in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert G.dimensions()==(17,17) and G.det()==948

R=G[:12,:12]
C=G[:12,12:]
Tail=G[12:,12:]
H=Tail-C.transpose()*R.inverse()*C
assert H.dimensions()==(5,5)

explicit_z=vector(ZZ,scan["explicit_mw_projection_A0_minus_R3"])
explicit_height=QQ(explicit_z*H*explicit_z)

def class_order(dual):
    o=ZZ(1)
    for x in dual:
        o=lcm(o,ZZ(QQ(x).denominator()))
    return o

def mod2(x):
    x=QQ(x)
    return x-2*(x/2).floor()

def frac_key(v):
    return tuple(QQ(x)-QQ(x).floor() for x in vector(QQ,v))

Rinv=R.inverse()
correction_by_class={frac_key(vector(QQ,[0]*12)):QQ(0)}
for i in range(12):
    weight=vector(QQ,Rinv.row(i))
    key=frac_key(weight)
    norm=QQ(weight*R*weight)
    if key not in correction_by_class or norm<correction_by_class[key]:
        correction_by_class[key]=norm
assert sorted(correction_by_class.values()) == [QQ(0),QQ(1),QQ(3),QQ(3)]

def profile(z):
    z=vector(ZZ,z)
    h=QQ(z*H*z)
    base=vector(ZZ,[0]*12+list(z))
    dual=vector(QQ,base*G[:,:12])*Rinv
    order=class_order(dual)
    key=frac_key(dual)
    if key not in correction_by_class:
        raise ArithmeticError(
            f"unknown D12 discriminant class for {tuple(z)}: {key}"
        )
    corr=correction_by_class[key]
    po=(h+corr-4)/2
    if po not in ZZ or po<0:
        raise ArithmeticError(
            f"invalid D12 profile for {tuple(z)}: h={h}, corr={corr}, P.O={po}"
        )
    return {
        "height":h,
        "class_order":order,
        "local_correction":corr,
        "P_dot_O":po,
        "denominator_Z_degree":po,
        "max_X_degree":2*po+4,
        "max_Y_degree":3*po+6,
    }

qprof=profile(explicit_z)

print(
    "Q24POINTED_MARK|"
    f"lattice_mw={','.join(map(str,explicit_z))}|"
    f"height={explicit_height}|PdotO={qprof['P_dot_O']}|"
    "quartic_opposite_section=plus_or_minus_this_direction|status=PASS",
    flush=True,
)

search=json.loads(search_path.read_text())
assert search["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
a11=[
    rec for rec in search["neighbors"]
    if tuple(rec["child_root_data"])==(11,132,12)
    and rec["child_ade"]=="A11"
    and int(rec["child_mw_rank"])==6
]
assert len(a11)==4

targets=[]
for rec in a11:
    z=vector(ZZ,rec["mw_projection"])
    pr=profile(z)
    row={
        "orbit_index":int(rec["orbit_index"]),
        "mw_projection":list(map(int,z)),
        "dominant_labels":rec["dominant_labels"],
        "witness":rec["witness"],
        "height":str(pr["height"]),
        "class_order":int(pr["class_order"]),
        "local_correction":str(pr["local_correction"]),
        "P_dot_O":int(pr["P_dot_O"]),
        "denominator_Z_degree":int(pr["denominator_Z_degree"]),
        "max_X_degree":int(pr["max_X_degree"]),
        "max_Y_degree":int(pr["max_Y_degree"]),
    }
    targets.append(row)
    print(
        "Q24POINTED_A11_PROFILE|"
        f"orbit={rec['orbit_index']}|mw={','.join(map(str,z))}|"
        f"height={pr['height']}|order={pr['class_order']}|"
        f"corr={pr['local_correction']}|PdotO={pr['P_dot_O']}|"
        f"degrees=Z{pr['denominator_Z_degree']},X<={pr['max_X_degree']},Y<={pr['max_Y_degree']}|"
        "status=PASS",
        flush=True,
    )

targets.sort(key=lambda row:(
    row["P_dot_O"],
    sum(abs(x) for x in row["mw_projection"]),
    row["orbit_index"],
))
best=targets[0]

route=[
    {"from":"H3","to":"E8+E6","q":6},
    {"from":"E8+E6","to":"D13","q":8},
    {"from":"D13","to":"D24eq/D12","q":24},
    {"from":"D12","to":"A11","q":6},
    {"from":"A11","to":"2A5","q":8},
    {"from":"2A5","to":"3A3","q":4},
    {"from":"3A3","to":"A3+2A2","q":4},
    {"from":"A3+2A2","to":"5A1","q":4},
    {"from":"5A1","to":"4A1","q":4},
    {"from":"4A1","to":"3A1","q":4},
    {"from":"3A1","to":"2A1","q":4},
    {"from":"2A1","to":"A1","q":4},
    {"from":"A1","to":"rootless/MW17","q":6},
    {"from":"rootless/MW17","to":"R17","q":None,"kind":"pinned lattice isometry"},
]

payload={
    "schema":"elkies-k3.h3-q24-pointed-d12-a11-profile.v1",
    "status":"PASS_Q24_POINTED_D12_A11_PROFILE",
    "prime":int(p),
    "route":route,
    "pointed_quartic":{
        "chosen_zero":"one of {R3,A0}; sign choice only negates the displayed opposite section",
        "old_I9star_root_QQ":str(alphaQ),
        "old_I9star_root_mod_p":int(alpha),
        "spinor_sqrt":norm_rf(r),
        "opposite_section":{
            "x":norm_rf(xj),
            "y":norm_rf(yj),
        },
        "compiler_jacobian_A":norm_rf(jacA),
        "compiler_jacobian_B":norm_rf(jacB),
        "sign_swap_negates_section":True,
    },
    "R3_zero_lattice_marking":{
        "explicit_A0_minus_R3_mw":list(map(int,explicit_z)),
        "height":str(explicit_height),
        "profile":{
            k:(str(v) if k in ("height","local_correction") else int(v))
            for k,v in qprof.items()
        },
    },
    "A11_targets":targets,
    "selected_smallest_profile":best,
    "next_solver_contract":{
        "target_orbit":best["orbit_index"],
        "P_dot_O":best["P_dot_O"],
        "Z_degree":best["denominator_Z_degree"],
        "X_degree_max":best["max_X_degree"],
        "Y_degree_max":best["max_Y_degree"],
        "require_child_root_data":[11,132,12],
        "require_child":"A11/MW6",
    },
    "proof_boundary":(
        "This fixes the q24 quartic-to-D12 group law by a pointed quartic "
        "conversion and profiles the four native A11 marked sections exactly "
        "in the R3-zero D12 lattice. It does not yet solve the target section "
        "coordinates or compile the q6 A11 pencil."
    ),
}

OUT=args.output.resolve() if args.output else OUTDIR/f"pointed-d12-a11-profile-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q24POINTED_A11_BEST|"
    f"orbit={best['orbit_index']}|PdotO={best['P_dot_O']}|"
    f"Zdeg={best['denominator_Z_degree']}|"
    f"Xmax={best['max_X_degree']}|Ymax={best['max_Y_degree']}|status=SELECTED",
    flush=True,
)
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24POINTED_RESULT|"
    "route_end=R17|A11_profiles=4|"
    "status=PASS_Q24_POINTED_D12_A11_PROFILE",
    flush=True,
)
