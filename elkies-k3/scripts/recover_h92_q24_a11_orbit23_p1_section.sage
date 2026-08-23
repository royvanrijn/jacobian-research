#!/usr/bin/env sage -python
"""
Recover the orbit-23 D12 -> A11 marked section on the correctly pointed q24
D12 model over GF(p).

Use the cheap translate
    S = P23 - Q,
where Q=A0-R3 is the explicit zero-pole section.  Lattice profile:
    S.O=1, height=5, D12 correction=1.

On the I8*-at-infinity model this forces
    Z=u-z0, deg X=5, deg Y=7,
with the u^15 term cancelling by
    c^3 + A6*c + B9 = 0
for the leading coefficient c of X.

Enumerate those small rational sections, identify their exact MW classes by
group-law pole signatures relative to Q, select S=P23-Q, then form P23=S+Q.
"""
import argparse, json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, IntegralLattice, PolynomialRing, QQ, ZZ,
    lcm, matrix, pari, vector
)

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)
F=GF(p)

EASY=OUTDIR/f"a11-easy-spinor-shift-p{p}.json"
PROF=OUTDIR/f"pointed-d12-a11-profile-p{p}.json"
SCAN=OUTDIR/f"d12-to-a11-equation-friendly-p{p}.json"
for path in (EASY,PROF,SCAN):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

easy=json.loads(EASY.read_text())
prof=json.loads(PROF.read_text())
scan=json.loads(SCAN.read_text())
assert easy["status"]=="Q24_A11_NO_ZERO_POLE_TRANSLATION_BY_KNOWN_Q"
assert prof["status"]=="PASS_Q24_POINTED_D12_A11_PROFILE"

# -------------------------------------------------------------------------
# 1. Correct pointed polynomial D12 model with I8* at infinity.
# -------------------------------------------------------------------------
model=easy["I8star_infinity"]
RU=PolynomialRing(F,"u")
u=RU.gen()
K=RU.fraction_field()

def poly(rec):
    return RU([F(v) for v in rec["coefficients_low_to_high"]])

A=poly(model["A"])
B=poly(model["B"])
Qx=poly(model["known_Q_x"])
Qy=poly(model["known_Q_y"])
assert A.degree()==6 and B.degree()==9
assert Qx.degree()<=4 and Qy.degree()<=6
assert Qy**2==Qx**3+A*Qx+B

E=EllipticCurve(K,[0,0,0,K(A),K(B)])
Q=E(K(Qx),K(Qy))

print(
    "Q24A11P1_MODEL|"
    f"prime={p}|Adeg={A.degree()}|Bdeg={B.degree()}|"
    f"Qxdeg={Qx.degree()}|Qydeg={Qy.degree()}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Exact R3-zero D12 lattice and exact local-correction computation.
# -------------------------------------------------------------------------
frame_path=ROOT/scan["frame"]
G=matrix(ZZ,[
    [ZZ(v) for v in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert G.dimensions()==(17,17) and G.det()==948

root=G[:12,:12]
coupling=G[:12,12:]
tail=G[12:,12:]
H=tail-coupling.transpose()*root.inverse()*coupling
assert H.dimensions()==(5,5)

q=vector(ZZ,prof["R3_zero_lattice_marking"]["explicit_A0_minus_R3_mw"])
t23=vector(
    ZZ,
    next(row for row in prof["A11_targets"] if int(row["orbit_index"])==23)["mw_projection"]
)
target=t23-q
assert tuple(target)==(0,-1,1,-1,0)

# D12 discriminant-class correction table.
#
# `dual = pairing * root^{-1}` is expressed in the SIMPLE-ROOT coefficient
# basis.  The previous version incorrectly fed this coefficient vector to
# IntegralLattice.enumerate_close_vectors(), whose vectors live in a different
# ambient realization.
#
# For D12 the discriminant group has four classes with minimal corrections
# 0,1,3,3.  Recover those classes intrinsically from the inverse Cartan matrix;
# this is invariant under permutation of the simple-root numbering.
root_inv=root.inverse()

def frac_key(v):
    return tuple(
        QQ(x)-QQ(x).floor()
        for x in vector(QQ,v)
    )

correction_by_class={frac_key(vector(QQ,[0]*12)):QQ(0)}
for i in range(12):
    weight=vector(QQ,root_inv.row(i))
    key=frac_key(weight)
    norm=QQ(weight*root*weight)
    if key not in correction_by_class or norm<correction_by_class[key]:
        correction_by_class[key]=norm

if sorted(correction_by_class.values()) != [QQ(0),QQ(1),QQ(3),QQ(3)]:
    raise ArithmeticError(
        f"unexpected D12 discriminant corrections: {correction_by_class}"
    )

def exact_correction(z):
    z=vector(ZZ,z)
    base=vector(ZZ,[0]*12+list(z))
    pair=vector(QQ,base*G[:,:12])
    dual=pair*root_inv
    key=frac_key(dual)
    try:
        return correction_by_class[key]
    except KeyError as error:
        raise ArithmeticError(
            f"unrecognized D12 discriminant class for {tuple(z)}: {key}"
        ) from error

def lattice_profile(z):
    z=vector(ZZ,z)
    h=QQ(z*H*z)
    corr=exact_correction(z)
    po=(h+corr-4)/2
    if po not in ZZ or po<0:
        raise ArithmeticError((tuple(z),h,corr,po))
    return h,corr,ZZ(po)

th,tc,tpo=lattice_profile(target)
print(
    "Q24A11P1_TARGET_CHECK|"
    f"height={th}|corr={tc}|PdotO={tpo}|status="
    f"{'PASS' if (th,tc,tpo)==(QQ(5),QQ(1),ZZ(1)) else 'MISMATCH'}",
    flush=True,
)
assert (th,tc,tpo)==(QQ(5),QQ(1),ZZ(1))

print(
    "Q24A11P1_TARGET|"
    f"mw={','.join(map(str,target))}|height={th}|corr={tc}|PdotO={tpo}|"
    "status=PASS",
    flush=True,
)

# Enumerate all abstract height-5/correction-1 MW vectors for signature matching.
scale=ZZ(1)
for value in H.list():
    scale=lcm(scale,ZZ(QQ(value).denominator()))
IH=(scale*H).change_ring(ZZ)
qf=pari(IH).qfminim(ZZ(scale*5))

abstract={}
seen=set()
for col in matrix(ZZ,qf[2]).columns():
    for sign in (1,-1):
        z=sign*vector(ZZ,col)
        key=tuple(map(int,z))
        if key in seen:
            continue
        seen.add(key)
        h=QQ(z*H*z)
        if h!=5:
            continue
        corr=exact_correction(z)
        if corr!=1:
            continue
        _,_,po=lattice_profile(z)
        if po==1:
            abstract[key]=z

assert tuple(target) in abstract
print(
    "Q24A11P1_ABSTRACT|"
    f"height5_corr1_PdotO1={len(abstract)}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 3. Pole degree / signature helpers in the pointed function-field group.
# -------------------------------------------------------------------------
def monic_power_root(P,e):
    P=RU(P)
    if not P:
        return RU.zero()
    lc=P.leading_coefficient()
    P=P/lc
    out=RU.one()
    for fac,m in P.factor():
        if int(m)%e:
            return None
        out*=fac.monic()**(int(m)//e)
    return out.monic()

def point_pole_degree(P):
    if P.is_zero():
        return None
    x,y=P.xy()
    dx=RU(x.denominator())
    dy=RU(y.denominator())
    zx=monic_power_root(dx,2)
    zy=monic_power_root(dy,3)
    if zx is None or zy is None or zx!=zy:
        raise ArithmeticError("section denominators are not a common Z^2/Z^3")
    return ZZ(zx.degree())

def lattice_signature(z):
    z=vector(ZZ,z)
    values=[]
    for k in (1,2):
        for n in range(-3,4):
            _,_,po=lattice_profile(k*z+n*q)
            values.append(int(po))
    return tuple(values)

def point_signature(P):
    values=[]
    for k in (1,2):
        for n in range(-3,4):
            R=k*P+n*Q
            po=point_pole_degree(R)
            if po is None:
                # Zero section corresponds to no denominator; use sentinel -1.
                values.append(-1)
            else:
                values.append(int(po))
    return tuple(values)

abstract_by_sig={}
for z in abstract.values():
    sig=lattice_signature(z)
    abstract_by_sig.setdefault(sig,[]).append(z)

target_sig=lattice_signature(target)

# -------------------------------------------------------------------------
# 4. Solve all correction-1 P.O=1 rational sections.
#
# Z=u-z0, X=c*u^5+x4*u^4+...+x0, Y degree 7.
# The u^15 coefficient of the cleared equation must vanish.
# -------------------------------------------------------------------------
lead_poly=PolynomialRing(F,"c")
cc=lead_poly.gen()
lead_eq=cc**3+F(A[6])*cc+F(B[9])
lead_roots=[]
for fac,e in lead_eq.factor():
    if fac.degree()==1:
        lead_roots.append(F(-fac[0]/fac[1]))
lead_roots=sorted(set(lead_roots),key=int)
if not lead_roots:
    raise SystemExit(f"no rational correction-1 leading branch: {lead_eq.factor()}")

print(
    "Q24A11P1_LEADS|"
    f"count={len(lead_roots)}|values={','.join(str(int(v)) for v in lead_roots)}|"
    "status=PASS",
    flush=True,
)

def solve_branch(clead):
    names=("z0","ell","x4","x3","x2","x1","x0","inv")
    S=PolynomialRing(F,names=names,order="degrevlex")
    z0,ell,x4,x3,x2,x1,x0,inv=S.gens()
    KF=S.fraction_field()
    UR=PolynomialRing(KF,"U")
    U=UR.gen()

    AA=UR([KF(v) for v in A.list()])
    BB=UR([KF(v) for v in B.list()])
    Z=U-KF(z0)
    X=KF(clead)*U**5+KF(x4)*U**4+KF(x3)*U**3+KF(x2)*U**2+KF(x1)*U+KF(x0)
    rhs=X**3+AA*X*Z**4+BB*Z**6
    assert rhs[15]==0
    assert rhs.degree()<=14

    ys={7:KF(ell)}
    equations=[S((KF(ell)**2-KF(rhs[14])).numerator())]
    for degree in range(13,6,-1):
        j=degree-7
        known=sum(
            ys[i]*ys[degree-i]
            for i in ys
            if (degree-i) in ys and i!=7 and (degree-i)!=7
        )
        ys[j]=(KF(rhs[degree])-known)/(KF(2)*ys[7])

    Y=sum(ys[i]*U**i for i in range(8))
    residual=Y**2-rhs
    for degree in range(7):
        equations.append(S(KF(residual[degree]).numerator()))
    equations.append(inv*ell-1)

    I=S.ideal(equations)
    print(
        "Q24A11P1_BRANCH|"
        f"lead={int(clead)}|vars={S.ngens()}|eqs={len(equations)}|"
        "status=GROEBNER_START",
        flush=True,
    )
    sols=I.variety()
    print(
        "Q24A11P1_BRANCH|"
        f"lead={int(clead)}|solutions={len(sols)}|status=GROEBNER_PASS",
        flush=True,
    )

    out=[]
    for sol in sols:
        vals={g:F(sol[g]) for g in S.gens()}
        if not vals[ell]:
            continue
        Zp=RU.base_ring() if False else None
        Z0=RU  # silence lint-like readers
        # Rebuild over F[u].
        Zpoly=RU  # replaced immediately below
        # Numeric polynomial ring is the outer RU from the script.
        Znum=u-vals[z0]
        Xnum=RU([
            vals[x0],vals[x1],vals[x2],vals[x3],vals[x4],F(clead)
        ])
        rhsnum=Xnum**3+A*Xnum*Znum**4+B*Znum**6
        yn=[F(0)]*8
        yn[7]=vals[ell]
        for degree in range(13,6,-1):
            j=degree-7
            known=sum(
                yn[i]*yn[degree-i]
                for i in range(8)
                if 0<=degree-i<8 and i!=7 and (degree-i)!=7
            )
            yn[j]=(rhsnum[degree]-known)/(F(2)*yn[7])
        Ynum=RU(yn)
        assert Ynum**2==rhsnum
        P=E(K(Xnum)/K(Znum**2),K(Ynum)/K(Znum**3))
        assert point_pole_degree(P)==1
        out.append((P,Znum,Xnum,Ynum))
    return out

sections=[]
for clead in lead_roots:
    sections.extend(solve_branch(clead))

# Deduplicate by exact rational coordinates.
unique={}
for P,Z,X,Y in sections:
    x,y=P.xy()
    key=(str(x),str(y))
    unique[key]=(P,Z,X,Y)
sections=list(unique.values())

print(
    "Q24A11P1_COUNT|"
    f"sections={len(sections)}|abstract={len(abstract)}|status="
    f"{'PASS_COUNT' if len(sections)==len(abstract) else 'PARTIAL_COUNT'}",
    flush=True,
)

# -------------------------------------------------------------------------
# 5. Match section signatures to exact MW vectors.
# -------------------------------------------------------------------------
matched=[]
unmatched=[]
for P,Z,X,Y in sections:
    sig=point_signature(P)
    candidates=abstract_by_sig.get(sig,[])
    row={
        "signature":list(sig),
        "candidate_mw":[list(map(int,z)) for z in candidates],
        "Z":[int(v) for v in Z.list()],
        "X":[int(v) for v in X.list()],
        "Y":[int(v) for v in Y.list()],
    }
    if len(candidates)==1:
        z=candidates[0]
        row["mw"]=list(map(int,z))
        matched.append((z,P,Z,X,Y,row))
    else:
        unmatched.append(row)

print(
    "Q24A11P1_MATCH|"
    f"unique={len(matched)}|ambiguous={len(unmatched)}|"
    f"target_signature_candidates={len(abstract_by_sig.get(target_sig,[]))}|"
    "status=PASS_SIGNATURE_AUDIT",
    flush=True,
)

target_hits=[item for item in matched if item[0]==target]
if len(target_hits)!=1:
    # Try direct signature if its abstract signature is unique and one equation
    # point has that signature, even if other abstract classes collide elsewhere.
    abs_target=abstract_by_sig.get(target_sig,[])
    pts=[item for item in sections if point_signature(item[0])==target_sig]
    if len(abs_target)==1 and abs_target[0]==target and len(pts)==1:
        P,Z,X,Y=pts[0]
        target_hits=[(target,P,Z,X,Y,{
            "mw":list(map(int,target)),
            "signature":list(target_sig),
            "Z":[int(v) for v in Z.list()],
            "X":[int(v) for v in X.list()],
            "Y":[int(v) for v in Y.list()],
        })]

if len(target_hits)!=1:
    payload={
        "schema":"elkies-k3.h3-q24-a11-orbit23-p1-section.v1",
        "status":"Q24_A11_P1_TARGET_NOT_UNIQUELY_IDENTIFIED",
        "prime":int(p),
        "target_mw":list(map(int,target)),
        "abstract_count":len(abstract),
        "section_count":len(sections),
        "matched_count":len(matched),
        "ambiguous":unmatched,
    }
    OUT=args.output.resolve() if args.output else OUTDIR/f"a11-orbit23-p1-section-p{p}.json"
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print(
        "Q24A11P1_RESULT|"
        f"sections={len(sections)}|matched={len(matched)}|target=0|"
        f"status={payload['status']}",
        flush=True,
    )
    raise SystemExit(0)

_,Spt,SZ,SX,SY,target_row=target_hits[0]

# Actual orbit-23 section P23=S+Q.
P23=Spt+Q
assert not P23.is_zero()
assert point_pole_degree(P23)==2
x23,y23=P23.xy()
assert y23**2==x23**3+K(A)*x23+K(B)

def rf_record(v):
    v=K(v)
    n=RU(v.numerator())
    d=RU(v.denominator())
    lc=d.leading_coefficient()
    n=n/lc
    d=d/lc
    return {
        "num_degree":int(n.degree()),
        "den_degree":int(d.degree()),
        "num":[int(x) for x in n.list()],
        "den":[int(x) for x in d.list()],
    }

payload={
    "schema":"elkies-k3.h3-q24-a11-orbit23-p1-section.v1",
    "status":"PASS_Q24_A11_ORBIT23_MARKED_SECTION_MODP",
    "prime":int(p),
    "route_end":"R17",
    "known_Q_mw":list(map(int,q)),
    "translated_target":{
        "mw":list(map(int,target)),
        "height":"5",
        "local_correction":"1",
        "P_dot_O":1,
        "Z":[int(v) for v in SZ.list()],
        "X":[int(v) for v in SX.list()],
        "Y":[int(v) for v in SY.list()],
        "x":rf_record(Spt.xy()[0]),
        "y":rf_record(Spt.xy()[1]),
    },
    "orbit23":{
        "orbit_index":23,
        "mw":list(map(int,t23)),
        "construction":"P23=(P23-Q)+Q",
        "P_dot_O":2,
        "x":rf_record(x23),
        "y":rf_record(y23),
        "expected_child":"A11/MW6",
        "expected_root_data":[11,132,12],
    },
    "enumeration":{
        "abstract_height5_corr1_PdotO1":len(abstract),
        "equation_sections":len(sections),
        "uniquely_signature_matched":len(matched),
    },
    "proof_boundary":(
        "The orbit-23 marked section is recovered and identified modulo p on "
        "the correctly pointed q24 D12 model. The q6 resolved RR pencil and "
        "A11 child equation are not yet compiled."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/f"a11-orbit23-p1-section-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q24A11P1_TARGET_FOUND|"
    f"mw={','.join(map(str,target))}|Zdeg={SZ.degree()}|"
    f"orbit23_xden={x23.denominator().degree()}|"
    "status=PASS",
    flush=True,
)
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24A11P1_RESULT|"
    f"sections={len(sections)}|matched={len(matched)}|target=1|"
    "status=PASS_Q24_A11_ORBIT23_MARKED_SECTION_MODP",
    flush=True,
)
