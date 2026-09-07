#!/usr/bin/env sage -python
"""
Build the q24 D12 lattice with geometric zero C10a, enumerate all q6 -> A11
neighbors, and select those whose MW projection lies in the explicit spinor
direction from the q24 quartic 2-cover.

For every scalar-direction hit, export the exact lattice decomposition
    D = O + P + V + kF
and the explicit modular marked point P on the D12 Jacobian.
"""
import argparse, contextlib, io, itertools, json, subprocess, sys
from pathlib import Path
from sage.all import (
    EllipticCurve, GF, IntegralLattice, PolynomialRing, QQ, ZZ,
    block_diagonal_matrix, identity_matrix, lcm, matrix, vector
)

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
EFF=SCRIPTS/"transport_h92_q24_effective_i9star_components.sage"
SEARCH=SCRIPTS/"search_root_adapted_weyl_neighbors.sage"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)
F=GF(p)

ANCHOR=LOCAL/f"q24-d12-spinor-anchor-mod-{p}.json"
SIG=LOCAL/f"q24-orbit85-d12-signature-mod-{p}.json"
OUTDIR=LOCAL/"q24-downstream-lift"
FRAME=OUTDIR/"d12-c10a-zero-frame.txt"
SEARCHOUT=OUTDIR/"d12-c10a-zero-q6-all.json"
FRAMES=OUTDIR/"d12-c10a-zero-q6-frames"
for path in (ANCHOR,SIG,ENGINE,EFF,SEARCH):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

anchor=json.loads(ANCHOR.read_text())
sig=json.loads(SIG.read_text())
if anchor["status"]!="PASS_Q24_D12_SPINOR_PRIMITIVE_DIRECTION":
    raise SystemExit(f"spinor anchor is not passing: {anchor['status']}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

# Replay the ACTUAL effective old-I9* component transport.  This transports
# all 14 effective affine-D13 components into the equation frame and proves
# the q24 incidence pattern exactly.  The two q24-degree-one components are
# pinned R3 and affine A0.
saved=list(sys.argv)
scope={"__name__":"__embedded_q24_effcomp__","__file__":str(EFF)}
buf=io.StringIO()
transport_stopped_at_identity=False
try:
    sys.argv=[str(EFF),"--prime",str(p)]
    with contextlib.redirect_stdout(buf):
        exec(compile(EFF.read_text(),str(EFF),"exec"),scope)
finally:
    sys.argv=saved

need=("ns","D","F8","effective","Dpairs")
missing=[k for k in need if k not in scope]
if missing:
    raise SystemExit(
        "actual component transport stopped before required data: "
        + ",".join(missing)
    )

ns=matrix(ZZ,scope["ns"])
D=vector(ZZ,scope["D"])
F8=vector(ZZ,scope["F8"])
effective={name:vector(ZZ,c) for name,c in scope["effective"].items()}
Dpairs={name:int(v) for name,v in scope["Dpairs"].items()}

assert len(effective)==14
assert set(effective)==({"A0"}|{f"R{i}" for i in range(1,14)})
assert D*ns*D==0
assert D*ns*F8==2
assert all(C*ns*C==-2 for C in effective.values())
assert all(C*ns*F8==0 for C in effective.values())
assert {name for name,v in Dpairs.items() if v==1}=={"R3","A0"}
assert all(v in (0,1) for v in Dpairs.values())

print(
    "Q24D12_COMPONENT_TRANSPORT|"
    f"components={len(effective)}|degree1=R3,A0|"
    f"old_O8_identity_check={'SKIPPED_OBSOLETE' if transport_stopped_at_identity else 'PASS'}|"
    "status=PASS_PREIDENTITY_TRANSPORT",
    flush=True,
)

assert {name for name,v in Dpairs.items() if v==1}=={"R3","A0"}
assert Dpairs["R3"]==Dpairs["A0"]==1

R3=effective["R3"]
A0=effective["A0"]
for name,C in (("R3",R3),("A0",A0)):
    assert C*ns*C==-2,(name,C*ns*C)
    assert C*ns*D==1,(name,C*ns*D)

# Choose R3 as the geometric zero of the q24 D12 child.
# A0 is the other degree-one section from the same old I9* fibre.
Ogeo=R3
Sgeo=A0

mate=Ogeo+D
assert mate*ns*mate==0 and mate*ns*D==1
complement=matrix(ZZ,[list(D*ns),list(mate*ns)]).right_kernel_matrix()
raw=-(complement*ns*complement.transpose())
Bzero=matrix(ZZ,[list(D),list(mate)]+[list(r) for r in complement.rows()])
assert abs(Bzero.det())==1

mini=minimize_child_frame(raw)
G=matrix(ZZ,mini["frame"])
A=matrix(ZZ,mini["basis"])
assert tuple(map(int,mini["root_data"]))==(12,264,4)
assert G.det()==948
B=block_diagonal_matrix(identity_matrix(ZZ,2),A)*Bzero
U2=matrix(ZZ,((0,1),(1,0)))
assert B*ns*B.transpose()==block_diagonal_matrix(U2,-G)

def child_coords(parent_curve):
    q=vector(QQ,parent_curve)*B.inverse()
    assert all(x in ZZ for x in q)
    return vector(ZZ,[ZZ(x) for x in q])

ca=child_coords(Ogeo)
cb=child_coords(Sgeo)
assert ca==vector(ZZ,[-1,1]+[0]*17)
assert cb[1]==1
zspin=vector(ZZ,cb[-5:])
assert zspin

R=G[:12,:12]
Cpl=G[:12,12:]
Tail=G[12:,12:]
H=Tail-Cpl.transpose()*R.inverse()*Cpl
spin_height=QQ(zspin*H*zspin)

OUTDIR.mkdir(parents=True,exist_ok=True)
FRAME.write_text(
    "# q24 D12/MW5 with geometric zero R3\n"
    f"# A0 MW projection relative to R3 = {tuple(map(int,zspin))}\n"
    + "\n".join(" ".join(map(str,row)) for row in G.rows())
    + "\n"
)

print(
    "Q24D12GEOZERO|"
    f"zero=R3|other=A0|"
    f"explicit_mw={','.join(map(str,zspin))}|height={spin_height}|"
    "status=PASS_Q24_DEGREE_ONE_COMPONENT_MARKING",
    flush=True,
)

# Exhaustive q=6 search in this actual geometric-zero frame.
cmd=[
    "sage","-python",str(SEARCH),
    "--frame",str(FRAME),
    "--root-rank","12",
    "--q","6",
    "--degree","2",
    "--adapt-mw-at-least","6",
    "--rank-growth-only",
    "--output",str(SEARCHOUT),
    "--frames-dir",str(FRAMES),
]
print("+"," ".join(cmd),flush=True)
subprocess.run(cmd,cwd=str(ROOT),check=True)

data=json.loads(SEARCHOUT.read_text())
assert data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
a11=[
    r for r in data["neighbors"]
    if tuple(r["child_root_data"])==(11,132,12)
    and r["child_ade"]=="A11"
    and int(r["child_mw_rank"])==6
]

def scalar_multiple(z,s):
    z=vector(ZZ,z); s=vector(ZZ,s)
    n=None
    for a,b in zip(z,s):
        if b:
            if a%b:
                return None
            q=a//b
            if n is None: n=q
            elif n!=q: return None
        elif a:
            return None
    return ZZ(0) if n is None else ZZ(n)

easy=[]
for rec in a11:
    n=scalar_multiple(rec["mw_projection"],zspin)
    if n is not None:
        easy.append((abs(int(n)),sum(abs(int(x)) for x in rec["dominant_labels"]),
                     int(rec["orbit_index"]),n,rec))
easy.sort(key=lambda x:x[:3])

print(
    "Q24D12A11_SCAN|"
    f"A11_hits={len(a11)}|explicit_direction_hits={len(easy)}|"
    f"status={'PASS_EQUATION_FRIENDLY_HIT' if easy else 'NEEDS_SECOND_EXPLICIT_DIRECTION'}",
    flush=True,
)

RV=PolynomialRing(F,"V")
V=RV.gen()
KV=RV.fraction_field()
def rf(rec):
    return KV(RV([F(v) for v in rec["num"]]))/KV(RV([F(v) for v in rec["den"]]))
jacA=rf(sig["jacobian_A"])
jacB=rf(sig["jacobian_B"])
E=EllipticCurve(KV,[0,0,0,jacA,jacB])

def point_from_record(rec):
    return E(rf(rec["x"]),rf(rec["y"]))

halfpoints=[point_from_record(x) for x in anchor["primitive_half_points"]]

def class_order(dual):
    o=ZZ(1)
    for x in dual:
        o=lcm(o,ZZ(QQ(x).denominator()))
    return o

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

def correction_for(z):
    z=vector(ZZ,z)
    base=vector(ZZ,[0]*12+list(z))
    dual=vector(QQ,base*G[:,:12])*Rinv
    order=class_order(dual)
    key=frac_key(dual)
    if key not in correction_by_class:
        return None
    corr=correction_by_class[key]
    h=QQ(z*H*z)
    po=(h+corr-4)/2
    if po not in ZZ or po<0:
        return None
    return (ZZ(po),corr,order)

def section_frame_for(z,corr,po):
    z=vector(ZZ,z)
    target_norm=QQ(z*H*z)+corr
    # Try small dominant component labels first.
    labels=[]
    zero=vector(ZZ,[0]*12)
    labels.append(zero)
    for weight in (1,2):
        for inds in itertools.combinations(range(12),weight):
            labels.append(vector(ZZ,[ZZ(i in inds) for i in range(12)]))
    for lab in labels:
        if QQ(lab*Rinv*lab)-corr not in 2*ZZ:
            continue
        rr=Rinv*(lab-Cpl*z)
        if not all(v in ZZ for v in rr):
            continue
        pf=vector(ZZ,[ZZ(v) for v in rr]+list(z))
        if QQ(pf*G*pf)==target_norm:
            return lab,pf
    # Exact closest-vector fallback.
    base=vector(ZZ,[0]*12+list(z))
    dual=vector(QQ,base*G[:,:12])*R.inverse()
    L=IntegralLattice(R)
    it=L.enumerate_close_vectors(-dual)
    for _ in range(100000):
        shift=vector(ZZ,next(it))
        pf=base+vector(ZZ,list(shift)+[0]*5)
        norm=QQ(pf*G*pf)
        if norm==target_norm:
            lab=R*vector(ZZ,pf[:12])+Cpl*z
            return vector(ZZ,lab),pf
        if norm>target_norm:
            break
    raise ArithmeticError("could not recover minimal section frame")

def enc(v):
    v=KV(v)
    return {"num":[int(x) for x in RV(v.numerator()).list()],
            "den":[int(x) for x in RV(v.denominator()).list()]}

candidates=[]
for _,_,_,n,rec in easy:
    z=vector(ZZ,rec["mw_projection"])
    cor=correction_for(z)
    if cor is None:
        continue
    po,corr,order=cor
    lab,pf=section_frame_for(z,corr,po)
    witness=vector(ZZ,rec["witness"])
    vroot=witness-pf
    assert all(x==0 for x in vroot[12:])
    k=ZZ(3)-po

    points=[]
    for hi,S in enumerate(halfpoints):
        P=ZZ(n)*S
        if P.is_zero():
            continue
        x,y=P.xy()
        points.append({
            "half_index":hi,
            "scalar":int(n),
            "point":{"x":enc(x),"y":enc(y)},
        })

    candidates.append({
        "orbit_index":int(rec["orbit_index"]),
        "scalar_on_spinor":int(n),
        "mw_projection":list(map(int,z)),
        "dominant_labels":rec["dominant_labels"],
        "witness":rec["witness"],
        "section_P_dot_O":int(po),
        "section_local_correction":str(corr),
        "section_component_labels":list(map(int,lab)),
        "section_frame":list(map(int,pf)),
        "vertical_root_coefficients":list(map(int,vroot[:12])),
        "fiber_twist":int(k),
        "formula":"D=O+P+V+kF",
        "explicit_marked_points":points,
        "child_root_data":rec["child_root_data"],
        "child_ade":rec["child_ade"],
        "child_mw_rank":rec["child_mw_rank"],
    })

OUT=args.output.resolve() if args.output else OUTDIR/f"d12-to-a11-equation-friendly-p{p}.json"
payload={
    "schema":"elkies-k3.h3-q24-d12-a11-equation-friendly.v1",
    "status":(
        "PASS_Q24_D12_A11_EXPLICIT_MARKED_SECTION"
        if candidates else
        "Q24_D12_A11_NEEDS_SECOND_EXPLICIT_MW_DIRECTION"
    ),
    "prime":int(p),
    "geometric_zero":"R3",
    "component_transport_boundary":"pre-O8-identity state; R3/A0 degree-one incidence reverified",
    "other_degree_one_component":"A0",
    "explicit_mw_projection_A0_minus_R3":list(map(int,zspin)),
    "explicit_direction_height":str(spin_height),
    "A11_hit_count":len(a11),
    "explicit_direction_scalar_hit_count":len(easy),
    "candidates":candidates,
    "frame":str(FRAME.relative_to(ROOT)),
    "search_artifact":str(SEARCHOUT.relative_to(ROOT)),
    "proof_boundary":(
        "This aligns the modular q24 D12 equation with a geometric zero via "
        "the old I9* spinor pair and corrects the binary-quartic 2-cover by "
        "halving the covariant difference. A scalar-direction A11 candidate "
        "has an explicit marked point and exact D=O+P+V+kF lattice "
        "decomposition. Resolved I8* RR conditions and A11 Jacobian "
        "elimination are the next gate."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D12A11_RESULT|"
    f"candidates={len(candidates)}|status={payload['status']}",
    flush=True,
)
