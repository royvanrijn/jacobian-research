#!/usr/bin/env sage -python
"""
Profile the R17-directed q24 D12 -> A11 orbit42 divisor under the two explicit
geometric q24 zeros R3 and A0.

Inputs:
  * exact current-equation-D13 -> historical D12 -> pinned R17 certificate;
  * exact orbit42 current-equation bridge;
  * actual effective old-I9* component transport.

For each zero:
  * rebuild the D12 child basis with that geometric section as zero;
  * transport the exact historical orbit42 A11 fibre into that D12 frame;
  * compute the exact D12 discriminant correction and P.O;
  * recover a minimal section representative and D=O+P+V+kF;
  * identify the matching record in the previous exhaustive q6 search when
    available;
  * test whether P lies in the explicit opposite-spinor MW direction.

No section enumeration or Groebner basis is used.
"""
import argparse, contextlib, io, json, sys
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector
)

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
EFF=SCRIPTS/"transport_h92_q24_effective_i9star_components.sage"
BRIDGE=LOCAL/"q24-orbit42-current-equation-bridge.json"
SEARCHOUT=OUTDIR/"d12-c10a-zero-q6-all.json"

for path in (ENGINE,EFF,BRIDGE):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

bridge=json.loads(BRIDGE.read_text())
assert bridge["status"]=="PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE"

# Replay effective component transport only up through its valid pre-identity state.
saved=list(sys.argv)
scope={"__name__":"__embedded_q24_effcomp__","__file__":str(EFF)}
buf=io.StringIO()
stopped=False
try:
    sys.argv=[str(EFF),"--prime",str(p)]
    with contextlib.redirect_stdout(buf):
        try:
            exec(compile(EFF.read_text(),str(EFF),"exec"),scope)
        except AssertionError:
            stopped=True
finally:
    sys.argv=saved

need=("ns","D","effective","Dpairs")
missing=[x for x in need if x not in scope]
if missing:
    raise SystemExit("effective transport missing: "+",".join(missing))

U2=matrix(ZZ,((0,1),(1,0)))

ns=matrix(ZZ,scope["ns"])
F24=vector(ZZ,scope["D"])
effective={k:vector(ZZ,v) for k,v in scope["effective"].items()}
Dpairs={k:int(v) for k,v in scope["Dpairs"].items()}

assert {name for name,v in Dpairs.items() if v==1}=={"R3","A0"}
R3=effective["R3"]
A0=effective["A0"]
assert R3*ns*F24==A0*ns*F24==1
assert R3*ns*R3==A0*ns*A0==-2

# Coordinate bridge:
# q24-orbit42-current-equation-bridge.json stores classes in the canonical
# equation-D13 BASIS coordinates.  The effective component transport stores
# R3/A0/D24eq in ambient/source-NS coordinates.  The embedded EFF script has
# already run close_h92_q8_q24_by_q6_translation.sage and retained its `cl`
# scope; Badapt has rows equal to the equation-D13 basis in ambient coords.
if "cl" not in scope or "Badapt" not in scope["cl"]:
    raise SystemExit("effective transport scope does not expose close/Badapt")
Badapt=matrix(ZZ,scope["cl"]["Badapt"])
adapted=matrix(ZZ,scope["cl"]["adapted"])
assert Badapt*ns*Badapt.transpose()==block_diagonal_matrix(U2,-adapted)

# Hard coordinate-system regression: the q24 fibre from the R17 bridge must
# become exactly the effective transport's q24 divisor in ambient coordinates.
F24eqcoords=vector(ZZ,bridge["current_equation_D13"]["q24_fibre"])
F24frombridge=vector(ZZ,F24eqcoords*Badapt)
assert F24frombridge==F24, (F24frombridge,F24)

# Now map the exact historical/R17-directed orbit42 A11 fibre into the same
# ambient coordinates as R3/A0.
D42coords=vector(ZZ,bridge["current_equation_D13"]["A11_fibre"])
D42eq=vector(ZZ,D42coords*Badapt)

print(
    "Q24O42SPINOR_COORDS|"
    f"q24_bridge_match={int(F24frombridge==F24)}|"
    f"D42_square={D42eq*ns*D42eq}|D42_degree={D42eq*ns*F24}|"
    "status=PASS",
    flush=True,
)

assert D42eq*ns*D42eq==0
assert D42eq*ns*F24==2



def scalar_multiple(z,s):
    z=vector(ZZ,z); s=vector(ZZ,s)
    n=None
    for a,b in zip(z,s):
        if b:
            if a%b:
                return None
            q=a//b
            if n is None:
                n=q
            elif n!=q:
                return None
        elif a:
            return None
    return ZZ(0) if n is None else ZZ(n)

def exact_profile(G,Dcoords):
    root_rank=12
    R=G[:root_rank,:root_rank]
    C=G[:root_rank,root_rank:]
    Tail=G[root_rank:,root_rank:]
    H=Tail-C.transpose()*R.inverse()*C

    w=vector(ZZ,Dcoords[2:])
    z=vector(ZZ,w[-5:])
    h=QQ(z*H*z)

    base=vector(ZZ,[0]*12+list(z))
    pair=vector(ZZ,base*G[:,:12])
    dual=vector(QQ,pair)*R.inverse()

    # Find the exact minimal dual norm in this discriminant coset.
    if all(x in ZZ for x in dual):
        corr=QQ(0)
        lam=pair
        r=vector(ZZ,[0]*12)
    else:
        detR=ZZ(R.det())
        Radj=(detR*R.inverse()).change_ring(ZZ)
        # D12 corrections are <=3, so norm <= 3*detR.
        qf=pari(Radj).qfminim(ZZ(3)*detR)
        found=[]
        M=matrix(ZZ,qf[2])
        seen=set()
        for col in M.columns():
            for sg in (1,-1):
                ll=sg*vector(ZZ,col)
                key=tuple(map(int,ll))
                if key in seen:
                    continue
                seen.add(key)
                rQ=vector(QQ,ll-pair)*R.inverse()
                if not all(x in ZZ for x in rQ):
                    continue
                rr=vector(ZZ,[ZZ(x) for x in rQ])
                ss=vector(QQ,rr)+dual
                cc=QQ(ss*R*ss)
                found.append((cc,sum(abs(int(x)) for x in rr),key,ll,rr))
        if not found:
            raise ArithmeticError(f"no dual representative for {tuple(z)}")
        found.sort(key=lambda x:(x[0],x[1],x[2]))
        corr,_,_,lam,r=found[0]

    po=(h+corr-4)/2
    if po not in ZZ or po<0:
        raise ArithmeticError((tuple(z),h,corr,po))
    po=ZZ(po)

    pframe=vector(ZZ,list(r)+list(z))
    assert QQ(pframe*G*pframe)==h+corr
    P=vector(ZZ,[po+1,1]+list(pframe))
    O=vector(ZZ,[-1,1]+[0]*17)
    F=vector(ZZ,[1,0]+[0]*17)

    Gns=block_diagonal_matrix(U2,-G)
    assert P*Gns*P==-2
    assert P*Gns*F==1
    assert P*Gns*O==po

    residual=vector(ZZ,Dcoords)-O-P
    assert residual[1]==0
    assert all(x==0 for x in residual[14:])
    k=ZZ(residual[0])
    vr=vector(ZZ,residual[2:14])

    return {
        "mw":z,
        "height":h,
        "correction":corr,
        "PdotO":po,
        "section":P,
        "section_frame":pframe,
        "dual_pairing":vector(ZZ,lam),
        "fiber_twist":k,
        "vertical_root":vr,
        "vertical_L1":sum(abs(int(x)) for x in vr),
        "vertical_support":sum(bool(x) for x in vr),
    }

search=None
if SEARCHOUT.exists():
    search=json.loads(SEARCHOUT.read_text())

rows=[]
for zero_name,Ogeo,other_name,Sgeo in (
    ("R3",R3,"A0",A0),
    ("A0",A0,"R3",R3),
):
    mate=Ogeo+F24
    assert mate*ns*mate==0 and mate*ns*F24==1
    complement=matrix(ZZ,[list(F24*ns),list(mate*ns)]).right_kernel_matrix()
    raw=-(complement*ns*complement.transpose())
    Bzero=matrix(ZZ,[list(F24),list(mate)]+[list(r) for r in complement.rows()])
    assert abs(Bzero.det())==1

    mini=minimize_child_frame(raw)
    G=matrix(ZZ,mini["frame"])
    A=matrix(ZZ,mini["basis"])
    assert tuple(map(int,mini["root_data"]))==(12,264,4)
    B=block_diagonal_matrix(identity_matrix(ZZ,2),A)*Bzero
    assert abs(B.det())==1
    assert B*ns*B.transpose()==block_diagonal_matrix(U2,-G)

    Binv=B.inverse()
    assert Binv.change_ring(ZZ)==Binv
    Binv=Binv.change_ring(ZZ)

    def child(v):
        c=vector(ZZ,v)*Binv
        assert all(x in ZZ for x in c)
        return vector(ZZ,c)

    czero=child(Ogeo)
    cother=child(Sgeo)
    cD42=child(D42eq)

    assert czero==vector(ZZ,[-1,1]+[0]*17)
    assert cother[1]==1
    assert cD42[1]==2
    assert cD42[0]*2*2 - vector(ZZ,cD42[2:])*G*vector(ZZ,cD42[2:])==0

    zspin=vector(ZZ,cother[-5:])
    pr=exact_profile(G,cD42)
    scalar=scalar_multiple(pr["mw"],zspin)

    match_orbits=[]
    if search is not None and zero_name=="R3":
        for rec in search["neighbors"]:
            if tuple(rec["child_root_data"])!=(11,132,12):
                continue
            f=vector(ZZ,rec["fiber"])
            if f==cD42:
                match_orbits.append(int(rec["orbit_index"]))

    row={
        "zero":zero_name,
        "other":other_name,
        "frame":[[int(x) for x in rr] for rr in G.rows()],
        "parent_to_child_basis":[[int(x) for x in rr] for rr in B.rows()],
        "D42_child_coordinates":list(map(int,cD42)),
        "opposite_spinor_mw":list(map(int,zspin)),
        "target_mw":list(map(int,pr["mw"])),
        "height":str(pr["height"]),
        "correction":str(pr["correction"]),
        "P_dot_O":int(pr["PdotO"]),
        "fiber_twist":int(pr["fiber_twist"]),
        "vertical_L1":int(pr["vertical_L1"]),
        "vertical_support":int(pr["vertical_support"]),
        "dual_pairing":list(map(int,pr["dual_pairing"])),
        "section_class":list(map(int,pr["section"])),
        "scalar_on_opposite_spinor":None if scalar is None else int(scalar),
        "matching_R3_zero_search_orbits":match_orbits,
    }
    rows.append(row)

    print(
        "Q24O42SPINOR_PROFILE|"
        f"zero={zero_name}|other={other_name}|"
        f"mw={','.join(map(str,pr['mw']))}|height={pr['height']}|"
        f"corr={pr['correction']}|PdotO={pr['PdotO']}|"
        f"vertical_F={pr['fiber_twist']}|vertical_L1={pr['vertical_L1']}|"
        f"scalar_on_other={scalar if scalar is not None else 'NO'}|"
        f"search_orbits={','.join(map(str,match_orbits)) if match_orbits else 'none'}|"
        "status=PASS",
        flush=True,
    )

# Select the equation-cheapest explicit zero.
rows.sort(
    key=lambda r:(
        r["P_dot_O"],
        r["vertical_L1"],
        0 if r["scalar_on_opposite_spinor"] is not None else 1,
        r["zero"],
    )
)
best=rows[0]

print(
    "Q24O42SPINOR_BEST|"
    f"zero={best['zero']}|PdotO={best['P_dot_O']}|"
    f"corr={best['correction']}|height={best['height']}|"
    f"vertical_F={best['fiber_twist']}|vertical_L1={best['vertical_L1']}|"
    f"scalar_on_other={best['scalar_on_opposite_spinor'] if best['scalar_on_opposite_spinor'] is not None else 'NO'}|"
    "status=SELECTED",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-orbit42-spinor-zero-profiles.v1",
    "status":"PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES",
    "prime":int(p),
    "source_bridge":str(BRIDGE.relative_to(ROOT)),
    "profiles":rows,
    "selected_zero":best["zero"],
    "route_end":"pinned R17",
    "proof_boundary":(
        "Exact integral lattice/marking and divisor decomposition under the two "
        "explicit q24 degree-one spinor zeros. No equation-level section has "
        "yet been recovered unless scalar_on_opposite_spinor is non-null."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/"q24-orbit42-spinor-zero-profiles.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42SPINOR_RESULT|"
    f"selected_zero={best['zero']}|PdotO={best['P_dot_O']}|"
    "status=PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES",
    flush=True,
)
