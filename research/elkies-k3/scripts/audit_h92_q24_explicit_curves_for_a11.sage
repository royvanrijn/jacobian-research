#!/usr/bin/env sage -python
"""
Audit already-explicit curves against the q24 fibration.

Goal: find a second explicit D12/MW5 direction without any section search.
All curve classes live in the same H3 NS frame.  We:
  * recover D24eq and the exact transported degree-one component R3;
  * use R3 as the q24/D12 zero and build its root-adapted child frame;
  * collect explicit (-2)-curves already present in the q8 equation certifier;
  * keep those with D24eq.C = 1;
  * compute their MW projections in the R3-zero D12 frame;
  * test the four q6 -> A11 targets against their integral span.
"""
import argparse
import contextlib
import io
import itertools
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

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
CERT=SCRIPTS/"certify_h92_q8_equation_ns_divisor.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"
SEARCHOUT=OUTDIR/"d12-c10a-zero-q6-all.json"

for path in (ENGINE,EFF,CERT,CLOSE,SEARCHOUT):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

def run_scope(path,argv=(),allow_late_assert=False):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__","__file__":str(path)}
    buf=io.StringIO()
    late=False
    try:
        sys.argv=[str(path)]+list(argv)
        with contextlib.redirect_stdout(buf):
            try:
                exec(compile(path.read_text(),str(path),"exec"),scope)
            except AssertionError:
                if not allow_late_assert:
                    raise
                late=True
    finally:
        sys.argv=saved
    return scope,buf.getvalue(),late

# q24 effective component transport: accept only its structurally reverified
# pre-O8-identity state.
eff,effout,late=run_scope(EFF,("--prime",str(p)),allow_late_assert=True)
need=("ns","D","F8","effective","Dpairs")
missing=[k for k in need if k not in eff]
if missing:
    raise SystemExit("effective transport stopped too early: "+",".join(missing))

ns=matrix(ZZ,eff["ns"])
D=vector(ZZ,eff["D"])
F8=vector(ZZ,eff["F8"])
effective={name:vector(ZZ,c) for name,c in eff["effective"].items()}
Dpairs={name:int(v) for name,v in eff["Dpairs"].items()}

assert D*ns*D==0 and D*ns*F8==2
assert len(effective)==14
assert {name for name,v in Dpairs.items() if v==1}=={"R3","A0"}
assert all(C*ns*C==-2 and C*ns*F8==0 for C in effective.values())

R3=effective["R3"]
A0=effective["A0"]
assert R3*ns*D==A0*ns*D==1

# Build q24 D12 with R3 as zero.
mate=R3+D
assert mate*ns*mate==0 and mate*ns*D==1
orth=matrix(ZZ,[list(D*ns),list(mate*ns)]).right_kernel_matrix()
raw=-(orth*ns*orth.transpose())
Bzero=matrix(ZZ,[list(D),list(mate)]+[list(r) for r in orth.rows()])
assert abs(Bzero.det())==1
mini=minimize_child_frame(raw)
G=matrix(ZZ,mini["frame"])
A=matrix(ZZ,mini["basis"])
assert tuple(map(int,mini["root_data"]))==(12,264,4)
B=block_diagonal_matrix(identity_matrix(ZZ,2),A)*Bzero
U2=matrix(ZZ,((0,1),(1,0)))
assert B*ns*B.transpose()==block_diagonal_matrix(U2,-G)

def child_coords(C):
    q=vector(QQ,C)*B.inverse()
    if not all(x in ZZ for x in q):
        raise ArithmeticError("nonintegral child coordinates")
    return vector(ZZ,[ZZ(x) for x in q])

def mw_projection(C):
    q=child_coords(C)
    assert q[1]==1
    return vector(ZZ,q[-5:]),q

zA0,_=mw_projection(A0)
assert zA0

# Load q8 equation certifier: it exposes many exact curve classes in same NS.
cert,certout,_=run_scope(CERT)
assert matrix(ZZ,cert["ns"])==ns

curves={}
def add(name,C,origin):
    C=vector(ZZ,C)
    if len(C)!=19:
        return
    curves[name]={"class":C,"origin":origin}

for name in ("Oold","Ostd","Smarked","S3std","F6","F8eq"):
    if name in cert:
        add(name,cert[name],"q8_equation_certifier")

for i,C in enumerate(cert.get("basis_sections",()),1):
    add(f"q6_basis_{i}",C,"q8_equation_certifier:basis_sections")

for i,C in enumerate(cert.get("simple",()),1):
    add(f"H3_simple_{i}",C,"H3_source_simple_roots")

# Current q8 component-nef root vectors are also exact curve classes.
selected=cert.get("selected")
if selected:
    for family in ("E6","E8"):
        for i,row in enumerate(selected[family]["simple_root_vectors_in_source_h3_ns"],1):
            add(f"q8_{family}_{i}",row,"q8_component_nef_roots")

# q24 data exposed by the close replay can contribute P24/O8.
close,closeout,_=run_scope(CLOSE)
assert matrix(ZZ,close["ns"])==ns
for name in ("O8","P24","S3std","Smarked","Ostd","Oold"):
    if name in close:
        add(f"close_{name}",close[name],"q24_close_replay")

# Actual affine-D13 components, including A0 and R3.
for name,C in effective.items():
    add(f"oldI9_{name}",C,"actual_effective_I9star")

records=[]
degree1=[]
seen_classes={}
for name,item in curves.items():
    C=item["class"]
    key=tuple(map(int,C))
    if key in seen_classes:
        continue
    seen_classes[key]=name
    square=int(C*ns*C)
    degree=int(D*ns*C)
    rec={
        "name":name,
        "origin":item["origin"],
        "square":square,
        "q24_degree":degree,
        "class":list(map(int,C)),
    }
    if square==-2 and degree==1:
        z,q=mw_projection(C)
        rec["d12_mw_projection"]=list(map(int,z))
        rec["d12_full_coordinates"]=list(map(int,q))
        degree1.append((name,C,z,q,item["origin"]))
    records.append(rec)

print(
    "Q24EXPLICIT_AUDIT|"
    f"curves={len(records)}|degree1={len(degree1)}|"
    f"names={','.join(name for name,_,_,_,_ in degree1)}|status=PASS",
    flush=True,
)

# Deduplicate MW directions up to sign and zero.
direction_map={}
for name,C,z,q,origin in degree1:
    if z==0:
        continue
    canonical=min(tuple(map(int,z)),tuple(map(int,-z)))
    direction_map.setdefault(canonical,[]).append(name)

directions=[vector(ZZ,key) for key in sorted(direction_map)]
if directions:
    M=matrix(ZZ,[list(z) for z in directions])
    L=M.row_module()
    basis=L.basis_matrix()
    rank=L.rank()
else:
    M=matrix(ZZ,0,5)
    L=None
    basis=matrix(ZZ,0,5)
    rank=0

print(
    "Q24EXPLICIT_MW|"
    f"directions={len(directions)}|rank={rank}|"
    f"A0_direction={','.join(map(str,zA0))}|status=PASS",
    flush=True,
)

search=json.loads(SEARCHOUT.read_text())
assert search["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
a11=[
    r for r in search["neighbors"]
    if tuple(r["child_root_data"])==(11,132,12)
    and r["child_ade"]=="A11"
    and int(r["child_mw_rank"])==6
]
assert len(a11)==4

matches=[]
for rec in a11:
    target=vector(ZZ,rec["mw_projection"])
    in_span=bool(L is not None and target in L)
    item={
        "orbit_index":int(rec["orbit_index"]),
        "mw_projection":list(map(int,target)),
        "dominant_labels":rec["dominant_labels"],
        "in_explicit_integral_span":in_span,
    }
    if in_span:
        coeff=basis.transpose().solve_right(target)
        assert all(c in ZZ for c in coeff)
        coeff=vector(ZZ,[ZZ(c) for c in coeff])
        assert coeff*basis==target
        item["explicit_span_basis"]=[list(map(int,row)) for row in basis.rows()]
        item["coefficients_in_span_basis"]=list(map(int,coeff))
        matches.append(item)
    print(
        "Q24EXPLICIT_A11|"
        f"orbit={rec['orbit_index']}|mw={','.join(map(str,target))}|"
        f"in_span={int(in_span)}|status={'PASS_MATCH' if in_span else 'NO_MATCH'}",
        flush=True,
    )

payload={
    "schema":"elkies-k3.h3-q24-explicit-curve-a11-span.v1",
    "status":(
        "PASS_Q24_A11_TARGET_IN_EXPLICIT_CURVE_SPAN"
        if matches else
        "Q24_A11_TARGET_NOT_IN_CURRENT_EXPLICIT_CURVE_SPAN"
    ),
    "prime":int(p),
    "geometric_zero":"R3",
    "A0_mw_projection":list(map(int,zA0)),
    "explicit_curve_records":records,
    "degree_one_sections":[
        {
            "name":name,
            "origin":origin,
            "mw_projection":list(map(int,z)),
            "full_child_coordinates":list(map(int,q)),
        }
        for name,C,z,q,origin in degree1
    ],
    "distinct_mw_directions":[
        {
            "canonical_mw":list(key),
            "curve_names":direction_map[key],
        }
        for key in sorted(direction_map)
    ],
    "explicit_span_rank":int(rank),
    "explicit_span_basis":[list(map(int,row)) for row in basis.rows()],
    "A11_targets":[
        {
            "orbit_index":int(rec["orbit_index"]),
            "mw_projection":rec["mw_projection"],
            "dominant_labels":rec["dominant_labels"],
        }
        for rec in a11
    ],
    "matching_A11_targets":matches,
    "boundary":(
        "This is purely an exact NS/lattice audit of already-explicit curves. "
        "It does not yet produce q24-quartic coordinates for a newly found "
        "degree-one curve. If an A11 target lies in the explicit span, the next "
        "step is to transport the corresponding explicit curve(s) through the "
        "q24 pencil/covariant and build the marked point."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/f"explicit-curves-a11-span-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24EXPLICIT_RESULT|"
    f"degree1={len(degree1)}|span_rank={rank}|"
    f"A11_matches={len(matches)}|status={payload['status']}",
    flush=True,
)
