#!/usr/bin/env sage -python
"""
Identify the equation-frame correction section R in

    Q_q24 = AJ(S3) + 2 R,

where the exhaustive structural D13 isometry audit proved

    Q_q24 = +/-(-2,1,-1,1)
    AJ(S3) = (0,-1,-1,1),

so for the displayed sign

    R = (-1,1,0,0).

Recover R as an actual (-2)-curve in the current equation q8 frame and test
whether it is an old q6 fibre component or a small q6 MW section.
"""

import json
import sys
from itertools import product
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix,
    lcm, matrix, pari, vector
)


def locate_repo():
    cwd=Path.cwd().resolve()
    candidates=[cwd,*cwd.parents]
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
    raise SystemExit("Could not locate jacobian-research")


ROOT=locate_repo()
EQ=ROOT/"elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUT=LOCAL/"q8-q24-correction-section-R.json"


def run(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


eq=run(EQ)
ns=eq["ns"]
F8=vector(ZZ,eq["F8eq"])
O8=vector(ZZ,eq["selected"]["E8"]["simple_root_vectors_in_source_h3_ns"][0])
F6=vector(ZZ,eq["F6"])
O6=vector(ZZ,eq["Oold"])
Fh3=vector(ZZ,eq["Fh3"])
roots_and_data=eq["roots_and_data"]

assert O8*ns*F8==1
assert O6*ns*F6==1


def child_frame_with_zero(ns,fibre,zero):
    mate=zero+fibre
    assert mate*ns*mate==0 and mate*ns*fibre==1
    orth=matrix(ZZ,[list(fibre*ns),list(mate*ns)]).right_kernel_matrix()
    B=matrix(ZZ,[list(fibre),list(mate)]+[list(r) for r in orth.rows()])
    assert abs(B.det())==1
    child=-(orth*ns*orth.transpose())
    return child,B


def deterministic_simple_roots(gram):
    roots,unused,data=roots_and_data(gram)
    rank=data[0]
    regular=None
    for shift in range(1,1000):
        c=vector(ZZ,[(i+1)**2+shift*(i+1)+1 for i in range(gram.nrows())])
        if all(c*r!=0 for r in roots):
            regular=c
            break
    assert regular is not None
    positive=[r for r in roots if regular*r>0]
    pset={tuple(r) for r in positive}
    simple=[r for r in positive if not any(tuple(r-left) in pset for left in positive)]
    simple=matrix(ZZ,[list(r) for r in simple])
    assert simple.nrows()==simple.rank()==rank
    return simple,simple*gram*simple.transpose()


def d13_root_adaptation(child):
    unused,rb,inv=roots_and_data(child)
    assert inv==(13,312,4)
    simple,cartan=deterministic_simple_roots(child)
    smith,left,right=rb.smith_form()
    assert tuple(abs(smith[i,i]) for i in range(13))==(1,)*13
    completion=right.inverse()
    initial=simple.stack(completion[13:])
    assert abs(initial.det())==1
    G=initial*child*initial.transpose()
    C=G[:13,13:]
    T=G[13:,13:]
    H=T-C.transpose()*cartan.inverse()*C
    scale=ZZ(1)
    for x in H.list():
        scale=lcm(scale,ZZ(QQ(x).denominator()))
    lll=matrix(ZZ,pari((scale*H).change_ring(ZZ)).qflllgram())
    change=block_diagonal_matrix(identity_matrix(ZZ,13),lll.transpose())
    A=change*initial
    G=A*child*A.transpose()
    R=G[:13,:13]
    C=G[:13,13:]
    T=G[13:,13:]
    H=T-C.transpose()*R.inverse()*C
    assert H.det()==237
    return A,G,H


def class_order(dual):
    o=ZZ(1)
    for x in dual:
        o=lcm(o,ZZ(QQ(x).denominator()))
    return o


def correction_for_dual(dual,root):
    o=class_order(dual)
    expected={ZZ(1):QQ(0),ZZ(2):QQ(1),ZZ(4):QQ(13)/4}
    assert o in expected
    return o,expected[o]


def section_for_mw(G,z,cap=100000):
    R=G[:13,:13]
    C=G[:13,13:]
    T=G[13:,13:]
    H=T-C.transpose()*R.inverse()*C
    z=vector(ZZ,z)
    h=QQ(z*H*z)
    base=vector(ZZ,[0]*13+list(z))
    pairing=vector(QQ,base*G[:,:13])
    dual=pairing*R.inverse()
    order,corr=correction_for_dual(dual,R)
    target=h+corr
    assert target in ZZ and target>=4 and target%2==0
    it=IntegralLattice(R).enumerate_close_vectors(-dual)
    chosen=None
    for _ in range(cap):
        shift=vector(ZZ,next(it))
        cand=base+vector(ZZ,list(shift)+[0]*4)
        norm=ZZ(cand*G*cand)
        if norm==target:
            chosen=cand
            break
        if norm>target:
            break
    assert chosen is not None
    pole=ZZ((target-4)/2)
    a=ZZ((target-2)/2)
    return vector(ZZ,[a,1]+list(chosen)),h,order,corr,pole


child8,Bzero8=child_frame_with_zero(ns,F8,O8)
A8,G8,H8=d13_root_adaptation(child8)
Badapt8=block_diagonal_matrix(identity_matrix(ZZ,2),A8)*Bzero8


def q8coords(C):
    c=vector(QQ,C)*Badapt8.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ,c)


A=vector(ZZ,(0,-1,-1,1))
Q=vector(ZZ,(-2,1,-1,1))
assert all((Q[i]-A[i])%2==0 for i in range(4))
Rmw=vector(ZZ,[(Q[i]-A[i])//2 for i in range(4)])
assert Rmw==vector(ZZ,(-1,1,0,0))

Pcoords,hR,orderR,corrR,poleR=section_for_mw(G8,Rmw)
PR=Pcoords*Badapt8
assert PR*ns*PR==-2 and PR*ns*F8==1
assert vector(ZZ,q8coords(PR)[-4:])==Rmw

hA=QQ(A*H8*A)
hQ=QQ(Q*H8*Q)
pairAR=QQ(A*H8*Rmw)

print(
    "Q8Q24R_PROFILE|"
    f"R={','.join(map(str,Rmw))}|height={hR}|class_order={orderR}|corr={corrR}|"
    f"PdotO8={poleR}|q6_degree={PR*ns*F6}|q6_Oold={PR*ns*O6}|H3_degree={PR*ns*Fh3}|"
    f"hAJ={hA}|hQ={hQ}|AJdotR={pairAR}|"
    "status=PASS",
    flush=True,
)

print(
    "Q8Q24R_HEIGHT_GRAM|rows="
    +";".join(",".join(str(x) for x in row) for row in H8.rows()),
    flush=True,
)

# ---------------------------------------------------------------------------
# Reconstruct actual q6 IV*/II* components in THIS equation ambient marking.
# ---------------------------------------------------------------------------

E6_CARTAN=matrix(ZZ,[
    [2,-1,0,0,0,0],
    [-1,2,-1,0,0,0],
    [0,-1,2,-1,0,-1],
    [0,0,-1,2,-1,0],
    [0,0,0,-1,2,0],
    [0,0,-1,0,0,2],
])
E8_CARTAN=matrix(ZZ,[
    [2,0,-1,0,0,0,0,0],
    [0,2,0,-1,0,0,0,0],
    [-1,0,2,-1,0,0,0,0],
    [0,-1,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,0],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,-1],
    [0,0,0,0,0,0,-1,2],
])
E6_QF_INDICES=(0,1,2,3,12,13)
E6_SIMPLE_IN_QF=(
    (-1,-1,-1,-1,0,0),
    (0,0,0,1,0,0),
    (0,0,1,0,0,0),
    (0,1,0,0,0,0),
    (1,0,0,0,0,1),
    (2,1,0,0,-1,1),
)


def highest_root(C):
    half=matrix(ZZ,pari(C).qfminim(2)[2]).transpose().rows()
    rr=[vector(ZZ,r) for r in half]
    rr += [-vector(ZZ,r) for r in half]
    return max((r for r in rr if all(x>=0 for x in r)),key=lambda r:sum(r))


q6orth=matrix(
    ZZ,[list(F6*ns),list((O6+F6)*ns)]
).right_kernel_matrix()
q6child=-(q6orth*ns*q6orth.transpose())
qfbasis=matrix(ZZ,pari(q6child).qfminim(2)[2]).transpose().row_module().basis_matrix()
assert qfbasis.rank()==14

e6qf=matrix(ZZ,[list(qfbasis[i]) for i in E6_QF_INDICES])
e6=matrix(ZZ,[list(vector(ZZ,row)*e6qf) for row in E6_SIMPLE_IN_QF])*q6orth
allphysical=qfbasis*q6orth
e8=allphysical[4:12,:]

assert -e6*ns*e6.transpose()==E6_CARTAN
assert -e8*ns*e8.transpose()==E8_CARTAN

e6aff=F6-highest_root(E6_CARTAN)*e6
e8aff=F6-highest_root(E8_CARTAN)*e8

components=[]
for i,row in enumerate(e6.rows(),1):
    components.append((f"IV*_E6_{i}",vector(ZZ,row)))
components.append(("IV*_affine",vector(ZZ,e6aff)))
for i,row in enumerate(e8.rows(),1):
    components.append((f"II*_E8_{i}",vector(ZZ,row)))
components.append(("II*_affine",vector(ZZ,e8aff)))

component_hits=[]
section_components=[]
for name,C in components:
    assert C*ns*C==-2 and C*ns*F6==0
    d=ZZ(C*ns*F8)
    if d==1:
        z=vector(ZZ,q8coords(C)[-4:])
        section_components.append((name,z,C))
        if C==PR or z==Rmw:
            component_hits.append((name,z,C==PR))

print(
    "Q8Q24R_Q6_COMPONENTS|q8_sections={}|hits={}|status=PASS".format(
        len(section_components),
        ";".join(
            f"{name}:{','.join(map(str,z))}:exact={int(exact)}"
            for name,z,exact in component_hits
        ) or "none"
    ),
    flush=True,
)

for name,z,C in section_components:
    print(
        "Q8Q24R_COMPONENT|"
        f"name={name}|mw={','.join(map(str,z))}|"
        f"exact_R={int(C==PR)}|status=CANDIDATE",
        flush=True,
    )

# ---------------------------------------------------------------------------
# Search small q6 MW words handled by the certified identity-component
# reconstruction. This is finite and cheap.
# ---------------------------------------------------------------------------

mw_hits=[]
for word in product(range(-5,6),repeat=3):
    if word==(0,0,0):
        continue
    if (word[0]-word[1])%3:
        continue
    try:
        C,pole,h=eq["section_from_old_mw"](vector(ZZ,word))
    except Exception:
        continue
    if C*ns*F8 != 1:
        continue
    z=vector(ZZ,q8coords(C)[-4:])
    if z==Rmw or C==PR:
        mw_hits.append((word,z,C==PR,pole,h))

print(
    "Q8Q24R_Q6_MW_SEARCH|"
    f"box=5|hits={len(mw_hits)}|"
    + (
        "records="+ ";".join(
            f"{','.join(map(str,w))}:{','.join(map(str,z))}:exact={int(exact)}:pole={pole}:h={h}"
            for w,z,exact,pole,h in mw_hits
        )
        if mw_hits else "records=none"
    )
    +"|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h92-q8-q24-correction-section-R.v1",
    "status":"PASS_Q24_CORRECTION_SECTION_PROFILE",
    "relation":{
        "AJ_S3_mw":list(map(int,A)),
        "q24_mw_up_to_sign":list(map(int,Q)),
        "R_mw":list(map(int,Rmw)),
        "formula":"Q_q24 = AJ(S3) + 2 R (up to global inversion)",
    },
    "R":{
        "height":str(hR),
        "class_order":int(orderR),
        "correction":str(corrR),
        "P_dot_O8":int(poleR),
        "q6_degree":int(PR*ns*F6),
        "q6_Oold_intersection":int(PR*ns*O6),
        "H3_degree":int(PR*ns*Fh3),
        "section_source_h3_ns":list(map(int,PR)),
    },
    "height_gram":[[str(x) for x in row] for row in H8.rows()],
    "q6_component_hits":[
        {"name":name,"mw":list(map(int,z)),"exact_curve":bool(exact)}
        for name,z,exact in component_hits
    ],
    "q6_mw_hits":[
        {
            "old_q6_mw":list(map(int,w)),
            "q8_mw":list(map(int,z)),
            "exact_curve":bool(exact),
            "q6_P_dot_O":int(pole),
            "q6_height":str(h),
        }
        for w,z,exact,pole,h in mw_hits
    ],
}

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q8Q24R_RESULT|"
    f"R={','.join(map(str,Rmw))}|height={hR}|corr={corrR}|PdotO8={poleR}|"
    f"component_hits={len(component_hits)}|q6mw_hits={len(mw_hits)}|"
    "status=PASS_Q24_CORRECTION_SECTION_IDENTIFICATION",
    flush=True,
)
