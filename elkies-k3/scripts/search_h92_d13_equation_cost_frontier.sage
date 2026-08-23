#!/usr/bin/env sage -python
"""
Equation-cost frontier from the exact H3 D13/MW4 frame.

Search degree-two elliptic neighbours for even q up to --qmax, but only build
children whose horizontal MW section has P.O <= --max-pole.  This optimizes
what matters for characteristic-zero equation construction, rather than
forcing the historical D12 child.

For a section MW class z on D13:
    P.O = (height(z) + local_correction(z) - 4)/2.
Thus small P.O means a small exact section:
    x = X/Z^2, y = Y/Z^3
with deg Z=P.O, deg X<=2P.O+4, deg Y<=3P.O+6.

Report every primitive rank-growing child (MW >= 5), sorted by:
    P.O, q, -MW(child), witness L1.
"""

import argparse, json, math
from pathlib import Path
from collections import defaultdict
from sage.all import (
    GF, QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix,
    pari, vector, xgcd
)

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
DATA=ROOT/"elkies-k3/data/fibrations"
FRAME=DATA/"h3_q6_q8_d13_mw4_root_adapted_frame.txt"

parser=argparse.ArgumentParser()
parser.add_argument("--qmax",type=int,default=40)
parser.add_argument("--max-pole",type=int,default=8)
parser.add_argument("--output",type=Path)
parser.add_argument("--frames-dir",type=Path)
args=parser.parse_args()
if args.qmax<4:
    raise ValueError("qmax must be >=4")
if args.max_pole<0:
    raise ValueError("max-pole must be >=0")

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

G=load_gram(FRAME)
assert G.dimensions()==(17,17) and G.det()==948
root_rank=13
C=G[:13,:13]
coupling=G[:13,13:]
tail=G[13:,13:]
H=tail-coupling.transpose()*C.inverse()*coupling
assert H.dimensions()==(4,4) and H.det()==237

U=matrix(ZZ,((0,1),(1,0)))
NS=block_diagonal_matrix(U,-G)
DET=abs(ZZ(G.det()))

def class_order(dual):
    o=ZZ(1)
    for v in dual:
        o=lcm(o,ZZ(QQ(v).denominator()))
    return o

def d13_correction(z):
    z=vector(ZZ,z)
    base=vector(ZZ,[0]*13+list(z))
    pair=vector(QQ,base*G[:,:13])
    dual=pair*C.inverse()
    order=class_order(dual)
    corr={ZZ(1):QQ(0),ZZ(2):QQ(1),ZZ(4):QQ(13)/4}[order]
    raw=QQ(dual*C*dual)
    mod2=lambda x: QQ(x)-2*(QQ(x)/2).floor()
    assert mod2(raw)==mod2(corr)
    return corr,order

def bezout_vector_for_pairing(ns,fiber):
    current=ZZ(0)
    result=[ZZ(0)]*ns.nrows()
    for i,value in enumerate(ns*fiber):
        if not value:
            continue
        g,a,b=xgcd(current,ZZ(value))
        result=[a*x for x in result]
        result[i]+=b
        current=g
    if abs(current)!=1:
        return None
    if current==-1:
        result=[-x for x in result]
    return vector(ZZ,result)

def child_frame(fiber):
    mate=bezout_vector_for_pairing(NS,fiber)
    if mate is None:
        return None
    ms=ZZ(mate*NS*mate)
    assert ms%2==0
    mate-=(ms//2)*fiber
    ker=matrix(ZZ,[list(fiber*NS),list(mate*NS)]).right_kernel_matrix()
    child=-(ker*NS*ker.transpose())
    assert child.is_positive_definite() and abs(child.det())==DET
    basis=matrix(ZZ,[list(fiber),list(mate)]+[list(r) for r in ker.rows()])
    assert abs(basis.det())==1
    return child,basis

def roots_and_data(gram):
    qf=pari(gram).qfminim(2)
    count=ZZ(qf[0])
    if not count:
        return (),matrix(ZZ,0,gram.nrows()),(0,0,1)
    roots0=[vector(ZZ,c) for c in matrix(ZZ,qf[2]).columns()]
    roots=tuple(roots0+[-r for r in roots0])
    rb=matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg=rb*gram*rb.transpose()
    return roots,rb,(int(rb.rank()),int(count),abs(int(rg.det())))

def root_rank_count(gram):
    qf=pari(gram).qfminim(2)
    count=ZZ(qf[0])
    if not count:
        return 0,0
    return int(matrix(ZZ,qf[2]).rank()),int(count)

def deterministic_simple_roots(gram):
    roots,unused,data=roots_and_data(gram)
    rr=data[0]
    pos=[r for r in roots if next(v for v in r if v)!=0 and next(v for v in r if v)>0]
    pset={tuple(r) for r in pos}
    simple=[
        r for r in pos
        if not any(tuple(r-left) in pset for left in pos)
    ]
    M=matrix(ZZ,[list(r) for r in simple])
    assert M.nrows()==M.rank()==rr
    return M,M*gram*M.transpose()

def connected_components(cartan):
    unseen=set(range(cartan.nrows()))
    out=[]
    while unseen:
        a=min(unseen); unseen.remove(a)
        todo=[a]; comp=[]
        while todo:
            i=todo.pop(); comp.append(i)
            for j in list(unseen):
                if cartan[i,j]:
                    unseen.remove(j); todo.append(j)
        out.append(tuple(sorted(comp)))
    return tuple(sorted(out,key=lambda c:(len(c),c)))

def component_name(cartan,comp):
    B=cartan.matrix_from_rows_and_columns(comp,comp)
    r=B.nrows(); det=abs(ZZ(B.det())); cnt=ZZ(pari(B).qfminim(2)[0])
    if det==r+1 and cnt==r*(r+1): return f"A{r}"
    if r>=4 and det==4 and cnt==2*r*(r-1): return f"D{r}"
    return {(6,3,72):"E6",(7,2,126):"E7",(8,1,240):"E8"}.get(
        (r,det,cnt),f"R{r}d{det}n{cnt}"
    )

def ade_name(cartan):
    return "+".join(component_name(cartan,c) for c in connected_components(cartan))

def root_adaptation(child):
    roots,rb,data=roots_and_data(child)
    rr=data[0]
    if rr==0:
        L=matrix(ZZ,pari(child).qflllgram()).transpose()
        return L*child*L.transpose(),L
    smith,left,right=rb.smith_form()
    if tuple(abs(smith[i,i]) for i in range(rr))!=(1,)*rr:
        return None
    simple,cartan=deterministic_simple_roots(child)
    completion=right.inverse()
    initial=simple.stack(completion[rr:])
    if abs(initial.det())!=1:
        return None
    adapted=initial*child*initial.transpose()
    cp=adapted[:rr,rr:]
    tt=adapted[rr:,rr:]
    hh=tt-cp.transpose()*cartan.inverse()*cp
    scale=lcm(x.denominator() for x in hh.list())
    L=matrix(ZZ,pari((scale*hh).change_ring(ZZ)).qflllgram())
    change=block_diagonal_matrix(identity_matrix(ZZ,rr),L.transpose())
    basis=change*initial
    adapted=basis*child*basis.transpose()
    return adapted,basis

def dominant_weights_d13(bound):
    inv=C.inverse()
    out=defaultdict(list)
    # Recursive enumeration of nonnegative Dynkin labels l with l C^-1 l <= bound.
    def rec(prefix,norm):
        i=len(prefix)
        if i==13:
            out[QQ(norm)].append(tuple(prefix))
            return
        v=0
        while True:
            added=inv[i,i]*v*v
            added+=2*v*sum(inv[i,j]*prefix[j] for j in range(i))
            nn=norm+added
            if nn>bound:
                break
            rec(prefix+[v],nn)
            v+=1
    rec([],QQ(0))
    return out

qvals=[q for q in range(4,args.qmax+1,2)]
targetmax=2*max(qvals)

# MW vectors only up to the largest shell.
scale=lcm(x.denominator() for x in H.list())
IH=(scale*H).change_ring(ZZ)
qf=pari(IH).qfminim(scale*targetmax)
mwmap={}
for col in matrix(ZZ,qf[2]).columns():
    for sign in (1,-1):
        z=sign*vector(ZZ,col)
        if z==0 or z*H*z>targetmax:
            continue
        mwmap[tuple(z)]=z

mw=[]
for z in mwmap.values():
    h=QQ(z*H*z)
    corr,order=d13_correction(z)
    po=(h+corr-4)/2
    if po in ZZ and 0<=po<=args.max_pole:
        mw.append((ZZ(po),h,corr,order,z))
mw.sort(key=lambda x:(x[0],x[1],sum(abs(int(v)) for v in x[4]),tuple(x[4])))

print(
    "D13COST_MW|"
    f"qmax={args.qmax}|max_pole={args.max_pole}|"
    f"all_mw_vectors={len(mwmap)}|low_cost_mw={len(mw)}|status=PASS",
    flush=True,
)
for po,h,c,o,z in mw[:30]:
    print(
        "D13COST_MWCLASS|"
        f"PdotO={po}|height={h}|corr={c}|order={o}|"
        f"mw={','.join(map(str,z))}|status=PASS",
        flush=True,
    )

weights=dominant_weights_d13(QQ(targetmax))
print(
    "D13COST_ROOTWEIGHTS|"
    f"norm_classes={len(weights)}|weights={sum(map(len,weights.values()))}|"
    f"bound={targetmax}|status=PASS",
    flush=True,
)

hits=[]
tested=0
primitive=0
for q in qvals:
    target=ZZ(2*q)
    a=ZZ(q//2)
    for po,h,corr,order,z in mw:
        if h>target:
            continue
        need=QQ(target-h)
        for labels_tuple in weights.get(need,()):
            labels=vector(ZZ,labels_tuple)
            rc=C.inverse()*(labels-coupling*z)
            if not all(v in ZZ for v in rc):
                continue
            w=vector(ZZ,list(map(ZZ,rc))+list(z))
            assert w*G*w==target
            tested+=1
            fiber=vector(ZZ,[a,2]+list(w))
            result=child_frame(fiber)
            if result is None:
                continue
            primitive+=1
            child,nb=result
            rr,rcount=root_rank_count(child)
            mw_rank=17-rr
            if mw_rank<5:
                continue
            roots,rb,rdata=roots_and_data(child)
            adapted_data=root_adaptation(child)
            ade="not-adapted"
            adapted=None
            abasis=None
            if adapted_data is not None:
                adapted,abasis=adapted_data
                if rr:
                    ade=ade_name(adapted[:rr,:rr])
                else:
                    ade="rootless"
            rec={
                "q":int(q),
                "D_dot_O":int(a-2),
                "P_dot_O":int(po),
                "P_height":str(h),
                "P_correction":str(corr),
                "P_class_order":int(order),
                "P_mw":list(map(int,z)),
                "section_profile":{
                    "Z_degree":int(po),
                    "X_degree_bound":int(2*po+4),
                    "Y_degree_bound":int(3*po+6),
                },
                "dominant_labels":list(map(int,labels)),
                "witness":list(map(int,w)),
                "witness_L1":sum(abs(int(v)) for v in w),
                "child_root_data":list(map(int,rdata)),
                "child_mw_rank":int(mw_rank),
                "child_ade":ade,
                "child_frame":[list(map(int,row)) for row in child.rows()],
                "neighbor_basis":[list(map(int,row)) for row in nb.rows()],
            }
            if adapted is not None:
                rec["child_root_adapted_frame"]=[
                    list(map(int,row)) for row in adapted.rows()
                ]
                rec["child_root_adapted_basis"]=[
                    list(map(int,row)) for row in abasis.rows()
                ]
            hits.append(rec)
            print(
                "D13COST_HIT|"
                f"q={q}|PdotO={po}|Pprofile={po},{2*po+4},{3*po+6}|"
                f"child={ade}|root={rdata[0]},{rdata[1]},{rdata[2]}|"
                f"MW={mw_rank}|mw={','.join(map(str,z))}|"
                f"labels={','.join(map(str,labels))}|status=PASS",
                flush=True,
            )

hits.sort(key=lambda r:(
    r["P_dot_O"],r["q"],-r["child_mw_rank"],r["witness_L1"],
    tuple(r["P_mw"]),tuple(r["dominant_labels"])
))

print(
    "D13COST_RESULT|"
    f"tested={tested}|primitive={primitive}|rank_growth_hits={len(hits)}|"
    f"best="+(
        "NONE" if not hits else
        f"q{hits[0]['q']}:P.O{hits[0]['P_dot_O']}:{hits[0]['child_ade']}:MW{hits[0]['child_mw_rank']}"
    )+
    "|status="+("PASS_LOW_COST_FRONTIER" if hits else "NO_LOW_COST_GROWTH"),
    flush=True,
)

for i,r in enumerate(hits[:20]):
    print(
        "D13COST_BEST|"
        f"rank={i}|q={r['q']}|PdotO={r['P_dot_O']}|"
        f"profile={r['section_profile']['X_degree_bound']}/"
        f"{2*r['P_dot_O']},"
        f"{r['section_profile']['Y_degree_bound']}/"
        f"{3*r['P_dot_O']}|"
        f"child={r['child_ade']}|MW={r['child_mw_rank']}|"
        f"root={','.join(map(str,r['child_root_data']))}|"
        f"mw={','.join(map(str,r['P_mw']))}|status=PASS",
        flush=True,
    )

OUT=args.output.resolve() if args.output else LOCAL/"d13-equation-cost-frontier.json"
OUT.write_text(json.dumps({
    "schema":"elkies-k3.h3-d13-equation-cost-frontier.v1",
    "status":"PASS_D13_EQUATION_COST_FRONTIER",
    "qmax":args.qmax,
    "max_pole":args.max_pole,
    "tested_orbits":tested,
    "primitive_neighbors":primitive,
    "rank_growth_hits":hits,
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)

if args.frames_dir:
    outdir=args.frames_dir.resolve()
    outdir.mkdir(parents=True,exist_ok=True)
    for i,r in enumerate(hits[:20]):
        if "child_root_adapted_frame" not in r:
            continue
        path=outdir/f"cost{i:02d}-q{r['q']}-po{r['P_dot_O']}-{r['child_ade'].replace('+','_')}.txt"
        path.write_text(
            "\n".join(
                ["# "+json.dumps({
                    "q":r["q"],"P_dot_O":r["P_dot_O"],
                    "child_ade":r["child_ade"],
                    "child_mw_rank":r["child_mw_rank"],
                    "P_mw":r["P_mw"],
                },sort_keys=True)]
                +[" ".join(map(str,row)) for row in r["child_root_adapted_frame"]]
            )+"\n"
        )
