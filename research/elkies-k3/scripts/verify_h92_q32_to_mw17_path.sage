#!/usr/bin/env sage -python
import contextlib, io, json, sys
from pathlib import Path
from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector, xgcd

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
SCOUT=LOCAL/"route-scout"
OUT=SCOUT/"q32-to-mw17-exact-lattice-certificate.json"

U=matrix(ZZ,((0,1),(1,0)))

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

def root_data(gram):
    qf=pari(gram).qfminim(2)
    count=ZZ(qf[0])
    if count==0:
        return (0,0,1)
    roots=matrix(ZZ,qf[2]).transpose()
    basis=roots.row_module().basis_matrix()
    rg=basis*gram*basis.transpose()
    return (int(basis.rank()),int(count),abs(int(rg.det())))

def bezout(ns,fiber):
    cur=ZZ(0)
    out=[ZZ(0)]*ns.nrows()
    for i,v in enumerate(ns*fiber):
        if not v:
            continue
        g,a,b=xgcd(cur,ZZ(v))
        out=[a*x for x in out]
        out[i]+=b
        cur=g
    assert abs(cur)==1
    if cur==-1:
        out=[-x for x in out]
    return vector(ZZ,out)

def neighbor_child(frame,fiber):
    ns=block_diagonal_matrix(U,-frame)
    mate=bezout(ns,fiber)
    mate-=(mate*ns*mate//2)*fiber
    ker=matrix(ZZ,[list(fiber*ns),list(mate*ns)]).right_kernel_matrix()
    child=-(ker*ns*ker.transpose())
    assert child.det()==948 and child.is_positive_definite()
    return child

def choose_record(path,q,target_root,target_ade,target_mw):
    data=json.loads(path.read_text())
    assert data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    hits=[
        r for r in data["neighbors"]
        if int(r["q"])==q
        and tuple(r["child_root_data"])==target_root
        and r["child_ade"]==target_ade
        and int(r["child_mw_rank"])==target_mw
    ]
    assert hits,(path,target_root,target_ade,target_mw)
    hits.sort(key=lambda r:(sum(abs(int(x)) for x in r["witness"]),int(r["orbit_index"])))
    return hits[0]

# Exact q32 prefix.
Q32TEST=SCRIPTS/"test_h92_q24_geometrically_clean_d12_neighbor.sage"
saved=list(sys.argv)
scope={"__name__":"__embedded__"}
buf=io.StringIO()
try:
    sys.argv=[str(Q32TEST)]
    with contextlib.redirect_stdout(buf):
        exec(compile(Q32TEST.read_text(),str(Q32TEST),"exec"),scope)
finally:
    sys.argv=saved

assert tuple(scope["root_data"])==(12,264,4)
assert int(scope["mw_rank"])==5
assert int(scope["norm"])==64
assert int(scope["DdotO"])==14

current=load_gram(SCOUT/"q32-d12-mw5-root-adapted-frame.txt")
assert root_data(current)==(12,264,4)

ledger=[{
    "step":0,"q":32,"source":"D13/MW4","child":"D12/MW5",
    "root_data":[12,264,4],"mw_rank":5,"old_fibre_degree":2,
    "D_dot_O":14,"status":"PASS_EXACT_Q32_PREFIX"
}]
print("Q32MW17CERT|step=0|q=32|source=D13/MW4|child=D12/MW5|root_data=12,264,4|MW=5|degree=2|status=PASS",flush=True)

steps=[
    (SCOUT/"q32-d12-next-q4q6q8.json",6,(11,132,12),"A11",6,"D12/MW5","A11/MW6"),
    (SCOUT/"q32-a11-q8-degree2.json",8,(10,60,36),"A5+A5",7,"A11/MW6","2A5/MW7"),
    (SCOUT/"q32-corridor-step01-q4-3a3.json",4,(9,36,64),"A3+A3+A3",8,"2A5/MW7","3A3/MW8"),
    (SCOUT/"q32-corridor-step02-q4-a3-2a2.json",4,(7,24,36),"A2+A2+A3",10,"3A3/MW8","A3+2A2/MW10"),
    (SCOUT/"q32-corridor-step03-q4-5a1.json",4,(5,10,32),"A1+A1+A1+A1+A1",12,"A3+2A2/MW10","5A1/MW12"),
    (SCOUT/"q32-corridor-step04-q4-4a1.json",4,(4,8,16),"A1+A1+A1+A1",13,"5A1/MW12","4A1/MW13"),
    (SCOUT/"q32-fast-step05-q4-3a1.json",4,(3,6,8),"A1+A1+A1",14,"4A1/MW13","3A1/MW14"),
    (SCOUT/"q32-fast-step06-q4-2a1.json",4,(2,4,4),"A1+A1",15,"3A1/MW14","2A1/MW15"),
    (SCOUT/"q32-fast-step07-q4-a1.json",4,(1,2,2),"A1",16,"2A1/MW15","A1/MW16"),
]

for idx,(path,q,target_root,target_ade,target_mw,source_name,child_name) in enumerate(steps,1):
    r=choose_record(path,q,target_root,target_ade,target_mw)
    witness=vector(ZZ,r["witness"])
    assert witness*current*witness==2*q
    fiber=vector(ZZ,[q//2,2]+list(witness))
    child=neighbor_child(current,fiber)
    assert child==matrix(ZZ,r["child_frame"])
    A=matrix(ZZ,r["child_root_adapted_basis"])
    adapted=A*child*A.transpose()
    assert adapted==matrix(ZZ,r["child_root_adapted_frame"])
    assert root_data(adapted)==target_root
    current=adapted

    ledger.append({
        "step":idx,"q":q,"source":source_name,"child":child_name,
        "orbit":int(r["orbit_index"]),"root_data":list(target_root),
        "mw_rank":target_mw,"old_fibre_degree":2,
        "mw_projection":r["mw_projection"],"dominant_labels":r["dominant_labels"],
        "witness":r["witness"],"status":"PASS_EXACT_NEIGHBOR_REPLAY"
    })
    print(
        f"Q32MW17CERT|step={idx}|q={q}|source={source_name}|child={child_name}|"
        f"orbit={r['orbit_index']}|root_data={','.join(map(str,target_root))}|MW={target_mw}|degree=2|status=PASS",
        flush=True,
    )

PARITY=SCOUT/"q32-final-rootless-parity-scan.json"
pdat=json.loads(PARITY.read_text())
assert pdat["status"]=="PASS_Q32_ROOTLESS_BY_PARITY"
hit=pdat["hit"]
w=vector(ZZ,hit["witness"])
assert w*current*w==12
fiber=vector(ZZ,[3,2]+list(w))
child=neighbor_child(current,fiber)
assert root_data(child)==(0,0,1)

ledger.append({
    "step":10,"q":6,"source":"A1/MW16","child":"rootless/MW17",
    "root_data":[0,0,1],"mw_rank":17,"old_fibre_degree":2,
    "mw_projection":hit["mw"],"witness":hit["witness"],
    "parity":hit["parity"],"parity_scan_column":hit["stored_column_index"],
    "status":"PASS_EXACT_ROOTLESS_CHILD"
})
print(
    f"Q32MW17CERT|step=10|q=6|source=A1/MW16|child=rootless/MW17|"
    f"column={hit['stored_column_index']}|root_data=0,0,1|MW=17|degree=2|status=PASS",
    flush=True,
)

route=(
    "D13/MW4-q32-D12/MW5-q6-A11/MW6-q8-2A5/MW7-"
    "q4-3A3/MW8-q4-A3+2A2/MW10-q4-5A1/MW12-"
    "q4-4A1/MW13-q4-3A1/MW14-q4-2A1/MW15-"
    "q4-A1/MW16-q6-rootless/MW17"
)

payload={
    "schema":"elkies-k3.h3-q32-to-mw17-exact-lattice-certificate.v1",
    "status":"PASS_EXACT_Q32_TO_ROOTLESS_MW17_LATTICE_ROUTE",
    "route":route,
    "all_old_fibre_degree":2,
    "steps":ledger,
    "equation_boundary":{
        "q32_mod100003_resolved_RR":"PASS",
        "q32_mod100003_child":"I8* + 10 I1 / D12",
        "note":"Suffix is exact integral lattice; characteristic-zero equations remain to be constructed."
    },
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print("Q32MW17CERT_RESULT|steps=11|all_degree2=1|status=PASS_EXACT_Q32_TO_ROOTLESS_MW17_LATTICE_ROUTE",flush=True)
