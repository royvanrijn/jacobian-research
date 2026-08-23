#!/usr/bin/env sage -python
"""
Recover the exact abstract D13-frame isometry between the native q24 frame and
the current equation q8 frame structurally.

Avoids a 17-dimensional qfisom:
  - root block: the two D13 Dynkin graph identifications;
  - MW block: enumerate all exact GL4(Z) isometries of the two height forms;
  - coupling block then forces the 4x13 root shift.

Filter the resulting exact 17x17 isometries by canonical q8 anchors:
  G1 = IV*_E6_1,
  AJ(old_E7_7) = G1+G3.

Then determine where native q24=(...) lands in the equation MW basis and replay
the complete native q24 divisor in the equation ambient frame.
"""

import json
import sys
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, block_matrix, gcd, identity_matrix,
    lcm, matrix, pari, vector, zero_matrix
)


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
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
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
EQ=SCRIPTS/"certify_h92_q8_equation_ns_divisor.sage"
Q24=SCRIPTS/"audit_h92_q8_q24_effective_zero_choices.sage"
OUT=LOCAL/"q8-q24-structural-d13-isometry.json"
TMP=LOCAL/"q8-q24-structural-native-temp.json"


def run(path, argv):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]+list(argv)
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


print("Q8Q24ISO|stage=equation",flush=True)
eq=run(EQ,[])
print("Q8Q24ISO|stage=native",flush=True)
qp=run(Q24,["--output",str(TMP)])

child_frame_with_zero=qp["child_frame_with_zero"]
d13_root_adaptation=qp["d13_root_adaptation"]

# ---------------------------------------------------------------------------
# Build the two D13 frames.
# ---------------------------------------------------------------------------

nsE=eq["ns"]
FE=vector(ZZ,eq["F8eq"])
OE=vector(ZZ,eq["selected"]["E8"]["simple_root_vectors_in_source_h3_ns"][0])
S3E=vector(ZZ,eq["S3std"])
assert OE*nsE*FE==1
childE,BzeroE=child_frame_with_zero(nsE,FE,OE)
AE,GEpos,HE=d13_root_adaptation(childE)
BadaptE=block_diagonal_matrix(identity_matrix(ZZ,2),AE)*BzeroE
assert BadaptE*nsE*BadaptE.transpose()==block_diagonal_matrix(
    matrix(ZZ,((0,1),(1,0))),-GEpos
)

nsN=qp["source_ns"]
FN=vector(ZZ,qp["F_actual"])
DN=vector(ZZ,qp["D_actual"])
profileN=next(p for p in qp["profiles"] if p["zero"]=="II*_E8_1")
ON=vector(ZZ,profileN["zero_source_h3_ns"])
PN=vector(ZZ,profileN["section_source_h3_ns"])
assert ON*nsN*FN==1
childN,BzeroN=child_frame_with_zero(nsN,FN,ON)
AN,GNpos,HN=d13_root_adaptation(childN)
BadaptN=block_diagonal_matrix(identity_matrix(ZZ,2),AN)*BzeroN
assert BadaptN*nsN*BadaptN.transpose()==block_diagonal_matrix(
    matrix(ZZ,((0,1),(1,0))),-GNpos
)


def coords(C,B):
    c=vector(QQ,C)*B.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ,c)


cS3E=coords(S3E,BadaptE)
zS3E=vector(ZZ,cS3E[-4:])
cDN=coords(DN,BadaptN)
zDN=vector(ZZ,cDN[-4:])
cPN=coords(PN,BadaptN)
zPN=vector(ZZ,cPN[-4:])

print(
    "Q8Q24ISO_INPUT|"
    f"AJ_eq={','.join(map(str,zS3E))}|"
    f"q24_native={','.join(map(str,zDN))}|"
    f"q24_section_native={','.join(map(str,zPN))}|"
    f"height_det_eq={HE.det()}|height_det_native={HN.det()}|status=PASS",
    flush=True,
)
assert zDN==zPN

# ---------------------------------------------------------------------------
# D13 graph isometries: exactly two, swapping the two spinor leaves.
# Rows of A express the native root basis in the equation root basis:
#     A * R_E * A^T = R_N.
# ---------------------------------------------------------------------------

RE=GEpos[:13,:13]
RN=GNpos[:13,:13]
CE=GEpos[:13,13:]
CN=GNpos[:13,13:]

def d13_structure(R):
    n=R.nrows()
    assert n==13
    assert all(R[i,i]==2 for i in range(n))
    adj=[]
    for i in range(n):
        nbr=[]
        for j in range(n):
            if i==j: continue
            assert R[i,j] in (0,-1)
            if R[i,j]==-1:
                nbr.append(j)
        adj.append(nbr)
    branch=[i for i,a in enumerate(adj) if len(a)==3]
    assert len(branch)==1
    b=branch[0]
    leaves=[j for j in adj[b] if len(adj[j])==1]
    chain_start=[j for j in adj[b] if len(adj[j])==2]
    assert len(leaves)==2 and len(chain_start)==1
    chain=[]
    prev=b
    cur=chain_start[0]
    while True:
        chain.append(cur)
        nxt=[j for j in adj[cur] if j!=prev]
        if not nxt:
            break
        assert len(nxt)==1
        prev,cur=cur,nxt[0]
    assert len(chain)==10
    assert len(adj[chain[-1]])==1
    return b,tuple(leaves),tuple(chain)

bN,lN,chN=d13_structure(RN)
bE,lE,chE=d13_structure(RE)

root_maps=[]
for swap in (False,True):
    dest={}
    dest[bN]=bE
    for i,j in zip(chN,chE):
        dest[i]=j
    leavesE=(lE[1],lE[0]) if swap else lE
    for i,j in zip(lN,leavesE):
        dest[i]=j
    assert len(dest)==13
    A=zero_matrix(ZZ,13,13)
    for i,j in dest.items():
        A[i,j]=1
    assert abs(A.det())==1
    assert A*RE*A.transpose()==RN
    root_maps.append((swap,A))

print("Q8Q24ISO_ROOT|maps=2|status=PASS_D13_GRAPH",flush=True)

# ---------------------------------------------------------------------------
# Enumerate ALL 4D MW isometries Q with Q*H_E*Q^T = H_N.
# Rank 4 makes exact bounded enumeration tiny compared with qfisom in rank 17.
# ---------------------------------------------------------------------------

scale=ZZ(1)
for x in HE.list() + HN.list():
    scale=lcm(scale,ZZ(QQ(x).denominator()))
QE=(scale*HE).change_ring(ZZ)
QN=(scale*HN).change_ring(ZZ)

target_norms=[ZZ(QN[i,i]) for i in range(4)]
maxnorm=max(target_norms)

res=pari(QE).qfminim(maxnorm)
half=[vector(ZZ,c) for c in matrix(ZZ,res[2]).columns()]
vectors=[]
for v in half:
    vectors.append(v)
    vectors.append(-v)

by_norm={}
for v in vectors:
    n=ZZ(v*QE*v)
    by_norm.setdefault(n,[]).append(v)

for n in target_norms:
    if n not in by_norm:
        raise ArithmeticError(f"no equation MW vectors of required norm {n}")

mw_isos=[]
rows=[]

def backtrack(i):
    if i==4:
        Q=matrix(ZZ,[list(r) for r in rows])
        if abs(Q.det())!=1:
            return
        if Q*QE*Q.transpose()!=QN:
            return
        mw_isos.append(Q)
        return
    for v in by_norm[target_norms[i]]:
        ok=True
        for j,u in enumerate(rows):
            if v*QE*u != QN[i,j]:
                ok=False
                break
        if not ok:
            continue
        rows.append(v)
        backtrack(i+1)
        rows.pop()

backtrack(0)

# Deduplicate, though qfminim returns each +/- pair once and we explicitly add signs.
uniq={}
for Q in mw_isos:
    uniq[tuple(Q.list())]=Q
mw_isos=list(uniq.values())

print(
    "Q8Q24ISO_MW|"
    f"scale={scale}|qfminim_half={len(half)}|"
    f"mw_isometries={len(mw_isos)}|status=PASS",
    flush=True,
)
if not mw_isos:
    raise ArithmeticError("no MW height-form isometry found")

# ---------------------------------------------------------------------------
# Marking-independent scalar anchors.
#
# Do NOT identify selected E6/E8 component vectors across the two ambient
# markings.  Instead use scalar intersection data of the native q24 horizontal
# section, which survives any basis change.
# ---------------------------------------------------------------------------

Paju,haju,ordaju,corraju,poleaju = qp["section_for_mw"](GEpos,zS3E,32768)
PAJ = vector(ZZ,Paju) * BadaptE
assert PAJ*nsE*PAJ == -2 and PAJ*nsE*FE == 1

native_scalar = {
    "q6_degree": int(profileN["P_q6_degree"]),
    "q6_oldzero": int(profileN["P_q6_oldzero_intersection"]),
    "h3_degree": int(profileN["P_H3_degree"]),
    "PdotO": int(profileN["P_dot_O"]),
    "height": str(profileN["height"]),
    "correction": str(profileN["D13_correction"]),
}

equation_AJ_scalar = {
    "q6_degree": int(PAJ * nsE * eq["F6"]),
    "q6_oldzero": int(PAJ * nsE * eq["Oold"]),
    "h3_degree": int(PAJ * nsE * eq["Fh3"]),
    "PdotO": int(PAJ * nsE * OE),
    "height": str(haju),
    "correction": str(corraju),
}

print(
    "Q8Q24ISO_SCALAR_AJ|"
    f"native_q6_degree={native_scalar['q6_degree']}|"
    f"equation_q6_degree={equation_AJ_scalar['q6_degree']}|"
    f"native_q6_oldzero={native_scalar['q6_oldzero']}|"
    f"equation_q6_oldzero={equation_AJ_scalar['q6_oldzero']}|"
    f"native_h3_degree={native_scalar['h3_degree']}|"
    f"equation_h3_degree={equation_AJ_scalar['h3_degree']}|"
    f"native_PdotO={native_scalar['PdotO']}|"
    f"equation_PdotO={equation_AJ_scalar['PdotO']}|"
    f"native_height={native_scalar['height']}|equation_height={equation_AJ_scalar['height']}|"
    f"native_corr={native_scalar['correction']}|equation_corr={equation_AJ_scalar['correction']}|"
    f"all_match={int(native_scalar == equation_AJ_scalar)}|"
    "status=PASS_DIAGNOSTIC",
    flush=True,
)

# ---------------------------------------------------------------------------
# Coupling forces B for each (A,Q):
#
# T = [ A  0 ]
#     [ B  Q ]
#
# and T G_E T^T = G_N.
# ---------------------------------------------------------------------------

full=[]
for swap,A in root_maps:
    Ainv=A.inverse()
    for qi,Q in enumerate(mw_isos):
        Bt=RE.inverse()*(Ainv*CN - CE*Q.transpose())
        if not all(x in ZZ for x in Bt.list()):
            continue
        B=matrix(ZZ,Bt.transpose())
        T=block_matrix(ZZ,[
            [A,zero_matrix(ZZ,13,4)],
            [B,Q],
        ])
        if abs(T.det())!=1:
            continue
        if T*GEpos*T.transpose()!=GNpos:
            continue

        full.append({
            "swap":swap,
            "mw_index":qi,
            "A":A,
            "Q":Q,
            "B":B,
            "T":T,
            "mapped_q24":vector(ZZ,zDN*Q),
            "maps_q24_to_AJ":bool(zDN*Q == zS3E),
        })

print(
    "Q8Q24ISO_FULL|"
    f"integral_full={len(full)}|"
    f"q24AJ_any={sum(r['maps_q24_to_AJ'] for r in full)}|"
    f"distinct_q24_images={len(set(tuple(r['mapped_q24']) for r in full))}|"
    "status=PASS",
    flush=True,
)

mapped_all=sorted(set(tuple(r["mapped_q24"]) for r in full))
print(
    "Q8Q24ISO_MAPPED|all_q24_images="
    +";".join(",".join(map(str,x)) for x in mapped_all)
    +"|AJ="+",".join(map(str,zS3E)),
    flush=True,
)

# ---------------------------------------------------------------------------
# Exhaustive exact replay.
#
# No extra ambient anchor is imposed. For every exact full D13 isometry:
#   * replay the native q24 horizontal section;
#   * compare the FULL (-2)-curve with the independently reconstructed
#     equation AJ(S3) effective section;
#   * replay the entire native q24 isotropic divisor;
#   * classify its orthogonal child.
#
# This distinguishes:
#   MW_match       : only Pic^0 coordinates agree;
#   section_match  : the exact effective section curve agrees in NS;
#   D12            : the full q24 divisor remains the rank-growing neighbour.
# ---------------------------------------------------------------------------

exact_replays=[]

for index,r in enumerate(full,1):
    T=r["T"]

    # Native q24 horizontal section -> equation frame.
    wPN=vector(ZZ,cPN[2:])
    wPE=vector(ZZ,wPN*T)
    cPE=vector(ZZ,[cPN[0],cPN[1]]+list(wPE))
    PE=cPE*BadaptE
    assert PE*nsE*PE==-2
    assert PE*nsE*FE==1

    zPE=vector(ZZ,cPE[-4:])
    mw_match=bool(zPE==zS3E)
    section_match=bool(PE==PAJ)

    # Native full q24 divisor -> equation frame.
    wDN=vector(ZZ,cDN[2:])
    wDE=vector(ZZ,wDN*T)
    cDE=vector(ZZ,[cDN[0],cDN[1]]+list(wDE))
    DE=cDE*BadaptE

    assert DE*nsE*DE==0
    assert DE*nsE*FE==2
    assert gcd(tuple(DE))==1

    mate=eq["isotropic_mate"](nsE,DE)
    orth=matrix(
        ZZ,[list(DE*nsE),list(mate*nsE)]
    ).right_kernel_matrix()
    child=-(orth*nsE*orth.transpose())
    rd=eq["roots_and_data"](child)[2]
    d12=bool(rd==(12,264,4))

    V=DE-OE-PE

    rec={
        "index":index,
        "swap":bool(r["swap"]),
        "mw_index":int(r["mw_index"]),
        "Q":r["Q"],
        "B":r["B"],
        "T":T,
        "mapped_q24":vector(ZZ,r["mapped_q24"]),
        "cPE":cPE,
        "PE":PE,
        "mw_match":mw_match,
        "section_match":section_match,
        "cDE":cDE,
        "DE":DE,
        "V":V,
        "rd":rd,
        "d12":d12,
    }
    exact_replays.append(rec)

    print(
        "Q8Q24ISO_EXACT_REPLAY|"
        f"index={index}|spin_swap={int(r['swap'])}|mw_iso={r['mw_index']}|"
        f"mapped_q24={','.join(map(str,r['mapped_q24']))}|"
        f"mw_match={int(mw_match)}|section_match={int(section_match)}|"
        f"child={rd[0]},{rd[1]},{rd[2]}|D12={int(d12)}|"
        f"V2={V*nsE*V}|VO={V*nsE*OE}|VP={V*nsE*PE}|"
        "status=PASS",
        flush=True,
    )

mw_matches=[r for r in exact_replays if r["mw_match"]]
section_matches=[r for r in exact_replays if r["section_match"]]
d12_matches=[r for r in exact_replays if r["d12"]]
closure=[r for r in exact_replays if r["section_match"] and r["d12"]]

print(
    "Q8Q24ISO_EXHAUSTIVE|"
    f"full={len(exact_replays)}|"
    f"mw_matches={len(mw_matches)}|"
    f"section_matches={len(section_matches)}|"
    f"D12={len(d12_matches)}|"
    f"section_and_D12={len(closure)}|"
    "status=PASS",
    flush=True,
)

# Summarize the MW maps explicitly: with only a handful, this is more useful
# than hiding them behind another anchor heuristic.
for r in exact_replays:
    print(
        "Q8Q24ISO_MW_MAP|"
        f"index={r['index']}|"
        f"Q_rows={';'.join(','.join(map(str,row)) for row in r['Q'].rows())}|"
        f"native_q24={','.join(map(str,zDN))}|"
        f"image={','.join(map(str,r['mapped_q24']))}|"
        f"AJ={','.join(map(str,zS3E))}|"
        f"section_match={int(r['section_match'])}|D12={int(r['d12'])}",
        flush=True,
    )

if not closure:
    payload={
        "schema":"elkies-k3.h92-q8-q24-structural-d13-isometry.v2",
        "status":"NO_EXACT_SECTION_D12_REPLAY",
        "counts":{
            "root_graph_maps":2,
            "mw_isometries":len(mw_isos),
            "full_integral_isometries":len(full),
            "mw_matches":len(mw_matches),
            "exact_section_matches":len(section_matches),
            "D12_replays":len(d12_matches),
            "section_and_D12":0,
        },
        "horizontal":{
            "native_q24_mw":list(map(int,zDN)),
            "equation_AJ_S3_mw":list(map(int,zS3E)),
            "images":[list(map(int,x)) for x in mapped_all],
        },
        "replays":[
            {
                "index":r["index"],
                "spin_swap":r["swap"],
                "mw_index":r["mw_index"],
                "mw_matrix_rows":[list(map(int,row)) for row in r["Q"].rows()],
                "mapped_q24":list(map(int,r["mapped_q24"])),
                "mw_match":r["mw_match"],
                "section_match":r["section_match"],
                "child_root_data":list(map(int,r["rd"])),
                "D12":r["d12"],
            }
            for r in exact_replays
        ],
        "interpretation":(
            "All exact structural D13 isometries were exhausted. No q6/H3 "
            "ambient scalar was used as an anchor."
        ),
    }
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print(
        "Q8Q24ISO_RESULT|"
        f"full={len(full)}|mw_match={len(mw_matches)}|"
        f"section_match={len(section_matches)}|D12={len(d12_matches)}|"
        "status=NO_EXACT_SECTION_D12_REPLAY",
        flush=True,
    )
    raise SystemExit(2)

# Multiple closure maps can differ only by automorphisms invisible on the
# selected q24 section/divisor. Deduplicate by the actual equation NS classes.
classes={}
for r in closure:
    key=(tuple(r["PE"]),tuple(r["DE"]))
    classes.setdefault(key,[]).append(r)

print(
    "Q8Q24ISO_CLOSURE_CLASSES|"
    f"maps={len(closure)}|distinct_NS_pairs={len(classes)}|status=PASS",
    flush=True,
)

# If every closure map gives the same actual section/divisor pair, the residual
# frame automorphism ambiguity is harmless. If there are several pairs, retain
# all rather than silently choose one.
canonical=list(classes.values())[0][0]
unique_pair=(len(classes)==1)

best=canonical

payload={
    "schema":"elkies-k3.h92-q8-q24-structural-d13-isometry.v2",
    "status":(
        "PASS_EXACT_Q24_AJ_STRUCTURAL_EQUIVALENCE"
        if unique_pair
        else "PASS_MULTIPLE_Q24_AJ_STRUCTURAL_EQUIVALENCES"
    ),
    "counts":{
        "root_graph_maps":2,
        "mw_isometries":len(mw_isos),
        "full_integral_isometries":len(full),
        "mw_matches":len(mw_matches),
        "exact_section_matches":len(section_matches),
        "D12_replays":len(d12_matches),
        "section_and_D12":len(closure),
        "distinct_equation_NS_pairs":len(classes),
    },
    "horizontal":{
        "native_q24_mw":list(map(int,zDN)),
        "equation_AJ_S3_mw":list(map(int,zS3E)),
        "mapped_equal":True,
        "exact_section_equal":True,
    },
    "representative_isometry":{
        "spin_swap":bool(best["swap"]),
        "mw_index":int(best["mw_index"]),
        "mw_matrix_rows":[list(map(int,row)) for row in best["Q"].rows()],
        "root_shift_rows":[list(map(int,row)) for row in best["B"].rows()],
        "full_positive_frame_matrix_rows":[list(map(int,row)) for row in best["T"].rows()],
    },
    "q24_equation":{
        "adapted_coordinates":list(map(int,best["cDE"])),
        "divisor_source_h3_ns":list(map(int,best["DE"])),
        "section_adapted_coordinates":list(map(int,best["cPE"])),
        "section_source_h3_ns":list(map(int,best["PE"])),
        "section_equals_equation_AJ_effective_lift":True,
        "vertical_source_h3_ns":list(map(int,best["V"])),
        "vertical_square":int(best["V"]*nsE*best["V"]),
        "vertical_dot_O":int(best["V"]*nsE*OE),
        "vertical_dot_P":int(best["V"]*nsE*best["PE"]),
        "child_root_data":list(map(int,best["rd"])),
        "child_root_lattice":"D12",
        "MW_rank_if_rho19":5,
    },
    "all_closure_maps":[
        {
            "index":r["index"],
            "spin_swap":r["swap"],
            "mw_index":r["mw_index"],
            "mw_matrix_rows":[list(map(int,row)) for row in r["Q"].rows()],
            "equation_divisor_source_h3_ns":list(map(int,r["DE"])),
            "equation_section_source_h3_ns":list(map(int,r["PE"])),
        }
        for r in closure
    ],
    "boundary":(
        "This proves existence of an exact abstract D13 lattice isometry taking "
        "the native q24 horizontal section to the independently reconstructed "
        "equation AJ(S3) section and the full q24 divisor to a D12-producing "
        "degree-two isotropic divisor. It does not identify the historically "
        "mismatched ambient source-H3 markings themselves."
    ),
}

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q8Q24ISO_RESULT|"
    f"full={len(full)}|mw_match={len(mw_matches)}|"
    f"section_match={len(section_matches)}|D12={len(d12_matches)}|"
    f"closure={len(closure)}|distinct_NS={len(classes)}|"
    f"status={payload['status']}",
    flush=True,
)
