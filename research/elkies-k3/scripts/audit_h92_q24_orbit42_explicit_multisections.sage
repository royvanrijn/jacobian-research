#!/usr/bin/env sage -python
"""
Audit ALL already-explicit q8 curves as q24 multisections for the exact
R17-directed D12 -> A11 orbit42 target, using A0 as the q24 zero.

Unlike the older audit, this does not discard q24-degree > 1 curves.
For any curve C of q24 degree d, the generic divisor C|F has Abel--Jacobi
class [C-d*O]; in a root-adapted D12 frame its MW direction is the five-entry
MW tail of C's child coordinates.

This script:
  * loads the exact A0-zero orbit42 profile;
  * reuses the 42 exact explicit curve classes already collected by the old
    explicit-curve audit;
  * computes every positive-degree curve's A0-zero MW direction;
  * tests the full and equation-recoverable integral spans against the exact
    orbit42 target;
  * searches for a sparse small-coefficient combination (<=3 curves) whose
    Abel--Jacobi directions equal the target.

No Groebner basis and no section solving.
"""
import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--coeff-bound",type=int,default=8)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)
BOUND=int(args.coeff_bound)

SPIN=LOCAL/"q24-orbit42-spinor-zero-profiles.json"
OLD=OUTDIR/f"explicit-curves-a11-span-p{p}.json"
for path in (SPIN,OLD):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

spin=json.loads(SPIN.read_text())
old=json.loads(OLD.read_text())
assert spin["status"]=="PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES"
assert "explicit_curve_records" in old

profiles=[r for r in spin["profiles"] if r["zero"]=="A0"]
assert len(profiles)==1
ap=profiles[0]
assert ap["target_mw"]==[1,0,1,1,0],ap["target_mw"]
target=vector(ZZ,ap["target_mw"])

G=matrix(ZZ,ap["frame"])
B=matrix(ZZ,ap["parent_to_child_basis"])
U2=matrix(ZZ,((0,1),(1,0)))
Gchild=block_diagonal_matrix(U2,-G)
assert abs(ZZ(B.det()))==1
Binv=B.inverse()
assert Binv.change_ring(ZZ)==Binv
Binv=Binv.change_ring(ZZ)

# Recover the ambient NS pairing from B*ns*B^t = Gchild.
ns=Binv*Gchild*Binv.transpose()
assert ns.change_ring(ZZ)==ns
ns=ns.change_ring(ZZ)

F24=vector(ZZ,B.row(0))
Ochild=vector(ZZ,[-1,1]+[0]*17)
Oamb=vector(ZZ,Ochild*B)
assert F24*ns*F24==0
assert Oamb*ns*Oamb==-2 and Oamb*ns*F24==1

# Verify the zero is the explicit old-I9* A0 curve from the old audit.
a0_records=[
    rec for rec in old["explicit_curve_records"]
    if rec["name"]=="oldI9_A0"
]
assert len(a0_records)==1
assert vector(ZZ,a0_records[0]["class"])==Oamb

# Verify the exact D42 profile supplied by the spinor profiler.
D42=vector(ZZ,ap["D42_child_coordinates"])
assert D42[1]==2
assert vector(ZZ,D42[-5:])==target
assert D42*Gchild*D42==0

print(
    "Q24O42MULTI_INPUT|"
    f"prime={p}|zero=A0|target_mw={','.join(map(str,target))}|"
    f"target_height={ap['height']}|target_corr={ap['correction']}|"
    f"target_PdotO={ap['P_dot_O']}|curves={len(old['explicit_curve_records'])}|"
    "status=PASS",
    flush=True,
)

def recoverability(origin,name):
    # These classes come with equation/geometric constructions that can be
    # replayed to actual curves without inventing a new NS representative.
    if origin.startswith("q8_equation_certifier"):
        return 0
    if origin=="q24_close_replay":
        return 0
    if origin=="actual_effective_I9star":
        return 1
    if origin=="q8_component_nef_roots":
        return 2
    if origin=="H3_source_simple_roots":
        return 3
    return 4

records=[]
for rec in old["explicit_curve_records"]:
    C=vector(ZZ,rec["class"])
    assert len(C)==19
    square=ZZ(C*ns*C)
    q=vector(ZZ,C*Binv)
    degree=ZZ(q[1])
    # Old artifact degree was computed against the same q24 fibre in ambient
    # coordinates; enforce this coordinate regression.
    assert degree==ZZ(C*ns*F24)
    assert degree==ZZ(rec["q24_degree"])

    z=vector(ZZ,q[-5:])
    item={
        "name":rec["name"],
        "origin":rec["origin"],
        "square":int(square),
        "q24_degree":int(degree),
        "child_coordinates":list(map(int,q)),
        "mw_direction":list(map(int,z)),
        "mw_L1":sum(abs(int(x)) for x in z),
        "recoverability":recoverability(rec["origin"],rec["name"]),
        "class":list(map(int,C)),
    }
    records.append(item)

positive=[
    r for r in records
    if r["square"]==-2 and r["q24_degree"]>0
]
nonzero=[
    r for r in positive
    if any(r["mw_direction"])
]

# Deduplicate exact MW directions, retaining the cheapest geometric
# representative for each direction.
by_direction={}
for r in nonzero:
    key=tuple(r["mw_direction"])
    score=(r["recoverability"],r["q24_degree"],r["name"])
    if key not in by_direction or score < by_direction[key][0]:
        by_direction[key]=(score,r)

options=[entry[1] for key,entry in sorted(by_direction.items())]
direct=[r for r in options if r["recoverability"]<=2]

def span_info(rows):
    if not rows:
        return 0,False,[],[]
    M=matrix(ZZ,[r["mw_direction"] for r in rows])
    L=M.row_module()
    inside=bool(target in L)
    basis=L.basis_matrix()
    coords=[]
    if inside:
        coords=list(map(int,L.coordinate_vector(target)))
        assert vector(ZZ,coords)*basis==target
    return int(L.rank()),inside,[
        list(map(int,row)) for row in basis.rows()
    ],coords

rank_all,in_all,basis_all,coords_all=span_info(options)
rank_direct,in_direct,basis_direct,coords_direct=span_info(direct)

print(
    "Q24O42MULTI_SPAN|"
    f"positive_curves={len(positive)}|nonzero={len(nonzero)}|"
    f"distinct_directions={len(options)}|rank={rank_all}|"
    f"target_in_span={int(in_all)}|"
    f"direct_directions={len(direct)}|direct_rank={rank_direct}|"
    f"target_in_direct_span={int(in_direct)}|status=PASS",
    flush=True,
)

# Surface the cheapest individual AJ directions.
for r in sorted(
    options,
    key=lambda x:(x["recoverability"],x["q24_degree"],x["mw_L1"],x["name"])
)[:16]:
    z=vector(ZZ,r["mw_direction"])
    scalar=None
    for c in range(-BOUND,BOUND+1):
        if c and c*z==target:
            scalar=c
            break
    print(
        "Q24O42MULTI_CURVE|"
        f"name={r['name']}|origin={r['origin']}|degree={r['q24_degree']}|"
        f"mw={','.join(map(str,z))}|recoverability={r['recoverability']}|"
        f"scalar_target={scalar if scalar is not None else 'NO'}|status=PROFILE",
        flush=True,
    )

def combo_score(terms):
    # terms = [(option_index, coefficient), ...]
    weighted=sum(
        abs(c)*options[i]["q24_degree"] for i,c in terms
    )
    coeff_l1=sum(abs(c) for i,c in terms)
    maxdeg=max(options[i]["q24_degree"] for i,c in terms)
    recov=max(options[i]["recoverability"] for i,c in terms)
    return (recov,weighted,len(terms),coeff_l1,maxdeg)

def search_sparse(pool_indices):
    """
    Search exact combinations with 1, 2, or 3 distinct MW directions and
    coefficients in [-BOUND,BOUND]\\{0}.  Pair sums are hashed, so the cubic
    search remains small.
    """
    coeffs=[c for c in range(-BOUND,BOUND+1) if c]
    best=None

    # one term
    for i in pool_indices:
        z=vector(ZZ,options[i]["mw_direction"])
        for c in coeffs:
            if c*z==target:
                terms=[(i,c)]
                score=combo_score(terms)
                if best is None or score<best[0]:
                    best=(score,terms)

    # two-term sums, and cache for triples
    pair_best={}
    for ai,i in enumerate(pool_indices):
        zi=vector(ZZ,options[i]["mw_direction"])
        for j in pool_indices[ai+1:]:
            zj=vector(ZZ,options[j]["mw_direction"])
            for ci in coeffs:
                vi=ci*zi
                for cj in coeffs:
                    val=vi+cj*zj
                    key=tuple(map(int,val))
                    terms=[(i,ci),(j,cj)]
                    score=combo_score(terms)
                    oldp=pair_best.get(key)
                    if oldp is None or score<oldp[0]:
                        pair_best[key]=(score,terms)
                    if val==target and (best is None or score<best[0]):
                        best=(score,terms)

    # three terms via target - c*z lookup
    for k in pool_indices:
        zk=vector(ZZ,options[k]["mw_direction"])
        for ck in coeffs:
            need=tuple(map(int,target-ck*zk))
            hit=pair_best.get(need)
            if hit is None:
                continue
            _,pair_terms=hit
            if any(i==k for i,c in pair_terms):
                continue
            terms=pair_terms+[(k,ck)]
            score=combo_score(terms)
            if best is None or score<best[0]:
                best=(score,terms)

    return best

direct_indices=[
    i for i,r in enumerate(options) if r["recoverability"]<=2
]
all_indices=list(range(len(options)))

best_direct=search_sparse(direct_indices)
best_all=best_direct if best_direct is not None else search_sparse(all_indices)

def combo_payload(found):
    if found is None:
        return None
    score,terms=found
    check=vector(ZZ,[0]*5)
    out=[]
    for i,c in terms:
        r=options[i]
        z=vector(ZZ,r["mw_direction"])
        check += c*z
        out.append({
            "coefficient":int(c),
            "name":r["name"],
            "origin":r["origin"],
            "q24_degree":int(r["q24_degree"]),
            "recoverability":int(r["recoverability"]),
            "mw_direction":list(map(int,z)),
        })
    assert check==target
    return {
        "score":list(map(int,score)),
        "terms":out,
        "coefficient_bound":BOUND,
    }

combo_direct=combo_payload(best_direct)
combo_all=combo_payload(best_all)

if combo_direct:
    chosen=combo_direct
    chosen_kind="DIRECT"
elif combo_all:
    chosen=combo_all
    chosen_kind="ALL_EXPLICIT"
else:
    chosen=None
    chosen_kind="NONE"

if chosen:
    print(
        "Q24O42MULTI_COMBO|"
        f"kind={chosen_kind}|terms={len(chosen['terms'])}|"
        f"weighted_degree={chosen['score'][1]}|"
        f"coeff_L1={chosen['score'][3]}|"
        + "|".join(
            f"term{j+1}={t['coefficient']}*{t['name']}"
            f"(d{t['q24_degree']},mw={','.join(map(str,t['mw_direction']))})"
            for j,t in enumerate(chosen["terms"])
        )
        + "|status=PASS_TARGET_COMBINATION",
        flush=True,
    )
else:
    print(
        "Q24O42MULTI_COMBO|"
        f"kind=NONE|coeff_bound={BOUND}|"
        f"target_in_span={int(in_all)}|target_in_direct_span={int(in_direct)}|"
        "status=NO_SPARSE_COMBINATION",
        flush=True,
    )

if combo_direct:
    terminal="PASS_Q24_ORBIT42_TARGET_IN_EXPLICIT_MULTISECTION_AJ_COMBINATION"
elif in_direct:
    terminal="Q24_ORBIT42_TARGET_IN_DIRECT_MULTISECTION_SPAN_NEEDS_INTEGER_COMBINATION"
elif combo_all:
    terminal="PASS_Q24_ORBIT42_TARGET_IN_BROAD_EXPLICIT_MULTISECTION_AJ_COMBINATION"
elif in_all:
    terminal="Q24_ORBIT42_TARGET_IN_BROAD_MULTISECTION_SPAN_ONLY"
else:
    terminal="Q24_ORBIT42_TARGET_NOT_IN_EXPLICIT_MULTISECTION_SPAN"

payload={
    "schema":"elkies-k3.h3-q24-orbit42-explicit-multisection-span.v1",
    "status":terminal,
    "prime":int(p),
    "zero":"A0",
    "target":{
        "mw":list(map(int,target)),
        "height":ap["height"],
        "correction":ap["correction"],
        "P_dot_O":int(ap["P_dot_O"]),
    },
    "counts":{
        "all_curve_records":len(records),
        "positive_minus2_curves":len(positive),
        "nonzero_AJ_curves":len(nonzero),
        "distinct_AJ_directions":len(options),
        "direct_AJ_directions":len(direct),
    },
    "full_span":{
        "rank":rank_all,
        "target_in_span":bool(in_all),
        "basis":basis_all,
        "target_coordinates":coords_all,
    },
    "direct_span":{
        "rank":rank_direct,
        "target_in_span":bool(in_direct),
        "basis":basis_direct,
        "target_coordinates":coords_direct,
    },
    "curve_profiles":records,
    "direction_representatives":options,
    "best_direct_sparse_combination":combo_direct,
    "best_sparse_combination":combo_all,
    "next":(
        "Recover fibrewise Abel-Jacobi sections for the chosen explicit "
        "multisections and combine them on the pointed q24 Jacobian."
        if combo_direct else
        "If the direct span fails, compile the exact D42 degree-two divisor "
        "by resolved D12/I8* Riemann-Roch rather than solving a P.O=3 section."
    ),
    "proof_boundary":(
        "Exact NS/Picard-group span audit only. A passing sparse combination "
        "identifies explicit curve classes whose fibrewise Abel-Jacobi sum has "
        "the exact orbit42 MW class. It does not yet compute those AJ rational "
        "sections or the q6 A11 equation."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/f"orbit42-explicit-multisection-span-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42MULTI_RESULT|"
    f"rank={rank_all}|direct_rank={rank_direct}|"
    f"target_in_span={int(in_all)}|target_in_direct_span={int(in_direct)}|"
    f"sparse={int(chosen is not None)}|status={terminal}",
    flush=True,
)
