#!/usr/bin/env sage -python
"""
Find the minimum-q24-degree exact integer combination of directly recoverable
explicit q24 multisection Abel--Jacobi directions that equals the exact
R17-directed orbit42 target.

This consumes the passing multisection-span audit.  It does NOT bound the
number of terms or coefficients.  Instead it solves the exact integer linear
problem

    sum_i c_i z_i = z_target

over the direct (recoverability <= 2) explicit directions, minimizing

    sum_i |c_i| * deg_q24(C_i)

with coefficient L1 as a secondary tie-breaker.

It also reports the degree-threshold filtration: the smallest maximum
individual q24 degree at which the direct directions generate the target.

No curve interpolation and no Groebner basis.
"""
import argparse, json
from pathlib import Path

from sage.all import (
    MixedIntegerLinearProgram, QQ, ZZ, matrix, vector
)

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)

IN=OUTDIR/f"orbit42-explicit-multisection-span-p{p}.json"
if not IN.exists():
    raise SystemExit(f"missing prerequisite: {IN}")

d=json.loads(IN.read_text())
assert d["status"] in (
    "PASS_Q24_ORBIT42_TARGET_IN_EXPLICIT_MULTISECTION_AJ_COMBINATION",
    "Q24_ORBIT42_TARGET_IN_DIRECT_MULTISECTION_SPAN_NEEDS_INTEGER_COMBINATION",
    "PASS_Q24_ORBIT42_TARGET_IN_BROAD_EXPLICIT_MULTISECTION_AJ_COMBINATION",
)
assert d["direct_span"]["target_in_span"] is True

target=vector(ZZ,d["target"]["mw"])
rows=[
    r for r in d["direction_representatives"]
    if int(r["recoverability"])<=2
]
assert rows

# Deduplicate again defensively.
uniq={}
for r in rows:
    key=tuple(map(int,r["mw_direction"]))
    score=(int(r["recoverability"]),int(r["q24_degree"]),str(r["name"]))
    if key not in uniq or score<uniq[key][0]:
        uniq[key]=(score,r)
rows=[v[1] for k,v in sorted(uniq.items())]

M=matrix(ZZ,[r["mw_direction"] for r in rows])
L=M.row_module()
assert target in L

print(
    "Q24O42DIRECT_INPUT|"
    f"prime={p}|directions={len(rows)}|rank={L.rank()}|"
    f"target={','.join(map(str,target))}|status=PASS",
    flush=True,
)

# ----------------------------------------------------------------------
# 1. Exact maximum-degree filtration.
# ----------------------------------------------------------------------
degrees=sorted(set(int(r["q24_degree"]) for r in rows))
filtration=[]
first_threshold=None
for D in degrees:
    sub=[r for r in rows if int(r["q24_degree"])<=D]
    S=matrix(ZZ,[r["mw_direction"] for r in sub])
    LS=S.row_module()
    inside=bool(target in LS)
    rec={
        "max_degree":int(D),
        "directions":len(sub),
        "rank":int(LS.rank()),
        "target_in_span":inside,
    }
    filtration.append(rec)
    print(
        "Q24O42DIRECT_FILTER|"
        f"max_degree={D}|directions={len(sub)}|rank={LS.rank()}|"
        f"target_in_span={int(inside)}|status=PASS",
        flush=True,
    )
    if inside and first_threshold is None:
        first_threshold=int(D)

assert first_threshold is not None

# ----------------------------------------------------------------------
# 2. Exact integer optimization over direct directions.
# ----------------------------------------------------------------------
#
# Use one objective with a large integer multiplier so weighted q24 degree
# dominates coefficient L1 lexicographically:
#
#   BIG * sum deg_i |c_i| + sum |c_i|.
#
# Since the secondary contribution is at most unbounded in principle, solve
# in two stages: first exact weighted degree, then pin it and minimize L1.
#
n=len(rows)
deg=[ZZ(r["q24_degree"]) for r in rows]
vec=[vector(ZZ,r["mw_direction"]) for r in rows]

def make_program():
    P=MixedIntegerLinearProgram(maximization=False,solver="GLPK")
    c=P.new_variable(integer=True,nonnegative=False)
    a=P.new_variable(integer=True,nonnegative=True)

    for j in range(5):
        P.add_constraint(
            sum(ZZ(vec[i][j])*c[i] for i in range(n)) == ZZ(target[j])
        )
    for i in range(n):
        P.add_constraint(a[i]-c[i],min=0)
        P.add_constraint(a[i]+c[i],min=0)
    return P,c,a

P,c,a=make_program()
weighted=sum(deg[i]*a[i] for i in range(n))
P.set_objective(weighted)
opt_weight=ZZ(round(P.solve()))
vals1=P.get_values(c)
coeff1=[ZZ(round(vals1[i])) for i in range(n)]
assert sum(coeff1[i]*vec[i] for i in range(n))==target
assert sum(abs(coeff1[i])*deg[i] for i in range(n))==opt_weight

# Secondary optimum: fix weighted degree, minimize coefficient L1.
P2,c2,a2=make_program()
weighted2=sum(deg[i]*a2[i] for i in range(n))
P2.add_constraint(weighted2==opt_weight)
l1expr=sum(a2[i] for i in range(n))
P2.set_objective(l1expr)
opt_l1=ZZ(round(P2.solve()))
vals2=P2.get_values(c2)
coeff=[ZZ(round(vals2[i])) for i in range(n)]

assert sum(coeff[i]*vec[i] for i in range(n))==target
assert sum(abs(coeff[i])*deg[i] for i in range(n))==opt_weight
assert sum(abs(coeff[i]) for i in range(n))==opt_l1

terms=[]
for i,ci in enumerate(coeff):
    if not ci:
        continue
    r=rows[i]
    terms.append({
        "coefficient":int(ci),
        "name":r["name"],
        "origin":r["origin"],
        "q24_degree":int(r["q24_degree"]),
        "recoverability":int(r["recoverability"]),
        "mw_direction":list(map(int,vec[i])),
        "weighted_degree_contribution":int(abs(ci)*deg[i]),
    })

assert terms

print(
    "Q24O42DIRECT_OPTIMUM|"
    f"weighted_degree={opt_weight}|coeff_L1={opt_l1}|"
    f"terms={len(terms)}|max_individual_degree={max(t['q24_degree'] for t in terms)}|"
    + "|".join(
        f"term{k+1}={t['coefficient']}*{t['name']}"
        f"(d{t['q24_degree']},mw={','.join(map(str,t['mw_direction']))})"
        for k,t in enumerate(terms)
    )
    + "|status=PASS_EXACT_INTEGER_OPTIMUM",
    flush=True,
)

# ----------------------------------------------------------------------
# 3. Practical route classification.
# ----------------------------------------------------------------------
#
# AJ interpolation cost is dominated by the degree of the involved
# multisections and the weighted multiplicity of repeated additions.
# We don't claim a hard runtime model here; expose the exact invariants.
#
maxdeg=max(t["q24_degree"] for t in terms)
if maxdeg<=100 and opt_weight<=500:
    recommendation="PURSUE_EXPLICIT_AJ"
elif maxdeg<=500 and opt_weight<=2000:
    recommendation="AJ_PLAUSIBLE_COMPARE_WITH_RESOLVED_RR"
else:
    recommendation="RESOLVED_RR_PREFERRED"

print(
    "Q24O42DIRECT_DECISION|"
    f"first_span_degree={first_threshold}|"
    f"opt_weighted_degree={opt_weight}|opt_coeff_L1={opt_l1}|"
    f"max_used_degree={maxdeg}|recommendation={recommendation}|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-orbit42-direct-aj-integer-optimum.v1",
    "status":"PASS_Q24_ORBIT42_DIRECT_AJ_INTEGER_OPTIMUM",
    "prime":int(p),
    "target_mw":list(map(int,target)),
    "direct_direction_count":len(rows),
    "direct_rank":int(L.rank()),
    "degree_filtration":filtration,
    "first_target_span_max_degree":first_threshold,
    "objective":{
        "primary":"sum(abs(coefficient)*q24_degree)",
        "secondary":"sum(abs(coefficient))",
        "optimal_weighted_degree":int(opt_weight),
        "optimal_coefficient_L1":int(opt_l1),
    },
    "optimal_terms":terms,
    "recommendation":recommendation,
    "proof_boundary":(
        "Exact integer Picard/MW combination optimization only.  It proves "
        "the displayed direct explicit curve combination has the orbit42 "
        "Abel--Jacobi MW class and is optimal for the stated weighted-degree "
        "objective.  It does not yet compute the fibrewise AJ rational points."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/f"orbit42-direct-aj-optimum-p{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42DIRECT_RESULT|"
    f"weighted_degree={opt_weight}|coeff_L1={opt_l1}|"
    f"recommendation={recommendation}|"
    "status=PASS_Q24_ORBIT42_DIRECT_AJ_INTEGER_OPTIMUM",
    flush=True,
)
