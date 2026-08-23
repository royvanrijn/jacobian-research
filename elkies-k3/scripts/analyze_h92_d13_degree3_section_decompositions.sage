#!/usr/bin/env sage -python
"""
Rank completed degree-3 D13 neighbours by actual cubic-divisor equation cost.

For a degree-3 divisor D on the old elliptic fibre, its MW projection z is
the Pic^0 class of D-3O.  If
    z = z1 + z2 + z3
with zi represented by low-pole sections, then D can be built from three
small section points instead of one enormous horizontal point.

We enumerate all D13 section classes with P.O <= --max-pole, then find the
best exact 3-section decomposition for every degree-3 rank-growing neighbour.
"""

import argparse, json
from pathlib import Path
from sage.all import QQ, ZZ, lcm, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
SCOUT=LOCAL/"route-scout"
FRAME=ROOT/"elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
INPUT=SCOUT/"d13-degree3-qle40-rank-growth.json"

parser=argparse.ArgumentParser()
parser.add_argument("--max-pole",type=int,default=8)
parser.add_argument("--top",type=int,default=40)
args=parser.parse_args()

if not INPUT.exists():
    raise SystemExit(f"missing completed degree-3 result {INPUT}")

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

G=load_gram(FRAME)
C=G[:13,:13]
coupling=G[:13,13:]
tail=G[13:,13:]
H=tail-coupling.transpose()*C.inverse()*coupling
assert H.det()==237

def class_order(dual):
    o=ZZ(1)
    for v in dual:
        o=lcm(o,ZZ(QQ(v).denominator()))
    return o

def section_profile(z):
    z=vector(ZZ,z)
    if z==0:
        # Zero section is available for free as a degree-one point.
        return {
            "po":ZZ(0),"height":QQ(0),"corr":QQ(0),
            "order":ZZ(1),"explicit_span":True,
        }
    base=vector(ZZ,[0]*13+list(z))
    dual=vector(QQ,base*G[:,:13])*C.inverse()
    order=class_order(dual)
    corr={ZZ(1):QQ(0),ZZ(2):QQ(1),ZZ(4):QQ(13)/4}[order]
    raw=QQ(dual*C*dual)
    mod2=lambda x: QQ(x)-2*(QQ(x)/2).floor()
    assert mod2(raw)==mod2(corr)
    h=QQ(z*H*z)
    po=(h+corr-4)/2
    if po not in ZZ or po<0:
        return None
    return {
        "po":ZZ(po),"height":h,"corr":corr,"order":order,
        # G1=(1,0,0,0), G3=(0,0,1,0) are exact over QQ(U).
        "explicit_span":bool(z[1]==0 and z[3]==0),
    }

# Height ceiling safely implied by P.O <= max-pole:
# h + corr = 2*P.O+4, corr >= 0.
height_bound=ZZ(2*args.max_pole+4)
scale=lcm(v.denominator() for v in H.list())
IH=(scale*H).change_ring(ZZ)
qf=pari(IH).qfminim(scale*height_bound)

cand={}
zero=vector(ZZ,[0,0,0,0])
cand[tuple(zero)]=(zero,section_profile(zero))
for col in matrix(ZZ,qf[2]).columns():
    for sign in (1,-1):
        z=sign*vector(ZZ,col)
        prof=section_profile(z)
        if prof is None or prof["po"]>args.max_pole:
            continue
        cand[tuple(z)]=(z,prof)

sections=list(cand.values())
sections.sort(key=lambda t:(
    t[1]["po"],
    0 if t[1]["explicit_span"] else 1,
    t[1]["height"],
    sum(abs(int(x)) for x in t[0]),
    tuple(t[0]),
))

print(
    "D13CUBIC_SECTIONS|"
    f"max_pole={args.max_pole}|count={len(sections)}|"
    f"explicit_G1G3_span={sum(p['explicit_span'] for _,p in sections)}|"
    "status=PASS",
    flush=True,
)

# For every pair sum retain Pareto-useful pair(s).  Candidate set is small,
# but retaining only a few metrics keeps the triple search trivial.
pairs={}
for i,(a,pa) in enumerate(sections):
    for b,pb in sections[i:]:
        s=tuple(a+b)
        metric=(
            max(int(pa["po"]),int(pb["po"])),
            int(pa["po"]+pb["po"]),
            int(not pa["explicit_span"])+int(not pb["explicit_span"]),
            sum(abs(int(x)) for x in a)+sum(abs(int(x)) for x in b),
        )
        old=pairs.get(s)
        if old is None or metric<old[0]:
            pairs[s]=(metric,a,pa,b,pb)

def triple_for(target, explicit_first=False):
    best=None
    for c,pc in sections:
        rem=tuple(target-c)
        rec=pairs.get(rem)
        if rec is None:
            continue
        pm,a,pa,b,pb=rec
        maxpo=max(pm[0],int(pc["po"]))
        sumpo=pm[1]+int(pc["po"])
        unknown=pm[2]+int(not pc["explicit_span"])
        l1=pm[3]+sum(abs(int(x)) for x in c)
        metric=(
            (unknown,maxpo,sumpo,l1)
            if explicit_first
            else (maxpo,sumpo,unknown,l1)
        )
        item=(metric,(a,pa),(b,pb),(c,pc))
        if best is None or metric<best[0]:
            best=item
    return best

data=json.loads(INPUT.read_text())
hits=[
    r for r in data.get("neighbors",[])
    if int(r.get("child_mw_rank",0))>=5
]

ranked=[]
for r in hits:
    z=vector(ZZ,r["mw_projection"])
    cheap=triple_for(z,False)
    explicit=triple_for(z,True)
    if cheap is None:
        continue
    cm=cheap[0]
    em=explicit[0] if explicit is not None else None
    ranked.append((cm,em,r,cheap,explicit))

ranked.sort(key=lambda x:(
    x[0][0],x[0][1],x[0][2],
    -int(x[2]["child_mw_rank"]),
    int(x[2]["q"]),int(x[2]["orbit_index"])
))

print(
    "D13CUBIC_INPUT|"
    f"rank_growth_neighbors={len(hits)}|"
    f"decomposable={len(ranked)}|status=PASS",
    flush=True,
)

def fmt_sec(item):
    z,p=item
    return "{}@po{}{}".format(
        ",".join(map(str,z)),
        p["po"],
        "*" if p["explicit_span"] else "",
    )

for rank,(cm,em,r,cheap,explicit) in enumerate(ranked[:args.top]):
    parts=";".join(fmt_sec(x) for x in cheap[1:])
    eparts="NA" if explicit is None else ";".join(fmt_sec(x) for x in explicit[1:])
    print(
        "D13CUBIC_BEST|"
        f"rank={rank}|q={r['q']}|orbit={r['orbit_index']}|"
        f"child={r['child_ade']}|MW={r['child_mw_rank']}|"
        f"target_mw={','.join(map(str,r['mw_projection']))}|"
        f"maxPO={cm[0]}|sumPO={cm[1]}|unknown={cm[2]}|"
        f"sections={parts}|"
        f"explicit_first_unknown={'NA' if em is None else em[0]}|"
        f"explicit_first_maxPO={'NA' if em is None else em[1]}|"
        f"explicit_sections={eparts}|status=PASS",
        flush=True,
    )

# Separate shortlist: maximize direct QQ constructibility first.
eranked=[x for x in ranked if x[4] is not None]
eranked.sort(key=lambda x:(
    x[1][0],x[1][1],x[1][2],
    -int(x[2]["child_mw_rank"]),
    int(x[2]["q"]),int(x[2]["orbit_index"])
))
for rank,(cm,em,r,cheap,explicit) in enumerate(eranked[:20]):
    print(
        "D13CUBIC_EXPLICIT|"
        f"rank={rank}|q={r['q']}|orbit={r['orbit_index']}|"
        f"child={r['child_ade']}|MW={r['child_mw_rank']}|"
        f"unknown={em[0]}|maxPO={em[1]}|sumPO={em[2]}|"
        f"sections={';'.join(fmt_sec(x) for x in explicit[1:])}|status=PASS",
        flush=True,
    )

if ranked:
    b=ranked[0]
    r=b[2]
    print(
        "D13CUBIC_RESULT|"
        f"best=q{r['q']}:o{r['orbit_index']}:{r['child_ade']}:"
        f"MW{r['child_mw_rank']}:maxPO{b[0][0]}:sumPO{b[0][1]}:"
        f"unknown{b[0][2]}|status=PASS_CUBIC_SECTION_FRONTIER",
        flush=True,
    )
else:
    print(
        "D13CUBIC_RESULT|best=NONE|status=NO_THREE_SECTION_DECOMPOSITION",
        flush=True,
    )

out=SCOUT/"d13-degree3-section-decomposition-frontier.json"
def serial_item(x):
    cm,em,r,cheap,explicit=x
    return {
        "q":int(r["q"]),
        "orbit_index":int(r["orbit_index"]),
        "child_ade":r["child_ade"],
        "child_mw_rank":int(r["child_mw_rank"]),
        "target_mw":list(map(int,r["mw_projection"])),
        "cheap_metric":list(map(int,cm)),
        "cheap_sections":[
            {
                "mw":list(map(int,z)),
                "P_dot_O":int(p["po"]),
                "height":str(p["height"]),
                "explicit_G1G3_span":bool(p["explicit_span"]),
            }
            for z,p in cheap[1:]
        ],
        "explicit_first_metric":None if em is None else list(map(int,em)),
        "explicit_first_sections":None if explicit is None else [
            {
                "mw":list(map(int,z)),
                "P_dot_O":int(p["po"]),
                "height":str(p["height"]),
                "explicit_G1G3_span":bool(p["explicit_span"]),
            }
            for z,p in explicit[1:]
        ],
    }

out.write_text(json.dumps({
    "schema":"elkies-k3.h3-d13-degree3-section-decomposition-frontier.v1",
    "status":"PASS_D13_DEGREE3_SECTION_DECOMPOSITION_FRONTIER",
    "max_section_pole":args.max_pole,
    "rank_growth_neighbors":len(hits),
    "decomposable_neighbors":len(ranked),
    "ranked":[serial_item(x) for x in ranked],
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
