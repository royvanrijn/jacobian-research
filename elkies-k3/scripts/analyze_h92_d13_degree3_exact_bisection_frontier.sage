#!/usr/bin/env sage -python
"""
Find degree-3 D13 rank-growth neighbours whose horizontal divisor can be
realized as:

    exact old q8 bisection + cheap q8 section.

The exact old-curve audit stores for every curve C:
    degree d = C.F8
    z(C) = AJ(C-d O8) in MW(D13).

Thus for d(C)=2 and a section P:
    AJ(C + P - 3 O8) = z(C) + z(P).

This is a much weaker/smaller requirement than splitting the cubic divisor
into three rational sections.
"""

import argparse, json
from pathlib import Path
from sage.all import QQ, ZZ, lcm, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
SCOUT=LOCAL/"route-scout"
FRAME=ROOT/"elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
OLD=LOCAL/"q8-explicit-old-curves.json"
INPUT=SCOUT/"d13-degree3-qle40-rank-growth.json"

parser=argparse.ArgumentParser()
parser.add_argument("--max-pole",type=int,default=12)
parser.add_argument("--top",type=int,default=50)
args=parser.parse_args()

for p in (FRAME,OLD,INPUT):
    if not p.exists():
        raise SystemExit(f"missing {p}")

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
        return {
            "po":ZZ(0),"height":QQ(0),"corr":QQ(0),"order":ZZ(1),
            "explicit_G1G3":True,
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
        # Exact equation-level points currently available over QQ(U):
        # G1=(1,0,0,0), G3=(0,0,1,0).
        "explicit_G1G3":bool(z[1]==0 and z[3]==0),
    }

# Enumerate every section with P.O <= max-pole.
height_bound=ZZ(2*args.max_pole+4)
scale=lcm(v.denominator() for v in H.list())
IH=(scale*H).change_ring(ZZ)
qf=pari(IH).qfminim(scale*height_bound)

sections={}
z0=vector(ZZ,[0,0,0,0])
sections[tuple(z0)]=(z0,section_profile(z0))
for col in matrix(ZZ,qf[2]).columns():
    for sign in (1,-1):
        z=sign*vector(ZZ,col)
        prof=section_profile(z)
        if prof is None or prof["po"]>args.max_pole:
            continue
        sections[tuple(z)]=(z,prof)

old=json.loads(OLD.read_text())
assert old["status"]=="PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE"

exact_curves=[]
for rec in old["curves"]:
    d=int(rec["q8_degree"])
    z=vector(ZZ,rec["mw_coordinates"])
    exact_curves.append({
        "name":rec["curve"],
        "degree":d,
        "mw":z,
        "height":rec["height"],
        "source_h3_ns":rec.get("source_h3_ns"),
    })
    print(
        "D13CUBIC_OLDCURVE|"
        f"curve={rec['curve']}|degree={d}|"
        f"mw={','.join(map(str,z))}|height={rec['height']}|status=PASS",
        flush=True,
    )

bisections=[c for c in exact_curves if c["degree"]==2]
direct_cubics=[c for c in exact_curves if c["degree"]==3]

print(
    "D13CUBIC_BISECTION_INPUT|"
    f"old_curves={len(exact_curves)}|bisections={len(bisections)}|"
    f"direct_cubics={len(direct_cubics)}|sections_le_po{args.max_pole}={len(sections)}|"
    "status=PASS",
    flush=True,
)

data=json.loads(INPUT.read_text())
neighbors=[
    r for r in data.get("neighbors",[])
    if int(r.get("child_mw_rank",0))>=5
]

hits=[]
for r in neighbors:
    target=vector(ZZ,r["mw_projection"])

    # A single already-explicit cubic multisection is the ideal case.
    for c in direct_cubics:
        if c["mw"]==target:
            hits.append({
                "mode":"exact_cubic",
                "section_po":-1,
                "section_explicit_G1G3":True,
                "curve":c,
                "section":None,
                "neighbor":r,
            })

    # Exact bisection + one section.
    for c in bisections:
        rem=target-c["mw"]
        item=sections.get(tuple(rem))
        if item is None:
            continue
        z,p=item
        hits.append({
            "mode":"bisection_plus_section",
            "section_po":int(p["po"]),
            "section_explicit_G1G3":bool(p["explicit_G1G3"]),
            "curve":c,
            "section":{
                "mw":list(map(int,z)),
                "P_dot_O":int(p["po"]),
                "height":str(p["height"]),
                "correction":str(p["corr"]),
                "explicit_G1G3":bool(p["explicit_G1G3"]),
            },
            "neighbor":r,
        })

def key(h):
    r=h["neighbor"]
    return (
        0 if h["mode"]=="exact_cubic" else 1,
        h["section_po"],
        0 if h["section_explicit_G1G3"] else 1,
        -int(r["child_mw_rank"]),
        int(r["q"]),
        int(r["orbit_index"]),
        h["curve"]["name"],
    )

hits.sort(key=key)

print(
    "D13CUBIC_BISECTION_GROUP|"
    f"rank_growth_neighbors={len(neighbors)}|hits={len(hits)}|"
    f"qq_direct={sum(h['mode']=='exact_cubic' or h['section_explicit_G1G3'] for h in hits)}|"
    "status=PASS",
    flush=True,
)

for rank,h in enumerate(hits[:args.top]):
    r=h["neighbor"]
    s=h["section"]
    print(
        "D13CUBIC_BISECTION_HIT|"
        f"rank={rank}|mode={h['mode']}|q={r['q']}|orbit={r['orbit_index']}|"
        f"child={r['child_ade']}|MW={r['child_mw_rank']}|"
        f"target_mw={','.join(map(str,r['mw_projection']))}|"
        f"curve={h['curve']['name']}|curve_mw={','.join(map(str,h['curve']['mw']))}|"
        f"section_mw={'NA' if s is None else ','.join(map(str,s['mw']))}|"
        f"section_PO={'NA' if s is None else s['P_dot_O']}|"
        f"section_G1G3={int(h['section_explicit_G1G3'])}|"
        "status=PASS",
        flush=True,
    )

if hits:
    b=hits[0]
    r=b["neighbor"]
    print(
        "D13CUBIC_BISECTION_RESULT|"
        f"best={b['mode']}:q{r['q']}:o{r['orbit_index']}:"
        f"{r['child_ade']}:MW{r['child_mw_rank']}:"
        f"curve={b['curve']['name']}:"
        f"sectionPO={'NA' if b['section'] is None else b['section']['P_dot_O']}:"
        f"G1G3={int(b['section_explicit_G1G3'])}|"
        "status=PASS_EXACT_MULTISECTION_FRONTIER",
        flush=True,
    )
else:
    print(
        "D13CUBIC_BISECTION_RESULT|best=NONE|"
        "status=NO_EXACT_BISECTION_PLUS_CHEAP_SECTION",
        flush=True,
    )

out=SCOUT/"d13-degree3-exact-bisection-frontier.json"
def ser(h):
    r=h["neighbor"]
    return {
        "mode":h["mode"],
        "q":int(r["q"]),
        "orbit_index":int(r["orbit_index"]),
        "child_ade":r["child_ade"],
        "child_mw_rank":int(r["child_mw_rank"]),
        "target_mw":list(map(int,r["mw_projection"])),
        "curve":{
            "name":h["curve"]["name"],
            "degree":h["curve"]["degree"],
            "mw":list(map(int,h["curve"]["mw"])),
            "height":h["curve"]["height"],
            "source_h3_ns":h["curve"]["source_h3_ns"],
        },
        "section":h["section"],
    }

out.write_text(json.dumps({
    "schema":"elkies-k3.h3-d13-degree3-exact-bisection-frontier.v1",
    "status":"PASS_D13_DEGREE3_EXACT_BISECTION_FRONTIER",
    "max_section_pole":args.max_pole,
    "rank_growth_neighbors":len(neighbors),
    "exact_old_curves":[
        {
            "name":c["name"],"degree":c["degree"],
            "mw":list(map(int,c["mw"])),"height":c["height"],
        } for c in exact_curves
    ],
    "hits":[ser(h) for h in hits],
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
