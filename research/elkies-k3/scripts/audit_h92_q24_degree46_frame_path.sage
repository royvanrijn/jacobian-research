#!/usr/bin/env sage -python
"""
Locate the exact frame path for the stored raw-q6 degree-46 q24 bridge.

We deliberately test both q8-target interpretations:
  A. stored target is already in the Weyl-transported q6 chamber;
  B. stored target must first be taken through the 22 q6 Weyl reflections.

For the bridge P_raw itself there is no ambiguity: it comes from
search_h92_q6_mw_for_q24_bridge_raw_v2.sage and must satisfy P_raw.F6_raw=1.

After identifying the unique exact composition that sends the physical q8
fibre to F8_eq, apply the same composition to P_raw and certify degree 46 on
F8_eq.  Only then read q6 MW coordinates.
"""

import argparse
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, matrix, vector


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
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
    raise SystemExit("Could not locate repo")


def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
CERT=ROOT/"elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
CHEAP=LOCAL/"q24-cheapest-bridge-current.json"
TARGET=LOCAL/"q8-target-component-nef.json"
OUT=args.output.resolve() if args.output else LOCAL/"q24-degree46-frame-path.json"

for p in (CERT,CHEAP,TARGET):
    if not p.exists():
        raise SystemExit(f"missing {p}")

eq=run_scope(CERT)
cheap=json.loads(CHEAP.read_text())
target=json.loads(TARGET.read_text())

ns=eq["ns"]
F6eq=vector(ZZ,eq["F6"])
F8eq=vector(ZZ,eq["F8eq"])
Oold=vector(ZZ,eq["Oold"])
Ostd=vector(ZZ,eq["Ostd"])
rawF=vector(ZZ,eq["raw_fiber"])
weyl=eq["weyl_transport"]

sel=cheap["selected"]
Praw=vector(ZZ,sel["ambient_ns_vector"])
F8stored=vector(ZZ,target["selected_q8"]["source_h3_ns_vector"])

assert Praw*ns*Praw==-2
assert Praw*ns*rawF==1
assert Praw*ns*F8stored==46
assert rawF*ns*rawF==0
assert weyl(rawF)==F6eq

Pweyl=weyl(Praw)
F8weyl=weyl(F8stored)

assert Pweyl*ns*Pweyl==-2
assert Pweyl*ns*F6eq==1

print(
    "Q24FRAME_RAW|"
    f"Praw_Fraw={Praw*ns*rawF}|Praw_F8stored={Praw*ns*F8stored}|"
    f"F8stored_Fraw={F8stored*ns*rawF}|"
    f"F8stored_F6eq={F8stored*ns*F6eq}|"
    f"F8weyl_F6eq={F8weyl*ns*F6eq}|"
    "status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Exact Eichler translation old zero -> standard zero in the equation q6 frame.
# ---------------------------------------------------------------------------

ov=ZZ(Ostd*ns*Oold)
v=Ostd-Oold-(ov+2)*F6eq
assert v*ns*F6eq==0
assert v*ns*v==-12

def tau(x):
    x=vector(ZZ,x)
    return (
        x
        +(x*ns*F6eq)*v
        -(x*ns*v)*F6eq
        -(v*ns*v//2)*(x*ns*F6eq)*F6eq
    )

assert tau(F6eq)==F6eq
assert tau(Oold)==Ostd

# Four exact candidates for the target path.  We demand exact equality with
# F8eq, not just Gram/intersection agreement.
target_candidates=[
    ("stored",F8stored),
    ("weyl",F8weyl),
    ("tau_stored",tau(F8stored)),
    ("tau_weyl",tau(F8weyl)),
]

for name,F in target_candidates:
    print(
        "Q24FRAME_TARGET|"
        f"name={name}|square={F*ns*F}|q6_degree={F*ns*F6eq}|"
        f"equals_eq={int(F==F8eq)}|status=DIAGNOSTIC",
        flush=True,
    )

paths=[]
# identity, Weyl, tau, tau∘Weyl applied to raw objects.
def identity(x): return vector(ZZ,x)
compositions=[
    ("identity",identity),
    ("weyl",weyl),
    ("tau",tau),
    ("tau_after_weyl",lambda x: tau(weyl(x))),
]

for name,fun in compositions:
    try:
        Ft=fun(F8stored)
        Pt=fun(Praw)
    except Exception as exc:
        print(
            f"Q24FRAME_PATH|name={name}|status=ERROR|type={type(exc).__name__}|message={exc}",
            flush=True,
        )
        continue
    fibre_match=(Ft==F8eq)
    p_is_section=(Pt*ns*F6eq==1 and Pt*ns*Pt==-2)
    degree=Pt*ns*F8eq
    print(
        "Q24FRAME_PATH|"
        f"name={name}|fibre_match={int(fibre_match)}|"
        f"P_section_q6eq={int(p_is_section)}|degree_on_eq={degree}|"
        "status=DIAGNOSTIC",
        flush=True,
    )
    if fibre_match and p_is_section and degree==46:
        paths.append((name,fun,Pt))

if len(paths)!=1:
    # If target artifact itself is not the right fibre representative for
    # selecting the path, retain every composition that at least sends the
    # raw bridge to a q6 section of equation degree 46. This is diagnostic,
    # not authorization to choose silently.
    bridge_only=[]
    for name,fun in compositions:
        try:
            Pt=fun(Praw)
        except Exception:
            continue
        if Pt*ns*Pt==-2 and Pt*ns*F6eq==1 and Pt*ns*F8eq==46:
            bridge_only.append((name,fun,Pt))
    print(
        "Q24FRAME_SELECTION|"
        f"exact_paths={','.join(x[0] for x in paths) or 'none'}|"
        f"bridge_only={','.join(x[0] for x in bridge_only) or 'none'}|"
        "status=NEEDS_INTERPRETATION",
        flush=True,
    )
    chosen = paths[0] if len(paths)==1 else (bridge_only[0] if len(bridge_only)==1 else None)
else:
    chosen=paths[0]

payload={
    "schema":"elkies-k3.h92-q24-degree46-frame-path.v1",
    "status":"PASS_DIAGNOSTIC_FRAME_PATH",
    "raw":{
        "P_Fraw":int(Praw*ns*rawF),
        "P_F8stored":int(Praw*ns*F8stored),
        "F8stored_Fraw":int(F8stored*ns*rawF),
    },
    "target_candidates":{
        name:{
            "equals_equation":bool(F==F8eq),
            "q6_degree":int(F*ns*F6eq),
        }
        for name,F in target_candidates
    },
}

if chosen is None:
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print("Q24FRAME_RESULT|status=NO_UNIQUE_FRAME_PATH",flush=True)
    raise SystemExit(2)

name,fun,Peq=chosen

# ---------------------------------------------------------------------------
# Direct old-zero MW coordinates in equation q6 frame.
# ---------------------------------------------------------------------------

phis=[vector(QQ,p) for p in eq["phis"]]
H=matrix(QQ,eq["EXPECTED_HEIGHT"])
projection=eq["projection"]

horizontal=Peq-Oold-(Peq*ns*Oold+2)*F6eq
assert horizontal*ns*F6eq==0 and horizontal*ns*Oold==0
phiP=vector(QQ,horizontal)*projection
pair=vector(QQ,[-phiP*ns*q for q in phis])
mw_old_q=pair*H.inverse()
assert all(x in ZZ for x in mw_old_q)
mw_old=vector(ZZ,mw_old_q)

# Standard zero has old-group coordinate (-2,1,0).
mw_std=mw_old-vector(ZZ,(-2,1,0))

payload.update({
    "status":"PASS_UNIQUE_DEGREE46_FRAME_PATH",
    "selected_path":name,
    "equation_bridge_ns":list(map(int,Peq)),
    "equation_degree":int(Peq*ns*F8eq),
    "q6_old_mw":list(map(int,mw_old)),
    "q6_standard_mw":list(map(int,mw_std)),
    "P_dot_Oold":int(Peq*ns*Oold),
    "P_dot_Ostd":int(Peq*ns*Ostd),
})

print(
    "Q24FRAME_MW|"
    f"path={name}|old_mw={','.join(map(str,mw_old))}|"
    f"standard_mw={','.join(map(str,mw_std))}|"
    f"PdotOold={Peq*ns*Oold}|PdotOstd={Peq*ns*Ostd}|"
    "status=PASS_DIRECT_SHIODA",
    flush=True,
)

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24FRAME_RESULT|"
    f"path={name}|equation_degree=46|"
    f"standard_mw={','.join(map(str,mw_std))}|"
    "status=PASS_UNIQUE_DEGREE46_FRAME_PATH",
    flush=True,
)
