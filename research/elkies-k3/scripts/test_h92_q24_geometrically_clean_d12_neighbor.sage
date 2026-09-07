#!/usr/bin/env sage -python
"""
Test the unique geometrically clean degree-two D12 neighbour visible on the
CURRENT resolved I9* fibre.

The resolved affine I9* enumeration found exactly one degree pattern whose
twelve zero-degree old components form D12:
    D.E10a = D.E10b = 1,
    D.C = 0 on every other old fibre component.

In the finite D13 root graph E10a,E10b are the two spinor leaves.  We do not
need to know the geometric Weyl chamber inside NS to test the abstract class:
in the current deterministic D13 frame the same pattern is simply "Dynkin
labels 1 on both spinor leaves, 0 elsewhere".  This is invariant under the
only D13 graph automorphism (spinor swap).

Fix the SAME q24 horizontal MW class z=(-2,1,-1,1), solve the root coordinates
with those Dynkin labels, impose D^2=0, and map the class back to the exact
source NS.

Certify:
  * primitive isotropic, old q8-fibre degree 2;
  * D.O=14 (equivalently U coefficient a=16);
  * child root data and MW rank;
  * difference from transported orbit85 D24;
  * whether the difference is a single root reflection / short Weyl move.

If child=D12/MW5, this is an equation-friendly alternative neighbour whose
special D12 fibre is already the resolved old-component configuration.
"""

import argparse
import contextlib
import io
import sys
from pathlib import Path

from sage.all import QQ, ZZ, gcd, matrix, pari, vector


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
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    buf=io.StringIO()
    try:
        sys.argv=[str(path)]
        with contextlib.redirect_stdout(buf):
            exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


def graph_shape(R):
    n=R.nrows()
    graph={i:set() for i in range(n)}
    for i in range(n):
        assert R[i,i]==2
        for j in range(i+1,n):
            assert R[i,j] in (0,-1)
            if R[i,j]==-1:
                graph[i].add(j); graph[j].add(i)
    branch=[i for i in graph if len(graph[i])==3]
    assert len(branch)==1
    b=branch[0]
    spin=[j for j in graph[b] if len(graph[j])==1]
    assert len(spin)==2
    return graph,b,tuple(sorted(spin))


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
CLOSE=ROOT/"elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"
if not CLOSE.exists():
    raise SystemExit(f"Missing {CLOSE}")

cl=run_scope(CLOSE)

need=("ns","F8eq","O8","D24eq","Badapt","adapted","H13","isotropic_mate","roots_and_data")
missing=[x for x in need if x not in cl]
if missing:
    raise SystemExit("close scope missing: "+",".join(missing))

ns=cl["ns"]
F=vector(ZZ,cl["F8eq"])
O=vector(ZZ,cl["O8"])
D85=vector(ZZ,cl["D24eq"])
Badapt=cl["Badapt"]
G=cl["adapted"]
H=cl["H13"]
isotropic_mate=cl["isotropic_mate"]
roots_and_data=cl["roots_and_data"]

R=G[:13,:13]
graph,branch,spin=graph_shape(R)

z=vector(ZZ,(-2,1,-1,1))
assert z*H*z==52

labels=vector(ZZ,[ZZ(i in spin) for i in range(13)])
assert sum(labels)==2

# Solve [r,z] * G[:,root] = labels.
C=G[:13,13:]
# r*R + z*C^T = labels
rhs=vector(QQ,labels)-vector(QQ,z)*C.transpose()
r=rhs*R.inverse()
assert all(v in ZZ for v in r),r
r=vector(ZZ,r)

w=vector(ZZ,list(r)+list(z))
norm=ZZ(w*G*w)
assert norm%4==0,norm
a=norm//4
Dcoords=vector(ZZ,[a,2]+list(w))

Gns=cl["Gadapt"]
assert Dcoords*Gns*Dcoords==0
assert Dcoords[1]==2

Dgeo=vector(ZZ,Dcoords*Badapt)
assert Dgeo*ns*Dgeo==0
assert Dgeo*ns*F==2
assert gcd(tuple(Dgeo))==1

DdotO=ZZ(Dgeo*ns*O)
assert DdotO==a-2

# Pairings with deterministic effective simple system -e_i.
eff=[]
for i in range(13):
    c=vector(ZZ,[0,0]+[-ZZ(i==j) for j in range(17)])
    eff.append(c)
pair=[int(Dcoords*Gns*c) for c in eff]
assert pair==list(labels)

mate=isotropic_mate(ns,Dgeo)
orth=matrix(ZZ,[list(Dgeo*ns),list(mate*ns)]).right_kernel_matrix()
child=-(orth*ns*orth.transpose())
root_data=tuple(map(int,roots_and_data(child)[2]))
mw_rank=17-root_data[0]

delta=Dgeo-D85
delta_sq=ZZ(delta*ns*delta)
delta_F=ZZ(delta*ns*F)
delta_O=ZZ(delta*ns*O)

# Search whether one root reflection sends D85 to Dgeo.
roots,unused,unused_data=roots_and_data(G)
single=[]
for rr in roots:
    # Root curve in full adapted NS coordinates.
    rho=vector(ZZ,[0,0]+list(rr))
    # Reflection s_r(D)=D+(D.r)r in NS convention for rho^2=-2.
    image=D85 + ZZ(D85*ns*(rho*Badapt))*(rho*Badapt)
    if image==Dgeo:
        single.append(list(map(int,rr)))

# Breadth-2 search in D13 roots, bounded and exact.
short_weyl=None
root_source=[vector(ZZ,[0,0]+list(rr))*Badapt for rr in roots]
if not single:
    seen={}
    for idx,rho in enumerate(root_source):
        d1=D85+ZZ(D85*ns*rho)*rho
        seen[tuple(d1)]=idx
        if d1==Dgeo:
            short_weyl=[idx]
            break
    if short_weyl is None:
        for key,i in list(seen.items()):
            d1=vector(ZZ,key)
            for j,rho in enumerate(root_source):
                d2=d1+ZZ(d1*ns*rho)*rho
                if d2==Dgeo:
                    short_weyl=[i,j]
                    break
            if short_weyl is not None:
                break

print(
    "Q24GEOCLEAN_CLASS|"
    f"spinor_nodes={spin[0]+1},{spin[1]+1}|"
    f"labels={','.join(map(str,labels))}|"
    f"mw={','.join(map(str,z))}|root={','.join(map(str,r))}|"
    f"norm={norm}|a={a}|DdotO={DdotO}|"
    f"degree={Dgeo*ns*F}|primitive={int(gcd(tuple(Dgeo))==1)}|status=PASS",
    flush=True,
)

print(
    "Q24GEOCLEAN_CHILD|"
    f"root_data={root_data[0]},{root_data[1]},{root_data[2]}|"
    f"MW={mw_rank}|"
    f"child={'D12' if root_data==(12,264,4) else 'OTHER'}|status=PASS",
    flush=True,
)

print(
    "Q24GEOCLEAN_COMPARE|"
    f"delta_square={delta_sq}|delta_F={delta_F}|delta_O={delta_O}|"
    f"same={int(Dgeo==D85)}|"
    f"single_root_reflection={int(bool(single))}|"
    f"weyl_length_le2={-1 if short_weyl is None else len(short_weyl)}|"
    "status=PASS",
    flush=True,
)

print(
    "Q24GEOCLEAN_RESULT|"
    f"DdotO={DdotO}|child_root={root_data[0]},{root_data[1]},{root_data[2]}|"
    f"MW={mw_rank}|"
    f"status={'PASS_ALTERNATIVE_D12_MW5' if root_data==(12,264,4) else 'NOT_D12_MW5'}",
    flush=True,
)
