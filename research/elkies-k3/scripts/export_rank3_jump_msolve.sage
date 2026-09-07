from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Export the historical parity-obstructed all-IV E8+A2^3 section system.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--slices",type=int,default=1)
ap.add_argument("--out",required=True)
a=ap.parse_args()
if not is_prime(a.p) or a.p in (2,3,79): raise SystemExit("choose good prime !=2,3,79")
F=GF(a.p); random.seed(a.seed); set_random_seed(a.seed)

# Family:
#   E_{lam,mu}: y^2 = x^3 + [t(t-1)(t-lam)]^2 (t-mu)
#
# Fibers generically: II* at infinity, IV at t=0,1,lam, II at t=mu.
#
# Polynomial section parameterization:
#   x=q^2+r, y=q^3+s
# with deg q=2, deg r,s<=1. Then y^2-x^3 automatically has degree <=7.
names=["lam","mu","q0","q1","q2","r0","r1","s0","s1"]
R=PolynomialRing(F,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

lam,mu=d["lam"],d["mu"]
q=d["q0"]+d["q1"]*t+d["q2"]*t^2
r=d["r0"]+d["r1"]*t
s=d["s0"]+d["s1"]*t
x=q^2+r
y=q^3+s
P=t*(t-1)*(t-lam)
g=P^2*(t-mu)

Fpoly=y^2-x^3-g
eqs=[R(Fpoly[k]) for k in range(8)]
assert all(e!=0 for e in eqs)

# Generic affine slices only to get zero-dimensional samples from the expected curve.
vars=list(R.gens())
for _ in range(a.slices):
    form=R(F.random_element())
    for v in random.sample(vars,min(5,len(vars))):
        c=F.random_element()
        while c==0: c=F.random_element()
        form += R(c)*v
    eqs.append(form)

out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(names)+"\n"+str(a.p)+"\n")
    for i,e in enumerate(eqs):
        h.write(str(e).replace("**","^"))
        h.write(",\n" if i+1<len(eqs) else "\n")

print(f"R3JUMP|stage=export|p={a.p}|seed={a.seed}|vars={len(names)}|base_eqs=8|slices={a.slices}|eqs={len(eqs)}|out={out}",flush=True)
