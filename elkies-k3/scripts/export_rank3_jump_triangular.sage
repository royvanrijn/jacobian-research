from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Triangularly eliminate the historical parity-obstructed all-IV section system.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--slice",action="store_true")
ap.add_argument("--out",required=True)
a=ap.parse_args()

if not is_prime(a.p) or a.p in (2,3,79): raise SystemExit("bad prime")
F=GF(a.p); random.seed(a.seed); set_random_seed(a.seed)

# Remaining primary variables. lam is retained initially because k5 becomes quadratic.
names=["lam","q0","q1","q2","r0","r1","s0"]
R=PolynomialRing(F,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

lam=d["lam"]; q0=d["q0"]; q1=d["q1"]; q2=d["q2"]
r0=d["r0"]; r1=d["r1"]; s0=d["s0"]

# Generic chart q2 != 0. k7 solves s1.
# k7:
# -12*q1*q2^3*r0 -18*q1^2*q2^2*r1 -12*q0*q2^3*r1 +2*q2^3*s1 -1 = 0
s1 = (
    12*q1*q2^3*r0
    +18*q1^2*q2^2*r1
    +12*q0*q2^3*r1
    +1
) / (2*q2^3)

# k6 solves mu.
# -18 q1^2 q2^2 r0 -12 q0 q2^3 r0 -12 q1^3 q2 r1
# -36 q0 q1 q2^2 r1 -3 q2^2 r1^2 +2 q2^3 s0
# +6 q1 q2^2 s1 +2 lam + mu +2 =0
mu = (
    18*q1^2*q2^2*r0
    +12*q0*q2^3*r0
    +12*q1^3*q2*r1
    +36*q0*q1*q2^2*r1
    +3*q2^2*r1^2
    -2*q2^3*s0
    -6*q1*q2^2*s1
    -2*lam
    -2
)

q=q0+q1*t+q2*t^2
r=r0+r1*t
s=s0+s1*t
x=q^2+r
y=q^3+s
P=t*(t-1)*(t-lam)
g=P^2*(t-mu)
Fpoly=(y^2-x^3-g)

# k7,k6 vanish by construction. Keep k0..k5.
eqs=[]
for k in range(6):
    e=Fpoly[k]
    # Clear denominators by taking numerator.
    e=R(e.numerator())
    # Remove repeated content/scalar.
    if e!=0:
        eqs.append(e)

print(f"R3TRI|stage=build|p={a.p}|vars={R.ngens()}|eqs={len(eqs)}")
for i,e in enumerate(eqs):
    print(f"R3TRI_EQ|i={i}|degree={e.degree()}|terms={len(e.monomials())}")

# One equation should encode the remaining quadratic dependence on lambda.
# We expect a one-dimensional section locus before a slice.
if a.slice:
    form=R(F.random_element())
    # Avoid q2 so generic-chart localization isn't accidentally forced to q2=0.
    for v in random.sample([d[n] for n in names if n!="q2"],4):
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

print(f"R3TRI|stage=export|vars={R.ngens()}|eqs={len(eqs)}|slice={a.slice}|out={out}")
print("R3TRI|chart=q2_nonzero")
