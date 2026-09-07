from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Use the three IV fibers directly: x=q^2, y=q^3 modulo P^2.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--slice",action="store_true")
ap.add_argument("--out",required=True)
a=ap.parse_args()

F=GF(a.p)
random.seed(a.seed); set_random_seed(a.seed)

# q is an element of F[t]/P^2, represented in degree < 6.
names=["lam"]+[f"q{i}" for i in range(6)]+["c"]
R=PolynomialRing(F,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

lam=d["lam"]; c=d["c"]
q=sum(d[f"q{i}"]*t^i for i in range(6))
P=t*(t-1)*(t-lam)
P2=P^2

# At each root of P, y^2=x^3 to second order.
# On the generic chart where x is nonzero at those roots, q=y/x gives
# x == q^2 (mod P^2), y == q^3 (mod P^2).
X=(q^2).mod(P2)
Y0=(q^3).mod(P2)

# x must have degree <=4, so the t^5 coefficient of X vanishes.
eqs=[R(X[5])]

# y has degree <=6, hence the only freedom beyond its remainder mod P^2
# is adding c*P^2.
Y=Y0+c*P2

# Compute quotient Q=(y^2-x^3)/P^2 exactly.
Dpoly=Y^2-X^3
quo,rem=Dpoly.quo_rem(P2)
for k in range(rem.degree()+1 if rem else 0):
    assert R(rem[k])==0

# Our family requires quotient exactly t-mu.
# Eliminate mu: constant term is free; coefficient of t must be 1,
# all coefficients t^2.. vanish.
eqs.append(R(quo[1]-1))
for k in range(2,quo.degree()+1):
    eqs.append(R(quo[k]))

eqs=[e for e in eqs if e!=0]

print(f"IVMOD|stage=build|p={a.p}|vars={R.ngens()}|eqs={len(eqs)}|quot_degree={quo.degree()}",flush=True)
for i,e in enumerate(eqs):
    print(f"IVMOD_EQ|i={i}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

if a.slice:
    f=R(F.random_element())
    for v in random.sample(list(R.gens()),min(4,R.ngens())):
        cc=F.random_element()
        while cc==0: cc=F.random_element()
        f += R(cc)*v
    eqs.append(f)

out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(names)+"\n"+str(a.p)+"\n")
    for i,e in enumerate(eqs):
        h.write(str(e).replace("**","^"))
        h.write(",\n" if i+1<len(eqs) else "\n")

print(f"IVMOD|stage=export|vars={R.ngens()}|eqs={len(eqs)}|slice={a.slice}|out={out}",flush=True)
print("IVMOD|mu_formula=-constant_term_of_quotient",flush=True)
