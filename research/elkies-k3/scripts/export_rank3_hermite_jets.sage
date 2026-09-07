from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Hermite-jet formulation for the E8+A2^3 rank-jump section.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--slice",action="store_true")
ap.add_argument("--out",required=True)
a=ap.parse_args()

F=GF(a.p)
random.seed(a.seed); set_random_seed(a.seed)

# Variables:
# lam plus local cusp parameters (u,v) at t=0,1,lam, plus one free
# degree-6 y coefficient c.
names=["lam","u0","v0","u1","v1","ul","vl","c"]
R=PolynomialRing(F,names,order="degrevlex")
d=R.gens_dict()
K=FractionField(R)
Kt=PolynomialRing(K,"t"); t=Kt.gen()

lam=K(d["lam"]); c=K(d["c"])
u0,v0,u1,v1,ul,vl=[K(d[n]) for n in ["u0","v0","u1","v1","ul","vl"]]

nodes=[K(0),K(1),lam]
us=[u0,u1,ul]
vs=[v0,v1,vl]

# Hermite interpolation polynomial of degree <=5 through values and derivatives.
# We construct via CRT modulo product (t-a)^2.
P=t*(t-1)*(t-lam)
M=P^2

def hermite(values, derivs):
    # Solve six linear equations for degree <=5 coefficients.
    coeffs=vector(K,6,[K(0)]*6)
    A=[]; b=[]
    for aa,val,der in zip(nodes,values,derivs):
        A.append([aa^j for j in range(6)])
        b.append(val)
        A.append([K(0) if j==0 else j*aa^(j-1) for j in range(6)])
        b.append(der)
    A=matrix(K,A)
    sol=A.solve_right(vector(K,b))
    return sum(sol[j]*t^j for j in range(6))

X5=hermite([u^2 for u in us],[2*u*v for u,v in zip(us,vs)])
Y5=hermite([u^3 for u in us],[3*u^2*v for u,v in zip(us,vs)])

# x must have degree <=4: coefficient t^5 = 0.
# y can have degree <=6; adding c*M (degree 6) preserves all 3 double jets.
X=X5
Y=Y5+c*M

raw_eqs=[]
raw_eqs.append(X[5])

# D is automatically divisible by M by construction.
D=Y^2-X^3
Q,Rem=D.quo_rem(M)
assert Rem==0

# Need Q=t-mu. Eliminate mu: coeff t^1 = 1; all coeff t^k for k>=2 vanish.
raw_eqs.append(Q[1]-1)
for k in range(2,Q.degree()+1):
    raw_eqs.append(Q[k])

# Convert rational functions to polynomial numerators and primitive scalar normalization.
eqs=[]
for e in raw_eqs:
    num=R(e.numerator())
    if num==0: continue
    # divide by scalar leading coefficient only
    lc=num.leading_coefficient()
    if lc in F and lc!=0:
        num=R(num/lc)
    eqs.append(num)

print(f"HJET|stage=build|p={a.p}|vars={R.ngens()}|eqs={len(eqs)}|Qdeg={Q.degree()}",flush=True)
for i,e in enumerate(eqs):
    print(f"HJET_EQ|i={i}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

if a.slice:
    form=R(F.random_element())
    # one generic slice for an expected one-dimensional component
    for v in random.sample(list(R.gens()),min(5,R.ngens())):
        cc=F.random_element()
        while cc==0: cc=F.random_element()
        form += R(cc)*v
    eqs.append(form)

out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(names)+"\n"+str(a.p)+"\n")
    for i,e in enumerate(eqs):
        h.write(str(e).replace("**","^"))
        h.write(",\n" if i+1<len(eqs) else "\n")

print(f"HJET|stage=export|vars={R.ngens()}|eqs={len(eqs)}|slice={a.slice}|out={out}",flush=True)
print("HJET|generic_chart=lam*(lam-1) != 0",flush=True)
