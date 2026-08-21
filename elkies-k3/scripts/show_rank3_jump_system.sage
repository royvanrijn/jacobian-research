from sage.all import *
from pathlib import Path

print("R3SYS|status=HISTORICAL_OBSTRUCTED_ALL_IV_CHART")

R=PolynomialRing(QQ,["lam","mu","q0","q1","q2","r0","r1","s0","s1"],order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()
lam,mu=d["lam"],d["mu"]
q=d["q0"]+d["q1"]*t+d["q2"]*t^2
r=d["r0"]+d["r1"]*t
s=d["s0"]+d["s1"]*t
P=t*(t-1)*(t-lam)
g=P^2*(t-mu)
F=(q^3+s)^2-(q^2+r)^3-g
eqs=[R(F[k]) for k in range(8)]

print("R3SYS|vars=9|eqs=8")
for k,e in enumerate(eqs):
    print(f"R3EQ|k={k}|terms={len(e.monomials())}|degree={e.degree()}|eq={e}")

# CM endpoint sanity: lam=mu=0. We do not require a generic MW section to remain
# nontrivial there; the rank disappears into the enhanced II* fiber.
cm=[e.subs({lam:0,mu:0}) for e in eqs]
print(f"R3SYS|cm_specialized_nonzero={sum(e!=0 for e in cm)}")
