from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6+A3^2+A1^2 MW3 fiber scaffold.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--show",action="store_true")
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("bad prime")

# Preferred fibers:
# IV* at infinity, I4 at 0, I4 at 1, I2 at lambda, I2 at mu.
#
# IV* at infinity in short Weierstrass:
# ord_inf A >=3 -> deg A<=5
# ord_inf B =4 -> deg B=8, B8 !=0.
# Normalize B8=1 by Weierstrass scaling on this construction chart.
names=[f"a{i}" for i in range(6)] + [f"b{i}" for i in range(8)]
names += ["lam","mu","s0","s1","sl","sm"]
R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

A=sum(d[f"a{i}"]*t**i for i in range(6))
B=sum(d[f"b{i}"]*t**i for i in range(8)) + t**8
Ap=A.derivative(t); App=Ap.derivative(t); A3=App.derivative(t)
Bp=B.derivative(t); Bpp=Bp.derivative(t); B3=Bpp.derivative(t)
Delta=-16*(4*A**3+27*B**2)

eqs=[]; tags=[]
def add(tag,e):
    e=R(e)
    if e!=0: tags.append(tag); eqs.append(e)

# Generic multiplicative local conditions.
def impose_I2(a,s,prefix):
    add(prefix+"_A",A(a)+3*s**2)
    add(prefix+"_B",B(a)-2*s**3)
    add(prefix+"_d1",Bp(a)+s*Ap(a))

# I4 needs Delta orders 0..3. First three local conditions via singular root,
# then impose Delta'''(a)=0 as the fourth vanishing condition.
def impose_I4(a,s,prefix):
    impose_I2(a,s,prefix)
    D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)
    add(prefix+"_d2",D2(a))
    add(prefix+"_d3",D3(a))

impose_I4(K(0),d["s0"],"I4_0")
impose_I4(K(1),d["s1"],"I4_1")
impose_I2(d["lam"],d["sl"],"I2_lam")
impose_I2(d["mu"],d["sm"],"I2_mu")

print(f"E6BUILD|stage=fibers|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)
print(f"E6BUILD|A={A}",flush=True)
print(f"E6BUILD|B={B}",flush=True)
for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6BUILD_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)
    if args.show:
        print(f"E6BUILD_FORMULA|i={i}|eq={e}",flush=True)

print("E6BUILD|saturate=B8*lam*(lam-1)*mu*(mu-1)*(lam-mu)*s0*s1*sl*sm and exact Kodaira residuals",flush=True)
