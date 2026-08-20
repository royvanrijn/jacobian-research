from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Build the semistable A10+A2+A1^2 elliptic-K3 fiber scaffold.")
ap.add_argument("--show-equations",action="store_true")
ap.add_argument("--out",default="artifacts/local/elkies-k3/mw3-fiber-scaffold.txt")
args=ap.parse_args()

names=[f"a{i}" for i in range(8)] + [f"b{i}" for i in range(12)] + ["lam"]
R=PolynomialRing(QQ,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t")
t=Rt.gen()

lam=d["lam"]
A=sum(d[f"a{i}"]*t**i for i in range(8)) - 3*t**8
B=sum(d[f"b{i}"]*t**i for i in range(12)) + 2*t**12
Delta=-16*(4*A**3+27*B**2)

eqs=[]
tags=[]

for k in range(14,24):
    e=R(Delta[k])
    if e:
        eqs.append(e)
        tags.append(f"inf_Delta_t{k}")

D1=Delta.derivative(t)
D2=D1.derivative(t)
for tag,e in [
    ("I3_0_D",Delta(0)),
    ("I3_0_D1",D1(0)),
    ("I3_0_D2",D2(0)),
    ("I2_1_D",Delta(1)),
    ("I2_1_D1",D1(1)),
    ("I2_lam_D",Delta(lam)),
    ("I2_lam_D1",D1(lam)),
]:
    eqs.append(R(e))
    tags.append(tag)

print(f"MW3BUILD|stage=fibers|vars={R.ngens()}|eqs={len(eqs)}|expected_dim={R.ngens()-len(eqs)}",flush=True)
print(f"MW3BUILD|A={A}",flush=True)
print(f"MW3BUILD|B={B}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"MW3BUILD_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)
    if args.show_equations:
        print(f"MW3BUILD_FORMULA|i={i}|eq={e}",flush=True)

D3=D2.derivative(t)
sat={
    "lambda":lam,
    "lambda_minus_1":lam-1,
    "c4_at_0":R(A(0)),
    "c4_at_1":R(A(1)),
    "c4_at_lambda":R(A(lam)),
    "Delta_t13":R(Delta[13]),
    "I3_exact":R(D3(0)),
    "I2_1_exact":R(D2(1)),
    "I2_lambda_exact":R(D2(lam)),
}
print("MW3BUILD|saturation="+",".join(sat.keys()),flush=True)

out=Path(args.out)
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write("vars="+repr(names)+"\n")
    h.write(f"equations={len(eqs)}\n")
    h.write("A="+str(A)+"\n")
    h.write("B="+str(B)+"\n")
    for tag,e in zip(tags,eqs):
        h.write(tag+" = "+str(e)+"\n")
    h.write("saturation:\n")
    for k,v in sat.items():
        h.write(f"{k} = {v}\n")
print(f"MW3BUILD|stage=done|out={out}",flush=True)
