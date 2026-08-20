from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Build fiber + three-section scaffold for the MW3 K3 node.")
ap.add_argument("--stage",choices=["p1","p12","all"],default="all")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
ap.add_argument("--show-first",type=int,default=0)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,11,79)):
    raise SystemExit("choose a good prime not 2,3,11,79")

names=[f"a{i}" for i in range(8)] + [f"b{i}" for i in range(12)] + ["lam"]
names += [f"x1_{i}" for i in range(5)] + [f"y1_{i}" for i in range(7)]

if args.stage in ("p12","all"):
    names += ["r2"]+[f"x2_{i}" for i in range(7)]+[f"y2_{i}" for i in range(10)]
if args.stage=="all":
    names += ["r3"]+[f"x3_{i}" for i in range(7)]+[f"y3_{i}" for i in range(10)]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t")
t=Rt.gen()

lam=d["lam"]
A=sum(d[f"a{i}"]*t**i for i in range(8)) - 3*t**8
B=sum(d[f"b{i}"]*t**i for i in range(12)) + 2*t**12
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t)
D2=D1.derivative(t)

eqs=[]
tags=[]

for k in range(14,24):
    e=R(Delta[k])
    if e:
        eqs.append(e)
        tags.append(f"F_inf_{k}")

for tag,e in [
    ("F_I3_0",Delta(0)),
    ("F_I3_0_d1",D1(0)),
    ("F_I3_0_d2",D2(0)),
    ("F_I2_1",Delta(1)),
    ("F_I2_1_d1",D1(1)),
    ("F_I2_lam",Delta(lam)),
    ("F_I2_lam_d1",D1(lam)),
]:
    eqs.append(R(e))
    tags.append(tag)

def coeff_equations(poly,prefix):
    out=[]
    if poly==0:
        return out
    for k in range(poly.degree()+1):
        e=R(poly[k])
        if e:
            out.append((f"{prefix}_{k}",e))
    return out

X1=sum(d[f"x1_{i}"]*t**i for i in range(5))
Y1=sum(d[f"y1_{i}"]*t**i for i in range(7))
S1=Y1**2-X1**3-A*X1-B
for tag,e in coeff_equations(S1,"P1"):
    tags.append(tag)
    eqs.append(e)

if args.stage in ("p12","all"):
    z2=t-d["r2"]
    X2=sum(d[f"x2_{i}"]*t**i for i in range(7))
    Y2=sum(d[f"y2_{i}"]*t**i for i in range(10))
    S2=Y2**2-X2**3-A*X2*z2**4-B*z2**6
    for tag,e in coeff_equations(S2,"P2"):
        tags.append(tag)
        eqs.append(e)

if args.stage=="all":
    z3=t-d["r3"]
    X3=sum(d[f"x3_{i}"]*t**i for i in range(7))
    Y3=sum(d[f"y3_{i}"]*t**i for i in range(10))
    S3=Y3**2-X3**3-A*X3*z3**4-B*z3**6
    for tag,e in coeff_equations(S3,"P3"):
        tags.append(tag)
        eqs.append(e)

field_name="QQ" if args.p==0 else f"GF({args.p})"
print(f"MW3SECT|stage={args.stage}|field={field_name}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    if i<args.show_first:
        print(f"MW3SECT_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}|eq={e}",flush=True)
    elif i<30 or i>=len(eqs)-5:
        print(f"MW3SECT_META|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

print("MW3SECT|component_target=P1:(2,1,0,1);P2:(6,2,1,1);P3:(10,2,0,1)",flush=True)
print("MW3SECT|intersection_target=P1O:0,P2O:1,P3O:1,P12:1,P13:2,P23:2",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"MW3SECT|export={out}",flush=True)
