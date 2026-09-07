from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Local-Tate MW3 construction scaffold with component-hit constraints.")
ap.add_argument("--stage",choices=["p1","p12","all"],default="all")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
ap.add_argument("--show-first",type=int,default=0)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,11,79)):
    raise SystemExit("choose a good prime not 2,3,11,79")

names=[f"a{i}" for i in range(8)] + [f"b{i}" for i in range(12)]
names += ["lam","s0","s1","sl"]
names += [f"x1_{i}" for i in range(5)] + [f"y1_{i}" for i in range(7)]
if args.stage in ("p12","all"):
    names += ["r2"]+[f"x2_{i}" for i in range(7)]+[f"y2_{i}" for i in range(10)]
if args.stage=="all":
    names += ["r3"]+[f"x3_{i}" for i in range(7)]+[f"y3_{i}" for i in range(10)]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t")
t=Rt.gen()

lam=d["lam"]; s0=d["s0"]; s1=d["s1"]; sl=d["sl"]
A=sum(d[f"a{i}"]*t**i for i in range(8)) - 3*t**8
B=sum(d[f"b{i}"]*t**i for i in range(12)) + 2*t**12
Delta=-16*(4*A**3+27*B**2)

Ap=A.derivative(t); App=Ap.derivative(t)
Bp=B.derivative(t); Bpp=Bp.derivative(t)

eqs=[]; tags=[]

def add(tag,e):
    e=R(e)
    if e!=0:
        tags.append(tag); eqs.append(e)

# I11 at infinity: keep the already-small high coefficient equations.
for k in range(14,24):
    add(f"INF_D{k}",Delta[k])

# Local singular cubic at I3, t=0.
add("I3_A0",A(0)+3*s0**2)
add("I3_B0",B(0)-2*s0**3)
add("I3_first",Bp(0)+s0*Ap(0))
add("I3_second",Ap(0)**2-6*s0*(s0*App(0)+Bpp(0)))

# I2, t=1.
add("I2_1_A",A(1)+3*s1**2)
add("I2_1_B",B(1)-2*s1**3)
add("I2_1_first",Bp(1)+s1*Ap(1))

# I2, t=lambda.
add("I2_lam_A",A(lam)+3*sl**2)
add("I2_lam_B",B(lam)-2*sl**3)
add("I2_lam_first",Bp(lam)+sl*Ap(lam))

def coeff_eqs(poly,prefix):
    if poly==0:return
    for k in range(poly.degree()+1):
        add(f"{prefix}_{k}",poly[k])

# P1 polynomial section.
X1=sum(d[f"x1_{i}"]*t**i for i in range(5))
Y1=sum(d[f"y1_{i}"]*t**i for i in range(7))
coeff_eqs(Y1**2-X1**3-A*X1-B,"P1")

# P1 components (2,1,0,1): nonidentity at infinity,0,lambda; identity at 1.
add("P1_inf_x",d["x1_4"]-1)
add("P1_inf_y",d["y1_6"])
add("P1_0_x",X1(0)-s0)
add("P1_0_y",Y1(0))
add("P1_lam_x",X1(lam)-sl)
add("P1_lam_y",Y1(lam))

if args.stage in ("p12","all"):
    r2=d["r2"]; z2=t-r2
    X2=sum(d[f"x2_{i}"]*t**i for i in range(7))
    Y2=sum(d[f"y2_{i}"]*t**i for i in range(10))
    coeff_eqs(Y2**2-X2**3-A*X2*z2**4-B*z2**6,"P2")

    # P2=(6,2,1,1): nonidentity at all four reducible fibers.
    add("P2_inf_x",d["x2_6"]-1)
    add("P2_inf_y",d["y2_9"])
    add("P2_0_x",X2(0)-s0*(0-r2)**2)
    add("P2_0_y",Y2(0))
    add("P2_1_x",X2(1)-s1*(1-r2)**2)
    add("P2_1_y",Y2(1))
    add("P2_lam_x",X2(lam)-sl*(lam-r2)**2)
    add("P2_lam_y",Y2(lam))

if args.stage=="all":
    r3=d["r3"]; z3=t-r3
    X3=sum(d[f"x3_{i}"]*t**i for i in range(7))
    Y3=sum(d[f"y3_{i}"]*t**i for i in range(10))
    coeff_eqs(Y3**2-X3**3-A*X3*z3**4-B*z3**6,"P3")

    # P3=(10,2,0,1): nonidentity at infinity,0,lambda; identity at 1.
    add("P3_inf_x",d["x3_6"]-1)
    add("P3_inf_y",d["y3_9"])
    add("P3_0_x",X3(0)-s0*(0-r3)**2)
    add("P3_0_y",Y3(0))
    add("P3_lam_x",X3(lam)-sl*(lam-r3)**2)
    add("P3_lam_y",Y3(lam))

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"MW3LOCAL|stage={args.stage}|field={field}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)

groups={}
for tag,e in zip(tags,eqs):
    group=tag.split("_")[0]
    groups[group]=groups.get(group,0)+1
print("MW3LOCAL|groups="+",".join(f"{k}:{v}" for k,v in groups.items()),flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    if i<args.show_first:
        print(f"MW3LOCAL_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}|eq={e}",flush=True)
    elif i<35 or i>=len(eqs)-8:
        print(f"MW3LOCAL_META|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

print("MW3LOCAL|saturate=lam*(lam-1)*s0*s1*sl*(r2-fibers)*(r3-fibers) as applicable",flush=True)
print("MW3LOCAL|note=I11 component indices 2,6,10 are not yet separated beyond nonidentity; deeper infinity blowups are next.",flush=True)

if args.export:
    if args.p==0: raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"MW3LOCAL|export={out}",flush=True)
