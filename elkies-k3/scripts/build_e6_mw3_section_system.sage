from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6 MW3 local-Tate section construction.")
ap.add_argument("--stage",choices=["p1","p13","all"],default="p13")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
ap.add_argument("--show-first",type=int,default=0)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose a good prime not 2,3,79")

# Fibers:
#   IV* @ infinity
#   I4  @ 0
#   I4  @ 1
#   I2  @ lam
#   I2  @ mu
#
# Canonical MW component labels:
# P1=(1,0,1,0,0), O-intersection 0
# P2=(1,1,2,1,1), O-intersection 1
# P3=(2,3,0,0,0), O-intersection 0
#
# Tuple order = (IV*, I4@0, I4@1, I2@lam, I2@mu)

names=[f"a{i}" for i in range(6)] + [f"b{i}" for i in range(8)]
names += ["lam","mu","s0","s1","sl","sm"]

# P1 polynomial x deg <=4, y deg <=6
names += [f"x1_{i}" for i in range(5)] + [f"y1_{i}" for i in range(7)]

if args.stage in ("p13","all"):
    # P3 is also polynomial.
    names += [f"x3_{i}" for i in range(5)] + [f"y3_{i}" for i in range(7)]

if args.stage=="all":
    # P2.O=1 -> one linear denominator.
    # c12 is the required simple good-fiber intersection P2=-P1, equivalently
    # (P1+P2).O=1 and <P1,P2>=-5/6 for the certified component profiles.
    names += ["r2", "c12"]
    names += [f"x2_{i}" for i in range(7)]+[f"y2_{i}" for i in range(10)]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

lam=d["lam"]; mu=d["mu"]
s0=d["s0"]; s1=d["s1"]; sl=d["sl"]; sm=d["sm"]

# IV* at infinity: deg A<=5, deg B=8. Normalize B8=1.
A=sum(d[f"a{i}"]*t^i for i in range(6))
B=sum(d[f"b{i}"]*t^i for i in range(8)) + t^8
Ap=A.derivative(t); Bp=B.derivative(t)
Delta=-16*(4*A^3+27*B^2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

eqs=[]; tags=[]
def add(tag,e):
    e=R(e)
    if e!=0:
        tags.append(tag); eqs.append(e)

def impose_I2(a,s,prefix):
    add(prefix+"_A",A(a)+3*s^2)
    add(prefix+"_B",B(a)-2*s^3)
    add(prefix+"_d1",Bp(a)+s*Ap(a))

def impose_I4(a,s,prefix):
    impose_I2(a,s,prefix)
    add(prefix+"_D2",D2(a))
    add(prefix+"_D3",D3(a))

impose_I4(K(0),s0,"I4_0")
impose_I4(K(1),s1,"I4_1")
impose_I2(lam,sl,"I2_lam")
impose_I2(mu,sm,"I2_mu")

def coeff_eqs(poly,prefix):
    if poly==0: return
    for k in range(poly.degree()+1):
        add(f"{prefix}_{k}",poly[k])

# P1
X1=sum(d[f"x1_{i}"]*t^i for i in range(5))
Y1=sum(d[f"y1_{i}"]*t^i for i in range(7))
coeff_eqs(Y1^2-X1^3-A*X1-B,"P1")

# IV* nonidentity: in local coordinates u=1/t, X=u^4 x, Y=u^6 y,
# the section specializes to the singular point (0,0).
add("P1_IVstar_x",d["x1_4"])
add("P1_IVstar_y",d["y1_6"])

# P1 label (1,0,1,0,0): nonidentity only at I4@1 among finite fibers.
add("P1_I4_1_x",X1(1)-s1)
add("P1_I4_1_y",Y1(1))

if args.stage in ("p13","all"):
    X3=sum(d[f"x3_{i}"]*t^i for i in range(5))
    Y3=sum(d[f"y3_{i}"]*t^i for i in range(7))
    coeff_eqs(Y3^2-X3^3-A*X3-B,"P3")
    add("P3_IVstar_x",d["x3_4"])
    add("P3_IVstar_y",d["y3_6"])
    # P3=(2,3,0,0,0): nonidentity only at I4@0 among finite fibers.
    add("P3_I4_0_x",X3(0)-s0)
    add("P3_I4_0_y",Y3(0))

if args.stage=="all":
    r2=d["r2"]; c12=d["c12"]; z2=t-r2
    X2=sum(d[f"x2_{i}"]*t^i for i in range(7))
    Y2=sum(d[f"y2_{i}"]*t^i for i in range(10))
    coeff_eqs(Y2^2-X2^3-A*X2*z2^4-B*z2^6,"P2")
    # IV* nonidentity
    add("P2_IVstar_x",d["x2_6"])
    add("P2_IVstar_y",d["y2_9"])
    # P2=(1,1,2,1,1): nonidentity at every finite reducible fiber.
    for name,a,s in [
        ("I4_0",K(0),s0),
        ("I4_1",K(1),s1),
        ("I2_lam",lam,sl),
        ("I2_mu",mu,sm),
    ]:
        add("P2_"+name+"_x",X2(a)-s*(a-r2)^2)
        add("P2_"+name+"_y",Y2(a))
    # At t=c12, P2=-P1.  Because the sum has a nonidentity class at every
    # reducible fiber, saturating c12 away from those fibers, r2, and Delta=0
    # makes this the unique required good-fiber intersection with O.
    add("P12_intersection_x", X2(c12)-X1(c12)*(c12-r2)^2)
    add("P12_intersection_y", Y2(c12)+Y1(c12)*(c12-r2)^3)

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6SECT|stage={args.stage}|field={field}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)

groups={}
for tag in tags:
    g=tag.split("_")[0]
    groups[g]=groups.get(g,0)+1
print("E6SECT|groups="+",".join(f"{k}:{v}" for k,v in groups.items()),flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    if i<args.show_first:
        print(f"E6SECT_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}|eq={e}",flush=True)
    elif i<35 or i>=len(eqs)-10:
        print(f"E6SECT_META|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

print("E6SECT|component_target=P1:(1,0,1,0,0);P2:(1,1,2,1,1);P3:(2,3,0,0,0)",flush=True)
print("E6SECT|intersection_target=P12=2,P13=2,P23=2; P1O=0,P2O=1,P3O=0",flush=True)
if args.stage=="all":
    print("E6SECT|pair_gate=(P1+P2).O=1|witness=c12|encoded=1",flush=True)
    print("E6SECT|pair_gate_open=c12*(c12-1)*(c12-lam)*(c12-mu)*(c12-r2)*Delta(c12)!=0",flush=True)
print("E6SECT|note=IV* labels 1 vs 2 and I4 labels 1/2/3 are not yet separated beyond first singular hit.",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6SECT|export={out}",flush=True)
