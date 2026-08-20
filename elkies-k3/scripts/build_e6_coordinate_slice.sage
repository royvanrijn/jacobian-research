from sage.all import *
from pathlib import Path
import argparse, random, time

ap=argparse.ArgumentParser(description="Build a correct 8x8 coordinate slice from the E6 P1 param0 system.")
ap.add_argument("--input",required=True)
ap.add_argument("--out",required=True)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--kill",default="r0,x1,a4")
args=ap.parse_args()

lines=[x.strip() for x in Path(args.input).read_text().splitlines() if x.strip()]
names=[x.strip() for x in lines[0].split(",") if x.strip()]
p=int(lines[1]); K=GF(p)
R=PolynomialRing(K,names,order="degrevlex"); d=R.gens_dict()
eqs=[R((x[:-1] if x.endswith(",") else x).replace("^","**")) for x in lines[2:]]
if len(eqs)!=11: raise RuntimeError(f"expected 11 param0 equations, got {len(eqs)}")

kill=[x.strip() for x in args.kill.split(",") if x.strip()]
if len(kill)!=3 or "r0" not in kill: raise RuntimeError("kill must contain exactly 3 variables including r0")
if "y1" in kill: raise RuntimeError("do not kill y1; it is eliminated from P1_1")
if any(n not in d for n in kill): raise RuntimeError(f"unknown kill variable in {kill}")

random.seed(int(args.seed))
vals={n:K(random.randrange(int(p))) for n in kill}
if vals["r0"]==0: vals["r0"]=K(1+random.randrange(int(p)-1))

keep=[n for n in names if n not in kill]
S=PolynomialRing(K,keep,order="degrevlex"); sd=S.gens_dict()
phi=R.hom([S(vals[n]) if n in kill else sd[n] for n in names],S)

print(f"E6COORD|stage=slice|p={p}|seed={args.seed}|kill={','.join(kill)}|vals="+",".join(f"{n}={int(vals[n])}" for n in kill),flush=True)
Es=[]
for i,e in enumerate(eqs):
    t0=time.time(); f=phi(e); Es.append(f)
    print(f"E6COORD|slice_eq={i+1}/11|seconds={time.time()-t0:.3f}|degree={f.degree()}|terms={len(f.monomials())}",flush=True)

# Param0 order: 8 fiber equations, then P1_1,P1_2,P1_3.
y=sd["y1"]; e1=Es[8]
if e1.degree(y)!=1: raise RuntimeError(f"P1_1 is not affine-linear in y1 (degree={e1.degree(y)})")
c=e1.derivative(y)
F=FractionField(S)
yexpr=F(-e1.subs({y:0}))/F(c)
if yexpr.denominator().degree()!=0:
    raise RuntimeError(f"y1 denominator is not constant after coordinate slice: {yexpr.denominator()}")
print(f"E6COORD|stage=y1|num_terms={len(S(yexpr.numerator()).monomials())}|den={yexpr.denominator()}",flush=True)

final_names=[n for n in keep if n!="y1"]
T=PolynomialRing(K,final_names,order="degrevlex"); td=T.gens_dict()
map0=S.hom([td[n] if n!="y1" else T(0) for n in keep],T)
y_num=S(yexpr.numerator()); y_den=K(yexpr.denominator())
if y_den==0: raise RuntimeError("zero y1 denominator")
yT=map0(y_num)/T(y_den)
psi=S.hom([yT if n=="y1" else td[n] for n in keep],T)

p12=psi(Es[9]); p13=psi(Es[10])
print(f"E6COORD|check=P1_2_zero|value={p12==0}",flush=True)
print(f"E6COORD|check=P1_3_zero|value={p13==0}",flush=True)
if p12!=0 or p13!=0:
    raise RuntimeError("P1_2/P1_3 did not vanish; refusing inconsistent reduction")

fiber=[]
for i,e in enumerate(Es[:8]):
    t0=time.time(); f=psi(e); fiber.append(f)
    print(f"E6COORD|final_eq={i+1}/8|seconds={time.time()-t0:.3f}|degree={f.degree()}|terms={len(f.monomials())}",flush=True)

out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(final_names)+"\n"+str(p)+"\n")
    for i,e in enumerate(fiber):
        h.write(str(e).replace("**","^")); h.write(",\n" if i+1<len(fiber) else "\n")
meta=out.with_suffix(".meta.txt")
meta.write_text("seed="+str(args.seed)+"\nkill="+repr(kill)+"\nvals="+repr({n:int(vals[n]) for n in kill})+"\ny1="+str(yexpr)+"\n")
print(f"E6COORD|stage=done|vars={T.ngens()}|eqs={len(fiber)}|out={out}|meta={meta}",flush=True)
