from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Reconstruct and verify the first explicit GF(31) point on the E6/P1 locus.")
ap.add_argument("--meta",default="artifacts/local/elkies-k3/e6-base.meta.txt")
args=ap.parse_args()

p=31
K=GF(p)

# Successful coordinate slice and reduced solution.
vals0={
    "r0":K(4), "s0":K(18), "x1":K(27),
    "a1":K(4), "a2":K(16), "a4":K(6), "s1":K(23),
    "lam":K(24), "mu":K(18), "sl":K(23), "sm":K(4),
}

r0=vals0["r0"]; s0=vals0["s0"]; x1=vals0["x1"]
x0=r0**2-2*s0
y0=r0*(r0**2-3*s0)
y1=(vals0["a1"]+3*(r0**2-s0)*x1)/(2*r0)

vals=dict(vals0)
vals.update({"x0":x0,"y0":y0,"y1":y1})

print("E6RECON|base|"+"|".join(f"{k}={int(v)}" for k,v in vals.items()),flush=True)

# Parse the four rational elimination expressions saved by export_e6_p1_sliced.sage.
meta=Path(args.meta)
if not meta.exists():
    raise SystemExit(f"missing meta file: {meta}")
md={}
for line in meta.read_text().splitlines():
    if "=" in line:
        k,v=line.split("=",1)
        md[k.strip()]=v.strip()

needed=["a3","y2","b5","b4"]
if any(k not in md for k in needed):
    raise RuntimeError(f"meta missing one of {needed}")

names=["a1","a2","a3","a4","b4","b5","lam","mu","s0","s1","sl","sm","x0","x1","y0","y1","y2"]
R=PolynomialRing(K,names)
d=R.gens_dict(); F=FractionField(R)

subs={d[k]:F(v) for k,v in vals.items() if k in d}

def eval_expr(text):
    e=F(text.replace("^","**"))
    for _ in range(8):
        old=e
        e=F(e.subs(subs))
        if e==old:
            break
    if e.denominator()==0:
        raise ZeroDivisionError(text)
    if any(str(g) in str(e) for g in R.gens()):
        # This is only a diagnostic; normal reconstruction should fully specialize.
        pass
    return K(e)

# Reverse dependency order: a3, y2, b5, b4.
for k in needed:
    v=eval_expr(md[k])
    vals[k]=v
    subs[d[k]]=F(v)
    print(f"E6RECON|elim|{k}={int(v)}",flush=True)

# Rebuild full Weierstrass coefficients.
a1=vals["a1"]; a2=vals["a2"]; a3=vals["a3"]; a4=vals["a4"]
s0=vals["s0"]; s1=vals["s1"]
b4=vals["b4"]; b5=vals["b5"]

a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*a3*s0**4)/(216*s0**3)
a5=-3*s1**2-(a0+a1+a2+a3+a4)

# Solve b6,b7 from the t=1 singular-root conditions, exactly as in the exporter.
U=PolynomialRing(K,["B6"]); B6=U.gen()
Ut=PolynomialRing(U,"t"); tt=Ut.gen()
AA=sum(U(c)*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
B7=2*U(s1)**3-(U(b0)+U(b1)+U(b2)+U(b3)+U(b4)+U(b5)+B6)-1
BB=U(b0)+U(b1)*tt+U(b2)*tt**2+U(b3)*tt**3+U(b4)*tt**4+U(b5)*tt**5+B6*tt**6+B7*tt**7+tt**8
EE=U(BB.derivative(tt)(1)+U(s1)*AA.derivative(tt)(1))
c6=EE.derivative(B6)
b6=K(-EE.subs({B6:0})/c6)
b7=K(B7.subs({B6:U(b6)}))

T=PolynomialRing(K,"t"); t=T.gen()
A=T(a0)+T(a1)*t+T(a2)*t**2+T(a3)*t**3+T(a4)*t**4+T(a5)*t**5
B=T(b0)+T(b1)*t+T(b2)*t**2+T(b3)*t**3+T(b4)*t**4+T(b5)*t**5+T(b6)*t**6+T(b7)*t**7+t**8

x0=vals["x0"]; x1=vals["x1"]; y0=vals["y0"]; y1=vals["y1"]; y2=vals["y2"]
x2=s1-x0-x1
y3=-1-y0-y1-y2
X=T(x0)+T(x1)*t+T(x2)*t**2
Y=T(y0)+T(y1)*t+T(y2)*t**2+T(y3)*t**3+t**4

sec=Y**2-X**3-A*X-B
Delta=-16*(4*A**3+27*B**2)
Dp=Delta.derivative(t)

print(f"E6RECON|A={A}",flush=True)
print(f"E6RECON|B={B}",flush=True)
print(f"E6RECON|X={X}",flush=True)
print(f"E6RECON|Y={Y}",flush=True)
print(f"E6RECON|section_zero={sec==0}",flush=True)

# Fiber checks: I4 at t=0,1 via Delta multiplicities; I2 at lambda,mu.
for label,z in [("0",K(0)),("1",K(1)),("lam",vals["lam"]),("mu",vals["mu"])]:
    mult=0
    D=Delta
    while D(z)==0 and D!=0:
        mult+=1
        D=D.derivative(t)
    print(f"E6RECON|fiber={label}|z={int(z)}|delta_mult={mult}",flush=True)

for label,z,s in [("lam",vals["lam"],vals["sl"]),("mu",vals["mu"],vals["sm"])]:
    okA=(A(z)==-3*s**2)
    okB=(B(z)==2*s**3)
    okD=(B.derivative(t)(z)+s*A.derivative(t)(z)==0)
    print(f"E6RECON|I2={label}|A={okA}|B={okB}|d1={okD}",flush=True)

print("E6RECON|full_zero="+str(sec==0),flush=True)
