from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Build/reduce semistable A10+A2+A1^2 MW3 P1 system.")
ap.add_argument("--p",type=int,default=31,help="finite field prime")
ap.add_argument("--out",default="/tmp/mw3-p1-reduced.ms")
ap.add_argument("--max-elim",type=int,default=40)
ap.add_argument("--show",action="store_true")
args=ap.parse_args()

if not is_prime(args.p) or args.p in (2,3,11,79):
    raise SystemExit("choose good prime != 2,3,11,79")

K=GF(args.p)

# Short Weierstrass K3:
# y^2 = x^3 + A(t)x + B(t), deg A<=8, deg B<=12.
#
# Fibers:
# I11 @ infinity
# I3  @ 0
# I2  @ 1
# I2  @ lambda
#
# Reduced P1 component profile:
# (I11,I3,I2@1,I2@lam) = (2,1,0,1), P1.O=0.
#
# P1.O=0 => x polynomial degree <=4, y degree <=6.
# Nonidentity at I3@0 forces P1 through its node.
# Nonidentity at I2@lambda forces P1 through its node.
# Identity at I2@1 is an OPEN condition, imposed only as verification.
#
# For I11 infinity, the discriminant has degree exactly 13.
# Thus coefficients t^24,...,t^14 vanish. D13 !=0 is open.
# P1 component class 2 at infinity implies in the minimal infinity chart
# X=u^4*x(1/u), Y=u^6*y(1/u), P1 hits the nodal point.
# We impose that first incidence exactly: x4=sinf, y6=0.
# Deeper class-2 separation is checked after solving, not guessed here.

names=[]
names += [f"a{i}" for i in range(9)]
names += [f"b{i}" for i in range(13)]
names += ["lam","sinf","s0","s1","sl"]
# P1: x0..x4, y0..y6; y0 is forced 0, x0=s0, x4=sinf, y6=0,
# but keep only genuinely free coefficients.
names += ["x1","x2","x3"]
names += ["y1","y2","y3","y4","y5"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t")
t=Rt.gen()

def V(n): return RF(d[n])

aa=[V(f"a{i}") for i in range(9)]
bb=[V(f"b{i}") for i in range(13)]
lam,sinf,s0,s1,sl=[V(n) for n in ["lam","sinf","s0","s1","sl"]]
x1,x2,x3=[V(n) for n in ["x1","x2","x3"]]
y1,y2,y3,y4,y5=[V(n) for n in ["y1","y2","y3","y4","y5"]]

A=sum(aa[i]*t**i for i in range(9))
B=sum(bb[i]*t**i for i in range(13))

# P1 finite t-chart. At t=0 it passes through node (s0,0).
# At infinity X=u^4*x(1/u) tends to x4=sinf and Y tends to y6=0.
X=s0 + x1*t + x2*t**2 + x3*t**3 + sinf*t**4
Y=y1*t + y2*t**2 + y3*t**3 + y4*t**4 + y5*t**5  # y0=y6=0

S=Y**2-X**3-A*X-B
Delta=-16*(4*A**3+27*B**2)

eqs=[]
tags=[]

def add(tag,e):
    e=RF(e)
    if e==0: return
    eqs.append(e)
    tags.append(tag)

# ---- I3 at t=0 ----
# Node (s0,0):
add("I3_0_A", A(0)+3*s0**2)
add("I3_0_B", B(0)-2*s0**3)
# Delta'(0)=0:
add("I3_0_d1", Delta[1])
# Delta'' condition is coefficient Delta_2=0.
add("I3_0_d2", Delta[2])

# ---- I2 at t=1 ----
add("I2_1_A", A(1)+3*s1**2)
add("I2_1_B", B(1)-2*s1**3)
add("I2_1_d1", Delta.derivative(t)(1))

# ---- I2 at t=lambda ----
add("I2_lam_A", A(lam)+3*sl**2)
add("I2_lam_B", B(lam)-2*sl**3)
add("I2_lam_d1", Delta.derivative(t)(lam))

# P1 is nonidentity at lambda: it passes through that node.
add("P1_lam_X", X(lam)-sl)
add("P1_lam_Y", Y(lam))

# ---- I11 at infinity ----
# deg Delta <= 13. Generic exact I11 additionally requires Delta_13 != 0.
for k in range(24,13,-1):
    add(f"I11_inf_D{k}", Delta[k])

# Infinity nodal point in u=1/t chart:
# Ainf(0)=a8=-3*sinf^2, Binf(0)=b12=2*sinf^3.
add("I11_inf_A", aa[8]+3*sinf**2)
add("I11_inf_B", bb[12]-2*sinf**3)
# P1 incidence x4=sinf,y6=0 is built directly into X,Y.

# ---- P1 section identity ----
for k in range(max(S.degree()+1,1)):
    if S[k] != 0:
        add(f"P1_{k}",S[k])

print(
    f"MW3P1|stage=raw|p={args.p}|vars={len(names)}|eqs={len(eqs)}"
    f"|section_degree={S.degree()}|delta_degree={Delta.degree()}",
    flush=True
)

# ---------- Greedy safe affine-linear elimination ----------
# Work in fraction field and repeatedly choose equations affine-linear in a
# variable. Prefer sparse equations, then variables with high occurrence count.
active_vars=[V(n) for n in names]
subs={}
history=[]

def settle(fr):
    fr=RF(fr)
    for _ in range(20):
        old=fr
        fr=RF(fr.subs(subs))
        if fr==old: break
    return fr

def live_equations():
    out=[]
    for tag,e in zip(tags,eqs):
        z=settle(e)
        if z!=0:
            out.append((tag,z))
    return out

for step in range(args.max_elim):
    live=live_equations()
    candidates=[]
    for ei,(tag,e) in enumerate(live):
        num=R(e.numerator())
        # denominator is an open condition after previous eliminations.
        for v in list(active_vars):
            vv=R(v)
            try:
                deg=num.degree(vv)
            except Exception:
                continue
            if deg==1:
                c=RF(num.derivative(vv))
                if c==0: continue
                # derivative must not still depend on v.
                if R(c.numerator()).degree(vv)!=0 or R(c.denominator()).degree(vv)!=0:
                    continue
                # Prefer sparse equations and later/high-index Weierstrass coeffs.
                score=(len(num.monomials()), num.total_degree(), names.index(str(v)))
                candidates.append((score,tag,e,v))
    if not candidates:
        break
    candidates.sort(key=lambda z:z[0])
    score,tag,e,v=candidates[0]
    ev=settle(e)
    vv=R(v)
    num=R(ev.numerator())
    c=RF(num.derivative(vv))
    rhs=RF(-num.subs({vv:0})/c)
    subs[v]=rhs
    active_vars.remove(v)
    history.append((str(v),tag,rhs))
    print(
        f"MW3P1|elim={len(history)}|var={v}|from={tag}"
        f"|terms={score[0]}|degree={score[1]}",
        flush=True
    )

live=live_equations()

print(
    f"MW3P1|stage=reduced|eliminated={len(history)}"
    f"|remaining_vars={len(active_vars)}|remaining_eqs={len(live)}",
    flush=True
)
print("MW3P1|remaining="+",".join(map(str,active_vars)),flush=True)

for i,(tag,e) in enumerate(live):
    num=R(e.numerator())
    print(
        f"MW3P1_EQ|i={i}|tag={tag}|degree={num.total_degree()}"
        f"|terms={len(num.monomials())}",
        flush=True
    )
    if args.show:
        print(f"MW3P1_EXPR|i={i}|{num}",flush=True)

# ---------- Export reduced numerator system to msolve ----------
rem_names=[str(v) for v in active_vars]
RR=PolynomialRing(K,rem_names,order="degrevlex")
rd=RR.gens_dict()

def convert_num(fr):
    fr=settle(fr)
    num=R(fr.numerator())
    out=RR(0)
    for mon,coef in num.dict().items():
        term=RR(coef)
        for i,exp in enumerate(mon):
            if not exp: continue
            nm=names[i]
            if nm in rd:
                term*=rd[nm]**exp
            else:
                # eliminated variables should have disappeared after settle()
                raise RuntimeError(f"eliminated variable survived conversion: {nm}")
        out+=term
    return out

red=[]
redtags=[]
for tag,e in live:
    try:
        q=convert_num(e)
    except RuntimeError:
        # Some substitutions depend on other eliminated vars and need more
        # compositional settling. Report clearly instead of silently corrupting.
        print(f"MW3P1|export_skip|tag={tag}|reason=unsettled_eliminated_variable",flush=True)
        continue
    if q!=0:
        red.append(q); redtags.append(tag)

out=Path(args.out)
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(rem_names)+"\n")
    h.write(str(args.p)+"\n")
    for i,q in enumerate(red):
        h.write(str(q).replace("**","^"))
        h.write(",\n" if i+1<len(red) else "\n")

meta=out.with_suffix(".meta.txt")
with meta.open("w") as h:
    h.write("prime="+str(args.p)+"\n")
    h.write("remaining="+repr(rem_names)+"\n")
    h.write("open_conditions=Delta13!=0, lambda!=0,1, P1(1) != singular node\n")
    h.write("profile=(I11,I3,I2@1,I2@lambda)=(2,1,0,1), deeper I11 class2 not yet imposed\n")
    h.write("\nELIMINATIONS\n")
    for v,tag,rhs in history:
        h.write(f"{v} <- {rhs}    # {tag}\n")
    h.write("\nEQUATIONS\n")
    for tag,q in zip(redtags,red):
        h.write(f"{tag}: {q}\n")

print(f"MW3P1|stage=export|out={out}|meta={meta}|export_eqs={len(red)}",flush=True)
print("MW3P1|NOTE=I11 component-2 deeper blowup condition intentionally deferred; first solve/reduce this exact incidence system.",flush=True)
