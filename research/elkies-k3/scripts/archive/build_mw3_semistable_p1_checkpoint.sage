from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--p",type=int,default=31)
ap.add_argument("--steps",type=int,default=11)
ap.add_argument("--checkpoint",default="/tmp/mw3-p1-checkpoint.txt")
args=ap.parse_args()

K=GF(args.p)
names=[]
names += [f"a{i}" for i in range(9)]
names += [f"b{i}" for i in range(13)]
names += ["lam","sinf","s0","s1","sl"]
names += ["x1","x2","x3"]
names += ["y1","y2","y3","y4","y5"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t"); t=Rt.gen()
V=lambda n: RF(d[n])

aa=[V(f"a{i}") for i in range(9)]
bb=[V(f"b{i}") for i in range(13)]
lam,sinf,s0,s1,sl=[V(n) for n in ["lam","sinf","s0","s1","sl"]]
x1,x2,x3=[V(n) for n in ["x1","x2","x3"]]
y1,y2,y3,y4,y5=[V(n) for n in ["y1","y2","y3","y4","y5"]]

A=sum(aa[i]*t**i for i in range(9))
B=sum(bb[i]*t**i for i in range(13))
X=s0+x1*t+x2*t**2+x3*t**3+sinf*t**4
Y=y1*t+y2*t**2+y3*t**3+y4*t**4+y5*t**5
S=Y**2-X**3-A*X-B
Delta=-16*(4*A**3+27*B**2)

eqs=[]; tags=[]
def add(tag,e):
    e=RF(e)
    if e!=0:
        tags.append(tag); eqs.append(e)

add("I3_0_A",A(0)+3*s0**2)
add("I3_0_B",B(0)-2*s0**3)
add("I3_0_d1",Delta[1])
add("I3_0_d2",Delta[2])

add("I2_1_A",A(1)+3*s1**2)
add("I2_1_B",B(1)-2*s1**3)
add("I2_1_d1",Delta.derivative(t)(1))

add("I2_lam_A",A(lam)+3*sl**2)
add("I2_lam_B",B(lam)-2*sl**3)
add("I2_lam_d1",Delta.derivative(t)(lam))
add("P1_lam_X",X(lam)-sl)
add("P1_lam_Y",Y(lam))

for k in range(24,13,-1):
    add(f"I11_inf_D{k}",Delta[k])
add("I11_inf_A",aa[8]+3*sinf**2)
add("I11_inf_B",bb[12]-2*sinf**3)

for k in range(S.degree()+1):
    if S[k]!=0:
        add(f"P1_{k}",S[k])

active=[V(n) for n in names]
subs={}
history=[]

def settle_limited(fr, passes=4):
    fr=RF(fr)
    for _ in range(passes):
        old=fr
        fr=RF(fr.subs(subs))
        if fr==old: break
    return fr

for step in range(args.steps):
    best=None
    for tag,e0 in zip(tags,eqs):
        e=settle_limited(e0)
        if e==0: continue
        num=R(e.numerator())
        for v in active:
            vv=R(v)
            try: deg=num.degree(vv)
            except: continue
            if deg!=1: continue
            c=RF(num.derivative(vv))
            if c==0: continue
            try:
                if R(c.numerator()).degree(vv)!=0 or R(c.denominator()).degree(vv)!=0:
                    continue
            except:
                continue
            score=(len(num.monomials()),num.total_degree(),names.index(str(v)))
            cand=(score,tag,e,v)
            if best is None or score<best[0]:
                best=cand
    if best is None:
        print(f"MW3CHK|stalled|step={step}",flush=True)
        break

    score,tag,e,v=best
    vv=R(v); num=R(e.numerator())
    c=RF(num.derivative(vv))
    rhs=RF(-num.subs({vv:0})/c)
    subs[v]=rhs
    active.remove(v)
    history.append((str(v),tag,str(rhs)))
    print(f"MW3CHK|elim={len(history)}|var={v}|from={tag}|terms={score[0]}|degree={score[1]}",flush=True)

print(f"MW3CHK|CHECKPOINT|eliminated={len(history)}|remaining={len(active)}",flush=True)
print("MW3CHK|remaining="+",".join(map(str,active)),flush=True)

cp=Path(args.checkpoint)
with cp.open("w") as h:
    h.write(f"p={args.p}\n")
    h.write("remaining="+",".join(map(str,active))+"\n")
    h.write("\nELIMINATIONS\n")
    for i,(v,tag,rhs) in enumerate(history,1):
        h.write(f"{i}|{v}|{tag}|{rhs}\n")
    h.write("\nALL_TAGS\n")
    for tag in tags:
        h.write(tag+"\n")

print(f"MW3CHK|saved={cp}",flush=True)
print("MW3CHK|done_without_global_expansion",flush=True)
