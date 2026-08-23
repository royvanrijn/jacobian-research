from sage.all import *
from pathlib import Path
import argparse,time

ap=argparse.ArgumentParser()
ap.add_argument("--p",type=int,default=31)
ap.add_argument("--checkpoint",default="/tmp/mw3-p1-checkpoint.txt")
ap.add_argument("--timeout",type=float,default=8.0,
                help="soft reporting threshold per equation; cannot interrupt Sage algebra")
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

eqmap={}
def add(tag,e):
    e=RF(e)
    if e!=0: eqmap[tag]=e

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
for k in range(24,13,-1): add(f"I11_inf_D{k}",Delta[k])
add("I11_inf_A",aa[8]+3*sinf**2)
add("I11_inf_B",bb[12]-2*sinf**3)
for k in range(S.degree()+1):
    if S[k]!=0: add(f"P1_{k}",S[k])

# Reproduce exactly the known 11 safe eliminations, but only compute the
# equation needed at each step.
plan=[
 ("a0","I3_0_A"),
 ("a1","P1_1"),
 ("a8","I11_inf_A"),
 ("a7","P1_11"),
 ("b0","I3_0_B"),
 ("b12","I11_inf_B"),
 ("a2","I3_0_d2"),
 ("a6","I11_inf_D22"),
 ("a5","I11_inf_D21"),
 ("y1","P1_lam_Y"),
 ("sinf","P1_lam_X"),
]
subs={}

def settle(fr,maxpasses=20):
    fr=RF(fr)
    for _ in range(maxpasses):
        old=fr
        fr=RF(fr.subs(subs))
        if fr==old: break
    return fr

for vname,tag in plan:
    v=V(vname); vv=R(v)
    e=settle(eqmap[tag])
    num=R(e.numerator())
    c=RF(num.derivative(vv))
    rhs=RF(-num.subs({vv:0})/c)
    subs[v]=rhs
    print(f"MW3PROBE|replay|var={vname}|from={tag}",flush=True)

remaining=[n for n in names if n not in [x[0] for x in plan]]
used=set(tag for _,tag in plan)

# Probe section equations first, then local finite fibers, then I11 discriminant.
order=[]
order += [f"P1_{k}" for k in range(13)]
order += ["I3_0_d1","I2_1_A","I2_1_B","I2_1_d1",
          "I2_lam_A","I2_lam_B","I2_lam_d1"]
order += [f"I11_inf_D{k}" for k in range(20,13,-1)]

seen=set()
order=[x for x in order if x in eqmap and x not in used and not (x in seen or seen.add(x))]

print("MW3PROBE|remaining="+",".join(remaining),flush=True)

for tag in order:
    t0=time.time()
    try:
        e=settle(eqmap[tag])
        dt=time.time()-t0
        if e==0:
            print(f"MW3PROBE|tag={tag}|zero=1|seconds={dt:.3f}",flush=True)
            continue
        num=R(e.numerator())
        vars_here=[]
        linear=[]
        for n in remaining:
            v=R(d[n])
            deg=num.degree(v)
            if deg>0:
                vars_here.append(f"{n}:{deg}")
                if deg==1:
                    c=num.derivative(v)
                    if c.degree(v)==0:
                        linear.append(n)
        print(
            f"MW3PROBE|tag={tag}|seconds={dt:.3f}|degree={num.total_degree()}"
            f"|terms={len(num.monomials())}|vars={','.join(vars_here)}"
            f"|linear={','.join(linear)}",
            flush=True
        )
    except Exception as ex:
        print(f"MW3PROBE|tag={tag}|ERROR={type(ex).__name__}:{ex}",flush=True)
