from sage.all import *
from itertools import product

# Reduced exact MW height Gram for A10+A2+A1+A1 branch.
G=matrix(QQ,[
    [QQ(79)/66, QQ(17)/66, -QQ(1)/66],
    [QQ(17)/66, QQ(106)/66, QQ(19)/66],
    [-QQ(1)/66, QQ(19)/66, QQ(259)/66],
])

# Fibers: I11, I3, I2, I2
ns=[11,3,2,2]
names=["I11","I3","I2a","I2b"]

def self_contr(n,c):
    return QQ(c*(n-c))/n

def pair_contr(n,a,b):
    a,b=min(a,b),max(a,b)
    return QQ(a*(n-b))/n

# For a K3, chi(O_X)=2:
# <P,P> = 4 + 2(P.O) - sum_v contr_v(P)
# <P,Q> = 2 + P.O + Q.O - P.Q - sum_v contr_v(P,Q)

opts=[]
for i in range(3):
    oi=[]
    for cs in product(*[range(n) for n in ns]):
        corr=sum(self_contr(n,c) for n,c in zip(ns,cs))
        for po in range(0,4):
            if QQ(4+2*po)-corr == G[i,i]:
                oi.append((tuple(cs),po,corr))
    opts.append(oi)
    print(f"MW3PROF|generator={i+1}|diagonal_options={len(oi)}")
    for cs,po,corr in oi:
        print(f"MW3PROF|diag|P={i+1}|components={cs}|P.O={po}|correction={corr}")

sol=[]
for a in opts[0]:
    for b in opts[1]:
        for c in opts[2]:
            arr=[a,b,c]
            ints={}
            good=True
            for i,j in [(0,1),(0,2),(1,2)]:
                ci,Oi,_=arr[i]
                cj,Oj,_=arr[j]
                corr=sum(pair_contr(n,x,y) for n,x,y in zip(ns,ci,cj))
                pij=QQ(2+Oi+Oj)-corr-G[i,j]
                if pij.denominator()!=1 or pij<0:
                    good=False
                    break
                ints[(i,j)]=ZZ(pij)
            if good:
                sol.append((arr,ints))

print(f"MW3PROF|full_solutions={len(sol)}")
for k,(arr,ints) in enumerate(sol,1):
    print(f"MW3PROF|solution={k}")
    for i,(cs,po,corr) in enumerate(arr,1):
        print(f"MW3PROF|P{i}|components={cs}|P.O={po}|selfcorr={corr}")
    print(
        f"MW3PROF|intersections|P1.P2={ints[(0,1)]}"
        f"|P1.P3={ints[(0,2)]}|P2.P3={ints[(1,2)]}"
    )

# Canonicalize by obvious symmetries:
# - negate all component labels of any generator
# - swap the two I2 fibers globally.
# The first printed solution is a convenient representative.
if sol:
    arr,ints=sol[0]
    print("MW3PROF|CANONICAL")
    for i,(cs,po,corr) in enumerate(arr,1):
        print(f"MW3PROF|CANONICAL|P{i}|components={cs}|P.O={po}")
    print(
        f"MW3PROF|CANONICAL|intersections="
        f"{ints[(0,1)]},{ints[(0,2)]},{ints[(1,2)]}"
    )
