from sage.all import *
from pathlib import Path
import itertools

# Reduced MW height Gram from exact glue computation.
H = matrix(QQ,[
    [QQ(79)/66, QQ(17)/66, QQ(-1)/66],
    [QQ(17)/66, QQ(106)/66, QQ(19)/66],
    [QQ(-1)/66, QQ(19)/66, QQ(259)/66],
])
assert H.det()==QQ(79)/11

# Semistable fibers: I11, I3, I2, I2.
ns=[11,3,2,2]

def contr_diag(n,k):
    if k==0: return QQ(0)
    return QQ(k*(n-k)) / QQ(n)

def contr_pair(n,i,j):
    # Inverse Cartan pairing for A_{n-1}, components numbered 0..n-1.
    if i==0 or j==0:
        return QQ(0)
    a=min(i,j); b=max(i,j)
    return QQ(a*(n-b)) / QQ(n)

labels=list(itertools.product(range(11),range(3),range(2),range(2)))
print(f"MW3COMP|stage=start|labels={len(labels)}|height_det={H.det()}",flush=True)

# For a K3, chi(O_X)=2. Shioda:
# <P,P> = 4 + 2(P.O) - sum_v contr_v(P)
# hence P.O = (<P,P>-4+sum contr)/2 must be a nonnegative integer.
diag_candidates=[]
for g in range(3):
    good=[]
    for lab in labels:
        c=sum(contr_diag(n,k) for n,k in zip(ns,lab))
        PO=(H[g,g]-4+c)/2
        if PO.denominator()==1 and PO>=0:
            good.append((lab,ZZ(PO),c))
    diag_candidates.append(good)
    print(f"MW3COMP|generator={g+1}|diag_candidates={len(good)}",flush=True)
    for row in good[:20]:
        print(f"MW3COMP_DIAG|g={g+1}|label={row[0]}|PdotO={row[1]}|local={row[2]}",flush=True)

def pair_intersection(g,h,A,B):
    labA,POA,_=A
    labB,POB,_=B
    c=sum(contr_pair(n,i,j) for n,i,j in zip(ns,labA,labB))
    # <P,Q> = 2 + P.O + Q.O - P.Q - sum contr(P,Q)
    PQ=QQ(2+POA+POB)-c-H[g,h]
    return PQ,c

# Pairwise-compatible candidate pairs.
pairs={}
for g,h in [(0,1),(0,2),(1,2)]:
    arr=[]
    for A in diag_candidates[g]:
        for B in diag_candidates[h]:
            PQ,c=pair_intersection(g,h,A,B)
            if PQ.denominator()==1 and PQ>=0:
                arr.append((A,B,ZZ(PQ),c))
    pairs[(g,h)]=arr
    print(f"MW3COMP|pair={g+1},{h+1}|compatible={len(arr)}",flush=True)

# Index pair compatibility for fast triple join.
ok01={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(0,1)]}
ok02={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(0,2)]}
ok12={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(1,2)]}

triples=[]
for A in diag_candidates[0]:
    ka=(A[0],A[1])
    for B in diag_candidates[1]:
        kb=(B[0],B[1])
        if (A[0],A[1],B[0],B[1]) not in ok01: continue
        for C in diag_candidates[2]:
            if (A[0],A[1],C[0],C[1]) not in ok02: continue
            if (B[0],B[1],C[0],C[1]) not in ok12: continue
            pq01,_=pair_intersection(0,1,A,B)
            pq02,_=pair_intersection(0,2,A,C)
            pq12,_=pair_intersection(1,2,B,C)
            triples.append((A,B,C,ZZ(pq01),ZZ(pq02),ZZ(pq12)))

print(f"MW3COMP|stage=triples|count={len(triples)}",flush=True)

# Score simpler geometric configurations first: low O intersections and low mutual intersections.
triples.sort(key=lambda z:(
    z[0][1]+z[1][1]+z[2][1]+z[3]+z[4]+z[5],
    z[0][1]+z[1][1]+z[2][1],
    z[3]+z[4]+z[5],
    z[0][0],z[1][0],z[2][0]
))

for i,T in enumerate(triples[:100],1):
    A,B,C,p12,p13,p23=T
    print(
        f"MW3COMP_BEST|rank={i}"
        f"|P1={A[0]}|P1O={A[1]}"
        f"|P2={B[0]}|P2O={B[1]}"
        f"|P3={C[0]}|P3O={C[1]}"
        f"|P1P2={p12}|P1P3={p13}|P2P3={p23}",
        flush=True
    )

out=Path("artifacts/local/elkies-k3/mw3-component-labels.txt")
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as f:
    f.write(f"H={H}\n")
    f.write(f"triples={len(triples)}\n")
    for i,T in enumerate(triples,1):
        A,B,C,p12,p13,p23=T
        f.write(
            f"{i}: P1={A[0]} P1O={A[1]} "
            f"P2={B[0]} P2O={B[1]} "
            f"P3={C[0]} P3O={C[1]} "
            f"P12={p12} P13={p13} P23={p23}\n"
        )

print(f"MW3COMP|stage=done|out={out}",flush=True)
