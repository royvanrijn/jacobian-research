from sage.all import *
import itertools
from pathlib import Path

# Exact reduced MW Gram for E6+A3^2+A1^2 node.
H=matrix(QQ,[
    [QQ(23)/12,QQ(-10)/12,QQ(-8)/12],
    [QQ(-10)/12,QQ(23)/12,QQ(1)/12],
    [QQ(-8)/12,QQ(1)/12,QQ(23)/12],
])
assert H.det()==QQ(79)/16

# Labels:
# IV*: class 0,1,2 in component group Z/3.
# I4a,I4b: classes 0..3.
# I2a,I2b: classes 0..1.
labels=list(itertools.product(range(3),range(4),range(4),range(2),range(2)))

def ivstar_diag(k):
    return QQ(0) if k==0 else QQ(4)/3

def ivstar_pair(i,j):
    if i==0 or j==0:
        return QQ(0)
    return QQ(4)/3 if i==j else QQ(2)/3

def In_diag(n,k):
    return QQ(0) if k==0 else QQ(k*(n-k))/QQ(n)

def In_pair(n,i,j):
    if i==0 or j==0:
        return QQ(0)
    a=min(i,j); b=max(i,j)
    return QQ(a*(n-b))/QQ(n)

def diag_contr(lab):
    e,a,b,c,d=lab
    return ivstar_diag(e)+In_diag(4,a)+In_diag(4,b)+In_diag(2,c)+In_diag(2,d)

def pair_contr(A,B):
    e,a,b,c,d=A
    f,g,h,i,j=B
    return ivstar_pair(e,f)+In_pair(4,a,g)+In_pair(4,b,h)+In_pair(2,c,i)+In_pair(2,d,j)

print(f"E6COMP|stage=start|labels={len(labels)}|det={H.det()}",flush=True)

diag=[]
for g in range(3):
    good=[]
    for lab in labels:
        c=diag_contr(lab)
        PO=(H[g,g]-4+c)/2
        if PO.denominator()==1 and PO>=0:
            good.append((lab,ZZ(PO),c))
    diag.append(good)
    print(f"E6COMP|generator={g+1}|diag_candidates={len(good)}",flush=True)
    for x in good[:30]:
        print(f"E6COMP_DIAG|g={g+1}|label={x[0]}|PdotO={x[1]}|local={x[2]}",flush=True)

def PQ(g,h,A,B):
    la,oa,_=A; lb,ob,_=B
    c=pair_contr(la,lb)
    pq=QQ(2+oa+ob)-c-H[g,h]
    return pq,c

pairs={}
for g,h in [(0,1),(0,2),(1,2)]:
    arr=[]
    for A in diag[g]:
        for B in diag[h]:
            pq,c=PQ(g,h,A,B)
            if pq.denominator()==1 and pq>=0:
                arr.append((A,B,ZZ(pq),c))
    pairs[(g,h)]=arr
    print(f"E6COMP|pair={g+1},{h+1}|compatible={len(arr)}",flush=True)

ok01={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(0,1)]}
ok02={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(0,2)]}
ok12={(A[0],A[1],B[0],B[1]) for A,B,_,_ in pairs[(1,2)]}

triples=[]
for A in diag[0]:
    for B in diag[1]:
        if (A[0],A[1],B[0],B[1]) not in ok01: continue
        for C in diag[2]:
            if (A[0],A[1],C[0],C[1]) not in ok02: continue
            if (B[0],B[1],C[0],C[1]) not in ok12: continue
            p12,_=PQ(0,1,A,B); p13,_=PQ(0,2,A,C); p23,_=PQ(1,2,B,C)
            triples.append((A,B,C,ZZ(p12),ZZ(p13),ZZ(p23)))

triples.sort(key=lambda z:(
    z[0][1]+z[1][1]+z[2][1]+z[3]+z[4]+z[5],
    z[0][1]+z[1][1]+z[2][1],
    z[3]+z[4]+z[5],
    z[0][0],z[1][0],z[2][0]
))
print(f"E6COMP|stage=triples|count={len(triples)}",flush=True)
for i,T in enumerate(triples[:100],1):
    A,B,C,p12,p13,p23=T
    print(f"E6COMP_BEST|rank={i}|P1={A[0]}|P1O={A[1]}|P2={B[0]}|P2O={B[1]}|P3={C[0]}|P3O={C[1]}|P12={p12}|P13={p13}|P23={p23}",flush=True)

out=Path("artifacts/local/elkies-k3/e6-component-labels.txt")
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as f:
    f.write(f"H={H}\ntriples={len(triples)}\n")
    for i,T in enumerate(triples,1):
        A,B,C,p12,p13,p23=T
        f.write(f"{i}: P1={A[0]} O={A[1]} P2={B[0]} O={B[1]} P3={C[0]} O={C[1]} P12={p12} P13={p13} P23={p23}\n")
print(f"E6COMP|stage=done|out={out}",flush=True)
