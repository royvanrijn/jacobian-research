from sage.all import *
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]
J=BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"

def readmat(p):
    rows=[]
    for line in p.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        line=line.replace("["," ").replace("]"," ")
        try:
            r=[ZZ(x) for x in line.split()]
        except:
            continue
        if r:
            rows.append(r)
    return matrix(ZZ,rows)

H=readmat(J/"rank17_gram.txt")
B=readmat(J/"short_vector_basis_coords.txt")

print("det B =",B.det())
assert abs(B.det())==1

n=17
coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(H[i,i]//2 if i==j else H[i,j])

Q=QuadraticForm(ZZ,n,coeff)
V=[vector(ZZ,v)
   for v in Q.short_vector_list_up_to_length(3,True)[2]]

Binv=B.inverse()
assert Binv.base_ring()==QQ

W=[]

for v in V:
    # rows use old coordinates; solve c*B=v.
    c=v*Binv
    assert all(x.denominator()==1 for x in c)
    W.append(vector(ZZ,c))

out=J/"all_1311_in_short_basis.txt"
out.write_text(
    "\n".join(" ".join(map(str,v)) for v in W)+"\n"
)

print("vectors =",len(W))
print("max coefficient =",max(abs(x) for v in W for x in v))

from collections import Counter
cc=Counter(abs(x) for v in W for x in v)
print("coefficient histogram:")
for x in sorted(cc):
    print(x,cc[x])
