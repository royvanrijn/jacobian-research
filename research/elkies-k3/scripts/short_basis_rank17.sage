from sage.all import *
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]
J=BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"
p=J/"rank17_gram.txt"

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

H=matrix(ZZ,rows)
n=17

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(H[i,i]//2 if i==j else H[i,j])

Q=QuadraticForm(ZZ,n,coeff)

S=Q.short_vector_list_up_to_length(3,True)
V=[vector(ZZ,v) for v in S[2]]

# Sage gives one per +/- pair with up_to_sign=True.
print("pairs norm4 =",len(V))

chosen=[]
M=matrix(ZZ,0,n)

for v in V:
    N=matrix(ZZ,chosen+[v])
    if N.rank()>len(chosen):
        chosen.append(v)
        M=N
        print("rank",len(chosen),"vector",v)
        if len(chosen)==n:
            break

assert len(chosen)==17

B=matrix(ZZ,chosen)
HB=B*H*B.T

print()
print("det original =",H.det())
print("det short-span Gram =",HB.det())

ratio=ZZ(HB.det()//H.det())
print("index^2 =",ratio)
print("index =",sqrt(ratio))

(J/"short_vector_basis_coords.txt").write_text(
    "\n".join(" ".join(map(str,v)) for v in B.rows())+"\n"
)

(J/"short_vector_basis_gram.txt").write_text(
    "\n".join(" ".join(map(str,v)) for v in HB.rows())+"\n"
)

print()
print("short basis Gram:")
print(HB)
