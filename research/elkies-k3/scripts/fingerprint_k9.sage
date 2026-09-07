from sage.all import *
from pathlib import Path
import json

BASE = Path(__file__).resolve().parents[1]
D = BASE / "data/lattice"
R = BASE / "results"

def readmat(path):
    rows=[]
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append([ZZ(x) for x in line.split()])
    return matrix(ZZ,rows)

G = readmat(D/"coxeter9_gram.txt")
n = G.nrows()

print("rank =", n)
print("det =", G.det())

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(G[i,i]//2 if i==j else G[i,j])

Q = QuadraticForm(ZZ,n,coeff)

S = Q.short_vector_list_up_to_length(11, True)

print("\nTHETA SHELLS")
theta={}
for q in range(1,len(S)):
    signed = 2*len(S[q])
    theta[2*q]=signed
    print("norm",2*q,"signed",signed)

# dual minimum.
#
# G^{-1} is rational. Clear denominators to get an integral
# Gram matrix; it need not be even, so use IntegralLattice
# directly rather than QuadraticForm's even-Hessian convention.
Gi = G.inverse()
den = lcm(x.denominator() for x in Gi.list())
Gdual = matrix(ZZ, den*Gi)

assert Gdual.is_positive_definite()

LD = IntegralLattice(Gdual)
dual_scaled_min = LD.minimum()
dual_min = QQ(dual_scaled_min) / den

print("\ndual clearing denominator =", den)
print("dual scaled minimum =", dual_scaled_min)
print("dual minimum =", dual_min)

minimum=4
gamma = QQ(minimum)/(QQ(G.det())**(QQ(1)/n))
print("Hermite invariant ~",N(gamma,20))

L=IntegralLattice(G)
A=L.orthogonal_group(is_finite=True)

print("\nAUT")
print("order =",A.order())
print("gens =",len(A.gens()))

try:
    print("structure =",A.structure_description())
except Exception as e:
    print("structure unavailable:",e)

# perfection:
mins=[vector(QQ,v) for v in S[2]]
mats=[]
for v in mins:
    M=v.column()*v.row()
    mats.append(vector(QQ,M.list()))

perfect_rank=matrix(QQ,mats).rank()

print("\nPERFECTION")
print("minimal pairs =",len(mins))
print("span rank of vv^T =",perfect_rank)
print("Sym^2 dimension =",n*(n+1)//2)
print("perfect =",perfect_rank == n*(n+1)//2)

(R/"k9-fingerprint.json").write_text(json.dumps({
    "rank":n,
    "det":int(G.det()),
    "minimum":int(4),
    "theta_signed":{int(k): int(v) for k,v in theta.items()},
    "aut_order":int(A.order()),
    "perfect_rank":int(perfect_rank),
}, indent=2)+"\n")
