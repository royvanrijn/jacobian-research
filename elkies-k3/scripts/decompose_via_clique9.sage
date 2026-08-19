from sage.all import *
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/lattice"

D.mkdir(parents=True,exist_ok=True)

def readmat(path):
    rows=[]
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        line=line.replace("["," ").replace("]"," ")
        try:
            row=[ZZ(x) for x in line.split()]
        except:
            continue
        if row:
            rows.append(row)
    return matrix(ZZ,rows)

H = readmat(J/"short_vector_basis_gram.txt")

signed = np.load(
    J/"all_2622_signed_short_basis.npy"
)

clique = [
    2127,2585,2431,2375,2051,
    1987,2115,1705,1483
]

Braw = matrix(ZZ,signed[clique].tolist())

# Saturated K9 inside Z^17.
Kmod = Braw.row_module(ZZ).saturation()
K = matrix(ZZ,Kmod.basis_matrix())

GK = K*H*K.T

print("========== K9 ==========")
print("rank =",K.rank())
print("det =",GK.det())
print("Gram:")
print(GK)

# ------------------------------------------------------------
# Orthogonal complement:
#
# x H K^T = 0
#
# i.e. rows x lie in left kernel of H K^T.
# ------------------------------------------------------------

A = H*K.T

Cmod = A.left_kernel()
C = matrix(ZZ,Cmod.basis_matrix())

# Saturate just in case.
C = matrix(ZZ,C.row_module(ZZ).saturation().basis_matrix())

print()
print("========== COMPLEMENT ==========")
print("shape =",C.nrows(),C.ncols())
print("rank =",C.rank())

assert C.rank()==8
assert (C*H*K.T).is_zero()

GC = C*H*C.T

print("Gram:")
print(GC)
print("det =",GC.det())
print("positive definite =",GC.is_positive_definite())
print("even =",all(GC[i,i]%2==0 for i in range(8)))

DS = GC.smith_form()[0]
smith=[
    abs(DS[i,i])
    for i in range(8)
    if DS[i,i]!=0
]

print("Smith =",smith)

# Shell.
coeff=[]
for i in range(8):
    for j in range(i,8):
        coeff.append(
            GC[i,i]//2 if i==j else GC[i,j]
        )

Q=QuadraticForm(ZZ,8,coeff)
S=Q.short_vector_list_up_to_length(5,True)

print()
print("shell pairs:")
for q in range(1,len(S)):
    print(
        "norm",2*q,
        "pairs",len(S[q]),
        "signed",2*len(S[q])
    )

# Automorphisms.
L=IntegralLattice(GC)
Agrp=L.orthogonal_group(is_finite=True)

print()
print("Aut(C8) order =",Agrp.order())
print("Aut(C8) generators =",len(Agrp.gens()))

# ------------------------------------------------------------
# Index of K (+) C in full rank17 lattice.
# ------------------------------------------------------------

BC = K.stack(C)

print()
print("========== GLUING ==========")
print("combined rank =",BC.rank())
print("det coordinate matrix =",abs(BC.det()))

index=abs(BC.det())

print("[MW : K+C] =",index)

block_det=GK.det()*GC.det()

print("det(K)*det(C) =",block_det)
print("det(MW) =",H.det())

print(
    "check index^2 =",
    block_det//H.det()
)

assert index^2 * H.det() == block_det

# Quotient invariants of Z17 / (K+C).
SM=BC.smith_form()[0]

coord_smith=[
    abs(SM[i,i])
    for i in range(17)
    if SM[i,i]!=0
]

print("coordinate Smith =",coord_smith)

# ------------------------------------------------------------
# Minimal vectors split by projection membership.
# ------------------------------------------------------------

KQ=K.row_space(QQ)
CQ=C.row_space(QQ)

kin=[]
cin=[]
mixed=[]

for i,row in enumerate(signed):
    v=matrix(QQ,[row.tolist()])

    rkK=K.stack(v).rank()==9
    rkC=C.stack(v).rank()==8

    if rkK:
        kin.append(i)
    elif rkC:
        cin.append(i)
    else:
        mixed.append(i)

print()
print("minimal signed vectors:")
print(" in K9 =",len(kin))
print(" in C8 =",len(cin))
print(" mixed =",len(mixed))

# Save.
for name,M in [
    ("coxeter9_basis.txt",K),
    ("coxeter9_gram.txt",GK),
    ("complement8_basis.txt",C),
    ("complement8_gram.txt",GC),
]:
    (D/name).write_text(
        "\n".join(
            " ".join(map(str,row))
            for row in M.rows()
        )+"\n"
    )
