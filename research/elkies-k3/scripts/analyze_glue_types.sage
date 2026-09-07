from sage.all import *
from pathlib import Path
from collections import Counter, defaultdict
from math import gcd, lcm
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/lattice"
OUT = BASE / "data/glue"

def readmat(path):
    rows=[]
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        line=line.replace("["," ").replace("]"," ")
        try:
            row=[ZZ(x) for x in line.split()]
        except Exception:
            continue
        if row:
            rows.append(row)
    return matrix(ZZ,rows)

H = readmat(J/"short_vector_basis_gram.txt")
K = readmat(D/"coxeter9_basis.txt")
C = readmat(D/"complement8_basis.txt")

signed = np.load(
    J/"all_2622_signed_short_basis.npy"
)

BC = K.stack(C)

S,U,V = BC.smith_form()

diag=[abs(ZZ(S[i,i])) for i in range(17)]
nontrivial=[
    (i,d)
    for i,d in enumerate(diag)
    if d>1
]

mods=[d for _,d in nontrivial]

def label(v):
    w=vector(ZZ,v)*V
    return tuple(
        int(w[i] % d)
        for i,d in nontrivial
    )

def element_order(lab):
    o=1
    for x,m in zip(lab,mods):
        if x:
            o=lcm(o,m//gcd(x,m))
    return o

GK=K*H*K.T
GC=C*H*C.T

# ------------------------------------------------------------
# First collect vector-level information.
# ------------------------------------------------------------

cosets=defaultdict(list)

for idx,row in enumerate(signed):
    v=vector(QQ,map(QQ,row.tolist()))

    sol=BC.T.solve_right(v.column())
    coeff=vector(
        QQ,
        [sol[j,0] for j in range(17)]
    )

    ak=coeff[:9]
    ac=coeff[9:]

    nk=(ak*GK*ak.column())[0]
    nc=(ac*GC*ac.column())[0]

    assert nk+nc==4

    lab=label(row.tolist())

    cosets[lab].append(
        (idx,nk,nc)
    )

# ------------------------------------------------------------
# Each occupied coset should usually have one projection type;
# verify and summarize.
# ------------------------------------------------------------

type_cosets=Counter()
type_vectors=Counter()
type_occupancies=defaultdict(Counter)

for lab,items in cosets.items():
    order=element_order(lab)

    ptypes=Counter(
        (nk,nc)
        for _,nk,nc in items
    )

    # Keep separate if a coset has multiple projection types
    # (the identity coset is one known example).
    for (nk,nc),n in ptypes.items():
        key=(order,nk,nc)

        type_vectors[key]+=n
        type_occupancies[key][n]+=1

    # Count coset once for each projection type represented.
    for ptype in ptypes:
        type_cosets[(order,*ptype)] += 1

print(
    "order | K norm | C norm | cosets | signed | pairs-ish | occupancy distribution"
)
print("-"*110)

grand=0

for key in sorted(
    type_vectors,
    key=lambda k:(k[0],k[1],k[2])
):
    order,nk,nc=key
    nv=type_vectors[key]
    ncsets=type_cosets[key]

    occ=", ".join(
        f"{mult}x{count}"
        for mult,count in sorted(
            type_occupancies[key].items()
        )
    )

    print(
        f"{order:5d} | "
        f"{str(nk):>6} | "
        f"{str(nc):>6} | "
        f"{ncsets:6d} | "
        f"{nv:6d} | "
        f"{float(nv)/2:9.1f} | "
        f"{occ}"
    )

    grand += nv

print("-"*110)
print("TOTAL SIGNED =",grand)
print("TOTAL PAIRS  =",grand//2)

# ------------------------------------------------------------
# Aggregate by projection type ignoring order.
# ------------------------------------------------------------

print()
print("========== PROJECTION TYPE BY ORDERS ==========")

by_projection=defaultdict(Counter)

for (order,nk,nc),nv in type_vectors.items():
    by_projection[(nk,nc)][order]+=nv

for (nk,nc),orders in sorted(
    by_projection.items(),
    key=lambda x:(x[0][0],x[0][1])
):
    print(
        f"K={str(nk):>6} C={str(nc):>6}"
        f" total={sum(orders.values()):4d}"
        f" orders={dict(sorted(orders.items()))}"
    )

# ------------------------------------------------------------
# Save compact table.
# ------------------------------------------------------------

with (OUT/"glue_type_table.txt").open("w") as f:
    f.write(
        "order K_norm C_norm cosets signed occupancy\n"
    )

    for key in sorted(
        type_vectors,
        key=lambda k:(k[0],k[1],k[2])
    ):
        order,nk,nc=key

        occ=",".join(
            f"{mult}:{count}"
            for mult,count in sorted(
                type_occupancies[key].items()
            )
        )

        f.write(
            f"{order} {nk} {nc} "
            f"{type_cosets[key]} "
            f"{type_vectors[key]} "
            f"{occ}\n"
        )
