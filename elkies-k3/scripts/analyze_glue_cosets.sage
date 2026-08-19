from sage.all import *
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/lattice"
OUT = BASE / "data/glue"

OUT.mkdir(parents=True, exist_ok=True)

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

assert BC.nrows() == 17
assert BC.ncols() == 17
assert abs(BC.det()) == 640
assert (K*H*C.T).is_zero()

# ------------------------------------------------------------
# Smith normal form
#
# U * BC * V = S
#
# Right multiplication by V is a unimodular change of ambient
# coordinates. Thus the row lattice becomes the diagonal row
# lattice S.
# ------------------------------------------------------------

snf = BC.smith_form()

print("smith_form tuple length =",len(snf))

S = snf[0]
U = snf[1]
V = snf[2]

assert U*BC*V == S

diag = [
    abs(ZZ(S[i,i]))
    for i in range(17)
]

print("SNF diagonal:")
print(diag)

nontrivial = [
    (i,d)
    for i,d in enumerate(diag)
    if d > 1
]

print("nontrivial quotient coordinates:")
print(nontrivial)

print(
    "quotient order =",
    prod(d for _,d in nontrivial)
)

# ------------------------------------------------------------
# Canonical quotient label.
# ------------------------------------------------------------

def coset_label(v):
    w = vector(ZZ,v) * V

    return tuple(
        int(w[i] % d)
        for i,d in nontrivial
    )

labels = defaultdict(list)

for i,row in enumerate(signed):
    lab = coset_label(row.tolist())
    labels[lab].append(i)

print()
print("========== MINIMAL VECTOR COSET DISTRIBUTION ==========")

print("occupied glue cosets =",len(labels))
print("empty glue cosets =",640-len(labels))

hist = Counter(
    len(indices)
    for indices in labels.values()
)

print("signed vectors per occupied coset:")
for n,c in sorted(hist.items()):
    print(
        f"  {n:3d} vectors : {c:3d} cosets"
    )

assert sum(
    n*c for n,c in hist.items()
) == 2622

# ------------------------------------------------------------
# Pair cosets under x -> -x.
# ------------------------------------------------------------

def neg_label(lab):
    return tuple(
        (-x) % d
        for x,(_,d) in zip(lab,nontrivial)
    )

seen=set()
coset_pairs=[]

for lab in sorted(labels):
    if lab in seen:
        continue

    neg=neg_label(lab)

    seen.add(lab)
    seen.add(neg)

    coset_pairs.append(
        (
            lab,
            neg,
            len(labels.get(lab,[])),
            len(labels.get(neg,[]))
        )
    )

self_inverse = [
    x
    for x in coset_pairs
    if x[0] == x[1]
]

print()
print("occupied cosets modulo +/- =",len(coset_pairs))
print("self-inverse occupied cosets =",len(self_inverse))

# ------------------------------------------------------------
# Orthogonal projections onto K and C.
#
# Since BC gives a rational basis of MW tensor Q:
#
#      v = a*K + b*C.
#
# K and C are orthogonal, so we can compute the two projection
# norms separately.
# ------------------------------------------------------------

GK = K*H*K.T
GC = C*H*C.T

projection_data=[]

proj_hist=Counter()
coset_proj=defaultdict(Counter)

for i,row in enumerate(signed):
    v=vector(QQ,map(QQ,row.tolist()))

    sol = BC.T.solve_right(v.column())
    coeff = vector(
        QQ,
        [sol[j,0] for j in range(17)]
    )

    ak=coeff[:9]
    ac=coeff[9:]

    nk = (ak * GK * ak.column())[0]
    nc = (ac * GC * ac.column())[0]

    assert nk+nc == 4

    lab=coset_label(row.tolist())

    projection_data.append(
        (i,lab,nk,nc)
    )

    proj_hist[(nk,nc)] += 1
    coset_proj[lab][(nk,nc)] += 1

print()
print("========== PROJECTION NORM DISTRIBUTION ==========")

for (nk,nc),count in sorted(
    proj_hist.items(),
    key=lambda x:(x[0][0],x[0][1])
):
    print(
        f"Knorm={nk!s:>8}"
        f"  Cnorm={nc!s:>8}"
        f"  count={count}"
    )

print()
print(
    "number of distinct projection norm pairs =",
    len(proj_hist)
)

# ------------------------------------------------------------
# Largest / most populated glue cosets.
# ------------------------------------------------------------

print()
print("========== MOST POPULATED COSETS ==========")

ordered=sorted(
    labels.items(),
    key=lambda x:(-len(x[1]),x[0])
)

for lab,inds in ordered[:40]:
    ph=coset_proj[lab]

    print(
        "coset",lab,
        "count",len(inds),
        "projection_types",
        sorted(
            (
                str(k),
                str(c),
                n
            )
            for (k,c),n in ph.items()
        )
    )

# ------------------------------------------------------------
# Save text tables.
# ------------------------------------------------------------

(OUT/"quotient_snf.txt").write_text(
    " ".join(map(str,diag))+"\n"
)

with (OUT/"minimal_vector_cosets.txt").open("w") as f:
    for lab in sorted(labels):
        f.write(
            " ".join(map(str,lab))
            + " : "
            + " ".join(map(str,labels[lab]))
            + "\n"
        )

with (OUT/"minimal_coset_histogram.txt").open("w") as f:
    for n,c in sorted(hist.items()):
        f.write(f"{n} {c}\n")

with (OUT/"projection_norm_histogram.txt").open("w") as f:
    for (nk,nc),count in sorted(
        proj_hist.items(),
        key=lambda x:(x[0][0],x[0][1])
    ):
        f.write(
            f"{nk} {nc} {count}\n"
        )

print()
print("saved under",OUT)
