from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "data/k3-model"
OUT.mkdir(parents=True, exist_ok=True)

H = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

signed = np.load(
    J / "all_2622_signed_short_basis.npy"
)

lookup = {
    tuple(map(int,row)): i
    for i,row in enumerate(signed)
}

# In the short basis, the actual MW basis vectors are e_0,...,e_16.
basis_idx = {}

for i in range(17):
    e = [0]*17
    e[i] = 1
    basis_idx[i] = lookup[tuple(e)]

print("BASIS SIGNED INDICES")
for i in range(17):
    print(i, basis_idx[i])

edges=[]

for i in range(17):
    for j in range(i):
        if H[i,j] != -2:
            continue

        e = [0]*17
        e[i] += 1
        e[j] += 1

        k = lookup.get(tuple(e))

        if k is None:
            raise RuntimeError(
                f"missing minimal sum e_{i}+e_{j}"
            )

        edges.append((i,j,basis_idx[i],basis_idx[j],k))

print()
print("-2 ADDITION TRIANGLES")

for i,j,pi,pj,s in edges:
    print(
        f"P{i} + P{j} = S{s}"
        f"   signed=({pi},{pj},{s})"
    )

print()
print("edges =",len(edges))

# Hub statistics
from collections import Counter
degree=Counter()

for i,j,*_ in edges:
    degree[i]+=1
    degree[j]+=1

print()
print("HUBS")
for i,d in degree.most_common():
    print(i,d)

with (OUT/"minus2_addition_triangles.tsv").open("w") as f:
    f.write(
        "i\tj\tPi_signed\tPj_signed\tsum_signed\n"
    )

    for row in edges:
        f.write("\t".join(map(str,row))+"\n")
