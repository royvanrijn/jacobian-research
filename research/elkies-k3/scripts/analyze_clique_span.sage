from sage.all import *
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

signed = np.load(J / "all_2622_signed_short_basis.npy")

clique = [2127, 2585, 2431, 2375, 2051, 1987, 2115, 1705, 1483]

B = matrix(QQ, signed[clique].tolist())
print("rank =", B.rank())

span_indices = []

for i,row in enumerate(signed):
    v = vector(QQ, row.tolist())

    if B.stack(matrix(QQ,[v])).rank() == B.rank():
        span_indices.append(i)

print("minimal signed vectors in rank-9 span =", len(span_indices))
print("pairs =", len(span_indices)//2)

# Gram and determinant of the clique lattice
H = matrix(
    ZZ,
    np.loadtxt(J / "short_vector_basis_gram.txt", dtype=np.int64).tolist()
)

BZ = matrix(ZZ, signed[clique].tolist())
G = BZ * H * BZ.T

print()
print("clique Gram:")
print(G)

print()
print("det =", G.det())

S = BZ.smith_form()[0]
smith = [
    abs(S[i,i])
    for i in range(min(S.nrows(),S.ncols()))
    if S[i,i] != 0
]

print("coordinate Smith =", smith)
print("primitive =", all(x == 1 for x in smith))

(J / "clique9_span_minimal_indices.txt").write_text(
    " ".join(map(str,span_indices)) + "\n"
)
