from sage.all import *
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/relations"
R = BASE / "results"

signed = np.load(J / "all_2622_signed_short_basis.npy")

motif = [
    int(x)
    for x in (D / "best_relation_motif_10_indices.txt").read_text().split()
]

seed = [1337, 736, 2474]

H = matrix(
    ZZ,
    np.loadtxt(
        J / "short_vector_basis_gram.txt",
        dtype=np.int64
    ).tolist()
)

V = {
    i: vector(ZZ, signed[i].tolist())
    for i in motif
}

B = matrix(ZZ, [V[i] for i in seed])

print("seed =", seed)
print("rank =", B.rank())

G = B * H * B.T

print()
print("SEED GRAM")
print(G)
print("det =", G.det())

# Smith data tells us whether this rank-3 row lattice is primitive
# in Z^17.
S = B.smith_form()[0]

smith = [
    abs(S[i,i])
    for i in range(min(S.nrows(), S.ncols()))
    if S[i,i] != 0
]

print()
print("coordinate Smith invariants =", smith)
print("primitive row lattice =", all(x == 1 for x in smith))

# ------------------------------------------------------------
# Solve every motif vector as c0*P + c1*Q + c2*R.
#
# Since B has rows P,Q,R, solve B^T*c = v^T.
# ------------------------------------------------------------

print()
print("MOTIF IN P,Q,R COORDINATES")

coords = {}

for idx in motif:
    v = V[idx]

    sol = B.T.solve_right(v.column())
    c = vector(QQ, sol)

    coords[idx] = c

    print(
        f"{idx:4d} = "
        f"({c[0]}) P + ({c[1]}) Q + ({c[2]}) R"
    )

# Verify.
for idx,c in coords.items():
    assert c * B == V[idx]

# Integral?
integral = all(
    x.denominator() == 1
    for c in coords.values()
    for x in c
)

print()
print("all coordinates integral =", integral)

# ------------------------------------------------------------
# Read and rewrite internal relations in P,Q,R coordinates.
# ------------------------------------------------------------

triples = np.load(J / "minimal_additive_triples.npy")
M = set(motif)

internal = [
    tuple(map(int,t))
    for t in triples
    if all(int(x) in M for x in t)
]

print()
print("RELATIONS IN P,Q,R COORDINATES")

for a,b,c in internal:
    ca,cb,cc = coords[a],coords[b],coords[c]

    assert ca + cb == cc

    print(
        f"{a} + {b} = {c}    "
        f"{tuple(ca)} + {tuple(cb)} = {tuple(cc)}"
    )

# ------------------------------------------------------------
# Save machine-readable simple coordinates.
# ------------------------------------------------------------

out = D / "rank3_motif_coordinates.txt"

out.write_text(
    "\n".join(
        f"{idx} " +
        " ".join(str(x) for x in coords[idx])
        for idx in motif
    ) + "\n"
)

(D / "rank3_motif_seed_indices.txt").write_text(
    " ".join(map(str,seed)) + "\n"
)

(D / "rank3_motif_seed_gram.txt").write_text(
    "\n".join(
        " ".join(map(str,row))
        for row in G.rows()
    ) + "\n"
)

print()
print("saved", out)
