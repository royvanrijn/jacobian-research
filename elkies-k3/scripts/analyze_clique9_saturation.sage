from sage.all import *
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "data/relations"

signed = np.load(J / "all_2622_signed_short_basis.npy")

H = matrix(
    ZZ,
    np.loadtxt(
        J / "short_vector_basis_gram.txt",
        dtype=np.int64
    ).tolist()
)

clique = [
    2127, 2585, 2431, 2375, 2051,
    1987, 2115, 1705, 1483
]

B = matrix(ZZ, signed[clique].tolist())

print("raw rank =", B.rank())

Graw = B * H * B.T

print()
print("RAW CLIQUE")
print("Gram:")
print(Graw)
print("det =", Graw.det())

# ------------------------------------------------------------
# Saturate the row lattice B inside Z^17.
#
# Row-space saturation = integer points in the rational row span.
# Sage's saturation() is available on the row module.
# ------------------------------------------------------------

M = B.row_module(ZZ)

print()
print("raw module index data:")
print(M)

S = M.saturation()
Bs = matrix(ZZ, S.basis_matrix())

print()
print("SATURATION")
print("basis shape =", Bs.nrows(), Bs.ncols())
print("rank =", Bs.rank())

Gs = Bs * H * Bs.T

print("Gram:")
print(Gs)

print("det =", Gs.det())

ratio = ZZ(Graw.det() // Gs.det())

print("det ratio =", ratio)
print("index =", sqrt(ratio))

assert Gs.det() == 1280
assert ratio == 4

# ------------------------------------------------------------
# Smith invariants of the saturated rank-9 Gram.
# ------------------------------------------------------------

D = Gs.smith_form()[0]

smith = [
    abs(D[i,i])
    for i in range(9)
    if D[i,i] != 0
]

print()
print("Gram Smith invariants:")
print(smith)

# ------------------------------------------------------------
# Shell of saturated rank-9 lattice.
# ------------------------------------------------------------

coeffs = []

for i in range(9):
    for j in range(i,9):
        coeffs.append(
            Gs[i,i]//2 if i == j
            else Gs[i,j]
        )

Q = QuadraticForm(ZZ,9,coeffs)

shell = Q.short_vector_list_up_to_length(3, True)

print()
print("norm-2 +/- pairs =", len(shell[1]))
print("norm-4 +/- pairs =", len(shell[2]))
print("norm-4 signed =", 2*len(shell[2]))

# ------------------------------------------------------------
# Exact automorphism group of this rank-9 sublattice.
# ------------------------------------------------------------

L = IntegralLattice(Gs)

print()
print("computing orthogonal group ...")

A = L.orthogonal_group(is_finite=True)

print("automorphism group order =", A.order())
print("automorphism generators =", len(A.gens()))

# ------------------------------------------------------------
# Find the actual glue vector:
#
# vector in saturation but not raw row lattice.
# ------------------------------------------------------------

print()
print("GLUE")

raw_Q = B.row_space(QQ)

for row in Bs.rows():
    v = vector(ZZ,row)

    # Try expressing v in the raw basis.
    c = B.T.solve_right(v.column())
    c = vector(QQ,[c[i,0] for i in range(c.nrows())])

    if any(x.denominator() != 1 for x in c):
        print("saturated basis vector outside raw lattice:")
        print(v)
        print("coordinates in 2A9 clique basis:")
        print(c)

        print("2*v clique coordinates:")
        print(2*c)

# ------------------------------------------------------------
# Save canonical saturated basis / Gram.
# ------------------------------------------------------------

OUT.mkdir(parents=True,exist_ok=True)

(OUT/"clique9_saturated_basis.txt").write_text(
    "\n".join(
        " ".join(map(str,row))
        for row in Bs.rows()
    ) + "\n"
)

(OUT/"clique9_saturated_gram.txt").write_text(
    "\n".join(
        " ".join(map(str,row))
        for row in Gs.rows()
    ) + "\n"
)
