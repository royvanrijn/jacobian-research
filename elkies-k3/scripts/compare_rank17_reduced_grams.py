from pathlib import Path

import numpy as np
from sage.all import Matrix, RDF

BASE = Path(__file__).resolve().parents[1]

J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/rank29"

HG_np = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=float
)

HS_np = np.loadtxt(
    D / "E29-height-gram.txt",
    dtype=float
)

assert HG_np.shape == (17, 17)
assert HS_np.shape == (29, 29)

SUBSET = [
    0, 1, 4, 5, 6, 8, 9, 11, 15,
    16, 17, 19, 22, 23, 24, 25, 26
]

M_np = HS_np[np.ix_(SUBSET, SUBSET)]

# Determinant-forced specialization scale.
sign_m, logdet_m = np.linalg.slogdet(M_np)
sign_g, logdet_g = np.linalg.slogdet(HG_np)

assert sign_m > 0
assert sign_g > 0

h_det = np.exp((logdet_m - logdet_g) / 17.0)

print(f"HDET|h={h_det:.17g}", flush=True)

# Compare M against h * HG.
TARGET_np = h_det * HG_np

M = Matrix(RDF, M_np.tolist())
TARGET = Matrix(RDF, TARGET_np.tolist())

print("LLL|target=start", flush=True)
UG = TARGET.LLL_gram()
print("LLL|target=done", flush=True)

print("LLL|candidate=start", flush=True)
UM = M.LLL_gram()
print("LLL|candidate=done", flush=True)

RG = UG.transpose() * TARGET * UG
RM = UM.transpose() * M * UM

RG_np = np.array(RG, dtype=float)
RM_np = np.array(RM, dtype=float)

# Direct comparison of the two LLL-reduced Gram matrices.
R = RM_np - RG_np

rel_frob = np.linalg.norm(R, "fro") / np.linalg.norm(RG_np, "fro")
max_rel = np.max(np.abs(R)) / np.max(np.abs(RG_np))

print()
print("REDUCED")
print(f"rel_frob={rel_frob:.17g}")
print(f"max_rel={max_rel:.17g}")
print(f"max_UG={max(abs(int(x)) for x in UG.list())}")
print(f"max_UM={max(abs(int(x)) for x in UM.list())}")

# If LLL happened to choose corresponding bases, this matrix converts
# candidate reduced coordinates toward generic reduced coordinates.
print()
print("TARGET LLL TRANSFORM UG")
print(UG)

print()
print("CANDIDATE LLL TRANSFORM UM")
print(UM)

print()
print("REDUCED TARGET GRAM")
print(RG)

print()
print("REDUCED CANDIDATE GRAM")
print(RM)

print()
print("RESIDUAL")
print(Matrix(RDF, R.tolist()))

# Basis-independent-ish diagnostics from reduced basis vector lengths.
# Sorted diagonals are lengths of the LLL basis vectors.
diag_g = sorted(float(RG[i, i]) for i in range(17))
diag_m = sorted(float(RM[i, i]) for i in range(17))

print()
print("SORTED REDUCED BASIS LENGTHS")
for i, (g, m) in enumerate(zip(diag_g, diag_m)):
    print(
        f"L|i={i}"
        f"|generic={g:.17g}"
        f"|candidate={m:.17g}"
        f"|ratio={m/g:.17g}"
    )

# Also compare successive-minima-like shape after normalising
# each list by its shortest reduced-basis vector.
g0 = diag_g[0]
m0 = diag_m[0]

shape_err = np.sqrt(
    np.mean([
        (np.log((m / m0) / (g / g0))) ** 2
        for g, m in zip(diag_g, diag_m)
    ])
)

print()
print(f"REDUCED_SHAPE|log_rms={shape_err:.17g}")

np.savetxt(
    BASE / "results/rank17-reduced-target-gram.txt",
    RG_np,
    fmt="%.17g"
)

np.savetxt(
    BASE / "results/rank17-reduced-candidate-gram.txt",
    RM_np,
    fmt="%.17g"
)

np.savetxt(
    BASE / "results/rank17-reduced-residual.txt",
    R,
    fmt="%.17g"
)

with open(BASE / "results/rank17-target-LLL-U.txt", "w") as f:
    for row in UG.rows():
        f.write(" ".join(str(int(x)) for x in row) + "\n")

with open(BASE / "results/rank17-candidate-LLL-U.txt", "w") as f:
    for row in UM.rows():
        f.write(" ".join(str(int(x)) for x in row) + "\n")
