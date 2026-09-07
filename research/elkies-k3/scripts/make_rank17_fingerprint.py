from pathlib import Path
import numpy as np
from sage.all import Matrix, RDF

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

H = np.loadtxt(J / "short_vector_basis_gram.txt", dtype=float)
G = Matrix(RDF, H.tolist())

U = G.LLL_gram()
R = U.transpose() * G * U
R = np.array(R, dtype=float)

# Our reduced Gram appears to have a common unit:
# diagonal ~= 4 * unit.
unit = np.median(np.diag(R)) / 4.0

Q = np.rint(R / unit).astype(int)

err = np.max(np.abs(R / unit - Q))

print("unit =", unit)
print("rounding_error =", err)
print()
print(Q)

assert err < 1e-8

np.savetxt(
    BASE / "results/rank17-fingerprint-Q.txt",
    Q,
    fmt="%d"
)
