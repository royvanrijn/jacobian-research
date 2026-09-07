from sage.all import *
from pathlib import Path
import json
import time

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "data/automorphisms"

OUT.mkdir(parents=True, exist_ok=True)

def readmat(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        line = line.replace("[", " ").replace("]", " ")
        try:
            row = [ZZ(x) for x in line.split()]
        except Exception:
            continue
        if row:
            rows.append(row)
    return matrix(ZZ, rows)

H = readmat(J / "short_vector_basis_gram.txt")

print("rank =", H.nrows(), flush=True)
print("det =", H.det(), flush=True)
print("positive definite =", H.is_positive_definite(), flush=True)

L = IntegralLattice(H)

print("AUT|start", flush=True)
t0 = time.time()

G = L.orthogonal_group(is_finite=True)

print("AUT|constructed|seconds=%.3f" % (time.time() - t0), flush=True)

try:
    order = G.order()
except Exception as e:
    print("AUT|order_failed|", repr(e), flush=True)
    order = None

gens = list(G.gens())

print("AUT|generators =", len(gens), flush=True)
print("AUT|order =", order, flush=True)

for i, g in enumerate(gens):
    print()
    print("generator", i)
    print(matrix(QQ, g))

# Save generators as integral matrices.
for i, g in enumerate(gens):
    M = matrix(QQ, g)

    if not all(x.denominator() == 1 for x in M.list()):
        raise RuntimeError("non-integral lattice automorphism")

    M = matrix(ZZ, M)

    assert M * H * M.transpose() == H

    (OUT / f"generator-{i:02d}.txt").write_text(
        "\n".join(
            " ".join(map(str, row))
            for row in M.rows()
        ) + "\n"
    )

summary = {
    "rank": int(H.nrows()),
    "determinant": int(H.det()),
    "number_of_generators": len(gens),
    "order": None if order is None else int(order),
}

(OUT / "summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)

print()
print("AUT|done")
