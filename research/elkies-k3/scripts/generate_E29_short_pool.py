from pathlib import Path
import argparse
import heapq
import random
import numpy as np
from sage.all import Matrix, RDF

parser = argparse.ArgumentParser()
parser.add_argument("gram")
parser.add_argument("--samples", type=int, default=5_000_000)
parser.add_argument("--keep", type=int, default=200_000)
parser.add_argument("--seed", type=int, default=290017)
parser.add_argument("--cutoff", type=float, default=300.0)
args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results"

HA = np.loadtxt(args.gram, dtype=float)
r = HA.shape[0]

G = Matrix(RDF, HA.tolist())

print("LLL start", flush=True)
U = G.LLL_gram()
print("LLL done", flush=True)

R = np.array(
    U.transpose() * G * U,
    dtype=float
)

U_np = np.array(U, dtype=np.int64)

rng = random.Random(args.seed)

# canonical vector -> norm
best = {}

def canonical(z):
    nz = np.flatnonzero(z)

    if not len(nz):
        return None

    z = z.copy()

    if z[nz[0]] < 0:
        z = -z

    return tuple(map(int, z))

def consider(z):
    key = canonical(z)

    if key is None or key in best:
        return

    z = np.asarray(key, dtype=np.int64)

    n = float(z @ R @ z)

    if n <= args.cutoff:
        best[key] = n


# ------------------------------------------------------------
# Deterministic seed pool
# ------------------------------------------------------------

for i in range(r):
    z = np.zeros(r, dtype=np.int64)
    z[i] = 1
    consider(z)

for i in range(r):
    for j in range(i + 1, r):
        for s in (-1, 1):
            z = np.zeros(r, dtype=np.int64)
            z[i] = 1
            z[j] = s
            consider(z)


seed_vectors = [
    np.array(k, dtype=np.int64)
    for k in best
]


# ------------------------------------------------------------
# Random local exploration
# ------------------------------------------------------------

current = seed_vectors[
    rng.randrange(len(seed_vectors))
].copy()

current_norm = float(current @ R @ current)

for step in range(args.samples):

    if rng.random() < 0.02:
        # restart from something already known
        keys = list(best.keys())

        current = np.array(
            keys[rng.randrange(len(keys))],
            dtype=np.int64
        )

        current_norm = float(
            current @ R @ current
        )

    p = current.copy()

    move = rng.randrange(100)

    if move < 70:
        i = rng.randrange(r)
        p[i] += rng.choice((-1, 1))

    elif move < 90:
        i = rng.randrange(r)
        p[i] += rng.choice((-2, 2))

    else:
        i, j = rng.sample(range(r), 2)
        p[i] += rng.choice((-1, 1)) * p[j]

    # Don't allow coefficient explosion.
    if np.max(np.abs(p)) > 12:
        continue

    pn = float(p @ R @ p)

    if pn <= args.cutoff:
        consider(p)

        # Bias walk toward short vectors but allow wandering.
        delta = pn - current_norm

        if (
            delta <= 0
            or rng.random() < np.exp(-delta / 20.0)
        ):
            current = p
            current_norm = pn

    if (step + 1) % 100000 == 0:
        print(
            f"STEP|n={step+1}"
            f"|unique={len(best)}"
            f"|current_norm={current_norm:.6g}",
            flush=True
        )


# ------------------------------------------------------------
# Keep shortest N
# ------------------------------------------------------------

items = sorted(
    best.items(),
    key=lambda x: x[1]
)[:args.keep]

Z = np.array(
    [
        np.array(k, dtype=np.int64)
        for k, n in items
    ],
    dtype=np.int64
)

N = np.array(
    [n for k, n in items],
    dtype=float
)

# Convert back to original E29 coordinates.
V = Z @ U_np.T

np.save(
    OUT / "E29-short-pool-vectors.npy",
    V
)

np.save(
    OUT / "E29-short-pool-norms.npy",
    N
)

print()
print("DONE")
print("unique =", len(best))
print("saved =", len(V))
print("min norm =", N.min())
print("max norm =", N.max())
