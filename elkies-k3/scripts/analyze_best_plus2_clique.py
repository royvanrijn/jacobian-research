from pathlib import Path
from collections import deque
import numpy as np
from sage.all import ZZ, matrix

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

signed = np.load(J / "all_2622_signed_short_basis.npy")
triples = np.load(J / "minimal_additive_triples.npy")

clique = [2127, 2585, 2431, 2375, 2051, 1987, 2115, 1705, 1483]

lookup = {
    tuple(map(int,row)): i
    for i,row in enumerate(signed)
}

print("clique size =", len(clique))

# rank
B = matrix(ZZ, signed[clique].tolist())
print("MW rank of clique =", B.rank())

# pairwise differences
diffs = set()

for a in range(len(clique)):
    for b in range(a):
        v = signed[clique[a]] - signed[clique[b]]
        k = lookup.get(tuple(map(int,v)))

        if k is None:
            raise RuntimeError("difference is not minimal")

        diffs.add(k)

print("pairwise differences =", len(diffs))
print("original + differences =", len(set(clique) | diffs))

# incidence structure for closure
N = len(signed)
M = len(triples)

flat = triples.reshape(-1)
tids = np.repeat(np.arange(M, dtype=np.int32), 3)

order = np.argsort(flat, kind="stable")
sv = flat[order]
si = tids[order]

offsets = np.searchsorted(
    sv,
    np.arange(N + 1, dtype=np.int32)
)

def incident(v):
    return si[offsets[v]:offsets[v+1]]

def closure(seed):
    known = np.zeros(N, dtype=bool)
    known[list(seed)] = True

    count = np.count_nonzero(
        known[triples],
        axis=1
    ).astype(np.uint8)

    ready = np.flatnonzero(count >= 2)

    queued = np.zeros(M, dtype=bool)
    queued[ready] = True

    q = deque(map(int, ready))

    while q:
        tid = q.popleft()

        for v in triples[tid]:
            v = int(v)

            if known[v]:
                continue

            known[v] = True

            ids = incident(v)
            old = count[ids].copy()
            count[ids] = old + 1

            nr = ids[
                (old < 2)
                & (count[ids] >= 2)
                & (~queued[ids])
            ]

            queued[nr] = True
            q.extend(map(int, nr))

    return np.flatnonzero(known)

C1 = closure(clique)
C2 = closure(set(clique) | diffs)

print("closure(clique) =", len(C1))
print("closure(clique + explicit differences) =", len(C2))

print()
print("difference indices:")
print(sorted(diffs))
