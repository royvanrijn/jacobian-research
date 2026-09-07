from pathlib import Path
from collections import deque
import numpy as np

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

signed = np.load(J / "all_2622_signed_short_basis.npy")
triples = np.load(J / "minimal_additive_triples.npy")

N = len(signed)
M = len(triples)

# ------------------------------------------------------------
# Incidence structure
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Locate +/- standard basis vectors.
#
# all_1311_in_short_basis.txt is expressed in the minimal-vector
# unimodular basis, so the MW basis vectors are e_1,...,e_17.
# ------------------------------------------------------------

lookup = {
    tuple(map(int,row)): i
    for i,row in enumerate(signed)
}

basis_plus = []
basis_minus = []

for i in range(17):
    e = [0] * 17
    e[i] = 1

    em = [0] * 17
    em[i] = -1

    basis_plus.append(lookup[tuple(e)])
    basis_minus.append(lookup[tuple(em)])

print("basis + indices:", basis_plus)
print("basis - indices:", basis_minus)

C = closure(basis_plus)

print()
print("17-BASIS CLOSURE:", len(C), "of", N)

# Signs of sections are available geometrically once P is known,
# so test explicitly seeding both +/- versions too.
Cpm = closure(basis_plus + basis_minus)

print("17-BASIS +/- CLOSURE:", len(Cpm), "of", N)

# ------------------------------------------------------------
# How fast does closure grow as basis elements are introduced?
# ------------------------------------------------------------

print()
print("PREFIX CLOSURES")

for k in range(2,18):
    c = closure(basis_plus[:k])
    print(f"k={k:2d} closure={len(c):4d}")

# ------------------------------------------------------------
# Greedily reorder the 17 actual basis vectors to maximize
# closure after each addition.
# ------------------------------------------------------------

remaining = set(basis_plus)
chosen = []

print()
print("GREEDY BASIS ORDER")

while remaining:
    best = None

    for v in remaining:
        c = closure(chosen + [v])
        score = len(c)

        if best is None or score > best[0]:
            best = (score,v)

    score,v = best

    chosen.append(v)
    remaining.remove(v)

    print(
        f"k={len(chosen):2d}"
        f" add={v:4d}"
        f" closure={score:4d}"
    )

# ------------------------------------------------------------
# Start from previous relation-rich 10-set, then ask which
# single vectors give the largest closure expansion.
# ------------------------------------------------------------

seed_path = J / "relation_seed_indices.txt"

if seed_path.exists():
    old_seed = [
        int(x)
        for x in seed_path.read_text().split()
    ]

    print()
    print("OLD 10-SEED:", old_seed)
    print("closure:", len(closure(old_seed)))

    seed = list(old_seed)
    remaining = set(range(N)) - set(seed)

    print()
    print("EXTENDING OLD SEED BY ACTUAL CLOSURE")

    for step in range(10):
        best = None

        # Evaluating all 2600 choices is now cheap enough.
        for v in remaining:
            c = closure(seed + [v])
            score = len(c)

            if best is None or score > best[0]:
                best = (score,v)

        score,v = best
        seed.append(v)
        remaining.remove(v)

        print(
            f"k={len(seed):2d}"
            f" add={v:4d}"
            f" closure={score:4d}"
        )

        if score == N:
            break

    (J / "closure_optimized_seed_indices.txt").write_text(
        " ".join(map(str,seed)) + "\n"
    )

    final = closure(seed)

    print()
    print("FINAL OPTIMIZED CLOSURE:",len(final),"of",N)
