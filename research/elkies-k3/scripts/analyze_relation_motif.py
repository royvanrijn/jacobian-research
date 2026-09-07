from pathlib import Path
from collections import defaultdict, deque
import numpy as np
from sage.all import ZZ, matrix

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/relations"

signed = np.load(J / "all_2622_signed_short_basis.npy")
triples = np.load(J / "minimal_additive_triples.npy")

# use final motif-10 result
idx_path = D / "best_relation_motif_10_indices.txt"
motif = [int(x) for x in idx_path.read_text().split()]

S = set(motif)

internal = []
for a,b,c in triples:
    a,b,c = map(int,(a,b,c))
    if a in S and b in S and c in S:
        internal.append((a,b,c))

print("motif =", motif)
print("size =", len(motif))

B = matrix(ZZ, signed[motif].tolist())
print("MW rank =", B.rank())

print()
print("INTERNAL RELATIONS")
for t in internal:
    print(t)

print()
print("internal count =", len(internal))

# ------------------------------------------------------------
# Can some motif vectors be derived from others purely using
# the internal relations?
# ------------------------------------------------------------

def closure(seed):
    known = set(seed)
    changed = True

    while changed:
        changed = False

        for a,b,c in internal:
            q = [a in known, b in known, c in known]

            if sum(q) >= 2:
                for x in (a,b,c):
                    if x not in known:
                        known.add(x)
                        changed = True

    return known

best = None

for mask in range(1,1 << len(motif)):
    seed = [
        motif[i]
        for i in range(len(motif))
        if (mask >> i) & 1
    ]

    C = closure(seed)

    if len(C) == len(motif):
        if best is None or len(seed) < len(best):
            best = seed

print()
print("MINIMUM INTERNAL GENERATING SEED")
print(best)
print("size =", len(best) if best else None)

# ------------------------------------------------------------
# For each non-seed vector, find a relation that derives it
# once possible.
# ------------------------------------------------------------

if best:
    known = set(best)
    derivations = []

    while len(known) < len(motif):
        progress = False

        for a,b,c in internal:
            arr = [a,b,c]
            present = [x in known for x in arr]

            if sum(present) == 2:
                missing = arr[present.index(False)]

                if missing not in known:
                    known.add(missing)
                    derivations.append((missing,a,b,c))
                    progress = True

        if not progress:
            break

    print()
    print("DERIVATION ORDER")

    for missing,a,b,c in derivations:
        print(
            f"{missing} derived from relation "
            f"{a} + {b} = {c}"
        )

    (D / "motif10_minimal_seed.txt").write_text(
        " ".join(map(str,best)) + "\n"
    )

    (D / "motif10_derivations.txt").write_text(
        "\n".join(
            f"{m} {a} {b} {c}"
            for m,a,b,c in derivations
        ) + "\n"
    )
