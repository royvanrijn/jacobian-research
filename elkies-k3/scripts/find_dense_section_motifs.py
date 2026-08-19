from pathlib import Path
import argparse
import random
import time
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--restarts", type=int, default=20000)
parser.add_argument("--seed", type=int, default=12001)
parser.add_argument("--motif-size", type=int, default=10)
args = parser.parse_args()

random.seed(args.seed)
rng = np.random.default_rng(args.seed)

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "data/relations"

OUT.mkdir(parents=True, exist_ok=True)

signed = np.load(J / "all_2622_signed_short_basis.npy")
triples = np.load(J / "minimal_additive_triples.npy")

# Gram in the same short basis as `signed`.
H = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64,
)

N = len(signed)
M = len(triples)

print("vectors =", N)
print("triples =", M)

# ============================================================
# Pairing matrix
# ============================================================

t0 = time.time()

P = (signed @ H) @ signed.T

assert np.all(np.diag(P) == 4)

print("pairing matrix seconds =", round(time.time() - t0, 3))

# ============================================================
# +2 graph.
#
# <a,b> = 2 => a-b is again minimal.
# ============================================================

plus2 = (P == 2)
np.fill_diagonal(plus2, False)

degrees = plus2.sum(axis=1)

print()
print("+2 degree:")
print(" min =", int(degrees.min()))
print(" max =", int(degrees.max()))
print(" mean =", float(degrees.mean()))

# ------------------------------------------------------------
# Greedy randomized maximal clique.
# ------------------------------------------------------------

def greedy_clique(start=None):
    if start is None:
        start = int(rng.integers(N))

    clique = [start]

    candidates = np.flatnonzero(plus2[start])

    while len(candidates):
        # Within current common-neighbor set, choose a vertex
        # with highest induced degree, with randomness among
        # the best handful.
        sub = plus2[np.ix_(candidates, candidates)]
        d = sub.sum(axis=1)

        order = np.argsort(-d)

        pick_window = min(5, len(order))
        k = int(order[int(rng.integers(pick_window))])
        v = int(candidates[k])

        clique.append(v)

        candidates = candidates[
            plus2[v, candidates]
        ]

    return clique


best_clique = []

# Seed a fraction of searches at high-degree vertices.
high = np.argsort(-degrees)[:200]

for r in range(args.restarts):

    if r < len(high):
        start = int(high[r])
    else:
        start = None

    C = greedy_clique(start)

    if len(C) > len(best_clique):
        best_clique = C

        print(
            f"CLIQUE|new_best"
            f"|restart={r}"
            f"|size={len(C)}"
            f"|vertices={C}",
            flush=True,
        )

print()
print("BEST +2 CLIQUE SIZE =", len(best_clique))
print("indices =", best_clique)

GC = P[np.ix_(best_clique, best_clique)]

print("clique Gram:")
print(GC)

# ============================================================
# Additive triple incidence
# ============================================================

incident = [[] for _ in range(N)]

for tid, (a,b,c) in enumerate(triples):
    incident[int(a)].append(tid)
    incident[int(b)].append(tid)
    incident[int(c)].append(tid)

incident = [
    np.asarray(x, dtype=np.int32)
    for x in incident
]

# Number of fully internal additive relations.
def internal_count(seed):
    inside = np.zeros(N, dtype=np.uint8)
    inside[seed] = 1

    return int(
        np.count_nonzero(
            inside[triples].sum(axis=1) == 3
        )
    )


# More useful score:
#
#   internal equations * large weight
#   + equations with exactly two known vertices
#
def motif_score(seed):
    inside = np.zeros(N, dtype=np.uint8)
    inside[seed] = 1

    counts = inside[triples].sum(axis=1)

    internal = int(np.count_nonzero(counts == 3))
    almost = int(np.count_nonzero(counts == 2))

    return internal, almost


# ------------------------------------------------------------
# Randomized local search for fixed-size relation-rich motif.
# ------------------------------------------------------------

K = args.motif_size

# Start from best +2 clique when possible.
if len(best_clique) >= K:
    current = list(best_clique[:K])
else:
    current = list(best_clique)

    remaining = [
        i for i in np.argsort(-degrees)
        if i not in current
    ]

    current += list(map(int, remaining[:K-len(current)]))

current_score = motif_score(current)
best_motif = list(current)
best_score = current_score

print()
print(
    "MOTIF|start",
    "size", K,
    "score", current_score,
    "vertices", current,
    flush=True,
)

temperature = 5.0

for step in range(200000):

    pos = random.randrange(K)

    candidate = random.randrange(N)

    if candidate in current:
        continue

    proposal = list(current)
    proposal[pos] = candidate

    score = motif_score(proposal)

    # Lexicographic objective:
    # internal triples matter much more.
    scalar_old = (
        100000 * current_score[0]
        + current_score[1]
    )
    scalar_new = (
        100000 * score[0]
        + score[1]
    )

    delta = scalar_new - scalar_old

    accept = delta >= 0

    if not accept:
        if random.random() < np.exp(
            delta / max(temperature, 0.001)
        ):
            accept = True

    if accept:
        current = proposal
        current_score = score

    if (
        score[0] > best_score[0]
        or (
            score[0] == best_score[0]
            and score[1] > best_score[1]
        )
    ):
        best_motif = list(proposal)
        best_score = score

        print(
            f"MOTIF|new_best"
            f"|step={step}"
            f"|internal={score[0]}"
            f"|almost={score[1]}"
            f"|vertices={best_motif}",
            flush=True,
        )

    temperature *= 0.99995

    if temperature < 0.05:
        temperature = 5.0

print()
print("BEST MOTIF")
print("size =", K)
print("internal =", best_score[0])
print("two-known =", best_score[1])
print("indices =", best_motif)

GM = P[np.ix_(best_motif, best_motif)]

print()
print("motif Gram:")
print(GM)

# List internal equations.
motif_set = set(best_motif)

internal_triples = [
    tuple(map(int,t))
    for t in triples
    if all(int(v) in motif_set for v in t)
]

print()
print("internal additive equations:")
for t in internal_triples:
    print(t)

# Save.
(OUT / "best_plus2_clique_indices.txt").write_text(
    " ".join(map(str,best_clique)) + "\n"
)

np.savetxt(
    OUT / "best_plus2_clique_vectors.txt",
    signed[best_clique],
    fmt="%d",
)

np.savetxt(
    OUT / "best_plus2_clique_gram.txt",
    GC,
    fmt="%d",
)

(OUT / f"best_relation_motif_{K}_indices.txt").write_text(
    " ".join(map(str,best_motif)) + "\n"
)

np.savetxt(
    OUT / f"best_relation_motif_{K}_vectors.txt",
    signed[best_motif],
    fmt="%d",
)

np.savetxt(
    OUT / f"best_relation_motif_{K}_gram.txt",
    GM,
    fmt="%d",
)

(OUT / f"best_relation_motif_{K}_triples.txt").write_text(
    "\n".join(
        " ".join(map(str,t))
        for t in internal_triples
    ) + "\n"
)

print()
print("done")
