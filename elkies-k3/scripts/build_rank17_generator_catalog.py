from pathlib import Path
from collections import Counter
import numpy as np
import csv

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
H_DIR = BASE / "data/holes"
OUT = BASE / "results/rank17-generator"

OUT.mkdir(parents=True, exist_ok=True)

H = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

signed = np.load(
    J / "all_2622_signed_short_basis.npy"
).astype(np.int64)

triples = np.load(
    J / "minimal_additive_triples.npy"
).astype(np.int64)

assert H.shape == (17,17)
assert signed.shape[1] == 17

# ------------------------------------------------------------
# Minimal-vector statistics
# ------------------------------------------------------------

P = (signed @ H) @ signed.T

assert np.all(np.diag(P) == 4)

degree2 = np.sum(P == 2, axis=1) - 1 * (np.diag(P) == 2)
degree1 = np.sum(np.abs(P) == 1, axis=1)

incident = np.zeros(len(signed), dtype=np.int64)

for t in triples:
    incident[t] += 1

# ------------------------------------------------------------
# Read ranked hole witnesses
# ------------------------------------------------------------

holes = []

with open(H_DIR / "candidate_half_holes.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        v = np.array(
            [int(x) for x in row["witness"].split()],
            dtype=np.int64
        )

        assert len(v) == 17

        norm = int(v @ H @ v)

        pair = signed @ H @ v

        holes.append({
            "priority": int(row["priority"]),
            "multiplicity": int(row["multiplicity"]),
            "mask": row["mask"],
            "hex": row["hex"],
            "v": v,
            "norm": norm,
            "pair0": int(np.sum(pair == 0)),
            "pair1": int(np.sum(np.abs(pair) == 1)),
            "pair2": int(np.sum(np.abs(pair) == 2)),
            "pair3plus": int(np.sum(np.abs(pair) >= 3)),
        })

# ------------------------------------------------------------
# Rank generator candidates.
#
# Prefer:
#   * high hole multiplicity
#   * simple pairings to minimal sections
#   * relatively short witness
# ------------------------------------------------------------

for x in holes:
    x["score"] = (
        1000 * x["multiplicity"]
        + 5 * x["pair2"]
        + 2 * x["pair1"]
        + x["pair0"]
        - 10 * x["pair3plus"]
        - x["norm"]
    )

holes.sort(
    key=lambda x: (
        -x["score"],
        x["norm"],
        x["priority"]
    )
)

print("minimal signed vectors =", len(signed))
print("additive triples =", len(triples))
print("hole candidates =", len(holes))
print()

print("TOP GENERATOR DIRECTIONS")

for i,x in enumerate(holes[:100]):
    print(
        f"GEN|rank={i+1}"
        f"|score={x['score']}"
        f"|priority={x['priority']}"
        f"|mult={x['multiplicity']}"
        f"|norm={x['norm']}"
        f"|p0={x['pair0']}"
        f"|p1={x['pair1']}"
        f"|p2={x['pair2']}"
        f"|p3plus={x['pair3plus']}"
        f"|v={' '.join(map(str,x['v']))}"
    )

with open(OUT / "generator_candidates.tsv", "w") as f:
    f.write(
        "rank\tscore\tpriority\tmultiplicity\tnorm\t"
        "pair0\tpair1\tpair2\tpair3plus\tvector\n"
    )

    for i,x in enumerate(holes):
        f.write(
            f"{i+1}\t{x['score']}\t{x['priority']}\t"
            f"{x['multiplicity']}\t{x['norm']}\t"
            f"{x['pair0']}\t{x['pair1']}\t"
            f"{x['pair2']}\t{x['pair3plus']}\t"
            + " ".join(map(str,x["v"]))
            + "\n"
        )

np.savetxt(
    OUT / "top100_vectors.txt",
    np.array([x["v"] for x in holes[:100]]),
    fmt="%d"
)

print()
print("saved under", OUT)
