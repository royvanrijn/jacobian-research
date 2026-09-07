from pathlib import Path
import argparse
import csv
import numpy as np
from sage.all import Matrix, RDF, ZZ

parser = argparse.ArgumentParser()
parser.add_argument("gram")
parser.add_argument("--templates", default=None)
parser.add_argument("--limit", type=int, default=5000)
parser.add_argument("--norm-factor", type=float, default=2.5)
args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

templates_path = (
    Path(args.templates)
    if args.templates
    else BASE / "results/rank17-generator/search_templates.tsv"
)

H17 = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=float
)

HA_np = np.loadtxt(args.gram, dtype=float)
r = HA_np.shape[0]

assert HA_np.shape == (r, r)

# ------------------------------------------------------------
# Load templates
# ------------------------------------------------------------

templates = []

with templates_path.open() as f:
    rd = csv.DictReader(f, delimiter="\t")

    for row in rd:
        templates.append({
            "id": int(row["id"]),
            "norm": float(row["norm"]),
            "coords": np.array(
                [int(x) for x in row["coords"].split(",")],
                dtype=np.int64
            ),
            "pairings": np.array(
                [int(x) for x in row["pairings"].split(",")],
                dtype=np.int64
            ),
        })

print(f"ambient_rank={r}")
print(f"templates={len(templates)}")

# ------------------------------------------------------------
# LLL-reduce ambient Gram
# ------------------------------------------------------------

HA = Matrix(RDF, HA_np.tolist())

print("LLL start", flush=True)
U = HA.LLL_gram()
print("LLL done", flush=True)

R = U.transpose() * HA * U
R_np = np.array(R, dtype=float)

# Ambient reduced basis vectors as coordinates in original basis.
U_np = np.array(U, dtype=np.int64)

diag = np.diag(R_np)

print(
    f"ambient_reduced_min={diag.min():.12g}"
    f"|median={np.median(diag):.12g}"
    f"|max={diag.max():.12g}"
)

# ------------------------------------------------------------
# Generate short integer combinations of reduced basis vectors.
#
# Start with:
#   ±e_i
#   ±e_i ± e_j
#   ±2e_i
#
# These are cheap and often enough to expose candidate shells.
# ------------------------------------------------------------

cand = {}

def add_candidate(z):
    z = np.asarray(z, dtype=np.int64)

    if not np.any(z):
        return

    # Canonical sign.
    nz = np.flatnonzero(z)
    if z[nz[0]] < 0:
        z = -z

    key = tuple(map(int, z))

    if key in cand:
        return

    norm = float(z @ R_np @ z)

    cand[key] = norm


for i in range(r):
    z = np.zeros(r, dtype=np.int64)
    z[i] = 1
    add_candidate(z)

    z2 = np.zeros(r, dtype=np.int64)
    z2[i] = 2
    add_candidate(z2)


for i in range(r):
    for j in range(i + 1, r):
        for a in (-1, 1):
            for b in (-1, 1):
                z = np.zeros(r, dtype=np.int64)
                z[i] = a
                z[j] = b
                add_candidate(z)


items = sorted(
    cand.items(),
    key=lambda kv: kv[1]
)[:args.limit]

print(f"short_candidates={len(items)}")

Z = np.array(
    [np.array(k, dtype=np.int64) for k, _ in items]
)

N = np.array(
    [n for _, n in items]
)

# Convert reduced coordinates back into original ambient basis.
V = Z @ U_np.T

# ------------------------------------------------------------
# Compare shells to rank17 template norms.
#
# Since specialization introduces an unknown scale h, each
# candidate/template pair implies:
#
#     h = ambient_norm / template_norm
#
# We cluster these implied scales.
# ------------------------------------------------------------

hits = []

for ti, t in enumerate(templates):

    tn = t["norm"]

    for ci, an in enumerate(N):

        h = an / tn

        if h <= 0:
            continue

        hits.append((
            h,
            ti,
            ci,
            an
        ))

hits.sort(key=lambda x: x[0])

# ------------------------------------------------------------
# Cluster nearby implied h-values.
# ------------------------------------------------------------

clusters = []

eps = 0.01  # 1% multiplicative window

for hit in hits:

    h = hit[0]

    if not clusters:
        clusters.append([hit])
        continue

    center = np.median(
        [x[0] for x in clusters[-1]]
    )

    if abs(h / center - 1.0) <= eps:
        clusters[-1].append(hit)
    else:
        clusters.append([hit])

clusters.sort(
    key=len,
    reverse=True
)

print()
print("TOP SCALE CLUSTERS")

for rank, cl in enumerate(clusters[:30], start=1):

    hs = [x[0] for x in cl]
    center = float(np.median(hs))

    unique_templates = len(set(x[1] for x in cl))
    unique_candidates = len(set(x[2] for x in cl))

    print(
        f"SCALE|rank={rank}"
        f"|h={center:.12g}"
        f"|hits={len(cl)}"
        f"|templates={unique_templates}"
        f"|candidates={unique_candidates}"
    )

    # Print a few representative candidate vectors.
    for h, ti, ci, an in cl[:8]:

        t = templates[ti]

        print(
            f" HIT"
            f"|template={t['id']}"
            f"|tnorm={t['norm']}"
            f"|anorm={an:.12g}"
            f"|h={h:.12g}"
            f"|ambient={V[ci].tolist()}"
        )

print()
print("DONE")
