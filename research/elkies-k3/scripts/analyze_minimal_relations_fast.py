from __future__ import annotations

from collections import Counter, deque
from pathlib import Path
from time import perf_counter
import argparse

import numpy as np
from sage.all import ZZ, matrix


def read_int_matrix(path: Path) -> np.ndarray:
    rows: list[list[int]] = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        line = line.replace("[", " ").replace("]", " ")
        row = [int(x) for x in line.split()]

        if row:
            rows.append(row)

    if not rows:
        raise RuntimeError(f"No matrix rows found in {path}")

    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise RuntimeError(f"Inconsistent row lengths in {path}")

    return np.asarray(rows, dtype=np.int64)


parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--seed-size", type=int, default=10)
args = parser.parse_args()

checkpoint = Path(args.checkpoint)

gram_path = checkpoint / "short_vector_basis_gram.txt"
vectors_path = checkpoint / "all_1311_in_short_basis.txt"

H = read_int_matrix(gram_path)
half = read_int_matrix(vectors_path)

if H.shape != (17, 17):
    raise RuntimeError(f"Expected 17x17 Gram, got {H.shape}")

if half.shape != (1311, 17):
    raise RuntimeError(
        f"Expected 1311x17 minimal vectors, got {half.shape}"
    )

m, dimension = half.shape
t0 = perf_counter()

# ----------------------------------------------------------------------
# All pairings among the 1311 representatives in one compiled operation.
# ----------------------------------------------------------------------

pairings = (half @ H) @ half.T

if not np.all(np.diag(pairings) == 4):
    raise RuntimeError("Minimal vectors do not all have norm 4")

upper = pairings[np.triu_indices(m, k=1)]
values, counts = np.unique(upper, return_counts=True)

print("minimal +/- pairs:", m)
print()
print("PAIRING DISTRIBUTION")

for value, count in zip(values, counts):
    print(int(value), int(count))

print(f"\npairing matrix: {perf_counter() - t0:.2f}s", flush=True)

# ----------------------------------------------------------------------
# Per-vector absolute-pairing profiles.
# ----------------------------------------------------------------------

absolute = np.abs(pairings)

profiles = np.column_stack([
    np.count_nonzero(absolute == 0, axis=1),
    np.count_nonzero(absolute == 1, axis=1),
    np.count_nonzero(absolute == 2, axis=1),
]).astype(np.int32)

profile_counter = Counter(map(tuple, profiles.tolist()))

print()
print("NUMBER OF DISTINCT ABS-PAIRING PROFILES:", len(profile_counter))

for profile, count in profile_counter.most_common(20):
    print("count", count, "profile", tuple(zip((0, 1, 2), profile)))

rare = sorted(
    range(m),
    key=lambda i: (
        profile_counter[tuple(profiles[i])],
        tuple(profiles[i]),
        i,
    ),
)

print()
print("RARE PROFILE VECTORS")

for i in rare[:30]:
    profile = tuple(zip((0, 1, 2), map(int, profiles[i])))
    print(
        "i", i,
        "multiplicity", profile_counter[tuple(profiles[i])],
        "v", tuple(map(int, half[i])),
        "profile", profile,
    )

# ----------------------------------------------------------------------
# Create signed vectors in the same alternating order as before:
#
#   2*i     =  v_i
#   2*i + 1 = -v_i
#
# We do NOT construct a 2622x2622 pairing matrix.
#
# For representatives i<j:
#
#   <v_i,v_j> = -2 gives pairs
#       ( v_i,  v_j)
#       (-v_i, -v_j)
#
#   <v_i,v_j> = +2 gives pairs
#       ( v_i, -v_j)
#       (-v_i,  v_j)
#
# Those are exactly all signed pairs whose sum has norm 4.
# ----------------------------------------------------------------------

signed = np.empty((2 * m, dimension), dtype=np.int64)
signed[0::2] = half
signed[1::2] = -half

minus_i, minus_j = np.nonzero(np.triu(pairings == -2, k=1))
plus_i, plus_j = np.nonzero(np.triu(pairings == 2, k=1))

left = np.concatenate([
    2 * minus_i,
    2 * minus_i + 1,
    2 * plus_i,
    2 * plus_i + 1,
]).astype(np.int32)

right = np.concatenate([
    2 * minus_j,
    2 * minus_j + 1,
    2 * plus_j + 1,
    2 * plus_j,
]).astype(np.int32)

if not np.all(left < right):
    raise RuntimeError("Signed pair ordering invariant failed")

# Map the sum to its signed minimal-vector index.
signed_lookup = {
    row.tobytes(): index
    for index, row in enumerate(signed)
}

sums = signed[left] + signed[right]

result = np.fromiter(
    (signed_lookup.get(row.tobytes(), -1) for row in sums),
    dtype=np.int32,
    count=len(sums),
)

missing = np.flatnonzero(result < 0)

if len(missing):
    i = int(missing[0])
    raise RuntimeError(
        f"Minimal-vector sum missing: "
        f"{signed[left[i]]} + {signed[right[i]]}"
    )

triples = np.column_stack((left, right, result)).astype(
    np.int32,
    copy=False,
)

triple_count = len(triples)
signed_count = len(signed)

print()
print("signed vectors =", signed_count)
print("additive equations =", triple_count)
print(f"relation extraction: {perf_counter() - t0:.2f}s", flush=True)

# ----------------------------------------------------------------------
# Degree and incidence index.
# ----------------------------------------------------------------------

degree = np.bincount(
    triples.reshape(-1),
    minlength=signed_count,
).astype(np.int32)

top = sorted(
    range(signed_count),
    key=lambda i: (-int(degree[i]), i),
)

print()
print("TOP RELATION-RICH VECTORS")

for i in top[:30]:
    print(
        "i", i,
        "degree", int(degree[i]),
        "v", tuple(map(int, signed[i])),
    )

# Compact incidence representation:
#
# incident triple IDs for vertex v are:
#
#   incident_ids[offsets[v] : offsets[v+1]]
#
flat_vertices = triples.reshape(-1)
flat_triple_ids = np.repeat(
    np.arange(triple_count, dtype=np.int32),
    3,
)

order = np.argsort(flat_vertices, kind="stable")
sorted_vertices = flat_vertices[order]
incident_ids = flat_triple_ids[order]

offsets = np.searchsorted(
    sorted_vertices,
    np.arange(signed_count + 1, dtype=np.int32),
)


def incident(vertex: int) -> np.ndarray:
    return incident_ids[offsets[vertex]:offsets[vertex + 1]]


# ----------------------------------------------------------------------
# Incremental greedy seed construction.
#
# The original code rescanned every triple for every candidate.
#
# Here, adding candidate v only affects triples incident to v.
# ----------------------------------------------------------------------

selected = np.zeros(signed_count, dtype=bool)
selected_per_triple = np.zeros(triple_count, dtype=np.uint8)

chosen: list[int] = [top[0]]
selected[top[0]] = True
selected_per_triple[incident(top[0])] += 1

print()
print("GREEDY RELATION SEED")
print(
    "choose 1",
    "index", top[0],
    "degree", int(degree[top[0]]),
    "internal", 0,
    "two-known", 0,
)

while len(chosen) < args.seed_size:
    scored: list[tuple[int, int, int, int]] = []

    for vertex in np.flatnonzero(~selected):
        vertex = int(vertex)
        ids = incident(vertex)
        existing = selected_per_triple[ids]

        # Adding vertex:
        #   count 2 -> count 3: a newly internal relation
        #   count 1 -> count 2: a newly almost-closed relation
        gain_internal = int(np.count_nonzero(existing == 2))
        gain_two = (
            int(np.count_nonzero(existing == 1))
            - gain_internal
        )

        scored.append((
            -gain_internal,
            -gain_two,
            -int(degree[vertex]),
            vertex,
        ))

    scored.sort()

    picked = None

    # Relation score first; exact independence check only for the
    # handful of top-scoring candidates.
    for _, _, _, vertex in scored:
        candidate = matrix(
            ZZ,
            signed[chosen + [vertex]].tolist(),
        )

        if candidate.rank() == len(chosen) + 1:
            picked = vertex
            break

    if picked is None:
        raise RuntimeError("Could not extend independent seed")

    chosen.append(picked)
    selected[picked] = True
    selected_per_triple[incident(picked)] += 1

    internal = int(np.count_nonzero(selected_per_triple == 3))
    two_known = int(np.count_nonzero(selected_per_triple == 2))

    print(
        "choose", len(chosen),
        "index", picked,
        "degree", int(degree[picked]),
        "internal", internal,
        "two-known", two_known,
    )

print(f"greedy seed: {perf_counter() - t0:.2f}s", flush=True)

# ----------------------------------------------------------------------
# Queue-based additive closure.
#
# Whenever any equation a+b=c has two known vertices, the third is
# forced. Each incidence is processed only a small number of times.
# ----------------------------------------------------------------------

closure = np.zeros(signed_count, dtype=bool)
closure[chosen] = True

closure_count = np.count_nonzero(
    closure[triples],
    axis=1,
).astype(np.uint8)

initial = np.flatnonzero(closure_count >= 2)

queued = np.zeros(triple_count, dtype=bool)
queued[initial] = True

queue = deque(map(int, initial))

while queue:
    triple_id = queue.popleft()
    vertices = triples[triple_id]

    missing_vertices = [
        int(v)
        for v in vertices
        if not closure[v]
    ]

    for vertex in missing_vertices:
        closure[vertex] = True

        ids = incident(vertex)
        old = closure_count[ids].copy()
        closure_count[ids] = old + 1

        newly_ready = ids[
            (old < 2)
            & (closure_count[ids] >= 2)
            & (~queued[ids])
        ]

        queued[newly_ready] = True
        queue.extend(map(int, newly_ready))

closure_indices = np.flatnonzero(closure)

print()
print("closure size =", len(closure_indices), "of", signed_count)
print(f"closure: {perf_counter() - t0:.2f}s", flush=True)

# ----------------------------------------------------------------------
# Save binary relation data and human-readable seed data.
# ----------------------------------------------------------------------

np.save(
    checkpoint / "all_2622_signed_short_basis.npy",
    signed,
)

np.save(
    checkpoint / "minimal_additive_triples.npy",
    triples,
)

np.save(
    checkpoint / "minimal_relation_degrees.npy",
    degree,
)

np.save(
    checkpoint / "minimal_profiles.npy",
    profiles,
)

np.savetxt(
    checkpoint / "relation_seed_indices.txt",
    np.asarray(chosen, dtype=np.int32),
    fmt="%d",
)

np.savetxt(
    checkpoint / "relation_seed_vectors.txt",
    signed[chosen],
    fmt="%d",
)

np.savetxt(
    checkpoint / "relation_closure_indices.txt",
    closure_indices,
    fmt="%d",
)

H_sage = matrix(ZZ, H.tolist())
B_sage = matrix(ZZ, signed[chosen].tolist())
seed_gram = B_sage * H_sage * B_sage.transpose()

(checkpoint / "relation_seed_gram.txt").write_text(
    "\n".join(
        " ".join(map(str, row))
        for row in seed_gram.rows()
    ) + "\n"
)

print()
print("CHOSEN INDICES")
print(chosen)

print()
print("CHOSEN GRAM")
print(seed_gram)

print()
print("saved:")
print(checkpoint / "minimal_additive_triples.npy")
print(checkpoint / "minimal_relation_degrees.npy")
print(checkpoint / "relation_seed_vectors.txt")
print(checkpoint / "relation_seed_gram.txt")
print(checkpoint / "relation_closure_indices.txt")

print()
print(f"TOTAL TIME: {perf_counter() - t0:.2f}s")
