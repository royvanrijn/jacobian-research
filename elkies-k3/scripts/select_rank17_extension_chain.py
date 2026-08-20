#!/usr/bin/env python3
"""Choose eight minimal sections extending the Coxeter-9 scaffold to rank 17.

The numerical reconstruction gives explicit coordinates for the nine oriented
Coxeter-clique sections.  The recovered rank-17 Mordell--Weil lattice tells us
which additional minimal vectors must exist, but there are many choices.

This script chooses an extension chain combinatorially.  At each step it:

* requires the new minimal vector to increase Q-rank;
* rewards |height pairing|=2 links to already selected generators;
* rewards additive-triple frontier incidence with the current closed shell;
* evaluates exact additive closure for the best candidates;
* at the final step prefers a determinant-1 17-vector coordinate basis when
  one is available.

The resulting chain is intended for staged numerical continuation
rank 9 -> 10 -> ... -> 17.
"""

from __future__ import annotations

from pathlib import Path
from collections import deque
import argparse
import csv

import numpy as np


BASE = Path(__file__).resolve().parents[1]
REL = BASE / "data" / "relations"
LAT = BASE / "data" / "lattice"


def rank_mod_rows(rows: np.ndarray, p: int) -> int:
    a = [[int(x) % p for x in row] for row in np.asarray(rows)]
    if not a:
        return 0
    nr = len(a)
    nc = len(a[0])
    r = 0
    for c in range(nc):
        pivot = next((i for i in range(r, nr) if a[i][c] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c] % p, -1, p)
        a[r] = [(v * inv) % p for v in a[r]]
        for i in range(nr):
            if i == r:
                continue
            f = a[i][c] % p
            if f:
                a[i] = [(u - f * v) % p for u, v in zip(a[i], a[r])]
        r += 1
        if r == nr:
            break
    return r


def independent(selected_vectors: list[np.ndarray], v: np.ndarray) -> bool:
    rows = np.vstack(selected_vectors + [v])
    target = len(selected_vectors) + 1
    for p in (1000003, 1000033, 1000037):
        if rank_mod_rows(rows, p) == target:
            return True
    return False


def det_bareiss(matrix: np.ndarray) -> int:
    a = [[int(x) for x in row] for row in np.asarray(matrix)]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = a[i][j] * pivot_value - a[i][k] * a[k][j]
                if k:
                    num //= prev
                a[i][j] = num
        prev = pivot_value
        for i in range(k + 1, n):
            a[i][k] = 0
        for j in range(k + 1, n):
            a[k][j] = 0
    return sign * a[n - 1][n - 1]


parser = argparse.ArgumentParser()
parser.add_argument("--top-candidates", type=int, default=64)
parser.add_argument(
    "--out",
    type=Path,
    default=BASE / "results" / "rank17-extension-chain-v1",
)
args = parser.parse_args()

signed = np.load(REL / "all_2622_signed_short_basis.npy").astype(np.int64, copy=False)
triples = np.load(REL / "minimal_additive_triples.npy").astype(np.int32, copy=False)
H = np.loadtxt(LAT / "short_vector_basis_gram.txt", dtype=np.int64)
clique = np.loadtxt(REL / "best_plus2_clique_indices.txt", dtype=np.int64).reshape(-1)

if len(clique) != 9:
    raise SystemExit(f"expected 9 clique indices, got {len(clique)}")

N = len(signed)
M = len(triples)
lookup = {tuple(map(int, row)): i for i, row in enumerate(signed)}
neg = np.empty(N, dtype=np.int32)
for i, row in enumerate(signed):
    neg[i] = lookup[tuple(map(int, -row))]
if not np.all(neg[neg] == np.arange(N)):
    raise RuntimeError("negation lookup is not involutive")

# One representative per +/- pair.
reps = [i for i in range(N) if i < int(neg[i])]

# Incidence index for additive closure.
flat = triples.reshape(-1)
tids = np.repeat(np.arange(M, dtype=np.int32), 3)
order = np.argsort(flat, kind="stable")
sv = flat[order]
si = tids[order]
offsets = np.searchsorted(sv, np.arange(N + 1, dtype=np.int32))


def incident(v: int) -> np.ndarray:
    return si[offsets[v] : offsets[v + 1]]


def additive_closure(seed) -> np.ndarray:
    known = np.zeros(N, dtype=bool)
    known[list(seed)] = True
    count = np.count_nonzero(known[triples], axis=1).astype(np.uint8)
    ready = np.flatnonzero(count >= 2)
    queued = np.zeros(M, dtype=bool)
    queued[ready] = True
    q = deque(map(int, ready))

    while q:
        tid = q.popleft()
        for vv in triples[tid]:
            v = int(vv)
            if known[v]:
                continue
            known[v] = True
            ids = incident(v)
            old = count[ids].copy()
            count[ids] = old + 1
            nr = ids[(old < 2) & (count[ids] >= 2) & (~queued[ids])]
            queued[nr] = True
            q.extend(map(int, nr))
    return known


selected_indices = list(map(int, clique))
selected_vectors = [signed[i].copy() for i in selected_indices]
if rank_mod_rows(np.vstack(selected_vectors), 1000003) != 9:
    raise RuntimeError("Coxeter clique does not have rank 9")

seed_pm = selected_indices + [int(neg[i]) for i in selected_indices]
known = additive_closure(seed_pm)
steps = []

print("START")
print("selected_rank = 9")
print("additive_closure =", int(np.count_nonzero(known)), "of", N)
print()

for step in range(1, 9):
    rank_now = len(selected_vectors)
    known_count = np.count_nonzero(known[triples], axis=1).astype(np.uint8)
    B = np.vstack(selected_vectors)

    scored = []
    for rep in reps:
        ni = int(neg[rep])
        if known[rep] or known[ni]:
            continue
        v = signed[rep]
        if not independent(selected_vectors, v):
            continue

        pairings = (v @ H) @ B.T
        direct_abs2 = int(np.count_nonzero(np.abs(pairings) == 2))

        ids = np.unique(np.concatenate([incident(rep), incident(ni)]))
        frontier = int(np.count_nonzero(known_count[ids] == 1))
        relation_degree = int(len(incident(rep)) + len(incident(ni)))

        det_abs = None
        if rank_now == 16:
            det_abs = abs(det_bareiss(np.vstack(selected_vectors + [v])))

        scored.append((rep, direct_abs2, frontier, relation_degree, det_abs))

    if not scored:
        raise RuntimeError(f"no independent extension candidate at step {step}")

    # On the final step, determinant one is a major structural win.  If such
    # candidates exist, only evaluate those; otherwise prefer the smallest
    # available determinant before closure scoring.
    if rank_now == 16:
        det1 = [s for s in scored if s[4] == 1]
        if det1:
            scored = det1
            print(f"STEP {step}: determinant-1 candidates = {len(scored)}", flush=True)
        else:
            min_det = min(int(s[4]) for s in scored if s[4] is not None and s[4] > 0)
            scored = [s for s in scored if s[4] == min_det]
            print(f"STEP {step}: no determinant-1 candidate; minimum determinant = {min_det}", flush=True)

    scored.sort(key=lambda s: (s[2], 8 * s[1], s[3]), reverse=True)
    finalists = scored[: max(1, args.top_candidates)]

    evaluated = []
    current_indices = np.flatnonzero(known).tolist()
    for rep, direct_abs2, frontier, relation_degree, det_abs in finalists:
        ni = int(neg[rep])
        c = additive_closure(current_indices + [rep, ni])
        closure_size = int(np.count_nonzero(c))
        evaluated.append((closure_size, direct_abs2, frontier, relation_degree, rep, det_abs, c))

    evaluated.sort(key=lambda z: (z[0], z[1], z[2], z[3]), reverse=True)
    closure_size, direct_abs2, frontier, relation_degree, rep, det_abs, new_known = evaluated[0]
    ni = int(neg[rep])

    # Choose the orientation that turns as many direct links as possible into
    # pairing -2 (sum-minimal) edges.  This is only a convenience for the
    # continuation equations; the x-coordinate class is unchanged.
    pv = (signed[rep] @ H) @ B.T
    minus_for_rep = int(np.count_nonzero(pv == -2))
    minus_for_neg = int(np.count_nonzero(pv == 2))
    chosen = rep if minus_for_rep >= minus_for_neg else ni

    selected_indices.append(int(chosen))
    selected_vectors.append(signed[chosen].copy())
    known = new_known

    step_row = {
        "step": step,
        "rank": len(selected_vectors),
        "signed_index": int(chosen),
        "negative_index": int(neg[chosen]),
        "direct_abs2": direct_abs2,
        "frontier": frontier,
        "relation_degree": relation_degree,
        "closure_size": closure_size,
        "det_abs_if_full": "" if det_abs is None else int(det_abs),
    }
    steps.append(step_row)

    print(
        f"EXT|step={step}"
        f"|rank={len(selected_vectors)}"
        f"|signed={chosen}"
        f"|neg={int(neg[chosen])}"
        f"|abs2={direct_abs2}"
        f"|frontier={frontier}"
        f"|closure={closure_size}"
        + (f"|det={det_abs}" if det_abs is not None else ""),
        flush=True,
    )

B = np.vstack(selected_vectors)
rank17 = rank_mod_rows(B, 1000003)
det_abs = abs(det_bareiss(B)) if B.shape == (17, 17) else 0
Gsel = (B @ H) @ B.T
final_closure = int(np.count_nonzero(known))

out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)
np.savetxt(out / "selected-vectors.txt", B, fmt="%d")
np.savetxt(out / "selected-gram.txt", Gsel, fmt="%d")
(out / "selected-signed-indices.txt").write_text(" ".join(map(str, selected_indices)) + "\n")

with (out / "extension-steps.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(steps[0].keys()))
    writer.writeheader()
    writer.writerows(steps)

# All |pairing|=2 links among the selected oriented generators.
with (out / "abs2-edges.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["a", "b", "signed_a", "signed_b", "pairing", "minimal_operation"])
    for a in range(17):
        for b in range(a + 1, 17):
            p = int(Gsel[a, b])
            if abs(p) != 2:
                continue
            writer.writerow([
                a,
                b,
                selected_indices[a],
                selected_indices[b],
                p,
                "sum" if p == -2 else "difference",
            ])

summary = "\n".join([
    "Coxeter-9 -> rank-17 minimal-section extension chain",
    "====================================================",
    f"initial_clique={' '.join(map(str, clique))}",
    f"final_indices={' '.join(map(str, selected_indices))}",
    f"final_rank={rank17}",
    f"coordinate_determinant_abs={det_abs}",
    f"final_additive_closure={final_closure}/{N}",
    f"selected_abs2_edges={sum(1 for a in range(17) for b in range(a+1,17) if abs(int(Gsel[a,b])) == 2)}",
    "",
    "Interpretation:",
    "  rank 9 is the Coxeter scaffold reconstructed numerically;",
    "  each following row is one independent height-4 section to impose;",
    "  determinant 1 means the final 17 vectors are a unimodular MW basis;",
    "  closure 2622 means their minimal-vector additions generate the full shell.",
]) + "\n"
(out / "summary.txt").write_text(summary)

print()
print(summary, end="")
print("saved =", out)
