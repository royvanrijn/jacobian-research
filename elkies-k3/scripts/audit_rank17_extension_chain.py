#!/usr/bin/env python3
"""Audit the Coxeter-9 -> rank-17 extension chain.

The raw 9-vector Coxeter clique is not a primitive basis of its saturated
rank-9 sublattice.  Therefore a determinant-1 17x17 matrix containing those
nine rows may be impossible even when the chosen 17 vectors are optimal.

This script:

* computes the exact saturation index of the raw Coxeter clique from Gram
  determinant ratios;
* checks the determinant of the selected 17-vector chain;
* for an index-2 final lattice, recovers the unique mod-2 parity obstruction;
* finds a minimal-vector representative of the missing saturation coset;
* measures additive closure after adjoining that bridge;
* exports all |height pairing|=2 links usable as continuation constraints;
* selects up to three earlier anchors for each of the eight independent
  extension sections.

No floating-point lattice decisions are used.
"""

from __future__ import annotations

from pathlib import Path
from collections import deque
import argparse
import csv
import math

import numpy as np


BASE = Path(__file__).resolve().parents[1]
REL = BASE / "data" / "relations"
LAT = BASE / "data" / "lattice"


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
        pv = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = a[i][j] * pv - a[i][k] * a[k][j]
                if k:
                    num //= prev
                a[i][j] = num
        prev = pv
        for i in range(k + 1, n):
            a[i][k] = 0
        for j in range(k + 1, n):
            a[k][j] = 0
    return sign * a[-1][-1]


def rank_mod2(rows: np.ndarray):
    a = (np.asarray(rows, dtype=np.int64) & 1).copy()
    nr, nc = a.shape
    pivots = []
    r = 0
    for c in range(nc):
        candidates = np.flatnonzero(a[r:, c])
        if not len(candidates):
            continue
        p = r + int(candidates[0])
        if p != r:
            a[[r, p]] = a[[p, r]]
        for i in range(nr):
            if i != r and a[i, c]:
                a[i] ^= a[r]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    return r, a, pivots


def right_null_vector_mod2(rows: np.ndarray) -> np.ndarray:
    rank, rref, pivots = rank_mod2(rows)
    n = rows.shape[1]
    free = [c for c in range(n) if c not in pivots]
    if len(free) != 1:
        raise RuntimeError(
            f"expected one-dimensional right nullspace mod 2, got {len(free)}"
        )
    f = free[0]
    h = np.zeros(n, dtype=np.int64)
    h[f] = 1
    for row, pivot in enumerate(pivots):
        h[pivot] = rref[row, f]
    if np.any((np.asarray(rows, dtype=np.int64) @ h) & 1):
        raise RuntimeError("computed parity vector is not in the mod-2 nullspace")
    return h


parser = argparse.ArgumentParser()
parser.add_argument(
    "--chain",
    type=Path,
    default=BASE / "results" / "rank17-extension-chain-v1",
)
parser.add_argument("--top-bridge-candidates", type=int, default=96)
parser.add_argument("--anchors-per-extension", type=int, default=3)
parser.add_argument(
    "--out",
    type=Path,
    default=BASE / "results" / "rank17-extension-audit-v1",
)
args = parser.parse_args()

chain = args.chain.resolve()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)

required = [
    chain / "selected-vectors.txt",
    chain / "selected-signed-indices.txt",
    REL / "all_2622_signed_short_basis.npy",
    REL / "minimal_additive_triples.npy",
    REL / "best_plus2_clique_gram.txt",
    REL / "clique9_saturated_gram.txt",
    LAT / "short_vector_basis_gram.txt",
]
for p in required:
    if not p.exists():
        raise SystemExit(f"missing required input: {p}")

B = np.loadtxt(chain / "selected-vectors.txt", dtype=np.int64)
selected_indices = [
    int(x) for x in (chain / "selected-signed-indices.txt").read_text().split()
]
signed = np.load(REL / "all_2622_signed_short_basis.npy").astype(np.int64, copy=False)
triples = np.load(REL / "minimal_additive_triples.npy").astype(np.int32, copy=False)
H = np.loadtxt(LAT / "short_vector_basis_gram.txt", dtype=np.int64)
Graw = np.loadtxt(REL / "best_plus2_clique_gram.txt", dtype=np.int64)
Gsat = np.loadtxt(REL / "clique9_saturated_gram.txt", dtype=np.int64)

if B.shape != (17, 17):
    raise SystemExit(f"expected selected-vectors shape 17x17, got {B.shape}")
if len(selected_indices) != 17:
    raise SystemExit(f"expected 17 selected signed indices, got {len(selected_indices)}")

raw_det = abs(det_bareiss(Graw))
sat_det = abs(det_bareiss(Gsat))
if raw_det % sat_det:
    raise RuntimeError("raw/saturated Coxeter Gram determinants are inconsistent")
ratio = raw_det // sat_det
clique_index = math.isqrt(ratio)
if clique_index * clique_index != ratio:
    raise RuntimeError("Coxeter Gram determinant ratio is not a square")

final_det = abs(det_bareiss(B))
optimal_given_fixed_clique = final_det == clique_index

print("COXETER SATURATION")
print("raw_clique_gram_det =", raw_det)
print("saturated_clique_gram_det =", sat_det)
print("raw_clique_saturation_index =", clique_index)
print()
print("FINAL CHAIN")
print("coordinate_determinant_abs =", final_det)
print("minimum_possible_with_fixed_raw_clique =", clique_index)
print("optimal_given_fixed_clique =", optimal_given_fixed_clique)

# ---------------------------------------------------------------------------
# Signed-vector negation and additive closure.
# ---------------------------------------------------------------------------

N = len(signed)
lookup = {tuple(map(int, row)): i for i, row in enumerate(signed)}
neg = np.empty(N, dtype=np.int32)
for i, row in enumerate(signed):
    neg[i] = lookup[tuple(map(int, -row))]

flat = triples.reshape(-1)
tids = np.repeat(np.arange(len(triples), dtype=np.int32), 3)
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
    queued = np.zeros(len(triples), dtype=bool)
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

seed_pm = selected_indices + [int(neg[i]) for i in selected_indices]
base_closure = additive_closure(seed_pm)
base_closure_size = int(np.count_nonzero(base_closure))
print("additive_closure_without_bridge =", f"{base_closure_size}/{N}")

# ---------------------------------------------------------------------------
# If final index is 2, membership in the selected row lattice is equivalent
# to one parity equation h.v = 0 mod 2.  This lets us find the missing coset
# exactly without rational matrix inversion.
# ---------------------------------------------------------------------------

bridge_index = None
bridge_neg = None
bridge_closure_size = base_closure_size
parity_vector = None
bridge_pairings = None

if final_det == 2:
    mod2_rank, _, _ = rank_mod2(B)
    if mod2_rank != 16:
        raise RuntimeError(f"determinant 2 but mod-2 rank is {mod2_rank}, expected 16")
    parity_vector = right_null_vector_mod2(B)
    print("index2_parity_vector =", " ".join(map(str, parity_vector)))

    reps = [i for i in range(N) if i < int(neg[i])]
    candidates = []
    for rep in reps:
        v = signed[rep]
        if int(np.dot(v, parity_vector)) % 2 == 0:
            continue
        pairings = (v @ H) @ B.T
        abs2 = int(np.count_nonzero(np.abs(pairings) == 2))
        relation_degree = int(len(incident(rep)) + len(incident(int(neg[rep]))))
        candidates.append((abs2, relation_degree, rep, pairings))

    candidates.sort(key=lambda z: (z[0], z[1]), reverse=True)
    finalists = candidates[: max(1, args.top_bridge_candidates)]

    evaluated = []
    current = np.flatnonzero(base_closure).tolist()
    for abs2, relation_degree, rep, pairings in finalists:
        c = additive_closure(current + [rep, int(neg[rep])])
        evaluated.append(
            (
                int(np.count_nonzero(c)),
                abs2,
                relation_degree,
                rep,
                pairings,
            )
        )
    evaluated.sort(key=lambda z: (z[0], z[1], z[2]), reverse=True)
    bridge_closure_size, bridge_abs2, bridge_degree, bridge_index, bridge_pairings = evaluated[0]
    bridge_neg = int(neg[bridge_index])

    # Pick the orientation with more -2 (sum-minimal) links for convenience.
    plus = signed[bridge_index]
    pp = (plus @ H) @ B.T
    if np.count_nonzero(pp == 2) > np.count_nonzero(pp == -2):
        bridge_index, bridge_neg = bridge_neg, bridge_index
        bridge_pairings = -bridge_pairings

    print("saturation_bridge_signed_index =", bridge_index)
    print("saturation_bridge_negative_index =", bridge_neg)
    print("saturation_bridge_abs2_links =", bridge_abs2)
    print("additive_closure_with_bridge =", f"{bridge_closure_size}/{N}")

    np.savetxt(out / "saturation-bridge-vector.txt", signed[bridge_index].reshape(1, -1), fmt="%d")
    (out / "saturation-bridge-index.txt").write_text(f"{bridge_index}\n")
    np.savetxt(out / "index2-parity-vector.txt", parity_vector.reshape(1, -1), fmt="%d")

# ---------------------------------------------------------------------------
# Export all abs-2 links and a compact continuation anchor plan.
# ---------------------------------------------------------------------------

Gsel = (B @ H) @ B.T
all_edges = []
for qpos in range(9, 17):
    for apos in range(qpos):
        p = int(Gsel[qpos, apos])
        if abs(p) != 2:
            continue
        all_edges.append(
            {
                "stage_rank": qpos + 1,
                "section_position": qpos,
                "section_signed_index": selected_indices[qpos],
                "anchor_position": apos,
                "anchor_signed_index": selected_indices[apos],
                "pairing": p,
                "operation": "sum" if p == -2 else "difference",
                "anchor_is_coxeter": int(apos < 9),
            }
        )

with (out / "all-extension-abs2-anchors.tsv").open("w", newline="") as handle:
    fields = list(all_edges[0].keys())
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(all_edges)

plan = []
for qpos in range(9, 17):
    edges = [e for e in all_edges if int(e["section_position"]) == qpos]
    # Prefer anchors from the original Coxeter scaffold: those coordinates are
    # reconstructed directly from the slope system and tend to be best conditioned.
    edges.sort(
        key=lambda e: (
            -int(e["anchor_is_coxeter"]),
            int(e["anchor_position"]),
        )
    )
    chosen = edges[: args.anchors_per_extension]
    if not chosen:
        raise RuntimeError(f"extension position {qpos} has no abs-2 earlier anchor")
    for n, e in enumerate(chosen):
        row = dict(e)
        row["role"] = "primary" if n == 0 else "secondary"
        plan.append(row)

with (out / "continuation-anchors.tsv").open("w", newline="") as handle:
    fields = [
        "stage_rank",
        "section_position",
        "section_signed_index",
        "anchor_position",
        "anchor_signed_index",
        "pairing",
        "operation",
        "anchor_is_coxeter",
        "role",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(plan)

summary = "\n".join(
    [
        "Rank-17 extension-chain saturation audit",
        "=========================================",
        f"raw_clique_gram_det={raw_det}",
        f"saturated_clique_gram_det={sat_det}",
        f"raw_clique_saturation_index={clique_index}",
        f"final_coordinate_determinant_abs={final_det}",
        f"minimum_possible_with_fixed_raw_clique={clique_index}",
        f"optimal_given_fixed_clique={optimal_given_fixed_clique}",
        f"base_additive_closure={base_closure_size}/{N}",
        f"bridge_signed_index={bridge_index if bridge_index is not None else ''}",
        f"bridge_negative_index={bridge_neg if bridge_neg is not None else ''}",
        f"bridge_additive_closure={bridge_closure_size}/{N}",
        f"extension_abs2_edges={len(all_edges)}",
        f"continuation_anchor_constraints={len(plan)}",
        "",
        "Interpretation:",
        "  determinant 2 is optimal when the fixed raw Coxeter clique itself",
        "  has saturation index 2; determinant 1 is then impossible without",
        "  replacing at least one of those nine raw clique generators.",
        "  The bridge is a norm-4 vector in the missing index-2 coset and is",
        "  useful as a final gluing/certification constraint after rank 17.",
    ]
) + "\n"
(out / "summary.txt").write_text(summary)

print()
print(summary, end="")
print("ANCHOR PLAN")
for row in plan:
    print(
        f"ANCHOR|rank={row['stage_rank']}"
        f"|qpos={row['section_position']}"
        f"|q={row['section_signed_index']}"
        f"|apos={row['anchor_position']}"
        f"|a={row['anchor_signed_index']}"
        f"|pairing={int(row['pairing']):+d}"
        f"|role={row['role']}"
    )
print("saved =", out)
