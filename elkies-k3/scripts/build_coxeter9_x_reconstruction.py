#!/usr/bin/env python3
"""Build the reduced x-coordinate reconstruction problem for the rank-17 K3.

The recovered Mordell--Weil lattice contains a 9-element clique of minimal
vectors with pairwise height pairing +2.  If V_i,V_j are two clique sections,
then D_ij=V_i-V_j is again a minimal section.  Up to sign (which does not
change x-coordinates), the 9 V_i and 36 D_ij give 45 height-4 sections.

For every minimal additive triple P+Q=R on a short Weierstrass model

    y^2 = x^3 + A(T)x + B(T),

with deg(x)<=4, the chord slope has degree <=2 and

    x_P(T) + x_Q(T) + x_R(T) = m_PQR(T)^2.

The Coxeter-9 subsystem supplies

  * 36 relations V_i + (-V_j) = D_ij;
  * 84 relations D_ij + D_jk = D_ik.

Hence there are 120 square relations on 45 quartic x-polynomials.  For fixed
quadratic slopes, all x-coefficients occur linearly with the same 120x45
incidence matrix C.  Eliminating x says that, coefficient by coefficient, the
120-vector of slope-square coefficients must lie in col(C).

This script constructs that subsystem, verifies every oriented relation against
the 184242 saved additive triples, computes exact modular rank certificates for
C, and saves a numerically orthonormal basis of the left nullspace for the next
finite-field/numerical reconstruction stage.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv

import numpy as np


BASE = Path(__file__).resolve().parents[1]
DEFAULT_REL = BASE / "data" / "relations"


def canonical_unsigned(v: np.ndarray) -> tuple[int, ...]:
    """Canonical key for the x-coordinate class {v,-v}."""
    a = np.asarray(v, dtype=np.int64)
    nz = np.flatnonzero(a)
    if len(nz) == 0:
        raise ValueError("zero vector is not a section node")
    if int(a[int(nz[0])]) < 0:
        a = -a
    return tuple(map(int, a))


def rank_mod(matrix: np.ndarray, p: int) -> int:
    """Exact Gaussian-elimination rank over F_p."""
    a = np.asarray(matrix, dtype=object).copy()
    a = np.vectorize(lambda x: int(x) % p, otypes=[object])(a)
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        pivot = None
        for i in range(r, rows):
            if int(a[i, c]) % p:
                pivot = i
                break
        if pivot is None:
            continue
        if pivot != r:
            a[[r, pivot]] = a[[pivot, r]]
        inv = pow(int(a[r, c]) % p, -1, p)
        for j in range(c, cols):
            a[r, j] = (int(a[r, j]) * inv) % p
        for i in range(rows):
            if i == r:
                continue
            factor = int(a[i, c]) % p
            if not factor:
                continue
            for j in range(c, cols):
                a[i, j] = (int(a[i, j]) - factor * int(a[r, j])) % p
        r += 1
        if r == rows:
            break
    return r


def vector_text(v: tuple[int, ...]) -> str:
    return " ".join(map(str, v))


parser = argparse.ArgumentParser()
parser.add_argument(
    "--relations-dir",
    type=Path,
    default=DEFAULT_REL,
    help="directory containing signed vectors, additive triples and clique indices",
)
parser.add_argument(
    "--out",
    type=Path,
    default=BASE / "results" / "coxeter9-x-reconstruction-v1",
)
args = parser.parse_args()

rel = args.relations_dir.resolve()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)

signed_path = rel / "all_2622_signed_short_basis.npy"
triples_path = rel / "minimal_additive_triples.npy"
clique_path = rel / "best_plus2_clique_indices.txt"

for path in (signed_path, triples_path, clique_path):
    if not path.exists():
        raise SystemExit(f"missing required input: {path}")

signed = np.load(signed_path).astype(np.int64, copy=False)
triples = np.load(triples_path).astype(np.int32, copy=False)
clique_indices = np.loadtxt(clique_path, dtype=np.int64).reshape(-1)

if len(clique_indices) != 9:
    raise SystemExit(f"expected a 9-vector clique, got {len(clique_indices)}")
if signed.ndim != 2 or signed.shape[1] != 17:
    raise SystemExit(f"unexpected signed-vector shape {signed.shape}")
if triples.ndim != 2 or triples.shape[1] != 3:
    raise SystemExit(f"unexpected triple shape {triples.shape}")

signed_lookup = {
    tuple(map(int, row)): i
    for i, row in enumerate(signed)
}

# Saved triples have left < right; the result index remains oriented.
triple_lookup = {
    (int(a), int(b), int(c))
    for a, b, c in triples
}


def signed_index(v: np.ndarray) -> int:
    key = tuple(map(int, v))
    try:
        return int(signed_lookup[key])
    except KeyError as exc:
        raise RuntimeError(f"minimal vector missing from signed catalog: {key}") from exc


def verify_addition(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[int, int, int]:
    """Verify a+b=c but preserve participant orientation in the return value."""
    if not np.array_equal(a + b, c):
        raise RuntimeError("internal vector-addition identity failed")
    ia = signed_index(a)
    ib = signed_index(b)
    ic = signed_index(c)
    lookup_a, lookup_b = (ia, ib) if ia < ib else (ib, ia)
    if (lookup_a, lookup_b, ic) not in triple_lookup:
        raise RuntimeError(
            "relation absent from global additive catalog: "
            f"signed=({ia},{ib},{ic}), lookup=({lookup_a},{lookup_b},{ic})"
        )
    # Do not sort ia,ib here: the next reconstruction stage needs signed_a to
    # correspond to node_a (and signed_b to node_b) in order to recover y signs.
    return ia, ib, ic


V = [signed[int(i)].copy() for i in clique_indices]

gram_path = BASE / "data" / "lattice" / "short_vector_basis_gram.txt"
if gram_path.exists():
    H = np.loadtxt(gram_path, dtype=np.int64)
    CV = np.vstack(V)
    G = (CV @ H) @ CV.T
    expected = np.full((9, 9), 2, dtype=np.int64)
    np.fill_diagonal(expected, 4)
    if not np.array_equal(G, expected):
        raise RuntimeError(f"saved clique does not have Coxeter +2 Gram:\n{G}")

# ---------------------------------------------------------------------------
# 45 unsigned x-coordinate nodes.
# ---------------------------------------------------------------------------

nodes: list[dict[str, object]] = []
node_by_key: dict[tuple[int, ...], int] = {}


def add_node(label: str, vector: np.ndarray, kind: str) -> int:
    key = canonical_unsigned(vector)
    if key in node_by_key:
        return node_by_key[key]
    node_id = len(nodes)
    node_by_key[key] = node_id
    representative = np.asarray(key, dtype=np.int64)
    nodes.append(
        {
            "node_id": node_id,
            "label": label,
            "kind": kind,
            "vector": key,
            "signed_index": signed_index(representative),
        }
    )
    return node_id


vertex_node: dict[int, int] = {}
for i in range(9):
    vertex_node[i] = add_node(f"V{i}", V[i], "vertex")

diff_vector: dict[tuple[int, int], np.ndarray] = {}
diff_node: dict[tuple[int, int], int] = {}
for i in range(9):
    for j in range(i + 1, 9):
        d = V[i] - V[j]
        diff_vector[i, j] = d
        diff_node[i, j] = add_node(f"D{i}_{j}", d, "difference")

if len(nodes) != 45:
    raise RuntimeError(
        f"expected 45 distinct unsigned x-section classes, got {len(nodes)}"
    )

# ---------------------------------------------------------------------------
# 120 square relations.
# ---------------------------------------------------------------------------

relations: list[dict[str, object]] = []
relation_node_sets: set[tuple[int, int, int]] = set()


def add_relation(
    kind: str,
    node_ids: tuple[int, int, int],
    signed_identity: tuple[np.ndarray, np.ndarray, np.ndarray],
    i: int,
    j: int,
    k: int | None = None,
) -> None:
    key = tuple(sorted(map(int, node_ids)))
    if key in relation_node_sets:
        raise RuntimeError(f"duplicate unsigned square relation {key}")
    relation_node_sets.add(key)

    ia, ib, ic = verify_addition(*signed_identity)
    relations.append(
        {
            "relation_id": len(relations),
            "kind": kind,
            "i": i,
            "j": j,
            "k": "" if k is None else k,
            "node_a": int(node_ids[0]),
            "node_b": int(node_ids[1]),
            "node_c": int(node_ids[2]),
            "signed_a": ia,
            "signed_b": ib,
            "signed_c": ic,
        }
    )


# V_i + (-V_j) = V_i - V_j.
for i in range(9):
    for j in range(i + 1, 9):
        d = diff_vector[i, j]
        add_relation(
            "vertex_difference",
            (vertex_node[i], vertex_node[j], diff_node[i, j]),
            (V[i], -V[j], d),
            i,
            j,
        )

# D_ij + D_jk = D_ik.
for i in range(9):
    for j in range(i + 1, 9):
        for k in range(j + 1, 9):
            dij = diff_vector[i, j]
            djk = diff_vector[j, k]
            dik = diff_vector[i, k]
            add_relation(
                "difference_triangle",
                (diff_node[i, j], diff_node[j, k], diff_node[i, k]),
                (dij, djk, dik),
                i,
                j,
                k,
            )

if len(relations) != 120:
    raise RuntimeError(f"expected 120 square relations, got {len(relations)}")

# C[r,n]=1 iff x_n appears in relation r.
C = np.zeros((len(relations), len(nodes)), dtype=np.int64)
for row in relations:
    r = int(row["relation_id"])
    for field in ("node_a", "node_b", "node_c"):
        C[r, int(row[field])] += 1

if not np.all(C.sum(axis=1) == 3):
    raise RuntimeError("incidence row does not contain exactly three nodes")

primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
mod_ranks = {p: rank_mod(C, p) for p in primes}
max_mod_rank = max(mod_ranks.values())
numeric_rank = int(np.linalg.matrix_rank(C.astype(np.float64)))
exact_full_column_rank_certified = max_mod_rank == C.shape[1]

U, singular_values, _ = np.linalg.svd(C.astype(np.float64), full_matrices=True)
svd_rank = int(np.sum(singular_values > singular_values[0] * 1e-12))
left_null = U[:, svd_rank:].T
null_error = float(np.max(np.abs(left_null @ C))) if len(left_null) else 0.0

np.savetxt(out / "incidence.txt", C, fmt="%d")
np.save(out / "incidence.npy", C)
np.save(out / "left_nullspace_orthonormal.npy", left_null)

with (out / "nodes.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["node_id", "label", "kind", "signed_index", "mw_vector"])
    for node in nodes:
        writer.writerow(
            [
                node["node_id"],
                node["label"],
                node["kind"],
                node["signed_index"],
                vector_text(node["vector"]),
            ]
        )

with (out / "square_relations.tsv").open("w", newline="") as handle:
    fields = [
        "relation_id",
        "kind",
        "i",
        "j",
        "k",
        "node_a",
        "node_b",
        "node_c",
        "signed_a",
        "signed_b",
        "signed_c",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(relations)

left_dim = C.shape[0] - max_mod_rank
slope_variables = len(relations) * 3
x_variables = len(nodes) * 5
full_variables = slope_variables + x_variables
coefficient_equations = len(relations) * 5
eliminated_equations = left_dim * 5

summary = "\n".join(
    [
        "Coxeter-9 x-coordinate reconstruction reduction",
        "================================================",
        f"signed_vectors={len(signed)}",
        f"global_additive_triples={len(triples)}",
        f"clique_indices={' '.join(map(str, clique_indices))}",
        f"x_nodes={len(nodes)}",
        "vertex_nodes=9",
        f"difference_nodes={len(nodes)-9}",
        f"square_relations={len(relations)}",
        "vertex_difference_relations=36",
        "difference_triangle_relations=84",
        f"incidence_shape={C.shape[0]}x{C.shape[1]}",
        f"numeric_rank={numeric_rank}",
        f"mod_ranks={mod_ranks}",
        f"max_mod_rank={max_mod_rank}",
        f"exact_full_column_rank_certified={exact_full_column_rank_certified}",
        f"left_nullity={left_dim}",
        f"svd_rank={svd_rank}",
        f"left_nullspace_error={null_error:.3e}",
        "",
        f"x_quartic_variables={x_variables}",
        f"quadratic_slope_variables={slope_variables}",
        f"full_x_plus_slope_variables={full_variables}",
        f"raw_coefficient_equations={coefficient_equations}",
        f"after_linear_x_elimination_quadratic_equations={eliminated_equations}",
        "",
        "For slope m_r = a_r U^2 + b_r U V + c_r V^2,",
        "each left-kernel vector lambda of C imposes the binary-quartic identity",
        "",
        "    sum_r lambda_r * m_r(U,V)^2 = 0.",
        "",
        "Once slopes are recovered, every coefficient vector of the 45 quartic",
        "x-polynomials is obtained by solving C X_k = q_k(m).",
        "",
        "Expected geometric freedom before gauge fixing: 5 dimensions:",
        "  3 from PGL_2 on the base, 1 from Weierstrass scaling, 1 true moduli.",
    ]
) + "\n"

(out / "summary.txt").write_text(summary)

print(summary, end="")
print("saved =", out)
