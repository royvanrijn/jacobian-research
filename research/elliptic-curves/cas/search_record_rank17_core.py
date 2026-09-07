#!/usr/bin/env python3
"""Search short-vector clouds for an exact rank-17 concentration.

The specialization of the 1,311 unoriented minimal vectors of R17 must lie in
one exact 17-dimensional rational subspace of the displayed Mordell--Weil
group.  This script enumerates short vectors with PARI/GP and uses deterministic
RANSAC to look for such a concentration.  A hit is only a numerical provenance
fingerprint: it is not an isometry or a specialization certificate.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import random
import subprocess

import numpy as np

from compare_record_height_lattices import CURVES, CurveData, gp_vector


def gp_program(curve: CurveData, digits: int, bound: float) -> str:
    points = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    return f"""
default(parisizemax,2000000000);
default(parisize,500000000);
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});
P={points};
H=ellheightmatrix(E,P);
U=qflllgram(H);
R=U~*H*U;
Q=qfminim(R,{bound},100000,2);
V=U*Q[3];
print("COUNT_SIGNED|",Q[1]);
print("BEGIN_VECTORS");
for(j=1,matsize(V)[2],for(i=1,matsize(V)[1],if(i>1,print1("|"));print1(V[i,j]));print());
print("END_VECTORS");
print("BEGIN_HEIGHTS");
for(j=1,matsize(V)[2],print(V[,j]~*H*V[,j]));
print("END_HEIGHTS");
"""


def parse_block(lines: list[str], begin: str, end: str) -> list[str]:
    i = lines.index(begin) + 1
    j = lines.index(end, i)
    return lines[i:j]


def enumerate_vectors(
    curve: CurveData, digits: int, bound: float
) -> tuple[np.ndarray, np.ndarray]:
    completed = subprocess.run(
        ["gp", "-q"],
        input=gp_program(curve, digits, bound),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    vectors = np.array(
        [[int(x) for x in line.split("|")] for line in parse_block(lines, "BEGIN_VECTORS", "END_VECTORS")],
        dtype=np.int64,
    )
    heights = np.array(
        [float(x) for x in parse_block(lines, "BEGIN_HEIGHTS", "END_HEIGHTS")],
        dtype=float,
    )
    order = np.argsort(heights)
    return vectors[order], heights[order]


def modular_rank(rows: np.ndarray, prime: int = 1_000_003) -> int:
    a = np.asarray(rows, dtype=object).astype(object)
    a = np.array([[int(x) % prime for x in row] for row in a], dtype=np.int64)
    rank = 0
    for col in range(a.shape[1]):
        pivot = next((i for i in range(rank, a.shape[0]) if a[i, col]), None)
        if pivot is None:
            continue
        a[[rank, pivot]] = a[[pivot, rank]]
        inv = pow(int(a[rank, col]), -1, prime)
        a[rank] = (a[rank] * inv) % prime
        for i in range(a.shape[0]):
            if i != rank and a[i, col]:
                a[i] = (a[i] - a[i, col] * a[rank]) % prime
        rank += 1
        if rank == a.shape[0]:
            break
    return rank


def exact_kernel(basis: np.ndarray) -> np.ndarray:
    matrix = "[" + ";".join(
        ",".join(str(int(x)) for x in row) for row in basis
    ) + "]"
    program = f"B={matrix};K=matkerint(B);for(i=1,matsize(K)[1],for(j=1,matsize(K)[2],if(j>1,print1(\"|\"));print1(K[i,j]));print());\n"
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    rows = [line.split("|") for line in completed.stdout.splitlines() if line.strip()]
    if not rows:
        return np.empty((basis.shape[1], 0), dtype=np.int64)
    return np.array([[int(x) for x in row] for row in rows], dtype=np.int64)


def exact_in_span_mask(vectors: np.ndarray, basis: np.ndarray) -> np.ndarray:
    kernel = exact_kernel(basis)
    products = np.asarray(vectors, dtype=object) @ np.asarray(kernel, dtype=object)
    return np.all(products == 0, axis=1)


def in_span_mask(vectors: np.ndarray, basis: np.ndarray) -> np.ndarray:
    # The vectors are small integers and the singular-value gap between exact
    # rank 17 and the other 12--14 directions is enormous.  Candidates are
    # subsequently certified by modular rank, so QR here is only a fast filter.
    _, singular, vh = np.linalg.svd(basis.astype(float), full_matrices=True)
    dimension = len(basis)
    if len(singular) < dimension or singular[dimension - 1] < 1e-9:
        return np.zeros(len(vectors), dtype=bool)
    null = vh[dimension:].T
    residual = np.max(np.abs(vectors @ null), axis=1)
    tentative = residual < 1e-7
    indices = np.flatnonzero(tentative)
    exact = np.zeros(len(vectors), dtype=bool)
    for i in indices:
        exact[i] = modular_rank(np.vstack((basis, vectors[i]))) == dimension
    return exact


def ransac(
    vectors: np.ndarray,
    heights: np.ndarray,
    *,
    trials: int,
    pool: int,
    seed: int,
    dimension: int,
    sampling_order: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    rng = random.Random(seed)
    pool = min(pool, len(vectors))
    if sampling_order is None:
        sampling_order = np.arange(len(vectors))
    sampling_indices = [int(i) for i in sampling_order[:pool]]
    best_basis = np.empty((0, vectors.shape[1]), dtype=np.int64)
    best_mask = np.zeros(len(vectors), dtype=bool)
    for trial in range(trials):
        indices = rng.sample(sampling_indices, dimension)
        basis = vectors[indices]
        if modular_rank(basis) != dimension:
            continue
        mask = in_span_mask(vectors, basis)
        if mask.sum() > best_mask.sum():
            best_basis = basis.copy()
            best_mask = mask
            selected_heights = heights[mask]
            print(
                f"BEST|trial={trial}|count={mask.sum()}"
                f"|height_min={selected_heights.min():.12g}"
                f"|height_max={selected_heights.max():.12g}",
                flush=True,
            )
    return best_basis, best_mask


def canonical_line(vector: np.ndarray) -> tuple[int, ...]:
    values = [int(x) for x in vector]
    first = next((x for x in values if x), 0)
    if first < 0:
        values = [-x for x in values]
    return tuple(values)


def additive_degrees(vectors: np.ndarray, pair_limit: int) -> np.ndarray:
    """Count visible ``a +/- b = c`` incidences in the short-vector cloud."""

    line_index = {canonical_line(vector): i for i, vector in enumerate(vectors)}
    degrees = np.zeros(len(vectors), dtype=np.int64)
    pair_limit = min(pair_limit, len(vectors))
    for i in range(pair_limit):
        a = vectors[i]
        for j in range(i):
            b = vectors[j]
            for candidate in (a + b, a - b):
                k = line_index.get(canonical_line(candidate))
                if k is not None:
                    degrees[i] += 1
                    degrees[j] += 1
                    degrees[k] += 1
    return degrees


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("label", choices=[curve.label for curve in CURVES])
    parser.add_argument("--bound", type=float, required=True)
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--dimension", type=int, default=17)
    parser.add_argument(
        "--dimensions",
        help="comma-separated dimension scan (overrides --dimension)",
    )
    parser.add_argument("--pool", type=int, default=1800)
    parser.add_argument(
        "--additive-pair-limit",
        type=int,
        default=0,
        help="rank sampling pool by visible additive relations; test pairs among this many shortest lines",
    )
    parser.add_argument("--seed", type=int, default=1729302)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    curve = next(curve for curve in CURVES if curve.label == args.label)
    vectors, heights = enumerate_vectors(curve, args.digits, args.bound)
    print(
        f"ENUM|label={curve.label}|lines={len(vectors)}"
        f"|height_min={heights.min():.12g}|height_max={heights.max():.12g}",
        flush=True,
    )
    sampling_order = None
    if args.additive_pair_limit:
        degrees = additive_degrees(vectors, args.additive_pair_limit)
        sampling_order = np.argsort(-degrees, kind="stable")
        print(
            "ADDITIVE|top_degrees="
            + ",".join(str(int(degrees[i])) for i in sampling_order[:20]),
            flush=True,
        )
    dimensions = (
        [int(value) for value in args.dimensions.split(",")]
        if args.dimensions
        else [args.dimension]
    )
    results = []
    for dimension in dimensions:
        basis, mask = ransac(
            vectors,
            heights,
            trials=args.trials,
            pool=args.pool,
            seed=args.seed + dimension,
            dimension=dimension,
            sampling_order=sampling_order,
        )
        if len(basis):
            mask = exact_in_span_mask(vectors, basis)
        rank = modular_rank(vectors[mask]) if mask.any() else 0
        print(f"DIMFINAL|dimension={dimension}|count={mask.sum()}|rank={rank}")
        results.append((dimension, basis, mask))
    dimension, basis, mask = results[-1]
    print(f"FINAL|dimension={dimension}|count={mask.sum()}|rank={modular_rank(vectors[mask]) if mask.any() else 0}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w") as handle:
            handle.write(f"label={curve.label}\n")
            handle.write(f"bound={args.bound}\n")
            handle.write(f"short_vector_lines={len(vectors)}\n")
            handle.write(f"best_subspace_count={int(mask.sum())}\n")
            handle.write(f"dimension={dimension}\n")
            handle.write("basis_rows_in_public_point_coordinates:\n")
            for row in basis:
                handle.write(" ".join(str(int(x)) for x in row) + "\n")


if __name__ == "__main__":
    main()
