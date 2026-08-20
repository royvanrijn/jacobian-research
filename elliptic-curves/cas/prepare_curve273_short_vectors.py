#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations, product
from math import gcd
import argparse
import heapq
import numpy as np

from sage.all import ZZ, matrix


PROTOCOL = "R31SHORT"

ROOT = Path(__file__).resolve().parents[2]

DEFAULT_GRAM = (
    ROOT /
    "artifacts/local/elliptic-curves/curve273-rank30/height-gram.txt"
)

DEFAULT_OUT = (
    ROOT /
    "artifacts/local/elliptic-curves/curve273-rank30/"
    "short-coefficient-vectors.tsv"
)


def parse_height_gram(path):
    rows = []

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue

        if not (line.startswith("[") and line.endswith("]")):
            raise ValueError(f"unexpected Gram row: {line}")

        values = [
            float(x.strip())
            for x in line[1:-1].split(",")
        ]

        rows.append(values)

    H = np.asarray(rows, dtype=float)

    if H.shape != (30, 30):
        raise ValueError(f"expected 30x30 Gram, got {H.shape}")

    H = (H + H.T) / 2.0

    eig = np.linalg.eigvalsh(H)

    if eig[0] <= 0:
        raise ValueError(
            f"height Gram not positive definite: eig_min={eig[0]}"
        )

    print(
        f"{PROTOCOL}|stage=gram"
        f"|eig_min={eig[0]:.12g}"
        f"|eig_max={eig[-1]:.12g}"
        f"|condition={eig[-1]/eig[0]:.12g}",
        flush=True,
    )

    return H


def canonical(coeffs):
    coeffs = tuple(int(x) for x in coeffs)

    first = next(
        (x for x in coeffs if x != 0),
        None,
    )

    if first is None:
        return None

    if first < 0:
        coeffs = tuple(-x for x in coeffs)

    g = 0
    for x in coeffs:
        g = gcd(g, abs(x))

    if g > 1:
        coeffs = tuple(x // g for x in coeffs)

    return coeffs


def qheight(H, c):
    v = np.asarray(c, dtype=float)
    return float(v @ H @ v)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--gram", type=Path, default=DEFAULT_GRAM)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)

    ap.add_argument(
        "--count",
        type=int,
        default=1000,
        help="number of short vectors written",
    )

    ap.add_argument(
        "--working-pool",
        type=int,
        default=10000,
    )

    ap.add_argument(
        "--raw-weight",
        type=int,
        default=4,
        help="exhaust signed raw combinations through this weight",
    )

    ap.add_argument(
        "--lll-scale",
        type=int,
        default=1_000_000_000,
    )

    ap.add_argument(
        "--max-abs-coeff",
        type=int,
        default=32,
    )

    ap.add_argument(
        "--max-l1",
        type=int,
        default=128,
    )

    args = ap.parse_args()

    H = parse_height_gram(args.gram)
    n = H.shape[0]

    # --------------------------------------------------------
    # Approximate Euclidean embedding of the height lattice.
    #
    # H = L L^T, so
    #
    #     c H c^T = || c L ||^2.
    #
    # Round a large multiple of L to ZZ and use Sage LLL.
    # The returned transformation rows are coefficient vectors
    # in the ORIGINAL 30-point basis.
    # --------------------------------------------------------

    L = np.linalg.cholesky(H)

    embedded = [
        [
            int(round(args.lll_scale * float(L[i, j])))
            for j in range(n)
        ]
        for i in range(n)
    ]

    B = matrix(ZZ, embedded)

    reduced, transform = B.LLL(
        delta=0.999,
        transformation=True,
    )

    if transform * B != reduced:
        raise AssertionError("unexpected Sage LLL transformation convention")

    print(
        f"{PROTOCOL}|stage=lll"
        f"|det_transform={transform.det()}",
        flush=True,
    )

    # max-heap implemented with negative heights.
    heap = []
    seen = set()

    def retain(c):
        c = canonical(c)

        if c is None or c in seen:
            return

        # Do not retain original basis vectors; the previous
        # 1,935-chart pass already searched them.
        support = sum(x != 0 for x in c)
        if support <= 1:
            return

        if max(abs(x) for x in c) > args.max_abs_coeff:
            return

        if sum(abs(x) for x in c) > args.max_l1:
            return

        seen.add(c)

        h = qheight(H, c)

        item = (-h, c)

        if len(heap) < args.working_pool:
            heapq.heappush(heap, item)
        elif h < -heap[0][0]:
            heapq.heapreplace(heap, item)

    # --------------------------------------------------------
    # Exhaust small signed raw combinations.
    #
    # Global sign is canonicalized by fixing the first selected
    # coefficient to +1.
    # --------------------------------------------------------

    for weight in range(2, args.raw_weight + 1):
        tested = 0

        for indices in combinations(range(n), weight):
            for tails in product((-1, 1), repeat=weight - 1):
                c = [0] * n
                c[indices[0]] = 1

                for index, sign in zip(indices[1:], tails):
                    c[index] = sign

                retain(c)
                tested += 1

        print(
            f"{PROTOCOL}|stage=raw"
            f"|weight={weight}"
            f"|tested={tested}"
            f"|pool={len(heap)}",
            flush=True,
        )

    # --------------------------------------------------------
    # LLL rows and short combinations of the LLL basis.
    # --------------------------------------------------------

    U = [
        tuple(int(x) for x in transform.row(i))
        for i in range(n)
    ]

    for c in U:
        retain(c)

    # Pair combinations of all LLL rows.
    for i, j in combinations(range(n), 2):
        for sign in (-1, 1):
            retain(
                tuple(
                    U[i][k] + sign * U[j][k]
                    for k in range(n)
                )
            )

    # Triple combinations among the first 14 reduced rows.
    m = min(14, n)

    for i, j, k in combinations(range(m), 3):
        for sj, sk in product((-1, 1), repeat=2):
            retain(
                tuple(
                    U[i][q]
                    + sj * U[j][q]
                    + sk * U[k][q]
                    for q in range(n)
                )
            )

    result = sorted(
        ((-neg, c) for neg, c in heap),
        key=lambda item: (item[0], item[1]),
    )[:args.count]

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open("w") as f:
        f.write(
            "height\t"
            + "\t".join(f"c{i+1}" for i in range(n))
            + "\n"
        )

        for h, c in result:
            f.write(
                f"{h:.17g}\t"
                + "\t".join(str(x) for x in c)
                + "\n"
            )

    print(
        f"{PROTOCOL}|stage=done"
        f"|vectors={len(result)}"
        f"|best_height={result[0][0]:.17g}"
        f"|out={args.out}",
        flush=True,
    )

    for rank, (h, c) in enumerate(result[:20], 1):
        support = sum(x != 0 for x in c)

        print(
            f"{PROTOCOL}|short={rank}"
            f"|height={h:.12g}"
            f"|support={support}"
            f"|l1={sum(abs(x) for x in c)}"
            f"|coeffs={','.join(map(str,c))}",
            flush=True,
        )


if __name__ == "__main__":
    main()
