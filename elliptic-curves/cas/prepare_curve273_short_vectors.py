#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations, product
from fractions import Fraction
from math import gcd
import argparse
import heapq
import json
import numpy as np

from sage.all import ZZ, matrix

from search_extra_points import gp_rational, gp_vector, run_gp


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

    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(f"expected a square Gram matrix, got {H.shape}")

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


def generate_height_gram(curve_id, path, curve_json=None):
    if curve_json is not None:
        payload = json.loads(curve_json.read_text(encoding="utf-8"))
        if int(payload["id"]) != curve_id:
            raise ValueError("--curve-id does not match --curve-json")
        coefficients = tuple(Fraction(value) for value in payload["ainvs"])
        public_points = tuple(
            tuple(Fraction(value) for value in point)
            for point in payload["points"]
        )
    elif curve_id == 273:
        import icarm_curve273 as data
        coefficients = data.GENERAL_WEIERSTRASS_COEFFICIENTS
        public_points = data.POINTS
    elif curve_id == 245:
        import icarm_curve245 as data
        coefficients = data.GENERAL_WEIERSTRASS_COEFFICIENTS
        public_points = data.POINTS
    elif curve_id == 90:
        import icarm_curve90 as data
        coefficients = data.GENERAL_WEIERSTRASS_COEFFICIENTS
        public_points = data.SEARCH_POINTS
    else:
        raise ValueError("--curve-json is required for this ICARM curve id")

    curve = ",".join(
        gp_rational(value) for value in coefficients
    )
    points = ",".join(gp_vector(point) for point in public_points)
    program = "\n".join(
        (
            "default(realprecision,120);",
            f"E=ellinit([{curve}]);",
            f"P=[{points}];",
            "H=ellheightmatrix(E,P);",
            'print("HEIGHT_BEGIN");',
            "for(i=1,matsize(H)[1],print(Vec(H[i,])));",
            'print("HEIGHT_END");',
            "quit",
        )
    ) + "\n"
    output, _wall = run_gp(program, timeout=120.0, stack_bytes=1_000_000_000)
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    start = lines.index("HEIGHT_BEGIN") + 1
    end = lines.index("HEIGHT_END")
    rows = lines[start:end]
    if len(rows) != len(public_points):
        raise AssertionError("PARI returned the wrong number of height rows")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(
        f"{PROTOCOL}|stage=height_gram|curve_id={curve_id}"
        f"|rows={len(rows)}|output={path}",
        flush=True,
    )


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

    ap.add_argument("--curve-id", type=int, default=273)
    ap.add_argument("--curve-json", type=Path)
    ap.add_argument("--gram", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--generate-gram", action="store_true")

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

    if args.gram is None:
        if args.curve_id == 273:
            args.gram = DEFAULT_GRAM
        elif args.curve_id == 245:
            args.gram = ROOT / "artifacts/local/elliptic-curves/curve245-rank20/height-gram.txt"
        else:
            raise SystemExit("--gram is required for this ICARM curve id")
    if args.out is None:
        if args.curve_id == 273:
            args.out = DEFAULT_OUT
        elif args.curve_id == 245:
            args.out = (
                ROOT
                / "artifacts/local/elliptic-curves/curve245-rank20/"
                "short-coefficient-vectors.tsv"
            )
        else:
            raise SystemExit("--out is required for this ICARM curve id")
    if args.generate_gram:
        generate_height_gram(args.curve_id, args.gram, args.curve_json)

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
    # in the original point basis.
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
