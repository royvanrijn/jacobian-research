from __future__ import annotations

from pathlib import Path
from time import perf_counter
import argparse
import csv
import math

import numpy as np
from sage.all import ZZ, matrix


parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--trials", type=int, default=200)
parser.add_argument("--core-pool", type=int, default=50000)
parser.add_argument("--hidden-pool", type=int, default=100000)
parser.add_argument("--seed", type=int, default=210021)
parser.add_argument("--out", default=None)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]

if args.out is None:
    OUT = BASE / "results" / "rank21-17plus4-lab"
else:
    OUT = Path(args.out)

OUT.mkdir(parents=True, exist_ok=True)

H = np.loadtxt(args.gram, dtype=np.float64)

VALL = np.load(
    args.vectors,
    mmap_mode="r",
)

NALL = np.load(
    args.norms,
    mmap_mode="r",
)

if H.shape != (21, 21):
    raise RuntimeError(f"Expected 21x21 Gram, got {H.shape}")

if VALL.ndim != 2 or VALL.shape[1] != 21:
    raise RuntimeError(f"Expected Nx21 vectors, got {VALL.shape}")

if len(VALL) != len(NALL):
    raise RuntimeError("vector/norm count mismatch")

rng = np.random.default_rng(args.seed)

core_limit = min(
    args.core_pool,
    len(VALL),
)

hidden_limit = min(
    args.hidden_pool,
    len(VALL),
)

CORE_POOL = np.asarray(
    VALL[:core_limit],
    dtype=np.int64,
)

CORE_NORMS = np.asarray(
    NALL[:core_limit],
    dtype=np.float64,
)

HIDDEN_POOL = np.asarray(
    VALL[:hidden_limit],
    dtype=np.int64,
)

HIDDEN_NORMS = np.asarray(
    NALL[:hidden_limit],
    dtype=np.float64,
)

print("RANK21 17+4 LAB")
print("full pool       =", len(VALL))
print("core pool       =", len(CORE_POOL))
print("hidden pool     =", len(HIDDEN_POOL))
print("core norm range =", float(CORE_NORMS[0]), "..", float(CORE_NORMS[-1]))
print("hidden range    =", float(HIDDEN_NORMS[0]), "..", float(HIDDEN_NORMS[-1]))
print()


def exact_rank(rows: list[np.ndarray]) -> int:
    if not rows:
        return 0

    return matrix(
        ZZ,
        np.asarray(rows, dtype=np.int64).tolist(),
    ).rank()


def saturation_index(B_np: np.ndarray) -> int:
    """
    For a rank-r row lattice B <= Z^21, the product of the nonzero
    Smith invariants is [sat(B):B].

    1 means primitive.
    """

    B = matrix(
        ZZ,
        B_np.tolist(),
    )

    try:
        result = B.smith_form()

        if isinstance(result, tuple):
            D = result[0]
        else:
            D = result

        product = 1
        count = 0

        for i in range(min(D.nrows(), D.ncols())):
            x = abs(int(D[i, i]))

            if x:
                product *= x
                count += 1

        if count != B.rank():
            return -2

        return int(product)

    except Exception as e:
        print(
            "SMITH_ERROR",
            repr(e),
            flush=True,
        )
        return -1


def normalized_grid_score(c: np.ndarray, k: int) -> float:
    """
    RMS distance to nearest (1/k)Z^17 point, measured in units
    of one denominator-k grid cell.

    This makes scores for k=1,2,3,4 directly comparable.
    """

    x = k * c
    delta = x - np.rint(x)

    return float(
        np.sqrt(
            np.mean(
                delta * delta
            )
        )
    )


def grid_metrics(c: np.ndarray, G: np.ndarray, k: int):
    """
    Distance of projection coefficients c to nearest (1/k) Z^17
    coordinate grid.

    This is deliberately numerical: the canonical-height Gram is
    numerical, so these are fingerprints, not claims about exact
    discriminant-group denominators.
    """

    q = np.rint(k * c) / k
    d = c - q

    euclidean = float(
        np.linalg.norm(d)
    )

    metric2 = float(
        d @ G @ d
    )

    if metric2 < 0 and metric2 > -1e-8:
        metric2 = 0.0

    metric = (
        math.sqrt(metric2)
        if metric2 >= 0
        else float("nan")
    )

    return euclidean, metric


def build_core() -> np.ndarray:
    """
    Construct a random independent 17-vector core.

    Bias strongly toward shorter vectors while retaining diversity.
    """

    chosen: list[np.ndarray] = []

    # Draw ranks with a squared distribution:
    # lots of very short vectors, but not exclusively the first few.
    attempts = 0

    while len(chosen) < 17:
        attempts += 1

        if attempts > 20000:
            raise RuntimeError("Could not build rank-17 core")

        u = rng.random()
        idx = int(
            (u * u) * core_limit
        )

        idx = min(
            idx,
            core_limit - 1,
        )

        candidate = CORE_POOL[idx]

        if any(
            np.array_equal(candidate, x)
            for x in chosen
        ):
            continue

        trial = chosen + [candidate]

        if exact_rank(trial) == len(trial):
            chosen.append(
                candidate.copy()
            )

    return np.asarray(
        chosen,
        dtype=np.int64,
    )


def analyze_core(
    core_id: int,
    B: np.ndarray,
):
    """
    Find four short independent directions extending B.

    Candidates are ranked by orthogonal residual height relative
    to span(B). Thus we deliberately ask which ambient vectors
    sit closest "above" the 17-dimensional core.
    """

    G = (B @ H) @ B.T

    sign, logdet = np.linalg.slogdet(G)

    if sign <= 0:
        raise RuntimeError("non-positive core Gram")

    eig = np.linalg.eigvalsh(G)

    Ginv = np.linalg.inv(G)

    # Pair every hidden-pool vector against the core.
    #
    # C[j,i] = <candidate_j, B_i>
    #
    C = (HIDDEN_POOL @ H) @ B.T

    # Projection coefficient rows:
    #
    # c = G^-1 b
    #
    COEFF = C @ Ginv

    # projection norm = b^T G^-1 b
    projection_norm = np.einsum(
        "ij,ij->i",
        C,
        COEFF,
    )

    residual = (
        HIDDEN_NORMS
        -
        projection_norm
    )

    # Numerical roundoff.
    residual[
        (residual < 0)
        &
        (residual > -1e-7)
    ] = 0.0

    # Positive orthogonal component means outside real span(B).
    viable = np.flatnonzero(
        residual > 1e-7
    )

    # Distribution of relative orthogonal heights among all viable
    # ambient candidates. This tells us whether the selected missing
    # direction is genuinely exceptional or merely the result of
    # deliberately minimizing orthogonal height.
    ratio_all = np.full(
        len(residual),
        np.nan,
        dtype=np.float64,
    )

    positive_norm = HIDDEN_NORMS > 0

    ratio_all[positive_norm] = (
        residual[positive_norm]
        /
        HIDDEN_NORMS[positive_norm]
    )

    valid_ratio = ratio_all[
        (residual > 1e-7)
        &
        np.isfinite(ratio_all)
    ]

    ratio_sorted = np.sort(
        valid_ratio
    )

    # Favor low absolute residual first.
    order = viable[
        np.argsort(
            residual[viable]
        )
    ]

    chosen_hidden: list[np.ndarray] = []
    chosen_rows = []

    current_rows = [
        row.copy()
        for row in B
    ]

    current_rank = 17

    for idx in order:

        v = HIDDEN_POOL[idx]

        trial_rows = (
            current_rows
            +
            [v]
        )

        r = exact_rank(
            trial_rows
        )

        if r != current_rank + 1:
            continue

        slot = len(chosen_hidden)

        c = COEFF[idx]

        n = float(
            HIDDEN_NORMS[idx]
        )

        proj = float(
            projection_norm[idx]
        )

        res = float(
            residual[idx]
        )

        residual_ratio = (
            res / n
        )

        residual_percentile = (
            np.searchsorted(
                ratio_sorted,
                residual_ratio,
                side="right",
            )
            /
            len(ratio_sorted)
        )

        grid1 = normalized_grid_score(c, 1)
        grid2 = normalized_grid_score(c, 2)
        grid3 = normalized_grid_score(c, 3)
        grid4 = normalized_grid_score(c, 4)

        # Negative values mean denominator-k fits better than the
        # coarser comparison grid.
        half_advantage = grid2 - grid1
        third_advantage = grid3 - grid1
        quarter_advantage = grid4 - min(grid1, grid2)

        integer_e, integer_m = grid_metrics(
            c,
            G,
            1,
        )

        half_e, half_m = grid_metrics(
            c,
            G,
            2,
        )

        third_e, third_m = grid_metrics(
            c,
            G,
            3,
        )

        quarter_e, quarter_m = grid_metrics(
            c,
            G,
            4,
        )

        frac = (
            c
            -
            np.rint(c)
        )

        pair_abs = np.abs(
            C[idx]
        )

        normalized_pair = (
            pair_abs
            /
            np.sqrt(
                np.maximum(
                    n * np.diag(G),
                    1e-30,
                )
            )
        )

        chosen_rows.append({
            "core_id": core_id,
            "slot": slot + 1,
            "pool_index": int(idx),

            "norm": n,
            "projection_norm": proj,
            "orthogonal_height": res,
            "orthogonal_ratio": residual_ratio,
            "orthogonal_percentile": residual_percentile,

            "grid1_score": grid1,
            "grid2_score": grid2,
            "grid3_score": grid3,
            "grid4_score": grid4,

            "half_advantage": half_advantage,
            "third_advantage": third_advantage,
            "quarter_advantage": quarter_advantage,

            "coeff_l2": float(np.linalg.norm(c)),
            "frac_l2": float(np.linalg.norm(frac)),
            "frac_max": float(np.max(np.abs(frac))),

            "integer_grid_euclid": integer_e,
            "integer_grid_metric": integer_m,

            "half_grid_euclid": half_e,
            "half_grid_metric": half_m,

            "third_grid_euclid": third_e,
            "third_grid_metric": third_m,

            "quarter_grid_euclid": quarter_e,
            "quarter_grid_metric": quarter_m,

            "pair_abs_mean": float(np.mean(pair_abs)),
            "pair_abs_max": float(np.max(pair_abs)),

            "pair_corr_mean": float(np.mean(normalized_pair)),
            "pair_corr_max": float(np.max(normalized_pair)),

            "support": int(np.count_nonzero(v)),
            "max_coeff": int(np.max(np.abs(v))),

            "vector": " ".join(
                map(str, map(int, v))
            ),

            "projection_coeffs": " ".join(
                f"{x:.12g}"
                for x in c
            ),
        })

        chosen_hidden.append(
            v.copy()
        )

        current_rows.append(
            v.copy()
        )

        current_rank += 1

        if len(chosen_hidden) == 4:
            break

    if len(chosen_hidden) != 4:
        return None, []

    sat = saturation_index(B)

    core_norms = np.einsum(
        "ij,ij->i",
        B @ H,
        B,
    )

    summary = {
        "core_id": core_id,

        "saturation_index": sat,

        "logdet": float(logdet),

        "eig_min": float(eig[0]),
        "eig_max": float(eig[-1]),
        "condition": float(eig[-1] / eig[0]),

        "core_norm_min": float(np.min(core_norms)),
        "core_norm_mean": float(np.mean(core_norms)),
        "core_norm_max": float(np.max(core_norms)),

        "hidden_orthogonal_sum": float(
            sum(
                row["orthogonal_height"]
                for row in chosen_rows
            )
        ),

        "hidden_orthogonal_max": float(
            max(
                row["orthogonal_height"]
                for row in chosen_rows
            )
        ),

        "hidden_half_metric_min": float(
            min(
                row["half_grid_metric"]
                for row in chosen_rows
            )
        ),

        "hidden_half_metric_mean": float(
            np.mean([
                row["half_grid_metric"]
                for row in chosen_rows
            ])
        ),

        "hidden_integer_metric_mean": float(
            np.mean([
                row["integer_grid_metric"]
                for row in chosen_rows
            ])
        ),

        "core_vectors": " | ".join(
            " ".join(map(str, map(int, row)))
            for row in B
        ),
    }

    return summary, chosen_rows


CORE_FIELDS = [
    "core_id",
    "saturation_index",
    "logdet",
    "eig_min",
    "eig_max",
    "condition",
    "core_norm_min",
    "core_norm_mean",
    "core_norm_max",
    "hidden_orthogonal_sum",
    "hidden_orthogonal_max",
    "hidden_half_metric_min",
    "hidden_half_metric_mean",
    "hidden_integer_metric_mean",
    "core_vectors",
]

HIDDEN_FIELDS = [
    "core_id",
    "slot",
    "pool_index",

    "norm",
    "projection_norm",
    "orthogonal_height",
    "orthogonal_ratio",
    "orthogonal_percentile",

    "grid1_score",
    "grid2_score",
    "grid3_score",
    "grid4_score",

    "half_advantage",
    "third_advantage",
    "quarter_advantage",

    "coeff_l2",
    "frac_l2",
    "frac_max",

    "integer_grid_euclid",
    "integer_grid_metric",

    "half_grid_euclid",
    "half_grid_metric",

    "third_grid_euclid",
    "third_grid_metric",

    "quarter_grid_euclid",
    "quarter_grid_metric",

    "pair_abs_mean",
    "pair_abs_max",

    "pair_corr_mean",
    "pair_corr_max",

    "support",
    "max_coeff",

    "vector",
    "projection_coeffs",
]


core_path = (
    OUT /
    "core_summary.tsv"
)

hidden_path = (
    OUT /
    "hidden_directions.tsv"
)

best_path = (
    OUT /
    "best_hidden_directions.tsv"
)

t0 = perf_counter()

core_rows = []
hidden_rows = []

best_residual = float("inf")
best_core = None


for trial in range(args.trials):

    B = build_core()

    summary, hidden = analyze_core(
        trial,
        B,
    )

    if summary is None:
        print(
            f"TRIAL|{trial}"
            f"|status=no_completion",
            flush=True,
        )
        continue

    core_rows.append(
        summary
    )

    hidden_rows.extend(
        hidden
    )

    score = summary[
        "hidden_orthogonal_sum"
    ]

    if score < best_residual:

        best_residual = score
        best_core = (
            summary,
            hidden,
        )

        print(
            f"BEST"
            f"|trial={trial}"
            f"|sat={summary['saturation_index']}"
            f"|orth_sum={score:.12g}"
            f"|orth_max={summary['hidden_orthogonal_max']:.12g}"
            f"|half_mean={summary['hidden_half_metric_mean']:.12g}",
            flush=True,
        )

        for row in hidden:
            print(
                "  HIDDEN"
                f"|slot={row['slot']}"
                f"|norm={row['norm']:.8f}"
                f"|orth={row['orthogonal_height']:.8f}"
                f"|ratio={row['orthogonal_ratio']:.8f}"
                f"|pct={100*row['orthogonal_percentile']:.6f}%"
                f"|g1={row['grid1_score']:.6f}"
                f"|g2={row['grid2_score']:.6f}"
                f"|g3={row['grid3_score']:.6f}"
                f"|g4={row['grid4_score']:.6f}"
                f"|halfAdv={row['half_advantage']:.6f}"
                f"|quarterAdv={row['quarter_advantage']:.6f}"
                f"|support={row['support']}"
                f"|maxc={row['max_coeff']}",
                flush=True,
            )

    if (
        trial % 10 == 0
        or trial + 1 == args.trials
    ):
        print(
            f"PROGRESS"
            f"|trial={trial+1}/{args.trials}"
            f"|valid={len(core_rows)}"
            f"|seconds={perf_counter()-t0:.2f}",
            flush=True,
        )


with core_path.open(
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=CORE_FIELDS,
    )

    writer.writeheader()
    writer.writerows(
        core_rows
    )


with hidden_path.open(
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=HIDDEN_FIELDS,
    )

    writer.writeheader()
    writer.writerows(
        hidden_rows
    )


if best_core is not None:

    summary, hidden = best_core

    with best_path.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=HIDDEN_FIELDS,
        )

        writer.writeheader()
        writer.writerows(
            hidden
        )


print()
print("=" * 72)
print("FINAL")
print("valid cores =", len(core_rows))
print("hidden rows =", len(hidden_rows))
print("seconds =", perf_counter() - t0)

print()
print("saved:")
print(core_path)
print(hidden_path)
print(best_path)

if best_core is not None:

    summary, hidden = best_core

    print()
    print("BEST CORE")
    print("core_id =", summary["core_id"])
    print("saturation_index =", summary["saturation_index"])
    print("hidden orthogonal sum =", summary["hidden_orthogonal_sum"])
    print("hidden orthogonal max =", summary["hidden_orthogonal_max"])
    print("half-grid mean =", summary["hidden_half_metric_mean"])

    print()
    print("FOUR MISSING DIRECTIONS")

    for row in hidden:

        print(
            f"slot={row['slot']}"
            f" pool={row['pool_index']}"
            f" norm={row['norm']:.12g}"
            f" orth={row['orthogonal_height']:.12g}"
            f" ratio={row['orthogonal_ratio']:.12g}"
            f" pct={100*row['orthogonal_percentile']:.8g}%"
            f" g1={row['grid1_score']:.8g}"
            f" g2={row['grid2_score']:.8g}"
            f" g3={row['grid3_score']:.8g}"
            f" g4={row['grid4_score']:.8g}"
            f" halfAdv={row['half_advantage']:.8g}"
            f" thirdAdv={row['third_advantage']:.8g}"
            f" quarterAdv={row['quarter_advantage']:.8g}"
        )

