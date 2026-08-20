from __future__ import annotations

from pathlib import Path
from time import perf_counter
import argparse
import csv
import json
import math

import numpy as np
from sage.all import ZZ, matrix


parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--trials", type=int, default=200)
parser.add_argument("--core-pool", type=int, default=50000)
parser.add_argument("--candidate-pool", type=int, default=100000)
parser.add_argument("--top", type=int, default=100)
parser.add_argument("--seed", type=int, default=212121)

parser.add_argument(
    "--replay-cores",
    default=None,
    help="Replay exact rank-17 cores from a previous core_vectors.tsv"
)

parser.add_argument("--run-name", default="rank21-growth-motifs-v1")
parser.add_argument("--progress-every", type=int, default=10)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]

OUT = (
    BASE
    / "results"
    / args.run_name
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

H = np.loadtxt(
    args.gram,
    dtype=np.float64,
)

VALL = np.load(
    args.vectors,
    mmap_mode="r",
)

NALL = np.load(
    args.norms,
    mmap_mode="r",
)

if H.shape != (21, 21):
    raise RuntimeError(
        f"Expected 21x21 Gram, got {H.shape}"
    )

if VALL.ndim != 2 or VALL.shape[1] != 21:
    raise RuntimeError(
        f"Expected Nx21 vectors, got {VALL.shape}"
    )

CORE_LIMIT = min(
    args.core_pool,
    len(VALL),
)

CAND_LIMIT = min(
    args.candidate_pool,
    len(VALL),
)

CORE_POOL = np.asarray(
    VALL[:CORE_LIMIT],
    dtype=np.int64,
)

CAND = np.asarray(
    VALL[:CAND_LIMIT],
    dtype=np.int64,
)

CAND_NORM = np.asarray(
    NALL[:CAND_LIMIT],
    dtype=np.float64,
)

print("precomputing CAND @ H ...", flush=True)

CAND_H = (
    CAND @ H
)

print("ready", flush=True)


def exact_rank(A):
    return matrix(
        ZZ,
        A.tolist(),
    ).rank()


def normalized_grid_score(c, k):
    x = k * c
    d = (
        x
        -
        np.rint(x)
    )

    return float(
        np.sqrt(
            np.mean(
                d * d
            )
        )
    )


def build_core(rng):

    while True:

        chosen = []
        used = set()

        attempts = 0

        while len(chosen) < 17:

            attempts += 1

            if attempts > 50000:
                break

            u = rng.random()

            idx = int(
                u * u * CORE_LIMIT
            )

            idx = min(
                idx,
                CORE_LIMIT - 1,
            )

            if idx in used:
                continue

            v = CORE_POOL[idx]

            A = np.asarray(
                chosen + [v],
                dtype=np.float64,
            )

            if (
                np.linalg.matrix_rank(A)
                !=
                len(chosen) + 1
            ):
                continue

            chosen.append(
                v.copy()
            )

            used.add(idx)

        if len(chosen) != 17:
            continue

        B = np.asarray(
            chosen,
            dtype=np.int64,
        )

        if exact_rank(B) == 17:
            return B


def geometry(B):

    G = (
        B @ H
    ) @ B.T

    # G = L L^T.  For projection coefficients c, the vector
    #
    #     z = L^T c
    #
    # represents the projection in an orthonormal Euclidean
    # realization of the height metric:
    #
    #     ||z||^2 = c^T G c.
    #
    # Unlike ordinary coeff_l2, the resulting energy-distribution
    # statistics are insensitive to simple rescaling/conditioning
    # of the integral basis.
    L = np.linalg.cholesky(G)

    C = (
        CAND_H
        @ B.T
    )

    coeff = np.linalg.solve(
        G,
        C.T,
    ).T

    projection_norm = np.einsum(
        "ij,ij->i",
        C,
        coeff,
    )

    orth = (
        CAND_NORM
        -
        projection_norm
    )

    orth[
        (orth < 0)
        &
        (orth > -1e-7)
    ] = 0.0

    ratio = np.full(
        len(orth),
        np.nan,
        dtype=np.float64,
    )

    good = (
        CAND_NORM > 0
    )

    ratio[good] = (
        orth[good]
        /
        CAND_NORM[good]
    )

    viable = np.flatnonzero(
        (orth > 1e-7)
        &
        np.isfinite(ratio)
    )

    return (
        G,
        L,
        C,
        coeff,
        projection_norm,
        orth,
        ratio,
        viable,
    )


FIELDS = [
    "trial",
    "slot",
    "rank_before",

    "shallow_rank",
    "pool_index",

    "norm",
    "projection_norm",
    "orthogonal_height",
    "orthogonal_ratio",

    "height_percentile",
    "ratio_percentile",

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

    "projection_max_share",
    "projection_effective_dim",
    "projection_entropy",
    "projection_entropy_norm",

    "pair_abs_mean",
    "pair_abs_max",
    "pair_corr_mean",
    "pair_corr_max",

    "support",
    "max_coeff",

    "vector",
    "projection_coeffs",
]


def load_replay_cores(path):
    """
    Load exact rank-17 cores written by an earlier motif run.

    Expected TSV:
        trial    row    vector
    """

    import csv

    cores = {}

    with Path(path).open() as f:
        reader = csv.DictReader(
            f,
            delimiter="\t",
        )

        for row in reader:
            trial = int(row["trial"])
            r = int(row["row"])

            v = np.asarray(
                [
                    int(x)
                    for x in row["vector"].split()
                ],
                dtype=np.int64,
            )

            cores.setdefault(
                trial,
                []
            ).append(
                (r, v)
            )

    result = []

    for trial in sorted(cores):
        rows = sorted(
            cores[trial],
            key=lambda x: x[0],
        )

        B = np.asarray(
            [
                v
                for _, v in rows
            ],
            dtype=np.int64,
        )

        if B.shape != (17, 21):
            raise RuntimeError(
                f"Replay core {trial} has shape {B.shape}"
            )

        if exact_rank(B) != 17:
            raise RuntimeError(
                f"Replay core {trial} does not have exact rank 17"
            )

        result.append(
            (trial, B)
        )

    return result


REPLAY_CORES = None

if args.replay_cores is not None:
    REPLAY_CORES = load_replay_cores(
        args.replay_cores
    )

    print(
        "replay cores =",
        len(REPLAY_CORES),
        flush=True,
    )


out_path = (
    OUT
    / "growth_motifs.tsv"
)

core_path = (
    OUT
    / "core_vectors.tsv"
)

manifest = {
    "version": 1,
    "trials": args.trials,
    "core_pool": CORE_LIMIT,
    "candidate_pool": CAND_LIMIT,
    "top": args.top,
    "seed": args.seed,
    "replay_cores": (
        None
        if args.replay_cores is None
        else str(Path(args.replay_cores).resolve())
    ),
    "gram": str(Path(args.gram).resolve()),
    "vectors": str(Path(args.vectors).resolve()),
    "norms": str(Path(args.norms).resolve()),
}

(
    OUT
    / "manifest.json"
).write_text(
    json.dumps(
        manifest,
        indent=2,
        sort_keys=True,
    )
    + "\n"
)


t0 = perf_counter()

with (
    out_path.open(
        "w",
        newline="",
    ) as fout,
    core_path.open(
        "w",
        newline="",
    ) as fcore
):

    writer = csv.DictWriter(
        fout,
        delimiter="\t",
        fieldnames=FIELDS,
    )

    writer.writeheader()

    core_writer = csv.writer(
        fcore,
        delimiter="\t",
    )

    core_writer.writerow(
        [
            "trial",
            "row",
            "vector",
        ]
    )

    if REPLAY_CORES is None:
        trial_source = [
            (trial, None)
            for trial in range(args.trials)
        ]
    else:
        trial_source = REPLAY_CORES[:args.trials]

    actual_trials = len(trial_source)

    for trial_position, (trial, replay_B) in enumerate(
        trial_source
    ):

        trial_start = perf_counter()

        trial_seed = (
            args.seed
            +
            trial * 1_000_003
        )

        rng = np.random.default_rng(
            trial_seed
        )

        if replay_B is None:
            B = build_core(
                rng
            )
        else:
            B = replay_B.copy()

        for r, v in enumerate(B):
            core_writer.writerow(
                [
                    trial,
                    r,
                    " ".join(
                        map(
                            str,
                            map(
                                int,
                                v,
                            ),
                        )
                    ),
                ]
            )

        chosen = []

        for slot in range(
            1,
            5
        ):

            (
                G,
                L,
                C,
                coeff,
                projection_norm,
                orth,
                ratio,
                viable,
            ) = geometry(B)

            if not len(viable):
                raise RuntimeError(
                    f"trial={trial} slot={slot}: no viable candidates"
                )

            # ---------------------------------------------
            # Rank all viable vectors by absolute
            # transverse height.
            # ---------------------------------------------

            order = viable[
                np.argsort(
                    orth[viable]
                )
            ]

            take = min(
                args.top,
                len(order),
            )

            top_ids = (
                order[:take]
            )

            height_sorted = np.sort(
                orth[viable]
            )

            ratio_sorted = np.sort(
                ratio[viable]
            )

            rank_before = (
                B.shape[0]
            )

            diagG = np.diag(
                G
            )

            for shallow_rank, idx in enumerate(
                top_ids,
                start=1,
            ):

                idx = int(idx)

                c = (
                    coeff[idx]
                )

                n = float(
                    CAND_NORM[idx]
                )

                proj = float(
                    projection_norm[idx]
                )

                o = float(
                    orth[idx]
                )

                rratio = float(
                    ratio[idx]
                )

                hpct = float(
                    np.searchsorted(
                        height_sorted,
                        o,
                        side="right",
                    )
                    /
                    len(height_sorted)
                )

                rpct = float(
                    np.searchsorted(
                        ratio_sorted,
                        rratio,
                        side="right",
                    )
                    /
                    len(ratio_sorted)
                )

                g1 = normalized_grid_score(
                    c,
                    1,
                )

                g2 = normalized_grid_score(
                    c,
                    2,
                )

                g3 = normalized_grid_score(
                    c,
                    3,
                )

                g4 = normalized_grid_score(
                    c,
                    4,
                )

                frac = (
                    c
                    -
                    np.rint(c)
                )

                # --------------------------------------------------
                # Basis-metric projection distribution
                # --------------------------------------------------

                z = (
                    L.T
                    @ c
                )

                energy = (
                    z * z
                )

                energy_sum = float(
                    np.sum(energy)
                )

                if energy_sum > 1e-30:

                    weights = (
                        energy
                        /
                        energy_sum
                    )

                    projection_max_share = float(
                        np.max(weights)
                    )

                    projection_effective_dim = float(
                        1.0
                        /
                        np.sum(
                            weights * weights
                        )
                    )

                    nz = weights[
                        weights > 1e-300
                    ]

                    projection_entropy = float(
                        -np.sum(
                            nz
                            *
                            np.log(nz)
                        )
                    )

                    projection_entropy_norm = float(
                        projection_entropy
                        /
                        np.log(len(weights))
                    )

                else:

                    projection_max_share = 1.0
                    projection_effective_dim = 1.0
                    projection_entropy = 0.0
                    projection_entropy_norm = 0.0

                pair_abs = np.abs(
                    C[idx]
                )

                corr = (
                    pair_abs
                    /
                    np.sqrt(
                        np.maximum(
                            n * diagG,
                            1e-30,
                        )
                    )
                )

                v = (
                    CAND[idx]
                )

                writer.writerow({
                    "trial": trial,
                    "slot": slot,
                    "rank_before": rank_before,

                    "shallow_rank": shallow_rank,
                    "pool_index": idx,

                    "norm": n,
                    "projection_norm": proj,
                    "orthogonal_height": o,
                    "orthogonal_ratio": rratio,

                    "height_percentile": hpct,
                    "ratio_percentile": rpct,

                    "grid1_score": g1,
                    "grid2_score": g2,
                    "grid3_score": g3,
                    "grid4_score": g4,

                    "half_advantage": g2 - g1,
                    "third_advantage": g3 - g1,
                    "quarter_advantage": g4 - min(g1, g2),

                    "coeff_l2": float(
                        np.linalg.norm(c)
                    ),

                    "frac_l2": float(
                        np.linalg.norm(frac)
                    ),

                    "frac_max": float(
                        np.max(
                            np.abs(frac)
                        )
                    ),

                    "projection_max_share":
                        projection_max_share,

                    "projection_effective_dim":
                        projection_effective_dim,

                    "projection_entropy":
                        projection_entropy,

                    "projection_entropy_norm":
                        projection_entropy_norm,

                    "pair_abs_mean": float(
                        np.mean(pair_abs)
                    ),

                    "pair_abs_max": float(
                        np.max(pair_abs)
                    ),

                    "pair_corr_mean": float(
                        np.mean(corr)
                    ),

                    "pair_corr_max": float(
                        np.max(corr)
                    ),

                    "support": int(
                        np.count_nonzero(v)
                    ),

                    "max_coeff": int(
                        np.max(
                            np.abs(v)
                        )
                    ),

                    "vector": " ".join(
                        map(
                            str,
                            map(
                                int,
                                v,
                            ),
                        )
                    ),

                    "projection_coeffs": " ".join(
                        f"{x:.15g}"
                        for x in c
                    ),
                })

            # ---------------------------------------------
            # The rank-1 shallow escape is the actual
            # sequential extension.
            # ---------------------------------------------

            selected = int(
                top_ids[0]
            )

            chosen.append(
                selected
            )

            B = np.vstack(
                [
                    B,
                    CAND[selected],
                ]
            )

        if exact_rank(B) != 21:
            raise RuntimeError(
                f"trial={trial}: final exact rank != 21"
            )

        fout.flush()
        fcore.flush()

        if (
            trial_position % args.progress_every == 0
            or trial_position + 1 == actual_trials
        ):

            print(
                f"PROGRESS"
                f"|trial={trial_position+1}/{actual_trials}"
                f"|source_trial={trial}"
                f"|chosen={chosen}"
                f"|trial_s={perf_counter()-trial_start:.3f}"
                f"|total_s={perf_counter()-t0:.1f}",
                flush=True,
            )


print()
print("=" * 72)
print("DONE")
print("motifs =", out_path)
print("cores  =", core_path)
print(
    "rows   =",
    (
        len(REPLAY_CORES[:args.trials])
        if REPLAY_CORES is not None
        else args.trials
    )
    * 4
    * args.top,
)
print(
    "seconds=",
    perf_counter() - t0,
)
