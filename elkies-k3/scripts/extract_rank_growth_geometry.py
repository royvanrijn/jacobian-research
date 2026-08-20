from __future__ import annotations

from pathlib import Path
from time import perf_counter
import argparse
import csv
import math

import numpy as np
import csv
from sage.all import ZZ, matrix


parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("motifs")
parser.add_argument("cores")

parser.add_argument("--trials", type=int, default=100)
parser.add_argument("--top", type=int, default=20)
parser.add_argument("--local-radius", type=int, default=1)
parser.add_argument("--run-name", default="rank21-growth-geometry-v1")

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / args.run_name
OUT.mkdir(parents=True, exist_ok=True)

H = np.loadtxt(args.gram, dtype=np.float64)

def read_tsv(path):
    with Path(path).open() as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


motifs = read_tsv(args.motifs)
cores_rows = read_tsv(args.cores)

if H.shape != (21, 21):
    raise RuntimeError(H.shape)


def parse_vec(s):
    return np.asarray(
        [int(x) for x in str(s).split()],
        dtype=np.int64,
    )


def load_core(trial):
    rows = [
        row
        for row in cores_rows
        if int(row["trial"]) == trial
    ]

    rows.sort(
        key=lambda row: int(row["row"])
    )

    B = np.vstack(
        [
            parse_vec(row["vector"])
            for row in rows
        ]
    )

    if B.shape != (17,21):
        raise RuntimeError(
            f"trial {trial}: core shape {B.shape}"
        )

    return B


def exact_rank(B):
    return matrix(
        ZZ,
        B.tolist(),
    ).rank()


def nearest_grid_local(c, G, denominator, radius):
    """
    Search around componentwise nearest denominator-grid point.

    Returns invariant squared height distance.

        q in (1/k) Z^r
        d = c-q
        dist^2 = d^T G d
    """

    k = denominator

    center = np.rint(
        k * c
    ).astype(np.int64)

    r = len(c)

    best = None
    best_q = None

    # Full (2R+1)^r is impossible.
    #
    # Instead identify the coordinates whose rounding ambiguity
    # contributes most metric energy and branch only there.
    fractional = (
        k * c
        -
        center
    )

    # crude but deterministic importance proxy
    diag = np.diag(G)

    importance = (
        fractional * fractional
        * diag
    )

    order = np.argsort(
        -importance
    )

    # Branch over at most 8 ambiguous coordinates.
    active = order[:min(8,r)]

    offsets = [-radius, 0, radius]

    candidates = [
        center.copy()
    ]

    for j in active:
        old = candidates
        candidates = []

        for z in old:
            for off in offsets:
                zz = z.copy()
                zz[j] += off
                candidates.append(zz)

    for z in candidates:

        q = z.astype(np.float64) / k

        d = c - q

        d2 = float(
            d @ G @ d
        )

        if (
            best is None
            or d2 < best
        ):
            best = d2
            best_q = q

    return (
        float(best),
        best_q,
    )


def short_vectors_from_basis(B, Hcur, max_vectors=2000):
    """
    Build a canonical-ish set of short vectors in the current span.

    Start with:
      ±e_i
      ±(e_i±e_j)

    then keep the shortest distinct vectors by height.

    This is deliberately cheap and invariant evaluation happens
    through the height metric.
    """

    r = len(B)

    coeffs = []

    for i in range(r):
        e = np.zeros(r, dtype=np.int64)
        e[i] = 1
        coeffs.append(e)
        coeffs.append(-e)

    for i in range(r):
        for j in range(i+1, r):

            for sj in (-1,1):

                e = np.zeros(r, dtype=np.int64)
                e[i] = 1
                e[j] = sj

                coeffs.append(e)
                coeffs.append(-e)

    C = np.asarray(
        coeffs,
        dtype=np.int64,
    )

    norms = np.einsum(
        "ij,ij->i",
        C @ Hcur,
        C,
    )

    idx = np.argsort(
        norms
    )

    C = C[idx]
    norms = norms[idx]

    seen = set()
    out_c = []
    out_n = []

    for c,n in zip(C,norms):

        key = tuple(
            map(int,c)
        )

        if key in seen:
            continue

        seen.add(key)

        out_c.append(c)
        out_n.append(float(n))

        if len(out_c) >= max_vectors:
            break

    return (
        np.asarray(out_c),
        np.asarray(out_n),
    )


FIELDS = [
    "trial",
    "slot",
    "shallow_rank",

    "orthogonal_height",
    "orthogonal_ratio",

    "nearest_L_d2",
    "nearest_halfL_d2",
    "nearest_thirdL_d2",
    "nearest_quarterL_d2",

    "nearest_half_advantage",
    "nearest_third_advantage",
    "nearest_quarter_advantage",

    "short_corr_max",
    "short_corr_mean_top10",
    "short_corr_mean_top50",

    "short_corr_gt_050",
    "short_corr_gt_070",
    "short_corr_gt_090",

    "short_effective_count",
    "short_entropy_norm",

    "new_block_projection_share",
    "last_added_projection_share",
    "new_block_max_share",
    "new_block_coeff_l2",

    "nearest_short_index",
    "nearest_short_corr",
    "nearest_short_norm",
]

rows = []

t0 = perf_counter()

trial_ids = sorted(
    {
        int(row["trial"])
        for row in motifs
    }
)[:args.trials]

for pos, trial in enumerate(
    trial_ids,
    start=1,
):

    B = load_core(
        trial
    )

    if exact_rank(B) != 17:
        raise RuntimeError(
            f"trial {trial}: bad initial rank"
        )

    trial_motifs = [
        row
        for row in motifs
        if int(row["trial"]) == trial
    ]

    for slot in range(1,5):

        current_rank = B.shape[0]

        Hcur = (
            B @ H
        ) @ B.T

        short_c, short_norms = short_vectors_from_basis(
            B,
            Hcur,
        )

        x = [
            row
            for row in trial_motifs
            if (
                int(row["slot"]) == slot
                and
                int(row["shallow_rank"]) <= args.top
            )
        ]

        x.sort(
            key=lambda row: int(row["shallow_rank"])
        )

        if not len(x):
            raise RuntimeError(
                f"trial={trial} slot={slot}: missing motifs"
            )

        for row in x:

            c = np.asarray(
                [
                    float(v)
                    for v in row["projection_coeffs"].split()
                ],
                dtype=np.float64,
            )

            if len(c) != current_rank:
                raise RuntimeError(
                    f"trial={trial} slot={slot}: coeff dimension "
                    f"{len(c)} != {current_rank}"
                )

            proj_norm = float(
                c @ Hcur @ c
            )

            d1, q1 = nearest_grid_local(
                c,
                Hcur,
                1,
                args.local_radius,
            )

            d2, q2 = nearest_grid_local(
                c,
                Hcur,
                2,
                args.local_radius,
            )

            d3, q3 = nearest_grid_local(
                c,
                Hcur,
                3,
                args.local_radius,
            )

            d4, q4 = nearest_grid_local(
                c,
                Hcur,
                4,
                args.local_radius,
            )

            # ------------------------------------------------
            # Sequential rank-growth block geometry
            #
            # B[:17]  = original rank-17 core
            # B[17:]  = previously selected +1/+2/+3 directions
            #
            # We ask how much of this candidate's projection onto
            # the current span lies specifically in the NEW
            # directions after quotienting out the original core.
            #
            # This tests:
            #
            #   independent spokes
            #          versus
            #   18 -> 19 -> 20 -> 21 cascade
            #
            # using the Schur complement of the original 17-block.
            # ------------------------------------------------

            if current_rank > 17:

                G00 = Hcur[:17, :17]
                G01 = Hcur[:17, 17:]
                G10 = Hcur[17:, :17]
                G11 = Hcur[17:, 17:]

                # Gram matrix of the previously-added directions
                # modulo the original rank-17 span.
                solve00_G01 = np.linalg.solve(
                    G00,
                    G01,
                )

                S = (
                    G11
                    -
                    G10 @ solve00_G01
                )

                # Pairing of candidate projection with the
                # original core and added directions.
                b_all = (
                    Hcur
                    @ c
                )

                b0 = b_all[:17]
                b1 = b_all[17:]

                # Pairings with the residualized +1/+2/+3 block.
                solve00_b0 = np.linalg.solve(
                    G00,
                    b0,
                )

                residual_pair = (
                    b1
                    -
                    G10 @ solve00_b0
                )

                # Coordinates in the residualized new block.
                alpha = np.linalg.solve(
                    S,
                    residual_pair,
                )

                new_block_norm = float(
                    residual_pair
                    @ alpha
                )

                if (
                    new_block_norm < 0
                    and
                    new_block_norm > -1e-8
                ):
                    new_block_norm = 0.0

                new_block_projection_share = (
                    new_block_norm / proj_norm
                    if proj_norm > 1e-30
                    else 0.0
                )

                # ------------------------------------------------
                # Internal distribution across the previously
                # added directions.
                #
                # S = L L^T
                # z = L^T alpha
                # ||z||^2 = alpha^T S alpha
                #
                # Note: these within-block shares depend on the
                # sequential ordering of the added directions,
                # which is exactly what we want to test here.
                # ------------------------------------------------

                try:
                    LS = np.linalg.cholesky(
                        S
                    )

                    znew = (
                        LS.T
                        @ alpha
                    )

                    energy_new = (
                        znew * znew
                    )

                    energy_new_sum = float(
                        energy_new.sum()
                    )

                    if energy_new_sum > 1e-30:

                        shares_new = (
                            energy_new
                            /
                            energy_new_sum
                        )

                        new_block_max_share = float(
                            shares_new.max()
                        )

                        last_added_projection_share = float(
                            shares_new[-1]
                        )

                    else:

                        new_block_max_share = 0.0
                        last_added_projection_share = 0.0

                except np.linalg.LinAlgError:

                    # Should be rare; retain total block share,
                    # but do not invent an internal decomposition.
                    new_block_max_share = float("nan")
                    last_added_projection_share = float("nan")

                new_block_coeff_l2 = float(
                    np.linalg.norm(
                        alpha
                    )
                )

            else:

                new_block_projection_share = 0.0
                last_added_projection_share = 0.0
                new_block_max_share = 0.0
                new_block_coeff_l2 = 0.0

            # ------------------------------------------------
            # Intrinsic correlation with short lattice vectors
            # ------------------------------------------------

            pair = (
                short_c
                @ Hcur
                @ c
            )

            denom = np.sqrt(
                np.maximum(
                    short_norms
                    *
                    proj_norm,
                    1e-30,
                )
            )

            corr = np.abs(
                pair
            ) / denom

            corr = np.clip(
                corr,
                0,
                1,
            )

            sorted_corr = np.sort(
                corr
            )[::-1]

            # Convert squared correlation strengths into weights.
            strength = corr * corr
            ssum = float(
                strength.sum()
            )

            if ssum > 0:

                w = strength / ssum

                effective = float(
                    1.0
                    /
                    np.sum(
                        w*w
                    )
                )

                nz = w[w > 1e-300]

                entropy = float(
                    -np.sum(
                        nz*np.log(nz)
                    )
                )

                entropy_norm = float(
                    entropy
                    /
                    np.log(len(w))
                )

            else:
                effective = 0.0
                entropy_norm = 0.0

            nearest_idx = int(
                np.argmax(corr)
            )

            rows.append({
                "trial": trial,
                "slot": slot,
                "shallow_rank": int(row["shallow_rank"]),

                "orthogonal_height":
                    float(row["orthogonal_height"]),

                "orthogonal_ratio":
                    float(row["orthogonal_ratio"]),

                "nearest_L_d2": d1,
                "nearest_halfL_d2": d2,
                "nearest_thirdL_d2": d3,
                "nearest_quarterL_d2": d4,

                "nearest_half_advantage":
                    d2 - d1,

                "nearest_third_advantage":
                    d3 - d1,

                "nearest_quarter_advantage":
                    d4 - min(d1,d2),

                "short_corr_max":
                    float(sorted_corr[0]),

                "short_corr_mean_top10":
                    float(
                        sorted_corr[
                            :min(10,len(sorted_corr))
                        ].mean()
                    ),

                "short_corr_mean_top50":
                    float(
                        sorted_corr[
                            :min(50,len(sorted_corr))
                        ].mean()
                    ),

                "short_corr_gt_050":
                    int(
                        np.count_nonzero(
                            corr >= .50
                        )
                    ),

                "short_corr_gt_070":
                    int(
                        np.count_nonzero(
                            corr >= .70
                        )
                    ),

                "short_corr_gt_090":
                    int(
                        np.count_nonzero(
                            corr >= .90
                        )
                    ),

                "short_effective_count":
                    effective,

                "short_entropy_norm":
                    entropy_norm,

                "new_block_projection_share":
                    new_block_projection_share,

                "last_added_projection_share":
                    last_added_projection_share,

                "new_block_max_share":
                    new_block_max_share,

                "new_block_coeff_l2":
                    new_block_coeff_l2,

                "nearest_short_index":
                    nearest_idx,

                "nearest_short_corr":
                    float(
                        corr[nearest_idx]
                    ),

                "nearest_short_norm":
                    float(
                        short_norms[nearest_idx]
                    ),
            })

        # rank-1 vector is the sequential extension.
        winner = next(
            row
            for row in x
            if int(row["shallow_rank"]) == 1
        )

        winner_vec = parse_vec(
            winner["vector"]
        )

        B = np.vstack(
            [
                B,
                winner_vec,
            ]
        )

    if exact_rank(B) != 21:
        raise RuntimeError(
            f"trial={trial}: final rank !=21"
        )

    if (
        pos % 10 == 0
        or pos == len(trial_ids)
    ):
        print(
            f"PROGRESS"
            f"|trial={pos}/{len(trial_ids)}"
            f"|rows={len(rows)}"
            f"|seconds={perf_counter()-t0:.1f}",
            flush=True,
        )


path = (
    OUT
    / "growth_geometry.tsv"
)

if rows:
    with path.open(
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=list(rows[0].keys()),
        )

        writer.writeheader()
        writer.writerows(rows)

print()
print("="*72)
print("DONE")
print("rows =",len(rows))
print("saved =",path)
print("seconds =",perf_counter()-t0)
