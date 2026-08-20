from __future__ import annotations

from pathlib import Path
from fractions import Fraction
from time import perf_counter
import argparse
import csv
import math

import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("motifs")
parser.add_argument("cores")

parser.add_argument("--trials", type=int, default=200)
parser.add_argument("--max-denominator", type=int, default=12)
parser.add_argument("--run-name", default="rank21-growth-chain-v1")

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / args.run_name
OUT.mkdir(parents=True, exist_ok=True)

H = np.loadtxt(
    args.gram,
    dtype=np.float64,
)


def read_tsv(path):
    with Path(path).open() as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


motifs = read_tsv(
    args.motifs
)

cores_rows = read_tsv(
    args.cores
)


def parse_vec(s):
    return np.asarray(
        [
            int(x)
            for x in str(s).split()
        ],
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
            f"trial={trial}: core shape {B.shape}"
        )

    return B


def winner_for(trial, slot):
    rows = [
        row
        for row in motifs
        if (
            int(row["trial"]) == trial
            and
            int(row["slot"]) == slot
            and
            int(row["shallow_rank"]) == 1
        )
    ]

    if len(rows) != 1:
        raise RuntimeError(
            f"trial={trial} slot={slot}: "
            f"expected one winner, got {len(rows)}"
        )

    return rows[0]


def nearest_rational(x, max_den):
    f = Fraction(float(x)).limit_denominator(
        max_den
    )

    value = float(f)

    return (
        f.numerator,
        f.denominator,
        value,
        abs(float(x) - value),
    )


def schur_new_block(B):
    """
    B has rows:
        first 17 = original core
        remaining = already-added directions

    Returns the residualized Gram S of the added block modulo
    the original 17-dimensional span.
    """

    G = (
        B @ H
    ) @ B.T

    if B.shape[0] == 17:
        return G, None

    G00 = G[:17, :17]
    G01 = G[:17, 17:]
    G10 = G[17:, :17]
    G11 = G[17:, 17:]

    S = (
        G11
        -
        G10
        @ np.linalg.solve(
            G00,
            G01,
        )
    )

    return G, S


FIELDS = [
    "trial",
    "slot",
    "rank_before",

    "candidate_norm",
    "orthogonal_height",
    "orthogonal_ratio",

    "new_block_dim",
    "new_block_norm",
    "new_block_projection_share",

    "previous_coefficient",
    "previous_abs_coefficient",

    "previous_coeff_num",
    "previous_coeff_den",
    "previous_coeff_rational",
    "previous_coeff_rational_error",

    "previous_residual_corr",
    "previous_energy_share",
    "last_increment_corr",

    "new_block_max_energy_share",
    "new_block_last_energy_share",

    "alpha_l2",
    "alpha_max_abs",

    "alpha_coefficients",
    "alpha_rationals",
]


trial_ids = sorted(
    {
        int(row["trial"])
        for row in motifs
    }
)[:args.trials]

rows_out = []

t0 = perf_counter()

for pos, trial in enumerate(
    trial_ids,
    start=1,
):

    B0 = load_core(
        trial
    )

    B = B0.copy()

    selected = []

    for slot in range(1,5):

        winner = winner_for(
            trial,
            slot,
        )

        q = parse_vec(
            winner["vector"]
        )

        current_rank = B.shape[0]

        G = (
            B @ H
        ) @ B.T

        # Pair candidate against current basis.
        b = (
            B @ H @ q
        )

        coeff = np.linalg.solve(
            G,
            b,
        )

        projection_norm = float(
            b @ coeff
        )

        candidate_norm = float(
            q @ H @ q
        )

        orth = float(
            candidate_norm
            -
            projection_norm
        )

        if orth < 0 and orth > -1e-8:
            orth = 0.0

        orth_ratio = (
            orth / candidate_norm
            if candidate_norm > 0
            else float("nan")
        )

        if slot == 1:
            # No previous + directions yet.
            rows_out.append({
                "trial": trial,
                "slot": slot,
                "rank_before": current_rank,

                "candidate_norm": candidate_norm,
                "orthogonal_height": orth,
                "orthogonal_ratio": orth_ratio,

                "new_block_dim": 0,
                "new_block_norm": 0.0,
                "new_block_projection_share": 0.0,

                "previous_coefficient": float("nan"),
                "previous_abs_coefficient": float("nan"),

                "previous_coeff_num": "",
                "previous_coeff_den": "",
                "previous_coeff_rational": "",
                "previous_coeff_rational_error": "",

                "previous_residual_corr": float("nan"),
                "previous_energy_share": float("nan"),
                "last_increment_corr": float("nan"),

                "new_block_max_energy_share": float("nan"),
                "new_block_last_energy_share": float("nan"),

                "alpha_l2": 0.0,
                "alpha_max_abs": 0.0,

                "alpha_coefficients": "",
                "alpha_rationals": "",
            })

        else:
            G00 = G[:17, :17]
            G01 = G[:17, 17:]
            G10 = G[17:, :17]
            G11 = G[17:, 17:]

            b0 = b[:17]
            b1 = b[17:]

            S = (
                G11
                -
                G10
                @ np.linalg.solve(
                    G00,
                    G01,
                )
            )

            residual_pair = (
                b1
                -
                G10
                @ np.linalg.solve(
                    G00,
                    b0,
                )
            )

            alpha = np.linalg.solve(
                S,
                residual_pair,
            )

            new_block_norm = float(
                residual_pair @ alpha
            )

            new_block_projection_share = (
                new_block_norm / projection_norm
                if projection_norm > 1e-30
                else 0.0
            )

            # ------------------------------------------------
            # Correlation with immediately previous residualized
            # added generator.
            # ------------------------------------------------

            previous_index = len(alpha) - 1

            prev_norm = float(
                S[previous_index, previous_index]
            )

            prev_pair = float(
                residual_pair[previous_index]
            )

            previous_residual_corr = (
                abs(prev_pair)
                /
                math.sqrt(
                    max(
                        prev_norm * new_block_norm,
                        1e-30,
                    )
                )
            )

            previous_residual_corr = min(
                previous_residual_corr,
                1.0,
            )

            previous_coefficient = float(
                alpha[-1]
            )

            (
                rn,
                rd,
                rv,
                re,
            ) = nearest_rational(
                previous_coefficient,
                args.max_denominator,
            )

            # ------------------------------------------------
            # Energy distribution in sequential residual block.
            # ------------------------------------------------

            try:
                L = np.linalg.cholesky(
                    S
                )

                z = (
                    L.T
                    @ alpha
                )

                e = z * z
                esum = float(
                    e.sum()
                )

                if esum > 1e-30:
                    shares = e / esum

                    max_share = float(
                        shares.max()
                    )

                    last_share = float(
                        shares[-1]
                    )

                    last_increment_corr = float(
                        math.sqrt(
                            max(
                                0.0,
                                min(
                                    1.0,
                                    last_share,
                                )
                            )
                        )
                    )
                else:
                    max_share = 0.0
                    last_share = 0.0
                    last_increment_corr = 0.0

            except np.linalg.LinAlgError:
                max_share = float("nan")
                last_share = float("nan")
                last_increment_corr = float("nan")

            # ------------------------------------------------
            # Rational approximants for all chain coefficients.
            # ------------------------------------------------

            rational_parts = []

            for a in alpha:
                n, d, value, err = nearest_rational(
                    float(a),
                    args.max_denominator,
                )

                rational_parts.append(
                    f"{n}/{d}:{err:.6g}"
                )

            rows_out.append({
                "trial": trial,
                "slot": slot,
                "rank_before": current_rank,

                "candidate_norm": candidate_norm,
                "orthogonal_height": orth,
                "orthogonal_ratio": orth_ratio,

                "new_block_dim": len(alpha),
                "new_block_norm": new_block_norm,
                "new_block_projection_share":
                    new_block_projection_share,

                "previous_coefficient":
                    previous_coefficient,

                "previous_abs_coefficient":
                    abs(previous_coefficient),

                "previous_coeff_num": rn,
                "previous_coeff_den": rd,
                "previous_coeff_rational": rv,
                "previous_coeff_rational_error": re,

                "previous_residual_corr":
                    previous_residual_corr,

                "previous_energy_share":
                    last_share,

                "last_increment_corr":
                    last_increment_corr,

                "new_block_max_energy_share":
                    max_share,

                "new_block_last_energy_share":
                    last_share,

                "alpha_l2":
                    float(np.linalg.norm(alpha)),

                "alpha_max_abs":
                    float(np.max(np.abs(alpha))),

                "alpha_coefficients":
                    " ".join(
                        f"{x:.15g}"
                        for x in alpha
                    ),

                "alpha_rationals":
                    " ".join(
                        rational_parts
                    ),
            })

        selected.append(
            q.copy()
        )

        B = np.vstack(
            [
                B,
                q,
            ]
        )

    if (
        pos % 20 == 0
        or pos == len(trial_ids)
    ):
        print(
            f"PROGRESS"
            f"|trial={pos}/{len(trial_ids)}"
            f"|rows={len(rows_out)}"
            f"|seconds={perf_counter()-t0:.2f}",
            flush=True,
        )


out_path = (
    OUT
    / "growth_chain.tsv"
)

with out_path.open(
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=FIELDS,
    )

    writer.writeheader()
    writer.writerows(
        rows_out
    )


print()
print("=" * 72)
print("DONE")
print("rows =", len(rows_out))
print("saved =", out_path)
print("seconds =", perf_counter()-t0)
