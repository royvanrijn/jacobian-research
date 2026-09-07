from __future__ import annotations

from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import math
import sqlite3

import numpy as np
from sage.all import ZZ, matrix


parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--trials", type=int, default=1000)
parser.add_argument("--core-pool", type=int, default=50000)
parser.add_argument("--candidate-pool", type=int, default=100000)
parser.add_argument("--controls", type=int, default=500)

parser.add_argument(
    "--selection-metric",
    choices=["height", "ratio"],
    default="height",
)

parser.add_argument("--seed", type=int, default=210021)
parser.add_argument("--run-name", default="rank21-longrun-v3")
parser.add_argument("--resume", action="store_true")
parser.add_argument("--progress-every", type=int, default=10)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / args.run_name
OUT.mkdir(parents=True, exist_ok=True)

DB_PATH = OUT / "experiment.sqlite"
MANIFEST_PATH = OUT / "manifest.json"

H = np.loadtxt(args.gram, dtype=np.float64)
VALL = np.load(args.vectors, mmap_mode="r")
NALL = np.load(args.norms, mmap_mode="r")

if H.shape != (21, 21):
    raise RuntimeError(f"Expected 21x21 Gram, got {H.shape}")

if VALL.ndim != 2 or VALL.shape[1] != 21:
    raise RuntimeError(f"Expected Nx21 vectors, got {VALL.shape}")

if len(VALL) != len(NALL):
    raise RuntimeError("vector/norm mismatch")

CORE_LIMIT = min(args.core_pool, len(VALL))
CAND_LIMIT = min(args.candidate_pool, len(VALL))

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

# ------------------------------------------------------------
# Main performance improvement.
#
# Every slot repeatedly needs
#
#     <candidate, B_i>
#
# so calculate candidate * H only once.
# ------------------------------------------------------------

print("precomputing candidate * H ...", flush=True)

CAND_H = CAND @ H

print("ready", flush=True)


def config_dict():
    return {
        "gram": str(Path(args.gram).resolve()),
        "vectors": str(Path(args.vectors).resolve()),
        "norms": str(Path(args.norms).resolve()),
        "shape": list(VALL.shape),
        "core_pool": CORE_LIMIT,
        "candidate_pool": CAND_LIMIT,
        "controls": args.controls,
        "selection_metric": args.selection_metric,
        "seed": args.seed,
        "version": 3,
    }


CONFIG = config_dict()

CONFIG_HASH = hashlib.sha256(
    json.dumps(
        CONFIG,
        sort_keys=True,
    ).encode()
).hexdigest()


# ============================================================
# Database
# ============================================================

conn = sqlite3.connect(DB_PATH)

conn.execute("PRAGMA journal_mode=WAL")
conn.execute("PRAGMA synchronous=NORMAL")
conn.execute("PRAGMA temp_store=MEMORY")

conn.executescript("""
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS cores (
    trial INTEGER PRIMARY KEY,
    trial_seed INTEGER NOT NULL,

    saturation_index INTEGER,
    logdet REAL,
    eig_min REAL,
    eig_max REAL,
    condition_number REAL,

    core_norm_min REAL,
    core_norm_mean REAL,
    core_norm_max REAL,

    selected_1 INTEGER,
    selected_2 INTEGER,
    selected_3 INTEGER,
    selected_4 INTEGER,

    completed_at REAL
);

CREATE TABLE IF NOT EXISTS samples (
    trial INTEGER NOT NULL,
    slot INTEGER NOT NULL,

    label INTEGER NOT NULL,
    pool_index INTEGER NOT NULL,

    rank_before INTEGER NOT NULL,

    norm REAL,
    projection_norm REAL,

    orthogonal_height REAL,
    orthogonal_ratio REAL,

    height_percentile REAL,
    ratio_percentile REAL,

    grid1_score REAL,
    grid2_score REAL,
    grid3_score REAL,
    grid4_score REAL,

    half_advantage REAL,
    third_advantage REAL,
    quarter_advantage REAL,

    coeff_l2 REAL,
    frac_l2 REAL,
    frac_max REAL,

    pair_abs_mean REAL,
    pair_abs_max REAL,
    pair_corr_mean REAL,
    pair_corr_max REAL,

    support INTEGER,
    max_coeff INTEGER,

    PRIMARY KEY (
        trial,
        slot,
        label,
        pool_index
    )
);

CREATE TABLE IF NOT EXISTS selected_vectors (
    trial INTEGER NOT NULL,
    slot INTEGER NOT NULL,
    pool_index INTEGER NOT NULL,

    vector TEXT NOT NULL,
    projection_coeffs TEXT NOT NULL,

    PRIMARY KEY (trial, slot)
);

CREATE INDEX IF NOT EXISTS samples_label_idx
    ON samples(label);

CREATE INDEX IF NOT EXISTS samples_slot_idx
    ON samples(slot);

CREATE INDEX IF NOT EXISTS samples_trial_idx
    ON samples(trial);
""")

stored_hash = conn.execute(
    "SELECT value FROM meta WHERE key='config_hash'"
).fetchone()

if stored_hash is not None:
    if not args.resume:
        raise RuntimeError(
            f"{DB_PATH} already exists.\n"
            "Use --resume or choose another --run-name."
        )

    if stored_hash[0] != CONFIG_HASH:
        raise RuntimeError(
            "Resume configuration differs from existing run.\n"
            f"existing={stored_hash[0]}\n"
            f"current ={CONFIG_HASH}"
        )
else:
    conn.execute(
        "INSERT INTO meta(key,value) VALUES (?,?)",
        ("config_hash", CONFIG_HASH),
    )

    conn.execute(
        "INSERT INTO meta(key,value) VALUES (?,?)",
        ("config", json.dumps(CONFIG, sort_keys=True)),
    )

    conn.commit()

MANIFEST_PATH.write_text(
    json.dumps(
        {
            **CONFIG,
            "config_hash": CONFIG_HASH,
            "database": str(DB_PATH),
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)


# ============================================================
# Helpers
# ============================================================

def exact_rank(B_np):
    return matrix(
        ZZ,
        B_np.tolist(),
    ).rank()


def saturation_index(B_np):
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

        for i in range(
            min(
                D.nrows(),
                D.ncols(),
            )
        ):
            x = abs(
                int(D[i, i])
            )

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


def normalized_grid_score(c, k):
    x = k * c
    delta = x - np.rint(x)

    return float(
        np.sqrt(
            np.mean(
                delta * delta
            )
        )
    )


def percentile(sorted_values, value):
    return float(
        np.searchsorted(
            sorted_values,
            value,
            side="right",
        )
        /
        len(sorted_values)
    )


def build_core(rng):
    """
    Generate an independent short-biased rank-17 core.

    Numerical rank is used while constructing it.
    Sage checks the final result exactly.
    """

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
                (u * u)
                *
                CORE_LIMIT
            )

            idx = min(
                idx,
                CORE_LIMIT - 1,
            )

            if idx in used:
                continue

            candidate = CORE_POOL[idx]

            if chosen:
                A = np.asarray(
                    chosen + [candidate],
                    dtype=np.float64,
                )

                if np.linalg.matrix_rank(A) != len(chosen) + 1:
                    continue

            chosen.append(
                candidate.copy()
            )

            used.add(idx)

        if len(chosen) != 17:
            continue

        B = np.asarray(
            chosen,
            dtype=np.int64,
        )

        # Exact check once, outside the inner loop.
        if exact_rank(B) == 17:
            return B


def slot_geometry(B):
    """
    Geometry relative to the CURRENT span.

    This is recomputed after each selected generator, so slot 2
    genuinely means rank18 -> rank19, etc.
    """

    G = (B @ H) @ B.T

    C = CAND_H @ B.T

    # Solve G * x = b rather than explicitly computing inverse.
    coeff = np.linalg.solve(
        G,
        C.T,
    ).T

    projection_norm = np.einsum(
        "ij,ij->i",
        C,
        coeff,
    )

    residual = (
        CAND_NORM
        -
        projection_norm
    )

    residual[
        (residual < 0)
        &
        (residual > -1e-7)
    ] = 0.0

    ratios = np.full(
        len(residual),
        np.nan,
        dtype=np.float64,
    )

    good_norm = CAND_NORM > 0

    ratios[good_norm] = (
        residual[good_norm]
        /
        CAND_NORM[good_norm]
    )

    viable = np.flatnonzero(
        (residual > 1e-7)
        &
        np.isfinite(ratios)
    )

    return (
        G,
        C,
        coeff,
        projection_norm,
        residual,
        ratios,
        viable,
    )


def feature_row(
    trial,
    slot,
    label,
    idx,
    B,
    G,
    C,
    coeff,
    projection_norm,
    residual,
    ratios,
    height_sorted,
    ratio_sorted,
):
    c = coeff[idx]

    n = float(
        CAND_NORM[idx]
    )

    proj = float(
        projection_norm[idx]
    )

    orth = float(
        residual[idx]
    )

    ratio = float(
        ratios[idx]
    )

    g1 = normalized_grid_score(c, 1)
    g2 = normalized_grid_score(c, 2)
    g3 = normalized_grid_score(c, 3)
    g4 = normalized_grid_score(c, 4)

    frac = c - np.rint(c)

    pair_abs = np.abs(
        C[idx]
    )

    denom = np.sqrt(
        np.maximum(
            n * np.diag(G),
            1e-30,
        )
    )

    corr = pair_abs / denom

    v = CAND[idx]

    return (
        trial,
        slot,
        label,
        int(idx),

        B.shape[0],

        n,
        proj,

        orth,
        ratio,

        percentile(
            height_sorted,
            orth,
        ),

        percentile(
            ratio_sorted,
            ratio,
        ),

        g1,
        g2,
        g3,
        g4,

        g2 - g1,
        g3 - g1,
        g4 - min(g1, g2),

        float(np.linalg.norm(c)),
        float(np.linalg.norm(frac)),
        float(np.max(np.abs(frac))),

        float(np.mean(pair_abs)),
        float(np.max(pair_abs)),

        float(np.mean(corr)),
        float(np.max(corr)),

        int(np.count_nonzero(v)),
        int(np.max(np.abs(v))),
    )


SAMPLE_INSERT = """
INSERT OR REPLACE INTO samples VALUES (
    ?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?,?
)
"""


# ============================================================
# Resume state
# ============================================================

completed = {
    int(row[0])
    for row in conn.execute(
        "SELECT trial FROM cores"
    )
}

print()
print("RANK21 17+4 LONG RUN")
print("run              =", args.run_name)
print("database         =", DB_PATH)
print("full pool        =", len(VALL))
print("core pool        =", CORE_LIMIT)
print("candidate pool   =", CAND_LIMIT)
print("controls/slot    =", args.controls)
print("selection metric =", args.selection_metric)
print("requested trials =", args.trials)
print("already complete =", len(completed))
print("config hash      =", CONFIG_HASH[:16])
print()

global_start = perf_counter()


# ============================================================
# Main experiment
# ============================================================

for trial in range(args.trials):

    if trial in completed:
        continue

    trial_start = perf_counter()

    # Deterministic per-trial seed means resume is reproducible.
    trial_seed = (
        args.seed
        +
        trial * 1_000_003
    )

    rng = np.random.default_rng(
        trial_seed
    )

    B = build_core(
        rng
    )

    G17 = (B @ H) @ B.T

    eig = np.linalg.eigvalsh(
        G17
    )

    sign, logdet = np.linalg.slogdet(
        G17
    )

    if sign <= 0:
        raise RuntimeError(
            "non-positive rank17 core"
        )

    sat = saturation_index(
        B
    )

    core_norms = np.einsum(
        "ij,ij->i",
        B @ H,
        B,
    )

    selected_indices = []

    sample_rows = []
    selected_rows = []

    # --------------------------------------------------------
    # TRUE sequential:
    #
    # rank17 -> rank18
    # rank18 -> rank19
    # rank19 -> rank20
    # rank20 -> rank21
    # --------------------------------------------------------

    for slot in range(1, 5):

        (
            G,
            C,
            coeff,
            projection_norm,
            residual,
            ratios,
            viable,
        ) = slot_geometry(B)

        if not len(viable):
            raise RuntimeError(
                f"trial {trial}: no viable extensions at slot {slot}"
            )

        height_sorted = np.sort(
            residual[viable]
        )

        ratio_sorted = np.sort(
            ratios[viable]
        )

        if args.selection_metric == "height":
            selected_idx = int(
                viable[
                    np.argmin(
                        residual[viable]
                    )
                ]
            )
        else:
            selected_idx = int(
                viable[
                    np.argmin(
                        ratios[viable]
                    )
                ]
            )

        # -----------------------------
        # Positive
        # -----------------------------

        sample_rows.append(
            feature_row(
                trial,
                slot,
                1,
                selected_idx,
                B,
                G,
                C,
                coeff,
                projection_norm,
                residual,
                ratios,
                height_sorted,
                ratio_sorted,
            )
        )

        selected_rows.append(
            (
                trial,
                slot,
                selected_idx,

                " ".join(
                    map(
                        str,
                        map(
                            int,
                            CAND[selected_idx],
                        ),
                    )
                ),

                " ".join(
                    f"{x:.15g}"
                    for x in coeff[selected_idx]
                ),
            )
        )

        # -----------------------------
        # Negative controls
        # -----------------------------

        controls_pool = viable[
            viable != selected_idx
        ]

        ncontrols = min(
            args.controls,
            len(controls_pool),
        )

        if ncontrols:
            control_ids = rng.choice(
                controls_pool,
                size=ncontrols,
                replace=False,
            )

            for idx in control_ids:
                sample_rows.append(
                    feature_row(
                        trial,
                        slot,
                        0,
                        int(idx),
                        B,
                        G,
                        C,
                        coeff,
                        projection_norm,
                        residual,
                        ratios,
                        height_sorted,
                        ratio_sorted,
                    )
                )

        selected_indices.append(
            selected_idx
        )

        # Add selected direction.
        #
        # residual > 0 already certifies real independence in
        # numerical height geometry. We exact-check final rank21.
        B = np.vstack(
            [
                B,
                CAND[selected_idx],
            ]
        )

    if exact_rank(B) != 21:
        raise RuntimeError(
            f"trial {trial}: final exact rank != 21"
        )

    # --------------------------------------------------------
    # Atomic transaction:
    #
    # either entire trial is persisted or none of it is.
    # --------------------------------------------------------

    conn.execute("BEGIN")

    try:
        conn.executemany(
            SAMPLE_INSERT,
            sample_rows,
        )

        conn.executemany(
            """
            INSERT OR REPLACE INTO selected_vectors
            VALUES (?,?,?,?,?)
            """,
            selected_rows,
        )

        conn.execute(
            """
            INSERT OR REPLACE INTO cores VALUES (
                ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
            )
            """,
            (
                trial,
                int(trial_seed),

                int(sat),
                float(logdet),
                float(eig[0]),
                float(eig[-1]),
                float(eig[-1] / eig[0]),

                float(np.min(core_norms)),
                float(np.mean(core_norms)),
                float(np.max(core_norms)),

                int(selected_indices[0]),
                int(selected_indices[1]),
                int(selected_indices[2]),
                int(selected_indices[3]),

                float(perf_counter() - trial_start),
            ),
        )

        conn.commit()

    except Exception:
        conn.rollback()
        raise

    if (
        trial % args.progress_every == 0
        or trial + 1 == args.trials
    ):
        positives = [
            row
            for row in sample_rows
            if row[2] == 1
        ]

        fields = []

        for row in positives:
            slot = row[1]
            orth = row[7]
            ratio = row[8]
            hpct = 100 * row[9]
            rpct = 100 * row[10]

            fields.append(
                f"s{slot}:"
                f"orth={orth:.5g},"
                f"ratio={100*ratio:.4g}%,"
                f"hPct={hpct:.4g}%,"
                f"rPct={rpct:.4g}%"
            )

        print(
            f"PROGRESS"
            f"|trial={trial}"
            f"|sat={sat}"
            f"|{' | '.join(fields)}"
            f"|trial_s={perf_counter()-trial_start:.3f}"
            f"|total_s={perf_counter()-global_start:.1f}",
            flush=True,
        )


print()
print("=" * 76)
print("DONE")

ncores = conn.execute(
    "SELECT COUNT(*) FROM cores"
).fetchone()[0]

nsamples = conn.execute(
    "SELECT COUNT(*) FROM samples"
).fetchone()[0]

npos = conn.execute(
    "SELECT COUNT(*) FROM samples WHERE label=1"
).fetchone()[0]

nneg = conn.execute(
    "SELECT COUNT(*) FROM samples WHERE label=0"
).fetchone()[0]

print("cores    =", ncores)
print("samples  =", nsamples)
print("positive =", npos)
print("controls =", nneg)
print("database =", DB_PATH)

conn.close()
