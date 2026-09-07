from pathlib import Path
import argparse
import sqlite3

import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("database")
args = parser.parse_args()

db = Path(args.database)

conn = sqlite3.connect(
    f"file:{db}?mode=ro",
    uri=True,
)

cores = pd.read_sql_query(
    "SELECT * FROM cores ORDER BY trial",
    conn,
)

p = pd.read_sql_query(
    """
    SELECT *
    FROM samples
    WHERE label=1
    ORDER BY trial, slot
    """,
    conn,
)

print("cores =", len(cores))
print("positive rows =", len(p))

if not len(p):
    raise SystemExit

print()
print("SEQUENTIAL EXTENSION STATISTICS")

for slot in range(1,5):

    x = p[p.slot == slot]

    print()
    print("slot", slot, f"(rank {16+slot} -> {17+slot})")

    r = 100 * x.orthogonal_ratio
    hp = 100 * x.height_percentile
    rp = 100 * x.ratio_percentile

    print(
        "  orth ratio % :",
        f"q25={r.quantile(.25):.6g}",
        f"median={r.median():.6g}",
        f"q75={r.quantile(.75):.6g}",
    )

    print(
        "  height pct   :",
        f"q25={hp.quantile(.25):.6g}",
        f"median={hp.median():.6g}",
        f"q75={hp.quantile(.75):.6g}",
    )

    print(
        "  ratio pct    :",
        f"q25={rp.quantile(.25):.6g}",
        f"median={rp.median():.6g}",
        f"q75={rp.quantile(.75):.6g}",
    )

print()
print("SATURATION")

print(
    cores.saturation_index
    .value_counts()
    .sort_index()
)

print()
print("GRID WINNERS")

cols = [
    "grid1_score",
    "grid2_score",
    "grid3_score",
    "grid4_score",
]

winner = (
    p[cols]
    .to_numpy()
    .argmin(axis=1)
    + 1
)

for k in range(1,5):
    print(
        k,
        int((winner == k).sum())
    )

print()
print("CORRELATIONS WITH ORTHOGONAL RATIO")

features = [
    "support",
    "max_coeff",
    "pair_corr_mean",
    "pair_corr_max",
    "grid1_score",
    "grid2_score",
    "grid3_score",
    "grid4_score",
]

for f in features:
    print(
        f"{f:20s}",
        f"{p.orthogonal_ratio.corr(p[f]):+.6f}"
    )

conn.close()
