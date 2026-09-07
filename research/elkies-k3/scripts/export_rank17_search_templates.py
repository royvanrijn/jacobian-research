from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]

J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
GDIR = BASE / "results/rank17-generator"

H = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

V = np.loadtxt(
    GDIR / "top100_vectors.txt",
    dtype=np.int64
)

if V.ndim == 1:
    V = V.reshape(1,-1)

OUT = GDIR / "search_templates.tsv"

with OUT.open("w") as f:
    f.write(
        "id\tnorm\tcoords\tpairings\n"
    )

    for i,v in enumerate(V):

        b = H @ v
        norm = int(v @ H @ v)

        f.write(
            f"{i}\t{norm}\t"
            + ",".join(map(str,v))
            + "\t"
            + ",".join(map(str,b))
            + "\n"
        )

        print(
            f"TEMPLATE|id={i}"
            f"|norm={norm}"
            f"|maxpair={np.max(np.abs(b))}"
            f"|coords={v.tolist()}"
            f"|pairings={b.tolist()}"
        )

print("saved", OUT)
