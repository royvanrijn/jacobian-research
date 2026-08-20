from pathlib import Path
import re
import numpy as np

from certify_nagao_rank21_t6793 import (
    load_pinned_input,
    PARAMETER_T,
)
from nagao_1994 import RANK21_CONSTRUCTION, short_jacobian_coefficients
from search_extra_points import gp_rational, gp_vector, run_gp

ROOT = Path(__file__).resolve().parents[2]

src = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_rank21_unbiased.json"
)

outdir = ROOT / "elkies-k3/data/rank19-nagao-t6793"
outdir.mkdir(parents=True, exist_ok=True)

points, provenance = load_pinned_input(src)

assert len(points) == 19

coeffs = short_jacobian_coefficients(
    RANK21_CONSTRUCTION,
    PARAMETER_T,
)

curve = ",".join(
    gp_rational(x)
    for x in coeffs
)

point_vector = ",".join(
    gp_vector(P)
    for P in points
)

program = "\n".join([
    "default(realprecision,120);",
    f"E=ellinit([{curve}]);",
    f"P=[{point_vector}];",
    "H=ellheightmatrix(E,P);",
    'print("HEIGHT_BEGIN");',
    "for(i=1,matsize(H)[1],print(Vec(H[i,])));",
    'print("HEIGHT_END");',
    "quit",
]) + "\n"

output, wall = run_gp(
    program,
    timeout=120.0,
    stack_bytes=1_000_000_000,
)

lines = [
    x.strip()
    for x in output.splitlines()
]

start = lines.index("HEIGHT_BEGIN") + 1
end = lines.index("HEIGHT_END")

rows = []

for line in lines[start:end]:
    line = line.strip()

    if not line:
        continue

    assert line.startswith("[") and line.endswith("]"), line

    values = line[1:-1].split(",")

    rows.append([
        float(x.strip())
        for x in values
    ])

H = np.asarray(rows, dtype=float)

assert H.shape == (19, 19), H.shape

# Numerical sanity.
sym_err = np.max(np.abs(H - H.T))
eig = np.linalg.eigvalsh(H)

print("shape =", H.shape)
print("symmetry error =", sym_err)
print("eig min =", eig[0])
print("eig max =", eig[-1])
print("condition =", eig[-1] / eig[0])
print("wall =", wall)

if sym_err > 1e-8:
    raise RuntimeError("height Gram is unexpectedly asymmetric")

if eig[0] <= 0:
    raise RuntimeError(
        "height Gram is not numerically positive definite"
    )

np.savetxt(
    outdir / "height-gram.txt",
    H,
    fmt="%.17g",
)

with (outdir / "points.txt").open("w") as f:
    for x, y in points:
        f.write(f"{x} {y}\n")

print("saved", outdir / "height-gram.txt")
