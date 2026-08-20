from pathlib import Path
import numpy as np

from nagao_1994 import (
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_POINTS,
)
from search_extra_points import gp_rational, gp_vector, run_gp

ROOT = Path(__file__).resolve().parents[2]

OUT = ROOT / "elkies-k3/data/rank21-nagao-published"
OUT.mkdir(parents=True, exist_ok=True)

coeffs = RANK21_PUBLISHED_MODEL
points = RANK21_PUBLISHED_POINTS

assert len(points) == 21

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
    'print("ONCURVE ",vecsum(vector(#P,i,ellisoncurve(E,P[i]))));',
    "H=ellheightmatrix(E,P);",
    'print("RANK ",matrank(H));',
    'print("HEIGHT_BEGIN");',
    "for(i=1,matsize(H)[1],print(Vec(H[i,])));",
    'print("HEIGHT_END");',
    "quit",
]) + "\n"

output, wall = run_gp(
    program,
    timeout=180.0,
    stack_bytes=1_500_000_000,
)

lines = [
    line.strip()
    for line in output.splitlines()
    if line.strip()
]

oncurve = next(
    int(x.split()[1])
    for x in lines
    if x.startswith("ONCURVE ")
)

rank = next(
    int(x.split()[1])
    for x in lines
    if x.startswith("RANK ")
)

assert oncurve == 21, oncurve
assert rank == 21, rank

start = lines.index("HEIGHT_BEGIN") + 1
end = lines.index("HEIGHT_END")

rows = []

for line in lines[start:end]:
    assert line.startswith("[") and line.endswith("]"), line

    rows.append([
        float(x.strip())
        for x in line[1:-1].split(",")
    ])

H = np.asarray(rows, dtype=float)

assert H.shape == (21, 21)

sym_err = np.max(np.abs(H - H.T))
eig = np.linalg.eigvalsh(H)

print("points =", len(points))
print("numerical rank =", rank)
print("shape =", H.shape)
print("symmetry error =", sym_err)
print("eig min =", eig[0])
print("eig max =", eig[-1])
print("condition =", eig[-1] / eig[0])
print("diag min =", np.diag(H).min())
print("diag median =", np.median(np.diag(H)))
print("diag max =", np.diag(H).max())
print("wall =", wall)

assert sym_err < 1e-8
assert eig[0] > 0

np.savetxt(
    OUT / "height-gram.txt",
    H,
    fmt="%.17g",
)

with (OUT / "points.txt").open("w") as f:
    for x, y in points:
        f.write(f"{x} {y}\n")

print("saved", OUT / "height-gram.txt")
