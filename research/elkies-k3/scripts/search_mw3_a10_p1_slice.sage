from sage.all import *
from pathlib import Path
import argparse
import itertools
import numpy as np


ap = argparse.ArgumentParser(
    description="Exhaust an A10/MW3 zero-dimensional coordinate slice over a small prime field."
)
ap.add_argument("--input", required=True)
ap.add_argument("--open-input", default=None)
ap.add_argument("--nonzero", default="s1,y4")
ap.add_argument("--max-hits", type=int, default=100)
args = ap.parse_args()

lines = [line.strip() for line in Path(args.input).read_text().splitlines() if line.strip()]
names = [name.strip() for name in lines[0].split(",") if name.strip()]
p = int(lines[1])
if p > 101:
    raise SystemExit("this exhaustive scanner is intended for small prime fields")
R = PolynomialRing(GF(p), names, order="degrevlex")
equations = [
    R((line[:-1] if line.endswith(",") else line).replace("^", "**"))
    for line in lines[2:]
]

# Saturation variables introduced for msolve are not part of the geometric
# coordinate scan.  Drop their inverse equation and reconstruct only the chart
# variables.
geometry_names = [name for name in names if name != "sat"]
if "sat" in names:
    sat_index = names.index("sat")
    equations = [e for e in equations if e.degree(R.gen(sat_index)) == 0]
    S = PolynomialRing(GF(p), geometry_names, order="degrevlex")
    sd = S.gens_dict()
    phi = R.hom([S(0) if name == "sat" else sd[name] for name in names], S)
    equations = [phi(e) for e in equations]
    R = S
    names = geometry_names

nonzero = {name.strip() for name in args.nonzero.split(",") if name.strip()}
if not nonzero.issubset(names):
    raise RuntimeError(f"unknown nonzero coordinate: {sorted(nonzero-set(names))}")

domains = [range(1, p) if name in nonzero else range(p) for name in names]
count = prod(len(domain) for domain in domains)
if count > 5_000_000:
    raise SystemExit(f"refusing exhaustive grid of {count} points")

grid = np.array(list(itertools.product(*domains)), dtype=np.int64).T
print(
    f"MW3A10SCAN|p={p}|vars={','.join(names)}|points={grid.shape[1]}"
    f"|eqs={len(equations)}|nonzero={','.join(sorted(nonzero))}",
    flush=True,
)


def evaluate(poly, points):
    values = np.zeros(points.shape[1], dtype=np.int64)
    power_cache = {}
    for exponents, coefficient in poly.dict().items():
        term = np.full(points.shape[1], int(coefficient), dtype=np.int64)
        for i, exponent in enumerate(exponents):
            if exponent == 0:
                continue
            key = (i, exponent)
            if key not in power_cache:
                power_cache[key] = np.power(points[i], exponent, dtype=np.int64) % p
            term = (term * power_cache[key]) % p
        values += term
        values %= p
    return values


if args.open_input:
    open_lines = [
        line.strip()
        for line in Path(args.open_input).read_text().splitlines()
        if line.strip()
    ]
    if open_lines[0].split(",") != names or int(open_lines[1]) != p:
        raise RuntimeError("open-condition ring does not match scan ring")
    open_polynomials = [
        R((line[:-1] if line.endswith(",") else line).replace("^", "**"))
        for line in open_lines[2:]
    ]
    for i, open_polynomial in enumerate(open_polynomials):
        open_values = evaluate(open_polynomial, grid)
        grid = grid[:, open_values != 0]
        print(
            f"MW3A10SCAN_OPEN|i={i}|degree={open_polynomial.total_degree()}"
            f"|terms={len(open_polynomial.monomials())}|survivors={grid.shape[1]}",
            flush=True,
        )


# Sparse-first filtering minimizes total work and is exact over GF(p).
order = sorted(range(len(equations)), key=lambda i: len(equations[i].monomials()))
survivors = grid
for step, equation_index in enumerate(order, 1):
    equation = equations[equation_index]
    values = evaluate(equation, survivors)
    survivors = survivors[:, values == 0]
    print(
        f"MW3A10SCAN_EQ|step={step}|source={equation_index}"
        f"|degree={equation.total_degree()}|terms={len(equation.monomials())}"
        f"|survivors={survivors.shape[1]}",
        flush=True,
    )
    if survivors.shape[1] == 0:
        break

for column in range(min(survivors.shape[1], args.max_hits)):
    print(
        "MW3A10SCAN_HIT|" + ",".join(
            f"{name}={int(survivors[i, column])}" for i, name in enumerate(names)
        ),
        flush=True,
    )
print(f"MW3A10SCAN|hits={survivors.shape[1]}", flush=True)
