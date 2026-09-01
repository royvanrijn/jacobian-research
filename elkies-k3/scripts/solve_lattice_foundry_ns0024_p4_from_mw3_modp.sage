#!/usr/bin/env sage -python
"""Solve the fourth minimum-basis section on a fixed oriented MW3 seed."""

import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ


parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=Path, required=True)
parser.add_argument(
    "--seed-index",
    type=int,
    default=0,
    help="zero-based record index when the seed file contains multiple MW3SEED lines",
)
parser.add_argument("--export-msolve", type=Path, required=True)
parser.add_argument("--c-minpoly", help="optional monic quadratic coefficients u,v for c^2+u*c+v")
parser.add_argument("--c-value", type=int, help="optional base-field value of c")
parser.add_argument(
    "--recursive-y",
    action="store_true",
    help=(
        "eliminate y2,...,y8 recursively; the default retains them to keep "
        "the fixed-fibre system sparse and cubic"
    ),
)
parser.add_argument(
    "--sparse-y",
    action="store_true",
    help="retain y1,...,y7; useful as a low-degree comparison system",
)
args = parser.parse_args()
if args.recursive_y and args.sparse_y:
    raise SystemExit("choose at most one of --recursive-y and --sparse-y")
if args.c_minpoly and args.c_value is not None:
    raise SystemExit("choose at most one of --c-minpoly and --c-value")

seed_records = [line for line in args.seed.resolve().read_text().splitlines() if line.strip()]
if args.seed_index < 0 or args.seed_index >= len(seed_records):
    raise SystemExit("--seed-index is outside the seed file")
fields = {}
for item in seed_records[args.seed_index].strip().split("|")[1:]:
    key, value = item.split("=", 1)
    fields[key] = value
prime = ZZ(fields["p"])
field = GF(prime)


def values(key):
    return [field(int(value)) for value in fields[key].split(",")]


names = ["c", "c0_inverse", "c1_inverse"] + [f"x{index}" for index in range(1, 5)]
if args.recursive_y:
    names += ["y1", "y1_inverse"]
elif args.sparse_y:
    names += [f"y{index}" for index in range(1, 8)]
    names += ["y0_tangent_inverse", "y1_tangent_inverse", "yi_tangent_inverse"]
else:
    names += ["y1", "y1_inverse", "y8", "y8_inverse"]
ring = PolynomialRing(field, names, order="degrevlex")
v = ring.gens_dict()
polynomials = PolynomialRing(ring, "t")
t = polynomials.gen()
A = polynomials(values("A"))
B = polynomials(values("B"))
r1, ri = field(int(fields["r1"])), field(int(fields["ri"]))
c = v["c"]
h = t - c
xs = [v[f"x{index}"] for index in range(1, 5)]
x5 = r1 * (1 - c)**2 - c**2 - ri - sum(xs)
X = c**2 + sum(xs[index - 1] * t**index for index in range(1, 5)) + x5 * t**5 + ri * t**6
rhs = X**3 + A * X * h**4 + B * h**6
if args.recursive_y:
    y = [ring.zero(), v["y1"]]
    for degree in range(3, 10):
        known = sum(
            y[left] * y[degree - left]
            for left in range(2, degree - 1)
        )
        y.append((rhs[degree] - known) * v["y1_inverse"] / 2)
    Y = sum(y[index] * t**index for index in range(1, 9))
    equations = [
        rhs[0], rhs[1], y[1]**2 - rhs[2],
        y[1] * v["y1_inverse"] - 1,
        sum(y[1:]),
    ]
    equations += [(Y**2 - rhs)[degree] for degree in range(10, 19)]
elif args.sparse_y:
    y = [ring.zero()] + [v[f"y{index}"] for index in range(1, 8)]
    y.append(-sum(y[1:]))
    Y = sum(y[index] * t**index for index in range(1, 9))
    equations = list(Y**2 - rhs)
    # Exact component one, rather than the closure of the raw node chart, has
    # nonzero first exceptional tangent at all three multiplicative fibres.
    equations += [
        y[1] * v["y0_tangent_inverse"] - 1,
        sum(index * y[index] for index in range(1, 9))
        * v["y1_tangent_inverse"] - 1,
        y[8] * v["yi_tangent_inverse"] - 1,
    ]
else:
    # Resolve from both multiplicative ends.  Depth one at I7 gives the
    # leading coefficient y1; depth one at I4 gives y8 after reversing the
    # base.  Three triangular steps from each end keep substitution shallow
    # while eliminating six middle coefficients.
    y = {1: v["y1"], 8: v["y8"]}
    for index in range(2, 5):
        degree = index + 1
        known = sum(
            y[left] * y[degree - left]
            for left in range(2, index)
        )
        y[index] = (rhs[degree] - known) * v["y1_inverse"] / 2
    for index in range(7, 4, -1):
        degree = index + 8
        known = sum(
            y[left] * y[degree - left]
            for left in range(index + 1, 9)
            if degree - left in y and degree - left != index
        )
        y[index] = (rhs[degree] - known) * v["y8_inverse"] / 2
    Y = sum(y[index] * t**index for index in range(1, 9))
    residual = Y**2 - rhs
    solved_degrees = set(range(2, 6)) | set(range(13, 17))
    equations = [
        v["y1"]**2 - rhs[2],
        v["y1"] * v["y1_inverse"] - 1,
        v["y8"]**2 - rhs[16],
        v["y8"] * v["y8_inverse"] - 1,
        Y(t=1),
    ]
    equations += [
        residual[degree]
        for degree in range(19)
        if degree not in solved_degrees
    ]
equations += [c * v["c0_inverse"] - 1, (c - 1) * v["c1_inverse"] - 1]
equations = [entry for entry in equations if entry != 0]
if args.c_minpoly:
    u, constant = (field(int(value)) for value in args.c_minpoly.split(","))
    equations.append(c**2 + u*c + constant)
if args.c_value is not None:
    equations.append(c - field(args.c_value))

output = args.export_msolve.resolve()
output.write_text(
    ",".join(names) + "\n" + str(prime) + "\n"
    + ",\n".join(str(equation) for equation in equations) + "\n"
)
print(
    f"NS0024P4SYSTEM|p={prime}|variables={ring.ngens()}|equations={len(equations)}"
    f"|y_mode={'recursive' if args.recursive_y else 'sparse' if args.sparse_y else 'bidirectional'}"
    f"|terms={sum(len(item.dict()) for item in equations)}"
    f"|output={output}",
    flush=True,
)
