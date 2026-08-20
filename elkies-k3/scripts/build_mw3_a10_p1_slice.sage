from sage.all import *
from pathlib import Path
import argparse
import random


ap = argparse.ArgumentParser(
    description="Build a reproducible coordinate slice of the normalized A10/MW3 P1 chart."
)
ap.add_argument("--input", required=True)
ap.add_argument("--open-input", default=None)
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--kill", default="rho,r1,lam")
ap.add_argument(
    "--saturate",
    default="s1",
    help="comma-separated remaining variables whose product is inverted",
)
args = ap.parse_args()

lines = [line.strip() for line in Path(args.input).read_text().splitlines() if line.strip()]
names = [name.strip() for name in lines[0].split(",") if name.strip()]
p = int(lines[1])
K = GF(p)
R = PolynomialRing(K, names, order="degrevlex")
d = R.gens_dict()
equations = [
    R((line[:-1] if line.endswith(",") else line).replace("^", "**"))
    for line in lines[2:]
]

kill = [name.strip() for name in args.kill.split(",") if name.strip()]
if len(kill) != len(names) - len(equations):
    raise RuntimeError(
        f"expected {len(names)-len(equations)} sliced variables, got {len(kill)}"
    )
if any(name not in d for name in kill):
    raise RuntimeError(f"unknown sliced variable in {kill}")

rng = random.Random(args.seed)
values = {name: K(rng.randrange(p)) for name in kill}
if "rho" in values and values["rho"] == 0:
    values["rho"] = K(1 + rng.randrange(p - 1))
if "lam" in values:
    while values["lam"] in (0, 1):
        values["lam"] = K(rng.randrange(p))

keep = [name for name in names if name not in kill]
sat_names = [name.strip() for name in args.saturate.split(",") if name.strip()]
if any(name not in keep for name in sat_names):
    raise RuntimeError(f"saturation variable is not free after slicing: {sat_names}")
output_names = keep + (["sat"] if sat_names else [])
S = PolynomialRing(K, output_names, order="degrevlex")
sd = S.gens_dict()
phi = R.hom([S(values[name]) if name in values else sd[name] for name in names], S)
sliced = [phi(equation) for equation in equations]
if any(equation == 0 for equation in sliced):
    raise RuntimeError("slice made a residual equation identically zero")
if sat_names:
    sat_product = prod(sd[name] for name in sat_names)
    sliced.append(sd["sat"] * sat_product - 1)

open_sliced = None
if args.open_input:
    open_lines = [
        line.strip()
        for line in Path(args.open_input).read_text().splitlines()
        if line.strip()
    ]
    if open_lines[0].split(",") != names or int(open_lines[1]) != p:
        raise RuntimeError("open-condition ring does not match input ring")
    open_sliced = [
        phi(R((line[:-1] if line.endswith(",") else line).replace("^", "**")))
        for line in open_lines[2:]
    ]

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as handle:
    handle.write(",".join(output_names) + "\n")
    handle.write(str(p) + "\n")
    for i, equation in enumerate(sliced):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if i + 1 < len(sliced) else "\n")

meta = out.with_suffix(".meta.txt")
with meta.open("w") as handle:
    handle.write(f"seed={args.seed}\n")
    handle.write("kill=" + repr(kill) + "\n")
    handle.write("values=" + repr({name: int(value) for name, value in values.items()}) + "\n")
    handle.write("remaining=" + repr(keep) + "\n")
    handle.write("saturate=" + repr(sat_names) + "\n")
if open_sliced is not None:
    open_out = out.with_suffix(".open.ms")
    with open_out.open("w") as handle:
        handle.write(",".join(output_names) + "\n")
        handle.write(str(p) + "\n")
        for i, polynomial in enumerate(open_sliced):
            handle.write(str(polynomial).replace("**", "^"))
            handle.write(",\n" if i + 1 < len(open_sliced) else "\n")

print(
    f"MW3A10SLICE|p={p}|seed={args.seed}|kill={','.join(kill)}"
    "|values=" + ",".join(f"{name}={int(values[name])}" for name in kill),
    flush=True,
)
if open_sliced is not None:
    print(f"MW3A10SLICE|open_out={open_out}", flush=True)
for i, equation in enumerate(sliced):
    print(
        f"MW3A10SLICE_EQ|i={i}|degree={equation.total_degree()}"
        f"|terms={len(equation.monomials())}",
        flush=True,
    )
print(
    f"MW3A10SLICE|vars={len(output_names)}|eqs={len(sliced)}"
    f"|saturate={','.join(sat_names)}|out={out}|meta={meta}",
    flush=True,
)
