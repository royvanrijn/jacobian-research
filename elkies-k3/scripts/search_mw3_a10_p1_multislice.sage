from sage.all import *
from pathlib import Path
import argparse
import itertools
import random
import time
import numpy as np


ap = argparse.ArgumentParser(
    description=(
        "Exhaust many normalized A10/P1 coordinate slices in one Sage process. "
        "The output is compatible with run_mw3_a10_p1p2_batch.py."
    )
)
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument(
    "--input",
    default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.ms",
)
ap.add_argument(
    "--open-input",
    default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.open.ms",
)
ap.add_argument("--dir", default="artifacts/local/elkies-k3/mw3-a10-p1")
ap.add_argument("--max-hits", type=int, default=20)
ap.add_argument("--fixed-names", default="rho,r1,lam")
ap.add_argument("--nonzero-keep", default="s1,y4")
ap.add_argument("--prefix", default="component2-valid")
ap.add_argument(
    "--deterministic-fixed",
    action="store_true",
    help="enumerate (rho,r1,lambda) lexicographically instead of seeded random slices",
)
args = ap.parse_args()

if args.seed_end < args.seed_start:
    raise SystemExit("--seed-end must be at least --seed-start")


def read_system(path):
    lines = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    return lines[0].split(","), int(lines[1]), lines[2:]


names, p, equation_strings = read_system(args.input)
open_names, open_p, open_strings = read_system(args.open_input)
if names != open_names or p != open_p:
    raise RuntimeError("open-condition ring does not match input ring")
if p > 101:
    raise SystemExit("this exhaustive scanner is intended for small prime fields")

K = GF(p)
R = PolynomialRing(K, names, order="degrevlex")
equations = [
    R((line[:-1] if line.endswith(",") else line).replace("^", "**"))
    for line in equation_strings
]
open_polynomials = [
    R((line[:-1] if line.endswith(",") else line).replace("^", "**"))
    for line in open_strings
]

fixed_names = [name for name in args.fixed_names.split(",") if name]
keep_names = [name for name in names if name not in fixed_names]
nonzero_keep = {name for name in args.nonzero_keep.split(",") if name}
if not fixed_names or len(keep_names) != 4:
    raise RuntimeError("scanner requires exactly four enumerated variables")
if set(fixed_names + keep_names) != set(names) or not nonzero_keep <= set(keep_names):
    raise RuntimeError(f"invalid fixed/nonzero coordinates for parent variables: {names}")
fixed_indices = [names.index(name) for name in fixed_names]
keep_indices = [names.index(name) for name in keep_names]


def term_records(poly):
    records = []
    for exponents, coefficient in poly.dict().items():
        exponents = tuple(exponents)
        records.append(
            (
                tuple(exponents[i] for i in keep_indices),
                tuple(exponents[i] for i in fixed_indices),
                int(coefficient),
            )
        )
    return records


equation_records = [term_records(poly) for poly in equations]
open_records = [term_records(poly) for poly in open_polynomials]

# s1 and y4 are open coordinates, so omit zero from their domains.  This is
# the same 864,900-point grid used by search_mw3_a10_p1_slice.sage, built once.
domains = [range(1, p) if name in nonzero_keep else range(p) for name in keep_names]
grid = np.array(list(itertools.product(*domains)), dtype=np.uint8).T
base_indices = np.arange(grid.shape[1], dtype=np.int64)

max_powers = [0] * len(keep_names)
for records in equation_records + open_records:
    for keep_exponents, _, _ in records:
        for i, exponent in enumerate(keep_exponents):
            max_powers[i] = max(max_powers[i], exponent)

# uint8 is sufficient for field residues.  Products are promoted to uint16
# before reduction, so no intermediate exceeds 30^2.
powers = []
for variable_index, maximum in enumerate(max_powers):
    variable_powers = [np.ones(grid.shape[1], dtype=np.uint8)]
    for _ in range(maximum):
        variable_powers.append(
            ((variable_powers[-1].astype(np.uint16) * grid[variable_index]) % p).astype(
                np.uint8
            )
        )
    powers.append(variable_powers)


def specialize(records, fixed_values):
    coefficients = {}
    for keep_exponents, fixed_exponents, coefficient in records:
        value = coefficient
        for fixed_value, exponent in zip(fixed_values, fixed_exponents):
            if exponent:
                value = (value * pow(fixed_value, exponent, p)) % p
        coefficients[keep_exponents] = (coefficients.get(keep_exponents, 0) + value) % p
    return {exponents: coefficient for exponents, coefficient in coefficients.items() if coefficient}


def evaluate(poly, indices):
    values = np.zeros(indices.shape[0], dtype=np.uint16)
    for exponents, coefficient in poly.items():
        term = np.full(indices.shape[0], coefficient, dtype=np.uint16)
        for variable_index, exponent in enumerate(exponents):
            if exponent:
                term = (term * powers[variable_index][exponent][indices]) % p
        values += term
        values %= p
    return values


def seed_values(seed):
    if args.deterministic_fixed:
        fixed_domains = [
            list(range(1, p)) if name == "rho"
            else list(range(2, p)) if name == "lam"
            else list(range(p))
            for name in fixed_names
        ]
        total = prod(len(domain) for domain in fixed_domains)
        if not 1 <= seed <= total:
            raise RuntimeError(f"deterministic slice index must lie in 1..{total}")
        index = seed - 1
        values = {}
        for name, domain in reversed(list(zip(fixed_names, fixed_domains))):
            values[name] = domain[index % len(domain)]
            index //= len(domain)
        return values
    rng = random.Random(seed)
    values = {name: rng.randrange(p) for name in fixed_names}
    if "rho" in values and values["rho"] == 0:
        values["rho"] = 1 + rng.randrange(p - 1)
    while "lam" in values and values["lam"] in (0, 1):
        values["lam"] = rng.randrange(p)
    return values


root = Path(args.dir)
root.mkdir(parents=True, exist_ok=True)
total_hits = 0
started = time.monotonic()

for seed in range(args.seed_start, args.seed_end + 1):
    slice_started = time.monotonic()
    fixed = seed_values(seed)
    fixed_tuple = tuple(fixed[name] for name in fixed_names)
    sliced_equations = [specialize(records, fixed_tuple) for records in equation_records]
    sliced_open = [specialize(records, fixed_tuple) for records in open_records]

    indices = base_indices
    valid_slice = True
    for polynomial in sliced_open:
        if not polynomial:
            valid_slice = False
            indices = indices[:0]
            break
        indices = indices[evaluate(polynomial, indices) != 0]
        if indices.shape[0] == 0:
            break

    survivor_counts = []
    if valid_slice:
        # The first residual is always the sparse one after specialization;
        # retain a term-count sort as a guard against future chart changes.
        order = sorted(range(len(sliced_equations)), key=lambda i: len(sliced_equations[i]))
        for equation_index in order:
            polynomial = sliced_equations[equation_index]
            if not polynomial:
                raise RuntimeError(f"seed {seed} made residual equation {equation_index} zero")
            indices = indices[evaluate(polynomial, indices) == 0]
            survivor_counts.append((equation_index, indices.shape[0]))
            if indices.shape[0] == 0:
                break

    hits = []
    for index in indices[: args.max_hits]:
        hits.append({name: int(grid[i, index]) for i, name in enumerate(keep_names)})
    total_hits += indices.shape[0]

    meta = root / f"{args.prefix}-seed{seed}.meta.txt"
    meta.write_text(
        f"seed={seed}\n"
        f"kill={fixed_names!r}\n"
        f"values={fixed!r}\n"
        f"remaining={keep_names!r}\n"
        "saturate=[]\n"
    )
    log = root / f"{args.prefix}-seed{seed}.scan.log"
    log_lines = [
        f"MW3A10SCAN|p={p}|vars={','.join(keep_names)}|points={grid.shape[1]}"
        f"|eqs={len(equations)}|nonzero={','.join(sorted(nonzero_keep))}",
        f"MW3A10SCAN_OPEN|survivors={indices.shape[0] if not survivor_counts else 'filtered'}",
    ]
    log_lines.extend(
        f"MW3A10SCAN_EQ|step={step}|source={source}|survivors={count}"
        for step, (source, count) in enumerate(survivor_counts, 1)
    )
    log_lines.extend(
        "MW3A10SCAN_HIT|" + ",".join(f"{name}={hit[name]}" for name in keep_names)
        for hit in hits
    )
    log_lines.append(f"MW3A10SCAN|hits={indices.shape[0]}")
    log.write_text("\n".join(log_lines) + "\n")

    print(
        f"MW3A10MULTI|seed={seed}|fixed="
        + ",".join(f"{name}:{fixed[name]}" for name in fixed_names)
        + f"|hits={indices.shape[0]}"
        f"|seconds={time.monotonic()-slice_started:.3f}",
        flush=True,
    )

print(
    f"MW3A10MULTI|done=1|seeds={args.seed_end-args.seed_start+1}"
    f"|hits={total_hits}|seconds={time.monotonic()-started:.3f}",
    flush=True,
)
