#!/usr/bin/env sage -python
"""Build reduced finite-field systems for the missing A11 MW direction.

On the exact orbit42 A11 child a section with P.O=4 can be written

    x = X/Z^2,  y = Y/Z^3,

with monic deg(Z)=4, deg(X)<=12 and deg(Y)<=18.  After fixing a nonzero
point (c,l) on the fibre at infinity, leading coefficients of Y can be
determined recursively by the coefficients of X and Z.  Eliminating all 18
leaves 16 variables and 18 residual equations, but those expanded equations
can be very large.  The default eliminates only four leading Y coefficients,
retaining a sparse 30-variable, 32-equation triangular system for msolve.

The output is an msolve input and a JSON construction record.  It is only a
finite-field discovery system; a solution does not prove a characteristic-
zero section or the next A11-to-2A5 equation lift.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
RR_PATH = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=53)
parser.add_argument(
    "--branch",
    type=int,
    default=0,
    help="index among infinity x-values with nonzero square right-hand side",
)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=LOCAL / "q24-a11-missing-mw-p4-modp",
)
parser.add_argument(
    "--eliminate-y-count",
    type=int,
    default=4,
    help="number of leading Y coefficients to eliminate recursively (0..18)",
)
args = parser.parse_args()

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
if not RR_PATH.exists():
    raise SystemExit(f"missing exact A11 certificate: {RR_PATH}")

rr = json.loads(RR_PATH.read_text())
assert rr["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"

F = GF(p)
FT = PolynomialRing(F, "t")
t0 = FT.gen()


def reduce_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError(f"bad denominator at p={p}")
    return F(value.numerator()) / F(value.denominator())


A0 = FT([reduce_q(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]])
B0 = FT([reduce_q(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]])
Delta0 = FT(-16 * (4 * A0**3 + 27 * B0**2))
if A0.degree() != 8 or B0.degree() != 12:
    raise ArithmeticError("the reduced A11 model lost its K3 degrees")

delta_pattern = sorted((int(factor.degree()), int(exponent)) for factor, exponent in Delta0.factor())
if delta_pattern != [(1, 12), (12, 1)]:
    raise ArithmeticError(
        f"p={p} does not retain the direct I12 + degree-12 I1 pattern: {delta_pattern}"
    )

# One y sign suffices: changing l to -l negates the section without changing X,Z.
branches = []
zero_y_branches = []
for c in F:
    rhs = c**3 + A0[8] * c + B0[12]
    if rhs == 0:
        zero_y_branches.append(int(c))
    elif rhs.is_square():
        ell = rhs.sqrt()
        branches.append((int(c), int(ell)))

if args.branch < 0 or args.branch >= len(branches):
    raise SystemExit(f"--branch must lie in [0,{len(branches) - 1}]")
if args.eliminate_y_count < 0 or args.eliminate_y_count > 18:
    raise SystemExit("--eliminate-y-count must lie in [0,18]")
c_int, ell_int = branches[args.branch]
c = F(c_int)
ell = F(ell_int)
assert ell and ell**2 == c**3 + A0[8] * c + B0[12]

remaining_y_top = 17 - args.eliminate_y_count
names = [f"z{i}" for i in range(3, -1, -1)] + [
    f"x{i}" for i in range(11, -1, -1)
] + [f"y{i}" for i in range(remaining_y_top, -1, -1)]
P = PolynomialRing(F, names=names, order="degrevlex")
generators = P.gens_dict()
PT = PolynomialRing(P, "t")
t = PT.gen()

Z = t**4 + sum(generators[f"z{i}"] * t**i for i in range(4))
X = P(c) * t**12 + sum(generators[f"x{i}"] * t**i for i in range(12))
A = PT([P(value) for value in A0.list()])
B = PT([P(value) for value in B0.list()])
rhs = X**3 + A * X * Z**4 + B * Z**6

# Hensel-style expansion at infinity.  Retained lower Y variables keep the
# system sparse; only the requested leading block is substituted away.
Y = P(ell) * t**18 + sum(
    generators[f"y{i}"] * t**i for i in range(remaining_y_top + 1)
)
solved_y = {}
for j in range(17, remaining_y_top, -1):
    residual = PT(rhs - Y**2)
    yj = P(residual[18 + j] / (2 * ell))
    solved_y[j] = yj
    Y += yj * t**j

final_residual = PT(Y**2 - rhs)
highest_residual_degree = 35 - args.eliminate_y_count
equations = [P(final_residual[k]) for k in range(highest_residual_degree + 1)]
if any(final_residual[k] for k in range(highest_residual_degree + 1, 37)):
    raise ArithmeticError("infinity recursion did not clear the requested leading block")
if any(not equation for equation in equations):
    raise ArithmeticError("unexpected identically-zero residual equation")

args.output_dir.mkdir(parents=True, exist_ok=True)
stem = (
    f"p{p}-branch{args.branch:02d}-c{c_int}-l{ell_int}"
    f"-elimY{args.eliminate_y_count:02d}"
)
msolve_path = args.output_dir / f"{stem}.ms"
metadata_path = args.output_dir / f"{stem}.json"

with msolve_path.open("w") as handle:
    handle.write(",".join(names) + "\n")
    handle.write(str(p) + "\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index + 1 < len(equations) else "\n")

payload = {
    "schema": "elkies-k3.h3-q24-a11-missing-mw-p4-modp-system.v1",
    "status": "PASS_A11_MISSING_MW_P4_REDUCED_MODP_SYSTEM",
    "input": {
        "exact_A11_certificate": str(RR_PATH.relative_to(ROOT)),
        "sha256": hashlib.sha256(RR_PATH.read_bytes()).hexdigest(),
    },
    "prime": int(p),
    "good_reduction_pattern": delta_pattern,
    "target_lattice_profile": {
        "mw_coordinates": [0, 0, 0, 1, 0, 1],
        "height": 12,
        "I12_correction": 0,
        "P_dot_O": 4,
        "note": "The equations enumerate all sections in this pole chart; lattice marking is a later filter.",
    },
    "infinity": {
        "nonzero_y_branch_count_up_to_sign": len(branches),
        "branch_index": args.branch,
        "x_leading_coefficient": c_int,
        "y_leading_coefficient": ell_int,
        "zero_y_x_leading_coefficients_excluded_from_this_chart": zero_y_branches,
    },
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "eliminated_Y_coefficient_count": args.eliminate_y_count,
        "retained_Y_coefficient_count": remaining_y_top + 1,
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "Z": str(Z),
        "X": str(X),
        "eliminated_Y_coefficients_high_to_low": [
            str(solved_y[j]) for j in range(17, remaining_y_top, -1)
        ],
        "msolve_input": str(msolve_path.relative_to(ROOT)),
        "msolve_sha256": hashlib.sha256(msolve_path.read_bytes()).hexdigest(),
    },
    "proof_boundary": (
        "Exact algebra over GF(p) fixes the infinity point and eliminates the stated "
        "leading Y block from the P.O=4 section identity. The chart excludes "
        "the listed zero-y infinity points. A modular solution is discovery evidence only; "
        "it must be marked, lifted to characteristic zero, and checked on the exact A11 "
        "equation before use in a resolved q8 Riemann--Roch lift."
    ),
}
metadata_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11P4MOD|prime={}|branch={}/{}|c={}|ell={}|elimY={}|vars={}|eqs={}|max_degree={}|status={}".format(
        p,
        args.branch,
        len(branches),
        c_int,
        ell_int,
        args.eliminate_y_count,
        len(names),
        len(equations),
        max(equation.degree() for equation in equations),
        payload["status"],
    ),
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
