#!/usr/bin/env sage -python
"""Build construction-targeted A11 section systems modulo p.

The retained A11-to-2A5 lattice neighbors have height 13, I12 correction 3,
P.O=6, and meet the central component (component 6) of the I12 fibre.  In a
split multiplicative local coordinate s at that fibre this forces

    x - center(s) = O(s^6),    y = O(s^6).

Writing Z monic of degree 6, this script builds the conditions into

    X = (center * Z^2 mod s^6) + s^6 Q,   deg Q <= 10,
    Y = s^6 R,                            deg R <= 18.

After fixing the point at infinity, the uneliminated system has 34 variables
and 36 nontrivial Weierstrass coefficient equations.  A Hensel-style expansion
at infinity can recursively eliminate any leading block of the 18 unknown R
coefficients.  Eliminating all 18 leaves 16 variables and 18 equations.  It is
a finite-field discovery chart, not a characteristic-zero section or
equation-lift certificate.

The ``--pole-order 5 --component-depth 3`` chart targets the exact
construction bridge M=(1,0,0,0,0,1), whose group-law combination with three
already exact identity-shell sections gives the q8/orbit12 section.  The
default remains the direct pole-order-6/component-6 target.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, is_prime


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
    "--free-infinity",
    action="store_true",
    help=(
        "retain the leading x/y coefficients and their fibre equation instead "
        "of fixing one smooth-fibre point; covers every infinity branch at once"
    ),
)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=LOCAL / "q24-a11-q8-component6-modp",
)
parser.add_argument("--pole-order", type=int, default=6)
parser.add_argument(
    "--component-depth",
    type=int,
    default=6,
    help="min(k,12-k) for the required I12 component",
)
parser.add_argument(
    "--eliminate-r-count",
    type=int,
    default=0,
    help="number of leading R coefficients to eliminate at infinity (0..18)",
)
parser.add_argument(
    "--retain-x-coefficients",
    action="store_true",
    help=(
        "retain the coefficients of X and impose the six component conditions "
        "as separate equations; this enlarges the system but keeps its defining "
        "equations cubic instead of substituting quadratic expressions into X^3"
    ),
)
parser.add_argument(
    "--quadraticize-z-powers",
    action="store_true",
    help=(
        "introduce W=Z^2 coefficient variables, making the component equations "
        "linear and the Weierstrass equations cubic; implies "
        "--retain-x-coefficients"
    ),
)
args = parser.parse_args()

if args.quadraticize_z_powers:
    args.retain_x_coefficients = True

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
pole_order = ZZ(args.pole_order)
component_depth = ZZ(args.component_depth)
if pole_order < 0 or not 1 <= component_depth <= 6:
    raise SystemExit("require --pole-order >= 0 and --component-depth in [1,6]")
r_top = ZZ(3 * pole_order + 6 - component_depth)
max_x_degree = ZZ(2 * pole_order + 4)
max_equation_degree = ZZ(6 * pole_order + 12)
if args.eliminate_r_count < 0 or args.eliminate_r_count > r_top:
    raise SystemExit(f"--eliminate-r-count must lie in [0,{r_top}]")
if args.free_infinity and args.eliminate_r_count:
    raise SystemExit("--free-infinity currently requires --eliminate-r-count 0")
if not RR_PATH.exists():
    raise SystemExit(f"missing exact A11 certificate: {RR_PATH}")

rr = json.loads(RR_PATH.read_text())
assert rr["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"

F = GF(p)
FT = PolynomialRing(F, "s")
s0 = FT.gen()


def reduce_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError(f"bad denominator at p={p}")
    return F(value.numerator()) / F(value.denominator())


A_t = FT([reduce_q(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]])
B_t = FT([reduce_q(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]])
Delta_t = FT(-16 * (4 * A_t**3 + 27 * B_t**2))
delta_factors = list(Delta_t.factor())
i12 = [(factor, int(exponent)) for factor, exponent in delta_factors if int(exponent) == 12]
if len(i12) != 1 or i12[0][0].degree() != 1:
    raise ArithmeticError(f"p={p} does not retain one rational I12 fibre")
i12_factor = i12[0][0].monic()
alpha = -i12_factor[0]

# Shift the old base so the I12 fibre is s=0.
A = FT(A_t(s0 + alpha))
B = FT(B_t(s0 + alpha))
Delta = FT(-16 * (4 * A**3 + 27 * B**2))
if Delta.valuation() != 12:
    raise ArithmeticError("shifted discriminant does not have order 12")

x_node = -F(3) * B[0] / (F(2) * A[0])
if x_node**3 + A[0] * x_node + B[0] != 0 or 3 * x_node**2 + A[0] != 0:
    raise ArithmeticError("failed to recover the I12 nodal x-coordinate")

PS = PowerSeriesRing(F, "s", default_prec=14)
s_series = PS.gen()
A_series = PS(A)
B_series = PS(B)
center = PS(x_node)
for unused in range(6):
    center = (center + (-A_series / 3) / center) / 2
if (center**2 + A_series / 3).valuation() < 13:
    raise ArithmeticError("formal center square root did not converge")
g0 = center**3 + A_series * center + B_series
if g0.valuation() != 12:
    raise ArithmeticError("centered multiplicative unit does not start in order 12")
center_trunc = FT([F(center[i]) for i in range(component_depth)])

# With fixed infinity, one y sign suffices because the other negates the
# section.  The free chart instead retains both leading coefficients.
branches = None
zero_y_branches = None
c_int = None
ell_int = None
if not args.free_infinity:
    branches = []
    zero_y_branches = []
    for c_value in F:
        infinity_rhs = c_value**3 + A[8] * c_value + B[12]
        if infinity_rhs == 0:
            zero_y_branches.append(int(c_value))
        elif infinity_rhs.is_square():
            branches.append((int(c_value), int(infinity_rhs.sqrt())))
    if args.branch < 0 or args.branch >= len(branches):
        raise SystemExit(f"--branch must lie in [0,{len(branches) - 1}]")
    c_int, ell_int = branches[args.branch]

remaining_r_top = r_top - 1 - args.eliminate_r_count
names = [f"z{i}" for i in range(pole_order - 1, -1, -1)]
if args.free_infinity:
    names += ["c_inf", "ell_inf"]
if args.quadraticize_z_powers:
    names += [f"w{i}" for i in range(2 * pole_order - 1, -1, -1)]
if args.retain_x_coefficients:
    names += [f"x{i}" for i in range(max_x_degree - 1, -1, -1)]
else:
    q_top = max_x_degree - component_depth
    names += [f"q{i}" for i in range(q_top - 1, -1, -1)]
names += [f"r{i}" for i in range(remaining_r_top, -1, -1)]
P = PolynomialRing(F, names=names, order="degrevlex")
generators = P.gens_dict()
PT = PolynomialRing(P, "s")
s = PT.gen()
c = generators["c_inf"] if args.free_infinity else P(F(c_int))
ell = generators["ell_inf"] if args.free_infinity else P(F(ell_int))

Z = s**pole_order + sum(
    generators[f"z{i}"] * s**i for i in range(pole_order)
)
if args.quadraticize_z_powers:
    W = s ** (2 * pole_order) + sum(
        generators[f"w{i}"] * s**i for i in range(2 * pole_order)
    )
    square_equations = [P((W - Z**2)[i]) for i in range(2 * pole_order)]
else:
    W = Z**2
    square_equations = []
r_coefficients = {
    i: generators[f"r{i}"] for i in range(remaining_r_top + 1)
}
r_coefficients[r_top] = P(ell)
A_P = PT([P(value) for value in A.list()])
B_P = PT([P(value) for value in B.list()])
C_P = PT([P(value) for value in center_trunc.list()])

local_x_remainder = PT((C_P * W) % s**component_depth)
if args.retain_x_coefficients:
    Q = None
    X = P(c) * s**max_x_degree + sum(
        generators[f"x{i}"] * s**i for i in range(max_x_degree)
    )
    component_equations = [
        P((X - local_x_remainder)[i]) for i in range(component_depth)
    ]
    if any(not equation for equation in component_equations):
        raise ArithmeticError("unexpected zero among retained-X component equations")
else:
    Q = P(c) * s**q_top + sum(
        generators[f"q{i}"] * s**i for i in range(q_top)
    )
    X = local_x_remainder + s**component_depth * Q
    component_equations = []
rhs = X**3 + A_P * X * W**2 + B_P * W**3

# Recursively clear the requested leading block coefficientwise.  In R^2,
# r_j first occurs linearly in convolution degree r_top+j with coefficient
# 2*ell; the outer s^(2*component_depth) shift fixes the corresponding
# Weierstrass coefficient.
# Computing the convolution directly avoids repeatedly expanding the whole
# increasingly dense square R^2.
solved_r = {}
for j in range(r_top - 1, remaining_r_top, -1):
    convolution_degree = r_top + j
    other = P.zero()
    for left in range(r_top):
        right = convolution_degree - left
        if 0 <= right < r_top:
            other += r_coefficients[left] * r_coefficients[right]
    rj = P((rhs[2 * component_depth + convolution_degree] - other) / (2 * ell))
    solved_r[j] = rj
    r_coefficients[j] = rj

if not args.retain_x_coefficients and any(
    rhs[k] for k in range(2 * component_depth)
):
    raise ArithmeticError("component ansatz did not clear the required local coefficients")
leading_equation = P(ell**2 - rhs[max_equation_degree])
if not args.free_infinity and leading_equation:
    raise ArithmeticError("selected infinity branch does not clear the leading coefficient")
highest_residual_degree = max_equation_degree - 1 - args.eliminate_r_count
equations = list(square_equations) + list(component_equations)
for coefficient_degree in range(2 * component_depth, highest_residual_degree + 1):
    convolution_degree = coefficient_degree - 2 * component_depth
    y2_coefficient = P.zero()
    for left in range(r_top + 1):
        right = convolution_degree - left
        if 0 <= right <= r_top:
            y2_coefficient += r_coefficients[left] * r_coefficients[right]
    equations.append(P(y2_coefficient - rhs[coefficient_degree]))
if args.free_infinity:
    equations.append(leading_equation)
if any(not equation for equation in equations):
    raise ArithmeticError("unexpected identically-zero equation")

R = sum(r_coefficients[i] * s**i for i in range(r_top + 1))

args.output_dir.mkdir(parents=True, exist_ok=True)
stem = (
    (
        f"p{p}-freeInfinity"
        if args.free_infinity
        else f"p{p}-branch{args.branch:02d}-c{c_int}-l{ell_int}"
    )
    + f"-elimR{args.eliminate_r_count:02d}"
    + (
        ""
        if (pole_order, component_depth) == (6, 6)
        else f"-pO{pole_order}-cpt{component_depth}"
    )
    + ("-keepX" if args.retain_x_coefficients else "")
    + ("-quadZ" if args.quadraticize_z_powers else "")
)
msolve_path = args.output_dir / f"{stem}.ms"
metadata_path = args.output_dir / f"{stem}.json"
with msolve_path.open("w") as handle:
    handle.write(",".join(names) + "\n")
    handle.write(str(p) + "\n")
    for index, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if index + 1 < len(equations) else "\n")

direct_q8 = (pole_order, component_depth) == (6, 6)
target_coset_bridge = (pole_order, component_depth) == (5, 3)
status = (
    "PASS_A11_Q8_COMPONENT6_REDUCED_MODP_SYSTEM"
    if direct_q8
    else "PASS_A11_TARGET_COSET_COMPONENT3_REDUCED_MODP_SYSTEM"
    if target_coset_bridge
    else "PASS_A11_CONSTRUCTION_TARGET_REDUCED_MODP_SYSTEM"
)
payload = {
    "schema": (
        "elkies-k3.h3-q24-a11-q8-component6-section-modp-system.v1"
        if direct_q8
        else "elkies-k3.h3-q24-a11-construction-target-section-modp-system.v2"
    ),
    "status": status,
    "input": {
        "exact_A11_certificate": str(RR_PATH.relative_to(ROOT)),
        "sha256": hashlib.sha256(RR_PATH.read_bytes()).hexdigest(),
    },
    "prime": int(p),
    "I12": {
        "old_base_root": int(alpha),
        "nodal_x": int(x_node),
        "center_coefficients_mod_s_depth": [int(value) for value in center_trunc.list()],
        "center_coefficients_mod_s6": (
            [int(value) for value in center_trunc.list()] if direct_q8 else None
        ),
        "centered_cubic_order": int(g0.valuation()),
        "component_depth": int(component_depth),
        "target_component": 6 if direct_q8 else None,
        "components_up_to_negation": [
            int(component_depth), int(12 - component_depth)
        ],
        "target_correction": str(
            QQ(component_depth * (12 - component_depth)) / 12
        ),
    },
    "target_lattice_profiles": (
        [
            {"orbit": 12, "mw": [0, 0, -1, 0, 0, 1]},
            {"orbit": 113, "mw": [-1, 0, 0, 1, 0, 1]},
            {"orbit": 1535, "mw": [-1, 0, 0, -1, 0, -1]},
            {"orbit": 2162, "mw": [-1, 0, -1, 0, 0, -1]},
        ]
        if direct_q8
        else [{"role": "target_coset_bridge", "mw": [1, 0, 0, 0, 0, 1]}]
        if target_coset_bridge
        else []
    ),
    "target_common_profile": {
        "pole_order": int(pole_order),
        "height": "13" if direct_q8 else "47/4" if target_coset_bridge else None,
        "I12_correction": (
            "3" if direct_q8 else "9/4" if target_coset_bridge else None
        ),
        "required_use": (
            "direct q8 divisor section"
            if direct_q8
            else "M in P12=M+S6-2*S2-2*S8"
            if target_coset_bridge
            else "unspecified construction target"
        ),
    },
    "infinity": {
        "free": bool(args.free_infinity),
        "nonzero_y_branch_count_up_to_sign": (
            None if branches is None else len(branches)
        ),
        "branch_index": None if args.free_infinity else args.branch,
        "x_leading_coefficient": c_int,
        "y_leading_coefficient": ell_int,
        "zero_y_x_leading_coefficients_excluded_from_this_chart": zero_y_branches,
    },
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "eliminated_R_coefficient_count": args.eliminate_r_count,
        "retained_R_coefficient_count": int(remaining_r_top + 1),
        "retained_X_coefficients": bool(args.retain_x_coefficients),
        "quadraticized_Z_powers": bool(args.quadraticize_z_powers),
        "square_equation_count": len(square_equations),
        "component_equation_count": len(component_equations),
        "infinity_equation_count": int(args.free_infinity),
        "equation_coefficient_degrees": list(
            range(2 * component_depth, highest_residual_degree + 1)
        ),
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "Z": str(Z),
        "W": str(W),
        "Q": None if Q is None else str(Q),
        "X": str(X),
        "R": str(R),
        "eliminated_R_coefficients_high_to_low": [
            str(solved_r[j]) for j in range(r_top - 1, remaining_r_top, -1)
        ],
        "X_local_remainder": str(local_x_remainder),
        "msolve_input": str(msolve_path.relative_to(ROOT)),
        "msolve_sha256": hashlib.sha256(msolve_path.read_bytes()).hexdigest(),
    },
    "proof_boundary": (
        f"Exact algebra over GF(p) builds the P.O={pole_order}, "
        f"I12-component-depth-{component_depth} chart "
        "and its recursively reduced section identity. "
        + (
            "The free-infinity chart covers every leading smooth-fibre point. "
            if args.free_infinity
            else "The fixed nonzero-y infinity chart excludes the listed leading x-values. "
        )
        + "A modular solution must still "
        + (
            "be assigned to one of the four lattice markings, compiled to a 2A5 child, "
            if direct_q8
            else "be identified with the selected target-coset bridge and combined "
            "with the exact identity-shell points, "
        )
        + "lifted to characteristic zero, and checked on the exact A11 equation."
    ),
}
metadata_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11Q8C6MOD|prime={}|branch={}/{}|c={}|ell={}|elimR={}|vars={}|eqs={}|max_degree={}|status={}".format(
        p,
        "free" if args.free_infinity else args.branch,
        "all" if branches is None else len(branches),
        "free" if c_int is None else c_int,
        "free" if ell_int is None else ell_int,
        args.eliminate_r_count,
        len(names),
        len(equations),
        max(equation.degree() for equation in equations),
        payload["status"],
    ),
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
