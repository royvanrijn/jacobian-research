#!/usr/bin/env sage-python
"""Build the component-adapted NS0011 pole-two section system modulo p.

The input models have fibres ``I9`` at zero, ``I7`` at one, and ``I3`` at
infinity.  The exact NS0011 rank-one generator has component depths two, one,
and zero at those fibres.  For ``Z`` monic of degree two this permits

    X = (center_0 Z^2 mod t^2, center_1 Z^2 mod (t-1))
        + t^2(t-1) Q,
    Y = t^2(t-1) R.

Consequently the homogeneous cubic is divisible by ``t^4(t-1)^2`` and the
section identity reduces to ``R^2=H``.  The default system retains all nine
non-leading coefficients of R: 16 variables and 18 equations of degree at
most six.  A leading block of R coefficients can optionally be eliminated.

This is an exact finite-field equation system.  Solving it does not by itself
prove a characteristic-zero section or a rational NS0011 source family.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0011-source-ansatz-mod5.json"
)
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0011-pole2-section-modp"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--example", type=int, default=0)
parser.add_argument(
    "--branch",
    type=int,
    default=0,
    help="index among infinity x-values with nonzero square right-hand side",
)
parser.add_argument(
    "--zero-y-branch",
    type=int,
    help="index among smooth zero-y infinity points; uses the nontriangular chart",
)
parser.add_argument("--eliminate-r-count", type=int, default=0)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_path = args.input.resolve()
payload = json.loads(input_path.read_text())
prime = ZZ(payload["prime"])
if not is_prime(prime) or prime in (2, 3):
    raise SystemExit("input characteristic must be a prime other than 2 or 3")
examples = payload["examples"]
if not 0 <= args.example < len(examples):
    raise SystemExit(f"--example must lie in [0,{len(examples) - 1}]")
if not 0 <= args.eliminate_r_count <= 9:
    raise SystemExit("--eliminate-r-count must lie in [0,9]")
if args.zero_y_branch is not None and args.eliminate_r_count:
    raise SystemExit("--zero-y-branch requires --eliminate-r-count 0")
example = examples[args.example]

field = GF(prime)
base_ring = PolynomialRing(field, "t")
t0 = base_ring.gen()
A0 = base_ring(example["A_coefficients_low_to_high"])
B0 = base_ring(example["B_coefficients_low_to_high"])


def formal_center(A, B, point, fibre_order, precision):
    """Return the nodal-center jet in the local coordinate t-point."""
    shifted_A = base_ring(A(t0 + point))
    shifted_B = base_ring(B(t0 + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    if node**3 + shifted_A[0] * node + shifted_B[0]:
        raise ArithmeticError("nodal x-coordinate is not on the cubic")
    if 3 * node**2 + shifted_A[0]:
        raise ArithmeticError("nodal x-coordinate is not a double root")
    series_ring = PowerSeriesRing(field, "s", default_prec=fibre_order + 2)
    center = series_ring(node)
    shifted_A_series = series_ring(shifted_A)
    for unused in range(6):
        center = (center + (-shifted_A_series / 3) / center) / 2
    if (center**2 + shifted_A_series / 3).valuation() < fibre_order + 1:
        raise ArithmeticError("formal nodal center did not converge")
    centered_cubic = center**3 + shifted_A_series * center + series_ring(shifted_B)
    if centered_cubic.valuation() != fibre_order:
        raise ArithmeticError("centered cubic has the wrong discriminant order")
    return base_ring([field(center[index]) for index in range(precision)]), node


center_zero, node_zero = formal_center(A0, B0, field.zero(), 9, 2)
center_one_shifted, node_one = formal_center(A0, B0, field.one(), 7, 1)
center_one = center_one_shifted[0]

# Fix one leading smooth point on the I3 cubic.  One y-sign suffices because
# negating Y preserves X and the component data.  Zero-y points need a
# separate chart because the triangular R elimination divides by 2*ell.
branches = []
zero_y = []
node_infinity = -field(3) * B0[12] / (field(2) * A0[8])
for c_value in field:
    infinity_rhs = c_value**3 + A0[8] * c_value + B0[12]
    if infinity_rhs == 0:
        zero_y.append(
            {"x": int(c_value), "is_node": bool(c_value == node_infinity)}
        )
    elif infinity_rhs.is_square():
        branches.append((int(c_value), int(infinity_rhs.sqrt())))
smooth_zero_y = [point for point in zero_y if not point["is_node"]]
if args.zero_y_branch is None:
    if not 0 <= args.branch < len(branches):
        raise SystemExit(f"--branch must lie in [0,{len(branches) - 1}]")
    c_int, ell_int = branches[args.branch]
else:
    if not 0 <= args.zero_y_branch < len(smooth_zero_y):
        raise SystemExit(
            f"--zero-y-branch must lie in [0,{len(smooth_zero_y) - 1}]"
        )
    c_int = smooth_zero_y[args.zero_y_branch]["x"]
    ell_int = 0

remaining_top = 8 - args.eliminate_r_count
names = ["z1", "z0"] + [f"q{i}" for i in range(4, -1, -1)]
names += [f"r{i}" for i in range(remaining_top, -1, -1)]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
generators = coefficient_ring.gens_dict()
polynomial_ring = PolynomialRing(coefficient_ring, "t")
t = polynomial_ring.gen()

Z = t**2 + generators["z1"] * t + generators["z0"]
W = Z**2
center_zero_polynomial = polynomial_ring(
    [coefficient_ring(value) for value in center_zero.list()]
)
local_zero = polynomial_ring((center_zero_polynomial * W) % t**2)
target_one = coefficient_ring(center_one) * coefficient_ring(Z(1)) ** 2
crt_remainder = local_zero + (target_one - coefficient_ring(local_zero(1))) * t**2
Q = coefficient_ring(c_int) * t**5 + sum(
    generators[f"q{i}"] * t**i for i in range(5)
)
D = t**2 * (t - 1)
X = crt_remainder + D * Q
A = polynomial_ring([coefficient_ring(value) for value in A0.list()])
B = polynomial_ring([coefficient_ring(value) for value in B0.list()])
rhs = X**3 + A * X * W**2 + B * W**3
D2 = D**2
H, remainder = rhs.quo_rem(D2)
if remainder:
    raise ArithmeticError("component-adapted cubic is not divisible by t^4(t-1)^2")
expected_h_degree = 17 if args.zero_y_branch is not None else 18
if H.degree() > expected_h_degree:
    raise ArithmeticError("unexpected quotient degree")
ell = coefficient_ring(ell_int)
if coefficient_ring(H[18]) != ell**2:
    raise ArithmeticError("selected infinity branch does not clear the leading equation")

r_coefficients = {9: ell}
for index in range(remaining_top + 1):
    r_coefficients[index] = generators[f"r{index}"]
solved_r = {}
for index in range(8, remaining_top, -1):
    degree = 9 + index
    other = coefficient_ring.zero()
    for left, left_value in r_coefficients.items():
        right = degree - left
        if right in r_coefficients:
            other += left_value * r_coefficients[right]
    solved = coefficient_ring((H[degree] - other) / (2 * ell))
    solved_r[index] = solved
    r_coefficients[index] = solved

highest_equation_degree = 17 - args.eliminate_r_count
equations = []
for degree in range(highest_equation_degree + 1):
    square_coefficient = coefficient_ring.zero()
    for left, left_value in r_coefficients.items():
        right = degree - left
        if right in r_coefficients:
            square_coefficient += left_value * r_coefficients[right]
    equation = coefficient_ring(square_coefficient - H[degree])
    if not equation:
        raise ArithmeticError(f"unexpected zero equation in degree {degree}")
    equations.append(equation)

R = sum(r_coefficients[index] * t**index for index in range(10))
output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
branch_tag = (
    f"zeroY{args.zero_y_branch:02d}"
    if args.zero_y_branch is not None
    else f"branch{args.branch:02d}"
)
stem = (
    f"example{args.example:02d}-{branch_tag}"
    f"-c{c_int}-l{ell_int}-elimR{args.eliminate_r_count:02d}"
)
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations)
msolve_text += "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0011-pole2-section-modp-system.v1",
    "status": "PASS_EXACT_COMPONENT_ADAPTED_MODULAR_SECTION_SYSTEM",
    "input": {
        "artifact": relative(input_path),
        "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "example_index": args.example,
    },
    "prime": int(prime),
    "fibres": {
        "0": {
            "type": "I9",
            "component_depth": 2,
            "nodal_x": int(node_zero),
            "center_coefficients_mod_t2": [int(value) for value in center_zero],
        },
        "1": {
            "type": "I7",
            "component_depth": 1,
            "nodal_x": int(node_one),
            "center_constant": int(center_one),
        },
        "infinity": {"type": "I3", "component_depth": 0},
    },
    "infinity": {
        "nonzero_y_branches_up_to_sign": [list(branch) for branch in branches],
        "selected_nonzero_y_branch": (
            None if args.zero_y_branch is not None else args.branch
        ),
        "selected_smooth_zero_y_branch": args.zero_y_branch,
        "selected_x_y": [c_int, ell_int],
        "nodal_x": int(node_infinity),
        "zero_y_points_excluded_from_this_chart": zero_y,
    },
    "system": {
        "variables": names,
        "base_variables_before_R": names[:7],
        "variable_count": len(names),
        "equation_count": len(equations),
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "eliminated_R_coefficient_count": args.eliminate_r_count,
        "eliminated_R_coefficients_high_to_low": [
            str(solved_r[index])
            for index in range(8, remaining_top, -1)
        ],
        "Z": str(Z),
        "X": str(X),
        "R": str(R),
        "H": str(H),
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "Exact algebra over the displayed finite field builds the pole-two "
        "section chart with I9 depth two, I7 depth one, and I3 identity-component "
        "data. The chart fixes a nonzero-y leading point and excludes the listed "
        "zero-y points. A solution must still be found, lifted to characteristic "
        "zero, and tied to a rational one-parameter NS0011 source family."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text:
        raise SystemExit("modular section msolve input is stale")
    if metadata_path.read_text() != metadata_text:
        raise SystemExit("modular section metadata is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0011SECTIONBUILD|"
    f"example={args.example}|branch={branch_tag}|nonzero_branches={len(branches)}|"
    f"elimR={args.eliminate_r_count}|vars={len(names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
