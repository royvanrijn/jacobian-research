#!/usr/bin/env sage -python
"""Recover the marked NS0024 MW4 source component over a finite field.

The equation chart has split fibres I7, I5, and I4 at 0, 1, and infinity.
Three polynomial sections and one P.O=1 section are imposed using the exact
root-adapted MW4 frame.  The script deliberately works only in finite
characteristic; it is not a characteristic-zero reconstruction tool.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
)


def branch_jet(ring, coefficients, root, root_inverse, precision):
    """Return 2*(-A/3)^(3/2), with its square-root branch marked by root."""
    unit = [-value / ring(3) for value in coefficients[:precision]]
    square_root = [root]
    for degree in range(1, precision):
        known = sum(
            square_root[left] * square_root[degree - left]
            for left in range(1, degree)
        )
        square_root.append((unit[degree] - known) * root_inverse / ring(2))
    square = [
        sum(square_root[left] * square_root[degree - left] for left in range(degree + 1))
        for degree in range(precision)
    ]
    cube = [
        sum(square[left] * square_root[degree - left] for left in range(degree + 1))
        for degree in range(precision)
    ]
    return [2 * value for value in cube]


def coefficient_equations(polynomial):
    return list(polynomial)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=11)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument(
    "--slice-a1",
    type=int,
    help="optional transverse finite-field slice used only to find closed points",
)
parser.add_argument(
    "--groebner",
    action="store_true",
    help="compute the grevlex Groebner basis and report the component dimension",
)
parser.add_argument("--export-msolve", type=Path)
args = parser.parse_args()

source = json.loads(args.source.resolve().read_text())
assert source["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_MW4_SOURCE_AND_NIEMEIER_CERTIFICATE"
assert source["source"]["root_type"] == "A3+A4+A6"
assert source["source"]["mw_rank_for_rho_19"] == 4

prime = ZZ(args.prime)
assert prime.is_prime() and prime not in (2, 3, 5, 7)
field = GF(prime)

surface_names = [f"a{index}" for index in range(1, 9)] + [
    "r1", "ri", "s1", "si", "b7", "b8",
]
p1_names = ["p1x1", "p1x2"] + [f"p1y{index}" for index in range(1, 5)]
p2_names = [f"p2x{index}" for index in range(1, 5)] + [
    f"p2y{index}" for index in range(1, 7)
]
p3_names = [f"p3x{index}" for index in range(4)] + [
    f"p3y{index}" for index in range(6)
]
p4_names = ["c", "c0_inverse", "c1_inverse"] + [
    f"p4x{index}" for index in range(1, 7)
] + [f"p4y{index}" for index in range(1, 10)]
names = surface_names + p1_names + p2_names + p3_names + p4_names
coefficient_ring = PolynomialRing(field, names, order="degrevlex")
variable = coefficient_ring.gens_dict()
polynomial_ring = PolynomialRing(coefficient_ring, "t")
t = polynomial_ring.gen()

a = [coefficient_ring(-3)] + [variable[f"a{index}"] for index in range(1, 9)]
r1, ri, s1, si, b7, b8 = [
    variable[name] for name in ("r1", "ri", "s1", "si", "b7", "b8")
]

# Seven zero jets and four infinity jets determine eleven B coefficients.
# The two middle coefficients b7,b8 remain free.  Five jets at one then give
# two determining equations and three compatibility equations, exactly the
# five-dimensional fibre stratum described in the foundry report.
zero_branch = branch_jet(coefficient_ring, a, coefficient_ring.one(), coefficient_ring.one(), 7)
infinity_branch = branch_jet(coefficient_ring, list(reversed(a)), ri, si, 4)
b = zero_branch + [b7, b8] + list(reversed(infinity_branch))
assert len(b) == 13

a_at_one = [
    sum(a[index] * binomial(index, jet) for index in range(jet, 9))
    for jet in range(5)
]
b_at_one_target = branch_jet(coefficient_ring, a_at_one, r1, s1, 5)
b_at_one_actual = [
    sum(b[index] * binomial(index, jet) for index in range(jet, 13))
    for jet in range(5)
]

equations = [
    a[8] + 3 * ri**2,
    sum(a) + 3 * r1**2,
    r1 * s1 - 1,
    ri * si - 1,
] + [b_at_one_actual[jet] - b_at_one_target[jet] for jet in range(5)]

A = polynomial_ring(a)
B = polynomial_ring(b)

# Exact minimum-complexity marked basis in the root-adapted MW quotient:
#   P1=(-1,0,0,0), P2=(-1,0,1,0), P3=(-1,1,0,0) have P.O=0;
#   P4=(-1,0,0,1) has P.O=1.
# Their resolved component labels are respectively
#   (I5:2,I4:2,I7:6), (I7:5), (I4:1), (I7:1),
# in the simple-root orientations stored in the source certificate.
p1x1, p1x2 = variable["p1x1"], variable["p1x2"]
p1x3 = r1 - ri - 1 - p1x1 - p1x2
X1 = 1 + p1x1 * t + p1x2 * t**2 + p1x3 * t**3 + ri * t**4
p1y = [variable[f"p1y{index}"] for index in range(1, 5)]
Y1 = sum(p1y[index - 1] * t**index for index in range(1, 5)) - sum(p1y) * t**5

X2 = 1 + sum(variable[f"p2x{index}"] * t**index for index in range(1, 5))
Y2 = sum(variable[f"p2y{index}"] * t**index for index in range(1, 7))

X3 = sum(variable[f"p3x{index}"] * t**index for index in range(4)) + ri * t**4
Y3 = sum(variable[f"p3y{index}"] * t**index for index in range(6))

c = variable["c"]
h = t - c
X4 = c**2 + sum(variable[f"p4x{index}"] * t**index for index in range(1, 7))
Y4 = sum(variable[f"p4y{index}"] * t**index for index in range(1, 10))
equations += [c * variable["c0_inverse"] - 1, (c - 1) * variable["c1_inverse"] - 1]

equations += coefficient_equations(Y1**2 - X1**3 - A * X1 - B)
equations += coefficient_equations(Y2**2 - X2**3 - A * X2 - B)
equations += coefficient_equations(Y3**2 - X3**3 - A * X3 - B)
equations += coefficient_equations(Y4**2 - X4**3 - A * X4 * h**4 - B * h**6)

if args.slice_a1 is not None:
    equations.append(a[1] - field(args.slice_a1))

ideal = coefficient_ring.ideal(equations)
print(
    "NS0024MW4SYSTEM|p={}|variables={}|equations={}|slice_a1={}".format(
        prime, coefficient_ring.ngens(), len(equations), args.slice_a1
    ),
    flush=True,
)

if args.export_msolve is not None:
    output = args.export_msolve.resolve()
    output.write_text(
        ",".join(names) + "\n"
        + str(prime) + "\n"
        + ",\n".join(str(equation) for equation in equations) + "\n"
    )
    print(f"NS0024MW4MSOLVE|output={output}", flush=True)

if args.groebner:
    basis = ideal.groebner_basis(algorithm="libsingular:slimgb")
    print(
        "NS0024MW4GROEBNER|p={}|basis={}|dimension={}|unit={}".format(
            prime, len(basis), ideal.dimension(), int(coefficient_ring.one() in ideal)
        ),
        flush=True,
    )
