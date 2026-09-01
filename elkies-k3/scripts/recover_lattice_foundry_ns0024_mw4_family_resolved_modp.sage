#!/usr/bin/env sage -python
"""Recover the resolved NS0024 MW4 marked family over a finite field.

This is the component-chart version of
``recover_lattice_foundry_ns0024_mw4_family_modp.sage``.  It incorporates two
lessons from the successful Q80/H92 computations:

* parameterize the I7, I5 and I4 discriminant jets by their formal double-root
  branches instead of asking a Groebner basis to rediscover those branches;
* impose a component of depth ``d=min(k,n-k)`` by substituting
  ``X = center (mod s^d)`` before forming a section-square system.

The exact MW basis has component profiles

    P1=(1,0,0), P2=(2,1,3), P3=(2,1,1), P4=(1,1,1).

Thus the only new coefficient constraints beyond node incidence are the
depth-two congruences for P2 and P3 at the I7 fibre.  All Y coefficients are
then eliminated recursively.  The result is a 31-variable resolved incidence
system, compared with 48 variables in the raw singular chart.

This script is finite-field only.  Exact component orientations, section
intersections, fibre orders, and the q4/orbit1 child remain certification gates
after a family component has been recovered.
"""

import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", required=True, type=int)
parser.add_argument("--slice-a1", type=int)
parser.add_argument("--export-msolve", type=Path)
parser.add_argument("--groebner", action="store_true")
parser.add_argument("--surface-only", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3, 5, 7):
    raise SystemExit("require a prime away from 2,3,5,7")
field = GF(prime)

surface_names = (
    [f"a{index}" for index in range(1, 7)]
    + ["r1", "ri", "r1_inverse", "ri_inverse"]
)
p1_names = [f"p1x{index}" for index in range(1, 5)] + ["p1y1", "p1y1_inverse"]
p2_names = ["p2x2", "p2y2", "p2y2_inverse"]
p3_names = ["p3x2", "p3y2", "p3y2_inverse"]
p4_names = (
    ["c", "c_inverse", "c1_inverse"]
    + [f"p4x{index}" for index in range(1, 5)]
    + ["p4y1", "p4y1_inverse"]
)
names = surface_names + p1_names + p2_names + p3_names + p4_names
ring = PolynomialRing(field, names=names, order="degrevlex")
v = ring.gens_dict()
polynomials = PolynomialRing(ring, "t")
t = polynomials.gen()


def formal_center(a_jet, root, root_inverse, precision):
    """Return C with C(0)=root and C^2=-A/3 modulo s^precision."""
    center = [ring(root)]
    for degree in range(1, precision):
        known = sum(
            center[left] * center[degree - left]
            for left in range(1, degree)
        )
        center.append((-a_jet[degree] / 3 - known) * root_inverse / 2)
    return center


def cube_coefficients(values, precision):
    local = PolynomialRing(ring, "s")
    s = local.gen()
    cube = 2 * local(values) ** 3
    return [ring(cube[index]) for index in range(precision)]


def recursively_squared_section(prefix, rhs, initial_degree, maximum_y_degree):
    """Eliminate Y coefficientwise after its first nonzero coefficient."""
    lead = v[f"{prefix}y{initial_degree}"]
    lead_inverse = v[f"{prefix}y{initial_degree}_inverse"]
    y = {index: ring.zero() for index in range(initial_degree)}
    y[initial_degree] = lead
    equations = [rhs[index] for index in range(2 * initial_degree)]
    equations += [lead**2 - rhs[2 * initial_degree], lead * lead_inverse - 1]
    for index in range(initial_degree + 1, maximum_y_degree + 1):
        degree = initial_degree + index
        known = sum(
            y[left] * y[degree - left]
            for left in range(initial_degree + 1, index)
        )
        y[index] = (rhs[degree] - known) * lead_inverse / 2
    Y = sum(y[index] * t**index for index in range(maximum_y_degree + 1))
    residual = Y**2 - rhs
    for degree in range(initial_degree + maximum_y_degree + 1,
                        max(Y.degree(), rhs.degree()) + 1):
        equations.append(residual[degree])
    return Y, [ring(equation) for equation in equations if equation != 0]


r1, ri = v["r1"], v["ri"]
r1_inverse, ri_inverse = v["r1_inverse"], v["ri_inverse"]
a = [ring(-3)] + [v[f"a{index}"] for index in range(1, 7)]
a += [None, -3 * ri**2]
a[7] = -3 * r1**2 - sum(a[:7]) - a[8]
A = polynomials(a)

equations = [r1 * r1_inverse - 1, ri * ri_inverse - 1]

# Formal double-root branches at 0, 1, and infinity.  The associated B branch
# is 2*C^3.  Seven bottom and four top coefficients leave b7,b8, which the
# value and first derivative at t=1 determine uniquely.  The remaining three
# I5 jet equations cut out the five-dimensional fibre stratum.
center0 = formal_center(a[:7], ring.one(), ring.one(), 7)
a_at_one = [
    sum(a[index] * binomial(index, jet) for index in range(jet, 9))
    for jet in range(5)
]
center1 = formal_center(a_at_one, r1, r1_inverse, 5)
a_at_infinity = [a[8 - index] for index in range(4)]
center_infinity = formal_center(a_at_infinity, ri, ri_inverse, 4)
b_bottom = cube_coefficients(center0, 7)
b_at_one = cube_coefficients(center1, 5)
b_top_reversed = cube_coefficients(center_infinity, 4)
b = b_bottom + [None, None] + list(reversed(b_top_reversed))
known_value = sum(b[index] for index in range(13) if index not in (7, 8))
known_derivative = sum(
    index * b[index] for index in range(13) if index not in (7, 8)
)
value_rhs = b_at_one[0] - known_value
derivative_rhs = b_at_one[1] - known_derivative
b[8] = derivative_rhs - 7 * value_rhs
b[7] = 8 * value_rhs - derivative_rhs
B = polynomials(b)
for jet in range(2, 5):
    equations.append(
        sum(b[index] * binomial(index, jet) for index in range(jet, 13))
        - b_at_one[jet]
    )

if args.surface_only:
    surface_equations = [ring(equation) for equation in equations]
    print(
        f"NS0024FIBREPARAM|p={prime}|variables={len(surface_names)}"
        f"|equations={len(surface_equations)}|degrees="
        + ",".join(str(equation.degree()) for equation in surface_equations)
        + "|terms="
        + ",".join(str(len(equation.dict())) for equation in surface_equations),
        flush=True,
    )
    if args.export_msolve is not None:
        output = args.export_msolve.resolve()
        output.write_text(
            ",".join(surface_names) + "\n" + str(prime) + "\n"
            + ",\n".join(map(str, surface_equations)) + "\n"
        )
    raise SystemExit(0)

# P1: component depth one at I7, identity components at I5 and I4.
X1 = 1 + sum(v[f"p1x{index}"] * t**index for index in range(1, 5))
rhs1 = X1**3 + A * X1 + B
Y1, section_equations = recursively_squared_section("p1", rhs1, 1, 6)
equations += section_equations


def depth_two_all_node_section(prefix):
    # X=center0 mod t^2 enforces component +/-2 at I7.  Node incidence at
    # t=1 and infinity supplies the remaining two coefficients.
    x1 = center0[1]
    x2 = v[f"{prefix}x2"]
    x3 = r1 - 1 - ri - x1 - x2
    X = 1 + x1 * t + x2 * t**2 + x3 * t**3 + ri * t**4
    rhs = X**3 + A * X + B
    Y, local_equations = recursively_squared_section(prefix, rhs, 2, 6)
    local_equations += [Y(t=1), Y[6]]
    return X, Y, local_equations


X2, Y2, p2_equations = depth_two_all_node_section("p2")
X3, Y3, p3_equations = depth_two_all_node_section("p3")
equations += p2_equations + p3_equations

# P4: simple pole at t=c and component depth one at all three reducible
# fibres.  The three node values determine x0,x5,x6; recurrence eliminates Y.
c = v["c"]
h = t - c
p4x = [v[f"p4x{index}"] for index in range(1, 5)]
x5 = r1 * (1 - c)**2 - c**2 - ri - sum(p4x)
X4 = (
    c**2
    + sum(p4x[index - 1] * t**index for index in range(1, 5))
    + x5 * t**5
    + ri * t**6
)
rhs4 = X4**3 + A * X4 * h**4 + B * h**6
Y4, p4_equations = recursively_squared_section("p4", rhs4, 1, 8)
equations += p4_equations + [Y4(t=1)]
equations += [c * v["c_inverse"] - 1, (c - 1) * v["c1_inverse"] - 1]

if args.slice_a1 is not None:
    equations.append(v["a1"] - field(args.slice_a1))
equations = [ring(equation) for equation in equations if equation != 0]
ideal = ring.ideal(equations)

print(
    f"NS0024MW4RESOLVED|p={prime}|variables={ring.ngens()}"
    f"|equations={len(equations)}|slice_a1={args.slice_a1}",
    flush=True,
)
if args.export_msolve is not None:
    output = args.export_msolve.resolve()
    output.write_text(
        ",".join(names) + "\n" + str(prime) + "\n"
        + ",\n".join(map(str, equations)) + "\n"
    )
    print(f"NS0024MW4RESOLVEDMSOLVE|output={output}", flush=True)
if args.groebner:
    basis = ideal.groebner_basis(algorithm="libsingular:slimgb")
    print(
        f"NS0024MW4RESOLVEDGB|p={prime}|basis={len(basis)}"
        f"|dimension={ideal.dimension()}|unit={int(ring.one() in ideal)}",
        flush=True,
    )
