#!/usr/bin/env sage -python
"""Recover a resolved NS0024 MW3/MW4 marking over a finite field.

This is the component-chart version of
``recover_lattice_foundry_ns0024_mw4_family_modp.sage``.  It incorporates two
lessons from the successful Q80/H92 computations:

* parameterize the I7, I5 and I4 discriminant jets by their formal double-root
  branches instead of asking a Groebner basis to rediscover those branches;
* impose a component of depth ``d=min(k,n-k)`` by substituting
  ``X = center (mod s^d)`` before forming a section-square system.

The default resolved-depth basis has absolute component profiles

    Q1=(6,0,0), Q2=(2,1,1), Q3=(4,2,0), Q4=(6,4,3).

The first two polynomial sections use the familiar depth-one/depth-two
charts.  Q3 is forced by ``X=center_0 mod t^3`` and
``X=center_1 mod (t-1)^2``; its remaining infinity identity-component open is
covered by two saturation charts.  In particular the surface and all three
MW3 sections remain variables in the exported ideal.  A closed point of this
ideal may therefore have its entire MW3 marking over ``GF(p^d)``.  This is not
the older fixed-rational-MW3 extension search.

Use ``--mw3-only`` to study the positive-dimensional three-section locus, or
leave P4 enabled to solve all four sections jointly.  Repeatable ``--slice``
constraints cut reproducible zero-dimensional sections without changing the
source field.  Sparse Y coefficients keep the exported identities low-degree.

This script is finite-field only.  Exact component orientations, section
intersections, fibre orders, and the q4/orbit1 child remain certification gates
after a family component has been recovered.
"""

import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", required=True, type=int)
parser.add_argument(
    "--basis-marking",
    choices=("resolved-depth13", "original"),
    default="resolved-depth13",
)
parser.add_argument("--slice-a1", type=int, help="compatibility alias for --slice a1=VALUE")
parser.add_argument(
    "--slice",
    action="append",
    default=[],
    metavar="NAME=VALUE",
    help="fix any named coordinate; repeat to cut a zero-dimensional slice",
)
parser.add_argument(
    "--surface-hyperplane",
    action="append",
    default=[],
    metavar="C1,C2,C3,C4,C5,C6,CR1,CRI,C0",
    help=(
        "add C1*a1+...+C6*a6+CR1*r1+CRI*ri+C0=0; repeatable, "
        "with one generic hyperplane expected for the joint MW4 curve"
    ),
)
parser.add_argument(
    "--mw3-only",
    action="store_true",
    help="omit Q4 and recover the joint surface/MW3 locus",
)
parser.add_argument(
    "--q3-infinity-chart",
    choices=("x", "y"),
    default="y",
    help="resolved Q3 identity-component saturation at infinity",
)
parser.add_argument(
    "--q1-i5-identity-chart", choices=("x", "y"), default="y"
)
parser.add_argument(
    "--q1-i4-identity-chart", choices=("x", "y"), default="y"
)
parser.add_argument("--export-msolve", type=Path)
parser.add_argument(
    "--fixed-rur-anchor",
    action="store_true",
    help="append a prime-independent full-coordinate separating-element anchor",
)
parser.add_argument(
    "--explicit-formal-centers",
    action="store_true",
    help=(
        "retain the I5/I4 formal-center jet coefficients as sparse auxiliary "
        "variables instead of recursively substituting powers of r1^-1 and ri^-1"
    ),
)
parser.add_argument(
    "--split-section-opens",
    action="store_true",
    help=(
        "give each section-chart nonvanishing condition its own inverse; "
        "this adds variables but keeps the saturation equations low-degree"
    ),
)
parser.add_argument("--groebner", action="store_true")
parser.add_argument("--surface-only", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3, 5, 7):
    raise SystemExit("require a prime away from 2,3,5,7")
field = GF(prime)

open_inverse_names = [
    "exact_i7_order_inverse", "exact_i5_order_inverse", "exact_i4_order_inverse"
]
section_open_inverse_names = {}
if not args.surface_only:
    if args.split_section_opens:
        section_open_inverse_names = {
            "q1": (
                "q1_i7_tangent_inverse",
                "q1_i5_identity_inverse",
                "q1_i4_identity_inverse",
            ),
            "q2": (
                "q2_i7_tangent_inverse",
                "q2_i5_tangent_inverse",
                "q2_i4_tangent_inverse",
            ),
            "q3": (
                "q3_i7_tangent_inverse",
                "q3_i5_tangent_inverse",
                "q3_i4_identity_inverse",
            ),
        }
        if not args.mw3_only:
            section_open_inverse_names["q4"] = (
                "q4_i7_tangent_inverse",
                "q4_i5_tangent_inverse",
                "q4_i4_tangent_inverse",
                "q4_pole_inverse",
                "q4_c0_inverse",
                "q4_c1_inverse",
            )
    else:
        section_open_inverse_names = {
            "q1": ("q1_chart_inverse",) * 3,
            "q2": ("q2_chart_inverse",) * 3,
            "q3": ("q3_chart_inverse",) * 3,
        }
        if not args.mw3_only:
            section_open_inverse_names["q4"] = ("q4_chart_inverse",) * 6
    open_inverse_names += sorted(
        {
            inverse_name
            for inverse_names in section_open_inverse_names.values()
            for inverse_name in inverse_names
        }
    )
surface_names = (
    [f"a{index}" for index in range(1, 7)]
    + ["r1", "ri", "r1_inverse", "ri_inverse"]
    + (
        [f"center1_{index}" for index in range(1, 6)]
        + [f"center_infinity_{index}" for index in range(1, 5)]
        if args.explicit_formal_centers else []
    )
    + open_inverse_names
)
p1_names = (
    [f"p1x{index}" for index in range(1, 5)]
    + [f"p1y{index}" for index in range(1, 7)]
)
p2_names = ["p2x2"] + [f"p2y{index}" for index in range(2, 5)]
p3_names = (
    ["p3y3", "p3y4"]
    if args.basis_marking == "resolved-depth13"
    else ["p3x2"] + [f"p3y{index}" for index in range(2, 5)]
)
p4_names = [] if args.mw3_only else (
    ["c"]
    + [f"p4x{index}" for index in range(1, 5)]
    + [f"p4y{index}" for index in range(1, 8)]
)
anchor_names = ["rur_anchor"] if args.fixed_rur_anchor else []
names = surface_names + p1_names + p2_names + p3_names + p4_names + anchor_names
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


def retained_formal_center(prefix, a_jet, root, precision):
    """Return a sparse formal center and its exact quadratic recurrence."""
    center = [ring(root)] + [v[f"{prefix}_{degree}"] for degree in range(1, precision)]
    recurrence = []
    for degree in range(1, precision):
        square_coefficient = sum(
            center[left] * center[degree - left]
            for left in range(degree + 1)
        )
        recurrence.append(3 * square_coefficient + a_jet[degree])
    return center, recurrence


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


def coefficient_equations(polynomial):
    """Return the nonzero coefficient equations without recursive substitution."""
    return [ring(coefficient) for coefficient in polynomial if coefficient != 0]


r1, ri = v["r1"], v["ri"]
r1_inverse, ri_inverse = v["r1_inverse"], v["ri_inverse"]
a = [ring(-3)] + [v[f"a{index}"] for index in range(1, 7)]
a += [None, -3 * ri**2]
a[7] = -3 * r1**2 - sum(a[:7]) - a[8]
A = polynomials(a)

equations = [r1 * r1_inverse - 1, ri * ri_inverse - 1]
open_factors = []

# Formal double-root branches at 0, 1, and infinity.  The associated B branch
# is 2*C^3.  Seven bottom and four top coefficients leave b7,b8, which the
# value and first derivative at t=1 determine uniquely.  The remaining three
# I5 jet equations cut out the five-dimensional fibre stratum.
center0 = formal_center(a[:8], ring.one(), ring.one(), 8)
a_at_one = [
    sum(a[index] * binomial(index, jet) for index in range(jet, 9))
    for jet in range(6)
]
if args.explicit_formal_centers:
    center1, center1_equations = retained_formal_center(
        "center1", a_at_one, r1, 6
    )
    equations += center1_equations
else:
    center1 = formal_center(a_at_one, r1, r1_inverse, 6)
a_at_infinity = [a[8 - index] for index in range(5)]
if args.explicit_formal_centers:
    center_infinity, center_infinity_equations = retained_formal_center(
        "center_infinity", a_at_infinity, ri, 5
    )
    equations += center_infinity_equations
else:
    center_infinity = formal_center(a_at_infinity, ri, ri_inverse, 5)
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

# Exclude enhancements beyond the exact I7/I5/I4 orders.  Since
# A=-3*C^2 and B=2*C^3 to the imposed order, the first unmatched coefficient
# of B-2*C^3 is nonzero exactly when the corresponding leading discriminant
# coefficient is nonzero (the omitted scalar is a unit at these good primes).
center_cube0 = cube_coefficients(center0, 8)
center_cube1 = cube_coefficients(center1, 6)
center_cube_infinity = cube_coefficients(center_infinity, 5)
delta7_at_zero = b[7] - center_cube0[7]
delta5_at_one = sum(
    b[index] * binomial(index, 5) for index in range(5, 13)
) - center_cube1[5]
delta20_at_infinity = b[8] - center_cube_infinity[4]
open_factors += [
    ("exact_i7_order_inverse", delta7_at_zero),
    ("exact_i5_order_inverse", delta5_at_one),
    ("exact_i4_order_inverse", delta20_at_infinity),
]

if args.surface_only:
    surface_equations = [ring(equation) for equation in equations]
    surface_equations += [
        ring(delta7_at_zero * v["exact_i7_order_inverse"] - 1),
        ring(delta5_at_one * v["exact_i5_order_inverse"] - 1),
        ring(delta20_at_infinity * v["exact_i4_order_inverse"] - 1),
    ]
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
Y1 = sum(v[f"p1y{index}"] * t**index for index in range(1, 7))
section_equations = coefficient_equations(Y1**2 - rhs1)
equations += section_equations
q1_i5_open = X1(t=1) - r1 if args.q1_i5_identity_chart == "x" else Y1(t=1)
q1_i4_open = X1[4] - ri if args.q1_i4_identity_chart == "x" else Y1[6]
open_factors += [
    (section_open_inverse_names["q1"][0], v["p1y1"]),
    (section_open_inverse_names["q1"][1], q1_i5_open),
    (section_open_inverse_names["q1"][2], q1_i4_open),
]


def depth_two_all_node_section(prefix):
    # X=center0 mod t^2 enforces component +/-2 at I7.  Node incidence at
    # t=1 and infinity supplies the remaining two coefficients.
    x1 = center0[1]
    x2 = v[f"{prefix}x2"]
    x3 = r1 - 1 - ri - x1 - x2
    X = 1 + x1 * t + x2 * t**2 + x3 * t**3 + ri * t**4
    rhs = X**3 + A * X + B
    y2, y3, y4 = (v[f"{prefix}y{index}"] for index in range(2, 5))
    y5 = -y2 - y3 - y4
    Y = y2 * t**2 + y3 * t**3 + y4 * t**4 + y5 * t**5
    local_equations = coefficient_equations(Y**2 - rhs)
    return X, Y, local_equations


X2, Y2, p2_equations = depth_two_all_node_section("p2")
open_factors += [
    (section_open_inverse_names["q2"][0], v["p2y2"]),
    (section_open_inverse_names["q2"][1], Y2.derivative()(t=1)),
    (section_open_inverse_names["q2"][2], Y2[5]),
]
if args.basis_marking == "resolved-depth13":
    # Q3 has absolute profile (4,2,0).  Depth three at I7 fixes the first
    # three coefficients of X to the formal center.  Depth two at I5 fixes
    # its value and derivative there, hence the remaining two coefficients.
    q3_known_at_one = ring.one() + center0[1] + center0[2]
    q3_value_rhs = r1 - q3_known_at_one
    center1_derivative = center1[1]
    q3_derivative_rhs = center1_derivative - center0[1] - 2 * center0[2]
    q3x4 = q3_derivative_rhs - 3 * q3_value_rhs
    q3x3 = 4 * q3_value_rhs - q3_derivative_rhs
    X3 = 1 + center0[1] * t + center0[2] * t**2 + q3x3 * t**3 + q3x4 * t**4
    rhs3 = X3**3 + A * X3 + B
    # Keeping four sparse Y coefficients avoids the severe expression swell
    # caused by recursively substituting the surface-dependent y4,y5,y6.
    p3y3, p3y4 = v["p3y3"], v["p3y4"]
    p3y5 = -3 * p3y3 - 2 * p3y4
    p3y6 = 2 * p3y3 + p3y4
    Y3 = p3y3 * t**3 + p3y4 * t**4 + p3y5 * t**5 + p3y6 * t**6
    p3_equations = coefficient_equations(Y3**2 - rhs3)
    q3_i5_second_taylor = sum(
        binomial(index, 2) * Y3[index] for index in range(3, 7)
    )
    open_factors += [
        (section_open_inverse_names["q3"][0], v["p3y3"]),
        (section_open_inverse_names["q3"][1], q3_i5_second_taylor),
    ]
    q3_infinity_open = q3x4 - ri if args.q3_infinity_chart == "x" else Y3[6]
    open_factors.append((section_open_inverse_names["q3"][2], q3_infinity_open))
else:
    X3, Y3, p3_equations = depth_two_all_node_section("p3")
    open_factors += [
        (section_open_inverse_names["q3"][0], v["p3y2"]),
        (section_open_inverse_names["q3"][1], Y3.derivative()(t=1)),
        (section_open_inverse_names["q3"][2], Y3[5]),
    ]
equations += p2_equations + p3_equations

if not args.mw3_only:
    # Q4: simple pole at t=c and component depth one at all three reducible
    # fibres.  The three node values determine x0,x5,x6.
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
    # As for Q3, sparse coefficients keep the joint system cubic and avoid
    # expanding nested inverse substitutions in the surface parameters.
    p4y = [v[f"p4y{index}"] for index in range(1, 8)]
    p4y8 = -sum(p4y)
    Y4 = sum(p4y[index - 1] * t**index for index in range(1, 8)) + p4y8 * t**8
    p4_equations = coefficient_equations(Y4**2 - rhs4)
    open_factors += [
        (section_open_inverse_names["q4"][0], v["p4y1"]),
        (section_open_inverse_names["q4"][1], Y4.derivative()(t=1)),
        (section_open_inverse_names["q4"][2], Y4[8]),
        (section_open_inverse_names["q4"][3], X4(t=c)),
        (section_open_inverse_names["q4"][4], c),
        (section_open_inverse_names["q4"][5], c - 1),
    ]
    equations += p4_equations + [Y4(t=1)]

# Group only the small section-chart factors.  The three larger formal-branch
# exact-order factors stay separate; multiplying all opens together causes
# catastrophic expression swell.
open_groups = {}
for inverse_name, factor in open_factors:
    open_groups.setdefault(inverse_name, []).append(factor)
for inverse_name, factors in open_groups.items():
    product = ring.one()
    for factor in factors:
        product *= factor
    equations.append(product * v[inverse_name] - 1)

if args.fixed_rur_anchor:
    geometric_names = (
        [f"a{index}" for index in range(1, 7)] + ["r1", "ri"]
        + p1_names + p2_names + p3_names + p4_names
    )
    # Fixed integral coefficients are cycled only after the displayed list;
    # the ordered geometric-name list and this list together define the
    # prime-independent separator used for cross-prime common-producer work.
    anchor_coefficients = (
        2, -3, 5, 7, -11, 13, 17, -19, 23, 29, -31, 37,
        41, -43, 47, 53, -59, 61, 67, -71, 73, 79, -83, 89,
        97, -101, 103, 107, -109, 113, 127, -131, 137, 139, -149, 151,
    )
    if len(geometric_names) > len(anchor_coefficients):
        raise ArithmeticError("fixed RUR anchor coefficient list is too short")
    anchor_form = sum(
        field(anchor_coefficients[index]) * v[name]
        for index, name in enumerate(geometric_names)
    )
    equations.append(v["rur_anchor"] - anchor_form)

slice_items = list(args.slice)
if args.slice_a1 is not None:
    slice_items.append("a1={}".format(args.slice_a1))
slices = {}
for item in slice_items:
    if "=" not in item:
        raise SystemExit("--slice must have the form NAME=VALUE")
    name, raw_value = item.split("=", 1)
    if name not in v:
        raise SystemExit("unknown slice coordinate: {}".format(name))
    if name in slices:
        raise SystemExit("duplicate slice coordinate: {}".format(name))
    try:
        value = field(ZZ(raw_value))
    except (TypeError, ValueError):
        raise SystemExit("slice values must be integers reduced modulo p")
    slices[name] = value
    equations.append(v[name] - value)
hyperplanes = []
surface_slice_variables = [v[f"a{index}"] for index in range(1, 7)] + [r1, ri]
for item in args.surface_hyperplane:
    try:
        coefficients = [field(ZZ(entry)) for entry in item.split(",")]
    except (TypeError, ValueError):
        raise SystemExit("surface-hyperplane coefficients must be integers")
    if len(coefficients) != 9:
        raise SystemExit("--surface-hyperplane requires exactly nine coefficients")
    if not any(coefficients[:8]):
        raise SystemExit("surface hyperplane must involve a surface coordinate")
    equation = sum(
        coefficients[index] * surface_slice_variables[index] for index in range(8)
    ) + coefficients[8]
    equations.append(equation)
    hyperplanes.append(tuple(coefficients))
equations = [ring(equation) for equation in equations if equation != 0]
ideal = ring.ideal(equations)

print(
    f"NS0024JOINTRESOLVED|p={prime}|basis={args.basis_marking}"
    f"|stage={'MW3' if args.mw3_only else 'MW4'}|variables={ring.ngens()}"
    f"|equations={len(equations)}|q3_infinity_chart={args.q3_infinity_chart}"
    f"|q1_identity_charts={args.q1_i5_identity_chart},{args.q1_i4_identity_chart}"
    f"|slices={','.join('{}={}'.format(name, value) for name, value in slices.items())}"
    f"|hyperplanes={len(hyperplanes)}|fixed_rur_anchor={int(args.fixed_rur_anchor)}"
    f"|explicit_formal_centers={int(args.explicit_formal_centers)}"
    f"|split_section_opens={int(args.split_section_opens)}",
    flush=True,
)
if args.export_msolve is not None:
    output = args.export_msolve.resolve()
    output.write_text(
        ",".join(names) + "\n" + str(prime) + "\n"
        + ",\n".join(map(str, equations)) + "\n"
    )
    print(f"NS0024JOINTRESOLVEDMSOLVE|output={output}", flush=True)
if args.groebner:
    basis = ideal.groebner_basis(algorithm="libsingular:slimgb")
    print(
        f"NS0024JOINTRESOLVEDGB|p={prime}|basis={len(basis)}"
        f"|dimension={ideal.dimension()}|unit={int(ring.one() in ideal)}",
        flush=True,
    )
