#!/usr/bin/env sage -python
"""Rank five-equation I2 pivots after extending a two-norm MW3 state."""

from sage.all import *

import argparse
import time


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--state", required=True)
parser.add_argument("--solve-omit-a5", action="store_true")
parser.add_argument("--solve-block-omit-a5", action="store_true")
parser.add_argument("--save-state", default=None)
args = parser.parse_args()

started = time.time()
state = load(args.state)
old_ring = state["ring"]
field = old_ring.base_ring()
old_names = list(old_ring.variable_names())
if "r2" in old_names:
    raise SystemExit("state already contains r2")

ring = PolynomialRing(field, old_names + ["r2"], order="degrevlex")
d = ring.gens_dict()
fraction_field = ring.fraction_field()
embedding = old_ring.hom([ring(d[name]) for name in old_names], ring)


def transport(value):
    value = old_ring.fraction_field()(value)
    return fraction_field(embedding(value.numerator())) / fraction_field(
        embedding(value.denominator())
    )


equations = {tag: transport(value) for tag, value in state["equations"].items()}
substitutions = {
    fraction_field(embedding(old_ring(key))): transport(value)
    for key, value in state["substitutions"].items()
}


def settle(value, passes=30):
    value = fraction_field(value)
    for unused in range(passes):
        old = value
        value = fraction_field(value.subs(substitutions))
        if value == old:
            break
    return value


t_ring = PolynomialRing(fraction_field, "t")
t = t_ring.gen()
A = sum(fraction_field(d[f"a{i}"]) * t**i for i in range(9))
X = (
    fraction_field(d["s0"])
    + fraction_field(d["x1"])*t
    + fraction_field(d["x2"])*t**2
    + fraction_field(d["x3"])*t**3
    + fraction_field(d["sinf"])*t**4
)
Y = (
    fraction_field(d["y1"])*t
    + fraction_field(d["y2"])*t**2
    + fraction_field(d["y3"])*t**3
    + fraction_field(d["y4"])*t**4
    + fraction_field(d["y5"])*t**5
)
r2 = fraction_field(d["r2"])
s1 = fraction_field(d["s1"])

components = {
    "I2_1_A": settle(equations["I2_1_A"]),
    "I2_lam_A": settle(equations["I2_lam_A"]),
    "point_X": settle(X(1) - (r2**2 - 2*s1)),
    "point_Y": settle(Y(1) - r2*(r2**2 - 3*s1)),
    "A_prime": settle(A.derivative(t)(1)),
    "X_prime": settle(X.derivative(t)(1)),
    "Y_prime": settle(Y.derivative(t)(1)),
}
if "a6" in state["remaining"]:
    components["P1_9"] = settle(equations["P1_9"])
components["tate"] = fraction_field(
    2*r2*components["Y_prime"]
    - 3*(r2**2-s1)*components["X_prime"]
    - components["A_prime"]
)

for name, value in components.items():
    numerator = ring(value.numerator())
    print(
        f"MW3S4COMP|name={name}|degree={numerator.total_degree()}|"
        f"terms={len(numerator.monomials())}",
        flush=True,
    )

if "a6" in state["remaining"]:
    a_four_equations = [
        ring(components[name].numerator())
        for name in ("I2_1_A", "I2_lam_A", "tate", "P1_9")
    ]
    a_four_variables = [
        ring(d[name]) for name in ("a3", "a4", "a5", "a6", "b2")
    ]
    for omitted in a_four_variables:
        pivot_variables = [
            variable for variable in a_four_variables if variable != omitted
        ]
        coefficient_matrix = matrix(
            ring,
            [
                [equation.derivative(variable) for variable in pivot_variables]
                for equation in a_four_equations
            ],
        )
        nonlinear_dependencies = sorted(
            {
                str(variable)
                for coefficient in coefficient_matrix.list()
                for variable in pivot_variables
                if coefficient.degree(variable) != 0
            }
        )
        if nonlinear_dependencies:
            print(
                f"MW3S4AFOUR|omit={omitted}|nonlinear="
                f"{','.join(nonlinear_dependencies)}",
                flush=True,
            )
            continue
        constants = vector(
            ring,
            [
                equation.subs({variable: 0 for variable in pivot_variables})
                for equation in a_four_equations
            ],
        )
        determinant = coefficient_matrix.det()
        numerator_terms = []
        numerator_degrees = []
        for column in range(4):
            replaced = matrix(ring, coefficient_matrix)
            replaced.set_column(column, -constants)
            numerator = replaced.det()
            numerator_terms.append(len(numerator.monomials()))
            numerator_degrees.append(numerator.total_degree())
        print(
            f"MW3S4AFOUR|omit={omitted}|"
            f"det_degree={determinant.total_degree()}|"
            f"det_terms={len(determinant.monomials())}|"
            f"num_degrees={','.join(map(str,numerator_degrees))}|"
            f"num_terms={','.join(map(str,numerator_terms))}",
            flush=True,
        )

five_equations = [
    ring(components[name].numerator())
    for name in ("I2_1_A", "I2_lam_A", "point_X", "point_Y", "tate")
]
six_variables = [ring(d[name]) for name in ("a3", "a4", "a5", "b2", "x2", "y2")]

a_equations = [
    ring(components[name].numerator())
    for name in ("I2_1_A", "I2_lam_A", "tate")
]
a_variables = [ring(d[name]) for name in ("a3", "a4", "a5", "b2")]
a_solve_data = None
for omitted in a_variables:
    pivot_variables = [variable for variable in a_variables if variable != omitted]
    coefficient_matrix = matrix(
        ring,
        [
            [equation.derivative(variable) for variable in pivot_variables]
            for equation in a_equations
        ],
    )
    # These three equations are genuinely affine-linear as a block in the A
    # variables; point-coordinate bilinearity is deliberately left outside.
    for coefficient in coefficient_matrix.list():
        assert all(coefficient.degree(variable) == 0 for variable in pivot_variables)
    constants = vector(
        ring,
        [
            equation.subs({variable: 0 for variable in pivot_variables})
            for equation in a_equations
        ],
    )
    determinant = coefficient_matrix.det()
    numerator_terms = []
    numerator_degrees = []
    for column in range(3):
        replaced = matrix(ring, coefficient_matrix)
        replaced.set_column(column, -constants)
        numerator = replaced.det()
        numerator_terms.append(len(numerator.monomials()))
        numerator_degrees.append(numerator.total_degree())
    print(
        f"MW3S4ATHREE|omit={omitted}|det_degree={determinant.total_degree()}|"
        f"det_terms={len(determinant.monomials())}|"
        f"num_degrees={','.join(map(str,numerator_degrees))}|"
        f"num_terms={','.join(map(str,numerator_terms))}",
        flush=True,
    )
    if str(omitted) == "a5":
        a_solve_data = (
            pivot_variables,
            coefficient_matrix,
            constants,
            determinant,
        )

if args.solve_block_omit_a5:
    pivot_variables, coefficient_matrix, constants, determinant = a_solve_data
    solve_numerators = []
    for column in range(3):
        replaced = matrix(ring, coefficient_matrix)
        replaced.set_column(column, -constants)
        solve_numerators.append(replaced.det())
    solve_values = [
        fraction_field(numerator) / fraction_field(determinant)
        for numerator in solve_numerators
    ]
    simultaneous = {
        fraction_field(variable): value
        for variable, value in zip(pivot_variables, solve_values)
    }
    for equation in a_equations:
        if fraction_field(equation).subs(simultaneous) != 0:
            raise ArithmeticError("three-equation A block solve failed")

    point_substitutions = {}
    for variable_name, component_name in (
        ("x2", "point_X"),
        ("y2", "point_Y"),
    ):
        variable = ring(d[variable_name])
        equation = ring(components[component_name].numerator())
        coefficient = equation.derivative(variable)
        assert coefficient.degree(variable) == 0
        root = fraction_field(-equation.subs({variable: 0})) / fraction_field(
            coefficient
        )
        if fraction_field(equation).subs({fraction_field(variable): root}) != 0:
            raise ArithmeticError(f"{component_name} linear solve failed")
        point_substitutions[fraction_field(variable)] = root

    substitutions.update(simultaneous)
    substitutions.update(point_substitutions)
    remaining = list(state["remaining"]) + ["r2"]
    solved_variables = list(pivot_variables) + [ring(d["x2"]), ring(d["y2"])]
    for variable in solved_variables:
        remaining.remove(str(variable))
    equations.update(
        {
            "I2_1_point_X": components["point_X"],
            "I2_1_point_Y": components["point_Y"],
            "I2_1_tate_compact": components["tate"],
        }
    )
    print(
        f"MW3S4BLOCKSOLVE|omit=a5|vars="
        f"{','.join(map(str,solved_variables))}|"
        f"remaining={','.join(remaining)}|verified=1",
        flush=True,
    )
    if args.save_state:
        save(
            {
                "ring": ring,
                "equations": equations,
                "substitutions": substitutions,
                "remaining": remaining,
                "prime": state["prime"],
            },
            args.save_state,
        )
        print(f"MW3S4|state_saved={args.save_state}", flush=True)

solve_data = None
for omitted in six_variables:
    pivot_variables = [variable for variable in six_variables if variable != omitted]
    coefficient_matrix = matrix(
        ring,
        [
            [equation.derivative(variable) for variable in pivot_variables]
            for equation in five_equations
        ],
    )
    constants = vector(
        ring,
        [
            equation.subs({variable: 0 for variable in pivot_variables})
            for equation in five_equations
        ],
    )
    determinant = coefficient_matrix.det()
    numerator_terms = []
    numerator_degrees = []
    for column in range(5):
        replaced = matrix(ring, coefficient_matrix)
        replaced.set_column(column, -constants)
        numerator = replaced.det()
        numerator_terms.append(len(numerator.monomials()))
        numerator_degrees.append(numerator.total_degree())
    print(
        f"MW3S4FIVE|omit={omitted}|det_degree={determinant.total_degree()}|"
        f"det_terms={len(determinant.monomials())}|"
        f"num_degrees={','.join(map(str,numerator_degrees))}|"
        f"num_terms={','.join(map(str,numerator_terms))}",
        flush=True,
    )
    if str(omitted) == "a5":
        solve_data = (
            pivot_variables,
            coefficient_matrix,
            constants,
            determinant,
        )

if args.solve_omit_a5:
    pivot_variables, coefficient_matrix, constants, determinant = solve_data
    solve_numerators = []
    for column in range(5):
        replaced = matrix(ring, coefficient_matrix)
        replaced.set_column(column, -constants)
        solve_numerators.append(replaced.det())
    solve_values = [
        fraction_field(numerator) / fraction_field(determinant)
        for numerator in solve_numerators
    ]
    simultaneous = {
        fraction_field(variable): value
        for variable, value in zip(pivot_variables, solve_values)
    }
    for equation in five_equations:
        if fraction_field(equation).subs(simultaneous) != 0:
            raise ArithmeticError("five-equation simultaneous solve failed")
    substitutions.update(simultaneous)
    remaining = list(state["remaining"]) + ["r2"]
    for variable in pivot_variables:
        remaining.remove(str(variable))
    equations.update(
        {
            "I2_1_point_X": components["point_X"],
            "I2_1_point_Y": components["point_Y"],
            "I2_1_tate_compact": components["tate"],
        }
    )
    print(
        f"MW3S4SOLVE|omit=a5|vars={','.join(map(str,pivot_variables))}|"
        f"remaining={','.join(remaining)}|verified=1",
        flush=True,
    )
    if args.save_state:
        save(
            {
                "ring": ring,
                "equations": equations,
                "substitutions": substitutions,
                "remaining": remaining,
                "prime": state["prime"],
            },
            args.save_state,
        )
        print(f"MW3S4|state_saved={args.save_state}", flush=True)

print(f"MW3S4|seconds={time.time()-started:.3f}|status=PASS_RANKED", flush=True)
