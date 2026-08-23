#!/usr/bin/env sage -python
"""Probe I11 conditions through the formal double root, not Delta expansion."""

from sage.all import *

import argparse
import time


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--state", required=True)
parser.add_argument("--k", type=int, action="append", required=True)
parser.add_argument("--show", action="store_true")
parser.add_argument("--rank-k4-pivots", action="store_true")
parser.add_argument("--solve-k4-y3", action="store_true")
parser.add_argument("--save-state", default=None)
parser.add_argument("--save-equation-state", action="store_true")
args = parser.parse_args()

state = load(args.state)
ring = state["ring"]
fraction_field = ring.fraction_field()
d = ring.gens_dict()
substitutions = state["substitutions"]


def settle(value, passes=30):
    value = fraction_field(value)
    for unused in range(passes):
        old = value
        value = fraction_field(value.subs(substitutions))
        if value == old:
            break
    return value


series_ring = PowerSeriesRing(fraction_field, "u", default_prec=13)
u = series_ring.gen()
a_infinity = sum(fraction_field(d[f"a{8-index}"]) * u**index for index in range(9))
b_infinity = sum(fraction_field(d[f"b{12-index}"]) * u**index for index in range(13))

s_coefficients = [fraction_field(d["sinf"])]
for degree in range(1, 11):
    partial = sum(s_coefficients[index] * u**index for index in range(degree))
    residual_coefficient = (a_infinity + 3*partial**2)[degree]
    s_coefficients.append(
        fraction_field(-residual_coefficient / (6*s_coefficients[0]))
    )

formal_root = sum(
    s_coefficients[index] * u**index for index in range(len(s_coefficients))
)
root_residual = b_infinity - 2*formal_root**3
k4_solution = None
k4_reduced_value = None

for degree in args.k:
    started = time.time()
    value = settle(root_residual[degree])
    seconds = time.time() - started
    if value == 0:
        print(f"MW3I11ROOT|k={degree}|zero=1|seconds={seconds:.3f}", flush=True)
        continue
    numerator = ring(value.numerator())
    if degree == 4:
        k4_reduced_value = value
    variables = []
    linear = []
    for name in state["remaining"]:
        variable = ring(d[name])
        variable_degree = numerator.degree(variable)
        if variable_degree:
            variables.append(f"{name}:{variable_degree}")
            if variable_degree == 1 and numerator.derivative(variable).degree(variable) == 0:
                linear.append(name)
    print(
        f"MW3I11ROOT|k={degree}|seconds={seconds:.3f}|"
        f"degree={numerator.total_degree()}|terms={len(numerator.monomials())}|"
        f"linear={','.join(linear)}|vars={','.join(variables)}",
        flush=True,
    )
    if args.rank_k4_pivots and degree == 4:
        for name in ("a5", "y3", "y5"):
            variable = ring(d[name])
            coefficient = numerator.derivative(variable)
            constant = numerator.subs({variable: 0})
            print(
                f"MW3I11PIVOT|k=4|var={name}|"
                f"coefficient_degree={coefficient.total_degree()}|"
                f"coefficient_terms={len(coefficient.monomials())}|"
                f"constant_degree={constant.total_degree()}|"
                f"constant_terms={len(constant.monomials())}",
                flush=True,
            )
    if args.solve_k4_y3 and degree == 4:
        variable = ring(d["y3"])
        coefficient = numerator.derivative(variable)
        constant = numerator.subs({variable: 0})
        assert numerator == coefficient*variable + constant
        k4_solution = fraction_field(-constant) / fraction_field(coefficient)
        print(
            f"MW3I11SOLVE|k=4|var=y3|"
            f"num_terms={len(constant.monomials())}|"
            f"den_terms={len(coefficient.monomials())}|verified=1",
            flush=True,
        )
    if args.show:
        print(f"MW3I11ROOT_EXPR|k={degree}|expr={numerator}", flush=True)

if args.solve_k4_y3:
    if k4_solution is None:
        raise SystemExit("--solve-k4-y3 requires --k=4")
    substitutions[fraction_field(d["y3"])] = k4_solution
    remaining = list(state["remaining"])
    remaining.remove("y3")
    equations = dict(state["equations"])
    equations["I11_root_R4"] = fraction_field(root_residual[4])
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
        print(f"MW3I11ROOT|state_saved={args.save_state}", flush=True)

if args.save_equation_state:
    if k4_reduced_value is None or not args.save_state:
        raise SystemExit("--save-equation-state requires --k=4 and --save-state")
    equations = dict(state["equations"])
    equations["I11_root_R4_reduced"] = k4_reduced_value
    save(
        {
            "ring": ring,
            "equations": equations,
            "substitutions": substitutions,
            "remaining": list(state["remaining"]),
            "prime": state["prime"],
        },
        args.save_state,
    )
    print(f"MW3I11ROOT|equation_state_saved={args.save_state}", flush=True)

print("MW3I11ROOT|status=PASS", flush=True)
