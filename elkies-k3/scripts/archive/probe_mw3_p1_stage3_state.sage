#!/usr/bin/env sage -python
"""Probe one residual equation from a serialized MW3 P1 reduction state."""

from sage.all import *

import argparse
import time


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--state", required=True)
parser.add_argument("--tag", action="append", required=True)
parser.add_argument("--show", action="store_true")
parser.add_argument("--factor", action="store_true")
args = parser.parse_args()

state = load(args.state)
R = state["ring"]
RF = R.fraction_field()
equations = state["equations"]
substitutions = state["substitutions"]
remaining = state["remaining"]
d = R.gens_dict()


def settle(value, passes=30):
    value = RF(value)
    for unused in range(passes):
        old = value
        value = RF(value.subs(substitutions))
        if value == old:
            break
    return value


for tag in args.tag:
    started = time.time()
    if tag in (
        "I2_1_tate_compact",
        "I2_1_A_prime",
        "I2_1_X_prime",
        "I2_1_Y_prime",
    ):
        Rt = PolynomialRing(RF, "t")
        t = Rt.gen()
        A = sum(RF(d[f"a{i}"]) * t**i for i in range(9))
        X = (
            RF(d["s0"])
            + RF(d["x1"])*t
            + RF(d["x2"])*t**2
            + RF(d["x3"])*t**3
            + RF(d["sinf"])*t**4
        )
        Y = (
            RF(d["y1"])*t
            + RF(d["y2"])*t**2
            + RF(d["y3"])*t**3
            + RF(d["y4"])*t**4
            + RF(d["y5"])*t**5
        )
        s1 = RF(d["s1"])
        if tag == "I2_1_A_prime":
            raw_value = A.derivative(t)(1)
        elif tag == "I2_1_X_prime":
            raw_value = X.derivative(t)(1)
        elif tag == "I2_1_Y_prime":
            raw_value = Y.derivative(t)(1)
        else:
            r2 = RF(d["r2"])
            raw_value = (
                2*r2*Y.derivative(t)(1)
                - 3*(r2**2-s1)*X.derivative(t)(1)
                - A.derivative(t)(1)
            )
    else:
        raw_value = equations[tag]
    value = settle(raw_value)
    seconds = time.time() - started
    if value == 0:
        print(f"MW3S3PROBE|tag={tag}|zero=1|seconds={seconds:.3f}", flush=True)
        continue
    numerator = R(value.numerator())
    linear = []
    variables = []
    for name in remaining:
        variable = R(d[name])
        degree = numerator.degree(variable)
        if degree:
            variables.append(f"{name}:{degree}")
            if degree == 1 and numerator.derivative(variable).degree(variable) == 0:
                linear.append(name)
    print(
        f"MW3S3PROBE|tag={tag}|seconds={seconds:.3f}|"
        f"degree={numerator.total_degree()}|terms={len(numerator.monomials())}|"
        f"linear={','.join(linear)}|vars={','.join(variables)}",
        flush=True,
    )
    if args.factor:
        factorization = numerator.factor()
        print(
            f"MW3S3FACTOR|tag={tag}|unit={factorization.unit()}|factors="
            + ";".join(
                f"degree:{factor.total_degree()},terms:{len(factor.monomials())},"
                f"exp:{exponent}"
                for factor, exponent in factorization
            ),
            flush=True,
        )
    if args.show:
        print(f"MW3S3EXPR|tag={tag}|expr={numerator}", flush=True)

print("MW3S3|done", flush=True)
