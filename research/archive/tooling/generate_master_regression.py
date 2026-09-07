#!/usr/bin/env python3
"""Print exact data for a new (m,r) master-cancellation regression."""
import argparse

import sympy as sp

from master_cancellation import hensel_jet, parameter_polynomial


parser = argparse.ArgumentParser()
parser.add_argument("m", type=int)
parser.add_argument("r", type=int)
args = parser.parse_args()
if args.m < 1 or args.r < 1:
    parser.error("m and r must be positive")

A, q = sp.symbols("A q")
M = parameter_polynomial(args.m, args.r, q)
h = hensel_jet(args.m, args.r, A, q)

print(f"m = {args.m}")
print(f"r = {args.r}")
print(f"generic fiber degree = {args.r * (args.m + 1) + 1}")
print("M_m_r(q) =", M)
print("h(A) =", h)
