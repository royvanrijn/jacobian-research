#!/usr/bin/env sage -python
"""Extend an MW3 P1 state by the P1_10 norm parameter r1."""

from sage.all import *

import argparse
import time


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--state", required=True)
parser.add_argument("--save-state", required=True)
args = parser.parse_args()

started = time.time()
state = load(args.state)
old_ring = state["ring"]
field = old_ring.base_ring()
old_names = list(old_ring.variable_names())
if "r1" in old_names:
    raise SystemExit("state already contains r1")

ring = PolynomialRing(field, old_names + ["r1"], order="degrevlex")
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


lam = fraction_field(d["lam"])
r0 = fraction_field(d["r0"])
r1 = fraction_field(d["r1"])
x1 = fraction_field(d["x1"])
x2 = fraction_field(d["x2"])
x3 = fraction_field(d["x3"])
y5 = fraction_field(d["y5"])
sl = fraction_field(d["sl"])
b11 = ring(d["b11"])

substitutions[sl] = lam**3*x3 + lam**2*x2 + lam*x1 + (r0**2-r1**2)/3
p110 = ring(settle(equations["P1_10"]).numerator())
assert p110.degree(b11) == 2
quadratic_a = ring(p110.derivative(b11, 2) / 2)
quadratic_b = ring(p110.derivative(b11).subs({b11: 0}))
quadratic_c = ring(p110.subs({b11: 0}))
discriminant = ring(quadratic_b**2 - 4*quadratic_a*quadratic_c)
expected_square = ring(y5**2 * lam**20 * r1**6)
quotient, remainder = discriminant.quo_rem(expected_square)
assert remainder == 0 and quotient in field and field(quotient).is_square()
square_unit = field(quotient)
square_root_unit = square_unit.sqrt()
quadratic_root = fraction_field(
    -quadratic_b + square_root_unit*y5*lam**10*r1**3
) / fraction_field(2*quadratic_a)
substitutions[fraction_field(b11)] = quadratic_root

remaining = list(state["remaining"]) + ["r1"]
remaining.remove("sl")
remaining.remove("b11")

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

print(
    f"MW3P110EXTEND|square_unit={square_unit}|"
    f"p110_degree={p110.total_degree()}|p110_terms={len(p110.monomials())}|"
    f"remaining={','.join(remaining)}|"
    f"seconds={time.time()-started:.3f}|status=PASS",
    flush=True,
)
print(f"MW3P110EXTEND|state_saved={args.save_state}", flush=True)
