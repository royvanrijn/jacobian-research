#!/usr/bin/env python3
"""Exact unit certificate on a four-parameter SIC(2) anchor boundary.

The family lies on the normalized s0=1 chart with A=B=mu_2=0.  Exact
contraction constructs moments three through seven.  Singular then proves
that these five restricted moments generate the unit ideal over Q.
"""

from __future__ import annotations

from math import factorial
import shutil
import subprocess

import sympy as sp


PARAMETERS = (
    "s0", "s1", "s2", "s3", "s4", "s5", "s6",
    "t0", "t1", "t2", "t3", "t4",
)
SEXTIC_MAP = (
    (0, 0, 3, -1), (0, 1, 4, -3), (0, 2, 5, -3), (0, 3, 6, -1),
    (1, 0, 2, 3), (1, 1, 3, 9), (1, 2, 4, 9), (1, 3, 5, 3),
    (2, 0, 1, -3), (2, 1, 2, -9), (2, 2, 3, -9), (2, 3, 4, -3),
    (3, 0, 0, 1), (3, 1, 1, 3), (3, 2, 2, 3), (3, 3, 3, 1),
)
QUARTIC_MAP_LOCAL = (
    (0, 0, 2, 1), (0, 1, 3, 2), (0, 2, 4, 1),
    (1, 0, 1, -2), (1, 1, 2, -3), (1, 3, 4, 1),
    (2, 0, 0, 1), (2, 2, 2, -3), (2, 3, 3, -2),
    (3, 1, 0, 1), (3, 2, 1, 2), (3, 3, 2, 1),
)
QUARTIC_MAP = tuple(
    (i, j, parameter + 7, coefficient)
    for i, j, parameter, coefficient in QUARTIC_MAP_LOCAL
)
NORMALIZED_QUADRATIC = (
    (0, 0, -1), (1, 1, -1), (2, 2, 1), (3, 3, 1),
)

a, b, h, q, x, y = sp.symbols("a b h q x y")
values: dict[str, sp.Expr] = {
    parameter: sp.Integer(0) for parameter in PARAMETERS
}
values.update({
    "s0": sp.Integer(1),
    "s4": -4 * q**2,
    "s5": h,
    "s6": (14 * a * b - 168 * a * q + 70) / 3,
    "t0": a,
    "t1": q,
    "t3": 3 * a,
    "t4": b,
})

# Verify the two reduced mu_3 pivot coefficients A and B, and mu_2.
coefficient_a = -3 * values["t0"] + values["t3"]
coefficient_b = -3 * values["s4"] - 12 * values["t1"] ** 2
normalized_mu2 = (
    -3 * values["s0"] * values["s6"]
    + 18 * values["s1"] * values["s5"]
    - 45 * values["s2"] * values["s4"]
    + 30 * values["s3"] ** 2
    + 14 * values["t0"] * values["t4"]
    - 56 * values["t1"] * values["t3"]
    + 42 * values["t2"] ** 2
    + 70
)
assert coefficient_a == 0
assert coefficient_b == 0
assert sp.expand(normalized_mu2) == 0

grid: list[list[sp.Expr]] = [
    [sp.Integer(0)] * 4 for _ in range(4)
]
for i, j, parameter, coefficient in SEXTIC_MAP + QUARTIC_MAP:
    grid[i][j] += coefficient * values[PARAMETERS[parameter]]
for i, j, coefficient in NORMALIZED_QUADRATIC:
    grid[i][j] += coefficient

polynomial = sp.Poly(sum(
    grid[i][j] * x**i * y**j
    for i in range(4)
    for j in range(4)
), x, y)
power = sp.Poly(1, x, y)
moments: list[sp.Expr] = []
for order in range(1, 8):
    power *= polynomial
    if order < 3:
        continue
    moment = sum(
        factorial(3 * order - diagonal)
        * factorial(diagonal)
        * power.coeff_monomial(x**diagonal * y**diagonal)
        for diagonal in range(3 * order + 1)
    )
    primitive = sp.Poly(moment, a, b, h, q).primitive()[1]
    moments.append(sp.expand(primitive.as_expr()))

assert sp.factor(moments[0]) == 3 * a**3 + 4 * a * q**4 + 8 * q**3


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "Singular is required"
generators = [singular_expression(moment) for moment in moments]
program = f"""
ring r=0,(a,b,h,q),dp;
option(redSB);
ideal I34={generators[0]},{generators[1]};
ideal G34=std(I34);
if (dim(G34)!=2)
{{
  print("DIMENSION_FAILURE_34");
  exit(1);
}}
ideal I345=G34,{generators[2]};
ideal G345=std(I345);
if (dim(G345)!=1)
{{
  print("DIMENSION_FAILURE_345");
  exit(1);
}}
ideal I3456=G345,{generators[3]};
ideal G3456=std(I3456);
if (dim(G3456)!=0 || vdim(G3456)!=372)
{{
  print("FINITE_QUOTIENT_FAILURE");
  exit(1);
}}
ideal I34567=G3456,{generators[4]};
ideal G34567=std(I34567);
if (G34567[1]!=1)
{{
  print("UNIT_FAILURE");
  exit(1);
}}
print("CERTIFICATE 2 1 0 372 UNIT");
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=120,
)
assert completed.stdout.strip() == "CERTIFICATE 2 1 0 372 UNIT", (
    completed.stdout
)

print("PASS: the four-parameter family lies on A=B=mu_2=0")
print("PASS: moments 3 through 6 leave an exact quotient of length 372")
print("PASS: moment 7 makes the characteristic-zero ideal the unit ideal")
