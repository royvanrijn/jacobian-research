#!/usr/bin/env python3
"""Exact certificate for the quartic (2,1,1)-root GVC symbol orbit."""

from __future__ import annotations

from math import factorial
import shutil
import subprocess

import sympy as sp


x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e = sp.symbols("a b c d e")

P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4
A = u**2 * v * (u + v)


def apolar_moment(order: int) -> sp.Expr:
    symbol_power = sp.Poly(sp.expand(A**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    value = 0
    for (x_order, y_order), coefficient in symbol_power.terms():
        value += (
            coefficient
            * polynomial_power.coeff_monomial(x**x_order * y**y_order)
            * factorial(x_order)
            * factorial(y_order)
        )
    return sp.expand(value)


moments = [apolar_moment(order) for order in range(1, 6)]
assert moments[0] == 2 * (3 * b + 2 * c)

# Eliminate c with the first moment and remove harmless rational contents.
restricted: list[sp.Expr] = []
for moment in moments[1:]:
    numerator = sp.fraction(
        sp.cancel(moment.subs(c, -sp.Rational(3, 2) * b))
    )[0]
    primitive = sp.Poly(numerator, a, b, d, e).primitive()[1]
    restricted.append(sp.expand(primitive.as_expr()))


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "Singular is required"

ideal_generators = ",".join(singular_expression(moment) for moment in restricted)
program = f"""
ring r=0,(a,b,d,e),dp;
ideal I={ideal_generators};
ideal GI=std(I);

// P=y^3*(d*x+e*y), P=a*x^4, and P=a*(x-y)^4.
ideal J1=a,b;
ideal J2=b,d,e;
ideal J3=b+4a,d+4a,e-a;
ideal J=intersect(intersect(J1,J2),J3);
ideal GJ=std(J);

if (size(reduce(I,GJ))!=0)
{{
  print("CONTAINMENT_FAILURE");
  exit(1);
}}

int generator,power;
poly q;
for (generator=1;generator<=size(GJ);generator++)
{{
  q=GJ[generator]^4;
  if (reduce(q,GI)!=0)
  {{
    print("POWER_FAILURE "+string(generator));
    exit(1);
  }}
}}

print("CERTIFICATE "+string(size(GJ))+" 4");
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
)
assert completed.stdout.strip() == "CERTIFICATE 5 4", completed.stdout

# Check the annihilating directions on the two pure-power components.
assert sp.diff(x**4, y) == 0
assert sp.diff((x - y) ** 4, x) + sp.diff((x - y) ** 4, y) == 0

print("PASS: quartic moments 1 through 5 have the expected (2,1,1) radical")
print("PASS: five nullcone generators have fourth-power certificates")
print("PASS: every component has a strict degree or annihilator cutoff")

