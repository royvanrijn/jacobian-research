#!/usr/bin/env python3
"""Exact anchor fiber for the generic squarefree-quartic GVC theorem.

For A=u*v*(u-v)*(u-lambda*v), the four expected annihilator sections are
x^4, y^4, (x+y)^4, and (lambda*x+y)^4.  The checker verifies them
symbolically and proves that at lambda=2 the first six moments cut out
exactly these four reduced projective points.
"""

from __future__ import annotations

from math import factorial
import shutil
import subprocess

import sympy as sp


x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e, lam = sp.symbols("a b c d e lam")
P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4
A = u * v * (u - v) * (u - lam * v)


def apolar_moment(order: int) -> sp.Expr:
    symbol_power = sp.Poly(sp.expand(A**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    return sp.expand(sum(
        coefficient
        * polynomial_power.coeff_monomial(x**x_order * y**y_order)
        * factorial(x_order)
        * factorial(y_order)
        for (x_order, y_order), coefficient in symbol_power.terms()
    ))


moments = [apolar_moment(order) for order in range(1, 7)]
expected_sections = (
    {a: 1, b: 0, c: 0, d: 0, e: 0},
    {a: 0, b: 0, c: 0, d: 0, e: 1},
    {a: 1, b: 4, c: 6, d: 4, e: 1},
    {
        a: lam**4,
        b: 4 * lam**3,
        c: 6 * lam**2,
        d: 4 * lam,
        e: 1,
    },
)
for section in expected_sections:
    for moment in moments:
        assert sp.expand(moment.subs(section)) == 0

# At lambda=2, mu_1=0 gives b=2*c-2*d.  Work in the remaining homogeneous
# coordinate ring Q[a,c,d,e].
specialized = [sp.expand(moment.subs(lam, 2)) for moment in moments]
assert sp.expand(specialized[0] - 6 * (b - 2 * c + 2 * d)) == 0
restricted: list[sp.Expr] = []
for moment in specialized[1:]:
    numerator = sp.fraction(
        sp.cancel(moment.subs(b, 2 * c - 2 * d))
    )[0]
    primitive = sp.Poly(numerator, a, c, d, e).primitive()[1]
    restricted.append(sp.expand(primitive.as_expr()))


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "Singular is required"
ideal_generators = ",".join(
    singular_expression(moment) for moment in restricted
)
program = f"""
ring r=0,(a,c,d,e),dp;
option(redSB);
ideal I={ideal_generators};
ideal GI=std(I);

// x^4, y^4, (x+y)^4, and (2*x+y)^4 after eliminating b.
ideal J1=c,d,e;
ideal J2=a,c,d;
ideal J3=c-6a,d-4a,e-a;
ideal J4=2c-3a,2d-a,16e-a;
ideal J=intersect(intersect(J1,J2),intersect(J3,J4));
ideal GJ=std(J);

if (dim(GI)!=1 || mult(GI)!=4)
{{
  print("SPECIAL_FIBER_DEGREE_FAILURE");
  exit(1);
}}
if (dim(GJ)!=1 || mult(GJ)!=4 || size(GJ)!=6)
{{
  print("EXPECTED_IDEAL_FAILURE");
  exit(1);
}}
if (size(reduce(I,GJ))!=0)
{{
  print("CONTAINMENT_FAILURE");
  exit(1);
}}
int generator;
for (generator=1;generator<=size(GJ);generator++)
{{
  if (reduce(GJ[generator]^5,GI)!=0)
  {{
    print("POWER_FAILURE "+string(generator));
    exit(1);
  }}
}}
print("CERTIFICATE 4 6 5");
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
)
assert completed.stdout.strip() == "CERTIFICATE 4 6 5", completed.stdout

print("PASS: four annihilator sections vanish for symbolic lambda")
print("PASS: at lambda=2 the projective moment fiber has degree four")
print("PASS: six radical generators have exact fifth-power certificates")
