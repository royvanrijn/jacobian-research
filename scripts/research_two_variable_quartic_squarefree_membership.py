#!/usr/bin/env python3
"""Research the remaining squarefree quartic Rabinowitsch membership.

This script reconstructs the ideal

  I=(f3,f4,f5,f6,z*p*(8*c-3*d^2)-1)

over Q and asks Singular's modular standard-basis routine for a candidate
basis.  It then checks both ideal containments exactly:

* every generator of I reduces to zero by the candidate basis; and
* ``lift(I,G)`` expresses every candidate generator as an exact
  Q[z,c,d,lambda]-linear combination of the original generators.

Only the second check turns modular reconstruction into a characteristic-zero
membership certificate.
"""

from __future__ import annotations

from math import factorial
from pathlib import Path
import shutil
import subprocess
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e, lam = sp.symbols("a b c d e lam")
P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4
A = u * v * (u - v) * (u - lam * v)


def apolar_moment(order: int) -> sp.Expr:
    symbol_power = sp.Poly(sp.expand(A**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    return sp.expand(
        sum(
            coefficient
            * polynomial_power.coeff_monomial(x**x_order * y**y_order)
            * factorial(x_order)
            * factorial(y_order)
            for (x_order, y_order), coefficient in symbol_power.terms()
        )
    )


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def main() -> None:
    moments = [apolar_moment(order) for order in range(1, 7)]
    b_solution = sp.Rational(2, 3) * c * (lam + 1) - d * lam
    chart_moments = [
        sp.cancel(moment.subs(b, b_solution).subs(e, 1))
        for moment in moments[1:]
    ]
    second = sp.Poly(chart_moments[0], a)
    pivot = sp.factor(second.coeff_monomial(a) / 576)
    constant = second.coeff_monomial(1)
    a_solution = sp.cancel(-constant / second.coeff_monomial(a))
    eliminated_polynomials = [
        sp.Poly(
            sp.fraction(sp.cancel(moment.subs(a, a_solution)))[0],
            c,
            d,
            domain=sp.QQ[lam],
        ).primitive()[1].as_expr()
        for moment in chart_moments[1:]
    ]
    h = 8 * c - 3 * d**2
    target = lam**4 * (lam - 1) ** 4

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    if len(sys.argv) == 2 and sys.argv[1] == "--saturation":
        program = """
LIB "modstd.lib";
option(redSB);
ring r=0,(c,d,lam),dp;
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)};
poly multiplier=p*h;
ideal F=f3,f4,f5,f6;
timer=1;
ideal G=modStd(F);
int modular_seconds=timer;
ideal V=std(G);
int exponent=-1;
poly powered=target;
int candidate;
for (candidate=0;candidate<=32;candidate++)
{{
  if (reduce(powered,V)==0)
  {{
    exponent=candidate;
    break;
  }}
  powered=powered*multiplier;
}}
if (exponent<0)
{{
  print("NO_SATURATION_EXPONENT_THROUGH_32");
  exit(1);
}}
ideal T=powered;
timer=1;
matrix L=lift(F,T);
int lift_seconds=timer;
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{{
  reconstructed=reconstructed+F[row]*L[row,1];
}}
if (reconstructed-powered!=0)
{{
  print("BAD_TARGET_ONLY_LIFT");
  exit(1);
}}
print("EXACT_TARGET_SATURATION_CERTIFICATE");
print("exponent",exponent);
print("modular_seconds",modular_seconds);
print("lift_seconds",lift_seconds);
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=1200,
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        return

    program = """
LIB "modstd.lib";
option(redSB);
ring r=0,(z,c,d,lam),lp;
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        program += f"poly f{order}={singular_expression(polynomial)};\n"
    program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)};
ideal I=f3,f4,f5,f6,z*p*h-1;
timer=1;
ideal G=modStd(I);
int modular_seconds=timer;
ideal V=std(G);
if (size(reduce(I,V))!=0)
{{
  print("BAD_FORWARD_CONTAINMENT");
  exit(1);
}}
if (reduce(target,V)!=0)
{{
  print("BAD_TARGET_REDUCTION");
  exit(1);
}}
timer=1;
matrix T=lift(I,G);
int lift_seconds=timer;
int row;
int column;
poly reconstructed;
for (column=1;column<=size(G);column++)
{{
  reconstructed=0;
  for (row=1;row<=size(I);row++)
  {{
    reconstructed=reconstructed+I[row]*T[row,column];
  }}
  if (reconstructed-G[column]!=0)
  {{
    print("BAD_REVERSE_CONTAINMENT");
    exit(1);
  }}
}}
print("EXACT_BIDIRECTIONAL_CERTIFICATE");
print("basis_size",size(G));
print("modular_seconds",modular_seconds);
print("lift_seconds",lift_seconds);
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=1200,
    )
    print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip())


if __name__ == "__main__":
    main()
