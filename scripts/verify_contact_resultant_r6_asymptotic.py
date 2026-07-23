#!/usr/bin/env python3
"""Exact branch-at-infinity certificate for the open contact column r=6.

The degree-29 endpoint eliminant is not treated as a growing Sylvester
resultant here.  With

    y = 1 + c/m,

the two endpoint equations have fixed limiting equations in (c,z).  Their
resultant is c^7 times a squarefree degree-29 polynomial P(c).  The c^7 is
the already removed endpoint factor.  A linear subresultant reconstructs a
unique algebraic limiting z on every root of P.

The Singular replay also checks directly that P is the complete top Newton
edge of H_6(m,1+y): every one of the 29 genuine branches therefore has

    y = 1 + c/m + O(m^-2),    z = d + O(m^-1).

Since P(0) != 0, Lindemann--Weierstrass gives exp(c) transcendental whereas
d is algebraic.  Thus d != exp(c) on all 29 branches, proving eventual
nonvanishing in the r=6 column.  The checker verifies the exact algebraic
inputs; Lindemann--Weierstrass is the cited transcendence input.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


m, y, z, c = sp.symbols("m y z c")
R = 6


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


def integral_polynomial(expression: sp.Expr) -> sp.Expr:
    """Clear parameter and rational-number denominators."""

    over_parameter_ring = sp.Poly(
        expression, z, y, domain=QQ.frac_field(m)
    ).clear_denoms(convert=True)[1]
    return sp.Poly(
        over_parameter_ring.as_expr(), z, y, m, domain=QQ
    ).clear_denoms(convert=True)[1].as_expr()


def beta_limit(k: int) -> sp.Rational:
    """Limit of m^(k+1) beta_k for k >= 1."""

    return sp.factorial(k) / sp.Integer(k) ** (k + 1)


def tail_limit(k: int) -> sp.Expr:
    """Limit of m^(k+1) T_k(1+c/m) for k >= 1.

    The displayed sum is the exact evaluation of

        integral_0^infinity exp(-k*s) (s-c)^k ds.
    """

    return sp.expand(
        sum(
            sp.binomial(k, j)
            * (-c) ** (k - j)
            * sp.factorial(j)
            / sp.Integer(k) ** (j + 1)
            for j in range(k + 1)
        )
    )


# Limits of m^6 E_6(m,1+c/m,z) and m^7 F_6(m,1+c/m,z).
# The k=0 summand of E_6 is c^6*z^5; the remaining terms use the
# incomplete-beta boundary-layer limits above.
limit_E6 = c**6 * z**5
for k in range(1, R):
    limit_E6 += (
        (-1) ** k
        * sp.binomial(R, k)
        * z ** (R - k - 1)
        * (beta_limit(k) - z**k * tail_limit(k))
        * (-c) ** (R - k - 1)
    )
limit_E6 = sp.factor(limit_E6)
limit_F6 = sp.factor(beta_limit(R) - z**R * tail_limit(R))

assert sp.degree(limit_E6, z) == 5
assert sp.degree(limit_F6, z) == 6

limit_resultant = sp.Poly(sp.resultant(limit_E6, limit_F6, z), c)
assert limit_resultant.degree() == 36
assert all(
    limit_resultant.nth(j) == 0 for j in range(7)
), "the endpoint factor must have order at least seven"

branch_polynomial = sp.Poly(
    sp.cancel(limit_resultant.as_expr() / c**7), c, domain=QQ
)
branch_polynomial = branch_polynomial.clear_denoms(convert=True)[1].primitive()[1]
assert branch_polynomial.degree() == 29
assert branch_polynomial.eval(0) != 0
assert sp.gcd(branch_polynomial, branch_polynomial.diff()).degree() == 0
assert sp.gcd(
    branch_polynomial, sp.Poly(tail_limit(R), c, domain=QQ)
).degree() == 0

# The penultimate subresultant gives A(c)z+B(c).  Coprimality of A and P
# makes z a unique algebraic rational function on every limiting branch.
subresultants = sp.subresultants(limit_E6, limit_F6, z)
assert [sp.degree(member, z) for member in subresultants] == [6, 5, 4, 3, 2, 1, 0]
linear_member = sp.Poly(subresultants[-2], z)
linear_coefficient = sp.Poly(linear_member.coeff_monomial(z), c, domain=QQ)
assert sp.gcd(branch_polynomial, linear_coefficient).degree() == 0


# Connect the limiting system to the already certified H_6 exactly.  In
# H_6(m,1+y), weight m=1 and y=-1.  Its maximal weight is 61, and that edge
# is a scalar multiple of P(y).  Hence all 29 roots, with multiplicity, are
# precisely the y-1 ~ c/m branches described by P.
D = 1 - y
E6 = sum(
    (-1) ** k
    * sp.binomial(R, k)
    * z ** (R - k - 1)
    * (beta(k) - y * z**k * endpoint_tail(k))
    * D ** (R - k - 1)
    for k in range(R)
)
F6 = beta(R) - y * z**R * endpoint_tail(R)
integral_E6 = sp.expand(integral_polynomial(E6))
integral_F6 = sp.expand(integral_polynomial(F6))


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


SINGULAR = shutil.which("Singular")
assert SINGULAR is not None, "the exact r=6 asymptotic certificate requires Singular"

program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E6)};
poly F={singular_expression(integral_F6)};
poly H=subst(resultant(E,F,z)/(y-1)^7,y,y+1);
poly remaining=H;
poly edge=0;
int violations=0;
intvec exponent;
int mdegree;
while(remaining!=0) {{
    exponent=leadexp(remaining);
    mdegree=exponent[3];
    if(exponent[3]-exponent[2]>61) {{violations=violations+1;}}
    if(exponent[3]-exponent[2]==61) {{
        edge=edge+lead(remaining)/m^mdegree;
    }}
    remaining=remaining-lead(remaining);
}}
poly expected={singular_expression(branch_polynomial.as_expr().subs(c, y))};
poly mismatch=leadcoef(expected)*edge-leadcoef(edge)*expected;
violations;
deg(edge);
size(edge);
size(mismatch);
"""

with tempfile.TemporaryDirectory(prefix="contact-r6-asymptotic-") as temporary:
    certificate = Path(temporary) / "certificate.sing"
    certificate.write_text(program)
    completed = subprocess.run(
        [SINGULAR, "-q", str(certificate)],
        stdin=subprocess.DEVNULL,
        check=True,
        capture_output=True,
        text=True,
    )

observed = [int(line) for line in completed.stdout.splitlines() if line.strip()]
assert observed == [0, 29, 30, 0], (observed, completed.stderr)

print("PASS: all 29 r=6 branches have y=1+c/m+O(m^-2)")
print("PASS: the degree-29 limiting polynomial is squarefree and excludes c=0")
print("PASS: a linear subresultant reconstructs a unique algebraic limiting z")
print("PASS: Lindemann--Weierstrass separates that z from exp(c) on every branch")
print("PASS: the r=6 contact resultant is nonzero for all sufficiently large m")
