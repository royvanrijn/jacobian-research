#!/usr/bin/env python3
"""Exact reduction and branch-at-infinity certificate for r=7.

The seventh cancellation contact-resultant column reduces to a degree-42
endpoint eliminant.  In the chart y=1+c/m, all 42 branches have ordinary
asymptotics and a unique algebraic limiting z.  Lindemann--Weierstrass then
proves nonvanishing for all sufficiently large m.  This checker supplies no
effective threshold; it is the r=7 analogue of the former r=6 tail theorem.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ


m, y, z, c = sp.symbols("m y z c")
R = 7
D = 1 - y


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


def beta_limit(k: int) -> sp.Rational:
    return sp.factorial(k) / sp.Integer(k) ** (k + 1)


def tail_limit(k: int) -> sp.Expr:
    return sp.expand(
        sum(
            sp.binomial(k, j)
            * (-c) ** (k - j)
            * sp.factorial(j)
            / sp.Integer(k) ** (j + 1)
            for j in range(k + 1)
        )
    )


# The k=0 summand is (-c)^7*z^6.  The sign is invisible in the even r=6
# column but essential here.
limit_E7 = (-c) ** R * z ** (R - 1)
for k in range(1, R):
    limit_E7 += (
        (-1) ** k
        * sp.binomial(R, k)
        * z ** (R - k - 1)
        * (beta_limit(k) - z**k * tail_limit(k))
        * (-c) ** (R - k - 1)
    )
limit_E7 = sp.factor(limit_E7)
limit_F7 = sp.factor(beta_limit(R) - z**R * tail_limit(R))
assert sp.degree(limit_E7, z) == 6 and sp.degree(limit_F7, z) == 7

limit_resultant = sp.Poly(sp.resultant(limit_E7, limit_F7, z), c)
assert limit_resultant.degree() == 49
assert all(limit_resultant.nth(j) == 0 for j in range(7))
assert limit_resultant.nth(7) != 0
branch_polynomial = sp.Poly(
    sp.cancel(limit_resultant.as_expr() / c**7), c, domain=QQ
).clear_denoms(convert=True)[1].primitive()[1]
assert branch_polynomial.degree() == 42
assert branch_polynomial.eval(0) != 0
assert sp.gcd(branch_polynomial, branch_polynomial.diff()).degree() == 0
assert sp.gcd(
    branch_polynomial, sp.Poly(tail_limit(R), c, domain=QQ)
).degree() == 0
assert branch_polynomial.count_roots(-sp.oo, sp.oo) == 0

subresultants = sp.subresultants(limit_E7, limit_F7, z)
assert [sp.degree(member, z) for member in subresultants] == [
    7,
    6,
    5,
    4,
    3,
    2,
    1,
    0,
]
linear = sp.Poly(subresultants[-2], z)
linear_coefficient = sp.Poly(linear.coeff_monomial(z), c, domain=QQ)
assert sp.gcd(branch_polynomial, linear_coefficient).degree() == 0


E7 = sum(
    (-1) ** k
    * sp.binomial(R, k)
    * z ** (R - k - 1)
    * (beta(k) - y * z**k * endpoint_tail(k))
    * D ** (R - k - 1)
    for k in range(R)
)
F7 = beta(R) - y * z**R * endpoint_tail(R)


def integral_polynomial(expression: sp.Expr) -> sp.Expr:
    over_parameter_ring = sp.Poly(
        expression, z, y, domain=QQ.frac_field(m)
    ).clear_denoms(convert=True)[1]
    return sp.Poly(
        over_parameter_ring.as_expr(), z, y, m, domain=QQ
    ).clear_denoms(convert=True)[1].as_expr()


integral_E7 = sp.expand(integral_polynomial(E7))
integral_F7 = sp.expand(integral_polynomial(F7))


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "the exact r=7 certificate requires Singular"

program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E7)};
poly F={singular_expression(integral_F7)};
poly RR=resultant(E,F,z);
poly H=RR/(y-1)^7;
poly shifted=subst(H,y,y+1);
poly remaining=shifted;
poly edge=0;
int violations=0;
intvec exponent;
int mdegree;
while(remaining!=0) {{
    exponent=leadexp(remaining);
    mdegree=exponent[3];
    if(exponent[3]-exponent[2]>84) {{violations=violations+1;}}
    if(exponent[3]-exponent[2]==84) {{
        edge=edge+lead(remaining)/m^mdegree;
    }}
    remaining=remaining-lead(remaining);
}}
poly expected={singular_expression(branch_polynomial.as_expr().subs(c, y))};
poly mismatch=leadcoef(expected)*edge-leadcoef(edge)*expected;
size(RR);
deg(RR,intvec(0,1,0));
deg(RR,intvec(0,0,1));
size(reduce(RR,std(ideal((y-1)^7))));
size(reduce(RR,std(ideal((y-1)^8))));
size(H);
deg(H,intvec(0,1,0));
deg(H,intvec(0,0,1));
size(reduce(H,std(ideal(y-1))));
violations;
deg(edge);
size(edge);
size(mismatch);
"""

completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
)
assert completed.stderr == "", completed.stderr
observed = [int(line) for line in completed.stdout.splitlines() if line.strip()]
expected = [6350, 49, 126, 0, 680, 5461, 42, 126, 85, 0, 42, 43, 0]
assert observed == expected, (observed, completed.stderr)

print("PASS contact resultant r=7: residual eliminant has bidegree (42,126)")
print("PASS contact resultant r=7: all 42 branches have y=1+c/m+O(m^-2)")
print("PASS contact resultant r=7: the limiting polynomial is squarefree and has no real roots")
print("PASS contact resultant r=7: a linear subresultant reconstructs a unique algebraic z")
print("PASS contact resultant r=7: Lindemann--Weierstrass gives eventual nonvanishing")
