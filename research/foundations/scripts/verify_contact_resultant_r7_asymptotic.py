#!/usr/bin/env python3
"""Exact reduction and branch-at-infinity certificate for r=7 or r=8.

In the chart y=1+c/m, every residual endpoint-eliminant branch has ordinary
asymptotics and a unique algebraic limiting z.  Lindemann--Weierstrass then
proves nonvanishing for all sufficiently large m.  The default is r=7; set
``CONTACT_R=8`` for the next fixed column.  This checker supplies no effective
threshold.
"""

from __future__ import annotations

import os
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ


m, y, z, c = sp.symbols("m y z c")
R = int(os.environ.get("CONTACT_R", "7"))
CONFIG = {
    7: {
        "endpoint_factor": 7,
        "resultant_degree": 49,
        "branch_degree": 42,
        "m_degree": 126,
        "edge_weight": 84,
        "real_branches": 0,
        "singular_expected": [
            6350, 49, 126, 0, 680, 5461, 42, 126, 85, 0, 42, 43, 0
        ],
    },
    8: {
        "endpoint_factor": 9,
        "resultant_degree": 64,
        "branch_degree": 55,
        "m_degree": 200,
        "edge_weight": 145,
        "real_branches": 1,
        "singular_expected": [
            13057, 64, 200, 0, 1460, 11248, 55, 200, 146, 0, 55, 56, 0
        ],
    },
}
assert R in CONFIG
ENDPOINT_FACTOR = CONFIG[R]["endpoint_factor"]
RESULTANT_DEGREE = CONFIG[R]["resultant_degree"]
BRANCH_DEGREE = CONFIG[R]["branch_degree"]
M_DEGREE = CONFIG[R]["m_degree"]
EDGE_WEIGHT = CONFIG[R]["edge_weight"]
REAL_BRANCHES = CONFIG[R]["real_branches"]
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


# The sign of the k=0 summand is invisible in even columns but essential in
# odd columns.
limit_E = (-c) ** R * z ** (R - 1)
for k in range(1, R):
    limit_E += (
        (-1) ** k
        * sp.binomial(R, k)
        * z ** (R - k - 1)
        * (beta_limit(k) - z**k * tail_limit(k))
        * (-c) ** (R - k - 1)
    )
limit_E = sp.factor(limit_E)
limit_F = sp.factor(beta_limit(R) - z**R * tail_limit(R))
assert sp.degree(limit_E, z) == R - 1 and sp.degree(limit_F, z) == R

limit_resultant = sp.Poly(sp.resultant(limit_E, limit_F, z), c)
assert limit_resultant.degree() == RESULTANT_DEGREE
assert all(limit_resultant.nth(j) == 0 for j in range(ENDPOINT_FACTOR))
assert limit_resultant.nth(ENDPOINT_FACTOR) != 0
branch_polynomial = sp.Poly(
    sp.cancel(limit_resultant.as_expr() / c**ENDPOINT_FACTOR), c, domain=QQ
).clear_denoms(convert=True)[1].primitive()[1]
assert branch_polynomial.degree() == BRANCH_DEGREE
assert branch_polynomial.eval(0) != 0
assert sp.gcd(branch_polynomial, branch_polynomial.diff()).degree() == 0
assert sp.gcd(
    branch_polynomial, sp.Poly(tail_limit(R), c, domain=QQ)
).degree() == 0
assert branch_polynomial.count_roots(-sp.oo, sp.oo) == REAL_BRANCHES

subresultants = sp.subresultants(limit_E, limit_F, z)
assert [sp.degree(member, z) for member in subresultants] == list(
    range(R, -1, -1)
)
linear = sp.Poly(subresultants[-2], z)
linear_coefficient = sp.Poly(linear.coeff_monomial(z), c, domain=QQ)
assert sp.gcd(branch_polynomial, linear_coefficient).degree() == 0


E = sum(
    (-1) ** k
    * sp.binomial(R, k)
    * z ** (R - k - 1)
    * (beta(k) - y * z**k * endpoint_tail(k))
    * D ** (R - k - 1)
    for k in range(R)
)
F = beta(R) - y * z**R * endpoint_tail(R)


def integral_polynomial(expression: sp.Expr) -> sp.Expr:
    over_parameter_ring = sp.Poly(
        expression, z, y, domain=QQ.frac_field(m)
    ).clear_denoms(convert=True)[1]
    return sp.Poly(
        over_parameter_ring.as_expr(), z, y, m, domain=QQ
    ).clear_denoms(convert=True)[1].as_expr()


integral_E = sp.expand(integral_polynomial(E))
integral_F = sp.expand(integral_polynomial(F))


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


singular = shutil.which("Singular")
assert singular is not None, "the exact fixed-r certificate requires Singular"

program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E)};
poly F={singular_expression(integral_F)};
poly RR=resultant(E,F,z);
poly H=RR/(y-1)^{ENDPOINT_FACTOR};
poly shifted=subst(H,y,y+1);
poly remaining=shifted;
poly edge=0;
int violations=0;
intvec exponent;
int mdegree;
while(remaining!=0) {{
    exponent=leadexp(remaining);
    mdegree=exponent[3];
    if(exponent[3]-exponent[2]>{EDGE_WEIGHT}) {{violations=violations+1;}}
    if(exponent[3]-exponent[2]=={EDGE_WEIGHT}) {{
        edge=edge+lead(remaining)/m^mdegree;
    }}
    remaining=remaining-lead(remaining);
}}
poly expected={singular_expression(branch_polynomial.as_expr().subs(c, y))};
poly mismatch=leadcoef(expected)*edge-leadcoef(edge)*expected;
size(RR);
deg(RR,intvec(0,1,0));
deg(RR,intvec(0,0,1));
size(reduce(RR,std(ideal((y-1)^{ENDPOINT_FACTOR}))));
size(reduce(RR,std(ideal((y-1)^{ENDPOINT_FACTOR + 1}))));
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
assert observed == CONFIG[R]["singular_expected"], (observed, completed.stderr)

print(
    f"PASS contact resultant r={R}: residual eliminant has "
    f"bidegree ({BRANCH_DEGREE},{M_DEGREE})"
)
print(
    f"PASS contact resultant r={R}: all {BRANCH_DEGREE} branches have "
    "y=1+c/m+O(m^-2)"
)
print(
    f"PASS contact resultant r={R}: the limiting polynomial is squarefree "
    f"with {REAL_BRANCHES} real roots"
)
print(
    f"PASS contact resultant r={R}: a linear subresultant reconstructs "
    "a unique algebraic z"
)
print(
    f"PASS contact resultant r={R}: "
    "Lindemann--Weierstrass gives eventual nonvanishing"
)
