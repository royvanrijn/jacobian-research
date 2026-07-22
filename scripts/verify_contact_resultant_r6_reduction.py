#!/usr/bin/env python3
"""Exact endpoint-eliminant reduction for the first open column r=6.

This is deliberately a reduction certificate, not a proof that the contact
resultant is nonzero.  It constructs the degree-five and degree-six endpoint
equations over Q[m,y], eliminates the auxiliary endpoint variable z with
Singular, and verifies the exact residual eliminant invariants.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


m, y, z = sp.symbols("m y z")
D = 1 - y
R = 6


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


# At a common root M_6=L_6=0, y and z=y^m are nonzero.  Divide the
# triangular contact identity by z and clear (1-y)^6.  The moment equation
# M_6=0 is already binomial in z.
E6 = sum(
    (-1) ** k
    * sp.binomial(R, k)
    * z ** (R - k - 1)
    * (beta(k) - y * z**k * endpoint_tail(k))
    * D ** (R - k - 1)
    for k in range(R)
)
F6 = beta(R) - y * z**R * endpoint_tail(R)

assert sp.degree(E6, z) == 5
assert sp.degree(F6, z) == 6


def integral_polynomial(expression: sp.Expr) -> sp.Expr:
    """Clear both parameter-function and rational-number denominators."""

    parameter_field = QQ.frac_field(m)
    over_parameter_ring = sp.Poly(
        expression, z, y, domain=parameter_field
    ).clear_denoms(convert=True)[1]
    return sp.Poly(
        over_parameter_ring.as_expr(), z, y, m, domain=QQ
    ).clear_denoms(convert=True)[1].as_expr()


integral_E6 = sp.expand(integral_polynomial(E6))
integral_F6 = sp.expand(integral_polynomial(F6))

assert sp.degree(integral_E6, z) == 5
assert sp.degree(integral_F6, z) == 6


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


SINGULAR = shutil.which("Singular")
assert SINGULAR is not None, "the exact r=6 reduction requires Singular"

program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E6)};
poly F={singular_expression(integral_F6)};
poly R=resultant(E,F,z);
ideal I7=(y-1)^7;
poly rem7=reduce(R,std(I7));
poly H=R/(y-1)^7;
ideal I1=y-1;
poly rem1=reduce(H,std(I1));
size(R);
deg(R,intvec(0,1,0));
deg(R,intvec(0,0,1));
size(rem7);
size(H);
deg(H,intvec(0,1,0));
deg(H,intvec(0,0,1));
size(rem1);
"""

with tempfile.TemporaryDirectory(prefix="contact-r6-") as temporary:
    certificate = Path(temporary) / "certificate.sing"
    certificate.write_text(program)
    completed = subprocess.run(
        [SINGULAR, "-q", str(certificate)],
        check=True,
        capture_output=True,
        text=True,
    )

observed = [int(line) for line in completed.stdout.splitlines() if line.strip()]
expected = [3361, 36, 90, 0, 2724, 29, 90, 62]
assert observed == expected, (observed, completed.stderr)

print("PASS: exact r=6 endpoint equations have z-degrees 5 and 6")
print("PASS: their resultant has exact endpoint factor (y-1)^7")
print("PASS: the residual eliminant has y-degree 29 and m-degree 90")
print("PASS: the residual eliminant is nonzero at y=1")
