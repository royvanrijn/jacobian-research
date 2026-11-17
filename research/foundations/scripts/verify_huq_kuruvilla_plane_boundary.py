#!/usr/bin/env python3
"""Exact boundary audit for the Huq--Kuruvilla--Mondello plane map.

The dependency-light block verifies the normalized-coordinate relations,
the reconstruction open, the two primes over Q=0, and the generic wild
local equation.  Singular independently computes the integral closure of
the monic cubic order, its conductor, the discriminant components upstairs,
and the two conductor branches.

The written proof and status boundary are in
verified/HUQ_KURUVILLA_PLANE_BOUNDARY_NORMALIZATION.md.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


def mod2(expr: sp.Expr, *generators: sp.Symbol) -> sp.Expr:
    """Return the canonical polynomial representative over F_2."""

    return sp.Poly(sp.expand(expr), *generators, modulus=2).as_expr()


def assert_zero(expr: sp.Expr, *generators: sp.Symbol) -> None:
    assert mod2(expr, *generators) == 0


# Source realization and the integral normalization generator Z=P^2/T.
x, y = sp.symbols("x y")
r = 1 + x * y
u = 1 + x**3 * r
P_source = x * r * u
Q_source = y + x**5 * r**3
T_source = r * u**2
Z_source = x**2 * r

P, Q, T, Z = sp.symbols("P Q T Z")
H = T**3 + T**2 + (P * Q + P**3) * T + P**3
relations = (
    Z * P * T + Z * P + P * Q + T**2 + T,
    Z * T + P**2,
    Z * P**2 * T + Z * P**2 + Z * Q * T + P * T**2 + P * T,
    Z**2 + Z * P**2 + Z * Q + P * T + P,
    H,
)
source_substitution = {
    P: P_source,
    Q: Q_source,
    T: T_source,
    Z: Z_source,
}
for relation in relations:
    assert_zero(relation.subs(source_substitution), x, y)
print("PASS: the source realizes the explicit normalized cubic presentation")

# The primitive cubic order has doubled coordinate-cross discriminant.
discriminant = sp.discriminant(H, T)
assert_zero(discriminant - (P * Q) ** 2, P, Q)
print("PASS: the hidden cubic discriminant is (P*Q)^2")

# Three denominator identities define x on the complement of their common
# zero curve.  The second source coordinate already belongs to the
# normalization ring.
D = T + P * Z
R = P + Z**2
assert_zero(x * Q_source - (T_source + 1), x, y)
assert_zero(x * D.subs(source_substitution) - P_source, x, y)
assert_zero(x * R.subs(source_substitution) - Z_source, x, y)
assert_zero(Q_source + Z_source**2 * (P_source + Z_source**2) - y, x, y)
print("PASS: the three x-charts and the global normalized formula for y agree")

# The boundary E and the ordinary companion A over Q=0.
z, p = sp.symbols("z p")
boundary_substitution = {P: z**2, Q: 0, T: z**3, Z: z}
ordinary_substitution = {P: p, Q: 0, T: 1, Z: p**2}
for relation in relations:
    assert_zero(relation.subs(boundary_substitution), z)
    assert_zero(relation.subs(ordinary_substitution), p)

factor_at_q_zero = sp.expand(H.subs(Q, 0))
assert_zero(factor_at_q_zero - (T + 1) * (T**2 + P**3), P, T)
assert_zero(D.subs(boundary_substitution), z)
assert_zero(R.subs(boundary_substitution), z)
assert_zero(D.subs(ordinary_substitution) - (1 + p**3), p)
assert_zero(R.subs(ordinary_substitution) - (p + p**4), p)
print("PASS: Q=0 has the ordinary sheet and the Frobenius boundary component")

# Near the generic boundary put pdev=P+Z^2.  The normalized equation becomes
# pdev*Z*(1+Z^3)+pdev^3+Q*Z^2=0.  Since Z(1+Z^3) is a
# generic unit, Q and pdev both have order one.  Implicit differentiation
# gives dP/dZ=pdev/(Z(1+Z^3)+pdev^2), hence different exponent one.
pdev = sp.symbols("pdev")
normalized_cubic = Z**3 + (Q + P**2) * Z**2 + P * Z + P**3
local_equation = pdev * Z * (1 + Z**3) + pdev**3 + Q * Z**2
assert_zero(normalized_cubic.subs(P, Z**2 + pdev) - local_equation, Z, Q, pdev)
local_dp = sp.diff(local_equation, pdev)
local_dz = sp.diff(local_equation, Z)
assert_zero(local_dp - (Z * (1 + Z**3) + pdev**2), Z, pdev)
assert_zero(local_dz - pdev, Z, Q, pdev)
print("PASS: the generic boundary has e=1, inseparable residue degree 2, and different exponent 1")

# A and E meet exactly where z^3=1: on E, the ordinary equations T=1 and
# Z=P^2 both reduce to that single separable condition.
assert_zero((T + 1).subs(boundary_substitution) - (z**3 + 1), z)
assert_zero((Z + P**2).subs(boundary_substitution) - z * (z**3 + 1), z)
assert sp.gcd(
    sp.Poly(z**3 + 1, z, modulus=2),
    sp.Poly(sp.diff(z**3 + 1, z), z, modulus=2),
).degree() == 0
print("PASS: the two Q=0 components meet in the three reduced cube roots of unity")

# At an intersection Z^3=1, Z is a unit.  The exact local parameter
# a=(Z(1+Z^3)+pdev^2)/Z^2 gives Q=pdev*a.  Its Z-derivative at pdev=0,
# Z^3=1 is Z^(-2), hence a unit.  Thus (pdev,a) are regular parameters,
# A=(a), E=(pdev), and the reduced boundary is an ordinary node.
a_numerator = Z * (1 + Z**3) + pdev**2
assert_zero(
    local_equation - (pdev * a_numerator + Q * Z**2), Z, Q, pdev
)
# In characteristic two the local equation is equivalent to
# Q=pdev*a_numerator/Z^2.  Check the substitution identically.
substituted_numerator = sp.together(
    local_equation.subs(Q, pdev * a_numerator / Z**2)
).as_numer_denom()[0]
assert_zero(substituted_numerator, Z, pdev)
a_derivative_numerator = sp.diff(a_numerator, Z)
assert_zero(a_derivative_numerator - 1, Z, pdev)
print("PASS: every retained/missing intersection has completed local equation Q=p*a and node conductor (p,a)")


SINGULAR_PROGRAM = r'''
LIB "primdec.lib";
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{
  if (reduce(f,std(G)) != 0)
  {
    "FAIL: "+label;
    exit(1);
  }
}

proc assertIdealEqual(ideal A, ideal B, string label)
{
  int i;
  for (i=1; i<=size(A); i++)
  {
    assertReductionZero(A[i],B,label);
  }
  for (i=1; i<=size(B); i++)
  {
    assertReductionZero(B[i],A,label);
  }
}

ring r=2,(P,Q,T),dp;
poly H=T3+T2+(P*Q+P3)*T+P3;
ideal cubic=H;

ideal jac=H,diff(H,P),diff(H,Q),diff(H,T);
ideal singularRadical=radical(jac);
assertIdealEqual(singularRadical,ideal(P,T),"singular locus");

list N=normal(cubic,"isPrim","withGens");
intvec normalTest=norTest(cubic,N);
if (normalTest[1] != 1 || normalTest[2] != 1 || normalTest[3] != 1)
{
  "FAIL: normalization test";
  exit(1);
}

ideal conductor=normalConductor(cubic);
assertIdealEqual(conductor,ideal(P,T),"conductor");

def Rn=N[1][1]; setring Rn;
poly ZZ=var(1);
ideal normalizationRelations=norid;

ideal E=normalizationRelations+ideal(Q,P+ZZ2,T+ZZ3);
ideal A=normalizationRelations+ideal(Q,T+1,ZZ+P2);
ideal qRadical=radical(normalizationRelations+ideal(Q));
ideal expectedQ=intersect(E,A);
assertIdealEqual(qRadical,expectedQ,"Q divisor");

ideal commonDenominators=radical(
  normalizationRelations+ideal(Q,T+P*ZZ,P+ZZ2)
);
assertIdealEqual(commonDenominators,E,"reconstruction boundary");

ideal C0=normalizationRelations+ideal(P,T,ZZ);
ideal C1=normalizationRelations+ideal(P,T,ZZ+Q);
ideal conductorRadical=radical(normalizationRelations+ideal(P,T));
ideal expectedConductor=intersect(C0,C1);
assertIdealEqual(conductorRadical,expectedConductor,"upstairs conductor");

"PASS_SINGULAR: normalization, conductor, boundary, and component decompositions";
'''


singular = shutil.which("Singular")
if singular is None:
    raise SystemExit("Singular is required for the normalization certificate")

completed = subprocess.run(
    [singular, "-q"],
    input=SINGULAR_PROGRAM,
    text=True,
    capture_output=True,
    check=False,
)
if completed.returncode != 0 or "PASS_SINGULAR" not in completed.stdout:
    raise AssertionError(
        "Singular normalization certificate failed\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
print("PASS: Singular certifies the normalization, conductor, and all declared component ideals")
print("PASS: plane boundary normalization audit complete")
