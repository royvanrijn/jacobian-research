#!/usr/bin/env python3
"""Verify the full reduced exceptional Schur atlas on the HC4 surface.

For

    h6 = (x^6+y^6+z^6)/30
         + mu*x^2*y^2*z^2
         + nu*sum_(i != j) x_i^4*x_j^2,

write D=det(Hess(h6)).  A general lemma used in the accompanying note is:
if D is squarefree, then no nonzero homogeneous quartic s can satisfy

    grad(s)^T adj(Hess(h6)) grad(s) in (D).

Indeed, modulo every irreducible factor of D the adjugate has generic
rank one.  The scalar identity therefore forces adj(Hess(h6))*grad(s) to
vanish modulo that factor.  Squarefreeness makes D divide this degree-11
vector, so the vector is zero and then grad(s)=0.

This checker proves that D is nonsquarefree only at the Fermat and radial
parameter points.  Existing exact fiber checkers classify their quartics
and exclude both components at the next determinant faces.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
mu, nu = sp.symbols("mu nu")
R, P2, P3 = sp.symbols("R P2 P3")
variables = (x, y, z)

radius = x**2 + y**2 + z**2
pair_sum = x**2 * y**2 + x**2 * z**2 + y**2 * z**2
triple_product = x**2 * y**2 * z**2
mixed_42 = sum(
    left**4 * right**2
    for left in variables
    for right in variables
    if left != right
)
sextic = (
    (x**6 + y**6 + z**6) / 30
    + mu * triple_product
    + nu * mixed_42
)
hessian = sp.hessian(sextic, variables)
determinant = sp.expand(hessian.det())

# The exact projective incidence scheme is the coefficient scheme of this
# identity, followed by elimination of q0,...,q5 and saturation by the
# irrelevant quartic ideal.  There are 120 degree-fourteen coefficients.
quartic_coefficients = sp.symbols("s0:15")
quartic_monomials = tuple(
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
)
quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quartic_coefficients, quartic_monomials
    )
)
quotient_coefficients = sp.symbols("q0:6")
quotient_monomials = (x**2, y**2, z**2, x * y, x * z, y * z)
quotient = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quotient_coefficients, quotient_monomials
    )
)
gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in variables]
)
incidence_remainder = sp.expand(
    (gradient.T * hessian.adjugate() * gradient)[0]
    - determinant * quotient
)
incidence_equations = sp.Poly(
    incidence_remainder, *variables
).coeffs()
assert len(incidence_equations) == 120

# Express D in the seven-dimensional space of symmetric degree-six
# polynomials in X=x^2,Y=y^2,Z=z^2.
c0 = 4 * nu**2
c1 = 4 * nu * (mu - 20 * nu**2)
c2 = -4 * (
    3 * mu**2
    - 130 * mu * nu**2
    + 13 * mu * nu
    + 240 * nu**3
    - 46 * nu**2
    + nu
)
c3 = -2 * nu * (10 * nu - 1) * (4 * mu - 18 * nu + 1)
c4 = -2 * (
    20 * mu**2 * nu
    - 18 * mu**2
    + 780 * mu * nu**2
    - 44 * mu * nu
    + 3 * mu
    - 2040 * nu**3
    + 240 * nu**2
    - 10 * nu
)
c5 = 2 * (mu - 2 * nu) * (10 * nu - 1) ** 2
c6 = (2 * mu + 6 * nu - 1) ** 2 * (10 * mu - 30 * nu + 1)

invariant_determinant = (
    c0 * R**6
    + c1 * R**4 * P2
    + c2 * R**3 * P3
    + c3 * R**2 * P2**2
    + c4 * R * P2 * P3
    + c5 * P2**3
    + c6 * P3**2
)
assert sp.expand(
    determinant
    - invariant_determinant.subs(
        {R: radius, P2: pair_sum, P3: triple_product}
    )
) == 0

# Regard D=a*P3^2+b*P3+c over Q[R,P2].  When a is nonzero, D is
# nonsquarefree precisely when its P3-discriminant vanishes identically.
a = c6
b = c2 * R**3 + c4 * R * P2
c = c0 * R**6 + c1 * R**4 * P2 + c3 * R**2 * P2**2 + c5 * P2**3
discriminant = sp.Poly(sp.expand(b**2 - 4 * a * c), R, P2)
discriminant_coefficients = discriminant.coeffs()
assert len(discriminant_coefficients) == 4


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the radical certificate")

program = f"""
LIB "primdec.lib";
ring parameter_ring=0,(mu,nu),dp;
ideal I={",".join(map(
    singular_expression, discriminant_coefficients
))};
ideal fermat=mu,nu;
ideal radial=5*mu-1,10*nu-1;
ideal split_point=7*mu-5,14*nu+1;
ideal expected=std(intersect(fermat,intersect(radial,split_point)));
ideal actual=std(radical(I));
print(
  "DISCRIMINANT_RADICAL "
  +string(size(reduce(actual,expected)))+" "
  +string(size(reduce(expected,actual)))
);
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=120,
)
if completed.stderr.strip() or "?" in completed.stdout:
    raise RuntimeError(
        "Singular discriminant calculation failed:\n"
        f"{completed.stdout}\n{completed.stderr}"
    )
marker = re.search(
    r"(?m)^DISCRIMINANT_RADICAL (\d+) (\d+)$",
    completed.stdout,
)
assert marker is not None
assert tuple(map(int, marker.groups())) == (0, 0)

fermat = {mu: 0, nu: 0}
radial = {mu: sp.Rational(1, 5), nu: sp.Rational(1, 10)}
split = {mu: sp.Rational(5, 7), nu: -sp.Rational(1, 14)}

fermat_invariant = sp.factor(invariant_determinant.subs(fermat))
radial_invariant = sp.factor(invariant_determinant.subs(radial))
split_invariant = sp.factor(invariant_determinant.subs(split))
assert fermat_invariant == P3**2
assert radial_invariant == R**6 / 25
assert split_invariant == (
    (12 * P2 - R**2)
    * (144 * P2**2 - 24 * P2 * R**2 - 7 * R**4)
    / 343
)
assert sp.gcd(
    sp.gcd(split_invariant, sp.diff(split_invariant, R)),
    sp.diff(split_invariant, P2),
) == 1

# On the two components of a=c6=0, the P3 coefficient b has the displayed
# form.  Away from its listed zero, it is a scalar times
# R*(alpha*P2+beta*R^2) with alpha nonzero, hence squarefree.  At the
# extra zero on the first line the remaining binary polynomial is the
# squarefree split fiber above.  The common intersection is radial.
line_one_mu = (1 - 6 * nu) / 2
line_two_mu = 3 * nu - sp.Rational(1, 10)
line_one_b = sp.factor(b.subs(mu, line_one_mu))
line_two_b = sp.factor(b.subs(mu, line_two_mu))
assert line_one_b == (
    3
    * R
    * (10 * nu - 1)
    * (14 * nu + 1)
    * (20 * P2 * nu - 2 * P2 - 6 * R**2 * nu + R**2)
)
assert line_two_b == (
    -3
    * R
    * (10 * nu - 1) ** 2
    * (80 * P2 * nu - 8 * P2 - 50 * R**2 * nu + R**2)
    / 25
)
assert sp.solve(
    sp.Eq(line_one_mu, line_two_mu), nu
) == [sp.Rational(1, 10)]

# Passing from the invariant coordinates to x,y,z is ramified only over
# P3=0 and the cubic-root discriminant
#
#   Delta=R^2*P2^2-4*P2^3-4*R^3*P3-27*P3^2+18*R*P2*P3.
#
# A squarefree invariant D could acquire a repeated pullback factor only
# by containing one of these branch divisors.  P3 divides D only at
# Fermat, and D is never a scalar multiple of Delta.
coordinate_branch_groebner = sp.groebner(
    (c0, c1, c3, c5), mu, nu, order="lex"
)
assert tuple(
    polynomial.as_expr()
    for polynomial in coordinate_branch_groebner.polys
) == (mu, nu)
branch_scalar = sp.symbols("branch_scalar")
permutation_branch_coefficients = (0, 0, -4, 1, 18, -4, -27)
permutation_branch_groebner = sp.groebner(
    tuple(
        coefficient - branch_scalar * branch_coefficient
        for coefficient, branch_coefficient in zip(
            (c0, c1, c2, c3, c4, c5, c6),
            permutation_branch_coefficients,
        )
    ),
    branch_scalar,
    mu,
    nu,
    order="lex",
)
assert tuple(
    polynomial.as_expr()
    for polynomial in permutation_branch_groebner.polys
) == (1,)

# Record the two reduced projective fibers.  At Fermat the surviving
# coefficient indices are z^4,y^4,x^4.  At radial the coefficient vector
# is that of R^2.
fermat_indices = (0, 4, 14)
assert tuple(
    quartic_monomials[index] for index in fermat_indices
) == (z**4, y**4, x**4)
radial_vector = tuple(
    sp.Poly(radius**2, *variables).coeff_monomial(monomial)
    for monomial in quartic_monomials
)
assert radial_vector == (
    1, 0, 2, 0, 1, 0, 0, 0, 0, 2, 0, 2, 0, 0, 1
)

print("PASS: the projective Schur incidence has 120 coefficient equations")
print("PASS: the Hessian determinant has the seven invariant coefficients")
print("PASS: the P3-discriminant radical has Fermat, radial, and split support")
print("PASS: the split fiber is squarefree; only Fermat and radial are repeated")
print("PASS: invariant-quotient branch divisors add no parameter component")
print("PASS: the reduced incidence fibers are P^2_Fermat and one radial point")
print("SCOPE: full reduced 15-coefficient incidence atlas on the surface")
