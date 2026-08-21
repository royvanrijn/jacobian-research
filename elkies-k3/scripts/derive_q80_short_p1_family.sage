#!/usr/bin/env sage
"""Derive the full short-P1 family in the q=80 ambient.

The D5 spinor and nonzero E6 component conditions give

    X=T*(1+c*T),  Y=T^2*(T-1)*(u+v*T),  d=X(1)=1+c.

Rather than carry the four interpolated B coefficients, define B from the
section identity.  Its T^3 normalization, A(1)=-3*d^2, and the remaining
three I4 discriminant jets cut out the complete one-section family.
"""

from sage.all import *


ring = PolynomialRing(QQ, names=("c", "u", "v", "p", "q", "r"), order="degrevlex")
c, u, v, p, q, r = ring.gens()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()
d = 1+c
A = T**2*(-3+p*T+q*T**2+r*T**3)
X = T*(1+c*T)
Y = T**2*(T-1)*(u+v*T)
B = Y**2-X**3-A*X
discriminant = 4*A**3+27*B**2
equations = (
    B[3]-2,
    A(1)+3*d**2,
    discriminant.derivative()(1),
    discriminant.derivative(2)(1),
    discriminant.derivative(3)(1),
)
ideal = ring.ideal(equations)
groebner = ideal.groebner_basis()
print(
    f"Q80SHORTP1|equations={len(equations)}|basis={len(groebner)}|"
    f"dimension={ideal.dimension()}|groebner={tuple(groebner)}",
    flush=True,
)
print("Q80SHORTP1|status=COMPUTED", flush=True)
