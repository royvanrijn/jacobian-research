#!/usr/bin/env sage
"""Recover the short polynomial CM24 basis on the quadratic q=80 surface.

The CM24 lattice audit gives polynomial-basis profiles

    D1=(1,1,1,1), D2=(0,0,1,2), D3=(1,1,0,2)

in (A1,A3,D5,E6) order.  The nonzero E6 labels force degree at most (2,4)
for (x,y).  This script first solves the especially small D2 chart, whose D5
spinor condition fixes x=T*(1+c*T).
"""

from sage.all import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--sign", type=int, choices=(-1, 1), default=-1)
args = parser.parse_args()

quadratic.<s> = QuadraticField(-3)
p = 72+args.sign*27*s
q = 18-2*p
e = (p-42)**2/36

ring = PolynomialRing(quadratic, names=("c", "u0", "u1", "u2"), order="degrevlex")
c, u0, u1, u2 = ring.gens()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3 + (2*p+e-45)*T**4 + (-9*p-4*e+186)*T**5
    + (12*p+6*e-299)*T**6 + (-5*p-4*e+210)*T**7 + e*T**8
)
X = T*(1+c*T)
Y = T**2*(u0+u1*T+u2*T**2)
identity = Y**2-X**3-A*X-B
equations = tuple(ring(value) for value in identity.list() if value)
ideal = ring.ideal(equations)
groebner = ideal.groebner_basis()
print(
    f"Q80CM24SHORT|stage=D2_groebner|equations={len(equations)}|"
    f"basis={len(groebner)}|dimension={ideal.dimension()}|groebner={tuple(groebner)}",
    flush=True,
)
solutions = ideal.variety(ring=quadratic)
print(f"Q80CM24SHORT|stage=D2_solutions|count={len(solutions)}", flush=True)
for index, solution in enumerate(solutions, 1):
    X_value = X.subs(solution)
    Y_value = Y.subs(solution)
    assert Y_value**2 == X_value**3+A*X_value+B
    print(f"Q80CM24SHORT|D2={index}|X={X_value}|Y={Y_value}", flush=True)
print("Q80CM24SHORT|status=PASS", flush=True)
