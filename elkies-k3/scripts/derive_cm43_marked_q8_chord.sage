#!/usr/bin/env sage
"""Derive the exact degree-two chord pencil of the marked CM-43 q=8 class.

After fixed-component reduction the q=8 fiber is

    D8_nef = O + R - 2*F,   R=P1-2*P2.

Thus its old-fiber degree is two and the natural base function is the chord
slope ``W=(y+y(R))/(x-x(R))`` (up to a base-dependent normalization needed
to select the global two-dimensional subsystem).  This script performs the
exact line substitution, removes the known root x=x(R), and factors the
quadratic discriminant as a curve over QQ(W).  It is an equation-level
boundary calculation; the Jacobian/fiber signature is certified only after
the resulting square-free model is explicitly checked below.
"""

from pathlib import Path

from sage.all import *


load(str(Path(__file__).resolve().parent / "verify_cm43_q8_short_section.sage"))

KW = FunctionField(QQ, "W")
Wnew = KW.gen()
RTW = PolynomialRing(KW, "t")
t = RTW.gen()


def change_base(poly):
    poly = RT(poly)
    return RTW([KW(QQ(value)) for value in poly])


a4w = change_base(a4)
a6w = change_base(a6)
Nx = change_base(q8_x.numerator())
Ny = change_base(q8_y.numerator())
h = change_base(h6)
# Writing x=x(R)+u gives the quadratic
#
#   -u^2+(W^2-3*x(R))*u-(2*W*y(R)+3*x(R)^2+a4).
#
# Its discriminant is the compact expression below.  Clearing h^4 avoids
# the prohibitively expensive nested function-field resultant/gcd path.
numerator = RTW(
    h**4*Wnew**4
    - 6*Nx*h**2*Wnew**2
    - 8*Ny*h*Wnew
    - 3*Nx**2
    - 4*a4w*h**4
)
denominator = h**4
repeated_part = gcd(numerator, numerator.derivative())
squarefree_model = RTW(numerator//repeated_part)

print(
    f"CM43MARKEDQ8CHORD|raw_num_degree={numerator.degree()}"
    f"|raw_den_degree={denominator.degree()}|gcd_degree={repeated_part.degree()}"
    f"|denominator_is_h4={int(denominator == h**4)}",
    flush=True,
)
print(
    f"CM43MARKEDQ8CHORD|squarefree_degree={squarefree_model.degree()}"
    f"|squarefree={squarefree_model}",
    flush=True,
)

if squarefree_model.degree() in (3, 4):
    coefficients = [squarefree_model[index] for index in range(5)]
    coefficients.extend([KW(0)]*(5-len(coefficients)))
    e, d, c, b, aa = coefficients[:5]
    invariant_i = 12*aa*e-3*b*d+c**2
    invariant_j = 72*aa*c*e+9*b*c*d-27*aa*d**2-27*b**2*e-2*c**3
    jacobian = EllipticCurve(KW, [0, 0, 0, -27*invariant_i, -27*invariant_j])
    print(
        f"CM43MARKEDQ8CHORD|jacobian_a4={jacobian.a4()}"
        f"|jacobian_a6={jacobian.a6()}|jacobian_delta={jacobian.discriminant()}",
        flush=True,
    )
    print("CM43MARKEDQ8CHORD|status=GENUS_ONE_MODEL", flush=True)
else:
    print("CM43MARKEDQ8CHORD|status=NEEDS_BASE_NORMALIZATION", flush=True)
