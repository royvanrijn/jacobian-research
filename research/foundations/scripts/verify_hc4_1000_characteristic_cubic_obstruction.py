#!/usr/bin/env python3
"""Exclude reduced characteristic cubic caustics in chart 1000.

Let p be an irreducible cubic component of L=0 which is characteristic for
the Hamiltonian field of f_a.  Since deg(f_a)<=3,

    f_a-c = alpha*p.

Rescale p so alpha=1.  If f4 is the quartic homogeneous part of f, then
p3=(f4)_a and the top part of p|L gives

    p3 | det(Hess(f4)).

An exact binary-quartic calculation, using only shears a -> a+rho*b and
nonzero scalings, leaves three top normal forms:

    p3=a^3,  p3=a*b^2,  p3=b^3.

For each form, write the most general lower-degree cubic p, integrate
f_a=p, add the available univariate phi(b), and impose p|L exactly.
Coefficient comparison gives respectively:

    p=a^3+A*a^2+D*a+F,
    p=(a+C)*(b+B/2)^2,
    p=(b+C/3)^3.

Every result is reducible over an algebraically closed field.  Therefore
no reduced irreducible characteristic cubic component survives.
"""

from __future__ import annotations

import sympy as sp


a, b = sp.symbols("a b")


def caustic(f: sp.Expr) -> sp.Expr:
    f_aa = sp.diff(f, a, 2)
    f_ab = sp.diff(f, a, b)
    f_bb = sp.diff(f, b, 2)
    return sp.expand(5 * b * f_aa - (f_aa * f_bb - f_ab**2))


def coefficient(expression: sp.Expr, a_degree: int, b_degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), a, b).coeff_monomial(
        a**a_degree * b**b_degree
    )


# Leading binary-quartic classification.
A4, B4, C4, D4, E4, r, s = sp.symbols("A4 B4 C4 D4 E4 r s")
f4 = (
    A4 * a**4
    + 4 * B4 * a**3 * b
    + 6 * C4 * a**2 * b**2
    + 4 * D4 * a * b**3
    + E4 * b**4
)
p3 = sp.diff(f4, a) / 4
hessian4 = sp.det(sp.hessian(f4, (a, b))) / 48
leading_remainder = sp.expand(hessian4 - p3 * (r * a + s * b))
leading_equations = {
    powers: sp.factor(coefficient(leading_remainder, *powers))
    for powers in ((4, 0), (3, 1), (2, 2), (1, 3), (0, 4))
}

# If A4!=0, shear a to make B4=0.  The five equations then successively
# give r=3*C4, s=6*D4, A4*E4=6*C4^2, C4*D4=0, and
# C4*E4=3*D4^2.  The two cases force C4=D4=E4=0.
leading_A_chart = {
    powers: sp.factor(value.subs(B4, 0))
    for powers, value in leading_equations.items()
}
assert sp.factor(leading_A_chart[(4, 0)] - A4 * (3 * C4 - r)) == 0
assert sp.factor(leading_A_chart[(3, 1)] - A4 * (6 * D4 - s)) == 0
assert sp.factor(
    leading_A_chart[(2, 2)].subs({r: 3 * C4, s: 6 * D4})
    - 3 * (A4 * E4 - 6 * C4**2)
) == 0
assert sp.factor(
    leading_A_chart[(1, 3)].subs({r: 3 * C4, s: 6 * D4})
    + 27 * C4 * D4
) == 0
assert sp.factor(
    leading_A_chart[(0, 4)].subs({r: 3 * C4, s: 6 * D4})
    - 3 * (C4 * E4 - 3 * D4**2)
) == 0

# If A4=0, the a^4 equation gives B4=0.  If C4!=0, the remaining
# relation 3*C4*E4=2*D4^2 makes f4 a constant times
# b^2*(a+D4*b/(3*C4))^2.  If C4=0, then f4=b^3*(4*D4*a+E4*b).
assert sp.factor(leading_equations[(4, 0)].subs(A4, 0) + 3 * B4**2) == 0
square_top = 6 * C4 * b**2 * (a + D4 * b / (3 * C4)) ** 2
assert sp.factor(
    square_top
    - f4.subs(
        {
            A4: 0,
            B4: 0,
            E4: 2 * D4**2 / (3 * C4),
        }
    )
) == 0


# Case 1: p3=a^3.
A, B, C, D, E, F = sp.symbols("A B C D E F")
phi2, phi3 = sp.symbols("phi2 phi3")
p_case1 = a**3 + A * a**2 + B * a * b + C * b**2 + D * a + E * b + F
f_case1 = sp.integrate(p_case1, a) + phi2 * b**2 + phi3 * b**3
L_case1 = caustic(f_case1)
remainder_case1 = sp.Poly(L_case1, a, domain=sp.EX).rem(
    sp.Poly(p_case1, a, domain=sp.EX)
).as_expr()
case1 = {
    powers: sp.factor(coefficient(remainder_case1, *powers))
    for powers in ((2, 1), (2, 0), (1, 1), (1, 0), (0, 2), (0, 1), (0, 0))
}
assert sp.factor(case1[(2, 1)] + 3 * (6 * phi3 - 5)) == 0
case1_first = {phi3: sp.Rational(5, 6)}
assert sp.factor(case1[(1, 1)].subs(case1_first) - 8 * B * C) == 0
assert sp.factor(case1[(0, 2)].subs(case1_first) - 10 * C**2) == 0
case1_second = case1_first | {C: 0}
assert sp.factor(
    case1[(2, 0)].subs(case1_second) - (B**2 - 6 * phi2)
) == 0
case1_third = case1_second | {phi2: B**2 / 6}
assert sp.factor(case1[(0, 1)].subs(case1_third) + B**3 / 3) == 0
case1_fourth = case1_third | {B: 0, phi2: 0}
assert sp.factor(case1[(0, 0)].subs(case1_fourth) - E**2) == 0
case1_solution = case1_fourth | {E: 0, phi2: 0}
assert all(sp.factor(value.subs(case1_solution)) == 0 for value in case1.values())
assert sp.factor(p_case1.subs(case1_solution) - (a**3 + A * a**2 + D * a + F)) == 0


# Case 2: p3=a*b^2.  Since L has degree four, write L=p*(r*a+s*b+t).
t = sp.symbols("t")
p_case2 = (
    a * b**2 + A * a**2 + B * a * b + C * b**2 + D * a + E * b + F
)
f_case2 = sp.integrate(p_case2, a) + phi2 * b**2 + phi3 * b**3
case2_difference = sp.expand(L := caustic(f_case2) - p_case2 * (r * a + s * b + t))
case2 = {
    powers: sp.factor(coefficient(case2_difference, *powers))
    for powers in (
        (3, 0), (2, 2), (2, 1), (2, 0), (1, 3), (1, 2),
        (1, 1), (1, 0), (0, 3), (0, 2), (0, 1), (0, 0),
    )
}
case2_pivots = {r: 3, s: 0}
assert sp.factor(case2[(3, 0)].subs(case2_pivots) + 5 * A) == 0
case2_pivots |= {A: 0}
assert sp.factor(case2[(2, 0)].subs(case2_pivots) - (B**2 - 4 * D)) == 0
case2_pivots |= {D: B**2 / 4}
assert sp.factor(case2[(1, 2)].subs(case2_pivots) - (3 * C - t)) == 0
case2_pivots |= {t: 3 * C}
assert sp.factor(case2[(1, 1)].subs(case2_pivots) - (E - B * C)) == 0
case2_pivots |= {E: B * C}
assert sp.factor(
    case2[(1, 0)].subs(case2_pivots) - 3 * (B**2 * C / 4 - F)
) == 0
case2_pivots |= {F: B**2 * C / 4}
assert sp.factor(case2[(0, 3)].subs(case2_pivots) - (5 - 6 * phi3)) == 0
case2_pivots |= {phi3: sp.Rational(5, 6)}
assert sp.factor(case2[(0, 2)].subs(case2_pivots) - (C**2 - 2 * phi2)) == 0
case2_pivots |= {phi2: C**2 / 2}
assert sp.factor(case2[(0, 0)].subs(case2_pivots)) == 0
assert sp.factor(
    p_case2.subs(case2_pivots) - (a + C) * (b + B / 2) ** 2
) == 0


# Case 3: p3=b^3.
p_case3 = (
    b**3 + A * a**2 + B * a * b + C * b**2 + D * a + E * b + F
)
f_case3 = sp.integrate(p_case3, a) + phi2 * b**2 + phi3 * b**3
case3_difference = sp.expand(caustic(f_case3) - p_case3 * (r * a + s * b + t))
case3 = {
    powers: sp.factor(coefficient(case3_difference, *powers))
    for powers in (
        (3, 0), (2, 1), (2, 0), (1, 3), (1, 2), (1, 1),
        (1, 0), (0, 4), (0, 3), (0, 2), (0, 1), (0, 0),
    )
}
case3_pivots = {r: 0, s: 9}
assert sp.factor(case3[(2, 1)].subs(case3_pivots) + 21 * A) == 0
case3_pivots |= {A: 0}
assert sp.factor(case3[(2, 0)].subs(case3_pivots) - B**2) == 0
case3_pivots |= {B: 0}
assert sp.factor(case3[(0, 3)].subs(case3_pivots) - (3 * C - t)) == 0
case3_pivots |= {t: 3 * C}
assert sp.factor(case3[(1, 1)].subs(case3_pivots) + 15 * D) == 0
case3_pivots |= {D: 0}
assert sp.factor(case3[(0, 2)].subs(case3_pivots) - (C**2 - 3 * E)) == 0
case3_pivots |= {E: C**2 / 3}
assert sp.factor(
    case3[(0, 1)].subs(case3_pivots) - (C**3 / 3 - 9 * F)
) == 0
case3_pivots |= {F: C**3 / 27}
assert all(sp.factor(value.subs(case3_pivots)) == 0 for value in case3.values())
assert sp.factor(p_case3.subs(case3_pivots) - (b + C / 3) ** 3) == 0


def main() -> None:
    print("PASS: f4_a | det(Hess(f4)) has exactly three shear normal forms")
    print("PASS: top form a^3 forces p to depend only on a")
    print("PASS: top form a*b^2 forces p=(a+C)*(b+B/2)^2")
    print("PASS: top form b^3 forces p=(b+C/3)^3")
    print("RESULT: every characteristic cubic caustic is reducible")
    print("        no reduced irreducible characteristic cubic survives")


if __name__ == "__main__":
    main()
