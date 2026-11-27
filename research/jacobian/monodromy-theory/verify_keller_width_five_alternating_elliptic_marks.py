"""Check the degree-four integral section and shared-jet identities."""

from sympy import Poly, Rational, expand, symbols


t, a, b, c, d = symbols("t a b c d")
u3, u2, u1, u0 = symbols("u3 u2 u1 u0")
section = d*t**4 + u3*t**3 + u2*t**2 + u1*t + u0
cubic = a*t**3 + b*t + c
cross_derivative = 3*section.diff(t)*cubic - section*cubic.diff(t)
equation = Poly(d*cross_derivative**2 - 81*a**2*(section**3+cubic), t)

assert expand(equation.nth(11) + 135*a**2*d**2*u3) == 0
substitutions = {u3: 0}
assert expand(
    equation.nth(10).subs(substitutions) + 9*a*d**2*(21*a*u2-22*b*d)
) == 0
substitutions[u2] = 22*b*d/(21*a)
assert expand(
    equation.nth(9).subs(substitutions) + 27*a*d**2*(9*a*u1-8*c*d)
) == 0
substitutions[u1] = 8*c*d/(9*a)
assert expand(
    equation.nth(8).subs(substitutions)
    + Rational(33, 49)*d**2*(441*a**2*u0-41*b**2*d)
) == 0
substitutions[u0] = 41*b**2*d/(441*a**2)
assert expand(
    equation.nth(6).subs(substitutions) + 640*b**3*d**3/(343*a)
) == 0
assert expand(
    cross_derivative.subs(substitutions).subs(b, 0)
    - d*(27*a**2*t**6 + 36*a*c*t**3 + 8*c**2)/(3*a)
) == 0

u0, u1, v0, v1, alpha, beta, gamma = symbols("u0 u1 v0 v1 alpha beta gamma")
u = u0 + u1*t + alpha*t**2
v = v0 + v1*t + beta*t**2 + gamma*t**3
shared_jet = Poly(v**2-u**3, t)
assert expand(shared_jet.nth(5) - (2*beta*gamma-3*u1*alpha**2)) == 0
assert expand(shared_jet.nth(6) - (gamma**2-alpha**3)) == 0

print("Verified the width-five integral-section and shared-jet identities")
