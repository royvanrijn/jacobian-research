"""Check the Laurent and section coefficient identities in the three-simple exclusion."""

from sympy import Poly, Rational, cancel, expand, symbols


tau, S, U, V, L, sum_C2, sum_C2_u = symbols("tau S U V L Q2 QU")
phi = tau**-2 - S/tau + L - U*tau + V*tau**2
dq = tau**-4 - 2*S*tau**-3 - sum_C2*tau**-2 + sum_C2_u
difference = Poly(expand(tau**4*(phi**2 - dq)), tau)
assert difference.nth(2) == 2*L + S**2 + sum_C2
assert difference.nth(3) == -2*L*S - 2*U
assert difference.nth(4) == L**2 + 2*S*U + 2*V - sum_C2_u

p, a, c, b1, b2, e1, e2, e3 = symbols("p a c b1 b2 e1 e2 e3")
u1 = p**2 + a*p + b1
u2 = p**2 + a*p + b2
u3 = p**2 + c**2
v1 = -p**3 - Rational(3, 2)*a*p**2 - (Rational(3, 2)*b1 + Rational(3, 8)*a**2)*p + e1
v2 = -p**3 - Rational(3, 2)*a*p**2 - (Rational(3, 2)*b2 + Rational(3, 8)*a**2)*p + e2
v3 = -p**3 - Rational(3, 2)*c**2*p + e3
delta = b1 - b2
cleared = Poly(
    (u1-u3)**2*delta - c*(v1+v2)*(u1-u3) + p*(v1+v3)*delta,
    p,
)
assert expand(cleared.nth(4) - 2*(a*c-delta)) == 0
assert expand(
    cleared.nth(3) - c*(2*(b1-c**2)+3*a**2) + Rational(3, 2)*a*delta
) == 0

alpha, beta, slope, constant, leading = symbols("alpha beta slope constant leading")
J = 3*beta/(2*alpha)
K = 3*beta/(8*alpha**2)
u = alpha*p**2 + slope*p + constant
v = beta*p**3 + J*slope*p**2 + (J*constant+K*slope**2)*p + leading
section = Poly(v**2 - u**3, p)
assert cancel(section.nth(5).subs(beta**2, alpha**3)) == 0
assert cancel(section.nth(4).subs(beta**2, alpha**3)) == 0
exception = {alpha: Rational(9, 64), beta: -Rational(27, 512)}
assert (J - (beta-1)/(alpha-1) + alpha - 1).subs(exception) != 0

print("Verified the three-simple Laurent and section coefficient identities")
