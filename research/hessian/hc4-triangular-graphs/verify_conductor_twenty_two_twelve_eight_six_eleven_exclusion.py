"""Replay the quartic square-root coefficients and the curve's critical equations."""

from sympy import Matrix, Poly, Rational, diff, expand, rem, symbols


z, a, b, c, ell, j = symbols("z a b c ell j")
m = z**3 + a*z**2 + b*z + c
L = (z**4 + 3*a*z**3/2 + (3*b/2+3*a**2/8)*z**2
     + (3*c/2+3*a*b/4-a**3/16+j)*z + ell)
delta = a**2 - 4*b
remainder = Poly(rem(m**2, L, z), z)
constraints = Matrix([remainder.nth(3), remainder.nth(2)])
assert constraints.jacobian([c, ell]).det() != 0
solution = {
    c: -a*delta/8 + 2*j,
    ell: -(5*a**4-16*a**2*b-128*a*j-16*b**2)/64,
}
assert all(expand(e.subs(solution)) == 0 for e in constraints)
assert not any(v in constraints.jacobian([c, ell]).free_symbols for v in (c, ell))
remainder = Poly(rem(m.subs(solution)**2, L.subs(solution), z), z)
defect = Poly(m.subs(solution)**3-z*L.subs(solution)**2, z)
assert expand(remainder.nth(1)-delta*(a*delta-64*j)/128) == 0
assert expand(remainder.nth(0)-(delta**3+2*(a*delta-32*j)**2)/512) == 0
assert defect.degree() == 6 and defect.nth(6) == -2*j

X, Y, A, B, kappa, d = symbols("X Y A B kappa d", nonzero=True)
curve = X**9 - kappa*Y**3*X**2 + 2*Y**3*(A*Y+B)**2
critical = X**7 - 2*kappa*Y**3/9
assert expand(diff(curve, X) - 9*X*critical) == 0
assert rem(diff(curve, X, 2), critical, X) == 14*kappa*Y**3
assert expand(curve-X**2*critical-Y**3*(2*(A*Y+B)**2-7*kappa*X**2/9)) == 0
for sign in (1, -1):
    branch_value = (A*Y+B)**7-sign*d*Y**3
    assert expand(Y*diff(branch_value,Y)-3*branch_value-(A*Y+B)**6*(4*A*Y-3*B)) == 0
assert Rational(5, 9)*2+3 < 5
assert -18 + 8 + 12 > 0
assert -18 + 8 + 8 + 2 == 0

print("Verified the quartic square-root coefficients and the nine-sheeted curve calculations.")
