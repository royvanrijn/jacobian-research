"""Replay the square-derivative identities and the final four-equation ideal."""

from sympy import Poly, QQ, Rational, cancel, expand, groebner, rem, symbols


x, y, a, u, v, B, C, L = symbols("x y a u v B C L")
k, ell, m, n = symbols("k ell m n")
b = v**2 - u**3 - a*u
f = (y+v)/(x-u)
phi = y + B*x + C*f + L
q = k*x*y + ell*x**2 + m*y + n*x - C**2*f

numerator = cancel(
    (phi**2 - 2*y*q.diff(x) - (3*x**2+a)*q.diff(y))*(x-u)**2
)
numerator = rem(numerator, y**2-x**3-a*x-b, y)


def coefficient(expression, monomial):
    return Poly(expression, x, y).coeff_monomial(monomial)


assert coefficient(numerator, x**5) == 1-5*k
assert coefficient(numerator, x**3*y) == 2*B-4*ell
numerator = numerator.subs({k: Rational(1, 5), ell: B/2})
assert coefficient(numerator, x**4) == B**2+2*C-3*m
assert expand(coefficient(numerator, x**2*y) - 2*(B*C+L-n)) == 0
numerator = numerator.subs({m: (B**2+2*C)/3, n: B*C+L})
assert expand(coefficient(numerator, x*y) - 2*C*(B*u+L+v)) == 0

equation_one = u*B**2+v*B-C**2-C*u-a/5
equation_two = (
    5*a*B**2-15*B**2*u**2-30*B*C*v-30*B*u*v-15*C**2*u
    -20*C*a-30*C*u**2+9*a*u+9*u**3-24*v**2
)
assert expand(
    numerator.subs(L, -B*u-v)
    + (2*equation_one*x+equation_two/15)*(x-u)**2
) == 0

p, A, b0, h, p0 = symbols("p A b0 h p0")
iv_star = equation_one.subs(
    {a: 0, u: 1+A**2*p**2, v: A*p*(Rational(3, 2)+A**2*p**2),
     B: -A*p+b0, C: h*(p-p0)}
)
iv_star = Poly(rem(expand(iv_star), A**2+Rational(4, 3), A), p)
assert expand(iv_star.nth(3) - 4*(A*b0+h)/3) == 0
assert expand(iv_star.nth(2) + (4*b0**2+3*h**2+4*h*p0-2)/3) == 0
assert expand(iv_star.nth(1) - (-A*b0+4*h**2*p0-2*h)/2) == 0

c, d, e, z, w = symbols("c d e z w")
i2_substitution = {
    a: -3*(p**2-1), u: -2*p+c, v: d*p+e,
    B: d*z, C: h*p+c*w,
}
eq_one = Poly(equation_one.subs(i2_substitution), p)
eq_two = Poly(equation_two.subs(i2_substitution), p)


def reduce_section(expression):
    expression = expand(expression.subs(e, -20/d))
    expression = rem(expression, d**2-9*c, d)
    return rem(expression, c**2-Rational(16, 3), c)


relations = [
    5*h**2-10*h-3,
    -2*h*w-h+2*w-18*z**2+9*z,
    -80*w**2-80*w+720*z**2-300*z-9,
    h**2-4*h*w+18*h*z-8*h+4*w+45*z**2-36*z+9,
]
assert reduce_section(eq_one.nth(2)+relations[0]/5) == 0
assert reduce_section(eq_one.nth(1)-c*relations[1]) == 0
assert reduce_section(eq_one.nth(0)-relations[2]/15) == 0
assert reduce_section(eq_two.nth(2)+15*c*relations[3]) == 0
basis = groebner(relations, w, z, h, domain=QQ)
assert [value.as_expr() for value in basis.polys] == [1]

print("Verified the simple-pole differential identities and unit ideal")
