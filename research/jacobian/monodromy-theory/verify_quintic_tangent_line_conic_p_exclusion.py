"""Verify four determinant identities in the tangent pure-fifth normal form."""

from sympy import Matrix, Poly, Rational, expand, symbols


x, y, z = symbols("x y z")
A, B, C, D = symbols("A B C D")
a, b, e, f, i, k, j, n, c = symbols("a b e f i k j n c")
p2_coefficients = symbols("u0:6")
q2_coefficients = symbols("v0:6")
p1_coefficients = symbols("r0:3")
q1_coefficients = symbols("s0:3")

p = z**2
q = x**2 + z*y
L = A*x + B*y
V = C*x + D*y
quadratic_monomials = (x**2, x*y, y**2, z*x, z*y, z**2)
linear_monomials = (x, y, z)


def polynomial(coefficients, monomials):
    return sum(coefficient*monomial for coefficient, monomial in zip(coefficients, monomials))


P = (
    Rational(1, 5)*z**5
    + p*z*L + a*p**2 + b*p*q
    + p*V + 2*z*L**2 + 4*a*p*L + i*z**3 + k*z*q
    + polynomial(p2_coefficients, quadratic_monomials)
    + polynomial(p1_coefficients, linear_monomials)
)
Q = (
    z**3*q
    + q*z*L + e*p**2 + f*p*q
    + 4*e*p*L + q*V + j*z**3 + n*z*q
    + polynomial(q2_coefficients, quadratic_monomials)
    + polynomial(q1_coefficients, linear_monomials)
)
R = z + p*(q + c)

jacobian = Matrix([P, Q, R]).jacobian((x, y, z))
determinant = Poly(jacobian.det(method="domain-ge"), x, y, z)
assert expand(determinant.coeff_monomial(x**3*y**2*z**3) + 12*B**3) == 0

u1 = p2_coefficients[1]
expected_coefficients = {
    x**4*z**4: 6*A**3 - A*k + 2*u1,
    y**2*z**6: -A*k - u1,
    x**5*z**2: -4*A*u1,
}
for monomial, expected in expected_coefficients.items():
    coefficient = determinant.coeff_monomial(monomial).subs({B: 0, D: 0})
    assert expand(coefficient - expected) == 0

print("Verified the four incompatible tangent pure-fifth determinant coefficients")
