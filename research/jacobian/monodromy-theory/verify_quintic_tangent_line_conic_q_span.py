"""Verify the final degree-five coefficient difference for the q span."""

from sympy import Matrix, Poly, Rational, symbols


x, y, z = symbols("x y z")
a, b, d, e, f, g, h, i, k, ell, m, n, r, v, w, c = symbols(
    "a b d e f g h i k ell m n r v w c"
)
p = z**2
q = x**2 + z*y
P = (
    z**3*q + a*p**2 + b*p*q + p*x + i*z**3 + k*z*q
    + 4*f*z*x + g*p + h*q + d*z + r*x
)
Q = (
    -Rational(1, 3)*z*q**2 + e*p*q + f*q**2 + q*x
    - Rational(16, 3)*a**2*z**3 + n*z*q - 4*a*z*x
    + ell*p + m*q + 2*y + v*x + w*z
)
R = z + z**2*(q + c)
jacobian = Matrix([P, Q, R]).jacobian((x, y, z))
determinant = Poly(jacobian.det(method="domain-ge"), x, y, z)
coefficient_difference = (
    determinant.coeff_monomial(x**2*z**3)
    - determinant.coeff_monomial(y*z**4)
)
assert coefficient_difference == -10
print("Verified the incompatible tangent-conic degree-five coefficients")
