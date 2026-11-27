"""Exact compatibility of the first two equations in the ternary quintic pencil."""

from sympy import Matrix, Poly, groebner, symbols


x, y, z = symbols("x y z")
p, q = z**2, x**2 + y * z
u = symbols("u0:6")
v = symbols("v0:6")
a, b, c, d, e = symbols("a b c d e")
parameters = (a, b, c, d, e, *u, *v)
quadratic_basis = (x**2, x*y, x*z, y**2, y*z, z**2)
U = sum(coefficient * monomial for coefficient, monomial in zip(u, quadratic_basis, strict=True))
V = sum(coefficient * monomial for coefficient, monomial in zip(v, quadratic_basis, strict=True))
top = Matrix([z*p**2, z*p*q, z*q**2])
quartic = Matrix([
    2*p*U + a*q**2 + b*p*q + c*p**2,
    q*U + p*V + d*p**2,
    2*q*V + e*p**2,
])
J5, J4 = (col.jacobian((x, y, z)) for col in (top, quartic))
mixed = Poly(sum(
    Matrix.hstack(*(J5[:, j] if j == index else J4[:, j] for j in range(3))).det()
    for index in range(3)
), x, y, z)

degree_ten = [x**i*y**j*z**(10-i-j) for i in range(11) for j in range(11-i)]
degree_three = [x**i*y**j*z**(3-i-j) for i in range(4) for j in range(4-i)]
normal = (q**2, -2*p*q, p**2)

def linear_column(monomial, component):
    f = normal[component] * monomial
    operator = Poly(5*z**3*(-z*f.diff(x) + 2*x*f.diff(y)), x, y, z)
    return Matrix([operator.coeff_monomial(m) for m in degree_ten])

operator = Matrix.hstack(*(
    linear_column(m, index) for index in range(3) for m in degree_three
))
constant_column = Matrix([mixed.coeff_monomial(m) for m in degree_ten])
assert operator.rank() == 20
conditions = [w.dot(constant_column).expand() for w in operator.T.nullspace()]
ideal = groebner(conditions, *parameters, order="grevlex", domain="QQ")

# These ideal memberships imply the six stated equations over every field.
for consequence in [u[1]**3, u[3]**2, v[3]**2, a*u[2], a**2*v[1], a*(u[0]-u[4])]:
    assert ideal.reduce(consequence)[1] == 0

# The two components of the stated conditions each annihilate the full system.
always = {u[1]: 0, u[3]: 0, v[3]: 0}
for branch in ({a: 0}, {u[2]: 0, v[1]: 0, u[4]: u[0]}):
    assert all(eq.subs(always | branch).expand() == 0 for eq in conditions)

print("Verified rank twenty and the necessary and sufficient first-equation conditions.")
