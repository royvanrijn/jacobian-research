"""Mixed determinant obstruction for the normalized nonzero x^4 sector."""

from itertools import permutations

from sympy import Matrix, Poly, symbols


x, y, z, u, v, w, b, c, d, e = symbols("x y z u v w b c d e")
variables = (x, y, z)
p, q = z**2, x**2 + y*z
shifted = q + (w-v)*y
top = Matrix([z*p**2, z*p*q, z*q**2])
quartic = Matrix([
    q**2 + (2*u+b)*p*q + c*p**2,
    u*q**2 + p*(v*x**2+w*y*z) + d*p**2,
    2*q*(v*x**2+w*y*z) + e*p**2,
])
model = Matrix([
    z**5 + shifted**2 + (2*u+b)*p*shifted + c*p**2,
    z**3*shifted + u*shifted**2 + v*p*shifted + d*p**2,
    z*shifted**2 + 2*v*shifted**2 + e*p**2,
])

for entry, leading, following in zip(model, top, quartic, strict=True):
    assert Poly(entry-leading-following, *variables).total_degree() <= 3


def monomials(degree):
    return [
        x**i*y**j*z**(degree-i-j)
        for i in range(degree+1)
        for j in range(degree+1-i)
    ]


output_basis = monomials(9)


def column(expression):
    polynomial = Poly(expression, *variables)
    return Matrix([polynomial.coeff_monomial(monomial) for monomial in output_basis])


operator_columns = []
for coefficient in (q**2, -2*p*q, p**2):
    for monomial in monomials(2):
        contraction = coefficient*monomial
        operator_columns.append(column(5*z**3*(-z*contraction.diff(x)+2*x*contraction.diff(y))))
operator = Matrix.hstack(*operator_columns)
annihilators = operator.T.nullspace()
rows = Matrix.vstack(*(annihilators[i].T for i in (0, 6, 11, 17)))
assert rows*operator == Matrix.zeros(4, operator.cols)

top_jacobian, quartic_jacobian = (vector.jacobian(variables) for vector in (top, quartic))
tangent_directions = (
    Matrix([2*p*x, q*x, 0]),
    Matrix([2*p*y, q*y, 0]),
    Matrix([0, p*x, 2*q*x]),
    Matrix([0, p*y, 2*q*y]),
)
mixed_columns = []
for tangent in tangent_directions:
    jacobians = (top_jacobian, quartic_jacobian, tangent.jacobian(variables))
    mixed = sum(
        Matrix.hstack(*(jacobians[index][:, j] for j, index in enumerate(order))).det()
        for order in permutations(range(3))
    )
    mixed_columns.append(column(mixed))
minor = (rows*Matrix.hstack(*mixed_columns)).applyfunc(lambda coefficient: coefficient.expand())
assert minor == Matrix([
    [-2, 0, 0, 0],
    [0, 4, 0, 0],
    [9*b+30*u, 0, 30, 0],
    [0, -6*b-20*u, 0, -20],
])
assert minor.det() == 4800
print("Verified the matching model and the constant degree-nine compatibility minor.")
