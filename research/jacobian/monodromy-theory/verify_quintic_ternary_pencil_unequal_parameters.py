"""Verify the degree-nine witnesses and weighted form for unequal a=0 parameters."""

from itertools import permutations

from sympy import Matrix, Poly, expand, symbols

import verify_quintic_ternary_pencil_first_equations as first


x, y, z = first.x, first.y, first.z
u, v = first.u, first.v
a, b = first.a, first.b
normal = first.normal


def homogeneous_terms(degree):
    return [x**i * y**j * z**(degree-i-j)
            for i in range(degree+1) for j in range(degree+1-i)]


def binary_weight_terms(polynomial, weight):
    return sum(coefficient * x**i * y**j * z**k
               for (i, j, k), coefficient in Poly(polynomial, x, y, z).terms()
               if i + 2*j == weight)


def check():
    substitution = {a: 0, u[1]: 0, u[3]: 0, v[3]: 0}

    operator = first.operator
    pivots = operator.rref()[1]
    rows = operator[:, pivots].T.rref()[1]
    solution = -operator.extract(rows, pivots).inv() * first.constant_column.extract(rows, [0])
    particular = Matrix.zeros(30, 1)
    for index, coefficient in zip(pivots, solution, strict=True):
        particular[index] = coefficient.subs(substitution)
    free = symbols("f:10")
    coefficients = particular + sum(
        (parameter * column for parameter, column in zip(free, operator.nullspace(), strict=True)),
        Matrix.zeros(30, 1),
    )
    cubic = Matrix([
        sum(coefficients[10*i+j]*monomial for j, monomial in enumerate(first.degree_three))
        for i in range(3)
    ])
    quartic = first.quartic.subs(substitution)
    J5, J4, J3 = first.J5, quartic.jacobian((x, y, z)), cubic.jacobian((x, y, z))

    degree_nine = homogeneous_terms(9)
    mixed = Poly(J4.det() + sum(
        Matrix.hstack(*[columns[row][:, index] for index, row in enumerate(order)]).det()
        for columns in [(J5, J4, J3)] for order in permutations(range(3))
    ), x, y, z)
    right_side = Matrix([mixed.coeff_monomial(term) for term in degree_nine])

    def quadratic_column(monomial, component):
        expression = normal[component] * monomial
        polynomial = Poly(5*z**3*(-z*expression.diff(x)+2*x*expression.diff(y)), x, y, z)
        return Matrix([polynomial.coeff_monomial(term) for term in degree_nine])

    quadratic_operator = Matrix.hstack(*[
        quadratic_column(monomial, component)
        for component in range(3) for monomial in homogeneous_terms(2)
    ])
    equations = [expand(row.dot(right_side)) for row in quadratic_operator.T.nullspace()]
    assert expand(-30*v[1]**2*(3*v[0]-v[4])) in equations

    normalized = {v[1]: 0, u[5]: 0, v[2]: 0, v[5]: 0}
    normalized_equations = [expand(eq.subs(normalized)) for eq in equations]
    delta = u[0] - u[4]
    assert expand(4*delta*u[0]*(3*b+10*u[0])/5) in normalized_equations
    assert expand(4*delta*(6*b*u[0]-3*b*u[4]+10*u[0]**2)/5) in normalized_equations

    unequal = normalized | {u[0]: 0, b: 0, u[2]: 0}
    vector = (first.top + quartic + cubic).subs(unequal)
    for coordinate in vector:
        assert all(i + 2*j <= 4 for (i, j, _), _ in Poly(coordinate, x, y, z).terms())
    k, v0, h = u[4], v[0], v[4] - v[0]
    q = first.q
    expected = (
        k**2*z*y**2,
        (k*z+free[4])*q*y+k*h*z*y**2,
        (z+2*v0)*q**2+(2*h*z+free[8])*q*y+h**2*z*y**2,
    )
    for coordinate, leading in zip(vector, expected, strict=True):
        assert expand(binary_weight_terms(coordinate, 4) - leading) == 0


if __name__ == "__main__":
    check()
    print("Verified the unequal-parameter compatibility witnesses and weight-four form.")
