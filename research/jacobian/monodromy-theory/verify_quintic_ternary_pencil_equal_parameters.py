"""Verify the weighted forms used for the equal a=0 ternary parameters."""

from sympy import Matrix, Poly, expand, symbols

import verify_quintic_ternary_pencil_first_equations as first


def weighted_part(expression, weight):
    x, y, z = first.x, first.y, first.z
    return sum(coefficient*x**i*y**j*z**k
               for (i, j, k), coefficient in Poly(expression, x, y, z).terms()
               if i+2*j == weight)


def check():
    x, y, z, q = first.x, first.y, first.z, first.q
    a, u, v = first.a, first.u, first.v
    parameters = {a: 0, u[1]: 0, u[3]: 0, v[3]: 0, v[1]: 0,
                  u[5]: 0, v[2]: 0, v[5]: 0, u[4]: u[0]}
    operator = first.operator
    pivots = operator.rref()[1]
    rows = operator[:, pivots].T.rref()[1]
    solution = -operator.extract(rows, pivots).inv() * first.constant_column.extract(rows, [0])
    particular = Matrix.zeros(30, 1)
    for index, coefficient in zip(pivots, solution, strict=True):
        particular[index] = coefficient.subs(parameters, simultaneous=True)
    free = symbols("f:10")
    coefficients = particular + sum(
        (parameter*column for parameter, column in zip(free, operator.nullspace(), strict=True)),
        Matrix.zeros(30, 1),
    )
    cubic = Matrix([
        sum(coefficients[10*i+j]*monomial for j, monomial in enumerate(first.degree_three))
        for i in range(3)
    ])
    vector = first.top+first.quartic.subs(parameters, simultaneous=True)+cubic
    for coordinate in vector:
        assert all(i+2*j <= 4 for (i, j, _), _ in Poly(coordinate, x, y, z).terms())
    h = v[4]-v[0]
    expected = (
        0,
        u[0]*q**2+free[4]*q*y,
        (z+2*v[0])*q**2+(2*h*z+free[8])*q*y+h**2*z*y**2,
    )
    for coordinate, leading in zip(vector, expected, strict=True):
        assert expand(weighted_part(coordinate, 4)-leading) == 0


if __name__ == "__main__":
    check()
    print("Verified the equal-parameter weight-four forms.")
