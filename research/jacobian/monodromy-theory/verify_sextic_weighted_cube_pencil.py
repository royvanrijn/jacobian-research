"""Exact coefficient checks for the weighted cube pencil."""

from itertools import product

from sympy import Matrix, Poly, QQ, expand, groebner, symbols

x, y, z = coordinates = symbols("x y z")
p = x * z + y**3
q = x**3
normal = (q**2, -2 * p * q, p**2)


def source_basis(weight):
    """All monomials of a fixed weight in weights (1, 1, 2)."""
    return [
        x**i * y**(weight - i - 2 * k) * z**k
        for k in range(weight // 2 + 1)
        for i in range(weight - 2 * k + 1)
    ]


def jacobian_rows(vector):
    return tuple(tuple(entry.diff(variable) for variable in coordinates) for entry in vector)


def determinant(row1, row2, row3):
    a, b, c = row1, row2, row3
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def source_coefficients(expression):
    return Poly(expand(expression), *coordinates).as_dict()


def remainder(expression, basis):
    expression = expand(expression)
    return expression if basis is None else basis.reduce(expression)[1]


def descend(weight, jacobians, constraints):
    """Solve the equation linear in the next triple of weighted parts."""
    basis = source_basis(weight)
    columns = [
        source_coefficients(2 * determinant(*jacobian_rows((p, q, factor * monomial))))
        for factor in normal
        for monomial in basis
    ]
    right = {}
    for i, j, k in product(jacobians, repeat=3):
        if i + j + k != weight + 12:
            continue
        term = determinant(jacobians[i][0], jacobians[j][1], jacobians[k][2])
        for monomial, coefficient in source_coefficients(term).items():
            right[monomial] = right.get(monomial, 0) - coefficient

    support = sorted(set(right).union(*(column.keys() for column in columns)))
    left = Matrix([[column.get(monomial, 0) for column in columns] for monomial in support])
    right_vector = Matrix([remainder(right.get(monomial, 0), constraints) for monomial in support])
    reduced, pivots = left.row_join(Matrix.eye(len(support))).rref()
    size = len(columns)
    free_indices = [index for index in range(size) if index not in pivots]
    free_parameters = symbols(f"u{weight}_0:{len(free_indices)}")
    coordinates_of_part = dict(zip(free_indices, free_parameters))
    transformed = reduced[:, size:] * right_vector
    compatibility = []
    for row, pivot in enumerate(pivots):
        value = remainder(transformed[row], constraints)
        if pivot < size:
            coordinates_of_part[pivot] = expand(
                value - sum(reduced[row, index] * coordinates_of_part[index] for index in free_indices)
            )
        elif value:
            compatibility.append(value)

    part = [
        expand(sum(monomial * coordinates_of_part[i * len(basis) + j]
                   for j, monomial in enumerate(basis)))
        for i in range(3)
    ]
    return part, free_parameters, compatibility


def verify_chart(chart):
    a, b, c, d, e, f, v, w = symbols("a b c d e f v w")
    variables = [a, b, c, d, e, f, v, w]
    charts = {
        "y-squared": {v: 1},
        "xy": {v: 0, w: 1},
        "zero": {v: 0, w: 0},
    }
    fixed = charts[chart]
    variables = [variable for variable in variables if variable not in fixed]
    lift = v * y**2 + w * x * y
    degree_five = (
        a * x**5 + b * x**2 * p,
        c * x**5 + d * x**2 * p + p * lift,
        e * x**5 + f * x**2 * p + 2 * q * lift,
    )
    jacobians = {
        6: jacobian_rows((p**2, p * q, q**2)),
        5: jacobian_rows([component.subs(fixed) for component in degree_five]),
    }
    constraints = None
    expected_dimensions = {4: 8, 3: 6, 2: 3, 1: 3}
    for weight in range(4, 0, -1):
        part, free, compatibility = descend(weight, jacobians, constraints)
        assert len(free) == expected_dimensions[weight], (chart, weight, len(free))
        variables.extend(free)
        old_generators = [] if constraints is None else list(constraints)
        generators = [*old_generators, *compatibility]
        if generators:
            constraints = groebner(generators, *variables, order="grevlex", domain=QQ)
        jacobians[weight] = jacobian_rows(part)

    origin = dict.fromkeys(coordinates, 0)
    # The first two columns come from weight one and the last from weight two.
    linear_jacobian = [
        (jacobians[1][i][0].subs(origin), jacobians[1][i][1].subs(origin),
         jacobians[2][i][2].subs(origin))
        for i in range(3)
    ]
    assert remainder(determinant(*linear_jacobian), constraints) == 0, chart
    print(f"{chart}: Jacobian at the origin reduces to zero")


def verify_nonzero_chart():
    a, b, c, d, e, f, v, w = symbols("a b c d e f v w")
    lift = v * y**2 + w * x * y + z
    jacobians = {
        6: jacobian_rows((p**2, p * q, q**2)),
        5: jacobian_rows((a * x**5 + b * x**2 * p,
                          c * x**5 + d * x**2 * p + p * lift,
                          e * x**5 + f * x**2 * p + 2 * q * lift)),
    }
    weight_four, free, compatibility = descend(4, jacobians, None)
    assert len(free) == 8
    constraints = groebner(compatibility, a, b, c, d, e, f, v, w, *free,
                           order="grevlex", domain=QQ)
    assert remainder(f, constraints) == 0
    expected = (0, 0, 1)
    assert tuple(remainder(item.coeff(z, 2) - coefficient, constraints)
                 for item, coefficient in zip(weight_four, expected, strict=True)) == (0, 0, 0)
    print("nonzero: the weight-four quadratic coefficient is (0, 0, 1)")


if __name__ == "__main__":
    for name in ("y-squared", "xy", "zero"):
        verify_chart(name)
    verify_nonzero_chart()
