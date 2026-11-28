"""Exact finite algebra in the six-cell (2,2,2) gap proof.

This checks the Young convention, cubic field identity, and all coefficients
of the transported adjacent residual. The determinant inequalities and the
operator comparison are proved in the accompanying mathematical note.
"""
import json
from fractions import Fraction as Q
from itertools import combinations


def transpose(matrix):
    return tuple(zip(*matrix))


def multiply(left, right):
    columns = transpose(right)
    return tuple(
        tuple(sum(x * y for x, y in zip(row, column)) for column in columns)
        for row in left
    )


def apply(matrix, column):
    return tuple(sum(x * y for x, y in zip(row, column)) for row in matrix)


def identity(size):
    return tuple(tuple(int(row == col) for col in range(size)) for row in range(size))


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def add(left, right):
    return tuple(x + y for x, y in zip(left, right))


def subtract(left, right):
    return tuple(x - y for x, y in zip(left, right))


def exterior_square(matrix):
    pairs = tuple(combinations(range(3), 2))
    return tuple(
        tuple(
            matrix[i][a] * matrix[j][b] - matrix[j][a] * matrix[i][b]
            for a, b in pairs
        )
        for i, j in pairs
    )


H = ((1, -1, 0, 0), (1, 1, -2, 0), (1, 1, 1, -3))
D = ((Q(1, 2), 0, 0), (0, Q(1, 6), 0), (0, 0, Q(1, 12)))
X = ((2, 1, 0), (-Q(4, 3), Q(4, 3), Q(2, 3)), (1, -1, 1))
YOUNG = (
    ((-1, 0, 0), (0, -Q(1, 3), Q(8, 9)), (0, 1, Q(1, 3))),
    ((-Q(1, 2), Q(3, 4), 0), (1, Q(1, 2), 0), (0, 0, -1)),
    ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
)
FIRST_TRANSPORT = ((1, 0, 0), (0, 1, 0), (0, 0, -2))


def check_young():
    # Tableau entries are box coordinates in numerical label order.
    a, b, c, d = (0, 1), (1, 1), (2, 0), (2, 1)
    prefix = ((0, 0), (1, 0))
    tableaux = tuple(prefix + suffix for suffix in ((a, b, c, d), (a, c, b, d), (c, a, b, d)))
    inversions = {}
    for tableau in tableaux:
        reading = sorted(range(6), key=tableau.__getitem__)
        inversions[tableau] = sum(reading[i] > reading[j] for i, j in combinations(range(6), 2))

    for k, generator in enumerate(YOUNG):
        actual = [[0] * 3 for _ in range(3)]
        for column, tableau in enumerate(tableaux):
            index = k + 2
            r, c = tableau[index]
            s, d = tableau[index + 1]
            axial = (d - s) - (c - r)
            actual[column][column] = Q(1, axial)
            adjacent = list(tableau)
            adjacent[index], adjacent[index + 1] = adjacent[index + 1], adjacent[index]
            adjacent = tuple(adjacent)
            if adjacent in inversions:
                row = tableaux.index(adjacent)
                actual[row][column] = (
                    1 if inversions[adjacent] > inversions[tableau] else 1 - Q(1, axial**2)
                )
        assert tuple(map(tuple, actual)) == generator
        assert multiply(generator, generator) == identity(3)

        permutation = list(identity(4))
        permutation[k], permutation[k + 1] = permutation[k + 1], permutation[k]
        helmert = multiply(multiply(multiply(H, permutation), transpose(H)), D)
        assert multiply(X, exterior_square(helmert)) == multiply(generator, X)

    # Compression of 2 R((23)) to (12)-minus leaves only the diagonal:
    # the tableau obtained by an admissible swap leaves the minus space.
    compressed = []
    for tableau in tableaux:
        r, c = tableau[1]
        s, d = tableau[2]
        compressed.append(Q(2, (d - s) - (c - r)))
        swapped = (tableau[0], tableau[2], tableau[1], *tableau[3:])
        assert swapped not in tableaux
    assert tuple(compressed) == tuple(FIRST_TRANSPORT[i][i] for i in range(3))


def field(points):
    a, b = (apply(H, component) for component in transpose(points))
    wedge = tuple(a[i] * b[j] - a[j] * b[i] for i, j in combinations(range(3), 2))
    return apply(X, wedge)


def increment_formula(points):
    d1, d2, d3 = (subtract(points[i + 1], points[i]) for i in range(3))
    b12, b13, b23 = determinant(d1, d2), determinant(d1, d3), determinant(d2, d3)
    return 6 * b12 + 3 * b13, 6 * b13 + 4 * b23, 6 * b23


def adjacent_formula(weights, alphas, betas):
    """The expected summed residual, read directly from equation (M7)."""
    f1, f2, f3, f4 = weights
    d1, d2, d3 = f1 - f2, f2 - f3, f3 - f4
    a1, a2, a3 = alphas
    b12, b13, b23 = betas
    return (
        -3 * (2 * d2 + d3) * a1 + 6 * d1 * a2 + 3 * d1 * a3 + 12 * (f3 + f4) * b12 + 12 * f3 * b13,
        -6 * d3 * a1 - 4 * d3 * a2 + (6 * d1 + 4 * d2) * a3
        + 8 * d3 * b12 + (16 * f2 + 8 * f3) * b13 + 16 * f2 * b23,
        -6 * d3 * (a2 + b12) + 6 * d2 * (a3 + b13) + (18 * f1 + 6 * f2) * b23,
    )


def check_field_and_residual():
    for i in range(4):
        for j in range(4):
            points = tuple((int(l == i), int(l == j)) for l in range(4))
            assert field(points) == increment_formula(points)

    # Select one f_l and one determinant coefficient at a time. Taking
    # two elementary planar vectors and setting the others to zero selects
    # the coefficient of the pair from gamma,d1,d2,d3. In the direct
    # calculation replace the actual slot, without using sorted increments
    # or (M6), and apply the transported adjacent matrix.
    for l, weights in enumerate(identity(4)):
        forward = identity(3)
        inverse = identity(3)
        for generator in YOUNG[:l]:
            forward = multiply(generator, forward)
            inverse = multiply(inverse, generator)
        transport = multiply(multiply(forward, FIRST_TRANSPORT), inverse)
        for i, j in combinations(range(4), 2):
            vectors = [(0, 0)] * 4
            vectors[i], vectors[j] = (1, 0), (0, 1)
            gamma, d1, d2, d3 = vectors
            points = [(0, 0), d1, add(d1, d2), add(add(d1, d2), d3)]
            zero_replacement = points.copy()
            tilted_replacement = points.copy()
            zero_replacement[l] = (0, 0)
            tilted_replacement[l] = gamma
            direct = subtract(
                field(zero_replacement),
                apply(transport, field(tilted_replacement)),
            )
            alphas = tuple(determinant(gamma, d) for d in (d1, d2, d3))
            betas = tuple(determinant(a, b) for a, b in combinations((d1, d2, d3), 2))
            assert direct == adjacent_formula(weights, alphas, betas)

    # Exact known failure of the unregrouped pointwise conditional bound.
    original = field(tuple((x, x * x) for x in (1, 2, 3, 4)))
    replacement = field(tuple((x, x * x) for x in (0, 2, 3, 4)))
    witness = subtract(original, apply(FIRST_TRANSPORT, replacement))
    assert witness == (-42, -36, 36)
    return tuple(map(int, witness))


def main():
    check_young()
    witness = check_field_and_residual()
    print(json.dumps({
        "shape": [2, 2, 2],
        "young_and_helmert_generators": 3,
        "cubic_field_coefficients_per_coordinate": 16,
        "adjacent_residual_coefficients_per_coordinate": 24,
        "unregrouped_conditional_residual_limit": witness,
    }, indent=2))


if __name__ == "__main__":
    main()
