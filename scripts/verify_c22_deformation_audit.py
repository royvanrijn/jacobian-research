"""Dependency-free regressions for the deformation-theoretic audit of C22.

This deliberately avoids invariant rings and computer algebra.  Coefficients
in the ordered conductor ribbon are represented by squarefree bit masks.
"""

from fractions import Fraction
from itertools import combinations
from math import comb


def boolean_add(left, right):
    answer = left.copy()
    for mask, coefficient in right.items():
        answer[mask] = answer.get(mask, Fraction(0)) + coefficient
        if answer[mask] == 0:
            del answer[mask]
    return answer


def boolean_scale(coefficient, value):
    return {
        mask: coefficient * entry
        for mask, entry in value.items()
        if coefficient * entry
    }


def boolean_multiply(left, right):
    answer = {}
    for left_mask, left_coefficient in left.items():
        for right_mask, right_coefficient in right.items():
            if left_mask & right_mask:
                continue
            mask = left_mask | right_mask
            answer[mask] = (
                answer.get(mask, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {mask: coefficient for mask, coefficient in answer.items() if coefficient}


def polynomial_add(left, right):
    size = max(len(left), len(right))
    zero = {}
    answer = []
    for degree in range(size):
        left_coefficient = left[degree] if degree < len(left) else zero
        right_coefficient = right[degree] if degree < len(right) else zero
        answer.append(boolean_add(left_coefficient, right_coefficient))
    while answer and not answer[-1]:
        answer.pop()
    return answer


def polynomial_scale(coefficient, value):
    return [boolean_scale(coefficient, entry) for entry in value]


def polynomial_multiply(left, right):
    if not left or not right:
        return []
    answer = [{} for _ in range(len(left) + len(right) - 1)]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            product = boolean_multiply(left_coefficient, right_coefficient)
            answer[left_degree + right_degree] = boolean_add(
                answer[left_degree + right_degree], product
            )
    while answer and not answer[-1]:
        answer.pop()
    return answer


def polynomial_power(value, exponent):
    answer = [{0: Fraction(1)}]
    for _ in range(exponent):
        answer = polynomial_multiply(answer, value)
    return answer


def ordered_norms(root_count):
    """Return U and V for distinct rational roots in the Boolean ribbon."""
    U = [{0: Fraction(1)}]
    V = [{0: Fraction(1)}]
    for index in range(root_count):
        root = Fraction(index + 1)
        epsilon = {1 << index: Fraction(1)}
        q = [{0: -root}, {0: Fraction(1)}]
        q_squared = polynomial_power(q, 2)
        q_cubed = polynomial_multiply(q_squared, q)
        vertical_square = polynomial_add(q_squared, [epsilon])
        vertical_cube = polynomial_add(
            q_cubed,
            polynomial_scale(Fraction(3, 2), polynomial_multiply([epsilon], q)),
        )
        V = polynomial_multiply(V, vertical_square)
        U = polynomial_multiply(U, vertical_cube)
    return U, V


def determinant_mod(matrix, prime):
    matrix = [[entry % prime for entry in row] for row in matrix]
    size = len(matrix)
    answer = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            answer = -answer
        pivot_value = matrix[column][column]
        answer = answer * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        for row in range(column + 1, size):
            factor = matrix[row][column] * inverse % prime
            for entry in range(column, size):
                matrix[row][entry] = (
                    matrix[row][entry] - factor * matrix[column][entry]
                ) % prime
    return answer % prime


def compound_matrix(matrix, degree, prime):
    if degree == 0:
        return [[1]]
    subsets = list(combinations(range(len(matrix)), degree))
    return [
        [
            determinant_mod(
                [[matrix[row][column] for column in columns] for row in rows],
                prime,
            )
            for columns in subsets
        ]
        for rows in subsets
    ]


def normalized_compound_matrix(matrix, degree, roots, prime):
    """Evaluation matrix of the Schur/divided-difference normal symbols."""
    compound = compound_matrix(matrix, degree, prime)
    if degree < 2:
        return compound
    subsets = list(combinations(range(len(matrix)), degree))
    for column, subset in enumerate(subsets):
        internal_vandermonde = 1
        for left_index, left in enumerate(subset):
            for right in subset[left_index + 1 :]:
                internal_vandermonde = (
                    internal_vandermonde * (roots[right] - roots[left])
                ) % prime
        inverse = pow(internal_vandermonde, prime - 2, prime)
        for row in range(len(compound)):
            compound[row][column] = compound[row][column] * inverse % prime
    return compound


def evaluate(coefficients, value):
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def derivative(coefficients):
    return [degree * coefficient for degree, coefficient in enumerate(coefficients)][1:]


def divided_difference(coefficients, left, right):
    if left == right:
        return evaluate(derivative(coefficients), left)
    return (evaluate(coefficients, left) - evaluate(coefficients, right)) / (
        left - right
    )


# The unique one-root transverse jet and all of its ordered tensor products.
for k in range(1, 8):
    U, V = ordered_norms(k)
    assert polynomial_power(U, 2) == polynomial_power(V, 3)

# Exact divided differences extend evaluation through a collision and satisfy
# the product rule used by the straightening argument.
f = [Fraction(2), Fraction(-3), Fraction(5), Fraction(7)]
g = [Fraction(-1), Fraction(4), Fraction(0), Fraction(3)]
fg = [
    sum(
        (f[i] * g[degree - i] for i in range(len(f)) if 0 <= degree - i < len(g)),
        Fraction(0),
    )
    for degree in range(len(f) + len(g) - 1)
]
for left, right in ((Fraction(2), Fraction(5)), (Fraction(3), Fraction(3))):
    assert divided_difference(fg, left, right) == (
        evaluate(f, left) * divided_difference(g, left, right)
        + evaluate(g, right) * divided_difference(f, left, right)
    )

# The vertical-degree-d normal symbols are Schur polynomials, obtained by
# dividing each alternant by its internal Vandermonde. The normalized
# compound determinant checks their number on a separated specialization.
prime = 1_000_003
for k in range(1, 9):
    roots = list(range(1, k + 1))
    vandermonde = [[pow(root, exponent, prime) for root in roots] for exponent in range(k)]
    vandermonde_determinant = determinant_mod(vandermonde, prime)
    assert vandermonde_determinant
    ranks = []
    for degree in range(k + 1):
        compound = normalized_compound_matrix(vandermonde, degree, roots, prime)
        compound_determinant = determinant_mod(compound, prime)
        exponent = 0 if degree in (0, k) else comb(k - 2, degree - 1)
        expected = pow(vandermonde_determinant, exponent, prime)
        assert compound_determinant == expected
        ranks.append(len(compound))
    assert ranks == [comb(k, degree) for degree in range(k + 1)]
    assert sum(ranks) == 2**k

print("PASS: ordered cusp jets satisfy U^2=V^3 through seven tensor factors")
print("PASS: exact divided differences satisfy the confluent product rule")
print("PASS: normalized compound bases have ranks binomial(k,d) through k=8")
print("PASS: the collision filtration has total rank 2^k through k=8")
