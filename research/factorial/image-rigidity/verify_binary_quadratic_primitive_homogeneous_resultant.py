"""Check the primitive homogeneous resultant at starts 1 through 5.

Requires SymPy. The all-start proof uses the moment recurrence and
multiplicativity of the resultant.
"""

from math import comb, factorial

from sympy import zeros


def monomials(degree):
    for i in range(degree + 1):
        for j in range(degree - i + 1):
            yield i, j, degree - i - j


def primitive_coefficients(degree):
    coefficients = {}
    for i, j, k in monomials(degree):
        numerator = factorial(2 * i + j) * factorial(2 * k + j)
        denominator = factorial(degree) * factorial(i) * factorial(j) * factorial(k)
        coefficient, remainder = divmod(numerator, denominator)
        assert remainder == 0
        coefficients[i, j, k] = coefficient
    return coefficients


def macaulay_resultant(start):
    degrees = (start, start + 1, start + 2)
    terms = list(monomials(sum(degrees) - 2))
    columns = {term: index for index, term in enumerate(terms)}
    forms = [primitive_coefficients(degree) for degree in degrees]
    matrix = zeros(len(terms))
    nonreduced = []
    for row, term in enumerate(terms):
        divisors = [axis for axis in range(3) if term[axis] >= degrees[axis]]
        if len(divisors) > 1:
            nonreduced.append(row)
        axis = divisors[0]
        shift = list(term)
        shift[axis] -= degrees[axis]
        for exponent, coefficient in forms[axis].items():
            shifted = tuple(power + offset for power, offset in zip(exponent, shift))
            matrix[row, columns[shifted]] = coefficient
    numerator = matrix.det(method="domain-ge")
    denominator = matrix.extract(nonreduced, nonreduced).det(method="domain-ge")
    quotient, remainder = divmod(numerator, denominator)
    assert remainder == 0
    return abs(quotient)


def expected_resultant(n):
    return (
        2 ** (4 * n**2 + 6 * n)
        * comb(2 * n + 1, n) ** ((n - 1) * (n + 2))
        * comb(2 * n + 3, n + 1) ** (n**2)
    )


if __name__ == "__main__":
    for n in range(1, 6):
        assert macaulay_resultant(n) == expected_resultant(n)
        print(f"Exact Macaulay resultant agrees at start {n}", flush=True)
