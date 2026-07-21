#!/usr/bin/env python3
"""Dependency-free sparse-polynomial audit of the explicit C14 model."""
from fractions import Fraction
from itertools import permutations


Poly = dict[tuple[int, ...], Fraction]


def clean(poly: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def const(value, variables):
    value = Fraction(value)
    return {} if not value else {(0,) * variables: value}


def var(index, variables):
    exponent = [0] * variables
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
    return clean(result)


def scale(polynomial, scalar):
    scalar = Fraction(scalar)
    return clean({monomial: scalar * coefficient for monomial, coefficient in polynomial.items()})


def mul(left, right):
    result = {}
    for lm, lc in left.items():
        for rm, rc in right.items():
            monomial = tuple(a + b for a, b in zip(lm, rm))
            result[monomial] = result.get(monomial, Fraction(0)) + lc * rc
    return clean(result)


def power(polynomial, exponent):
    variables = len(next(iter(polynomial))) if polynomial else 1
    result = const(1, variables)
    for _ in range(exponent):
        result = mul(result, polynomial)
    return result


def derivative(polynomial, index):
    result = {}
    for monomial, coefficient in polynomial.items():
        if monomial[index]:
            exponent = list(monomial)
            factor = exponent[index]
            exponent[index] -= 1
            result[tuple(exponent)] = coefficient * factor
    return clean(result)


def divide_monomial(polynomial, index, exponent):
    assert all(monomial[index] >= exponent for monomial in polynomial)
    result = {}
    for monomial, coefficient in polynomial.items():
        reduced = list(monomial)
        reduced[index] -= exponent
        result[tuple(reduced)] = coefficient
    return result


def evaluate(polynomial, values):
    total = Fraction(0)
    for monomial, coefficient in polynomial.items():
        term = coefficient
        for exponent, value in zip(monomial, values):
            term *= Fraction(value) ** exponent
        total += term
    return total


def det3(matrix):
    result = {}
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = const(-1 if inversions % 2 else 1, 3)
        for row, column in enumerate(permutation):
            term = mul(term, matrix[row][column])
        result = add(result, term)
    return result


# Rebuild the map without importing its implementation.
x, y, z = (var(index, 3) for index in range(3))
one3 = const(1, 3)
u = add(one3, scale(mul(x, y), 3))
gamma = add(one3, scale(mul(x, y), -4), scale(mul(power(x, 2), z), -1))
num_a = add(scale(u, 2), power(u, 2), scale(mul(power(u, 4), power(gamma, 2)), -3))
num_b = add(one3, u, scale(mul(power(u, 3), power(gamma, 2)), -2))
A_map = divide_monomial(num_a, 0, 2)
B_map = divide_monomial(num_b, 0, 1)
C_map = mul(x, gamma)
mapping = (A_map, B_map, C_map)
jacobian = [[derivative(component, index) for index in range(3)] for component in mapping]
assert det3(jacobian) == const(-6, 3)
assert tuple(evaluate(component, (1, 0, 0)) for component in mapping) == (0, 0, 1)
assert tuple(evaluate(component, (-1, 0, 2)) for component in mapping) == (0, 0, 1)

# Source marking satisfies the quartic and its derivative relation.
W_map = mul(u, gamma)
E_map = add(
    power(W_map, 2),
    scale(power(W_map, 4), -1),
    scale(mul(mul(B_map, C_map), W_map), -2),
    mul(A_map, power(C_map, 2)),
)
dE_map = add(
    scale(W_map, 2),
    scale(power(W_map, 3), -4),
    scale(mul(B_map, C_map), -2),
    scale(gamma, 2),
)
assert not E_map and not dE_map

# Independently evaluate the simplified quartic discriminant formula in
# Q[A,B,C] for a=-1,b=0,c=1,d=-2BC,e=AC^2.
A, B, C = (var(index, 3) for index in range(3))
d = scale(mul(B, C), -2)
e = mul(A, power(C, 2))
discriminant = add(
    scale(power(e, 3), -256),
    scale(power(e, 2), -128),
    scale(mul(mul(power(d, 2), e), const(1, 3)), 144),
    scale(power(d, 4), -27),
    scale(e, -16),
    scale(power(d, 2), 4),
)
Q4 = add(
    A,
    scale(power(B, 2), -1),
    mul(power(C, 2), add(scale(power(B, 4), 27), scale(mul(A, power(B, 2)), -36), scale(power(A, 2), 8))),
    scale(mul(power(A, 3), power(C, 4)), 16),
)
assert discriminant == scale(mul(power(C, 2), Q4), -16)

# Repeated-root normalization is checked coefficientwise by direct rational
# evaluation at several values; both sides have total r-degree at most six,
# so seven points certify the identity.
for r in range(-3, 4):
    rr = Fraction(r)
    for w in range(-2, 3):
        ww = Fraction(w)
        left = ww**2 - ww**4 - 2*(rr-2*rr**3)*ww + rr**2 - 3*rr**4
        right = -(ww-rr)**2 * (ww**2 + 2*rr*ww + 3*rr**2 - 1)
        assert left == right

# C=0 chart identities and the omitted node factorization.
assert (3*Fraction(2)+2)**2 - 3*(Fraction(2)+1)*(3*Fraction(2)+1) == 1
for w in (Fraction(-2), Fraction(-1, 2), Fraction(0), Fraction(3, 2)):
    assert w*w - w**4 - Fraction(1, 4) == -(w*w-Fraction(1, 2))**2

print("PASS: clean-room quartic map, collision, incidence, and det=-6")
print("PASS: clean-room quartic discriminant and repeated-root normalization")
print("PASS: clean-room C=0 and omitted-node fiber identities")
