#!/usr/bin/env python3
"""Replay the exact coefficient witnesses for the tangent double-conic exclusion."""

from functools import cache
from itertools import product

from sympy import Matrix, Poly, Rational, expand, integrate, linear_eq_to_matrix, symbols


x, r = symbols("x r")
a, b, d, e, g, h = symbols("a b d e g h")
DEGREES = (5, 5, 4)
H0 = Matrix([r**2, r, r**2])
RADIAL = Matrix([a + b*r, d + e*r, g + h*r])


def determinant_piece(pieces, drop):
    """At z=1, use columns of homogeneous degrees and the r and x derivatives."""
    degree_columns = [
        Matrix([(degree - i)*entry for degree, entry in zip(DEGREES, piece)])
        for i, piece in enumerate(pieces)
    ]
    r_columns = [piece.diff(r) for piece in pieces]
    x_columns = [piece.diff(x) for piece in pieces]
    return expand(sum(
        Matrix.hstack(degree_columns[i], r_columns[j], x_columns[k]).det()
        for i, j, k in product(range(len(pieces)), repeat=3)
        if i + j + k == drop
    ))


def remove_radial(polynomial):
    return expand(polynomial - polynomial.subs(x, 0))


def coefficient_vector(polynomial, keys):
    polynomial = Poly(polynomial, x, r)
    return Matrix([polynomial.coeff_monomial(x**i*r**j) for i, j in keys])


def monomials(degree):
    # Original variables are x, y=r-x^2 at z=1.
    return [
        expand(x**i*(r - x**2)**j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    ]


@cache
def linear_system(drop):
    bases = [monomials(degree - drop) for degree in DEGREES]
    cofactors = (6*r**2, -2*r**3, -5*r**2)
    columns = [
        remove_radial(cofactor*term)
        for cofactor, basis in zip(cofactors, bases)
        for term in basis
    ]
    keys = sorted(
        {key for column in columns for key in Poly(column, x, r).monoms()},
        reverse=True,
    )
    matrix = Matrix.hstack(*(coefficient_vector(column, keys) for column in columns))
    pivots = matrix.rref()[1]
    rows = matrix[:, pivots].T.rref()[1]
    inverse = matrix.extract(rows, pivots).inv()
    return bases, keys, matrix, pivots, rows, inverse


def solve_piece(pieces, prefix):
    """Parametrize candidates and return necessary compatibility polynomials."""
    drop = len(pieces)
    # Integrate before removing the radial terms. In particular an x-independent
    # part of the determinant itself would give a constraint linear in x.
    known = remove_radial(integrate(determinant_piece(pieces, drop), x))
    bases, keys, matrix, pivots, rows, inverse = linear_system(drop)
    vector = coefficient_vector(known, keys)
    coefficients = Matrix.zeros(matrix.cols, 1)
    for pivot, value in zip(pivots, -inverse*vector.extract(rows, [0]), strict=True):
        coefficients[pivot] = value
    compatibility = list(matrix*coefficients + vector)
    compatibility.extend(
        value for key, value in Poly(known, x, r).terms() if key not in keys
    )

    kernel = matrix.nullspace()
    free = symbols(f"{prefix}:{len(kernel)}")
    for generator, parameter in zip(kernel, free, strict=True):
        coefficients += parameter*generator

    piece = []
    start = 0
    for basis in bases:
        piece.append(expand(sum(
            term*coefficients[start + i] for i, term in enumerate(basis)
        )))
        start += len(basis)
    return Matrix(piece), free, [expand(value) for value in compatibility]


def witness(conditions, expression):
    assert expand(expression) in conditions


def specialize(pieces, substitution):
    return [
        piece.subs(substitution, simultaneous=True).applyfunc(expand)
        for piece in pieces
    ]


def forced_solution(conditions, variables):
    """Solve selected rows with an invertible constant-determinant matrix."""
    matrix, vector = linear_eq_to_matrix(conditions, variables)
    rows = matrix.T.rref()[1]
    square = matrix.extract(rows, range(len(variables)))
    determinant = expand(square.det())
    assert determinant.is_Rational and determinant != 0
    values = square.inv()*vector.extract(rows, [0])
    substitution = {var: expand(value) for var, value in zip(variables, values)}
    assert all(expand(value.subs(substitution, simultaneous=True)) == 0
               for value in conditions)
    return substitution


def upper_correction(U):
    return RADIAL + Matrix([r*U, 3*U, 0])


def check_upper_corrections():
    aa, bb, alpha, beta = symbols("aa bb alpha beta")
    y = r - x**2
    U = aa*x*y + bb*y**2 + alpha*x + beta*y
    primitive = (
        6*r**2*U**2 + r*(15*g - 24*a)*U
        + r**2*((-12*b + 8*d + 5*h) + 4*e*r)*U
    )
    assert expand(determinant_piece([H0, upper_correction(U)], 2)
                  - primitive.diff(x)) == 0

    pieces = [H0, upper_correction(U.subs(bb, 0)).subs(g, 8*a/5)]
    H2, _, conditions = solve_piece(pieces, "f")
    assert not any(conditions)
    _, _, conditions = solve_piece([*pieces, H2], "l")
    witness(conditions, 10*aa**3)
    conditions = [expand(value.subs(aa, 0)) for value in conditions]
    witness(conditions, 12*a*beta**2)
    witness([expand(value.subs(beta, 0)) for value in conditions], 12*a*alpha**2)

    for correction in (x, y):
        pieces = [H0, upper_correction(correction).subs({a: 0, g: 0})]
        H2, f, conditions = solve_piece(pieces, "f")
        assert not any(conditions)
        H3, _, conditions = solve_piece([*pieces, H2], "l")
        substitution = forced_solution(conditions, [f[i] for i in (6, 8, 9, 10, 11)])
        _, _, conditions = solve_piece(specialize([*pieces, H2, H3], substitution), "m")
        witness(conditions, Rational(35, 3))


def check_radial_correction():
    pieces = [H0, RADIAL]
    H2, f, conditions = solve_piece(pieces, "f")
    assert not any(conditions)
    H3, ell, drop3 = solve_piece([*pieces, H2], "l")
    _, _, drop4 = solve_piece([*pieces, H2, H3], "m")

    for value, parameter, substitute in (
        (5*f[8]**2/12, f[8], 0),
        (5*f[10]**2/12, f[10], 0),
        (2*(f[3] - f[5])**2/3, f[3], f[5]),
        (5*(f[11] - f[7])**2/12, f[11], f[7]),
    ):
        witness(drop4, value)
        drop3 = [expand(expr.subs(parameter, substitute)) for expr in drop3]
        drop4 = [expand(expr.subs(parameter, substitute)) for expr in drop4]
        H2 = H2.subs(parameter, substitute)

    theta = 8*a - 5*g
    witness(drop3, f[9]*theta/2)
    witness(drop4, (5*f[9]**2 - 6*theta*ell[5])/12)
    H2 = H2.subs(f[9], 0).applyfunc(expand)
    assert H2 == Matrix([
        f[0] + f[1]*r + f[4]*r*x/3,
        f[2] + f[4]*x + f[5]*r,
        f[6] + f[7]*r,
    ])

    H3, ell, _ = solve_piece([*pieces, H2], "l")
    H4, _, _ = solve_piece([*pieces, H2, H3], "m")
    pieces = [*pieces, H2, H3, H4]
    degree_five = Poly(determinant_piece(pieces, 6), x, r)
    assert degree_five.coeff_monomial(x**3) == 5*ell[5]**2/3
    assert expand(degree_five.coeff_monomial(x**2*r)
                  - 10*f[4]*(3*e*ell[5] - f[4]**2)/9) == 0

    pieces = specialize(pieces, {ell[5]: 0, f[4]: 0})
    assert determinant_piece(pieces, 11) == 0


if __name__ == "__main__":
    check_upper_corrections()
    check_radial_correction()
    print("Tangent double conic: exact upper and lower determinant witnesses verified.")
