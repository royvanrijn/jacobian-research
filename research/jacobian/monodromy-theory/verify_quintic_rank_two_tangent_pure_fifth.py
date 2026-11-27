#!/usr/bin/env python3
"""Replay the coefficient witnesses excluding the tangent pure-fifth conic pencil."""

from functools import cache

from sympy import Matrix, Poly, Rational, expand, integrate, symbols

from verify_quintic_rank_two_tangent_double_conic import (
    DEGREES,
    RADIAL,
    b,
    e,
    coefficient_vector,
    determinant_piece,
    monomials,
    r,
    remove_radial,
    witness,
    x,
)


H0 = Matrix([1, r, r**2])


@cache
def linear_system(drop):
    bases = [monomials(degree - drop) for degree in DEGREES]
    cofactors = (6*r**2, -10*r, 5)
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
    drop = len(pieces)
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


def upper_pieces(U):
    return [H0, RADIAL + Matrix([U, 3*r*U/5, 0])]


def through_degree_seven(U):
    pieces = upper_pieces(U)
    H2, f, conditions = solve_piece(pieces, "f")
    assert not any(conditions)
    H3, _, drop3 = solve_piece([*pieces, H2], "l")
    _, _, drop4 = solve_piece([*pieces, H2, H3], "m")
    return pieces, H2, f, drop3, drop4


def check_upper_correction():
    A, B, alpha, beta = symbols("A B alpha beta")
    y = r - x**2
    U = A*x*y + B*y**2 + alpha*x + beta*y
    _, _, conditions = solve_piece(upper_pieces(U), "f")
    witness(conditions, 24*B**2/5)
    witness([expand(value.subs(B, 0)) for value in conditions], 12*A**2/5)

    _, _, _, degree_eight, _ = through_degree_seven(y)
    witness(degree_eight, Rational(-2, 25))

    _, _, f, degree_eight, degree_seven = through_degree_seven(x)
    witness(degree_eight, 2*(20*b*e - 50*b*f[5] - 15*f[1] + 4)/25)
    witness(degree_seven, 4*b*f[4]/15)
    quadratic = 625*f[4]**2 + 75*f[4] + 3
    second = 375*f[1]*f[4] + 1250*f[4]**2 + 3
    witness(degree_seven, -2*quadratic/375)
    witness(degree_seven, 2*second/375)
    assert expand(second.subs(f[1], Rational(4, 15)) - 2*quadratic) == -50*f[4] - 3
    assert quadratic.subs(f[4], Rational(-3, 50)) == Rational(3, 4)


def check_radial_correction():
    pieces, H2, f, _, drop4 = through_degree_seven(0)
    witness(drop4, -10*f[4]**2/3)
    pieces.append(H2.subs(f[4], 0).applyfunc(expand))
    H3, _, _ = solve_piece(pieces, "l")
    pieces.append(H3)
    H4, _, _ = solve_piece(pieces, "m")
    pieces.append(H4)
    degree_five = Poly(determinant_piece(pieces, 6), x, r)
    assert degree_five.coeff_monomial(x**2*r**2) == 10*f[5]**3/9
    pieces = [piece.subs(f[5], 0).applyfunc(expand) for piece in pieces]
    assert all(piece.diff(x).is_zero_matrix for piece in pieces)
    assert determinant_piece(pieces, 11) == 0


if __name__ == "__main__":
    check_upper_correction()
    check_radial_correction()
    print("Tangent pure fifth: exact determinant witnesses verified.")
