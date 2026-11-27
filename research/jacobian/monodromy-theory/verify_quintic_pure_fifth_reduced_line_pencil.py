#!/usr/bin/env python3
"""Replay the pure-fifth reduced-line determinant coefficients."""

from functools import cache

from sympy import Matrix, Poly, expand, symbols

from verify_quintic_rank_two_tangent_double_conic import determinant_piece, r, x


DEGREES = (5, 5, 4)
TOP = Matrix([1, r, r])
V = Matrix([0, 1, 1])
W = Matrix([1, r / 5, 0])
y = r / x
a, d, g, a2, c2, b1, d1, b2, d2, t0, t1 = symbols("a d g a2 c2 b1 d1 b2 d2 t0 t1")
FIRST = (
    Matrix([a + t0*r**2, d + t1*r**2, g])
    + (a2*x**2 + c2*y**2)*V
    + (b1*x + d1*y + b2*x**2 + d2*y**2)*W
)


def coefficients(expression):
    """Collect powers of x, allowing negative Laurent exponents."""
    modes = {}
    for term in expand(expression).as_ordered_terms():
        mode = int(term.as_powers_dict().get(x, 0))
        modes[mode] = modes.get(mode, 0) + term/x**mode
    return {
        (mode, exponent[0]): value
        for mode, coefficient in modes.items()
        for exponent, value in Poly(expand(coefficient), r).terms()
        if value != 0
    }


def determinant(pieces, drop):
    return coefficients(x * determinant_piece(pieces, drop))


@cache
def system(drop):
    bases = [
        [x**i * y**j for i in range(degree-drop+1) for j in range(degree-drop+1-i)]
        for degree in DEGREES
    ]
    columns = [
        coefficients(factor*term)
        for factor, basis in zip((r, -5, 5), bases, strict=True)
        for term in basis
    ]
    keys = sorted({key for col in columns for key in col if key[0]}, reverse=True)
    matrix = Matrix([[col.get(key, 0) for col in columns] for key in keys])
    pivots = matrix.rref()[1]
    rows = matrix[:, pivots].T.rref()[1]
    return bases, keys, matrix, pivots, rows, matrix.extract(rows, pivots).inv()


def next_piece(pieces, prefix):
    drop = len(pieces)
    known = determinant(pieces, drop)
    primitive = {key: val/key[0] for key, val in known.items() if key[0]}
    bases, keys, matrix, pivots, rows, inverse = system(drop)
    vector = Matrix([primitive.get(key, 0) for key in keys])
    solution = Matrix.zeros(matrix.cols, 1)
    for pivot, value in zip(pivots, -inverse*vector.extract(rows, [0]), strict=True):
        solution[pivot] = value
    compatibility = dict(zip(keys, matrix*solution+vector, strict=True))
    compatibility.update({key: val for key, val in primitive.items() if key not in keys})
    compatibility.update({key: val for key, val in known.items() if not key[0]})
    kernel = matrix.nullspace()
    free = symbols(f"{prefix}:{len(kernel)}")
    for direction, parameter in zip(kernel, free, strict=True):
        solution += parameter*direction
    piece = []
    offset = 0
    for basis in bases:
        piece.append(expand(sum(term*solution[offset+i] for i, term in enumerate(basis))))
        offset += len(basis)
    return Matrix(piece), free, {key: expand(val) for key, val in compatibility.items()}


def witness(compatibility, key, value):
    assert expand(compatibility.get(key, 0)-value) == 0, (key, value, compatibility.get(key))


def check_first_piece():
    _, _, eq = next_piece([TOP, FIRST], "f")
    witness(eq, (4, 1), -2*b2**2/5)
    first = FIRST.subs({b2: 0, d2: 0})
    second, _, eq2 = next_piece([TOP, first], "f")
    _, _, eq3 = next_piece([TOP, first, second], "g")
    witness(eq2, (2, 2), -2*a2*t0)
    witness(eq3, (4, 0), -5*a2**2*t1)
    witness(eq3, (3, 1), 2*b1*(50*a2*t1-b1**2)/25)


def check_radial_fourth_terms():
    first = FIRST.subs({a2: 0, c2: 0, b1: 0, d1: 0, b2: 0, d2: 0})
    second, f, _ = next_piece([TOP, first], "f")
    third, _, eq3 = next_piece([TOP, first, second], "g")
    _, _, eq4 = next_piece([TOP, first, second, third], "h")
    witness(eq3, (1, 3), 4*f[5]*t0)
    witness(eq3, (1, 2), -2*(10*f[5]*t1+f[9]*t0))
    witness(eq3, (2, 2), -2*f[11]*t0)
    witness(eq4, (2, 2), 40*f[11]*t1**2)
    known = {t0: 0, f[5]: 0, f[11]: 0}
    assert expand(eq4[(1, 2)].subs(known)) == 40*f[9]*t1**2
    second = second.subs({f[i]: 0 for i in (4, 5, 7, 8, 9, 11)})
    third, g3, _ = next_piece([TOP, first, second], "g")
    fourth, _, eq4 = next_piece([TOP, first, second, third], "h")
    witness(eq4, (1, 2), -2*g3[6]*t0)
    eq5 = determinant([TOP, first, second, third, fourth], 5)
    witness(eq5, (1, 2), 40*g3[6]*t1**2)


def check_bilinear_reduction():
    first = FIRST.subs({t0: 0, t1: 0, b1: 0, d1: 0, b2: 0, d2: 0})
    second, f, _ = next_piece([TOP, first], "f")
    third, _, eq3 = next_piece([TOP, first, second], "g")
    _, _, eq4 = next_piece([TOP, first, second, third], "h")
    witness(eq3, (3, 0), 5*a2*f[5])
    witness(eq4, (-3, 3), -5*c2**2*f[5]+(4*a*c2+5*f[8])*f[4])
    witness(eq4, (1, 2), f[1]*f[5])
    witness(eq4, (3, 0), (4*a*a2+5*f[11])*f[5]-5*a2**2*f[4])
    witness(eq4, (2, 1), -a2*symbols("g1")-f[1]*f[11]-10*f[5]**2)
    second = second.subs({f[4]: 0, f[5]: 0})
    third, g3, _ = next_piece([TOP, first, second], "g")
    assert expand(second[0]-f[0]-f[1]*r) == 0
    assert expand(third[0]-g3[0]-g3[1]*r-f[1]*(a2*x**2+c2*y**2)) == 0


if __name__ == "__main__":
    check_first_piece()
    check_radial_fourth_terms()
    check_bilinear_reduction()
    print("Pure-fifth reduced-line pencil: exact Laurent modes and bilinear reduction verified.")
