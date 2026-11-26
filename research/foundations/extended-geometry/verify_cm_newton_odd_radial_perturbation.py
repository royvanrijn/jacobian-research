#!/usr/bin/env python3
"""Replay the exact unmarked moments and coherent reconstruction for g_b."""

from math import comb, factorial, prod

from sage.all import GF, QQ, ZZ, PolynomialRing, identity_matrix, matrix


def height_integral(height, count):
    assert height >= 0 and height % 2 == 0
    return QQ(2**count * factorial(count)) / prod(
        height + 1 + 2*j for j in range(count + 1)
    )


def verify_two_moments(ring):
    X = ring.gen()
    moments = []
    for k in (2, 4):
        moments.append(sum(
            comb(k, h) * comb(4*h, h//2)
            * X**((k-h)//2) * height_integral(17*k-h, h)
            for h in range(0, k + 1, 2)
        ))
    second, fourth = moments
    assert second == X/35 + QQ(64)/42735
    assert fourth == X**2/69 + QQ(128)*X/109411 + QQ(3072)/103831039
    quotient, remainder = fourth.quo_rem(second)
    assert quotient == 35*X/69 + QQ(5754560)/400772493
    assert remainder == QQ(3753620480)/464386710041397
    print("Orders 2 and 4: exact nonzero Bézout remainder")


def normalized_row(k, ring):
    X = ring.gen()
    coherent = [
        sum(
            comb(k, h) * (-1)**h * comb(4*h+2, (h-j)//2)
            * X**((k-h)//2) * height_integral(17*k-h+12, h)
            for h in range(k % 2, k + 1, 2)
            if h >= j and (h-j) % 2 == 0
        )
        for j in range(2)
    ]
    corrections = [
        sum(
            comb(k, h) * (-1)**h * comb(4*h, (h-j)//2)
            * X**((k-h-1)//2) * height_integral(17*k-h+11, h)
            for h in range(1-k % 2, k, 2)
            if h >= j and (h-j) % 2 == 0
        )
        for j in range(6)
    ]
    return [ring(entry) for entry in coherent + corrections]


def odd_row(k, ring):
    row = normalized_row(k, ring)
    return row[1], row[2::2]


def verify_coherent_reconstruction(ring):
    X = ring.gen()
    solutions = []
    for k in (1, 3, 5):
        coherent, columns = odd_row(k, ring)
        j = len(solutions)
        assert columns[j].is_constant() and columns[j] and not any(columns[j+1:])
        residual = coherent + sum(col*value for col, value in zip(columns, solutions))
        solutions.append(-residual / QQ(columns[j]))
    assert solutions == [
        QQ(2)/31, 427*X/93-QQ(204)/2077,
        -707*X**2/6+551306*X/6231-QQ(360474)/213931,
    ]

    finite = PolynomialRing(GF(31), "X")
    x = finite.gen()
    reductions = []
    for k in (7, 9):
        coherent, columns = odd_row(k, ring)
        residual = coherent + sum(col*value for col, value in zip(columns, solutions, strict=True))
        cleared = (residual * residual.denominator()).change_ring(ZZ)
        integral = cleared.parent()([c // cleared.content() for c in cleared.list()])
        if integral.leading_coefficient() > 0:
            integral = -integral
        reductions.append(finite(integral))
        assert reductions[-1].degree() == residual.degree()
    cubic, quartic = reductions
    assert cubic == 8*x**3 + 27*x + 15
    assert quartic == 25*x**4 + 5*x**2 + 20*x
    assert quartic % cubic == 2*x**2 + 8*x
    assert cubic % (quartic % cubic) == 15
    print("Rows 1, 3, 5, 7, 9: coherent component separated exactly")


def verify_polynomial_reconstruction(ring):
    shifts = [0, 1, 2, 3, 4, 5, 6, 7, 9]
    N = matrix(ring, [normalized_row(k, ring) for k in shifts])
    coherent, corrections = N[:, :2], N[:, 2:]
    early = list(range(1, 7))
    block = corrections.matrix_from_rows(early)
    diagonal = [(j+1)*(-1)**j*height_integral(16*j+28, j) for j in range(6)]
    assert block.diagonal() == diagonal
    assert not any(block[i, j] for i in range(6) for j in range(i+1, 6))
    inverse = block.inverse().change_ring(ring)

    rows = identity_matrix(ring, len(shifts))
    a0 = 13*rows[0]
    early_coherent = coherent.matrix_from_rows(early)
    initial = inverse * (
        rows.matrix_from_rows(early) - early_coherent[:, :1]*matrix(ring, [a0])
    )
    r7, r9 = [
        rows[i] - coherent[i, 0]*a0 - corrections[i]*initial
        for i in (7, 8)
    ]
    s7, s9 = [row.dot_product(N.column(1)) for row in (r7, r9)]
    gcd, alpha, beta = s7.xgcd(s9)
    assert gcd == 1 and alpha*s7 + beta*s9 == 1
    a1 = alpha*r7 + beta*r9
    weighted_c = initial - inverse*early_coherent[:, 1:2]*matrix(ring, [a1])
    left_inverse = matrix(ring, [a0, a1, *weighted_c.rows()])
    assert left_inverse*N == identity_matrix(ring, 8)

    smith_rows = matrix(ring, [*left_inverse.rows(), -s9*r7+s7*r9])
    assert smith_rows.det() == 13/prod(diagonal)
    assert smith_rows*N == identity_matrix(ring, 8).stack(matrix(ring, 1, 8))
    print("Verified the polynomial inverse and unimodular Smith transformation")


if __name__ == "__main__":
    ring = PolynomialRing(QQ, "X")
    verify_two_moments(ring)
    verify_coherent_reconstruction(ring)
    verify_polynomial_reconstruction(ring)
