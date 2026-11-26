#!/usr/bin/env python3
"""Check finite critical orbits and the two sphere response matrices."""

from collections import defaultdict
from math import comb

from sage.all import GF, QQ, PolynomialRing, matrix
from sympy import I, Matrix, Poly, Rational, diff, expand, gcd, resultant, sqrt, symbols


def finite_critical_orbits():
    z, h = symbols("z h")
    v = z-1
    P = (1-z*z)**4
    cases = [
        (z*v**4, v*v, 11, 11, 10, 1, 2,
         [6627, 28341, -5474, -46, -6173, 1653, 156, 4]),
        (2*z*v*v*(4+5*v*v-v**3), 1+2*v*v, 73, 331, 11, -4, 110,
         [3654, -66528, 344629, -787520, 904709, -459752, 88183, 264, 9]),
    ]
    for V, W, irreducibility_prime, check_prime, power, scale, determinant, coefficients in cases:
        Phi = P+h*V+h*h*W
        a = 18*P+17*h*V+16*h*h*W
        b = z*diff(Phi, z)-Phi
        R = Poly.from_list(coefficients, z).as_expr()
        assert expand(resultant(a, b, h) - scale*(z+1)**3*v**power*R) == 0
        assert Poly(gcd(a.subs(z, -1), b.subs(z, -1)), h).monic().as_expr() == h
        if power == 10:
            assert a.subs(z, 1) == b.subs(z, 1) == 0 and Phi.subs(z, 1) == 0
            assert expand(V**2-4*P*W + v**6*(z*z+5*z+2)*(3*z*z+3*z+2)) == 0
        else:
            assert Poly(gcd(a.subs(z, 1), b.subs(z, 1)), h).monic().as_expr() == h*h
            assert expand(V**2-4*P*W - 4*(z-3)*v**9*(z*z-8*z-1)) == 0

        ring = PolynomialRing(GF(irreducibility_prime), "Z")
        modulus = ring(list(reversed(coefficients)))
        assert modulus.degree() == len(coefficients)-1 and modulus.is_irreducible()

        eliminated = Poly(16*W*b - (z*diff(W, z)-W)*a, h)
        assert eliminated.degree() == 1

        def quotient_solution(field):
            ring = PolynomialRing(field, "Z")
            modulus = ring(list(reversed(coefficients))).monic()
            assert modulus.degree() == len(coefficients)-1 and modulus.is_squarefree()
            quotient = ring.quotient(modulus, "q")
            Z = quotient.gen()

            def evaluate(polynomial, H=None):
                return sum(
                    quotient(field(int(coefficient)))*Z**j*(H**k if k else 1)
                    for (j, k), coefficient in Poly(polynomial, z, h).terms()
                )

            assert all(x.is_unit() for x in (Z, Z*Z-1, evaluate(W), evaluate(P)))
            H = -evaluate(eliminated.nth(0))/evaluate(eliminated.nth(1))
            assert H.is_unit()
            assert not evaluate(a, H) and not evaluate(b, H)
            return Z, H, evaluate

        # Check both zero identities in characteristic zero before the modular checks.
        quotient_solution(QQ)
        field = GF(check_prime)
        Z, H, evaluate = quotient_solution(field)
        jacobian = diff(a, z)*diff(b, h)-diff(a, h)*diff(b, z)
        assert evaluate(jacobian, H).is_unit()
        critical_value = evaluate(Phi, H)/(Z*H**18)
        assert critical_value.is_unit()
        degree = len(coefficients)-1
        powers = matrix(field, [
            [(critical_value**j).lift()[k] for k in range(degree)]
            for j in range(degree)
        ])
        assert powers.det() == determinant
        print(f"{degree} distinct nonzero finite Morse values in one Galois orbit"
              f" (power-basis determinant {determinant} modulo {check_prime})")


def finite_marker_rank():
    field = GF(1009)
    i, sqrt3 = field(469), field(149)
    assert i*i == -1 and sqrt3*sqrt3 == 3
    ring = PolynomialRing(field, "z")
    z = ring.gen()
    v = z-i
    tau = v/i
    top = (z*z+1)**4
    markers = [
        [(0, 0, 1), (2, 0, 2), (4, 0, 1)],
        [(1, 0, 1), (3, 0, 2), (5, 0, 1)],
    ] + [[(j, 1, 1)] for j in range(6)]

    phases = [
        (v**4, z-2*i, 349, 576),
        (-2/sqrt3*v*v*(4+5*tau*tau-tau**3), -4*i/3+2*z/3, 323, 52),
        (1+z**5, z, 3, None),
    ]
    for B, D, determinant, short_determinant in phases:
        phase = {(j-1, 0): c for j, c in top.dict().items()}
        phase.update({(j, 1): c for j, c in B.dict().items()})
        phase.update({(j, 2): c for j, c in D.dict().items()})
        phase[-1, 2] = -field.one()
        coefficients = {(0, 0): field.one()}
        rows = []
        for k in range(8):
            rows.append([
                sum(
                    weight*c/field(13-parity+18*k-deficit)
                    for j, parity, weight in column
                    for (angle, deficit), c in coefficients.items()
                    if angle == -j and deficit % 2 == parity
                )
                for column in markers
            ])
            if k < 7:
                product = defaultdict(field.zero)
                for (angle, deficit), c in coefficients.items():
                    for (step, loss), value in phase.items():
                        product[angle+step, deficit+loss] += c*value
                coefficients = product
        assert matrix(field, rows).det() == determinant
        if short_determinant is None:
            print(f"Generic eight-row witness: determinant {determinant} modulo 1009")
        else:
            assert matrix(field, [row[:5] for row in rows[:5]]).det() == short_determinant
            print(f"Verified marked ranks eight and five: determinants {determinant}, {short_determinant}")


def normalized_period_rows():
    z, a0, a1, c0, c1, c2 = symbols("z a0 a1 c0 c1 c2")
    inputs = (a0, a1, c0, c1, c2)
    C = c0+c1*z+c2*z*z
    g = -4*I*sqrt(3)
    numerator = (z+I)**2*(a0+a1*z) + (g+g*(z-I)/I)*C
    leading = expand(numerator.subs(z, I))
    resonant = expand(diff(numerator, z).subs(z, I)-4*leading/(3*I))
    rows = Matrix([a0, a1, resonant, leading]).jacobian(inputs)
    assert list(rows.row(2)) == [-Rational(4, 3)*I, -Rational(8, 3),
                                  -g/(3*I), 2*g/3, 5*I*g/3]
    assert list(rows.row(3)) == [-4, -4*I, g, I*g, -g]
    assert expand(rows[2:4, 2:4].det()) == -g*g
    assert (rows[:3, :].rank(), rows[:3, 2:].rank()) == (3, 1)
    assert (rows.rank(), rows[:, 2:].rank()) == (4, 2)
    print("Verified the resonant coefficient rows and their two rank differences")


def radial_leading_determinant():
    field = GF(1009)
    columns = list(zip((0, 1, 0, 1, 2, 3, 4, 5), (0, 0, 1, 1, 1, 1, 1, 1)))
    rows = []
    for k in range(8):
        row = []
        for j, parity in columns:
            if k < j or (k-j) % 2 != parity:
                row.append(field.zero())
                continue
            height = 17*k-j+12-parity
            integral = sum(
                field((-1)**(j-r)*comb(j, r))/field(height+2*r+1)
                for r in range(j+1)
            )
            row.append(comb(k, j)*integral)
        rows.append(row)
    assert matrix(field, rows).det() == 789
    print("Verified the nonzero degree-12 leading coefficient on every radial line")


if __name__ == "__main__":
    finite_critical_orbits()
    finite_marker_rank()
    normalized_period_rows()
    radial_leading_determinant()
