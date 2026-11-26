"""Exact symbolic and native coefficient checks for the prime-first lift."""

from math import comb, factorial

import sympy as sp

from verify_binary_even_quadratic_prime_front_zero_lifts import (
    factorial_polynomial,
    moment_jet,
    value_and_gradient,
)


X, G, H, z, W = sp.symbols("X G H z W")
D0 = 1 + 4 * W + 24 * W**2
Q1 = 120 * W**2 + 12 * W + 1
Q2 = 360 * W**2 + 24 * W + 1
I = 129600 * W**5 + 43200 * W**4 + 14400 * W**3 + 1800 * W**2 + 225 * W + 128
LOWER_H = (12 * G**3 - 3 * G - 6 * X * G + 2) / (6 * X)
LEADING_G = sp.Rational(8, 15) / D0


def derivative(poly):
    return sp.expand(
        X * sp.diff(poly, X)
        + (G - H) * sp.diff(poly, G) / 2
        + ((G - 1) / X - H) * sp.diff(poly, H) / 2
    )


def g_small(order, variable):
    return sum(comb(order, i) * factorial(2 * i) * variable**i for i in range(order + 1))


def exterior(order, layer):
    result = 0
    for u in range(order + 1):
        r = 2 * (order - u)
        polynomial = sum(
            comb(u, i) * sp.rf(r + layer + 1, 2 * i) * X**i for i in range(u + 1)
        )
        result += 2 * comb(order, u) * (-2)**u * factorial(r) * z ** (r + layer) * polynomial
    return result


def native_jets():
    g2 = 1 + 6 * X * G - 3 * X * H
    g3 = 1 + 20 * X * g2 - 15 * X * G
    g4 = 1 + 42 * X * g3 - 35 * X * g2
    jets = []
    for j, gj, next_g in zip(range(3), (G, g2, g3), (g2, g3, g4)):
        eta = -(-4)**(j + 1) * factorial(j) / sp.factorial2(2 * j + 1)
        nu = -sp.Rational((-2)**(j + 1), j + 1)
        small = g_small(j + 1, X)
        leading = exterior(j, 1) + eta * gj
        first = exterior(j, 2) + eta * z * (2 * derivative(gj) + gj) + nu * small
        shifted = 2 * derivative(gj) + 2 * gj
        second = (
            exterior(j, 3)
            + eta * z**2 * (2 * derivative(shifted) + shifted) / 2
            + nu * z * (2 * derivative(small) + small)
            + eta * next_g / (2 * j + 3)
        )
        jets.append(tuple(sp.expand(value) for value in (leading, first, second)))
    return jets


def rational_mod(expression, prime):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    return int(numerator) * pow(int(denominator), -1, prime) % prime


def verify_native(jets):
    for prime in (11, 13, 17, 19, 23, 37, 41, 47):
        h = (prime - 1) // 2
        sign = (-1)**h
        for j, layers in enumerate(jets):
            exponent = 2 * prime + 2 * j + 1
            coefficients, raw_min = moment_jet(prime + j, prime, exponent + 2)
            assert exponent + 2 < raw_min + prime
            assert min(power for power, _, _ in coefficients) == exponent
            for residue_x in (1, 2, prime - 1):
                changes = {X: residue_x, G: factorial_polynomial(h + 1, residue_x, prime),
                           H: factorial_polynomial(h, residue_x, prime)}
                for residue_z in (1, 2, prime - 1):
                    changes[z] = residue_z
                    for offset, layer in enumerate(layers):
                        power = exponent + offset
                        direct, _, _ = value_and_gradient(
                            coefficients, power, residue_x, sign * residue_z, prime
                        )
                        assert rational_mod(layer.subs(changes), prime) == sign**power * direct % prime
    print("Verified 648 native moment layer evaluations.")


def reduced_leading(expression):
    with_lower_h = sp.factor(expression.subs(H, LOWER_H))
    return sp.factor(with_lower_h.subs(X, W + 2 * G**2).subs(G, LEADING_G))


def residue(expression, modulus):
    a, b = sp.fraction(sp.cancel(expression))
    return sp.rem(sp.rem(a, modulus, W) * sp.invert(sp.rem(b, modulus, W), modulus, W), modulus, W)


def zero_mod(expression, modulus):
    numerator = sp.fraction(sp.cancel(expression))[0]
    return sp.rem(numerator, modulus, W) == 0


def eliminate_z(jets):
    argument = z**2 - 2 * X
    p = -8 * (1 - 2 * argument + 6 * argument**2)
    boundary = (sp.Rational(128, 15) - z * p, sp.Rational(8, 3) - z * sp.Rational(128, 15) + argument * p, 0)
    z0 = -2 * G
    _, first, second = jets[0]
    z1 = -first.subs(z, z0) / 2
    z2 = -(second.subs(z, z0) + z1 * sp.diff(first, z).subs(z, z0)) / 2
    equations = []
    for e, f, second in (*jets[1:], boundary):
        value = e.subs(z, z0)
        first_layer = f.subs(z, z0) + z1 * sp.diff(e, z).subs(z, z0)
        second_layer = (
            sp.sympify(second).subs(z, z0) + z1 * sp.diff(f, z).subs(z, z0)
            + z2 * sp.diff(e, z).subs(z, z0) + z1**2 * sp.diff(e, z, 2).subs(z, z0) / 2
        )
        equations.append(tuple(map(sp.expand, (value, first_layer, second_layer))))
    return equations


def verify_first(equations):
    e1, e2, eb = [entry[0] for entry in equations]
    assert sp.expand(e1 - sp.Rational(8, 3) * (-12 * G**3 + 3 * G + 6 * X * G - 2 + 6 * X * H)) == 0
    assert reduced_leading(e2) == reduced_leading(eb) == 0
    assert sp.cancel(reduced_leading(X) - I / (225 * D0**2)) == 0
    j1, j2, jb = [reduced_leading(derivative(value)) for value in (e1, e2, eb)]
    r1, r2, _ = [reduced_leading(entry[1]) for entry in equations]
    assert sp.cancel(j2 - 120 * D0 * Q1**2 / I) == 0
    assert sp.cancel(j1 * r2 - j2 * r1 - 80 * D0 * Q1 * Q2 / I) == 0
    assert all(zero_mod(value, Q1) for value in (j1, j2, jb))
    assert zero_mod(reduced_leading(derivative(derivative(eb))) - 448 * (10 * W + 1) / 25, Q1)
    assert zero_mod(reduced_leading(derivative(derivative(e2))), Q1)
    assert zero_mod(r2 - sp.Rational(8, 3), Q1)
    assert Q1.subs(W, -sp.Rational(1, 10)) == 1
    return (j1, j2), (r1, r2)


def second_data(equations):
    theta = (-32 * LEADING_G * g_small(3, W) + sp.Rational(512, 35)) / 720
    data = [
        tuple(reduced_leading(value) for value in (
            derivative(derivative(e)), derivative(r), s
        ))
        for e, r, s in equations[:2]
    ]
    d2, dr, second = data[1]
    data[1] = d2, dr, second + 24 * theta
    return data


def verify_second(equations, jacobian, first_layer):
    data = second_data(equations)
    j1, j2 = [residue(j, Q2) for j in jacobian]
    u = sp.rem(-residue(first_layer[1], Q2) * sp.invert(j2, Q2, W), Q2, W)
    assert u == sp.Rational(9967, 980) - sp.Rational(150, 49) * W
    shifts = [
        residue(d2, Q2) * u**2 / 2 + residue(dr, Q2) * u + residue(s, Q2)
        for d2, dr, s in data
    ]
    compatibility = sp.rem(j2 * shifts[0] - j1 * shifts[1], Q2, W)
    linear = 81133029162 * W + 3369377549
    assert sp.expand(compatibility - sp.Rational(16, 10409878395) * linear) == 0
    assert sp.resultant(linear, Q2, W) == 4 * 9 * 7 * 83 * 33047233 * 5944184339
    support = {2, 3, 5, 7, 83, 33047233}
    for polynomial in (D0, I, Q1):
        assert set(sp.factorint(abs(sp.resultant(polynomial, Q2, W)))) <= support
    assert set(sp.factorint(10409878395)) <= support
    return data


def verify_exceptional(jacobian, first_layer, data):
    cases = {
        83: {78: "scale", 16: 22},
        33047233: {26667925: "shape", 21801350: 31229009},
        5944184339: {5504500288: 2163193724, 43405095: 0},
    }
    for prime, outcomes in cases.items():
        assert sp.isprime(prime)
        assert len(outcomes) == 2
        for root, expected in outcomes.items():
            assert int(Q2.subs(W, root)) % prime == 0
            if expected == "scale":
                assert int(D0.subs(W, root)) % prime == 0
                continue
            scale = rational_mod(LEADING_G.subs(W, root), prime)
            shape = (root + 2 * scale**2) % prime
            if expected == "shape":
                assert shape == 0
                continue
            assert shape
            j1, j2 = [rational_mod(value.subs(W, root), prime) for value in jacobian]
            u = -rational_mod(first_layer[1].subs(W, root), prime) * pow(j2, -1, prime) % prime
            corrections = [
                (rational_mod(d2.subs(W, root), prime) * u**2 * pow(2, -1, prime)
                 + rational_mod(dr.subs(W, root), prime) * u
                 + rational_mod(s.subs(W, root), prime)) % prime
                for d2, dr, s in data
            ]
            actual = (j2 * corrections[0] - j1 * corrections[1]) % prime
            assert actual == expected
            if actual == 0:
                assert (prime, shape, scale) == (5944184339, 1951318268, 4250342174)
    print("Verified the first and second correction identities and all six exceptional residues.")


if __name__ == "__main__":
    jets = native_jets()
    verify_native(jets)
    equations = eliminate_z(jets)
    jacobian, first_layer = verify_first(equations)
    data = verify_second(equations, jacobian, first_layer)
    verify_exceptional(jacobian, first_layer, data)
