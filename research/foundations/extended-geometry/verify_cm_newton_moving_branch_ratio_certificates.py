#!/usr/bin/env python3
"""Verify the conductor-18 unit ideal for a moving branch ratio."""

from collections import defaultdict

from sage.all import GF, QQ, ZZ, PolynomialRing, gcd, matrix, vector


def residue_rows(odd_degree, last_shift):
    odd_names = tuple(f"l{j}" for j in range(1, odd_degree + 1))
    phase_ring = PolynomialRing(QQ, names=("u", *odd_names, "v", "eta"))
    u, *odd_coefficients, v, eta = phase_ring.gens()
    z = PolynomialRing(phase_ring, "z").gen()
    quadratic = z**2 - eta*z - 1
    phase = {(i - 1, 0): value for i, value in (quadratic**4).dict().items()}
    phase.update({(-1, 2): -1, (0, 1): u, (0, 2): v})
    phase.update({(j, 1): u * value for j, value in enumerate(odd_coefficients, start=1)})
    coherent_weights = {-i - 1: value for i, value in (quadratic**2).dict().items()}
    ring = PolynomialRing(QQ, names=("r", *odd_names, "v", "eta"))
    coefficients = {(0, 0): phase_ring.one()}

    def normalize(polynomial, parity):
        assert all(index[0] % 2 == parity for index in polynomial.dict())
        return ring({
            (degree // 2, *rest): value
            for (degree, *rest), value in polynomial.dict().items()
        })

    for k in range(1, last_shift + 1):
        product = defaultdict(phase_ring.zero)
        for (angle, deficit), value in coefficients.items():
            for (step, loss), weight in phase.items():
                # A remaining factor lowers the angle by at most one.
                if angle + step <= last_shift - k:
                    product[angle + step, deficit + loss] += weight * value
        coefficients = {key: value for key, value in product.items() if value}
        coherent = phase_ring.zero()
        poles = [phase_ring.zero(), phase_ring.zero()]
        for (angle, deficit), value in coefficients.items():
            if deficit % 2 == 0 and angle in coherent_weights:
                coherent += coherent_weights[angle] * value / QQ(13 + 18*k - deficit)
            if deficit % 2 and angle in (0, -1):
                poles[-angle] += value / QQ(12 + 18*k - deficit)
        yield k, normalize(coherent, 0), [normalize(entry, 1) for entry in poles]


def residual_rows(odd_degree, last_shift):
    solved_poles = []
    substitutions = {}
    stored = []
    for k, coherent, poles in residue_rows(odd_degree, last_shift):
        ring = coherent.parent()
        r, l, *_, v, eta = ring.gens()
        row = coherent + sum(entry * value for entry, value in zip(poles, solved_poles))
        if k <= 2:
            pivot = poles[k - 1]
            assert not any(poles[k:]) and pivot.is_constant() and pivot
            solved_poles.append(-row / QQ(pivot))
            continue
        row = ring(row.subs(substitutions))
        if k == 3 and odd_degree == 1:
            pivot = row.coefficient({r: 1})
            assert pivot and pivot.is_constant() and row.degree(r) == 1
            radial_square = (
                (-6*l**2 + 48*l*v) / 427
                + (10092*eta*l + 3936*eta*v) / 20923
                + (1309476 - 865554*eta**2) / 1401841
            )
            assert ring(row.subs({r: radial_square})) == 0
            substitutions[r] = radial_square
            continue
        if 3 <= k <= odd_degree + 1:
            variable = ring.gen(k - 1)
            pivot = row.coefficient({variable: 1})
            assert pivot and pivot.is_constant() and row.degree(variable) == 1
            substitutions[variable] = -(row - pivot * variable) / QQ(pivot)
            assert variable not in substitutions[variable].variables()
            continue
        assert not set(substitutions) & set(row.variables())
        kept = tuple(i for i, variable in enumerate(ring.gens()) if variable not in substitutions)
        assert all(sum(exponents[1:]) % 2 == (k - 1) % 2 for exponents in row.dict())
        integral = row * row.denominator()
        content = gcd([ZZ(value) for value in integral.coefficients()])
        stored.append({
            tuple(exponents[i] for i in kept): ZZ(value / content)
            for exponents, value in integral.dict().items()
        })
    assert solved_poles == [QQ(2) / 31, (-2*l + 16*v) / 31 + 3364*eta / 1519]
    return stored


def monomial_multipliers(weights, bound, parity):
    if len(weights) == 1:
        for exponent in range(parity, bound + 1, 2):
            yield (exponent,)
        return
    weight, *remaining = weights
    for exponent in range(bound // weight + 1):
        tail_bound = bound - weight * exponent
        tail_parity = (parity - weight * exponent) % 2
        for tail in monomial_multipliers(remaining, tail_bound, tail_parity):
            yield (exponent, *tail)


def verify(odd_degree, last_shift, degree_limit, matrix_size):
    monomials = {}
    multiples = []
    polynomials = residual_rows(odd_degree, last_shift)
    weights = (1, 1, 1) if odd_degree == 1 else (2, 1, 1, 1)
    degrees = [
        max(sum(a*b for a, b in zip(exponents, weights)) for exponents in polynomial)
        for polynomial in polynomials
    ]
    assert degrees == list(range(odd_degree + 1 if odd_degree > 1 else 3, last_shift))
    # Take every multiple within the weighted degree bound and with even parity.
    for polynomial, degree in zip(polynomials, degrees):
        for shift in monomial_multipliers(weights, degree_limit - degree, degree % 2):
            row = {}
            for exponents, value in polynomial.items():
                key = tuple(a + b for a, b in zip(shift, exponents))
                row[monomials.setdefault(key, len(monomials))] = value
            multiples.append(row)
    coefficients = {
        (i, j): value for i, row in enumerate(multiples) for j, value in row.items()
    }
    reduced = matrix(GF(65521), len(multiples), len(monomials), coefficients, sparse=False)
    pivot_rows = reduced.transpose().pivots()
    pivot_columns = reduced.pivots()
    assert (reduced.nrows(), reduced.ncols(), len(pivot_rows)) == matrix_size
    selected = matrix(ZZ, len(pivot_rows), len(monomials), {
        (i, j): value
        for i, row_index in enumerate(pivot_rows)
        for j, value in multiples[row_index].items()
    })
    constant = monomials[(0,) * len(weights)]
    square = selected.matrix_from_columns(pivot_columns).transpose().dense_matrix()
    target = vector(ZZ, [int(j == constant) for j in pivot_columns])
    # The modular calculation chooses a submatrix. The identity is checked over Q.
    multipliers = square.solve_right(target)
    expected = vector(ZZ, [int(j == constant) for j in range(len(monomials))])
    assert multipliers * selected == expected
    print(f"Moving branch ratio, odd degree {odd_degree}: exact rational identity")


if __name__ == "__main__":
    verify(1, 7, 12, (340, 251, 245))
    verify(3, 10, 16, (1912, 1484, 1427))
    verify(4, 11, 18, (2548, 2199, 2078))
