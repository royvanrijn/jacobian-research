#!/usr/bin/env python3
"""Verify conductor-18 reconstruction certificates over Q with SageMath."""

from collections import defaultdict
import json
from math import comb
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


def phase_terms(odd_terms, even_radial=0, even_linear=0):
    # Keys are the angular exponent and the weight deficit.
    terms = {(2 * j - 1, 0): (-1)**j * comb(4, j) for j in range(5)}
    terms[-1, 2] = -1
    terms.update({(j, 1): value for j, value in enumerate(odd_terms)})
    terms[0, 2] = even_radial
    terms[1, 2] = even_linear
    return {key: value for key, value in terms.items() if value}


def reconstruction_rows(phase, columns, final_shift, ring):
    coefficients = {(0, 0): ring.one()}
    for k in range(final_shift + 1):
        row = []
        for j, parity in columns:
            entry = sum(
                value / QQ(13 - parity + 18 * k - deficit)
                for (angle, deficit), value in coefficients.items()
                if angle == -j and deficit % 2 == parity
            )
            row.append(ring(entry))
        yield row
        if k == final_shift:
            break
        product = defaultdict(ring.zero)
        for (angle, deficit), value in coefficients.items():
            for (step, weight), coefficient in phase.items():
                # Remaining factors can each lower the angle by only one.
                if angle + step <= final_shift - k - 1:
                    product[angle + step, deficit + weight] += value * coefficient
        coefficients = {key: value for key, value in product.items() if value}


def verify_exceptional_radial_line():
    ring = PolynomialRing(QQ, "b")
    b = ring.gen()
    phase = phase_terms((b, 0, 0, 1))
    columns = [(j, 0) for j in range(7)] + [(j, 1) for j in range(6)]
    rows = list(reconstruction_rows(phase, columns, 13, ring))
    delta = matrix(ring, rows[:13]).det()
    replacement = matrix(ring, rows[:12] + rows[13:]).det()
    assert delta.degree() == 42 and delta(0) != 0
    assert delta.gcd(delta.derivative()) == 1
    assert replacement.degree() == 39
    gcd, left, right = delta.xgcd(replacement)
    assert gcd == 1 and left * delta + right * replacement == 1
    print("Exceptional radial line: 42 distinct roots, two minors with exact Bezout identity")


def normalize_radial_power(polynomial, odd, ring):
    terms = {}
    for (radial_degree, *others), value in polynomial.dict().items():
        assert radial_degree % 2 == odd
        terms[(radial_degree // 2, *others)] = value
    return ring(terms)


def solve_pivot(polynomial, variable):
    pivot = polynomial.coefficient({variable: 1})
    assert polynomial.degree(variable) == 1 and pivot and pivot.is_constant()
    solution = -(polynomial - pivot * variable) / QQ(pivot)
    assert variable not in solution.variables()
    return solution


def primitive_polynomial(polynomial):
    integral = polynomial * polynomial.denominator()
    return integral / QQ(integral.content())


def reduced_residual_rows(odd_degree, final_shift=None, *, even_linear=False, primitive=True):
    even_names = ("v", "w") if even_linear else ("v",)
    odd_names = tuple(f"l{j}" for j in range(1, odd_degree + 1))
    names = odd_names + even_names
    phase_ring = PolynomialRing(QQ, names=("u", *names))
    u, *others = phase_ring.gens()
    phase = phase_terms((u, *(u * value for value in others[:odd_degree])), *others[odd_degree:])

    ring = PolynomialRing(QQ, names=("r", *names), order="degrevlex")
    nuisance_count = 3 if even_linear else 2
    surviving_odd = ("l", "h") if even_linear else ("l",)
    subring = PolynomialRing(QQ, names=("r", *surviving_odd, *even_names), order="degrevlex")
    generators = subring.gens()
    first_even = 1 + len(surviving_odd)
    reduction = ring.hom(
        [*generators[:first_even], *([0] * (odd_degree - len(surviving_odd))),
         *generators[first_even:]],
        subring,
    )
    poles = []
    odd_coefficients = {}
    residuals = []
    columns = [(1, 0), (3, 0), (5, 0)] + [(j, 1) for j in range(nuisance_count)]
    rows = reconstruction_rows(phase, columns, final_shift or odd_degree + 5, phase_ring)
    for k, (a1, a3, a5, *pole_columns) in enumerate(rows):
        if k == 0:
            assert not any((a1, a3, a5, *pole_columns))
            continue
        coherent = normalize_radial_power(a1 - 2*a3 + a5, 0, ring)
        nuisance = [normalize_radial_power(column, 1, ring) for column in pole_columns]
        residual = coherent + sum(column * entry for column, entry in zip(nuisance, poles))
        if k <= nuisance_count:
            pivot = nuisance[k-1]
            assert pivot and pivot.is_constant() and not any(nuisance[k:])
            poles.append(-residual / QQ(pivot))
            continue

        residual = ring(residual.subs(odd_coefficients))
        if k <= odd_degree + 1:
            variable = ring.gen(k-1)
            odd_coefficients[variable] = solve_pivot(residual, variable)
            continue

        assert not set(residual.variables()) & odd_coefficients.keys()
        reduced = reduction(residual)
        assert reduced
        residuals.append(primitive_polynomial(reduced) if primitive else reduced)
    if even_linear:
        r, l, h, *_, v, w = ring.gens()
        assert poles == [
            QQ(2)/31,
            (16*v-2*l)/31,
            (2*l*l-16*l*v-2*h+34*w)/31 + QQ(427)*r/93 - QQ(8908)/2077,
        ]
    return residuals


def verify_odd_cubic_family():
    residuals = reduced_residual_rows(3)
    assert [(f.degree(), len(f.monomials())) for f in residuals] == [
        (4, 13), (5, 16), (6, 27), (7, 33),
    ]
    ring = residuals[0].parent()
    multipliers = ring.one().lift(ring.ideal(residuals))
    assert sum(g * f for g, f in zip(multipliers, residuals, strict=True)) == 1
    print("Odd cubic polynomial and even radial term: rows 5 through 8 certify the unit ideal")


def verify_odd_quartic_family():
    residuals = reduced_residual_rows(4)
    assert [(f.degree(), len(f.monomials())) for f in residuals] == [
        (5, 19), (6, 26), (7, 37), (8, 48),
    ]
    ring = residuals[0].parent()
    certificate = Path(__file__).with_name("cm_newton_conductor_18_quartic_bezout.json")
    multipliers = []
    for polynomial in json.loads(certificate.read_text(encoding="utf-8")):
        terms = {(i, j, k): QQ(numerator) / QQ(denominator)
                 for i, j, k, numerator, denominator in polynomial}
        assert len(terms) == len(polynomial)
        multipliers.append(ring(terms))
    assert sum(g * f for g, f in zip(multipliers, residuals, strict=True)) == 1
    print("Odd quartic polynomial and even radial term: rows 6 through 9 certify the unit ideal")


def weighted_exponents(weights, bound):
    if not weights:
        yield ()
        return
    first, *rest = weights
    for exponent in range(bound // first + 1):
        for tail in weighted_exponents(rest, bound - exponent * first):
            yield (exponent, *tail)


def coefficient_system(polynomials, weights, degree, parity_indices):
    order = (*parity_indices, *(i for i in range(len(weights)) if i not in parity_indices))
    ordered_weights = tuple(weights[i] for i in order)
    positions = tuple(order.index(i) for i in range(len(weights)))
    columns = {}
    coefficients = {}
    multiples = []
    for row in polynomials:
        terms = row.dict()
        weight = max(sum(d * w for d, w in zip(exponent, weights)) for exponent in terms)
        parities = {sum(exponent[i] for i in parity_indices) % 2 for exponent in terms}
        assert len(parities) == 1
        parity = next(iter(parities))
        for ordered in weighted_exponents(ordered_weights, degree - weight):
            if sum(ordered[:len(parity_indices)]) % 2 != parity:
                continue
            shift = tuple(ordered[i] for i in positions)
            index = len(multiples)
            multiples.append((row, shift))
            for exponent, coefficient in terms.items():
                lifted = tuple(a + b for a, b in zip(exponent, shift))
                column = columns.setdefault(lifted, len(columns))
                coefficients[index, column] = coefficient
    return coefficients, len(multiples), len(columns), columns[(0,) * len(weights)]


def verify_unit_combination(polynomials, weights, degree, parity_indices):
    # Pick an invertible minor modulo a prime. Solve and verify over Q.
    data, row_count, column_count, target = coefficient_system(
        polynomials, weights, degree, parity_indices,
    )
    coefficients = {position: ZZ(value) for position, value in data.items()}
    reduced = matrix(GF(65521), row_count, column_count, coefficients, sparse=False)
    rows, pivots = reduced.transpose().pivots(), reduced.pivots()
    assert target in pivots
    row_indices = {row: i for i, row in enumerate(rows)}
    exact = matrix(ZZ, len(rows), column_count, {
        (row_indices[i], j): value for (i, j), value in coefficients.items() if i in row_indices
    }, sparse=True)
    rhs = vector(ZZ, [int(j == target) for j in pivots])
    square = exact.matrix_from_columns(pivots).transpose().dense_matrix()
    multipliers = square._solve_right_nonsingular_square(rhs, check_rank=False, algorithm="flint")
    target_row = vector(QQ, [int(j == target) for j in range(column_count)])
    assert multipliers * exact == target_row
    return len(rows)


def verify_three_pole_coordinate_slices():
    rows = reduced_residual_rows(5, 16, even_linear=True)
    source = rows[0].parent()
    weights = (2, 1, 2, 1, 2)
    for index, label, bound, rank in [
        (4, "e", 18, 1085), (2, "b2", 18, 1071),
        (3, "c", 20, 957), (1, "b1", 20, 953),
    ]:
        names = list(source.variable_names())
        names.pop(index)
        ring = PolynomialRing(QQ, names=tuple(names), order="degrevlex")
        images = list(ring.gens())
        images.insert(index, 0)
        specialize = source.hom(images, ring)
        polynomials = [primitive_polynomial(specialize(row)) for row in rows]
        slice_weights = weights[:index] + weights[index + 1:]
        parity = tuple(i for i, w in enumerate(slice_weights) if w == 1)
        assert verify_unit_combination(polynomials, slice_weights, bound, parity) == rank
        print(f"Three poles, {label}=0: shifts 7 through 16 certify the unit ideal")


def verify_three_pole_prime_331():
    prime = 331
    assert prime % 36 == 7 and 13 + 18 * 16 < prime
    field = GF(prime)
    rows = reduced_residual_rows(5, 16, even_linear=True, primitive=False)
    ring = PolynomialRing(field, names=rows[0].parent().variable_names())
    reductions = [ring({e: field(c) for e, c in row.dict().items()}) for row in rows]
    data, row_count, column_count, target = coefficient_system(reductions, (2, 1, 2, 1, 2), 22, (1, 3))
    echelon = matrix(field, row_count, column_count, data, sparse=False).echelon_form()
    pivots = echelon.pivots()
    assert echelon[pivots.index(target)] == vector(field, [int(j == target) for j in range(column_count)])
    print("Three poles, unrestricted coefficients: shifts 7 through 16 certify 1 over F_331")


def verify_third_pole_example():
    prime, order = 151, 109
    ring = PolynomialRing(GF(prime), names=("x", "t"))
    phase = phase_terms((2, 3, 5, 7, 11, 13), 17, 19)
    # Form x*f in the literal x,t chart to avoid negative exponents.
    lifted = ring({(i + 1, 18 - deficit - 2*i): value for (i, deficit), value in phase.items()})
    coefficient = (lifted**order)[order + 2*prime, 8*prime - 1]
    assert coefficient == 97
    print("Third pole is attained: C2=97 modulo 151")


if __name__ == "__main__":
    verify_exceptional_radial_line()
    verify_odd_cubic_family()
    verify_odd_quartic_family()
    verify_third_pole_example()
    verify_three_pole_coordinate_slices()
    verify_three_pole_prime_331()

