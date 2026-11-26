#!/usr/bin/env python3
"""Verify conductor-18 reconstruction certificates over Q with SageMath."""

from collections import defaultdict
import json
from math import comb
from pathlib import Path

from sage.all import GF, LaurentPolynomialRing, PolynomialRing, QQ, ZZ, matrix, vector


def phase_terms(odd_terms, even_radial=0, even_linear=0, *, top_sign=-1):
    # Keys are the angular exponent and the weight deficit.
    terms = {(2 * j - 1, 0): top_sign**j * comb(4, j) for j in range(5)}
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
    verify_quintic_reduction_size(rows)
    verify_reduced_degree_twelve_algebra(rows)
    verify_one_normalized_pole(rows)


def coherent_quintic(b, b1, b2, b3, c, e, *, top_sign=1):
    return (155708*c*b**5 + 1407*b1*b**4 + top_sign*8880*c*b**3
            + (3417*b1*e + top_sign*13362*b1 + 1608*b2*c + 201*b3)*b**2
            - (1608*b1*b1*c + 402*b1*b2)*b + 201*b1**3)


def verify_coherent_quintic():
    ring = PolynomialRing(QQ, names=("b", "b1", "b2", "b3", "b4", "b5", "c", "e"))
    b, b1, b2, b3, b4, b5, c, e = ring.gens()
    phase = phase_terms((b, b1, b2, b3, b4, b5), c, e, top_sign=1)
    columns = [(1, 0), (3, 0), (5, 0)] + [(j, 1) for j in range(3)]
    rows = list(reconstruction_rows(phase, columns, 4, ring))
    assert not any(rows[0])
    reduced = matrix(ring, [[a1 + 2*a3 + a5, *poles] for a1, a3, a5, *poles in rows[1:]])
    scale = -QQ(4096)/433798734017054481075
    assert reduced.det() == scale*coherent_quintic(b, b1, b2, b3, c, e)
    print("The two-coherent, three-correction response matrix has the stated quintic minor over Q")


def verify_coherent_quintic_norm():
    # At b1=b2=c=e=1, b3=0, use the sign rotation to evaluate S7.
    field = GF(1009)
    i = field(469)
    assert i*i == -1
    ring = PolynomialRing(field, "b")
    b = ring.gen()
    quintic = coherent_quintic(b, 1, 1, 0, 1, 1)
    assert quintic.degree() == 5 and quintic[0]
    residual, = reduced_residual_rows(5, 7, even_linear=True, primitive=False)
    cleared = ring.zero()
    for (radial, linear, quadratic, even_radial, even_linear), coefficient in residual.dict().items():
        power = 6 + 2*radial - linear - quadratic
        assert power >= 0
        cleared += field(coefficient)*(-1)**(radial + quadratic + even_linear)*i**(linear + even_radial)*b**power
    algebra = ring.quotient(quintic.monic(), "u")
    assert algebra(b).matrix().det()**6 == 336
    assert algebra(cleared).matrix().det() == 801
    print("Row 7 is a unit in the rank-five algebra on a nonempty open set of the base")


def verify_quintic_reduction_size(rows):
    base = LaurentPolynomialRing(QQ, names=("b1", "b2", "b3", "c", "e"))
    b1, b2, b3, c, e = base.gens()
    ring = PolynomialRing(base, "b")
    b = ring.gen()
    quintic = coherent_quintic(b, b1, b2, b3, c, e, top_sign=-1)/base(155708*c)
    for row, expected in [(rows[0], (71, 6, 12, 353, 8)), (rows[-1], (1330, 13, 27, 12816, 22))]:
        terms = row.dict()
        clearing_power = max(linear + quadratic for _, linear, quadratic, _, _ in terms)
        cleared = sum(
            ring(coefficient*b1**linear*b2**quadratic*c**even_radial*e**even_linear)
            * b**(clearing_power + 2*radial - linear - quadratic)
            for (radial, linear, quadratic, even_radial, even_linear), coefficient in terms.items()
        )
        remainder = cleared % quintic
        assert remainder.degree() < 5
        exponents = [exponent for coefficient in remainder for exponent in coefficient.exponents()]
        assert (len(terms), clearing_power, cleared.degree(), len(exponents),
                max(-exponent[3] for exponent in exponents)) == expected
    print("Quintic normal forms: row 7 grows from 71 to 353 terms, row 16 from 1330 to 12816")


def verify_reduced_degree_twelve_algebra(rows):
    base = PolynomialRing(QQ, names=("b1", "b2", "c", "e"))
    b1, b2, c, e = base.gens()
    ring = PolynomialRing(base, "b")
    b = ring.gen()
    cleared = []
    for row in rows:
        terms = row.dict()
        power = max(linear + quadratic - 2*radial for radial, linear, quadratic, _, _ in terms)
        cleared.append(sum(
            ring(coefficient*b1**linear*b2**quadratic*c**even_radial*e**even_linear)
            * b**(power + 2*radial - linear - quadratic)
            for (radial, linear, quadratic, even_radial, even_linear), coefficient in terms.items()
        ))
    first = cleared[0]
    leading = -QQ(23055471104)/584486340375375
    assert first.degree() == 12 and first.leading_coefficient() == leading
    assert first[0] == -QQ(2048)/64942926708375*b1**6
    monic = first/leading
    for index, count in [(1, 105), (9, 2753)]:
        remainder = cleared[index] % monic
        assert remainder.degree() < 12
        assert sum(len(coefficient.dict()) for coefficient in remainder) == count
    print("Degree-12 normal forms on the four-parameter base: 105 and 2753 terms")


def verify_331_radial_leading_minor():
    field = GF(331)
    columns = list(zip((0, 1, 0, 1, 2, 3, 4, 5), (0, 0, 1, 1, 1, 1, 1, 1)))
    leading = matrix(field, [
        [sum(field(comb(k, j)*(-1)**(j-a)*comb(j, a))/field(17*k-j+13-parity+2*a)
             for a in range(j+1)) if k >= j and (k-j) % 2 == parity else 0
         for j, parity in columns]
        for k in range(8)
    ])
    assert leading.det() == 5
    print("The common degree-12 radial leading coefficient is 5 modulo 331")


def verify_one_normalized_pole(rows):
    """Exclude a single normalized pole at 331 using S_7 through S_11."""
    finite_field = GF(331)
    residue_ring = PolynomialRing(finite_field, names=("r", "l", "h", "v", "w"))
    r, l, h, v, w = residue_ring.gens()
    source = rows[0].parent()

    def leading(row, index):
        variable = source.gen(index)
        coefficient = row.coefficient({variable: row.degree(variable)})
        return residue_ring({exponent: finite_field(value)
                             for exponent, value in coefficient.dict().items()})

    for index, expected in enumerate((17, -121, -121)):
        assert leading(rows[0], index) == expected

    c_leads = [leading(row, 3) for row in rows[:2]]
    assert -102 * c_leads[0] - 41 * c_leads[1] == r
    assert c_leads[0].subs({r: 0}) == 109 * (l**2 - h - 34)
    e_leads = [leading(row, 4) for row in rows[:3]]
    assert 56 * e_leads[0] + 96 * v * e_leads[1] - 120 * e_leads[2] == r
    assert e_leads[1] == 81 * l
    assert e_leads[0].subs({r: 0, l: 0}) == 117 * (h - 137)

    # Put s = 1/w. Homogenize in w to the rational degree before
    # reducing modulo 331 or truncating at order two in s.
    base = PolynomialRing(QQ, "v")
    c = base.gen()
    jets = PolynomialRing(base, "s")
    s = jets.gen()
    coordinates = [
        QQ(342144) / 19747853293 * s**2,
        QQ(576) / 289 * c * s + QQ(257472) / 83521 * c * s**2,
        QQ(1466) / 1139 + QQ(54144) / 4913 * s
        + (QQ(216576) / 83521 * c**2
           + QQ(10771301327656320) / 235622922126379) * s**2,
        c,
    ]

    initial_value = source.hom([base.zero(), base.zero(), base(QQ(1466) / 1139),
                                c, base.zero()], base)

    e = source.gen(4)
    char_zero_leads = [row.coefficient({e: row.degree(e)}) for row in rows[:3]]
    jacobian = matrix(base, [[initial_value(f.derivative(x)) for x in source.gens()[:3]]
                             for f in char_zero_leads])
    determinant = jacobian.det()
    assert determinant.degree() == 0 and finite_field(determinant[0])
    # Reduction here also tests integrality of the matrix entries and jets.
    for entry in jacobian.list():
        for coefficient in entry.coefficients():
            finite_field(coefficient)
    for coordinate in coordinates:
        for block in jets(coordinate).coefficients():
            for coefficient in block.coefficients():
                finite_field(coefficient)

    def truncate(poly):
        return jets({j: entry for j, entry in poly.dict().items() if j <= 2})

    def evaluate(row):
        degree = row.degree(e)
        result = jets.zero()
        for (i, j, k, v_degree, w_degree), coefficient in row.dict().items():
            order = degree - w_degree
            if order <= 2:
                term = truncate(coordinates[0]**i * coordinates[1]**j)
                term = truncate(term * coordinates[2]**k)
                result += coefficient * c**v_degree * truncate(s**order * term)
        return result

    assert all(evaluate(row) == 0 for row in rows[:3])
    obstruction = -QQ(56321516765184) / 1239568379329485639766985107
    assert finite_field(obstruction) and evaluate(rows[4]) == obstruction * s**2
    # At a c-pole, keep l,w integral and solve for r,h using S_7,S_8.
    base = PolynomialRing(QQ, names=("l", "w"))
    ell, eps = base.gens()
    jets = PolynomialRing(base, "s")
    s = jets.gen()
    coordinates = [jets.zero(), jets(ell), jets(ell**2 + QQ(370)/67),
                   jets.zero(), jets(eps)]
    initial_value = source.hom([coordinate[0] for coordinate in coordinates], base)
    c = source.gen(3)
    char_zero_leads = [row.coefficient({c: row.degree(c)}) for row in rows[:2]]
    jacobian = matrix(QQ, [[QQ(initial_value(f.derivative(source.gen(index))))
                            for index in (0, 2)] for f in char_zero_leads])
    assert finite_field(jacobian.det())

    def check_integral(poly):
        for block in poly.coefficients():
            for coefficient in block.coefficients():
                finite_field(coefficient)

    def c_chart_row(row, bound):
        degree = row.degree(c)

        def truncation(poly, order):
            return jets({j: entry for j, entry in poly.dict().items() if j <= order})

        result = jets.zero()
        for (i, j, k, v_degree, w_degree), coefficient in row.dict().items():
            order = degree - v_degree
            if order <= bound:
                term = truncation(coordinates[0]**i, bound - order)
                term = truncation(term * coordinates[2]**k, bound - order)
                result += coefficient * ell**j * eps**w_degree * s**order * term
        return result

    inverse = jacobian.inverse()
    assert not any(c_chart_row(row, 0) for row in rows[:5])
    for order in range(1, 7):
        residual = vector(base, [c_chart_row(row, order)[order] for row in rows[:2]])
        corrections = -inverse * residual
        for index, correction in zip((0, 2), corrections, strict=True):
            coordinates[index] += correction * s**order
    for coordinate in coordinates:
        check_integral(coordinate)

    substituted = [c_chart_row(row, 6) for row in rows[:5]]
    assert not any(substituted[:2])
    assert all(not f[j] for f in substituted[2:] for j in range(4))
    pivot = substituted[2]
    leading_pivot = pivot[4]
    assert leading_pivot.degree(ell) == 6
    assert finite_field(leading_pivot.lc())
    reduced_ring = PolynomialRing(finite_field, names=base.variable_names())
    for row, slope, constant in zip(substituted[3:], (92, 244), (121, 243), strict=True):
        for j in (4, 5, 6):
            quotient, remainder = row[j].quo_rem(leading_pivot)
            check_integral(jets(quotient))
            if j < 6:
                assert remainder == 0
                row -= s**(j - 4) * quotient * pivot
            else:
                reduced = reduced_ring({e: finite_field(x)
                                        for e, x in remainder.dict().items()})
                assert reduced == slope * reduced_ring.gen(1) + constant
    assert finite_field(92 * 243 - 244 * 121)
    print("Excluded each of the five possible single normalized poles.")

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
    verify_coherent_quintic()
    verify_coherent_quintic_norm()
    verify_331_radial_leading_minor()
    verify_three_pole_coordinate_slices()
    verify_three_pole_prime_331()

