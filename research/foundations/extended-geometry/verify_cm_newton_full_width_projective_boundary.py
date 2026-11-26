#!/usr/bin/env python3
"""Check the projective boundary and local exclusions for the 331 test."""

from math import prod

from sage.all import (
    GF, QQ, PolynomialRing, PowerSeriesRing, binomial, gcd, identity_matrix, lcm, matrix, vector,
)

from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows


WEIGHTS = (2, 1, 2, 1, 2)


def weighted_degree(exponent, weights=WEIGHTS):
    return sum(a * b for a, b in zip(exponent, weights, strict=True))


def leading_forms(rows, weights):
    ring = rows[0].parent()
    forms = []
    for row in rows:
        terms = row.dict()
        degree = max(weighted_degree(exponent, weights) for exponent in terms)
        forms.append(ring({
            exponent: coefficient for exponent, coefficient in terms.items()
            if weighted_degree(exponent, weights) == degree
        }))
    return forms


def verify_boundary(rows):
    for k, row in enumerate(rows, 7):
        assert max(weighted_degree(e) for e in row.dict()) == k - 1
    source = rows[0].parent()
    r, l, h, v, w = source.gens()
    field = GF(331)
    target = PolynomialRing(field, names=source.variable_names(), order="degrevlex")
    rf, lf, hf, vf, wf = target.gens()

    ordinary, weighted = (leading_forms(rows, weights) for weights in [(1,) * 5, WEIGHTS])
    ordinary_components = [
        (r, 0, h, v, w),
        (r, 8*v, h, v, w),
    ]
    weighted_components = [
        (0, 0, 0, v, w),
        (0, l, 0, v, (8*l*v-l*l)/17),
        (0, 8*v, 17*w, v, w),
    ]
    for forms, components in [(ordinary, ordinary_components), (weighted, weighted_components)]:
        for images in components:
            substitute = source.hom(images, source)
            assert all(not substitute(form) for form in forms)

    def reduce_coefficients(forms):
        reductions = [target({e: field(c) for e, c in form.dict().items()}) for form in forms]
        assert all(reductions)
        return reductions

    ordinary_basis = target.ideal(reduce_coefficients(ordinary)).groebner_basis()
    assert not ((lf * (lf-8*vf))**4).reduce(ordinary_basis)
    weighted_basis = target.ideal(reduce_coefficients(weighted)).groebner_basis()
    for polynomial, power in [
        (rf, 9),
        (hf*(hf-17*wf), 7),
        (hf*(lf-8*vf), 7),
        ((lf*lf+hf)*(lf-8*vf)-lf*(hf-17*wf), 7),
    ]:
        assert not (polynomial**power).reduce(weighted_basis)
    print("At 331: two ordinary boundary components and three weighted boundary components")


def verify_first_weighted_component(rows):
    # On w=1, this component is r=l=h=0. The variable v is free.
    base = PolynomialRing(QQ, "v")
    v = base.gen()
    fractions = base.fraction_field()
    series = PowerSeriesRing(fractions, "tau", default_prec=4)
    tau = series.gen()
    terms = []
    for k, row in enumerate(rows, 7):
        expansion = []
        for exponent, coefficient in row.dict().items():
            r, l, h, v_degree, _ = exponent
            deficit = k - 1 - weighted_degree(exponent)
            assert deficit >= 0 and deficit % 2 == 0
            power = deficit // 2
            if r + l + h + power <= 3:
                expansion.append(((r, l, h, power), coefficient * v**v_degree))
        terms.append(expansion)

    normals = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]
    jacobian = matrix(base, [
        [sum(c for e, c in row if e == normal) for normal in normals]
        for row in terms
    ])

    def evaluate(row, coordinates, depth):
        r, l, h = coordinates
        values = (coefficient * r**a * l**b * h**c * tau**d
                  for (a, b, c, d), coefficient in row if a + b + c + d <= depth)
        return sum(values, series.zero()).add_bigoh(depth + 1)

    finite_base = PolynomialRing(GF(331), "v")
    finite_determinants = []
    for shifts in [(7, 8, 9), (7, 8, 10)]:
        pivots = [k - 7 for k in shifts]
        square = jacobian.matrix_from_rows(pivots)
        determinant = finite_base(square.det())
        assert determinant
        finite_determinants.append(determinant)
        inverse = square.change_ring(fractions).inverse()
        coordinates = [series.zero()] * 3
        for depth in range(1, 4):
            rhs = vector(fractions, [-evaluate(terms[i], coordinates, depth)[depth] for i in pivots])
            coordinates = [x + delta * tau**depth for x, delta in zip(coordinates, inverse * rhs)]
            coefficients = [evaluate(row, coordinates, depth)[depth] for row in terms]
            assert all(not coefficients[i] for i in pivots)
            if depth < 3:
                assert not any(coefficients)
        common = finite_base.zero()
        for coefficient in coefficients:
            if coefficient:
                numerator = finite_base(coefficient.numerator())
                denominator = finite_base(coefficient.denominator())
                assert denominator and not determinant**denominator.degree() % denominator
                common = common.gcd(numerator)
        assert common == 1
    assert finite_determinants[0].gcd(finite_determinants[1]) == 1
    print("At 331: exact transverse jets exclude the weighted component r=l=h=0 wherever w != 0")


def truncate(polynomial, weights, cap):
    return polynomial.parent()({
        e: c for e, c in polynomial.dict().items() if weighted_degree(e, weights) <= cap
    })


def truncated_powers(value, weights, cap, exponent):
    powers = [value.parent().one()]
    for _ in range(exponent):
        powers.append(truncate(powers[-1] * value, weights, cap))
    return powers


def substitute_truncated(terms, values, weights, cap):
    ring = values[0].parent()
    powers = [truncated_powers(x, weights, cap, max(e[i] for e in terms))
              for i, x in enumerate(values)]
    result = ring.zero()
    for exponents, coefficient in terms.items():
        term = ring(coefficient)
        for table, degree in zip(powers, exponents, strict=True):
            term = truncate(term * table[degree], weights, cap)
            if not term:
                break
        result += term
    return result


def verify_second_weighted_component(rows):
    # On l=1, its equations are r=h=0 and w=(8*v-1)/17.
    base = PolynomialRing(QQ, "v")
    v = base.gen()
    chart = PolynomialRing(base, names=("r", "h", "d", "tau"))
    r, h, d, tau = chart.gens()
    free = PolynomialRing(base, names=("x", "y"))
    x, y = free.gens()
    cap = 6
    free_weights = (1, 2)
    values = (r, chart.one(), h, chart(v), (8*v-1+d)/17, tau)
    charts = []
    for k, row in enumerate(rows[:4], 7):
        terms = {}
        for exponent, coefficient in row.dict().items():
            deficit = k - 1 - weighted_degree(exponent)
            assert deficit >= 0 and deficit % 2 == 0
            terms[tuple(exponent) + (deficit // 2,)] = coefficient
        chart_row = substitute_truncated(terms, values, (1, 1, 1, 1), cap)
        assert not chart_row.constant_coefficient()
        charts.append(chart_row)

    jacobian = matrix(base, [[f[e] for e in [(1, 0, 0, 0), (0, 0, 1, 0)]] for f in charts[:2]])
    determinant = jacobian.det()
    adjugate = jacobian.adjugate().change_ring(free)
    values = [free.zero(), determinant**2*x, free.zero(), determinant**2*y]
    for degree in range(1, cap + 1):
        errors = adjugate * vector(free, [
            substitute_truncated(f.dict(), values, free_weights, degree) for f in charts[:2]
        ])
        for i, error in enumerate(errors):
            corrections = {}
            for exponent, coefficient in error.dict().items():
                assert weighted_degree(exponent, free_weights) == degree
                quotient, remainder = coefficient.quo_rem(determinant)
                assert not remainder
                corrections[exponent] = quotient
            values[2*i] -= free(corrections)
    assert all(not substitute_truncated(f.dict(), values, free_weights, cap) for f in charts[:2])

    ninth, tenth = [substitute_truncated(f.dict(), values, free_weights, cap) for f in charts[2:]]
    assert not truncate(ninth, free_weights, 4) and not truncate(tenth, free_weights, 4)
    field = base.fraction_field()
    ratio_ring = PolynomialRing(field, "U")
    U = ratio_ring.gen()

    def leading_at(polynomial, degree):
        return ratio_ring(sum(field(c)*U**b for (a,b),c in polynomial.dict().items()
                              if a+2*b == degree))

    Q = leading_at(ninth, 5)
    gamma, remainder = leading_at(tenth, 5).quo_rem(Q)
    assert not remainder and gamma.degree() == 0
    H = leading_at(tenth, 6) - gamma*leading_at(ninth, 6)
    assert Q.degree() == 2 and H.degree() == 3
    finite = PolynomialRing(GF(331), "v")
    vf = finite.gen()
    exceptional = prod(vf - GF(331)(a) for a in [QQ(2)/17, QQ(5)/56, QQ(1)/8])
    for value in map(field, [determinant, Q[0], Q[2], H[3], Q.resultant(H)]):
        for part in (value.numerator(), value.denominator()):
            reduction = finite(part)
            assert reduction and all(exceptional % factor == 0 for factor, _ in reduction.factor())
    gamma_denominator = finite(gamma[0].denominator())
    assert gamma_denominator and all(exceptional % factor == 0 for factor, _ in gamma_denominator.factor())
    print("At 331: on the second weighted component with l != 0 only three v/l values survive")


def intersection_chart(rows, keep_l):
    """The weighted chart v=1, with q=l^2-8l+17w and tau=t^2."""
    L = PolynomialRing(QQ, "l")
    l = L.gen()
    x = L.fraction_field()
    w0 = (8*l-l*l)/17
    chart = []
    for row in rows:
        degree = max(weighted_degree(e, (2, 1, 2, 1, 2)) for e in row.dict())
        chart_row = {}
        for (r, j, h, v, w), coefficient in row.dict().items():
            tau, parity = divmod(degree-(2*r+j+2*h+v+2*w), 2)
            assert parity == 0
            for q in range(w+1):
                term = coefficient*l**j*binomial(w,q)*w0**(w-q)/17**q
                if keep_l:
                    exponent = (r, h, q, tau)
                    chart_row[exponent] = chart_row.get(exponent, x.zero()) + x(term)
                else:
                    for j_out, value in term.dict().items():
                        exponent = (r, h, j_out, q, tau)
                        chart_row[exponent] = chart_row.get(exponent, QQ.zero()) + value
        chart.append({e:c for e,c in chart_row.items() if c})
    return chart, (x if keep_l else QQ)


def implicit_intersection_jets(chart, field, weights, cap):
    """Use S7,S8 to solve r,h; retain exact weighted jets to the given cap."""
    basis = [(1, 0, *(0 for _ in weights)), (0, 1, *(0 for _ in weights))]
    square = matrix(field, [[row.get(e,0) for e in basis] for row in chart[:2]])
    determinant = square.det()
    assert determinant
    inverse = square.inverse()
    zero = (0,)*len(weights)

    def trunc_product(left, right, bound):
        result = {}
        for ex, a in left.items():
            for ey, b in right.items():
                ez = tuple(x + y for x, y in zip(ex, ey, strict=True))
                if weighted_degree(ez, weights) <= bound:
                    result[ez] = result.get(ez, 0) + a*b
        return {e:c for e,c in result.items() if c}

    def evaluate(polynomials, base, bound):
        min_weights = (1,1,*weights)
        terms = [[(e,c) for e,c in f.items() if weighted_degree(e,min_weights) <= bound]
                 for f in polynomials]
        powers = []
        for index, value in enumerate(base):
            largest = max(e[index] for f in terms for e,c in f)
            known = [{zero:field.one()}]
            for _ in range(largest):
                known.append(trunc_product(known[-1], value, bound))
            powers.append(known)
        evaluated = []
        for f in terms:
            result = {}
            for e,c in f:
                coeffs = {zero:c}
                for j,a in enumerate(e):
                    coeffs = trunc_product(coeffs, powers[j][a], bound)
                for exponent,value in coeffs.items():
                    result[exponent] = result.get(exponent, 0) + value
            evaluated.append({e:c for e,c in result.items() if c})
        return evaluated

    r, h = {}, {}
    independent = [{tuple(int(j == i) for j in range(len(weights))):field.one()}
                   for i in range(len(weights))]
    for degree in range(1, cap+1):
        error = evaluate(chart[:2], [r,h,*independent], degree)
        for j, series in enumerate((r,h)):
            for i, poly in enumerate(error):
                for e,c in poly.items():
                    assert weighted_degree(e,weights) == degree
                    series[e] = series.get(e,0) - inverse[j,i]*c
    values = evaluate(chart, [r,h,*independent], cap)
    assert not any(values[:2])
    return determinant, values[2:]


def verify_weighted_component_intersection(rows):
    p = 331
    k = GF(p)
    D = PolynomialRing(k, names=("l","q","tau"))
    l,q,tau = D.gens()
    chart, field = intersection_chart(rows, keep_l=False)
    # Along the q-axis with r=h=l=tau=0 the *full* equations vanish.
    assert all(e[0] or e[1] or e[2] or e[4] for row in chart for e in row)
    weights = (1,1,2)
    det, jets = implicit_intersection_jets(chart,field,weights,10)
    assert k(det)
    columns = sorted({e for f in jets for e in f}, key=lambda e:(weighted_degree(e,weights),e))
    A = matrix(QQ,[[f.get(e,0) for e in columns] for f in jets])
    combinations = A.augment(identity_matrix(QQ,A.nrows())).rref()
    leading = []
    for combination in combinations:
        value = {e:c for e,c in zip(columns,combination[:len(columns)],strict=True) if c}
        assert value
        degree = min(weighted_degree(e,weights) for e in value)
        form = {e:c for e,c in value.items() if weighted_degree(e,weights)==degree}
        denominator = lcm([c.denominator() for c in form.values()])
        scale = QQ(denominator)/gcd([QQ(c*denominator) for c in form.values()])
        for c in combination[len(columns):]:
            k(c*scale)  # Every multiplier must be integral at 331.
        leading.append(D({e:k(c*scale) for e,c in form.items()}))
    assert D.ideal(leading).radical() == D.ideal(l*q,tau)
    # On the q-dominant axis, put l=q L and tau=q^2 S.
    # S9/q^6 solves L. The displayed leading form then supplies a unit
    # cubic coefficient in S, while both lower jets vanish by the
    # generic A-chart calculation with S7,S8,S9.
    assert min(weighted_degree(e,weights) for e in jets[0]) == 6
    assert k(jets[0][(1,5,0)])
    assert (8*l+q)*tau**3 in leading

    chart, field = intersection_chart(rows, keep_l=True)
    weights = (1,2)
    determinant, jets = implicit_intersection_jets(chart,field,weights,6)

    def local_evaluation(function):
        # Normalize denominator to 1 at l=0 and check regularity over Z_(331).
        f = field(function)
        constant = f.denominator()(0)
        assert constant
        for c in (f.denominator()/constant).list() + (f.numerator()/constant).list():
            k(c)
        return k(f(0))

    assert local_evaluation(determinant)
    Q = PolynomialRing(k,names=("q","tau"))
    x = field.gen()
    assert all(weighted_degree(e,weights)>=5 for f in jets for e in f)
    first = jets[0]
    form = Q({e:local_evaluation(c/x) for e,c in first.items() if weighted_degree(e,weights)==5})
    leading = [form]
    assert form
    for f in jets[1:]:
        ratio = f[(5,0)]/first[(5,0)]
        local_evaluation(ratio)
        assert all(f.get(e,0)==ratio*c for e,c in first.items() if weighted_degree(e,weights)==5)
        assert all(not c or e in first for e,c in f.items() if weighted_degree(e,weights)==5)
        combo = {
            e:local_evaluation((f.get(e,0)-ratio*first.get(e,0))/x)
            for e in set(f)|set(first) if weighted_degree(e,weights)==6
        }
        leading.append(Q(combo))
    assert Q.ideal(leading).radical() == Q.ideal(Q.gens())
    print("At 331: exact jets exclude the intersection of the first two weighted components")


if __name__ == "__main__":
    residuals = reduced_residual_rows(5, 16, even_linear=True, primitive=False)
    verify_boundary(residuals)
    verify_first_weighted_component(residuals)
    verify_second_weighted_component(residuals)
    verify_weighted_component_intersection(residuals)
