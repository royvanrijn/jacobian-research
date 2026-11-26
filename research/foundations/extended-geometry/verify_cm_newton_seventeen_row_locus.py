"""Exact finite-algebra and lower-window checks for the three-pole system."""
from collections import defaultdict

from sage.all import GF, PolynomialRing, QQ, TermOrder, ZZ, matrix, vector

from verify_cm_newton_full_width_reconstruction_certificates import (
    coefficient_system,
    primitive_polynomial,
    reduced_residual_rows,
    verify_unit_combination,
)


def specialize_rows(rows, slice_name, ratio):
    if slice_name == "r=0":
        names = ("h", "l", "v", "w")
    elif slice_name == "h=17w":
        names = ("r", "l", "v", "w")
    elif slice_name == "ell=alpha*v":
        names = ("r", "u", "h", "w")
    else:
        names = ("r", "u", "w")
    ring = PolynomialRing(QQ, names, order="degrevlex")
    specialized = []
    for k, row in enumerate(rows, start=7):
        terms = defaultdict(lambda: QQ.zero())
        for (r, l, h, v, w), coefficient in row.dict().items():
            assert (l + v) % 2 == (k - 7) % 2
            if slice_name == "r=0":
                if r:
                    continue
                exponent = (h, l, v, w)
            elif slice_name == "h=17w":
                coefficient *= 17**h
                exponent = (r, l, v, h + w)
            else:
                coefficient *= ratio**l
                u = (l + v) // 2
                if slice_name == "ell=alpha*v":
                    exponent = (r, u, h, w)
                else:
                    coefficient *= 17**h
                    exponent = (r, u, h + w)
            terms[exponent] += coefficient
        poly = ring(terms)
        clearing = poly.denominator()
        integral = poly * clearing
        content = ZZ(integral.content())
        assert clearing % 331 and content % 331
        specialized.append(integral / content)
    return specialized


def exact_constant_span(polynomials, weights, parity, bound, expected):
    prime = 331
    table, nrows, ncols, target = coefficient_system(polynomials, weights, bound, parity)
    reduced = matrix(GF(prime), nrows, ncols, {key: GF(prime)(value) for key, value in table.items()})
    row_indices = reduced.transpose().pivots()
    col_indices = reduced.pivots()
    assert (nrows, ncols, len(col_indices)) == expected
    positions = {r: i for i, r in enumerate(row_indices)}
    exact = matrix(ZZ, len(row_indices), ncols, {
        (positions[i], j): ZZ(value) for (i, j), value in table.items() if i in positions
    }, sparse=True)
    rhs = vector(ZZ, [int(j == target) for j in col_indices])
    square = exact.matrix_from_columns(col_indices).transpose().dense_matrix()
    multiplier = square._solve_right_nonsingular_square(rhs, check_rank=False, algorithm="flint")
    assert all(QQ(value).denominator() % prime for value in multiplier)
    assert multiplier * exact == vector(QQ, [int(j == target) for j in range(ncols)])


def verify_affine_slices(rows):
    cases = [
        ("r=0", None, 16, 18, (2, 1, 1, 2), (1, 2), (1209, 1199, 792)),
        ("h=17w", None, 16, 18, (2, 1, 1, 2), (1, 2), (1209, 1203, 1080)),
        ("ell=alpha*v", QQ(8), 13, 8, (1, 1, 1, 1), (), (477, 484, 408)),
        ("ell=alpha*v", QQ(17)/2, 16, 9, (1, 1, 1, 1), (), (912, 715, 680)),
        ("ell=alpha*v", QQ(56)/5, 16, 9, (1, 1, 1, 1), (), (912, 715, 682)),
        ("ell=8v,h=17w", QQ(8), 11, 6, (1, 1, 1), (), (64, 68, 60)),
    ]
    for label, alpha, last, degree, weights, parity, shape in cases:
        transformed = specialize_rows(rows[:last - 6], label, alpha)
        exact_constant_span(transformed, weights, parity, degree, shape)
        print(f"Exact 331-integral unit identity on {label}, alpha {alpha}, through shift {last}")


def verify_rank_three_norm(rows):
    rows = rows[:3]
    # If wt(r,l,h,v,w)=(2,1,2,1,2), multiplication columns have
    # filtered weights 7,9,11 and determinant weight at most 21.
    weights = (2, 1, 2, 1, 2)
    assert all(max(sum(i*j for i,j in zip(e,weights,strict=True)) for e in g.dict()) == k-1
               for k,g in enumerate(rows,7))
    coefficient_of_r_cubed = next(c for e, c in rows[0].dict().items() if tuple(e) == (3, 0, 0, 0, 0))
    assert coefficient_of_r_cubed != 0

    # Specialize h=v=w=1. The characteristic-0 gcd checks below are
    # small univariate calculations after the rank-three reduction.
    base = PolynomialRing(QQ, "l")
    ell = base.gen()
    extension = PolynomialRing(base, "r")
    r = extension.gen()
    def specialize(row):
        terms = {}
        for (i, j, _, _, _), coefficient in row.dict().items():
            terms[i] = terms.get(i, base.zero()) + coefficient*ell**j
        return extension(terms)
    monic, second, third = map(specialize, rows)
    monic = monic.monic()
    second %= monic
    third %= monic
    def coordinates(poly):
        return vector(base, [poly[i] for i in range(3)])
    c0 = coordinates(second)
    c1 = coordinates(r*second % monic)
    c2 = coordinates(r*r*second % monic)
    norm = matrix([c0, c1, c2]).det()
    mixed = matrix([c0, c1, coordinates(third)]).det()
    assert (norm.degree(), mixed.degree()) == (21, 18)
    assert norm[0] and norm.gcd(norm.derivative()) == 1
    assert norm.gcd(monic[0]) == 1
    gcd, left, right = norm.xgcd(mixed)
    assert gcd == 1 and left*norm + right*mixed == 1

    residue = GF(65521)
    residue_base = PolynomialRing(residue, "l")
    reduced_norm = residue_base([residue(c) for c in norm])
    reduced_mixed = residue_base([residue(c) for c in mixed])
    assert list(map(int, [residue(coefficient_of_r_cubed),
                         reduced_norm.leading_coefficient(), reduced_mixed.leading_coefficient(),
                         reduced_norm[0], reduced_norm.resultant(reduced_mixed)])) == [
                             49912, 23322, 63396, 45123, 49856]
    print("Rank-three norm: 21 distinct common points for rows 7 and 8 at h=v=w=1, all excluded by row 9", flush=True)


def verify_rank_two_prefix(rows):
    ring = rows[0].parent()
    r, ell, j, v, w = ring.gens()
    s7, s8, s9, s10 = [row.subs({j: j + 17*w}) for row in rows[:4]]
    coefficient = lambda polynomial, exponent: polynomial.coefficient({w: exponent})
    k, L, C = [coefficient(s7, exponent) for exponent in (2, 1, 0)]
    assert s7.degree(w) == 2
    assert k == -QQ(12829696)/194828780125125 * r

    def linear_remainder(polynomial):
        assert polynomial.degree(w) <= 3
        for exponent in range(polynomial.degree(w), 1, -1):
            lead = coefficient(polynomial, exponent)
            if lead:
                polynomial -= ring(lead/k) * s7*w**(exponent-2)
        assert polynomial.degree(w) <= 1
        return coefficient(polynomial, 0), coefficient(polynomial, 1)

    a, d = linear_remainder(s8)
    x, y = linear_remainder(s9)
    norm = lambda constant, linear: k*constant**2 - L*constant*linear + C*linear**2
    n8, n9 = norm(a, d), norm(x, y)
    cross = a*y - d*x
    assert all(f.degree(w) <= 0 for f in (n8, n9, cross))
    assert (n8.degree(), len(n8.dict())) == (16, 744)
    assert (cross.degree(), len(cross.dict())) == (13, 383)
    degenerate = a.resultant(d, j)
    assert (degenerate.degree(), len(degenerate.dict())) == (17, 276)

    def reduction(prime):
        field = GF(prime)
        residue = PolynomialRing(field, names=ring.variable_names())
        return lambda polynomial: residue({
            exponent: field(c) for exponent, c in polynomial.dict().items()
        })

    reduce = reduction(65521)
    for polynomial in (n8, cross, degenerate):
        reduced = reduce(polynomial)
        assert reduced.degree() == polynomial.degree()
        assert len(reduced.factor()) == 1 and reduced.factor()[0][1] == 1
    jbar = reduce(j)
    second = reduce(d).resultant(reduce(n9), jbar)
    assert second.reduce([reduce(degenerate)])

    reduce = reduction(107)
    point = (1, 1, 37, 13, 0)
    mixed = norm(a+x, d+y)
    jacobian = n8.derivative(j)*n9.derivative(v) - n8.derivative(v)*n9.derivative(j)
    vals = [int(reduce(f)(*point)) for f in (n8, n9, cross, mixed, jacobian)]
    assert vals == [0, 0, 33, 55, 69]

    # Prefix elimination in j over the three-parameter base Q[r,ell,v].
    weights = (2, 1, 2, 1, 2)

    def weighted_degree(polynomial):
        return max(sum(weight*p for weight, p in zip(weights, exponent))
                   for exponent in polynomial.dict())

    n_lead = n8.coefficient({j: 7})
    d_lead = cross.coefficient({j: 6})
    assert n8.degree(j) == 7 and cross.degree(j) == 6
    assert n_lead and ring(n_lead/(ell*(ell-8*v))).is_constant()
    assert d_lead and ring(d_lead/ell).is_constant()
    assert (weighted_degree(n8), weighted_degree(cross),
            weighted_degree(degenerate)) == (16, 13, 17)

    E_top = ring({exponent: c for exponent, c in degenerate.dict().items()
                  if sum(weight*p for weight, p in zip(weights, exponent)) == 17})
    restricted = E_top.subs({r: 0})
    factor = QQ(restricted/(ell**2*(ell-8*v)**15))
    assert GF(331)(factor)

    # At this point of the degree-17 component mod 107 the common root
    # of a,d is unique and simple for D=a*y-d*x. Full degrees in a,d
    # and a nonzero slope in their linear Euclidean remainder certify
    # the same generic properties over Q.
    prime = 107
    reduce = reduction(prime)
    field = GF(prime)
    univariate = PolynomialRing(field, "u")
    u = univariate.gen()
    specialize = reduce(r).parent().hom([1, 1, u, 85, 0], univariate)
    aa, dd, delta = [specialize(reduce(f)) for f in (a, d, cross)]
    assert (aa.degree(), dd.degree(), delta.degree()) == (3, 2, 6)
    assert aa.gcd(dd) == u - 19 and delta.derivative()(19)
    reduced_E = reduce(degenerate)
    assert not reduced_E(1, 1, 19, 85, 0)
    assert reduced_E.derivative(reduce(v))(1, 1, 19, 85, 0)

    # The mod-p quotient exhibits the characteristic-zero factor F
    # without expanding a large rational resultant. Divisibility by ell
    # and E^2 over Q follows from the two leading coefficients and the
    # simple-root DVR argument in the accompanying proof.
    reduce = reduction(65521)
    jbar, lbar, wbar = [reduce(value) for value in (j, ell, w)]
    f = reduce(n8).resultant(reduce(cross), jbar)
    E = reduce(degenerate)
    F = f // (lbar*E**2)
    assert F*lbar*E**2 == f
    assert (weighted_degree(f), weighted_degree(F), F.degree(), len(F.dict())) == (103, 68, 63, 11742)
    assert len(F.factor()) == 1 and F.factor()[0][1] == 1
    assert reduce(cross).resultant(reduce(n9), jbar).reduce([E])

    r10 = reduce(s10)
    kbar = reduce(k)
    s7bar = reduce(s7)
    while r10.degree(wbar) > 1:
        d10 = r10.degree(wbar)
        r10 = (kbar*r10
               - r10.coefficient({wbar: d10}) * s7bar * wbar**(d10 - 2))
    a10, b10 = [r10.coefficient({wbar: i}) for i in (0, 1)]
    cross810 = reduce(a)*b10 - reduce(d)*a10
    assert reduce(n8).resultant(cross810, jbar).reduce([F])
    print("Rank-two prefixes: two exact component cuts and failure of separate norms")


def verify_finite_even_projection(rows):
    weights = (2, 1, 2)
    term_order = TermOrder("wdegrevlex", weights)
    for field in (QQ, GF(331)):
        ring = PolynomialRing(field, names=("r", "l", "h"), order=term_order)
        leading = [
            ring({tuple(e[:3]): field(c) for e, c in row.dict().items()
                  if not e[3] and not e[4]
                  and sum(i*j for i, j in zip(e[:3], weights, strict=True)) == k-1})
            for k, row in enumerate(rows[:3], start=7)
        ]
        assert all(f and max(sum(i*j for i, j in zip(e, weights, strict=True)) for e in f.dict()) == k-1 for k, f in enumerate(leading, start=7))
        ideal = ring.ideal(leading)
        assert ideal.dimension() == 0 and ideal.vector_space_dimension() == 84

    field = GF(331)
    specialized = [
        ring({tuple(e[:3]): field(c) for e, c in row.dict().items() if not e[3] and not e[4]})
        for row in (*rows[:3], rows[4])
    ]
    assert ring.ideal(specialized).groebner_basis() == [1]
    at_one = [
        sum(ring({tuple(e[:3]): field(c)}) for e, c in row.dict().items())
        for row in rows[:4]
    ]
    assert ring.ideal(at_one).groebner_basis() == [1]
    print("Finite projection: rank 84, and unit norms for rows 11 at the origin and 10 at (1,1)")


def verify_horizontal_chart():
    field = GF(331)
    slopes = [field(0), field(1), field(2)/17, field(5)/56]
    values = [field(1), field(0), *((8*v-1)/17 for v in slopes[2:])]
    g = [v*v+w for v, w in zip(slopes, values, strict=True)]
    assert all(g)
    assert [w/a for w, a in zip(values, g, strict=True)] == [1, 0, 110, 205]
    print("One horizontal chart separates the four residue clusters")


def verify_row_twelve_witness(rows):
    field = GF(23)
    point = (17, 5, 12, 6, 17)
    ring = PolynomialRing(field, names=rows[0].parent().variable_names())
    polynomials = [
        ring({e: field(c) for e, c in primitive_polynomial(row).dict().items()})
        for row in rows[:6]
    ]
    assert not any(f(*point) for f in polynomials[:5])
    jacobian = matrix(field, [[f.derivative(x)(*point) for x in ring.gens()]
                              for f in polynomials[:5]])
    assert jacobian.det() == 18 and polynomials[5](*point) == 21
    print("Simple algebraic blind branch through shift 11, detected by shift 12")


def verify_angular_zero_slice(rows):
    rows = rows[:7]
    ring = PolynomialRing(QQ, names=("r", "h", "w"))
    r, h, w = ring.gens()
    specialize = rows[0].parent().hom([r, 0, h, 0, w], ring)
    assert all(not specialize(row) for row in rows[1::2])
    even = [primitive_polynomial(specialize(row)) for row in rows[::2]]

    residue_field = GF(11)
    residue_ring = PolynomialRing(residue_field, names=ring.variable_names())
    reduced = [residue_ring(poly) for poly in even]
    point = (residue_field(2), residue_field(8), residue_field(2))
    assert all(not poly(*point) for poly in reduced[:3])
    jacobian = matrix(residue_field, [
        [poly.derivative(var)(*point) for var in residue_ring.gens()]
        for poly in reduced[:3]
    ])
    assert jacobian.det() == 7 and reduced[3](*point) == 10
    print("Hensel lift of (2,8,2) mod 11 on ell=v=0 is blind through row 12, with row 13 nonzero")
    assert verify_unit_combination(even, (1, 1, 1), 12, ()) == 441
    print("Exact unit on ell=v=0 from rows 7,9,11,13 at degree 12")


if __name__ == "__main__":
    residuals = reduced_residual_rows(5, 16, even_linear=True, primitive=False)
    verify_rank_three_norm(residuals)
    verify_rank_two_prefix(residuals)
    verify_finite_even_projection(residuals)
    verify_horizontal_chart()
    verify_row_twelve_witness(residuals)
    verify_angular_zero_slice(residuals)
    verify_affine_slices(residuals)

