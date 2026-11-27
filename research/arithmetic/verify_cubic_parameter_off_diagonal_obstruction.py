"""Verify the diagonal restriction and the exact singularity exclusion in CP6."""

from sage.all import Matrix, NumberField, PolynomialRing, QQ, gcd


def interpolation_basis():
    field_polynomials = PolynomialRing(QQ, "z")
    z = field_polynomials.gen()
    field = NumberField(z**4 - 2 * z**2 - QQ(1) / 3, "alpha")
    alpha = field.gen()
    pairs = PolynomialRing(field, ("u", "v"))
    u, v = pairs.gens()

    def branch_quotient(t):
        return t**3 + alpha * t**2 + (alpha**2 - 2) * t + alpha**3 - 2 * alpha

    plane = PolynomialRing(QQ, ("p", "q"))
    p, q = plane.gens()
    power_sums = (2, p, p**2 - 2 * q, p**3 - 3 * p * q)
    pair_product = branch_quotient(u) * branch_quotient(v)
    return [
        sum((QQ((value * alpha**weight).trace()) * q**j
             * (power_sums[i - j] if i > j else 1)
             for (i, j), value in pair_product.dict().items() if i >= j), plane.zero())
        for weight in range(4)
    ]


def check_off_diagonal_criterion():
    basis = interpolation_basis()
    single = PolynomialRing(QQ, "t")
    t = single.gen()
    diagonal = [single(form(2 * t, t**2)) for form in basis]
    coefficients = Matrix(QQ, [[form[i] for i in range(7)] for form in diagonal])
    relations = Matrix(QQ, [
        [1, 0, 0, 0, -QQ(1) / 9, 0, QQ(2) / 9],
        [0, 1, 0, 1, 0, -QQ(1) / 3, 0],
        [0, 0, 1, 0, QQ(2) / 3, 0, -QQ(7) / 3],
    ])
    assert coefficients.rank() == 4 and relations.rank() == 3
    assert coefficients * relations.transpose() == 0

    quadratics = PolynomialRing(QQ, ("a", "b", "c"))
    a, b, c = quadratics.gens()
    powers = PolynomialRing(quadratics, "t")
    t = powers.gen()
    H = a * t**2 + b * t + c
    f = t**4 - 2 * t**2 - QQ(1) / 3
    constraints = Matrix(quadratics, [
        [sum(relations[i, k] * (H**2 * t**j)[k] for k in range(7)) for j in range(3)]
        for i in range(3)
    ])
    assert constraints.det() == -QQ(1) / 3 * (a**2 + b**2 + 2 * a * c - 3 * c**2) * f.resultant(H)

    parameters = PolynomialRing(QQ, "r")
    r = parameters.gen()
    specialization = quadratics.hom([1 - 3 * r**2, 4 * r, 1 + r**2], parameters)
    conic_constraints = constraints.apply_map(specialization)
    assert conic_constraints * Matrix(parameters, 3, 1, [r**2, -2 * r, 1]) == 0
    assert gcd(conic_constraints.minors(2)) == 1
    at_infinity = constraints.apply_map(quadratics.hom([-3, 0, 1], QQ))
    assert at_infinity.rank() == 2 and at_infinity.column(0) == 0

    family_ring = PolynomialRing(parameters, "t")
    t = family_ring.gen()
    H_r = (1 - 3 * r**2) * t**2 + 4 * r * t + 1 + r**2
    diagonal_target = 3 * H_r**2 * (t - r)**2
    gamma = coefficients.change_ring(parameters).solve_left(
        Matrix(parameters, 1, 7, [diagonal_target[i] for i in range(7)])
    ).row(0)
    family_plane = PolynomialRing(parameters, ("p", "q"))
    p, q = family_plane.gens()
    cubic = sum((family_plane(form) * family_plane(g) for g, form in zip(gamma, basis)), family_plane.zero())
    derivative = family_ring(cubic.derivative(q)(2 * t, t**2))
    f = t**4 - 2 * t**2 - QQ(1) / 3
    f_plus = t**4 + 2 * t**2 - QQ(1) / 3
    assert derivative(r) == -9 * f(r)**2
    assert f.resultant(H_r) == 64 * f(r) * f_plus(r)
    assert derivative.resultant(H_r) == 11664 * (r**2 - QQ(1) / 3)**2 * f(r)**2 * f_plus(r)**2

    special = family_plane({exp: value % (r**2 - QQ(1) / 3) for exp, value in cubic.dict().items()})
    assert special == QQ(16) / 9 * (3 * p**2 * q - 3 * q**2 - 6 * q + 1)
    infinity = family_plane({exp: value[6] for exp, value in cubic.dict().items()})
    assert infinity == 3 * (3 * p**2 * q - 3 * q**2 - 6 * q + 1)


if __name__ == "__main__":
    check_off_diagonal_criterion()
    print("Verified the equianharmonic interpolation, diagonal square criterion, and singularity resultants.")
