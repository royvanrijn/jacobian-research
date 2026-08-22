#!/usr/bin/env sage
"""Exact generic group relation on the diameter-235 two-section component.

This is deliberately a Sage rational-function computation, rather than a
Grobner calculation or an expansion of the two-section residual.  It rebuilds
the Mestre square root recursively over QQ(p,T), maps the two triangular
affine sections and the relevant visible points to the covariant short
Jacobian, and checks the relation found at the p=-294,T=2 specialization.
"""

R = PolynomialRing(QQ, names=("p", "T"))
p, T = R.gens()
K = R.fraction_field()
p, T = K(p), K(T)
SX = PolynomialRing(K, "X")
X = SX.gen()


def polynomial(value, coefficients):
    answer = K.zero()
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def setup():
    b = (p - 66) * (p + 54) * (p^2 + 18 * p + 456)
    k = 3 * p^2 + 4 * p + 1068
    c1 = -polynomial(p, (-7245936, -246096, 4704, -876, 79)) / (2 * b)
    c2 = 5 * polynomial(p, (14277841000704, 888999349248, -30000699648,
        2783206656, -278543520, -14802624, 293232, -40272, 1849)) / (16 * b^2)
    c3 = -polynomial(p, (-32402742765539106816, -2677382748374243328,
        191800563007140864, -2244816893422080, 657002839322880,
        128514275004672, -3626525084544, 276573028032, -9211864320,
        -761052280, 16418244, -1903778, 59329)) / (16 * b^3)
    c4 = 25 * (p - 26) * (p - 6) * (p + 6) * (p + 14) * (p^2 - 12 * p + 276) * k * (7 * p^2 - 204 * p - 2628) * (29 * p^3 + 378 * p^2 + 3132 * p + 177336) * (37 * p^3 - 126 * p + 8316 * p - 269352) / (64 * b^4)
    roots = [K.zero(), K.one(),
        8 - (p + 294) * (3 * p^3 - 314 * p^2 - 7356 * p - 161208) / (4 * b),
        k * (7 * p^2 - 204 * p - 2628) / (2 * b),
        (p + 14) * (37 * p^3 - 126 * p^2 + 8316 * p - 269352) / (4 * b),
        QQ(235) / 17 - 45 * (p + 294) * (p^3 - 118 * p^2 - 2292 * p - 57416) / (34 * b)]
    lines = [
        (-(29 * p^3 + 378 * p^2 + 3132 * p + 177336) / (3 * (p - 66) * k), -(13 * p^2 + 204 * p + 1908) / (3 * k)),
        (polynomial(p, (-2491329312, -167485968, -6454080, -770760, 9870, -1053, 53)) / (3 * b * k), -2 * (p - 6) * (p + 54) / (3 * k))]
    q = X * (X - 1) * (X^4 + c1 * X^3 + c2 * X^2 + c3 * X + c4)
    product = q(X - T) * q(X + T)
    approximant = X^6
    for lower in range(5, -1, -1):
        approximant += (product[6 + lower] - (approximant^2)[6 + lower]) * X^lower / 2
    remainder = approximant^2 - product
    failures = [degree for degree in range(5, 13) if remainder[degree]]
    if failures:
        raise AssertionError(f"the recursive Mestre square root failed at {failures}")
    quartic = [remainder[degree] / T^2 for degree in range(5)]
    root_d = 15 * (p - 26) * (p - 6) * (p + 6) * (p + 14) * k / ((p - 66) * (p + 54) * (p^2 + 18 * p + 456)^2)
    return approximant, remainder, quartic, roots, lines, root_d


def triangular_ordinate(remainder, intercept, slope, root_d):
    substituted = remainder(intercept + slope * T) / T^2
    f = [substituted[degree] for degree in range(7)]
    y3 = (1 - slope^2) * root_d / 2
    y2 = f[5] / (2 * y3)
    y1 = (f[4] - y2^2) / (2 * y3)
    y0 = (f[3] - 2 * y1 * y2) / (2 * y3)
    values = [y0, y1, y2, y3]
    if any(f[degree] != sum(values[left] * values[degree - left] for left in range(max(0, degree - 3), min(3, degree) + 1)) for degree in range(7)):
        raise AssertionError("a triangular ordinate failed")
    return y0 + y1 * T + y2 * T^2 + y3 * T^3


def covariant_point(quartic, x_value, y_value):
    e, d, c, b, a = quartic
    g0 = b^2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    g2 = c^2 / 12 - b * d / 8 - a * e
    g3 = c * d / 12 - b * e / 2
    g4 = d^2 / 16 - c * e / 6
    g_value = g0 * x_value^4 + g1 * x_value^3 + g2 * x_value^2 + g3 * x_value + g4
    g_x = 4 * g0 * x_value^3 + 3 * g1 * x_value^2 + 2 * g2 * x_value + g3
    g_y = g1 * x_value^3 + 2 * g2 * x_value^2 + 3 * g3 * x_value + 4 * g4
    u_x = 4 * a * x_value^3 + 3 * b * x_value^2 + 2 * c * x_value + d
    u_y = b * x_value^3 + 2 * c * x_value^2 + 3 * d * x_value + 4 * e
    h_value = (u_x * g_y - u_y * g_x) / 8
    return 36 * g_value / y_value^2, 108 * h_value / y_value^3


def add(coefficient_a, left, right):
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    slope = (3 * x_left^2 + coefficient_a) / (2 * y_left) if x_left == x_right else (y_right - y_left) / (x_right - x_left)
    x_sum = slope^2 - x_left - x_right
    return x_sum, -y_left + slope * (x_left - x_sum)


def signed_sum(coefficient_a, terms):
    answer = None
    for point, sign in terms:
        answer = add(coefficient_a, answer, point if sign > 0 else (point[0], -point[1]))
    return answer


def replay():
    approximant, remainder, quartic, roots, lines, root_d = setup()
    visible = []
    for root in roots:
        for sign in (-1, 1):
            x_value = root + sign * T
            visible.append(covariant_point(quartic, x_value, approximant(x_value) / T))
    affine = []
    for intercept, slope in lines:
        ordinate = triangular_ordinate(remainder, intercept, slope, root_d)
        affine.append(covariant_point(quartic, intercept + slope * T, ordinate))
    invariant_i = 12 * quartic[4] * quartic[0] - 3 * quartic[3] * quartic[1] + quartic[2]^2
    coefficient_a = -27 * invariant_i
    # Orientations at p=-294,T=2 turn the discovered specialization vector
    # into this raw-covariant relation.
    relation = (-1, 0, 1, -1, -1, 1, 0, 1, 1, 0, 1, 0)
    left = add(coefficient_a, affine[0], affine[1])
    right = signed_sum(coefficient_a, [(point, sign) for point, sign in zip(visible, relation) if sign])
    if left != right:
        raise AssertionError("the generic diameter-235 relation failed")
    return {"status": "exact generic relation verified", "relation_vector": list(relation)}


if __name__ == "__main__":
    print(replay())
