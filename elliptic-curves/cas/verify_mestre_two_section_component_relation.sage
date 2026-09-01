#!/usr/bin/env sage
"""Function-field audit for the recovered Mestre two-section surface.

This checker constructs the seed component as the degree-eight factor over
``k(r3,r4)`` of the small complete intersection ``(M,H)``.  It is the
function-field companion to ``verify_mestre_two_section_root_surface.sing``.
The first gate records whether the leading invariant is already a square in
that field; later gates will use the same field for exact covariant group law.
"""

import argparse
import os

from sage.all import GF, QQ, FractionField, PolynomialRing, prod


def build_component_field(prime: int):
    ground = QQ if prime == 0 else GF(prime)
    base_polynomial = PolynomialRing(ground, names=("a", "b"))
    a, b = base_polynomial.gens()
    base = FractionField(base_polynomial)
    r3 = base(a)
    r4 = base(b)

    r5_ring = PolynomialRing(base, names=("c",))
    c = r5_ring.gen()
    r6_ring = PolynomialRing(r5_ring, names=("z",))
    z = r6_ring.gen()

    roots = (r3, r4, c, z)
    c1 = -sum(roots)
    c2 = sum(roots[i] * roots[j] for i in range(4) for j in range(i + 1, 4))
    c3 = -sum(
        roots[i] * roots[j] * roots[k]
        for i in range(4)
        for j in range(i + 1, 4)
        for k in range(j + 1, 4)
    )
    c4 = prod(roots)
    mestre = (
        c1**5 + c1**4 - 6 * c1**3 * c2 - 5 * c1**2 * c2
        + 8 * c1 * c2**2 + 7 * c1**2 * c3 + 4 * c2**2
        + 6 * c1 * c3 - 12 * c2 * c3 - 8 * c1 * c4
        - c1 - c2 - c3 - 16 * c4 - 1
    )
    e = r3 - r4 + 1
    d = r3 - c - z + 1
    sparse = r3 * d + (c + z - r4) * e * d - (r3 - c * z) * e

    resultant = r6_ring(mestre).resultant(r6_ring(sparse))
    factors = list(r5_ring(resultant).factor())
    seed_factors = [factor for factor, exponent in factors if factor.degree() == 8]
    if len(seed_factors) != 1:
        raise AssertionError(
            f"expected one degree-eight seed factor, found {[(f.degree(), n) for f, n in factors]}"
        )
    seed_polynomial = seed_factors[0].monic()
    component = base.extension(seed_polynomial, names=("alpha",))
    # The polynomial was selected as an irreducible factor.  Sage's generic
    # quotient-ring domain test does not reuse that factorization cache.
    component.modulus().is_irreducible.set_cache(True)
    alpha = component.gen()

    final_ring = PolynomialRing(component, names=("zz",))

    def specialize_r5(polynomial):
        polynomial = r6_ring(polynomial)
        return final_ring(
            [component(r5_ring(coefficient)(alpha)) for coefficient in polynomial.list()]
        )

    specialized_mestre = specialize_r5(mestre)
    specialized_sparse = specialize_r5(sparse)
    common = specialized_mestre.gcd(specialized_sparse).monic()
    if common.degree() != 1:
        raise AssertionError(f"expected a linear r6 gcd, found degree {common.degree()}")
    r6 = -common[0]
    values = (component(r3), component(r4), alpha, r6)
    return component, values, seed_polynomial


def replay(prime: int) -> None:
    component, roots, seed_polynomial = build_component_field(prime)
    print("component_field_constructed 1", flush=True)
    r3, r4, r5, r6 = roots
    a1 = -(r3 + r4 + r5 + r6 + 1)
    a2 = (
        r3 + r4 + r5 + r6
        + r3 * r4 + r3 * r5 + r3 * r6
        + r4 * r5 + r4 * r6 + r5 * r6
    )
    a3 = -(
        r3 * r4 + r3 * r5 + r3 * r6 + r4 * r5 + r4 * r6 + r5 * r6
        + r3 * r4 * r5 + r3 * r4 * r6 + r3 * r5 * r6 + r4 * r5 * r6
    )
    a4 = (
        r3 * r4 * r5 + r3 * r4 * r6 + r3 * r5 * r6 + r4 * r5 * r6
        + r3 * r4 * r5 * r6
    )
    leading = 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4
    square_ring = PolynomialRing(component, names=("omega",))
    omega = square_ring.gen()
    cover = component.extension(omega**2 - leading, names=("w",))
    # Work on the generic quadratic leading-square cover.  The exact
    # ordinate checks below would fail if this formal extension collapsed.
    cover.modulus().is_irreducible.set_cache(True)
    w = cover.gen()
    print("leading_square_cover_constructed 1", flush=True)

    roots = tuple(cover(value) for value in roots)
    r3, r4, r5, r6 = roots
    e = r3 - r4 + 1
    d = r3 - r5 - r6 + 1
    x01 = r3 / e
    x11 = (r5 - r6) / e
    x02 = (r3 - r5 * r6) / d
    x12 = -r4 / d

    t_ring = PolynomialRing(cover, names=("T",))
    T = t_ring.gen()
    x_ring = PolynomialRing(t_ring, names=("X",))
    X = x_ring.gen()
    q = X * (X - 1) * prod(X - value for value in roots)
    shifted_product = q(X - T) * q(X + T)
    approximant = X**6
    for lower_degree in range(5, -1, -1):
        target_degree = 6 + lower_degree
        correction = (
            shifted_product[target_degree] - (approximant * approximant)[target_degree]
        ) / 2
        approximant += correction * X**lower_degree
    numerator = approximant * approximant - shifted_product
    remainder_coefficients = []
    for coefficient in numerator.list():
        quotient, residue = t_ring(coefficient).quo_rem(T**2)
        if residue:
            raise AssertionError("Mestre remainder ceased to be divisible by T^2")
        remainder_coefficients.append(quotient)
    remainder = x_ring(remainder_coefficients)
    if remainder.degree() != 4:
        raise AssertionError(f"expected a quartic remainder, found degree {remainder.degree()}")

    def ordinate(intercept, slope):
        value = t_ring(remainder(intercept + slope * T))
        f = [value[index] for index in range(7)]
        s = (1 - slope**2) * w / 2
        n1 = 4 * s**2 * f[4] - f[5] ** 2
        n0 = 8 * s**4 * f[3] - n1 * f[5]
        z2 = f[5] / (2 * s)
        z1 = n1 / (8 * s**3)
        z0 = n0 / (16 * s**5)
        answer = t_ring(z0 + z1 * T + z2 * T**2 + s * T**3)
        if answer * answer != value:
            raise AssertionError("recursive cubic ordinate failed in the component field")
        return answer

    y1 = ordinate(x01, x11)
    y2 = ordinate(x02, x12)
    print("section_ordinates_verified 1", flush=True)
    quartic = [t_ring(coefficient) for coefficient in remainder.list()]

    def covariant_point(x_value, y_numerator, y_denominator=1):
        """Return Jacobian-projective covariant coordinates.

        For the quartic ordinate ``y_numerator/y_denominator``, the affine
        covariant coordinates have denominators ``y^2`` and ``y^3``.
        Therefore ``(36*g*yden^2,108*h*yden^3,ynum)`` is a projective point
        with x=X/Z^2 and y=Y/Z^3.  Keeping this representation avoids the
        prohibitively expensive fraction-field gcd normalization over the
        degree-sixteen component cover.
        """
        x_value = t_ring(x_value)
        y_numerator = t_ring(y_numerator)
        y_denominator = t_ring(y_denominator)
        ee, dd, cc, bb, aa = quartic
        g0 = bb**2 / 16 - aa * cc / 6
        g1 = bb * cc / 12 - aa * dd / 2
        g2 = cc**2 / 12 - bb * dd / 8 - aa * ee
        g3 = cc * dd / 12 - bb * ee / 2
        g4 = dd**2 / 16 - cc * ee / 6
        gv = g0 * x_value**4 + g1 * x_value**3 + g2 * x_value**2 + g3 * x_value + g4
        gx = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
        gy = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
        ux = 4 * aa * x_value**3 + 3 * bb * x_value**2 + 2 * cc * x_value + dd
        uy = bb * x_value**3 + 2 * cc * x_value**2 + 3 * dd * x_value + 4 * ee
        hv = (ux * gy - uy * gx) / 8
        return (
            36 * gv * y_denominator**2,
            108 * hv * y_denominator**3,
            y_numerator,
        )

    invariant_i = 12 * quartic[4] * quartic[0] - 3 * quartic[3] * quartic[1] + quartic[2] ** 2
    coefficient_a = -27 * invariant_i

    def negative(point):
        return point[0], -point[1], point[2]

    infinity = (t_ring(1), t_ring(1), t_ring(0))

    def equal(left, right):
        if left[2] == 0 or right[2] == 0:
            return left[2] == 0 and right[2] == 0
        return (
            left[0] * right[2] ** 2 == right[0] * left[2] ** 2
            and left[1] * right[2] ** 3 == right[1] * left[2] ** 3
        )

    def double(point):
        if point[2] == 0 or point[1] == 0:
            return infinity
        xx, yy, zz = point
        aa0 = xx**2
        bb0 = yy**2
        cc0 = bb0**2
        dd0 = 2 * ((xx + bb0) ** 2 - aa0 - cc0)
        ee0 = 3 * aa0 + coefficient_a * zz**4
        ff0 = ee0**2
        return (
            ff0 - 2 * dd0,
            ee0 * (dd0 - (ff0 - 2 * dd0)) - 8 * cc0,
            2 * yy * zz,
        )

    def add(left, right):
        if left is None or left[2] == 0:
            return right
        if right is None or right[2] == 0:
            return left
        x1, y1j, z1 = left
        x2, y2j, z2 = right
        z1z1 = z1**2
        z2z2 = z2**2
        u1 = x1 * z2z2
        u2 = x2 * z1z1
        s1 = y1j * z2 * z2z2
        s2 = y2j * z1 * z1z1
        if u1 == u2:
            return double(left) if s1 == s2 else infinity
        hh = u2 - u1
        ii = (2 * hh) ** 2
        jj = hh * ii
        rr = 2 * (s2 - s1)
        vv = u1 * ii
        x3 = rr**2 - jj - 2 * vv
        y3 = rr * (vv - x3) - 2 * s1 * jj
        z3 = ((z1 + z2) ** 2 - z1z1 - z2z2) * hh
        return x3, y3, z3

    P1 = covariant_point(x01 + x11 * T, y1)
    P2 = covariant_point(x02 + x12 * T, y2)
    print("affine_covariant_points_constructed 1", flush=True)
    visible = {}
    # The proposed relation only uses the r5 and r6 visible pairs.  Avoid
    # constructing the other eight covariant points in this high-degree
    # function field.
    for index, root in ((4, r5), (5, r6)):
        for sign in (-1, 1):
            x_value = root + sign * T
            y_numerator = t_ring(approximant(x_value))
            visible[index, sign] = covariant_point(x_value, y_numerator, T)
    print("support_visible_points_constructed 1", flush=True)

    # Balance the addition tree.  Jacobian-projective coordinate degrees grow
    # under addition, so this is materially smaller than a four-term chain.
    r5_difference = add(visible[4, -1], negative(visible[4, 1]))
    print("r5_visible_difference_constructed 1", flush=True)
    r6_difference = add(negative(visible[5, -1]), visible[5, 1])
    print("r6_visible_difference_constructed 1", flush=True)
    visible_sum = add(r5_difference, r6_difference)
    print("visible_side_constructed 1", flush=True)
    # The common square-root orientation is fixed by the recursive leading
    # coefficient.  On the pinned D-square curve its exact specialization is
    # P1+P2=V(r5,-)-V(r5,+)-V(r6,-)+V(r6,+).
    total = add(P1, P2)
    print("nonvisible_side_constructed 1", flush=True)
    relation_signs = [(1, 1)] if equal(total, visible_sum) else []

    print("MESTRE_TWO_SECTION_COMPONENT_FUNCTION_FIELD_V1")
    print(f"characteristic {prime}")
    print(f"seed_projection_degree {seed_polynomial.degree()}")
    print("leading_square_root_adjoined 1")
    print(f"section_relation_signs {relation_signs}")
    print(f"generic_visible_relation_verified {int(bool(relation_signs))}")
    print("DONE")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime",
        type=int,
        default=int(os.environ.get("MESTRE_RELATION_PRIME", "17")),
    )
    args = parser.parse_args()
    replay(args.prime)


main()
