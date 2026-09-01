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

from sage.all import (
    GF,
    QQ,
    FractionField,
    FunctionField,
    PolynomialRing,
    matrix,
    prod,
    vector,
)


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
    full_a1 = -(1 + sum(roots))
    full_a2 = sum(roots) + c2
    full_a3 = -(c2 + sum(
        roots[i] * roots[j] * roots[k]
        for i in range(4)
        for j in range(i + 1, 4)
        for k in range(j + 1, 4)
    ))
    full_a4 = (
        sum(
            roots[i] * roots[j] * roots[k]
            for i in range(4)
            for j in range(i + 1, 4)
            for k in range(j + 1, 4)
        )
        + c4
    )
    raw_leading = (
        5 * full_a1**4 - 24 * full_a1**2 * full_a2
        + 32 * full_a1 * full_a3 + 16 * full_a2**2 - 64 * full_a4
    )
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
    specialized_leading = specialize_r5(raw_leading)
    leading = specialized_leading(r6)
    return component, values, seed_polynomial, leading


def replay(prime: int) -> None:
    component, roots, seed_polynomial, leading = build_component_field(prime)
    print("component_field_constructed 1", flush=True)
    r3, r4, r5, r6 = roots
    if os.environ.get("MESTRE_QUADRATIC_TOWER") == "1":
        base = component.base_ring()
        alpha = component.gen()

        def component_coordinates(value):
            coefficients = list(component(value).lift())
            coefficients += [base(0)] * (8 - len(coefficients))
            return vector(base, coefficients[:8])

        delta_powers = [component(1)]
        for _ in range(3):
            delta_powers.append(delta_powers[-1] * leading)

        # Write alpha^2-S(delta)*alpha+P(delta)=0.  This is a linear
        # calculation in the alpha power basis, despite the quadratic tower
        # it uncovers.
        columns = [
            component_coordinates(-power * alpha) for power in delta_powers
        ] + [component_coordinates(power) for power in delta_powers]
        tower_matrix = matrix(
            base, 8, 8, lambda row, column: columns[column][row]
        )
        solution = tower_matrix.solve_right(-component_coordinates(alpha**2))
        trace = sum(solution[index] * delta_powers[index] for index in range(4))
        norm = sum(
            solution[4 + index] * delta_powers[index] for index in range(4)
        )
        if alpha**2 - trace * alpha + norm != 0:
            raise AssertionError("quadratic tower relation did not verify")

        discriminant = trace**2 - 4 * norm
        quotient = discriminant / leading
        subfield_columns = [component_coordinates(power) for power in delta_powers]
        subfield_matrix = matrix(
            base, 8, 4, lambda row, column: subfield_columns[column][row]
        )
        pivot_rows = subfield_matrix.transpose().pivots()
        square_submatrix = subfield_matrix.matrix_from_rows(pivot_rows)

        def subfield_coordinates(value):
            target = component_coordinates(value)
            answer = square_submatrix.solve_right(
                vector(base, [target[index] for index in pivot_rows])
            )
            if subfield_matrix * answer != target:
                raise AssertionError("element did not lie in the degree-four subfield")
            return answer

        quotient_coordinates = subfield_coordinates(quotient)
        leading_coordinates = subfield_coordinates(leading)
        minimal_polynomial = leading.minpoly()

        # If q has degree four, then q is a square in F(q) exactly when
        # minpoly_q(Y^2) has a degree-four factor.  For such a factor
        # Y^4+A Y^3+B Y^2+C Y+E, its root is recovered inside F(q) by
        #     sqrt(q)=-(q^2+Bq+E)/(Aq+C).
        # This avoids the much larger four-coefficient square-root system.
        y_ring = PolynomialRing(base, names=("Y",))
        Y = y_ring.gen()

        def square_via_even_minpoly(value):
            value_minpoly = value.minpoly()
            even_polynomial = sum(
                base(coefficient) * Y ** (2 * index)
                for index, coefficient in enumerate(value_minpoly.list())
            )
            factorization = list(even_polynomial.factor())
            degrees = [factor.degree() for factor, _ in factorization]
            candidate = None
            candidate_factor = None
            for factor, _ in factorization:
                if factor.degree() != value_minpoly.degree():
                    continue
                factor = factor.monic()
                if factor.degree() != 4:
                    continue
                aa = factor[3]
                bb = factor[2]
                cc = factor[1]
                ee = factor[0]
                denominator = aa * value + cc
                if denominator == 0:
                    continue
                proposed = -(value**2 + bb * value + ee) / denominator
                if proposed**2 == value:
                    candidate = proposed
                    candidate_factor = factor
                    break
            return value_minpoly, degrees, candidate, candidate_factor

        (
            quotient_minpoly,
            quotient_factor_degrees,
            quotient_square_root,
            quotient_square_factor,
        ) = square_via_even_minpoly(quotient)
        (
            leading_minpoly,
            leading_factor_degrees,
            leading_square_root,
            leading_square_factor,
        ) = square_via_even_minpoly(leading)
        quotient_is_square = quotient_square_root is not None
        leading_is_square = leading_square_root is not None
        component_square_root = leading_square_root
        if component_square_root is None and quotient_square_root is not None:
            component_square_root = (2 * alpha - trace) / quotient_square_root
        component_is_square = (
            component_square_root is not None
            and component_square_root**2 == leading
        )

        print("MESTRE_TWO_SECTION_QUADRATIC_TOWER_V1")
        print(f"characteristic {prime}")
        print(f"component_degree {seed_polynomial.degree()}")
        print(f"subfield_degree {minimal_polynomial.degree()}")
        print(f"tower_matrix_rank {tower_matrix.rank()}")
        print("quadratic_relation_verified 1")
        print(
            "discriminant_over_D_square_in_subfield",
            int(quotient_is_square),
        )
        print(
            "D_square_in_subfield",
            int(leading_is_square),
        )
        print("discriminant_over_D_minpoly_degree", quotient_minpoly.degree())
        print("quotient_even_factor_degrees", quotient_factor_degrees)
        print("D_minpoly_degree", leading_minpoly.degree())
        print("leading_even_factor_degrees", leading_factor_degrees)
        print("D_square_in_component", int(component_is_square))
        print("component_square_root_verified", int(component_is_square))
        if os.environ.get("MESTRE_TOWER_VERBOSE") == "1":
            print("TRACE_COEFFICIENTS")
            for value in solution[:4]:
                print(value)
            print("NORM_COEFFICIENTS")
            for value in solution[4:]:
                print(value)
            print("DISCRIMINANT_OVER_D_COORDINATES")
            for value in quotient_coordinates:
                print(value)
            if quotient_square_factor is not None:
                print("QUOTIENT_SQUARE_FACTOR")
                print(quotient_square_factor)
            if leading_square_factor is not None:
                print("LEADING_SQUARE_FACTOR")
                print(leading_square_factor)
        if os.environ.get("MESTRE_TOWER_FACTOR_ONLY") == "1":
            print("LEADING_SQUARE_FACTOR")
            print(leading_square_factor)
        print("DONE")
        return
    if os.environ.get("MESTRE_MINPOLY_D") == "1":
        minimal_polynomial = leading.minpoly()
        print("MESTRE_TWO_SECTION_LEADING_MINPOLY_V1")
        print(f"characteristic {prime}")
        print(f"degree {minimal_polynomial.degree()}")
        print(minimal_polynomial)
        print("DONE")
        return
    if os.environ.get("MESTRE_SAMPLE_SQRT") == "1":
        if prime == 0:
            raise ValueError("finite-extension square-root samples require a prime")
        ground = GF(prime)
        sample_ring = PolynomialRing(ground, names=("c",))
        c_sample = sample_ring.gen()

        def specialize_base(value, first, second):
            value = component.base_ring()(value)
            numerator = value.numerator()(ground(first), ground(second))
            denominator = value.denominator()(ground(first), ground(second))
            if denominator == 0:
                raise ZeroDivisionError
            return ground(numerator) / ground(denominator)

        lifted_leading = leading.lift()
        requested = int(os.environ.get("MESTRE_SAMPLE_COUNT", "12"))
        interpolate = os.environ.get("MESTRE_INTERPOLATE_SQRT") == "1"
        found = 0
        sample_records = []

        def interpolate_records():
            if not interpolate:
                return
            display_ring = PolynomialRing(ground, names=("A", "B"))
            A, B = display_ring.gens()
            declared_training_count = int(
                os.environ.get("MESTRE_TRAINING_COUNT", str(len(sample_records)))
            )
            training_count = min(declared_training_count, len(sample_records))
            training_records = sample_records[:training_count]
            holdout_records = sample_records[training_count:]
            summary_only = os.environ.get("MESTRE_INTERPOLATE_SUMMARY") == "1"
            for ratio_index in range(7):
                reconstructed = False
                for degree in range(1, 41):
                    exponents = [
                        (first_degree, total_degree - first_degree)
                        for total_degree in range(degree + 1)
                        for first_degree in range(total_degree + 1)
                    ]
                    if 2 * len(exponents) > len(training_records) + 1:
                        break
                    rows = []
                    for first, second, ratios in training_records:
                        values = [
                            ground(first) ** first_degree * ground(second) ** second_degree
                            for first_degree, second_degree in exponents
                        ]
                        ratio = ground(ratios[ratio_index])
                        rows.append(values + [-ratio * value for value in values])
                    kernel = matrix(ground, rows).right_kernel()
                    if kernel.dimension() != 1:
                        continue
                    vector = list(kernel.basis()[0])
                    count = len(exponents)
                    numerator_coefficients = vector[:count]
                    denominator_coefficients = vector[count:]
                    pivot = next(value for value in denominator_coefficients if value)
                    numerator_coefficients = [value / pivot for value in numerator_coefficients]
                    denominator_coefficients = [value / pivot for value in denominator_coefficients]
                    numerator = sum(
                        coefficient * A**first_degree * B**second_degree
                        for coefficient, (first_degree, second_degree)
                        in zip(numerator_coefficients, exponents)
                    )
                    denominator = sum(
                        coefficient * A**first_degree * B**second_degree
                        for coefficient, (first_degree, second_degree)
                        in zip(denominator_coefficients, exponents)
                    )
                    holdout_failures = 0
                    for first, second, ratios in holdout_records:
                        numerator_value = numerator(ground(first), ground(second))
                        denominator_value = denominator(ground(first), ground(second))
                        if (
                            denominator_value == 0
                            or numerator_value
                            != ground(ratios[ratio_index]) * denominator_value
                        ):
                            holdout_failures += 1
                    if holdout_failures:
                        continue
                    print(f"RATIO {ratio_index} DEGREE {degree}")
                    print(f"HOLDOUTS {len(holdout_records)} FAILURES 0")
                    if summary_only:
                        print(f"NUMERATOR_TERMS {len(numerator.dict())}")
                        print(f"DENOMINATOR_TERMS {len(denominator.dict())}")
                    else:
                        print(f"NUMERATOR {numerator}")
                        print(f"DENOMINATOR {denominator}")
                    reconstructed = True
                    break
                if not reconstructed:
                    print(f"RATIO {ratio_index} UNRESOLVED")

        print("MESTRE_TWO_SECTION_SQUARE_ROOT_SAMPLES_V1")
        print(f"characteristic {prime}")
        for first in range(prime):
            for second in range(prime):
                try:
                    specialized_modulus = sample_ring(
                        [
                            specialize_base(value, first, second)
                            for value in component.modulus().list()
                        ]
                    )
                    if specialized_modulus.degree() != 8 or not specialized_modulus.is_irreducible():
                        continue
                    specialized_leading = sample_ring(
                        [
                            specialize_base(value, first, second)
                            for value in lifted_leading.list()
                        ]
                    )
                except ZeroDivisionError:
                    continue
                extension = GF(prime**8, names=("theta",), modulus=specialized_modulus)
                theta = extension.gen()
                leading_value = extension(
                    sum(
                        extension(specialized_leading[index]) * theta**index
                        for index in range(len(specialized_leading.list()))
                    )
                )
                if not leading_value.is_square():
                    print(f"NONSQUARE {first} {second}")
                    raise AssertionError("an irreducible specialization made D nonsquare")
                root = leading_value.sqrt(all=True)[0]
                root_coefficients = list(root.polynomial())
                root_coefficients += [ground(0)] * (8 - len(root_coefficients))
                normalization_index = max(
                    index for index, value in enumerate(root_coefficients) if value
                )
                normalization = root_coefficients[normalization_index]
                ratios = [int(value / normalization) for value in root_coefficients]
                sample_records.append((first, second, ratios))
                if not interpolate:
                    print(
                        "SAMPLE",
                        first,
                        second,
                        normalization_index,
                        *ratios,
                    )
                found += 1
                if requested > 0 and found >= requested:
                    interpolate_records()
                    print(f"sample_count {found}")
                    print("DONE")
                    return
        if requested == 0:
            interpolate_records()
            print(f"sample_count {found}")
            print("DONE")
            return
        raise AssertionError(f"found only {found} irreducible specializations")
    if os.environ.get("MESTRE_SOLVE_SQRT") == "1":
        base = component.base_ring()
        coefficient_ring = PolynomialRing(
            base,
            names=tuple(f"s{index}" for index in range(8)),
            order="degrevlex",
        )
        coefficients = coefficient_ring.gens()
        alpha_ring = PolynomialRing(coefficient_ring, names=("c",))
        c = alpha_ring.gen()
        modulus = alpha_ring(
            [coefficient_ring(value) for value in component.modulus().list()]
        )
        lifted_leading = leading.lift()
        leading_polynomial = alpha_ring(
            [coefficient_ring(value) for value in lifted_leading.list()]
        )
        candidate = sum(coefficients[index] * c**index for index in range(8))
        remainder = (candidate**2 - leading_polynomial).quo_rem(modulus)[1]
        equations = [remainder[index] for index in range(8)]
        square_root_ideal = coefficient_ring.ideal(equations)
        print("MESTRE_TWO_SECTION_SQUARE_ROOT_SYSTEM_V1")
        print(f"characteristic {prime}")
        print(f"equation_count {len(equations)}")
        basis = square_root_ideal.groebner_basis()
        basis_ideal = coefficient_ring.ideal(basis)
        print(f"ideal_dimension {basis_ideal.dimension()}", flush=True)
        print(f"basis_size {len(basis)}")
        for value in basis:
            print(value)
        print("DONE")
        return
    if os.environ.get("MESTRE_COMPONENT_SQRT") == "1":
        split_a = 3 * r3**2 - 3 * r4**2 - 6 * r3 + 3
        split_b = (
            -r3**4 + 6 * r3**3 * r4 - 7 * r3**2 * r4**2
            + 6 * r3 * r4**3 - r4**4 - 8 * r3**3 - 6 * r3**2 * r4
            + 2 * r3 * r4**2 + 6 * r4**3 + 18 * r3**2
            - 6 * r3 * r4 - 7 * r4**2 - 8 * r3 + 6 * r4 - 1
        )
        split_c = (
            12 * r3**4 * r4**2 - 12 * r3**2 * r4**4
            - 48 * r3**3 * r4**2 + 24 * r3 * r4**4
            + 72 * r3**2 * r4**2 - 12 * r4**4
            - 48 * r3 * r4**2 + 12 * r4**2
        )
        split_e = (
            -144 * r3**4 * r4**3 + 144 * r3**3 * r4**4
            + 144 * r3**4 * r4**2 + 144 * r3**3 * r4**3
            - 288 * r3**2 * r4**4 - 288 * r3**3 * r4**2
            + 144 * r3**2 * r4**3 + 144 * r3 * r4**4
            + 144 * r3**2 * r4**2 - 144 * r3 * r4**3
        )
        w = -(leading**2 + split_b * leading + split_e) / (
            split_a * leading + split_c
        )
        if w**2 != leading:
            raise AssertionError("component square-root formula failed")
        cover = component
        print("leading_square_in_component_verified 1", flush=True)
    else:
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
