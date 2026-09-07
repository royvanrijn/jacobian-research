#!/usr/bin/env sage
"""Fail-closed local audit for the missing determinant-1236 double cover.

Let C be X_0^6(103)/<w_618>.  The currently recorded Atkin--Lehner
assignment identifies B=C/<w_2>, E=C/<w_2,w_3>, and a genus-three quotient
D=C/<w_3> whose Prym has the 618c1 and 618d1 factors.  Under that assignment,
C is the normalization of B x_E D, so a rational point of B can lift to C
only if its image on E lifts to D.

The independent exhaustive V4 audit confirms that this factor assignment is
locally compatible at p=5 and p=7 with the displayed B -> 618f1 squareclass.
This script remains a local reconstruction screen: none of its local
survivors should be promoted to a characteristic-zero cover without a global
squareclass reconstruction.

For each requested good prime this script enumerates an exhaustive superset
of the possible quadratic extensions of F_p(E) having:

* rational branch point +3G or -3G;
* three further geometric branch points;
* Frobenius power traces of 618c1 x 618d1 through degrees 1, 2, and 3.

It then reports the possible quadratic characters at the six non-fixed
rational images +/-G, +/-4G, +/-10G.  A singleton [-1] is a local
obstruction.  Sets containing +1 are inconclusive.  An empty candidate set
is also inconclusive.  The enumeration includes both even-pole charts
(``L(6O)`` and ``L(4O)``) and the charts where one branch point is ``O``
(``L(5O)`` and ``L(3O)``).

The value ``delta`` below is a quadratic CHARACTER (+1 for the current
twist, -1 for its nonsquare twist), not the literal field element -1.  This
distinction is essential when p == 1 (mod 4).

Quick replay:
    sage elkies-k3/scripts/explore_det1236_double_cover_local_gate.sage 5 7 11 13

Append the positional word ``details`` to emit every surviving local
coefficient vector and its chart metadata.

Append ``cmfield`` to require the residual branch factorization dictated by
``a^3-a^2+4*a+12``.  This is valid only after that cubic has independently
been certified as the residue field of the discriminant ``-1236`` orbit.

At every even zero, including the imposed double zero Q, character sums use
the leading local unit after removing the square uniformizer factor.  Merely
evaluating the function to zero there would give an incorrect point count.

The implementation is deliberately exhaustive rather than optimized.  Use
small primes for routine replays; the degree-three point enumeration becomes
expensive quickly.

Extended exploratory replay:
    sage elkies-k3/scripts/explore_det1236_double_cover_local_gate.sage \
        5 7 11 13 17 23 29 31 37 41 43 47
"""

import json
import sys
from itertools import product


EQ = EllipticCurve(QQ, [1, 0, 0, -185, 1401])
ELLIPTIC_FACTORS = {
    letter: EllipticCurve("618%s1" % letter) for letter in "abcd"
}
PRYM_FACTORS = {
    pair: tuple(ELLIPTIC_FACTORS[letter] for letter in pair)
    for pair in ("ab", "ac", "ad", "bc", "bd", "cd")
}
COMPLEMENTARY_PRYM_PAIR = {
    "ab": "cd",
    "ac": "bd",
    "ad": "bc",
    "bc": "ad",
    "bd": "ac",
    "cd": "ab",
}


def power_trace(ap, p, degree):
    if degree == 0:
        return ZZ(2)
    if degree == 1:
        return ZZ(ap)
    previous, current = ZZ(2), ZZ(ap)
    for _ in range(2, degree + 1):
        previous, current = current, ap * current - p * previous
    return ZZ(current)


def target_triple(p, prym_pair):
    return tuple(
        -power_trace(PRYM_FACTORS[prym_pair][0].ap(p), p, degree)
        - power_trace(PRYM_FACTORS[prym_pair][1].ap(p), p, degree)
        for degree in (1, 2, 3)
    )


def quadratic_character(value):
    if not value:
        return 0
    return 1 if value.is_square() else -1


def normalize(vector_value):
    for value in vector_value:
        if value:
            return tuple(value ** (-1) * entry for entry in vector_value)
    raise ArithmeticError("zero projective vector")


def projective_vectors(field, vector_dimension):
    """Enumerate each point of P^(vector_dimension-1)(field) once."""
    for pivot in range(vector_dimension):
        for tail in product(field, repeat=vector_dimension-pivot-1):
            yield vector(field, [0]*pivot + [1] + list(tail))


def evaluation_row(point, field):
    x_value, y_value = point[0], point[1]
    return vector(
        field,
        [1, x_value, y_value, x_value**2, x_value * y_value, x_value**3],
    )


def derivative_row(point, field):
    x_value, y_value = point[0], point[1]
    denominator = 2 * y_value + x_value
    numerator = 3 * x_value**2 - 185 - y_value
    if not denominator:
        # At a geometric 2-torsion point the tangent is vertical.  With y
        # as local parameter, x has no linear term.
        return vector(field, [0, 0, 1, 0, x_value, 0])
    return vector(
        field,
        [
            0,
            denominator,
            numerator,
            2 * x_value * denominator,
            y_value * denominator + x_value * numerator,
            3 * x_value**2 * denominator,
        ],
    )


def local_basis_expansion(point, field, precision=7):
    """Expand 1,x,y,x^2,xy,x^3 in a local uniformizer."""
    series_ring = PowerSeriesRing(field, "t", default_prec=precision)
    t_value = series_ring.gen()
    x_zero, y_zero = point[0], point[1]
    partial_y = 2*y_zero+x_zero
    if partial_y:
        x_series = series_ring(x_zero)+t_value
        y_series = series_ring(y_zero)
        derivative = partial_y
        solve_for_y = True
    else:
        y_series = series_ring(y_zero)+t_value
        x_series = series_ring(x_zero)
        derivative = y_zero-3*x_zero**2+field(185)
        assert derivative
        solve_for_y = False

    def equation_value(x_argument, y_argument):
        return (
            y_argument**2+x_argument*y_argument
            - x_argument**3+field(185)*x_argument-field(1401)
        )

    for degree in range(1, precision):
        error = equation_value(x_series, y_series)[degree]
        correction = -error/derivative
        if solve_for_y:
            y_series += correction*t_value**degree
        else:
            x_series += correction*t_value**degree
    final_error = equation_value(x_series, y_series)
    assert all(final_error[degree] == 0 for degree in range(precision))
    basis = [
        series_ring(1), x_series, y_series, x_series**2,
        x_series*y_series, x_series**3,
    ]
    return [
        [basis[index][degree] for index in range(6)]
        for degree in range(precision)
    ]


def local_order_and_leading(coefficients, expansion_rows):
    for order, row in enumerate(expansion_rows):
        leading = sum(
            coefficients[index]*row[index] for index in range(6)
        )
        if leading:
            return order, leading
    raise ArithmeticError(
        "nonzero element of L(6O) vanished to order greater than six"
    )


def local_character(coefficients, expansion_rows):
    order, leading = local_order_and_leading(coefficients, expansion_rows)
    return 0 if order % 2 else quadratic_character(leading)


def residual_branch_splitting(
    coefficients, field, rational_branch, q_point, pole_degree
):
    polynomial_ring = PolynomialRing(field, "x")
    x_value = polynomial_ring.gen()
    a_value = (
        coefficients[0]
        + coefficients[1] * x_value
        + coefficients[3] * x_value**2
        + coefficients[5] * x_value**3
    )
    b_value = coefficients[2] + coefficients[4] * x_value
    cubic = x_value**3 - 185 * x_value + 1401
    norm = a_value**2 - x_value * a_value * b_value - cubic * b_value**2
    prescribed = x_value - rational_branch[0]
    if q_point is not None:
        prescribed *= (x_value - q_point[0]) ** 2
    quotient, remainder = norm.quo_rem(prescribed)
    # An even pole leaves three residual finite branch points.  An odd pole
    # means that O itself is a branch point and leaves two.  Squarefreeness of
    # the x-polynomial is deliberately not required: distinct opposite points
    # on E share an x-coordinate.  This is a superset screen, so retaining a
    # degenerate class cannot create a false local exclusion.
    residual_degree = 3 if pole_degree % 2 == 0 else 2
    if remainder != 0 or quotient.degree() != residual_degree:
        return None
    return tuple(
        sorted(
            factor.degree()
            for factor, multiplicity in quotient.factor()
            for _ in range(multiplicity)
        )
    )


def character_sum(coefficients, cache, pole_degree):
    field, local_data = cache
    lifted = [field(value) for value in coefficients]
    leading_index = {6: 5, 5: 4, 4: 3, 3: 2}[pole_degree]
    leading = lifted[leading_index]
    if not leading:
        return None
    # At an even pole there are two or zero points above O according to the
    # leading squareclass.  At an odd pole O is ramified and contributes zero
    # to #D-#E.
    answer = quadratic_character(leading) if pole_degree % 2 == 0 else 0
    for _, expansion_rows, _, _ in local_data:
        answer += local_character(lifted, expansion_rows)
    return ZZ(answer)


def product_with_known_cover_character_sum(coefficients, cache, pole_degree):
    """Character sum after multiplying by (X-4)/54, the B -> E cover."""
    field, local_data = cache
    lifted = [field(value) for value in coefficients]
    leading_index = {6: 5, 5: 4, 4: 3, 3: 2}[pole_degree]
    leading = lifted[leading_index]
    if not leading:
        return None
    answer = (
        quadratic_character(leading / field(54))
        if pole_degree % 2 == 0
        else 0
    )
    for _, expansion_rows, h_order, h_leading in local_data:
        b_order, b_leading = local_order_and_leading(
            lifted, expansion_rows
        )
        if (b_order+h_order) % 2 == 0:
            answer += quadratic_character(b_leading*h_leading)
    return ZZ(answer)


def enumerate_prime(
    p,
    include_candidates=False,
    require_cm_field_splitting=False,
    prym_pair="cd",
    require_v4_compatibility=False,
    all_rational_branches=False,
):
    field = GF(p)
    elliptic_curve = EllipticCurve(field, [1, 0, 0, -185, 1401])
    generator = elliptic_curve(10, -29)
    target = target_triple(p, prym_pair)
    complementary_pair = COMPLEMENTARY_PRYM_PAIR[prym_pair]
    complementary_target = target_triple(p, complementary_pair)
    polynomial_ring = PolynomialRing(field, "u")
    u_value = polynomial_ring.gen()
    cm_field_polynomial = u_value**3-u_value**2+4*u_value+12
    cm_field_splitting = tuple(
        sorted(
            factor.degree()
            for factor, multiplicity in cm_field_polynomial.factor()
            for _ in range(multiplicity)
        )
    )

    caches = {}
    for degree in (1, 2, 3):
        extension = GF(p**degree, "z")
        extended_curve = EllipticCurve(extension, [1, 0, 0, -185, 1401])
        local_data = []
        for point in extended_curve:
            if not point.is_zero():
                expansion_rows = local_basis_expansion(point, extension)
                h_order, h_leading = local_order_and_leading(
                    [extension(-4), extension(1), 0, 0, 0, 0],
                    expansion_rows,
                )
                local_data.append(
                    (
                        point,
                        expansion_rows,
                        h_order,
                        h_leading/extension(54),
                    )
                )
        caches[degree] = (extension, local_data)

        # Literal squares contribute +1 at every place, including their
        # double zeros; multiplying by one leaves the h character unchanged.
        square_coefficients = vector(
            extension, [100, -20, 0, 1, 0, 0]
        )
        assert character_sum(
            square_coefficients, caches[degree], 4
        ) == extended_curve.cardinality()
        known_h_sum = (
            quadratic_character(extension(1)/extension(54))
            + sum(
                0 if h_order % 2 else quadratic_character(h_leading)
                for _, _, h_order, h_leading in local_data
            )
        )
        assert product_with_known_cover_character_sum(
            square_coefficients, caches[degree], 4
        ) == known_h_sum

    candidates = []
    tested = 0
    chart_compatible = 0
    trace_one_compatible = 0
    trace_two_compatible = 0
    branch_multiples = (
        range(generator.order()) if all_rational_branches else (3, -3)
    )
    for branch_multiple in branch_multiples:
        rational_branch = branch_multiple * generator
        local_candidates = set()

        # Generic representatives in L(6O): zeros at the rational branch
        # point and a double zero Q.
        q_point = elliptic_curve(0)
        for q_index in range(generator.order()):
            if not q_point.is_zero():
                matrix_value = matrix(
                    field,
                    [
                        evaluation_row(rational_branch, field),
                        evaluation_row(q_point, field),
                        derivative_row(q_point, field),
                    ],
                )
                kernel = matrix_value.right_kernel().basis()
                for projective_coefficients in projective_vectors(
                    field, len(kernel)
                ):
                    coefficient_vector = sum(
                        (
                            projective_coefficients[index] * kernel[index]
                            for index in range(len(kernel))
                        ),
                        vector(field, 6),
                    )
                    if coefficient_vector[5]:
                        local_candidates.add(
                            (
                                normalize(coefficient_vector),
                                q_point,
                                q_index,
                                6,
                            )
                        )
            q_point += generator

        # If Q=O, cancel 2O and use L(4O).
        matrix_value = matrix(
            field, [evaluation_row(rational_branch, field)[:4]]
        )
        kernel = matrix_value.right_kernel().basis()
        for projective_coefficients in projective_vectors(field, len(kernel)):
            coefficient_vector = sum(
                (
                    projective_coefficients[index] * kernel[index]
                    for index in range(3)
                ),
                vector(field, 4),
            )
            if coefficient_vector[3]:
                local_candidates.add(
                    (
                        normalize(vector(field, list(coefficient_vector) + [0, 0])),
                        None,
                        None,
                        4,
                    )
                )

        # If one branch point is O, use an odd pole.  The generic chart is
        # L(5O), with basis 1,x,y,x^2,xy and a finite double zero Q.
        q_point = elliptic_curve(0)
        for q_index in range(generator.order()):
            if not q_point.is_zero():
                matrix_value = matrix(
                    field,
                    [
                        evaluation_row(rational_branch, field)[:5],
                        evaluation_row(q_point, field)[:5],
                        derivative_row(q_point, field)[:5],
                    ],
                )
                kernel = matrix_value.right_kernel().basis()
                for projective_coefficients in projective_vectors(
                    field, len(kernel)
                ):
                    coefficient_vector = sum(
                        (
                            projective_coefficients[index] * kernel[index]
                            for index in range(len(kernel))
                        ),
                        vector(field, 5),
                    )
                    if coefficient_vector[4]:
                        local_candidates.add(
                            (
                                normalize(
                                    vector(field, list(coefficient_vector) + [0])
                                ),
                                q_point,
                                q_index,
                                5,
                            )
                        )
            q_point += generator

        # If the double zero is also O, cancel 2O and use L(3O).
        matrix_value = matrix(
            field, [evaluation_row(rational_branch, field)[:3]]
        )
        kernel = matrix_value.right_kernel().basis()
        for projective_coefficients in projective_vectors(field, len(kernel)):
            coefficient_vector = sum(
                (
                    projective_coefficients[index] * kernel[index]
                    for index in range(2)
                ),
                vector(field, 3),
            )
            if coefficient_vector[2]:
                local_candidates.add(
                    (
                        normalize(
                            vector(field, list(coefficient_vector) + [0, 0, 0])
                        ),
                        None,
                        None,
                        3,
                    )
                )

        for coefficients, q_point, q_index, pole_degree in local_candidates:
            tested += 1
            branch_splitting = residual_branch_splitting(
                coefficients, field, rational_branch, q_point, pole_degree
            )
            if branch_splitting is None:
                continue
            expected_splitting = list(cm_field_splitting)
            if pole_degree % 2:
                if 1 not in expected_splitting:
                    if require_cm_field_splitting:
                        continue
                else:
                    expected_splitting.remove(1)
            if (
                require_cm_field_splitting
                and branch_splitting != tuple(expected_splitting)
            ):
                continue
            chart_compatible += 1
            # Filter in increasing extension degree.  The degree-two trace is
            # unchanged by a quadratic twist, while the odd-degree traces pick
            # up the twist character.  Computing over F_(p^3) only for classes
            # that already pass degrees one and two makes larger-prime exact
            # replays practical without changing the enumerated set.
            trace_one = character_sum(coefficients, caches[1], pole_degree)
            possible_twists = [
                twist_character
                for twist_character in (1, -1)
                if trace_one == twist_character * target[0]
            ]
            if not possible_twists:
                continue
            trace_one_compatible += 1
            trace_two = character_sum(coefficients, caches[2], pole_degree)
            if trace_two != target[1]:
                continue
            trace_two_compatible += 1
            trace_three = character_sum(coefficients, caches[3], pole_degree)
            for twist_character in possible_twists:
                if trace_three == twist_character * target[2]:
                    if require_v4_compatibility:
                        product_traces = [
                            product_with_known_cover_character_sum(
                                coefficients, caches[degree], pole_degree
                            )
                            for degree in (1, 2, 3)
                        ]
                        expected_product_traces = [
                            twist_character * complementary_target[0],
                            complementary_target[1],
                            twist_character * complementary_target[2],
                        ]
                        if product_traces != expected_product_traces:
                            continue
                    candidates.append(
                        (
                            coefficients,
                            q_index,
                            twist_character,
                            branch_multiple,
                            pole_degree,
                            branch_splitting,
                        )
                    )

    possible_values = {}
    for multiple in (1, -1, 4, -4, 10, -10):
        point = multiple * generator
        expansion_rows = local_basis_expansion(point, field)
        values = set()
        for coefficients, _, twist_character, _, _, _ in candidates:
            local_value = local_character(
                coefficients, expansion_rows
            )
            values.add(
                0 if local_value == 0 else twist_character*local_value
            )
        possible_values[str(multiple)] = [int(value) for value in sorted(values)]

    result = {
        "prime": int(p),
        "cm_field_splitting": [int(value) for value in cm_field_splitting],
        "cm_field_splitting_required": bool(require_cm_field_splitting),
        "target_power_traces": [int(value) for value in target],
        "tested_projective_classes": int(tested),
        "chart_compatible_classes": int(chart_compatible),
        "trace_one_compatible_classes": int(trace_one_compatible),
        "trace_two_compatible_classes": int(trace_two_compatible),
        "surviving_twisted_classes": int(len(candidates)),
        "possible_fibre_characters": possible_values,
        "decisive_local_obstructions": [
            multiple
            for multiple, values in possible_values.items()
            if values == [-1]
        ],
        "boundary": (
            "empty candidate sets and sets containing 0 or +1 are not exclusions"
        ),
    }
    if include_candidates:
        result["surviving_candidates"] = [
            {
                "coefficients": [int(value) for value in coefficients],
                "double_zero_index": (
                    None if q_index is None else int(q_index)
                ),
                "twist_character": int(twist_character),
                "rational_branch_multiple": int(branch_multiple),
                "pole_degree": int(pole_degree),
                "branch_splitting": [int(value) for value in branch_splitting],
            }
            for (
                coefficients,
                q_index,
                twist_character,
                branch_multiple,
                pole_degree,
                branch_splitting,
            ) in sorted(
                candidates,
                key=lambda row: (
                    row[4],
                    row[3],
                    -1 if row[1] is None else row[1],
                    row[2],
                    tuple(row[0]),
                ),
            )
        ]
    return result


def main():
    arguments = list(sys.argv[1:])
    include_candidates = "details" in arguments
    require_cm_field_splitting = "cmfield" in arguments
    require_v4_compatibility = "v4" in arguments
    all_rational_branches = "allbranch" in arguments
    requested_pairs = [value for value in arguments if value in PRYM_FACTORS]
    if len(requested_pairs) > 1:
        raise ValueError("choose at most one Prym pair")
    prym_pair = requested_pairs[0] if requested_pairs else "cd"
    arguments = [
        value
        for value in arguments
        if value not in ("details", "cmfield", "v4", "allbranch")
        and value not in PRYM_FACTORS
    ]
    primes = [ZZ(value) for value in arguments] or [5, 7, 11, 13]
    bad = [p for p in primes if not p.is_prime() or p in (2, 3, 103)]
    if bad:
        raise ValueError("not good primes: %s" % bad)
    payload = {
        "status": "EXPLORATORY_LOCAL_SCREEN_ONLY",
        "cover": {
            "ab": "provisional D=C_1236/<w_6> -> E=618f1",
            "cd": "provisional D=C_1236/<w_3> -> E=618f1",
        }.get(
            prym_pair,
            "diagnostic hypothetical Prym pair %s over E=618f1" % prym_pair,
        ),
        "prym_factors": [
            curve.cremona_label() for curve in PRYM_FACTORS[prym_pair]
        ],
        "rational_branch_candidates": ["+3G", "-3G"],
        "v4_compatibility_required": bool(require_v4_compatibility),
        "all_rational_branches_tested": bool(all_rational_branches),
        "nonfixed_images": ["+G", "-G", "+4G", "-4G", "+10G", "-10G"],
        "results": [
            enumerate_prime(
                p,
                include_candidates=include_candidates,
                require_cm_field_splitting=require_cm_field_splitting,
                prym_pair=prym_pair,
                require_v4_compatibility=require_v4_compatibility,
                all_rational_branches=all_rational_branches,
            )
            for p in primes
        ],
    }
    payload["decisive"] = any(
        row["decisive_local_obstructions"] for row in payload["results"]
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


main()
