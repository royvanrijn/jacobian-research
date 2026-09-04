#!/usr/bin/env sage
"""Certify the exact determinant-1236 candidate double cover.

This script verifies the characteristic-zero algebra behind the candidate

    D: z^2 = 2*b0(X,Y)

over E=618f1, where

    b0 = -X^3-328*X^2-2772*X+66512+(32*X+1600)*Y.

It proves the divisor description, pulls the squareclass back to the
published genus-two quotient B, and evaluates its normalization fibres at
all twelve non-fixed rational points of B.  It also compares the candidate's
finite-field character sums through extension degree three with the asserted
618c1*618d1 and complementary 618a1*618b1 Prym factors.

Independently, the Gonzalez--Rotger residue-field formula and the complete
class-number-one/two lists prove that the actual marked curve has exactly ten
rational CM points: two of discriminant -3 and four each of discriminants
-43 and -67.  Consequently, if the displayed candidate is the Shimura cover,
its eight non-fixed rational lifts are all CM; they do not pass the rank-19
period gate.  The same count rules out the no-nonfixed-lift twist as the
actual descent once the candidate branch divisor is identified.

The finite-field matches do not identify this algebraic candidate with the
Shimura-theoretic cover C_1236 -> B.  The output therefore remains explicitly
fail-closed until a characteristic-zero modular/CM calculation proves that
the cubic branch orbit used here is the order-discriminant -1236 orbit on the
specified marked quotient.

Replay:
    sage elkies-k3/scripts/certify_det1236_candidate_double_cover.sage
    sage elkies-k3/scripts/certify_det1236_candidate_double_cover.sage check
"""

import argparse
import hashlib
import json
from pathlib import Path


# Sage preparses ``.sage`` files through a temporary directory, so ``__file__``
# is not a stable repository anchor here.  Replays are intentionally run from
# the repository root, as documented below.
ROOT = Path.cwd().resolve()
if not (ROOT / "elkies-k3").is_dir():
    raise RuntimeError("run this certificate from the repository root")
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-det1236-candidate-double-cover-v1.json"
)
TRACE_PRIMES = (5, 11, 13, 17, 19, 23, 29, 31)
LOCAL_PRECISION = 7


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quadratic_character(value):
    if not value:
        return 0
    return 1 if value.is_square() else -1


def power_trace(ap_value, prime, degree):
    if degree == 0:
        return ZZ(2)
    if degree == 1:
        return ZZ(ap_value)
    previous, current = ZZ(2), ZZ(ap_value)
    for _ in range(2, degree + 1):
        previous, current = current, ap_value*current-prime*previous
    return ZZ(current)


def prym_target(curve_labels, prime, degree):
    return -sum(
        power_trace(EllipticCurve(label).ap(prime), prime, degree)
        for label in curve_labels
    )


def evaluation(coefficients, point):
    x_value, y_value = point[0], point[1]
    basis = [1, x_value, y_value, x_value**2, x_value*y_value, x_value**3]
    return sum(coefficients[index]*basis[index] for index in range(6))


def local_basis_expansion(point, field, precision=LOCAL_PRECISION):
    """Expand 1,X,Y,X^2,XY,X^3 in a local uniformizer."""
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


def local_order_and_leading(coefficients, point, field):
    if evaluation(coefficients, point):
        return 0, evaluation(coefficients, point)
    rows = local_basis_expansion(point, field)
    for order, row in enumerate(rows):
        leading = sum(
            coefficients[index]*row[index] for index in range(6)
        )
        if leading:
            return order, leading
    raise ArithmeticError("local vanishing order exceeds L(6O) precision")


def character_sums(prime, degree, candidate_coefficients):
    field = GF(prime**degree, "z")
    curve = EllipticCurve(field, [1, 0, 0, -185, 1401])
    coefficients = [field(value) for value in candidate_coefficients]
    h_coefficients = [field(-4), field(1), 0, 0, 0, 0]
    # At O, b has pole order 6 with leading coefficient coefficients[5].
    # The known h=(X-4)/54 has pole order 2 with leading coefficient 1/54.
    sum_b = quadratic_character(coefficients[5])
    sum_product = quadratic_character(coefficients[5]/field(54))
    for point in curve:
        if point.is_zero():
            continue
        b_order, b_leading = local_order_and_leading(
            coefficients, point, field
        )
        h_order, h_leading = local_order_and_leading(
            h_coefficients, point, field
        )
        h_leading /= field(54)
        if b_order % 2 == 0:
            sum_b += quadratic_character(b_leading)
        if (b_order+h_order) % 2 == 0:
            sum_product += quadratic_character(b_leading*h_leading)

    # Regression: a literal square contributes +1 at every point of E,
    # including its double zeros, and does not alter h.
    square = [field(100), field(-20), 0, field(1), 0, 0]
    square_sum = quadratic_character(square[3])
    square_product_sum = quadratic_character(square[3]/field(54))
    for point in curve:
        if point.is_zero():
            continue
        square_order, square_leading = local_order_and_leading(
            square, point, field
        )
        h_order, h_leading = local_order_and_leading(
            h_coefficients, point, field
        )
        h_leading /= field(54)
        if square_order % 2 == 0:
            square_sum += quadratic_character(square_leading)
        if (square_order+h_order) % 2 == 0:
            square_product_sum += quadratic_character(
                square_leading*h_leading
            )
    h_sum = quadratic_character(field(1)/field(54))
    for point in curve:
        if point.is_zero():
            continue
        h_order, h_leading = local_order_and_leading(
            h_coefficients, point, field
        )
        if h_order % 2 == 0:
            h_sum += quadratic_character(h_leading/field(54))
    assert square_sum == curve.cardinality()
    assert square_product_sum == h_sum
    return ZZ(sum_b), ZZ(sum_product)


def rational_squareclass(value):
    value = QQ(value)
    if not value:
        return {"value": "0", "squareclass": "ramified_or_even_zero"}
    numerator = ZZ(value.numerator())
    denominator = ZZ(value.denominator())
    sign = -1 if numerator < 0 else 1
    squarefree = ZZ(sign)
    for prime, exponent in factor(abs(numerator*denominator)):
        if exponent % 2:
            squarefree *= prime
    return {
        "value": str(value),
        "squareclass": str(squarefree),
        "is_square": bool(value.is_square()),
    }


def cm_order_record(class_number, conductor, field_discriminant):
    """Apply the optimal-embedding and w_618 residue-degree criteria."""
    conductor = ZZ(conductor)
    field_discriminant = ZZ(field_discriminant)
    order_discriminant = conductor**2*field_discriminant

    def eichler_symbol(prime):
        return ZZ(1) if conductor % prime == 0 else ZZ(
            kronecker(field_discriminant, prime)
        )

    local_factors = {
        "2": int(1-eichler_symbol(2)),
        "3": int(1-eichler_symbol(3)),
        "103": int(1+eichler_symbol(103)),
    }
    top_cm_points = ZZ(class_number)*prod(local_factors.values())
    d_r = prod(
        prime for prime in (2, 3)
        if conductor % prime and kronecker(field_discriminant, prime) == -1
    )
    n_r = prod(
        prime for prime in (103,)
        if conductor % prime == 0
        or kronecker(field_discriminant, prime) == 1
    )
    nstar_r = prod(
        prime for prime in (103,)
        if conductor % prime
        and kronecker(field_discriminant, prime) == 1
    )
    m_r = gcd(
        ZZ(618), abs(order_discriminant)//gcd(ZZ(103), conductor)
    )
    quotient = ZZ(618)//m_r
    d_r_nstar_r = ZZ(d_r*nstar_r)

    # Corollary 5.14 fixes H_R by one involution in every case except
    # D(R)N*(R)=1 and m/m_r=1, where two involutions occur.  Since the
    # class of the norm-m_r ideal has order at most two, a rational fixed
    # field forces h(R)<=2.  For h=1 the following condition is exact.  For
    # h=2 the displayed condition is necessary; no class-number-two order
    # with a nonempty CM locus reaches it here, so no Artin-class choice is
    # needed.
    if class_number == 1:
        rational_image = bool(
            top_cm_points
            and (
                d_r_nstar_r == 1
                or quotient == d_r_nstar_r
            )
        )
    else:
        rational_image = bool(
            top_cm_points
            and d_r_nstar_r == 1
            and quotient == 1
        )

    return {
        "class_number": int(class_number),
        "conductor": int(conductor),
        "field_discriminant": int(field_discriminant),
        "order_discriminant": int(order_discriminant),
        "local_embedding_factors_p2_p3_p103": local_factors,
        "top_curve_cm_points": int(top_cm_points),
        "D_R": int(d_r),
        "N_R": int(n_r),
        "N_star_R": int(nstar_r),
        "m_R": int(m_r),
        "m_over_m_R": int(quotient),
        "rational_image_on_w618_quotient": rational_image,
    }


def build_payload():
    # Complete lists of imaginary quadratic orders of class number one and
    # two, represented as (conductor, fundamental discriminant).  These are
    # the pinned Watkins lists used by Padurariu--Saia.  Corollary 5.14 shows
    # that no order of larger class number can have a rational image on this
    # single Atkin--Lehner quotient.
    class_number_one_orders = [
        (1, -3), (2, -3), (3, -3), (1, -4), (2, -4),
        (1, -7), (2, -7), (1, -8), (1, -11), (1, -19),
        (1, -43), (1, -67), (1, -163),
    ]
    class_number_two_orders = [
        (4, -3), (5, -3), (7, -3), (3, -4), (4, -4),
        (5, -4), (4, -7), (2, -8), (3, -8), (3, -11),
        (1, -15), (2, -15), (1, -20), (1, -24), (1, -35),
        (1, -40), (1, -51), (1, -52), (1, -88), (1, -91),
        (1, -115), (1, -123), (1, -148), (1, -187),
        (1, -232), (1, -235), (1, -267), (1, -403),
        (1, -427),
    ]
    cm_order_rows = []
    for class_number, orders in (
        (1, class_number_one_orders),
        (2, class_number_two_orders),
    ):
        for conductor, field_discriminant in orders:
            assert len(BinaryQF_reduced_representatives(
                ZZ(conductor**2*field_discriminant),
                primitive_only=True,
            )) == class_number
            cm_order_rows.append(cm_order_record(
                class_number, conductor, field_discriminant
            ))
    rational_cm_rows = [
        row for row in cm_order_rows
        if row["rational_image_on_w618_quotient"]
    ]
    assert [row["order_discriminant"] for row in rational_cm_rows] == [
        -3, -43, -67
    ]
    for row in rational_cm_rows:
        # None lies in the order-discriminant -2472 fixed locus of w_618,
        # so the quotient pairs its top-curve CM points freely.
        assert row["order_discriminant"] != -2472
        row["marked_curve_rational_cm_points"] = (
            row["top_curve_cm_points"]//2
        )
    assert [
        row["marked_curve_rational_cm_points"] for row in rational_cm_rows
    ] == [2, 4, 4]
    rational_cm_point_count = sum(
        row["marked_curve_rational_cm_points"] for row in rational_cm_rows
    )
    assert rational_cm_point_count == 10

    polynomial_ring = PolynomialRing(QQ, "u")
    u_value = polynomial_ring.gen()
    cubic_field = NumberField(u_value**3-u_value**2+4*u_value+12, "a")
    a_value = cubic_field.gen()

    # The order discriminant -1236 is fundamental.  Its Hilbert class field
    # has generalized-dihedral Galois closure over QQ.  The three embedded
    # cubic subfields are conjugate and have one QQ-isomorphism class.  This
    # verifies that the candidate branch field is the exact cubic CM
    # residue-field class; it deliberately does not identify which K-point
    # on E is the Shimura CM image.
    cm_field = QuadraticField(-1236, "d")
    cm_class_group_invariants = tuple(
        ZZ(value) for value in cm_field.class_group().invariants()
    )
    assert cm_class_group_invariants == (6, 2)
    hilbert_class_field = cm_field.hilbert_class_field("h")
    assert hilbert_class_field.relative_degree() == 12
    assert hilbert_class_field.absolute_degree() == 24
    absolute_hilbert_class_field = hilbert_class_field.absolute_field("r")
    cubic_subfield_tuples = absolute_hilbert_class_field.subfields(3)
    cubic_subfields = [entry[0] for entry in cubic_subfield_tuples]
    assert len(cubic_subfields) == 3
    assert all(field.discriminant() == -1236 for field in cubic_subfields)
    assert all(field.is_isomorphic(cubic_field) for field in cubic_subfields)

    curve_k = EllipticCurve(cubic_field, [1, 0, 0, -185, 1401])
    generator_k = curve_k(10, -29)
    seed = curve_k(
        -3*a_value**2-9*a_value-20,
        15*a_value**2-57*a_value-155,
    )
    branch_point = 9*generator_k+2*seed
    expected_branch_point = curve_k(
        18*a_value**2-48*a_value+130,
        270*a_value**2-672*a_value+2059,
    )
    assert branch_point == expected_branch_point

    splitting_field = cubic_field.defining_polynomial().splitting_field("s")
    roots = cubic_field.defining_polynomial().change_ring(
        splitting_field
    ).roots(multiplicities=False)
    curve_split = EllipticCurve(
        splitting_field, [1, 0, 0, -185, 1401]
    )
    seed_trace = curve_split(0)
    branch_trace = curve_split(0)
    branch_orbit = []
    for root in roots:
        seed_conjugate = curve_split(
            -3*root**2-9*root-20,
            15*root**2-57*root-155,
        )
        branch_conjugate = 9*curve_split(10, -29)+2*seed_conjugate
        seed_trace += seed_conjugate
        branch_trace += branch_conjugate
        branch_orbit.append(branch_conjugate)
    assert seed_trace == -16*curve_split(10, -29)
    assert branch_trace == -5*curve_split(10, -29)

    curve_q = EllipticCurve(QQ, [1, 0, 0, -185, 1401])
    generator = curve_q(10, -29)
    rational_branch = -3*generator
    double_zero = 4*generator
    assert rational_branch == curve_q(4, -29)
    assert double_zero == curve_q(58, 403)
    assert branch_trace + curve_split(rational_branch) + 2*curve_split(double_zero) == curve_split(0)

    candidate_coefficients = [
        ZZ(2*66512), ZZ(-2*2772), ZZ(2*1600),
        ZZ(-2*328), ZZ(2*32), ZZ(-2),
    ]
    unscaled_coefficients = [value/2 for value in candidate_coefficients]
    assert evaluation(unscaled_coefficients, rational_branch) == 0
    assert evaluation(unscaled_coefficients, double_zero) == 0
    q_order, q_leading_unscaled = local_order_and_leading(
        unscaled_coefficients, double_zero, QQ
    )
    assert q_order == 2
    assert q_leading_unscaled == 18
    assert all(evaluation(unscaled_coefficients, point) == 0 for point in branch_orbit)

    x_ring = PolynomialRing(QQ, "X")
    X = x_ring.gen()
    right_side = X**3-185*X+1401
    a_polynomial = -X**3-328*X**2-2772*X+66512
    y_coefficient = 32*X+1600
    norm = a_polynomial**2-X*a_polynomial*y_coefficient-right_side*y_coefficient**2
    cubic_branch_polynomial = X**3-216*X**2-6924*X-62224
    expected_norm = (X-4)*(X-58)**2*cubic_branch_polynomial
    assert norm == expected_norm
    assert cubic_branch_polynomial(branch_point[0]) == 0
    assert cubic_branch_polynomial.discriminant() == -1236*ZZ(24192)**2

    xy_ring = PolynomialRing(QQ, names=("x", "y"))
    x_value, y_value = xy_ring.gens()
    X_pullback = 54*x_value**2+4
    Y_pullback = 9*y_value-27*x_value**2-2
    b0_pullback = (
        -X_pullback**3-328*X_pullback**2-2772*X_pullback+66512
        +(32*X_pullback+1600)*Y_pullback
    )
    q_polynomial = (
        -81*x_value**6-534*x_value**4+8*x_value**2*y_value
        -177*x_value**2+8*y_value+24
    )
    assert 2*b0_pullback == 3888*q_polynomial
    assert 3888 == 3*ZZ(36)**2
    b_right_side = (
        1944*x_value**6+441*x_value**4-90*x_value**2+9
    )
    q_even_part = -81*x_value**6-534*x_value**4-177*x_value**2+24
    q_y_part = 8*(x_value**2+1)
    q_norm = q_even_part**2-q_y_part**2*b_right_side
    branch_x_polynomial = (
        81*x_value**6-306*x_value**4-239*x_value**2-48
    )
    assert q_norm == (
        81*x_value**2*(x_value**2-1)**2*branch_x_polynomial
    )
    assert branch_x_polynomial.gcd(
        branch_x_polynomial.derivative(x_value)
    ) == 1
    assert cubic_branch_polynomial(54*x_value**2+4) == 1944*branch_x_polynomial
    candidate_branch_degree = branch_x_polynomial.degree()
    candidate_genus = 2*2-1+candidate_branch_degree//2
    assert candidate_branch_degree == 6
    assert candidate_genus == 6

    # The canonical plane-quartic model of the genus-three candidate uses
    # [U:V:W]=[X-58:Y-403:z].  Direct ideal reduction certifies the model;
    # this gives the next modular-identification calculation an exact target.
    quartic_ring = PolynomialRing(QQ, names=("U", "V", "W"))
    U, V, W = quartic_ring.gens()
    canonical_quartic = (
        87616*U**4-43784*U**3*V+7944*U**2*V**2-608*U*V**3
        +16*V**4+268*U**2*W**2-118*U*V*W**2
        +10*V**2*W**2+W**4
    )
    projective_plane = ProjectiveSpace(QQ, 2, names=("U", "V", "W"))
    plane_U, plane_V, plane_W = projective_plane.coordinate_ring().gens()
    plane_quartic = canonical_quartic(
        U=plane_U, V=plane_V, W=plane_W
    )
    assert Curve(plane_quartic).is_smooth()
    relation_ring = PolynomialRing(QQ, names=("XX", "YY"), order="lex")
    XX, YY = relation_ring.gens()
    relation_b0 = (
        -XX**3-328*XX**2-2772*XX+66512+(32*XX+1600)*YY
    )
    quartic_relation = canonical_quartic(
        U=XX-58, V=YY-403, W=relation_ring(1)
    )
    # Homogeneity says W=z, not W=1, before projectivizing.
    quartic_relation_with_z = canonical_quartic(
        U=XX-58, V=YY-403, W=relation_ring(1)
    )
    # Substitute W^2=2*b0 in the terms of W-degree two and four.
    quartic_after_cover = (
        87616*(XX-58)**4
        -43784*(XX-58)**3*(YY-403)
        +7944*(XX-58)**2*(YY-403)**2
        -608*(XX-58)*(YY-403)**3
        +16*(YY-403)**4
        +(268*(XX-58)**2-118*(XX-58)*(YY-403)
          +10*(YY-403)**2)*(2*relation_b0)
        +(2*relation_b0)**2
    )
    elliptic_relation = YY**2+XX*YY-XX**3+185*XX-1401
    assert relation_ring.ideal(elliptic_relation).reduce(
        quartic_after_cover
    ) == 0

    b_points = [
        (QQ(1), QQ(48)), (QQ(-1), QQ(48)),
        (QQ(1), QQ(-48)), (QQ(-1), QQ(-48)),
        (QQ(1)/3, QQ(8)/3), (QQ(-1)/3, QQ(8)/3),
        (QQ(1)/3, QQ(-8)/3), (QQ(-1)/3, QQ(-8)/3),
        (QQ(1)/5, QQ(312)/125), (QQ(-1)/5, QQ(312)/125),
        (QQ(1)/5, QQ(-312)/125), (QQ(-1)/5, QQ(-312)/125),
    ]
    point_evaluations = []
    rational_lifting_points = []
    for x_coordinate, y_coordinate in b_points:
        assert y_coordinate**2 == (
            1944*x_coordinate**6+441*x_coordinate**4
            -90*x_coordinate**2+9
        )
        e_point = curve_q(
            54*x_coordinate**2+4,
            9*y_coordinate-27*x_coordinate**2-2,
        )
        multiple = next(
            value for value in (1, -1, 4, -4, 10, -10)
            if value*generator == e_point
        )
        pullback_value = 3*q_polynomial(x=x_coordinate, y=y_coordinate)
        if pullback_value:
            square_data = rational_squareclass(pullback_value)
            lifts = bool(pullback_value.is_square())
            normalization_note = "ordinary_nonzero_fibre"
        else:
            order, leading = local_order_and_leading(
                candidate_coefficients, e_point, QQ
            )
            assert order == 2
            square_data = rational_squareclass(leading)
            lifts = bool(leading.is_square())
            normalization_note = (
                "even_zero_on_the_affine_model; normalization fibre is "
                "decided by the leading local unit"
            )
        record = {
            "B_point": [str(x_coordinate), str(y_coordinate)],
            "E_image_multiple_of_G": int(multiple),
            "affine_pullback_value_3q": str(pullback_value),
            "normalization_fibre_squareclass": square_data,
            "rational_normalization_fibre": lifts,
            "normalization_note": normalization_note,
        }
        point_evaluations.append(record)
        if lifts:
            rational_lifting_points.append(record["B_point"])
    assert rational_lifting_points == [
        ["1", "48"], ["-1", "48"], ["1/3", "8/3"], ["-1/3", "8/3"]
    ]

    alternative_twist_lifting_points = []
    for record in point_evaluations:
        representative = QQ(
            record["normalization_fibre_squareclass"]["value"]
        )
        if (-3*representative).is_square():
            alternative_twist_lifting_points.append(record["B_point"])
    assert alternative_twist_lifting_points == []

    # The two fixed fibres explain why the remaining constant twist matters.
    # For the displayed candidate (0,3) has value 144.  At (0,-3), write
    # y=-3+15*x^2+O(x^4); then 3*q=-243*x^2+O(x^4).
    fixed_fibres = [
        {
            "B_point": ["0", "3"],
            "candidate_local_representative": "144",
            "candidate_residue_field": "QQ x QQ",
            "minus_three_twist_residue_field": "QQ(sqrt(-3))",
        },
        {
            "B_point": ["0", "-3"],
            "candidate_local_representative": "-243",
            "candidate_residue_field": "QQ(sqrt(-3))",
            "minus_three_twist_residue_field": "QQ x QQ",
        },
    ]
    assert QQ(144).is_square()
    assert not QQ(-243).is_square()
    assert QQ((-3)*(-243)).is_square()

    rational_e_fibres = []
    for multiple in (1, -1, 4, -4, 10, -10, 3, -3):
        point = multiple*generator
        value = evaluation(candidate_coefficients, point)
        order, leading = local_order_and_leading(
            candidate_coefficients, point, QQ
        )
        rational_e_fibres.append(
            {
                "multiple_of_G": int(multiple),
                "point": [str(point[0]), str(point[1])],
                "value": str(value),
                "local_order": int(order),
                "leading_unit": str(leading),
                "leading_squareclass": rational_squareclass(leading),
            }
        )

    trace_rows = []
    for prime in TRACE_PRIMES:
        degree_rows = []
        for degree in (1, 2, 3):
            b_sum, product_sum = character_sums(
                prime, degree, candidate_coefficients
            )
            cd_target = prym_target(("618c1", "618d1"), prime, degree)
            ab_target = prym_target(("618a1", "618b1"), prime, degree)
            assert b_sum == cd_target
            assert product_sum == ab_target
            degree_rows.append(
                {
                    "degree": int(degree),
                    "candidate_b_character_sum": int(b_sum),
                    "target_618c1_618d1": int(cd_target),
                    "candidate_hb_character_sum": int(product_sum),
                    "target_618a1_618b1": int(ab_target),
                }
            )
        trace_rows.append({"prime": int(prime), "degrees": degree_rows})

    return {
        "schema": "elkies-k3.det1236-candidate-double-cover.v1",
        "status": "UNRESOLVED_FOR_EXPLICIT_REASON",
        "candidate": {
            "elliptic_base": "618f1: Y^2+X*Y=X^3-185*X+1401",
            "generator": "G=(10,-29)",
            "cover_squareclass": (
                "2*(-X^3-328*X^2-2772*X+66512+(32*X+1600)*Y)"
            ),
            "divisor": (
                "(-3G)+Orbit(P)+2*(4G)-6O, "
                "P=(18*a^2-48*a+130,270*a^2-672*a+2059)"
            ),
            "cubic_field": "a^3-a^2+4*a+12=0",
            "cubic_field_discriminant": int(cubic_field.discriminant()),
            "cm_residue_field_match": {
                "quadratic_cm_field_discriminant": -1236,
                "cm_class_group_invariants": [
                    int(value) for value in cm_class_group_invariants
                ],
                "hilbert_class_field_relative_degree": int(
                    hilbert_class_field.relative_degree()
                ),
                "hilbert_class_field_absolute_degree": int(
                    hilbert_class_field.absolute_degree()
                ),
                "embedded_cubic_subfield_count": len(cubic_subfields),
                "embedded_cubic_subfield_discriminants": [
                    int(field.discriminant()) for field in cubic_subfields
                ],
                "unique_cubic_QQ_isomorphism_class": True,
                "candidate_field_is_that_class": True,
                "boundary": (
                    "This identifies the exact CM residue-field class, not "
                    "the candidate point P with the Shimura CM image on E."
                ),
            },
            "seed_point": (
                "A=(-3*a^2-9*a-20,15*a^2-57*a-155)"
            ),
            "branch_point_relation": "P=9G+2A",
            "trace_relations": ["Tr(A)=-16G", "Tr(P)=-5G"],
            "norm_factorization": str(expected_norm.factor()),
            "cubic_x_polynomial": str(cubic_branch_polynomial),
            "cubic_x_discriminant": int(
                cubic_branch_polynomial.discriminant()
            ),
            "smooth_canonical_plane_quartic": str(canonical_quartic),
            "canonical_coordinates": "[U:V:W]=[X-58:Y-403:z]",
        },
        "pullback_to_B": {
            "B": "y^2=1944*x^6+441*x^4-90*x^2+9",
            "map_to_E": ["X=54*x^2+4", "Y=9*y-27*x^2-2"],
            "squareclass": (
                "3*(-81*x^6-534*x^4+8*x^2*y-177*x^2+8*y+24)"
            ),
            "exact_identity": "2*b0(X(x),Y(x,y))=3888*q=3*(36^2)*q",
            "norm_factorization": (
                "Norm_B(q)=81*x^2*(x^2-1)^2*"
                "(81*x^6-306*x^4-239*x^2-48)"
            ),
            "geometric_branch_degree": int(candidate_branch_degree),
            "normalization_genus": int(candidate_genus),
        },
        "nonfixed_rational_point_evaluations": point_evaluations,
        "conditional_rational_lifting_points": rational_lifting_points,
        "conditional_rational_point_count_above_nonfixed_points": int(
            2*len(rational_lifting_points)
        ),
        "fixed_cm_fibre_twist_fork": {
            "fibres": fixed_fibres,
            "alternative_squareclass": "-3 times the displayed candidate",
            "alternative_nonfixed_rational_lifting_points": (
                alternative_twist_lifting_points
            ),
            "interpretation": (
                "The genus-two model identifies the rational and quadratic "
                "fixed CM fibres only up to y-sign. Thus the fixed-fibre "
                "data alone leave precisely the displayed twist and its "
                "-3 twist as this reconstruction's two arithmetic outcomes."
            ),
        },
        "actual_marked_curve_rational_cm_locus": {
            "rational_order_discriminants": [-3, -43, -67],
            "rows": rational_cm_rows,
            "rational_cm_point_count": int(rational_cm_point_count),
            "completeness_argument": (
                "For C_1236=X_0^6(103)/<w_618>, Gonzalez--Rotger "
                "Corollary 5.14 makes the residue degree at least h(R), "
                "except in its two-involution case where it is at least "
                "h(R)/2. Hence a rational CM image forces h(R)<=2. "
                "The complete class-number-one/two order lists and exact "
                "local optimal-embedding factors leave precisely -3, -43, "
                "and -67, contributing 2, 4, and 4 rational marked points."
            ),
            "conditional_cover_consequence": (
                "If the displayed candidate is C_1236, its total rational "
                "point count is 2 fixed plus 8 non-fixed. This equals the "
                "complete ten-point rational CM locus, so every conditional "
                "lift is CM. The -3 twist has only the two fixed rational "
                "points and therefore cannot be the actual marked descent."
            ),
        },
        "elliptic_base_local_evaluations": rational_e_fibres,
        "finite_field_prym_trace_matches": {
            "primes": [int(prime) for prime in TRACE_PRIMES],
            "extension_degrees": [1, 2, 3],
            "rows": trace_rows,
            "interpretation": (
                "At every listed prime the candidate cover has Prym traces "
                "618c1+618d1 and its product with h=(X-4)/54 has Prym "
                "traces 618a1+618b1. These exact finite computations are "
                "strong identification evidence, not a characteristic-zero "
                "isomorphism or modular-function proof."
            ),
        },
        "arithmetic_decision": {
            "status": "UNRESOLVED_FOR_EXPLICIT_REASON",
            "precise_obstruction": (
                "The candidate cubic field is already certified as the "
                "unique QQ-isomorphism class of cubic subfields in the "
                "discriminant -1236 Hilbert class field. Prove in "
                "characteristic zero that the particular point P=9G+2A "
                "is the image of that CM orbit on the specified "
                "Atkin-Lehner quotient. The exact "
                "rational CM count then forces the displayed constant twist, "
                "but until the branch-orbit identification is certified the "
                "displayed cover cannot be promoted to C_1236."
            ),
            "conditional_consequence": (
                "If the candidate is the marked Shimura cover, the four "
                "non-fixed B-points (+/-1,48) and (+/-1/3,8/3) each have "
                "two rational normalization lifts. Together with the two "
                "fixed points these exhaust the complete ten-point rational "
                "CM locus. Thus the candidate would prove "
                "ARITHMETICALLY_EXCLUDED, not realizable."
            ),
        },
        "reproduce": (
            "sage elkies-k3/scripts/"
            "certify_det1236_candidate_double_cover.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", choices=("write", "check"), default="write")
    args = parser.parse_args()
    output = OUTPUT
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n"
    if args.mode == "check":
        if not output.is_file():
            raise FileNotFoundError(output)
        if output.read_text() != rendered:
            raise AssertionError("generated artifact changed: %s" % output)
        print(json.dumps({
            "status": "PASS_DET1236_CANDIDATE_DOUBLE_COVER_CHECK",
            "output": relative(output),
            "sha256": digest(output),
            "arithmetic_decision": payload["arithmetic_decision"]["status"],
            "conditional_lifting_points": len(
                payload["conditional_rational_lifting_points"]
            ),
        }, sort_keys=True))
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
    print(json.dumps({
        "status": "WROTE_DET1236_CANDIDATE_DOUBLE_COVER_CERTIFICATE",
        "output": relative(output),
        "sha256": digest(output),
        "arithmetic_decision": payload["arithmetic_decision"]["status"],
        "conditional_lifting_points": len(
            payload["conditional_rational_lifting_points"]
        ),
    }, sort_keys=True))


main()
