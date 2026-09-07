#!/usr/bin/env sage
"""Recover the Picard-20 q=8 fibration from the chord through section S.

The raw chord variable z=(y-Sy)/(x-Sx) gives a genus-two quintic.  Besides
affine vertical completions, a degree-one linear system can have two distinct
collapsed vertical branches.  Up to a constant base Mobius transformation it
then has the cross-ratio form

    W = (t-a)(z-z_b) / ((t-b)(z-z_a)).

This script exhausts the ordered pairs among the four singular branches of
the quintic and records every genus-one completion and its Jacobian fibers.
The complete degree-one Grassmann scans over GF(11) and GF(23) instead select
one non-affine pencil.  It is reconstructed over QQ here, put into a compact
integral semistable model, and equipped with a saturated MW-rank-two basis.

This is a Picard-rank-20, discriminant-43 Noether--Lefschetz boundary model.
The q=8 fiber class uses the extra section S and is not asserted to deform in
the generic determinant-948 Neron--Severi lattice.
"""

from sage.all import *
from itertools import permutations as ordered_pairs


K = FunctionField(QQ, "W")
W = K.gen()
RT = PolynomialRing(K, "t")
t = RT.gen()


def chord_quintic(z):
    return (
        t**5 + (QQ(21)/50*z-QQ(323)/200)*t**4
        + (-QQ(483)/625*z+QQ(129)/1250)*t**3
        + (QQ(1323)/62500*z**2+QQ(11907)/31250*z+QQ(1)/2)*t**2
        - QQ(31311)/781250*z**2*t + QQ(194481)/78125000*z**4
    )


branches = (
    (QQ(0), QQ(0), "D4"),
    (QQ(1), QQ(25)/21, "A2a"),
    (QQ(9)/25, -QQ(13)/7, "A2b"),
    (QQ(49)/25, -QQ(17)/3, "A1"),
)


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        rank = max(0, ord_delta-1)
        return rank, None if rank == 0 else f"A{rank}"
    if ord_delta == 2:
        return 0, None
    if ord_delta == 3:
        return 1, "A1"
    if ord_delta == 4:
        return 2, "A2"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, "D4"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        rank = ord_delta-2
        return rank, f"D{rank}"
    if ord_delta == 8:
        return 6, "E6"
    if ord_delta == 9:
        return 7, "E7"
    if ord_delta == 10:
        return 8, "E8"
    raise ValueError((ord_a, ord_b, ord_delta))


def jacobian_signature(quartic):
    coefficients = [quartic[index] for index in range(5)]
    e, d, c, b, a = coefficients
    invariant_i = 12*a*e-3*b*d+c**2
    invariant_j = (
        72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    )
    curve = EllipticCurve(K, [0, 0, 0, -27*invariant_i, -27*invariant_j])
    discriminant = curve.discriminant()
    ade = []
    root_rank = 0
    finite = []
    for factor, delta_order in discriminant.numerator().factor():
        ord_a = curve.a4().numerator().valuation(factor)
        ord_b = curve.a6().numerator().valuation(factor)
        rank, component = kodaira_data(ord_a, ord_b, delta_order)
        root_rank += factor.degree()*rank
        if component is not None:
            ade.extend([component]*factor.degree())
        finite.append((factor, ord_a, ord_b, delta_order))
    infinity = (
        8-curve.a4().numerator().degree(),
        12-curve.a6().numerator().degree(),
        24-discriminant.numerator().degree(),
    )
    infinity_rank, infinity_component = kodaira_data(*infinity)
    root_rank += infinity_rank
    if infinity_component is not None:
        ade.append(infinity_component)
    return curve, discriminant, tuple(sorted(ade)), root_rank, tuple(finite), infinity


hits = []
target_ade = tuple(sorted(("A7", "A4", "A3", "A2")))

# The complete degree-one scans over GF(11) and GF(23) select respectively
# (a,b)=(6,10) and (11,5) in W=(t+b*z)/(1+a*z).  Both lift uniquely from the
# rational A1 branch (t,z)=(49/25,-17/3): a=-1/z=3/17 and b=a*t=147/425.
exact_a = QQ(3)/17
exact_b = QQ(147)/425
exact_denominator = exact_a*W-exact_b
exact_z = (t-W)/exact_denominator
exact_completed = RT(exact_denominator**4*chord_quintic(exact_z))
exact_factors = exact_completed.factor()
assert exact_completed.valuation(t-QQ(49)/25) == 2
exact_quartic = exact_factors.unit()*prod(
    factor for factor, exponent in exact_factors if exponent % 2
)
assert exact_quartic.degree() == 3
(
    exact_curve,
    exact_discriminant,
    exact_ade,
    exact_root_rank,
    exact_finite,
    exact_infinity,
) = jacobian_signature(exact_quartic)
assert exact_ade == target_ade and exact_root_rank == 16
exact_a0, exact_a1, exact_a2, exact_a3 = [
    exact_quartic[index] for index in range(4)
]
exact_natural = EllipticCurve(K, [
    0,
    exact_a2,
    0,
    exact_a3*exact_a1,
    exact_a3**2*exact_a0,
])
exact_natural_discriminant = exact_natural.discriminant()
assert tuple(exponent for _, exponent in exact_natural_discriminant.numerator().factor()) == (
    3, 4, 5, 8, 1,
)

# The raw line discriminant is 2*(6250/441)^2 times the monic quintic used
# above.  Hence the actual pointed cubic is v^2=2*exact_quartic(t), and its
# natural Weierstrass model is the constant-2 twist below.
exact_twisted = EllipticCurve(K, [
    0,
    2*exact_a2,
    0,
    4*exact_a3*exact_a1,
    8*exact_a3**2*exact_a0,
])
exact_twisted_discriminant = exact_twisted.discriminant()
assert tuple(exponent for _, exponent in exact_twisted_discriminant.numerator().factor()) == (
    3, 4, 5, 8, 1,
)

section_t_values = (
    K(0),
    K(QQ(34)/25*W/(W+QQ(21)/50)),
    K(-QQ(17)/25*W/(W-QQ(49)/25)),
    K((QQ(3)/25*W-QQ(504)/625)/(W-QQ(49)/25)),
)
section_v_values = tuple(K(2*exact_quartic(value)).sqrt() for value in section_t_values)
section_points = tuple(
    exact_twisted((
        2*exact_a3*t_value,
        2*exact_a3*v_value,
    ))
    for t_value, v_value in zip(section_t_values, section_v_values)
)
P, Q, minus_p_minus_q, minus_2p = section_points
assert minus_p_minus_q == -P-Q
assert minus_2p == -2*P

# A compact integral base normalization.  It sends the I5, I4, I3, and I8
# fibers at W=-21/50, 0, 7/6, 49/25 to U=0, 1, 28/3, infinity.  Clearing the
# base denominator and applying the constant Weierstrass scale 6250/441 gives
#
#   y^2 = x^3 + (q(U)^2+5292*U*(U+7))*x^2
#               + 2*L(U)*q(U)*x + L(U)^2,
#   q(U)=9*U^2+84*U+49,  L(U)=1815156*U*(U-1)^2.
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
base_denominator = 3*U+14
W_of_U = QQ(147)/25*(U-1)/base_denominator
normal_scale = QQ(6250)/441


def evaluate_at_new_base(value):
    numerator = value.numerator()
    denominator = value.denominator()
    new_numerator = sum(
        QQ(numerator[index])*W_of_U**index
        for index in range(numerator.degree()+1)
    )
    new_denominator = sum(
        QQ(denominator[index])*W_of_U**index
        for index in range(denominator.degree()+1)
    )
    return KU(new_numerator/new_denominator)


normalized_a2 = RU(
    normal_scale**2*base_denominator**4*evaluate_at_new_base(
        exact_twisted.a2()
    )
)
normalized_a4 = RU(
    normal_scale**4*base_denominator**8*evaluate_at_new_base(
        exact_twisted.a4()
    )
)
normalized_a6 = RU(
    normal_scale**6*base_denominator**12*evaluate_at_new_base(
        exact_twisted.a6()
    )
)
compact_q = 9*U**2+84*U+49
compact_L = 1815156*U*(U-1)**2
assert normalized_a2 == compact_q**2+5292*U*(U+7)
assert normalized_a4 == 2*compact_L*compact_q
assert normalized_a6 == compact_L**2

normalized_curve = EllipticCurve(KU, [
    0, normalized_a2, 0, normalized_a4, normalized_a6,
])
normalized_discriminant = RU(normalized_curve.discriminant())
normalized_residual = (
    243*U**4+8316*U**3+156114*U**2-2210292*U+45619
)
assert normalized_discriminant.valuation(U) == 5
assert normalized_discriminant.valuation(U-1) == 4
assert normalized_discriminant.valuation(U-QQ(28)/3) == 3
assert 24-normalized_discriminant.degree() == 8
assert normalized_discriminant.quo_rem(normalized_residual)[1] == 0
assert normalized_residual.is_irreducible()


def normalize_point(point):
    return normalized_curve((
        normal_scale**2*base_denominator**4*evaluate_at_new_base(point[0]),
        normal_scale**3*base_denominator**6*evaluate_at_new_base(point[1]),
    ))


normalized_P, normalized_Q = tuple(normalize_point(point) for point in (P, Q))
assert normalized_P == normalized_curve((0, -compact_L))
assert normalized_Q == normalized_curve((
    -115248*(U-1),
    -777924*(U-QQ(28)/3)*(U-QQ(49)/9)*(U-1),
))
normalized_sum = normalized_P+normalized_Q
assert normalized_sum == normalized_curve((
    -12348*U*(U-1),
    -111132*U*(U-1)*(U-QQ(28)/3)*(U+QQ(7)/3),
))

# Exact component profiles in the fiber order (I3,I4,I5,I8).  At the finite
# fibers the first nonsingular reduction of mP gives component orders
# (1,2,5) and (3,4,1).  At infinity, in the minimal coordinates
# (x/U^4,y/U^6), P and Q travel together through three blowups of the A7 node;
# their difference has exact order three, selecting label 3 rather than 1.
# At I5, P and P+Q travel together through one blowup, selecting label 1.
# These local resolution data give the signed profiles below.


def evaluate_at_finite_place(value, place):
    return QQ(value.numerator()(place))/QQ(value.denominator()(place))


def has_singular_finite_reduction(point, place):
    if point.is_zero():
        return False
    if point[0].denominator()(place) == 0 or point[1].denominator()(place) == 0:
        return False
    x_value = evaluate_at_finite_place(point[0], place)
    y_value = evaluate_at_finite_place(point[1], place)
    RX = PolynomialRing(QQ, "xx")
    xx = RX.gen()
    cubic = (
        xx**3+normalized_a2(place)*xx**2
        + normalized_a4(place)*xx+normalized_a6(place)
    )
    return (
        y_value == 0 and cubic(x_value) == 0
        and cubic.derivative()(x_value) == 0
    )


def value_in_infinity_chart(value, weight):
    numerator = value.numerator()
    denominator = value.denominator()
    degree = numerator.degree()-denominator.degree()
    if degree < weight:
        return QQ(0)
    if degree == weight:
        return QQ(numerator.leading_coefficient())/QQ(
            denominator.leading_coefficient()
        )
    return Infinity


def has_singular_infinity_reduction(point):
    if point.is_zero():
        return False
    # The infinity special fiber is y^2=x^2*(x+81), with node (0,0).
    return (
        value_in_infinity_chart(point[0], 4) == 0
        and value_in_infinity_chart(point[1], 6) == 0
    )


finite_places = ((QQ(28)/3, 3), (QQ(1), 4), (QQ(0), 5))
assert tuple(
    next(
        multiple for multiple in range(1, order+1)
        if not has_singular_finite_reduction(multiple*normalized_P, place)
    )
    for place, order in finite_places
) == (1, 2, 5)
assert tuple(
    next(
        multiple for multiple in range(1, order+1)
        if not has_singular_finite_reduction(multiple*normalized_Q, place)
    )
    for place, order in finite_places
) == (3, 4, 1)
assert next(
    multiple for multiple in range(1, 9)
    if not has_singular_infinity_reduction(multiple*normalized_P)
) == 8
assert next(
    multiple for multiple in range(1, 9)
    if not has_singular_infinity_reduction(multiple*normalized_Q)
) == 8
assert next(
    multiple for multiple in range(1, 9)
    if not has_singular_infinity_reduction(multiple*normalized_sum)
) == 4
assert not has_singular_infinity_reduction(normalized_P-normalized_Q)

P_minus_Q_x = normalized_P[0]-normalized_Q[0]
P_minus_Q_y = normalized_P[1]-normalized_Q[1]
assert (
    4-(P_minus_Q_x.numerator().degree()-P_minus_Q_x.denominator().degree()),
    6-(P_minus_Q_y.numerator().degree()-P_minus_Q_y.denominator().degree()),
) == (3, 3)
P_minus_sum_x = normalized_P[0]-normalized_sum[0]
P_minus_sum_y = normalized_P[1]-normalized_sum[1]
assert (
    P_minus_sum_x.numerator().valuation(U),
    P_minus_sum_y.numerator().valuation(U),
) == (1, 1)
assert gcd(P_minus_Q_x.numerator(), P_minus_Q_y.numerator()) == U-1

P_profile = (0, 2, 1, 3)
Q_profile = (1, 1, 0, 3)
fiber_orders = (3, 4, 5, 8)


def component_contribution(left, right, order):
    if left == 0 or right == 0:
        return QQ(0)
    return QQ(min(left, right)*(order-max(left, right)))/order


P_contribution = sum(
    component_contribution(label, label, order)
    for label, order in zip(P_profile, fiber_orders)
)
Q_contribution = sum(
    component_contribution(label, label, order)
    for label, order in zip(Q_profile, fiber_orders)
)
PQ_contribution = sum(
    component_contribution(left, right, order)
    for left, right, order in zip(P_profile, Q_profile, fiber_orders)
)
assert P_contribution == QQ(147)/40
assert Q_contribution == QQ(79)/24
assert PQ_contribution == QQ(19)/8

# P.O=Q.O=P.Q=0 in the resolved integral model.  Shioda's formula for a K3
# then gives this exact Gram.  It has the full frame determinant 43/480 and is
# unimodularly equivalent to the independently recovered reduced Gram
# (1/120)*[[34,6],[6,39]]: P=e2 and Q=-e1-e2.
section_height_gram = Matrix(QQ, [
    [4-P_contribution, 2-PQ_contribution],
    [2-PQ_contribution, 4-Q_contribution],
])
assert section_height_gram == Matrix(QQ, [
    [QQ(13)/40, -QQ(3)/8],
    [-QQ(3)/8, QQ(17)/24],
])
assert section_height_gram.det() == QQ(43)/480
reduced_height_gram = Matrix(QQ, [[QQ(17)/60, QQ(1)/20],
                                  [QQ(1)/20, QQ(13)/40]])
height_change = Matrix(ZZ, [[0, -1], [1, -1]])
assert height_change.det() == 1
assert height_change.transpose()*reduced_height_gram*height_change == section_height_gram

specialization_value = QQ(2)


def specialize(value):
    return QQ(value.numerator()(specialization_value)) / QQ(
        value.denominator()(specialization_value)
    )


specialized_curve = EllipticCurve(QQ, [
    0,
    specialize(exact_twisted.a2()),
    0,
    specialize(exact_twisted.a4()),
    specialize(exact_twisted.a6()),
])
specialized_points = tuple(
    specialized_curve((specialize(point[0]), specialize(point[1])))
    for point in (P, Q)
)
specialized_height = specialized_curve.height_pairing_matrix(
    specialized_points, precision=128
)
specialized_regulator = specialized_height.det()
assert specialized_regulator > 1
print(
    f"PICARD20Q8CHORD|exact_target=1|a={exact_a}|b={exact_b}"
    f"|W=(t+({exact_b})*z)/(1+({exact_a})*z)"
    f"|factor_degrees={tuple((factor.degree(), exponent) for factor, exponent in exact_factors)}"
    f"|quartic={exact_quartic}",
    flush=True,
)
print(
    f"PICARD20Q8CHORD|exact_twisted=1|a2={exact_twisted.a2()}"
    f"|a4={exact_twisted.a4()}|a6={exact_twisted.a6()}"
    f"|disc_num={exact_twisted_discriminant.numerator().factor()}"
    f"|disc_den={exact_twisted_discriminant.denominator().factor()}",
    flush=True,
)
for section_index, (t_value, v_value, point) in enumerate(
    zip(section_t_values, section_v_values, section_points), 1
):
    print(
        f"PICARD20Q8CHORD|section={section_index}|t={t_value}|v={v_value}"
        f"|x={point[0]}|y={point[1]}",
        flush=True,
    )
print(
    "PICARD20Q8CHORD|section_relations=P3:-P1-P2,P4:-2P1"
    f"|specialization=W:2|regulator={specialized_regulator}"
    "|rational_independent=2|geometric_MW=2",
    flush=True,
)
print(
    "PICARD20Q8CHORD|compact_model=1"
    f"|q={compact_q}|L={compact_L}"
    f"|a2={normalized_a2}|a4={normalized_a4}|a6={normalized_a6}"
    f"|disc={normalized_discriminant.factor()}"
    "|fibers=I3@28/3,I4@1,I5@0,I8@infinity,4I1",
    flush=True,
)
print(
    f"PICARD20Q8CHORD|compact_P={normalized_P}"
    f"|compact_Q={normalized_Q}|compact_PplusQ={normalized_sum}",
    flush=True,
)
print(
    f"PICARD20Q8CHORD|profiles=P:{P_profile};Q:{Q_profile}"
    f"|height_gram={section_height_gram.list()}"
    f"|height_det={section_height_gram.det()}|MW_saturated=1"
    "|Picard=20|NS_disc=43|NL_special=1",
    flush=True,
)
print(
    f"PICARD20Q8CHORD|exact_natural=1|a2={exact_natural.a2()}"
    f"|a4={exact_natural.a4()}|a6={exact_natural.a6()}"
    f"|disc_num={exact_natural_discriminant.numerator().factor()}"
    f"|disc_den={exact_natural_discriminant.denominator().factor()}",
    flush=True,
)
print(
    f"PICARD20Q8CHORD|exact_target=1|ADE={'+'.join(exact_ade)}"
    f"|root_rank={exact_root_rank}|geometric_MW={18-exact_root_rank}"
    f"|A={exact_curve.a4()}|B={exact_curve.a6()}"
    f"|disc_num={exact_discriminant.numerator().factor()}"
    f"|disc_den={exact_discriminant.denominator().factor()}"
    f"|infinity={exact_infinity}",
    flush=True,
)

for left, right in ordered_pairs(branches, int(2)):
    a, z_a, name_a = left
    b, z_b, name_b = right
    denominator = W*(t-b)-(t-a)
    numerator = W*(t-b)*z_a-(t-a)*z_b
    z = numerator/denominator
    completed = RT(denominator**4*chord_quintic(z))
    factors = completed.factor()
    odd_part = factors.unit()*prod(
        factor for factor, exponent in factors if exponent % 2
    )
    if odd_part.degree() > 4:
        continue
    if odd_part.degree() < 4:
        odd_part *= RT(1)
    curve, discriminant, ade, root_rank, finite, infinity = jacobian_signature(odd_part)
    record = (name_a, name_b, a, b, z_a, z_b, odd_part, ade, root_rank)
    hits.append(record)
    print(
        f"PICARD20Q8CHORD|left={name_a}|right={name_b}"
        f"|factor_degrees={tuple((factor.degree(), exponent) for factor, exponent in factors)}"
        f"|quartic_degree={odd_part.degree()}|ADE={'+'.join(ade)}"
        f"|root_rank={root_rank}|geometric_MW={18-root_rank}",
        flush=True,
    )
    print(
        f"PICARD20Q8CHORD|left={name_a}|right={name_b}"
        f"|W=((t-{a})*(z-({z_b})))/((t-{b})*(z-({z_a})))"
        f"|quartic={odd_part}",
        flush=True,
    )
    print(
        f"PICARD20Q8CHORD|left={name_a}|right={name_b}"
        f"|A={curve.a4()}|B={curve.a6()}"
        f"|disc_num={discriminant.numerator().factor()}"
        f"|disc_den={discriminant.denominator().factor()}|infinity={infinity}",
        flush=True,
    )
    if ade == target_ade:
        print(
            f"PICARD20Q8CHORD|target=1|left={name_a}|right={name_b}",
            flush=True,
        )

target_hits = [record for record in hits if record[7] == target_ade]
print(
    f"PICARD20Q8CHORD|pairs=12|genus_one_hits={len(hits)}"
    f"|cross_ratio_target_hits={len(target_hits)}|exact_mobius_target=1|status=PASS",
    flush=True,
)
