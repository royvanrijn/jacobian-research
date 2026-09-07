from sage.all import *
from itertools import combinations, product
import argparse


parser = argparse.ArgumentParser(
    description="Search CRT vertical completions of the MW2 degree-two chord map."
)
parser.add_argument("--p", type=int, default=23)
parser.add_argument("--chord", choices=("p1", "s"), default="p1")
parser.add_argument("--q-degree", type=int, default=1, choices=(0, 1, 2))
parser.add_argument(
    "--denominator-degree", type=int, default=0, choices=(0, 1, 2, 3),
    help="allow d(t) supported at old special fibers in z=(g*U+h)/d",
)
parser.add_argument(
    "--mobius-linear", action="store_true",
    help="use z=((t*h+g)*U+h)/(t*U+1), the first translated-pencil chart",
)
parser.add_argument(
    "--mobius-quadratic-middle", type=int,
    help="use c=t^2+a*t (a modulo p) in z=((c*h+g)*U+h)/(c*U+1)",
)
args = parser.parse_args()

k = GF(args.p)
assert args.p == 23, "the pinned local branch data below are for GF(23)"
K = FunctionField(k, "U")
U = K.gen()
RT = PolynomialRing(K, "t")
t = RT.gen()


def chord_discriminant(z):
    if args.chord == "p1":
        return (
            t**5 + (11*z + 20)*t**4 + (4*z + 19)*t**3
            + (19*z**2 + 16*z + 20)*t**2
            + (2*z**2 + 4*z + 2)*t
            + 13*z**4 + 14*z**2 + 11*z + 7
        )
    return (
        t**5 + (11*z + 10)*t**4 + 19*t**3
        + (9*z**2 + z + 12)*t**2 + 17*z**2*t + 13*z**4
    )


if args.chord == "p1":
    local_options = (
        (k(0), (k(3), k(22))),
        (k(1), (k(0),)),
        (k(16), (k(21), k(2))),
    )
    denominator_support = tuple(map(k, (0, 1, 16, 13, 8)))
    target_ade = tuple(sorted(("A5", "D4", "A2", "A2", "A1")))
else:
    # Reductions of (t,z)=(0,0),(1,25/21),(9/25,-13/7),
    # (49/25,-17/3), the four singular branches of the S-chord quintic.
    local_options = (
        (k(0), (k(0),)),
        (k(1), (k(22),)),
        (k(16), (k(8),)),
        (k(15), (k(2),)),
    )
    denominator_support = tuple(map(k, (0, 1, 16, 15, 8)))
    target_ade = tuple(sorted(("A7", "A4", "A3", "A2")))


def interpolate(points):
    answer = RT.zero()
    for index, (abscissa, ordinate) in enumerate(points):
        numerator = RT.one()
        denominator = k.one()
        for other, (other_abscissa, _) in enumerate(points):
            if other == index:
                continue
            numerator *= t - other_abscissa
            denominator *= abscissa - other_abscissa
        answer += K(ordinate / denominator) * numerator
    return answer


hits = []
choices = []
if args.chord == "s":
    # The infinity-only correction u=U-(50/21)t has no finite branch support.
    choices.append([])
for states in product(*[(None,) + options for _, options in local_options]):
    points = [
        (local_options[index][0], value)
        for index, value in enumerate(states) if value is not None
    ]
    if points:
        choices.append(points)

denominator_choices = []
for degree in range(args.denominator_degree + 1):
    denominator_choices.extend(combinations(denominator_support, degree))

target_hits = []
tested = 0
for denominator_points in denominator_choices:
  denominator_point_set = set(denominator_points)
  dpoly = RT(prod(t - abscissa for abscissa in denominator_points))
  for points in choices:
    # A shared zero would cancel from g,h,d and is represented by a smaller
    # denominator support, so omit it from this normalized enumeration.
    if any(abscissa in denominator_point_set for abscissa, _ in points):
        continue
    g = prod(t - abscissa for abscissa, _ in points)
    h0 = interpolate([
        (abscissa, ordinate * dpoly(abscissa))
        for abscissa, ordinate in points
    ])
    coefficient_ranges = [k] * args.q_degree
    for coefficients in product(*coefficient_ranges):
        tested += 1
        # In the affine chart a constant multiple of g only translates U.
        # In the t-dependent Mobius chart it is not a constant base change and
        # must be retained.
        first_power = (
            0 if args.mobius_linear or args.mobius_quadratic_middle is not None
            else 1
        )
        q = sum(
            K(coefficients[index]) * t**(index + first_power)
            for index in range(len(coefficients))
        )
        h = h0 + q * g
        if args.mobius_linear or args.mobius_quadratic_middle is not None:
            if dpoly != 1:
                raise ValueError("Mobius charts currently require denominator degree 0")
            if args.mobius_quadratic_middle is None:
                mobius_c = t
            else:
                mobius_c = t**2 + k(args.mobius_quadratic_middle)*t
            mobius_denominator = mobius_c*U + 1
            mobius_numerator = (mobius_c*h + g)*U + h
            rational_polynomial = (
                mobius_denominator**4
                * chord_discriminant(mobius_numerator / mobius_denominator)
            )
        else:
            rational_polynomial = dpoly**4 * chord_discriminant((g * U + h) / dpoly)
        polynomial = RT(rational_polynomial)
        factors = polynomial.factor()
        odd_part = prod(factor for factor, exponent in factors if exponent % 2)
        if odd_part.degree() <= 4:
            binary_quartic = factors.unit() * odd_part
            coefficients4 = [binary_quartic[index] for index in range(5)]
            e, d, c, b, a = coefficients4
            invariant_i = 12 * a * e - 3 * b * d + c**2
            invariant_j = (
                72 * a * c * e + 9 * b * c * d - 27 * a * d**2
                - 27 * b**2 * e - 2 * c**3
            )
            jacobian = EllipticCurve(K, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
            jacobian_discriminant = jacobian.discriminant()

            def kodaira_data(ord_a, ord_b, ord_delta):
                if ord_a == 0 or ord_b == 0:
                    rank = max(0, ord_delta - 1)
                    return rank, None if rank == 0 else f"A{rank}"  # I_n
                if ord_delta == 2:
                    return 0, None  # II
                if ord_delta == 3:
                    return 1, "A1"  # III
                if ord_delta == 4:
                    return 2, "A2"  # IV
                if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
                    return 4, "D4"  # I0*
                if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
                    rank = ord_delta - 2
                    return rank, f"D{rank}"  # I_(n)*
                if ord_delta == 8:
                    return 6, "E6"  # IV*
                if ord_delta == 9:
                    return 7, "E7"  # III*
                if ord_delta == 10:
                    return 8, "E8"  # II*
                raise ValueError((ord_a, ord_b, ord_delta))

            root_rank = 0
            ade = []
            local_valuations = []
            for delta_factor, delta_order in jacobian_discriminant.numerator().factor():
                ord_a = jacobian.a4().numerator().valuation(delta_factor)
                ord_b = jacobian.a6().numerator().valuation(delta_factor)
                local_rank, local_ade = kodaira_data(ord_a, ord_b, delta_order)
                root_rank += delta_factor.degree() * local_rank
                if local_ade is not None:
                    ade.extend([local_ade] * delta_factor.degree())
                local_valuations.append((delta_factor.degree(), ord_a, ord_b, delta_order, local_rank))
            infinity_orders = (
                8 - jacobian.a4().numerator().degree(),
                12 - jacobian.a6().numerator().degree(),
                24 - jacobian_discriminant.numerator().degree(),
            )
            infinity_rank, infinity_ade = kodaira_data(*infinity_orders)
            root_rank += infinity_rank
            if infinity_ade is not None:
                ade.append(infinity_ade)
            jacobian_signature = (
                tuple((factor.degree(), exponent) for factor, exponent in jacobian_discriminant.numerator().factor()),
                tuple((factor.degree(), exponent) for factor, exponent in jacobian_discriminant.denominator().factor()),
                jacobian.a4().numerator().degree(),
                jacobian.a4().denominator().degree(),
                jacobian.a6().numerator().degree(),
                jacobian.a6().denominator().degree(),
                tuple(local_valuations),
                infinity_orders,
                infinity_rank,
                root_rank,
                17 - root_rank,
                tuple(sorted(ade)),
            )
            record = (
                tuple(map(int, denominator_points)),
                tuple((int(a), int(value)) for a, value in points),
                tuple(int(value) for value in coefficients),
                polynomial.degree(),
                tuple((factor.degree(), exponent) for factor, exponent in factors),
                odd_part,
                jacobian_signature,
            )
            hits.append(record)
            if tuple(sorted(ade)) == target_ade:
                target_hits.append(record)
            print(
                f"MW2CRT|denominator_points={record[0]}|points={record[1]}"
                f"|q={record[2]}|degree={record[3]}"
                f"|factors={record[4]}|squarefree_degree={odd_part.degree()}"
                f"|squarefree={odd_part}|jacobian_signature={jacobian_signature}",
                flush=True,
            )

print(
    f"MW2CRT|p={args.p}|chord={args.chord}|q_degree={args.q_degree}|cases={len(choices)}"
    f"|denominators={len(denominator_choices)}|tested={tested}|hits={len(hits)}"
    f"|mobius_linear={int(args.mobius_linear)}"
    f"|mobius_quadratic_middle={args.mobius_quadratic_middle}"
    f"|target_hits={len(target_hits)}|status=PASS",
    flush=True,
)
