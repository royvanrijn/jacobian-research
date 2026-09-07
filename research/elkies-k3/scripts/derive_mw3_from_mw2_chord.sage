from sage.all import *
import argparse


parser = argparse.ArgumentParser(
    description="Test the degree-two chord map for the MW2 -> MW3 inverse neighbor."
)
parser.add_argument("--p", type=int, default=23, help="prime field; use 0 for QQ")
parser.add_argument(
    "--chord", choices=("p1", "raw-q", "s"), default="p1",
    help=(
        "use O+P1 after translation, the unreduced O+(3P1-4P2) divisor, "
        "or the Picard-20 third-section chord O+(-S)"
    ),
)
parser.add_argument(
    "--translated-pole", action="store_true",
    help="use w=1/(z_P1-k(t)) with k fixed by the exact P1-2P2 translation",
)
parser.add_argument(
    "--search-affine", action="store_true",
    help="over GF(p), search z=u+a*t+b corrections that lower the t-genus",
)
parser.add_argument(
    "--corrected", action="store_true",
    help="construct the corrected quartic and its binary-quartic Jacobian",
)
parser.add_argument(
    "--mw3-crt", action="store_true",
    help="construct the exact cubic MW3 completion through t=9/25, z=2",
)
parser.add_argument(
    "--mw3-target", action="store_true",
    help="construct the exact D5+A3+2A2+2A1 MW3 completion",
)
parser.add_argument("--section-scan-bound", type=int, default=0)
args = parser.parse_args()

F = QQ if args.p == 0 else GF(args.p)
Kt = FunctionField(F, "t")
t = Kt.gen()

if args.p == 0:
    A = (
        -F(32447500) / F(583443) * t**2
        - F(906250) / F(194481) * t**3
        + F(31250000) / F(194481) * t**4
        - F(19531250) / F(194481) * t**5
    )
    B = (
        F(300827000000) / F(2315685267) * t**3
        + F(340001171875) / F(1029193452) * t**4
        - F(498857421875) / F(257298363) * t**5
        + F(29541015625) / F(10501974) * t**6
        - F(152587890625) / F(85766121) * t**7
        + F(152587890625) / F(343064484) * t**8
    )
    p1x = 1 - F(800) / F(1323) * t + F(625) / F(147) * t**2
    p1y = (
        1 - F(400) / F(441) * t - F(394375) / F(18522) * t**2
        + F(484375) / F(9261) * t**3 - F(390625) / F(18522) * t**4
    )
    pole = t - F(16) / F(25)
    numerator_x = (
        F(77824) / F(33075) * t + F(12400) / F(1323) * t**2
        - F(30500) / F(1323) * t**3 + F(5000) / F(441) * t**4
    )
    numerator_y = (
        F(4096) / F(343) * t**2 + F(281408) / F(9261) * t**3
        - F(1517000) / F(9261) * t**4 + F(1296875) / F(6174) * t**5
        - F(1015625) / F(9261) * t**6 + F(390625) / F(18522) * t**7
    )
else:
    A = 13 * t**2 + 17 * t**3 + 11 * t**4 + 19 * t**5
    B = 20 * t**3 + 7 * t**4 + 3 * t**6 + 19 * t**7 + t**8
    p1x = 1 + 10 * t + 3 * t**2
    p1y = 1 + 15 * t + 14 * t**2 + 15 * t**3 + t**4
    pole = t - 8
    numerator_x = 15 * t + 6 * t**2 + 19 * t**3 + 8 * t**4
    numerator_y = 22 * t**2 + 14 * t**3 + 13 * t**4 + 2 * t**5 + 19 * t**6 + 22 * t**7

p2x = numerator_x / pole**2
p2y = numerator_y / pole**3
sx = -F(625)/F(441)*t**2 + F(3550)/F(1323)*t
sy = (
    F(390625)/F(18522)*t**4 - F(359375)/F(9261)*t**3
    + F(1875)/F(98)*t**2
)
E = EllipticCurve(Kt, [0, 0, 0, A, B])
P1 = E((p1x, p1y))
P2 = E((p2x, p2y))
S = E((sx, sy))
Q = {
    "p1": P1,
    "raw-q": 3 * P1 - 4 * P2,
    # The common line convention below passes through -Q.
    "s": -S,
}[args.chord]
qx, qy = Q[0], Q[1]
translation_pole = None
if args.translated_pole:
    if args.chord != "p1":
        raise ValueError("--translated-pole requires --chord p1")
    translation = P1 - 2*P2
    raw_horizontal = 3*P1 - 4*P2
    translated_left = -translation
    translated_right = raw_horizontal - translation

    def p1_chord_value(point):
        return (point[1] + p1y) / (point[0] - p1x)

    translation_pole = p1_chord_value(translated_left)
    assert translation_pole == p1_chord_value(translated_right)
    print(
        "MW2CHORD|translation=S=P1-2P2|pole_value_degree="
        f"{translation_pole.numerator().degree()},"
        f"{translation_pole.denominator().degree()}"
        f"|pole_value={translation_pole}",
        flush=True,
    )

print(
    "MW2CHORD|field={}|chord={}|x_num_deg={}|x_den_deg={}|y_num_deg={}|y_den_deg={}".format(
        "QQ" if args.p == 0 else f"GF({args.p})",
        args.chord,
        qx.numerator().degree(), qx.denominator().degree(),
        qy.numerator().degree(), qy.denominator().degree(),
    ),
    flush=True,
)

# Rebase from F(t) to F(u)(t), so the discriminant can be factored in t.
Fu = FractionField(PolynomialRing(F, "u"))
u = Fu.gen()
T = PolynomialRing(Fu, "T")
Tvar = T.gen()
KT = FractionField(T)


def rebase(value):
    numerator = value.numerator()
    denominator = value.denominator()
    top = sum(Fu(numerator[index]) * Tvar**index for index in range(numerator.degree() + 1))
    bottom = sum(Fu(denominator[index]) * Tvar**index for index in range(denominator.degree() + 1))
    return KT(top) / KT(bottom)


Rx = PolynomialRing(KT, "x")
x = Rx.gen()
A2, B2, qx2, qy2 = map(rebase, (A, B, qx, qy))
if translation_pole is None:
    chord_slope = KT(u)
else:
    chord_slope = rebase(translation_pole) + 1/KT(u)
line_y = chord_slope * (x - qx2) - qy2
cubic = line_y**2 - (x**3 + A2 * x + B2)
quadratic, remainder = cubic.quo_rem(x - qx2)
assert remainder == 0 and quadratic.degree() == 2
discriminant = quadratic.discriminant()
numerator = discriminant.numerator()
factorization = numerator.factor()
factorization_unit = factorization.unit()
factor_degrees = [(factor.degree(), exponent) for factor, exponent in factorization]
squarefree = prod(factor for factor, exponent in factorization if exponent % 2)
full_squarefree = factorization_unit * squarefree
print(
    f"MW2CHORD|quadratic=PASS|disc_num_degree={numerator.degree()}"
    f"|disc_den_degree={discriminant.denominator().degree()}"
    f"|factor_degrees={factor_degrees}|squarefree_degree={squarefree.degree()}",
    flush=True,
)
print(f"MW2CHORD|factorization_unit={factorization_unit}", flush=True)
print(f"MW2CHORD|discriminant_denominator={discriminant.denominator()}", flush=True)
if args.p != 0 or squarefree.degree() <= 6:
    print(f"MW2CHORD|squarefree={squarefree}", flush=True)

if args.search_affine:
    if args.p == 0:
        raise ValueError("--search-affine currently requires a prime field")
    BU = PolynomialRing(F, names=("U", "T"))
    Uvar, BT = BU.gens()

    def coefficient_at(value, substitution):
        top = value.numerator()
        bottom = value.denominator()
        assert bottom.degree() == 0
        return sum(BU(top[index]) * substitution**index for index in range(top.degree() + 1)) / F(bottom[0])

    affine_hits = []
    for slope in F:
        substitution = Uvar + slope * BT
        transformed = sum(
            coefficient_at(squarefree[index], substitution) * BT**index
            for index in range(squarefree.degree() + 1)
        )
        degree = transformed.degree(BT)
        if degree <= 4:
            affine_hits.append((slope, degree))
    print(f"MW2CHORD|affine_slope_hits={affine_hits}", flush=True)

    # Translation by b does not alter the degree at infinity, but record the
    # squarefree t-degree for every surviving slope/intercept pair.
    FU2 = FractionField(PolynomialRing(F, "U2"))
    U2 = FU2.gen()
    TT = PolynomialRing(FU2, "TT")
    TTvar = TT.gen()
    for slope, _ in affine_hits:
        for intercept in F:
            substitution = U2 + FU2(slope) * TTvar + FU2(intercept)
            transformed = sum(
                sum(
                    FU2(squarefree[index].numerator()[j]) * substitution**j
                    for j in range(squarefree[index].numerator().degree() + 1)
                ) / FU2(squarefree[index].denominator()[0]) * TTvar**index
                for index in range(squarefree.degree() + 1)
            )
            factors = transformed.factor()
            odd_part = prod(factor for factor, exponent in factors if exponent % 2)
            if odd_part.degree() <= 4:
                print(
                    f"MW2CHORD|affine_target|slope={slope}|intercept={intercept}"
                    f"|degree={transformed.degree()}|factor_degrees="
                    f"{[(factor.degree(), exponent) for factor, exponent in factors]}"
                    f"|squarefree_degree={odd_part.degree()}|squarefree={odd_part}",
                    flush=True,
                )

if args.corrected:
    slope = F(-50) / F(21)
    FU3 = FractionField(PolynomialRing(F, "U3"))
    U3 = FU3.gen()
    T3ring = PolynomialRing(FU3, "T3")
    T3 = T3ring.gen()

    def rebase_coefficient(value, substitution):
        top = value.numerator()
        bottom = value.denominator()
        assert bottom.degree() == 0
        return sum(FU3(top[index]) * substitution**index for index in range(top.degree() + 1)) / FU3(bottom[0])

    substitution = U3 + FU3(slope) * T3
    quartic = sum(
        rebase_coefficient(squarefree[index], substitution) * T3**index
        for index in range(squarefree.degree() + 1)
    )
    assert quartic.degree() == 4
    a, b, c, d, e = [quartic[index] for index in (4, 3, 2, 1, 0)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    jacobian = EllipticCurve(FU3, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    jacobian_discriminant = jacobian.discriminant()
    print(f"MW2CHORD|corrected_slope={slope}|quartic={quartic}", flush=True)
    print(
        f"MW2CHORD|jacobian|A={jacobian.a4()}|B={jacobian.a6()}"
        f"|disc_num_factor={jacobian_discriminant.numerator().factor()}"
        f"|disc_den_factor={jacobian_discriminant.denominator().factor()}",
        flush=True,
    )

if args.mw3_crt:
    if args.p != 0:
        raise ValueError("--mw3-crt is the exact QQ reconstruction")
    KU = FractionField(PolynomialRing(QQ, "V"))
    V = KU.gen()
    TVring = PolynomialRing(KU, "T4")
    T4 = TVring.gen()

    def rebase_at_v(value, substitution):
        top = value.numerator()
        bottom = value.denominator()
        assert bottom.degree() == 0
        return sum(KU(top[index]) * substitution**index for index in range(top.degree() + 1)) / KU(bottom[0])

    chord_value = (T4 - QQ(9) / 25) * V + 2
    completed = sum(
        rebase_at_v(squarefree[index], chord_value) * T4**index
        for index in range(squarefree.degree() + 1)
    )
    cubic, remainder = completed.quo_rem((T4 - QQ(9) / 25)**2)
    assert remainder == 0 and cubic.degree() == 3
    e, d, c, b, a = [cubic[index] for index in range(5)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    mw3 = EllipticCurve(KU, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    mw3_discriminant = mw3.discriminant()
    print(
        f"MW2CHORD|mw3_base=V=(z-2)/(t-9/25)|cubic={cubic}",
        flush=True,
    )
    print(
        f"MW2CHORD|mw3_jacobian|A={mw3.a4()}|B={mw3.a6()}"
        f"|disc_num_factor={mw3_discriminant.numerator().factor()}"
        f"|disc_den_factor={mw3_discriminant.denominator().factor()}",
        flush=True,
    )
    for factor, delta_order in mw3_discriminant.numerator().factor():
        a_order = mw3.a4().numerator().valuation(factor)
        b_order = mw3.a6().numerator().valuation(factor)
        print(
            f"MW2CHORD|mw3_finite_factor={factor}|ordA={a_order}"
            f"|ordB={b_order}|ordDelta={delta_order}",
            flush=True,
        )
    print(
        f"MW2CHORD|mw3_infinity|ordA={8 - mw3.a4().numerator().degree()}"
        f"|ordB={12 - mw3.a6().numerator().degree()}"
        f"|ordDelta={24 - mw3_discriminant.numerator().degree()}",
        flush=True,
    )
    for test_t in (QQ(0), QQ(1), QQ(9) / 25, QQ(49) / 25):
        value = KU(cubic(test_t))
        print(
            f"MW2CHORD|mw3_point_test|t={test_t}|square={int(value.is_square())}"
            + (f"|y={value.sqrt()}" if value.is_square() else ""),
            flush=True,
        )
    for m in range(-args.section_scan_bound, args.section_scan_bound + 1):
        for n in range(-args.section_scan_bound, args.section_scan_bound + 1):
            old_section = m * P1 + n * P2
            if old_section.is_zero() or old_section == P1:
                continue
            xr, yr = old_section[0], old_section[1]
            denominator = xr - p1x
            if denominator == 0:
                continue
            section_base = ((yr + p1y) / denominator - 2) / (t - QQ(9) / 25)
            degree = max(section_base.numerator().degree(), section_base.denominator().degree())
            if degree <= 2:
                print(
                    f"MW2CHORD|old_section=({m},{n})|new_base_degree={degree}"
                    f"|new_base={section_base}",
                    flush=True,
                )

if args.mw3_target:
    if args.p != 0:
        raise ValueError("--mw3-target is the exact QQ reconstruction")
    KW = FractionField(PolynomialRing(QQ, "W"))
    W = KW.gen()
    TWring = PolynomialRing(KW, "T5")
    T5 = TWring.gen()

    def rebase_at_w(value, substitution):
        top = value.numerator()
        bottom = value.denominator()
        assert bottom.degree() == 0
        return sum(KW(top[index]) * substitution**index for index in range(top.degree() + 1)) / KW(bottom[0])

    chord_value = T5 * (T5 - 1) * W - 1 - QQ(25) / 21 * T5
    completed = sum(
        rebase_at_w(squarefree[index], chord_value) * T5**index
        for index in range(squarefree.degree() + 1)
    )
    quartic, remainder = completed.quo_rem(T5**2 * (T5 - 1)**2)
    assert remainder == 0 and quartic.degree() == 4 and quartic[0] == 0
    e, d, c, b, a = [quartic[index] for index in range(5)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    target = EllipticCurve(KW, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    target_discriminant = target.discriminant()
    print(
        f"MW2CHORD|target_base=W=(z+1+(25/21)t)/(t*(t-1))|quartic={quartic}",
        flush=True,
    )
    print(
        f"MW2CHORD|target_jacobian|A={target.a4()}|B={target.a6()}"
        f"|disc_num_factor={target_discriminant.numerator().factor()}"
        f"|disc_den_factor={target_discriminant.denominator().factor()}",
        flush=True,
    )
    for factor, delta_order in target_discriminant.numerator().factor():
        print(
            f"MW2CHORD|target_finite_factor={factor}"
            f"|ordA={target.a4().numerator().valuation(factor)}"
            f"|ordB={target.a6().numerator().valuation(factor)}"
            f"|ordDelta={delta_order}",
            flush=True,
        )
    print(
        f"MW2CHORD|target_infinity|ordA={8 - target.a4().numerator().degree()}"
        f"|ordB={12 - target.a6().numerator().degree()}"
        f"|ordDelta={24 - target_discriminant.numerator().degree()}",
        flush=True,
    )
    natural = EllipticCurve(KW, [0, c, 0, b * d, a * d**2])
    natural_discriminant = natural.discriminant()
    print(
        f"MW2CHORD|target_natural|a2={natural.a2()}|a4={natural.a4()}"
        f"|a6={natural.a6()}|disc_num_factor={natural_discriminant.numerator().factor()}"
        f"|disc_den_factor={natural_discriminant.denominator().factor()}",
        flush=True,
    )
    twisted = EllipticCurve(KW, [0, 2 * c, 0, 4 * b * d, 8 * a * d**2])
    twist_points = []
    pulled_back_points = []
    for index, alpha in enumerate((QQ(1), -QQ(3) / 25, QQ(7) / 25), 1):
        t_candidate = KW(alpha) - KW(QQ(25) / 21) / W
        half_value = KW(quartic(t_candidate)) / 2
        assert half_value.is_square()
        rational_v = half_value.sqrt()
        point = twisted((2 * d / t_candidate, 4 * d * rational_v / t_candidate**2))
        twist_points.append(point)
        print(
            f"MW2CHORD|twist_section={index}|t={t_candidate}"
            f"|x={point[0]}|y={point[1]}",
            flush=True,
        )

        # Pull the section back through the chord quadratic and identify it in
        # the original MW2 group.  The original radicand differs from the
        # monic quartic by 2*(6250/441)^2.
        old_w = -Kt(QQ(25) / 21) / (t - Kt(alpha))
        old_z = t * (t - 1) * old_w - 1 - Kt(QQ(25) / 21) * t
        OldX = PolynomialRing(Kt, "old_x")
        old_x = OldX.gen()
        old_line_y = old_z * (old_x - p1x) - p1y
        old_cubic = old_line_y**2 - (old_x**3 + A*old_x + B)
        old_quadratic, old_remainder = old_cubic.quo_rem(old_x - p1x)
        assert old_remainder == 0 and old_quadratic.degree() == 2
        old_discriminant = old_quadratic.discriminant()
        assert old_discriminant.is_square()
        old_sqrt = old_discriminant.sqrt()
        old_a, old_b = old_quadratic[2], old_quadratic[1]
        old_roots = ((-old_b + old_sqrt)/(2*old_a), (-old_b - old_sqrt)/(2*old_a))
        identifications = []
        for root in old_roots:
            old_y = old_z * (root-p1x) - p1y
            old_point = E((root, old_y))
            pulled_back_points.append(old_point)
            print(
                f"MW2CHORD|pulled_back_section={len(pulled_back_points)}"
                f"|x_num_degree={old_point[0].numerator().degree()}"
                f"|x_den_degree={old_point[0].denominator().degree()}"
                f"|y_num_degree={old_point[1].numerator().degree()}"
                f"|y_den_degree={old_point[1].denominator().degree()}"
                f"|x={old_point[0]}|y={old_point[1]}",
                flush=True,
            )
            match = None
            for m in range(-6, 7):
                for n in range(-6, 7):
                    if old_point == m*P1 + n*P2:
                        match = (m, n)
            identifications.append(match)
        print(
            f"MW2CHORD|twist_section={index}|old_MW_identifications={tuple(identifications)}",
            flush=True,
        )
    print(
        f"MW2CHORD|target_twist|a2={twisted.a2()}|a4={twisted.a4()}"
        f"|a6={twisted.a6()}|sections=3",
        flush=True,
    )
    old_specialization_value = QQ(2)

    def specialize_kt(value):
        return QQ(value.numerator()(old_specialization_value)) / QQ(value.denominator()(old_specialization_value))

    old_specialized_curve = EllipticCurve(
        QQ, [0, 0, 0, specialize_kt(A), specialize_kt(B)]
    )
    old_specialized_points = [
        old_specialized_curve((specialize_kt(point[0]), specialize_kt(point[1])))
        for point in (P1, P2, pulled_back_points[0])
    ]
    old_height_matrix = old_specialized_curve.height_pairing_matrix(
        old_specialized_points, precision=128
    )
    old_regulator = old_height_matrix.det()
    assert old_regulator > 1e-20
    print(
        f"MW2CHORD|old_specialization=t:2|height_matrix={old_height_matrix}"
        f"|regulator={old_regulator}|independent=3",
        flush=True,
    )
    specialization_value = -QQ(25) / 21

    def specialize_kw(value):
        return QQ(value.numerator()(specialization_value)) / QQ(value.denominator()(specialization_value))

    specialized_curve = EllipticCurve(
        QQ,
        [0, specialize_kw(twisted.a2()), 0, specialize_kw(twisted.a4()), specialize_kw(twisted.a6())],
    )
    specialized_points = [
        specialized_curve((specialize_kw(point[0]), specialize_kw(point[1])))
        for point in twist_points
    ]
    height_matrix = specialized_curve.height_pairing_matrix(specialized_points, precision=128)
    regulator = height_matrix.det()
    assert regulator > 1e-20
    print(
        f"MW2CHORD|twist_specialization=W:-25/21|curve={specialized_curve.global_minimal_model()}"
        f"|height_matrix={height_matrix}|regulator={regulator}"
        f"|independent=3",
        flush=True,
    )
    for label, t_candidate in (
        ("mobius_1", (W + 1) / W),
        ("mobius_10", (10 * W + 1) / W),
        ("mobius_15", (15 * W + 1) / W),
    ):
        quartic_value = KW(quartic(t_candidate))
        print(
            f"MW2CHORD|target_point_test={label}|square={int(quartic_value.is_square())}"
            + (f"|t={t_candidate}|v={quartic_value.sqrt()}" if quartic_value.is_square() else ""),
            flush=True,
        )
    RU = PolynomialRing(QQ, "U6")
    U6 = RU.gen()

    def substitute_w(polynomial):
        numerator = polynomial.numerator()
        denominator = polynomial.denominator()
        assert denominator.degree() == 0
        return sum(
            QQ(numerator[index]) * (QQ(25) / 21 * U6)**index
            for index in range(numerator.degree() + 1)
        ) / QQ(denominator[0])

    normalized_a2 = RU(substitute_w(natural.a2()))
    normalized_a4 = RU(substitute_w(natural.a4()))
    normalized_a6 = RU(substitute_w(natural.a6()))
    print(
        f"MW2CHORD|target_base_normalized=U6=21W/25"
        f"|a2_factor={normalized_a2.factor()}"
        f"|a4_factor={normalized_a4.factor()}"
        f"|a6_factor={normalized_a6.factor()}",
        flush=True,
    )
    for index, point in enumerate(twist_points, 1):
        normalized_x = RU(substitute_w(point[0]))
        normalized_y = RU(substitute_w(point[1]))
        print(
            f"MW2CHORD|normalized_twist_section={index}"
            f"|x_factor={normalized_x.factor()}|y_factor={normalized_y.factor()}",
            flush=True,
        )
print("MW2CHORD|status=PASS", flush=True)
