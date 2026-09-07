#!/usr/bin/env sage
"""Search a compact first-q4 Riemann--Roch coordinate at the CM24 anchor.

The exact chamber certificate puts the first neighbor in a degree-two,
zero-MW-projection relative RR space.  This bounded equation diagnostic tests

    U = (x + f(T)) / T^k

on the exact rational CM24 q80 surface.  For each ``k<=--max-k`` it imposes
that, after clearing a common even power ``T^(2m)``, all four coefficients of
the cubic in U have T-degree at most four.  Any surviving solution therefore
gives a binary quartic genus-one model immediately.

This monomial-denominator ansatz is not a completeness theorem for the full
RR space: an empty scan leaves denominators supported at several old fibers
or unresolved-component compensators open.
"""

import argparse
from collections import Counter

from sage.all import GF, PolynomialRing, QQ, ZZ


parser = argparse.ArgumentParser()
parser.add_argument("--max-k", type=int, default=8)
parser.add_argument(
    "--general-degree",
    type=int,
    default=0,
    help="also solve the monic general-denominator chart of this degree",
)
parser.add_argument(
    "--lift-pair-survivors",
    action="store_true",
    help="solve surviving discriminant-pair charts exactly over QQ",
)
parser.add_argument(
    "--discriminant-pairs-mod",
    type=int,
    default=0,
    help="screen all (g,H) discriminant-divisor pairs over the given prime",
)
parser.add_argument(
    "--discriminant-divisors",
    action="store_true",
    help="scan all degree-four monic divisors of the CM24 discriminant",
)
args = parser.parse_args()

base = PolynomialRing(QQ, "T")
T0 = base.gen()
A0 = T0**2 * (-3 + QQ(9) / 4 * T0 - QQ(9) / 4 * T0**2 + QQ(9) / 4 * T0**3)
B0 = T0**3 * (
    2
    - QQ(315) / 32 * T0
    + 9 * T0**2
    - QQ(9) / 16 * T0**3
    - QQ(27) / 32 * T0**5
)


def outside_coefficients(polynomial, low, high):
    return tuple(
        coefficient
        for degree, coefficient in enumerate(polynomial.list())
        if coefficient and not low <= degree <= high
    )


hits = []
charts = 0
for k in range(args.max_k + 1):
    # The U^3 coefficient is T^(3k), so only these square clearings can
    # possibly leave a quartic.
    for m in range(max(0, (3 * k - 4 + 1) // 2), 3 * k // 2 + 1):
        low = 2 * m
        high = low + 4
        if not low <= 3 * k <= high:
            continue
        support = tuple(
            degree
            for degree in range(max(0, low - 2 * k), high - 2 * k + 1)
        )
        # f -> f+c*T^k only translates U. Fix that gauge when present.
        free_support = tuple(degree for degree in support if degree != k)
        names = tuple(f"c{degree}" for degree in free_support)
        coefficient_ring = PolynomialRing(QQ, names=names or ("dummy",))
        coefficients = dict(zip(free_support, coefficient_ring.gens()))
        if not names:
            coefficients = {}
        local = PolynomialRing(coefficient_ring, "T")
        T = local.gen()
        A = local(A0)
        B = local(B0)
        f = sum(coefficients.get(degree, coefficient_ring(0)) * T**degree for degree in support)
        g = T**k
        u3 = g**3
        u2 = -3 * g**2 * f
        u1 = 3 * g * f**2 + A * g
        u0 = -f**3 - A * f + B
        equations = []
        for polynomial in (u3, u2, u1, u0):
            equations.extend(outside_coefficients(polynomial, low, high))
        # Remove constants already known to vanish and reject contradictions.
        if any(not equation.variables() and equation for equation in equations):
            continue
        equations = tuple(equation for equation in equations if equation)
        ideal = coefficient_ring.ideal(equations)
        charts += 1
        if ideal.is_one():
            continue
        dimension = ideal.dimension()
        solutions = ()
        if dimension == 0:
            try:
                solutions = tuple(ideal.variety(ring=QQ))
            except (ValueError, TypeError):
                solutions = ()
        print(
            f"Q80FIRSTQ4RR|k={k}|m={m}|support={support}|"
            f"variables={len(names)}|equations={len(equations)}|"
            f"dimension={dimension}|rational_solutions={len(solutions)}",
            flush=True,
        )
        for solution in solutions:
            f_solution = base(
                sum(
                    QQ(solution[coefficients[degree]]) * T0**degree
                    for degree in free_support
                )
            )
            hits.append((k, m, f_solution))
            print(
                f"Q80FIRSTQ4RR|hit|k={k}|m={m}|f={f_solution}",
                flush=True,
            )


def binary_quartic_jacobian(k, m, f):
    coefficients = PolynomialRing(QQ, "U")
    U = coefficients.gen()
    local = PolynomialRing(coefficients, "T")
    T = local.gen()
    A = local(A0)
    B = local(B0)
    f = local(f)
    quartic, remainder = (
        (T**k * U - f) ** 3 + A * (T**k * U - f) + B
    ).quo_rem(T ** (2 * m))
    assert remainder == 0 and quartic.degree() <= 4
    e, d, c, b, a = [quartic[index] for index in range(5)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    jacobian_a = -27 * invariant_i
    jacobian_b = -27 * invariant_j
    discriminant = 4 * jacobian_a**3 + 27 * jacobian_b**2
    return quartic, jacobian_a, jacobian_b, discriminant


def geometric_root_rank(jacobian_a, jacobian_b, discriminant):
    rank = ZZ(0)
    signature = []
    for factor, delta_order in discriminant.factor():
        a_order = jacobian_a.valuation(factor)
        b_order = jacobian_b.valuation(factor)
        degree = factor.degree()
        if a_order == 0 or b_order == 0:
            local_rank = max(ZZ(0), ZZ(delta_order) - 1)
            fiber = f"I{delta_order}"
        elif a_order >= 2 and b_order >= 3 and delta_order >= 6:
            local_rank = ZZ(delta_order) - 2
            fiber = f"I{delta_order - 6}*"
        else:
            raise RuntimeError(
                f"unclassified finite fiber valuations {(a_order, b_order, delta_order)}"
            )
        rank += degree * local_rank
        signature.append((str(factor), degree, a_order, b_order, ZZ(delta_order), fiber))
    infinity = (
        ZZ(8 - jacobian_a.degree()),
        ZZ(12 - jacobian_b.degree()),
        ZZ(24 - discriminant.degree()),
    )
    a_order, b_order, delta_order = infinity
    if a_order == 0 or b_order == 0:
        infinity_rank = max(ZZ(0), delta_order - 1)
        infinity_fiber = f"I{delta_order}"
    elif a_order >= 2 and b_order >= 3 and delta_order >= 6:
        infinity_rank = delta_order - 2
        infinity_fiber = f"I{delta_order - 6}*"
    else:
        raise RuntimeError(f"unclassified infinity valuations {infinity}")
    rank += infinity_rank
    signature.append(("infinity", 1, a_order, b_order, delta_order, infinity_fiber))
    return rank, tuple(signature)


def is_target_cm_signature(signature):
    fibers = Counter()
    for _, degree, _, _, _, fiber in signature:
        fibers[fiber] += degree
    # The rank-19 family has I5*+I5+8I1.  At the CM24 boundary the finite I5
    # becomes I6, while two conjugate pairs of residual roots become I2.  The
    # root rank therefore rises by three (and the child MW rank drops by two),
    # even though the old q80 fibration acquires only one new NS direction.
    return fibers == Counter({"I5*": 1, "I6": 1, "I2": 2, "I1": 3})


target_signature_hits = 0
for k, m, f in hits:
    quartic, jacobian_a, jacobian_b, discriminant = binary_quartic_jacobian(k, m, f)
    root_rank, signature = geometric_root_rank(
        jacobian_a, jacobian_b, discriminant
    )
    # The generic target has D9+A4 and eight I1 fibers.  At CM24 it specializes
    # to D9+A5+2A1.  Root rank alone is not sufficient: other compact CM
    # pencils can have a different ADE type.
    target_signature = is_target_cm_signature(signature)
    target_signature_hits += ZZ(target_signature)
    print(
        f"Q80FIRSTQ4RR|classified|k={k}|m={m}|f={f}|quartic={quartic}|"
        f"jacobian_A={jacobian_a.factor()}|jacobian_B={jacobian_b.factor()}|"
        f"Delta={discriminant.factor()}|root_rank={root_rank}|fibers={signature}|"
        f"cm_root_rank16={ZZ(root_rank == 16)}|"
        f"target_cm_signature={ZZ(target_signature)}",
        flush=True,
    )

general_hits = []
if args.general_degree:
    degree = ZZ(args.general_degree)
    assert 1 <= degree <= 4
    names = tuple(f"g{i}" for i in range(degree)) + tuple(
        f"f{i}" for i in range(5) if i != degree
    )
    coefficient_ring = PolynomialRing(QQ, names=names, order="degrevlex")
    generators = dict(zip(names, coefficient_ring.gens()))
    local = PolynomialRing(coefficient_ring, "T")
    T = local.gen()
    g = T**degree + sum(generators[f"g{i}"] * T**i for i in range(degree))
    # f -> f+c*g translates U, so set the degree-d coefficient of f to zero.
    f = sum(
        generators[f"f{i}"] * T**i
        for i in range(5)
        if i != degree
    )
    A = local(A0)
    B = local(B0)
    linear_numerator = 3 * f**2 + A
    constant_numerator = -f**3 - A * f + B
    linear_quotient, linear_remainder = linear_numerator.quo_rem(g)
    constant_quotient, constant_remainder = constant_numerator.quo_rem(g**2)
    equations = tuple(
        coefficient
        for polynomial in (linear_remainder, constant_remainder)
        for coefficient in polynomial.list()
        if coefficient
    )
    if linear_quotient.degree() > 4:
        equations += tuple(linear_quotient[index] for index in range(5, linear_quotient.degree() + 1))
    if constant_quotient.degree() > 4:
        equations += tuple(constant_quotient[index] for index in range(5, constant_quotient.degree() + 1))
    ideal = coefficient_ring.ideal(equations)
    dimension = ideal.dimension()
    solutions = ()
    if dimension == 0:
        try:
            solutions = tuple(ideal.variety(ring=QQ))
        except (ValueError, TypeError):
            solutions = ()
    print(
        f"Q80FIRSTQ4RR|general_degree={degree}|variables={len(names)}|"
        f"equations={len(equations)}|dimension={dimension}|"
        f"quotient_degrees={linear_quotient.degree()},{constant_quotient.degree()}|"
        f"rational_solutions={len(solutions)}",
        flush=True,
    )
    rational = PolynomialRing(QQ, "T")
    TT = rational.gen()
    for solution in solutions:
        g_solution = rational(
            TT**degree
            + sum(QQ(solution[generators[f"g{i}"]]) * TT**i for i in range(degree))
        )
        f_solution = rational(
            sum(
                QQ(solution[generators[f"f{i}"]]) * TT**i
                for i in range(5)
                if i != degree
            )
        )
        general_hits.append((g_solution, f_solution))
        print(
            f"Q80FIRSTQ4RR|general_hit|g={g_solution}|f={f_solution}",
            flush=True,
        )

divisor_hits = []
if args.discriminant_divisors:
    from itertools import product

    discriminant_old = 4 * A0**3 + 27 * B0**2
    factors = tuple(discriminant_old.factor())
    divisors = []
    for exponents in product(*(range(exponent + 1) for _, exponent in factors)):
        candidate = base.prod(
            factor**chosen
            for (factor, _), chosen in zip(factors, exponents)
        )
        if candidate.degree() == 4:
            divisors.append(candidate.monic())
    divisors = sorted(set(divisors), key=str)
    assert len(divisors) == 15

    prime = ZZ(101)
    finite = PolynomialRing(GF(prime), names=("f0", "f1", "f2", "f3"), order="degrevlex")
    finite_T = PolynomialRing(finite, "T")
    Tf = finite_T.gen()
    f_finite = sum(finite.gen(i) * Tf**i for i in range(4))

    for divisor_index, g_rational in enumerate(divisors):
        g_finite = finite_T(
            [GF(prime)(coefficient) for coefficient in g_rational]
        )
        A_finite = finite_T([GF(prime)(coefficient) for coefficient in A0])
        B_finite = finite_T([GF(prime)(coefficient) for coefficient in B0])
        finite_equations = tuple(
            coefficient
            for polynomial in (
                (3 * f_finite**2 + A_finite) % g_finite,
                (-f_finite**3 - A_finite * f_finite + B_finite) % (g_finite**2),
            )
            for coefficient in polynomial.list()
            if coefficient
        )
        finite_ideal = finite.ideal(finite_equations)
        finite_solutions = () if finite_ideal.is_one() else tuple(finite_ideal.variety())
        print(
            f"Q80FIRSTQ4RR|divisor={divisor_index}|g={g_rational}|"
            f"mod101_solutions={len(finite_solutions)}",
            flush=True,
        )
        if not finite_solutions:
            continue

        rational_ring = PolynomialRing(
            QQ, names=("f0", "f1", "f2", "f3"), order="degrevlex"
        )
        rational_T = PolynomialRing(rational_ring, "T")
        Tr = rational_T.gen()
        f_rational = sum(rational_ring.gen(i) * Tr**i for i in range(4))
        g_local = rational_T(g_rational)
        equations = tuple(
            coefficient
            for polynomial in (
                (3 * f_rational**2 + rational_T(A0)) % g_local,
                (-f_rational**3 - rational_T(A0) * f_rational + rational_T(B0))
                % (g_local**2),
            )
            for coefficient in polynomial.list()
            if coefficient
        )
        ideal = rational_ring.ideal(equations)
        solutions = () if ideal.is_one() else tuple(ideal.variety(ring=QQ))
        for solution in solutions:
            f_solution = base(
                sum(QQ(solution[rational_ring.gen(i)]) * T0**i for i in range(4))
            )
            quartic_ring = PolynomialRing(QQ["U"], "T")
            Tq = quartic_ring.gen()
            Uq = quartic_ring.base_ring().gen()
            quartic, remainder = (
                (quartic_ring(g_rational) * Uq - quartic_ring(f_solution)) ** 3
                + quartic_ring(A0)
                * (quartic_ring(g_rational) * Uq - quartic_ring(f_solution))
                + quartic_ring(B0)
            ).quo_rem(quartic_ring(g_rational) ** 2)
            assert remainder == 0 and quartic.degree() <= 4
            e, d, c, b, a = [quartic[index] for index in range(5)]
            invariant_i = 12 * a * e - 3 * b * d + c**2
            invariant_j = (
                72 * a * c * e
                + 9 * b * c * d
                - 27 * a * d**2
                - 27 * b**2 * e
                - 2 * c**3
            )
            jacobian_a = -27 * invariant_i
            jacobian_b = -27 * invariant_j
            discriminant = 4 * jacobian_a**3 + 27 * jacobian_b**2
            root_rank, signature = geometric_root_rank(
                jacobian_a, jacobian_b, discriminant
            )
            divisor_hits.append((g_rational, f_solution, root_rank))
            print(
                f"Q80FIRSTQ4RR|divisor_hit|g={g_rational}|f={f_solution}|"
                f"quartic={quartic}|Delta={discriminant.factor()}|"
                f"root_rank={root_rank}|fibers={signature}|"
                f"cm_root_rank16={ZZ(root_rank == 16)}|"
                f"target_cm_signature={ZZ(is_target_cm_signature(signature))}",
                flush=True,
            )

pair_survivors = []
if args.discriminant_pairs_mod:
    from itertools import product

    prime = ZZ(args.discriminant_pairs_mod)
    assert prime.is_prime() and prime not in (2, 3)
    discriminant_old = 4 * A0**3 + 27 * B0**2
    factors = tuple(discriminant_old.factor())
    pairs = []
    for g_exponents in product(*(range(exponent + 1) for _, exponent in factors)):
        g = base.prod(
            factor**chosen
            for (factor, _), chosen in zip(factors, g_exponents)
        ).monic()
        if not 1 <= g.degree() <= 4:
            continue
        for h_exponents in product(
            *(range((3 * chosen) // 2 + 1) for chosen in g_exponents)
        ):
            h = base.prod(
                factor**chosen
                for (factor, _), chosen in zip(factors, h_exponents)
            ).monic()
            quotient_degree = 3 * g.degree() - 2 * h.degree()
            if 0 <= quotient_degree <= 4:
                pairs.append((g, h))
    pairs = sorted(set(pairs), key=lambda pair: (str(pair[0]), str(pair[1])))
    assert len(pairs) == 99

    for pair_index, (g_rational, h_rational) in enumerate(pairs):
        degree = g_rational.degree()
        names = tuple(f"f{i}" for i in range(5) if i != degree)
        finite = PolynomialRing(GF(prime), names=names, order="degrevlex")
        finite_T = PolynomialRing(finite, "T")
        Tf = finite_T.gen()
        f_finite = sum(
            finite.gen(index) * Tf**power
            for index, power in enumerate(i for i in range(5) if i != degree)
        )
        g_finite = finite_T([GF(prime)(coefficient) for coefficient in g_rational])
        h_finite = finite_T([GF(prime)(coefficient) for coefficient in h_rational])
        A_finite = finite_T([GF(prime)(coefficient) for coefficient in A0])
        B_finite = finite_T([GF(prime)(coefficient) for coefficient in B0])
        denominator = h_finite**2
        polynomials = (
            g_finite**3,
            -3 * g_finite**2 * f_finite,
            3 * g_finite * f_finite**2 + A_finite * g_finite,
            -f_finite**3 - A_finite * f_finite + B_finite,
        )
        equations = []
        for polynomial in polynomials:
            quotient, remainder = polynomial.quo_rem(denominator)
            equations.extend(coefficient for coefficient in remainder if coefficient)
            if quotient.degree() > 4:
                equations.extend(
                    quotient[index]
                    for index in range(5, quotient.degree() + 1)
                    if quotient[index]
                )
        ideal = finite.ideal(tuple(equations))
        if ideal.is_one():
            continue
        dimension = ideal.dimension()
        if dimension != 0:
            pair_survivors.append((g_rational, h_rational, dimension, None))
            print(
                f"Q80FIRSTQ4RR|pair_survivor={pair_index}|g={g_rational}|"
                f"H={h_rational}|mod={prime}|dimension={dimension}|solutions=NA",
                flush=True,
            )
            continue
        solutions = tuple(ideal.variety())
        if solutions:
            pair_survivors.append((g_rational, h_rational, dimension, len(solutions)))
            print(
                f"Q80FIRSTQ4RR|pair_survivor={pair_index}|g={g_rational}|"
                f"H={h_rational}|mod={prime}|dimension=0|solutions={len(solutions)}",
                flush=True,
            )
    print(
        f"Q80FIRSTQ4RR|pair_charts={len(pairs)}|mod={prime}|"
        f"pair_survivors={len(pair_survivors)}",
        flush=True,
    )

pair_exact_hits = []
pair_target_hits = 0
if args.lift_pair_survivors:
    assert args.discriminant_pairs_mod and pair_survivors
    for survivor_index, (g_rational, h_rational, _, _) in enumerate(pair_survivors):
        degree = g_rational.degree()
        powers = tuple(i for i in range(5) if i != degree)
        names = tuple(f"f{i}" for i in powers)
        rational_ring = PolynomialRing(QQ, names=names, order="degrevlex")
        rational_T = PolynomialRing(rational_ring, "T")
        T = rational_T.gen()
        f = sum(
            rational_ring.gen(index) * T**power
            for index, power in enumerate(powers)
        )
        g = rational_T(g_rational)
        h = rational_T(h_rational)
        A = rational_T(A0)
        B = rational_T(B0)
        equations = []
        for polynomial in (
            g**3,
            -3 * g**2 * f,
            3 * g * f**2 + A * g,
            -f**3 - A * f + B,
        ):
            quotient, remainder = polynomial.quo_rem(h**2)
            equations.extend(coefficient for coefficient in remainder if coefficient)
            if quotient.degree() > 4:
                equations.extend(
                    quotient[index]
                    for index in range(5, quotient.degree() + 1)
                    if quotient[index]
                )
        ideal = rational_ring.ideal(tuple(equations))
        solutions = () if ideal.is_one() else tuple(ideal.variety(ring=QQ))
        print(
            f"Q80FIRSTQ4RR|pair_exact={survivor_index}|g={g_rational}|"
            f"H={h_rational}|rational_solutions={len(solutions)}",
            flush=True,
        )
        for solution in solutions:
            f_solution = base(
                sum(
                    QQ(solution[rational_ring.gen(index)]) * T0**power
                    for index, power in enumerate(powers)
                )
            )
            quartic_ring = PolynomialRing(QQ["U"], "T")
            Tq = quartic_ring.gen()
            Uq = quartic_ring.base_ring().gen()
            quartic, remainder = (
                (quartic_ring(g_rational) * Uq - quartic_ring(f_solution)) ** 3
                + quartic_ring(A0)
                * (quartic_ring(g_rational) * Uq - quartic_ring(f_solution))
                + quartic_ring(B0)
            ).quo_rem(quartic_ring(h_rational) ** 2)
            assert remainder == 0 and quartic.degree() <= 4
            e, d, c, b, a = [quartic[index] for index in range(5)]
            invariant_i = 12 * a * e - 3 * b * d + c**2
            invariant_j = (
                72 * a * c * e
                + 9 * b * c * d
                - 27 * a * d**2
                - 27 * b**2 * e
                - 2 * c**3
            )
            jacobian_a = -27 * invariant_i
            jacobian_b = -27 * invariant_j
            discriminant = 4 * jacobian_a**3 + 27 * jacobian_b**2
            root_rank, signature = geometric_root_rank(
                jacobian_a, jacobian_b, discriminant
            )
            target_signature = is_target_cm_signature(signature)
            pair_target_hits += ZZ(target_signature)
            pair_exact_hits.append(
                (g_rational, h_rational, f_solution, root_rank, target_signature)
            )
            print(
                f"Q80FIRSTQ4RR|pair_exact_hit|g={g_rational}|H={h_rational}|"
                f"f={f_solution}|quartic={quartic}|Delta={discriminant.factor()}|"
                f"root_rank={root_rank}|fibers={signature}|"
                f"cm_root_rank16={ZZ(root_rank == 16)}|"
                f"target_cm_signature={ZZ(target_signature)}",
                flush=True,
            )

print(
    f"Q80FIRSTQ4RR|charts={charts}|hits={len(hits)}|"
    f"target_signature_hits={target_signature_hits}|general_hits={len(general_hits)}|"
    f"divisor_hits={len(divisor_hits)}|"
    f"pair_survivors={len(pair_survivors)}|"
    f"pair_exact_hits={len(pair_exact_hits)}|"
    f"pair_target_hits={pair_target_hits}|"
    f"max_k={args.max_k}|"
    f"status=PASS_BOUNDED_ANSATZ",
    flush=True,
)
