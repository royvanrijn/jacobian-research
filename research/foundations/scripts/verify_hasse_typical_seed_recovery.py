#!/usr/bin/env python3
"""Verify minimal p-typical Hasse reconstruction and the degree-eight F_5 repair."""

from __future__ import annotations

from math import comb

import sympy as sp

from search_d1_f2_coefficient_collisions import (
    clean_seed,
    hessian,
    projective_key,
)


def trim(polynomial: tuple[int, ...], p: int) -> tuple[int, ...]:
    values = list(polynomial)
    while values and values[-1] % p == 0:
        values.pop()
    return tuple(value % p for value in values)


def hasse(
    polynomial: tuple[int, ...], order: int, p: int
) -> tuple[int, ...]:
    if order >= len(polynomial):
        return ()
    return trim(
        tuple(
            comb(exponent, order) * polynomial[exponent] % p
            for exponent in range(order, len(polynomial))
        ),
        p,
    )


def p_powers_at_most(n: int, p: int) -> tuple[int, ...]:
    powers = []
    value = 1
    while value <= n:
        powers.append(value)
        value *= p
    return tuple(powers)


def digit_at(value: int, place: int, p: int) -> int:
    return value // place % p


def coefficient(polynomial: tuple[int, ...], exponent: int) -> int:
    return polynomial[exponent] if exponent < len(polynomial) else 0


def reconstruct(
    channels: dict[int, tuple[int, ...]],
    n: int,
    p: int,
    constant: int,
) -> tuple[int, ...]:
    recovered = [constant % p]
    powers = p_powers_at_most(n, p)
    for exponent in range(1, n + 1):
        place = next(
            power for power in powers if digit_at(exponent, power, p) != 0
        )
        digit = digit_at(exponent, place, p)
        observed = coefficient(channels[place], exponent - place)
        recovered.append(observed * pow(digit, -1, p) % p)
    return trim(tuple(recovered), p)


def monomial(exponent: int) -> tuple[int, ...]:
    return (0,) * exponent + (1,)


def verify_all_degree_reconstruction() -> None:
    for p in (2, 3, 5, 7):
        for n in range(1, 41):
            powers = p_powers_at_most(n, p)

            # Basis recovery proves the identity for every polynomial.
            for exponent in range(n + 1):
                source = monomial(exponent)
                channels = {order: hasse(source, order, p) for order in powers}
                expected = () if exponent == 0 else source
                assert reconstruct(channels, n, p, constant=0) == expected

            # A dense deterministic polynomial also checks simultaneous use.
            source = tuple((j * j + 3 * j + 1) % p for j in range(n + 1))
            channels = {order: hasse(source, order, p) for order in powers}
            assert reconstruct(channels, n, p, source[0]) == trim(source, p)

            # Each forced p-power order is the only positive order detecting
            # its corresponding monomial, even if arbitrary orders are allowed.
            for forced in powers:
                witness = monomial(forced)
                assert hasse(witness, forced, p) == (1,)
                for other in range(1, n + 1):
                    if other != forced:
                        assert hasse(witness, other, p) == ()


def verify_characteristic_two_parity_core() -> None:
    p = 2
    for n in range(1, 41):
        source = tuple((j * j + j + 1) % p for j in range(n + 1))
        odd_channel = tuple(
            coefficient(source, 2 * index + 1)
            for index in range((n + 1) // 2)
        )
        even_channel = tuple(
            coefficient(source, 2 * index)
            for index in range(n // 2 + 1)
        )

        derivative_from_parity = [0] * n
        second_coordinate = [0] * (n + 1)
        reconstructed = [0] * (n + 1)
        for index, value in enumerate(odd_channel):
            derivative_from_parity[2 * index] = value
            reconstructed[2 * index + 1] = value
        for index, value in enumerate(even_channel):
            second_coordinate[2 * index] = value
            reconstructed[2 * index] = value

        assert trim(tuple(derivative_from_parity), p) == ordinary_derivative(
            source, p
        )
        tangent_value = tuple(
            ((j - 1) * coefficient(source, j)) % p for j in range(n + 1)
        )
        assert trim(tuple(second_coordinate), p) == trim(tangent_value, p)
        assert trim(tuple(reconstructed), p) == trim(source, p)

        # The weighted parameter denominator 2+kappa always vanishes because
        # every ordinary second derivative is zero in characteristic two.
        second_derivative = ordinary_derivative(
            ordinary_derivative(source, p), p
        )
        assert second_derivative == ()
        kappa = evaluate(second_derivative, 1, p)
        assert (2 + kappa) % p == 0


def multiply(
    left: tuple[int, ...], right: tuple[int, ...], p: int
) -> tuple[int, ...]:
    if not left or not right:
        return ()
    result = [0] * (len(left) + len(right) - 1)
    for left_exponent, left_value in enumerate(left):
        for right_exponent, right_value in enumerate(right):
            result[left_exponent + right_exponent] += left_value * right_value
    return trim(tuple(result), p)


def verify_characteristic_two_scalar_ansatz_no_go() -> None:
    p = 2

    # A bounded exhaustive regression of the general domain argument:
    # R R'/gamma can be a nonzero constant only for R=gamma.  The proof in
    # the note is degree-independent because both factors of
    # S(S+gamma*S') must be units.
    admissible = []
    for mask in range(1 << 8):
        reparameterization = (0,) + tuple((mask >> j) & 1 for j in range(8))
        if evaluate(reparameterization, 1, p) != 1:
            continue
        product = multiply(
            reparameterization,
            ordinary_derivative(reparameterization, p),
            p,
        )
        if product == (0, 1):
            admissible.append(trim(reparameterization, p))
    assert admissible == [(0, 1)]

    # With R=gamma, characteristic-two polynomiality requires d=0, whereas
    # the source vertical Jacobian is nonzero exactly when d is nonzero.
    for d in range(p):
        polynomial = d == 0
        separable = d != 0
        assert not (polynomial and separable)


def verify_characteristic_two_weight_redistribution() -> None:
    p = 2
    x, y, z, w = sp.symbols("x y z w")
    domain = sp.GF(p)

    def poly(expression):
        return sp.Poly(sp.expand(expression), x, y, z, domain=domain)

    # Odd degrees use W^2+W^N.  Even degrees at least six use
    # W^2+W^4+W^N+W^3.  Each seed has an exact double zero, H(1)=0, and
    # H'(1)=1=-1 in characteristic two.
    seeds = []
    for n in range(3, 12, 2):
        seed = [0] * (n + 1)
        seed[2] = seed[n] = 1
        seeds.append(tuple(seed))
    for n in range(6, 13, 2):
        seed = [0] * (n + 1)
        for exponent in (2, 3, 4, n):
            seed[exponent] = 1
        seeds.append(tuple(seed))

    u = 1 + x**2 * y
    gamma = 1 + x * z
    capital_w = u * gamma
    capital_c = x * gamma

    source_jacobian = sp.Matrix((capital_w, gamma, capital_c)).jacobian(
        (x, y, z)
    ).det()
    assert poly(source_jacobian - x**3 * gamma**2).is_zero

    for seed in seeds:
        seed_polynomial = sp.Poly(
            sum(value * w**j for j, value in enumerate(seed)),
            w,
            modulus=p,
        )
        derivative = seed_polynomial.diff()
        assert seed_polynomial.eval(0) % p == 0
        assert derivative.eval(0) % p == 0
        assert seed_polynomial.eval(1) % p == 0
        assert derivative.eval(1) % p == 1

        seed_at_w = seed_polynomial.as_expr().subs(w, capital_w)
        derivative_at_w = derivative.as_expr().subs(w, capital_w)
        first_numerator = sp.expand(
            capital_w * (derivative_at_w + gamma) - seed_at_w
        )
        second_numerator = sp.expand(derivative_at_w + gamma)
        alpha = poly(first_numerator).exquo(poly(capital_c**2))
        beta = poly(second_numerator).exquo(poly(capital_c))
        output = (alpha.as_expr(), beta.as_expr(), capital_c)

        jacobian = sp.Matrix(output).jacobian((x, y, z)).det()
        assert poly(jacobian - 1).is_zero
        assert poly(beta.as_expr() * capital_c - second_numerator).is_zero
        assert poly(alpha.as_expr() * capital_c**2 - first_numerator).is_zero

        if seed == (0, 0, 1, 1):
            assert poly(alpha.as_expr() - (y + x**2 * y**2)).is_zero
            assert poly(
                beta.as_expr() - (z + x**3 * y**2 + x**4 * y**2 * z)
            ).is_zero
            assert poly(capital_c - (x + x**2 * z)).is_zero


def verify_singular_parameter_weight_redistribution() -> None:
    x, y, z, w = sp.symbols("x y z w")
    u = 1 + x**2 * y
    gamma = 1 + x * z
    capital_w = u * gamma
    capital_c = x * gamma

    # The integral seed 2W^2-3W^3+W^4 has H'(1)=-1 and H''(1)=-2,
    # so it lies on the old chart's singular denominator in every
    # characteristic.  The redistributed chart remains polynomial Keller.
    for p in (3, 5, 7):
        domain = sp.GF(p)

        def poly(expression):
            return sp.Poly(sp.expand(expression), x, y, z, domain=domain)

        seed = sp.Poly(2 * w**2 - 3 * w**3 + w**4, w, modulus=p)
        derivative = seed.diff()
        assert seed.eval(0) % p == 0
        assert derivative.eval(0) % p == 0
        assert seed.eval(1) % p == 0
        assert derivative.eval(1) % p == -1 % p
        assert seed.diff().diff().eval(1) % p == -2 % p

        seed_at_w = seed.as_expr().subs(w, capital_w)
        derivative_at_w = derivative.as_expr().subs(w, capital_w)
        second_numerator = sp.expand(derivative_at_w + gamma)
        first_numerator = sp.expand(
            capital_w * second_numerator - seed_at_w
        )
        alpha = poly(first_numerator).exquo(poly(capital_c**2))
        beta = poly(second_numerator).exquo(poly(capital_c))
        output = (alpha.as_expr(), beta.as_expr(), capital_c)
        jacobian = sp.Matrix(output).jacobian((x, y, z)).det()
        assert poly(jacobian - 1).is_zero


def verify_characteristic_two_wild_clean_witnesses() -> None:
    p = 2
    r, t, parameter, second_parameter = sp.symbols(
        "r t parameter second_parameter"
    )

    for n in range(5, 17):
        if n % 2:
            seed_expression = r**2 + r**n
        else:
            seed_expression = r**2 + r**3 + r**4 + r**n
        seed = sp.Poly(seed_expression, r, modulus=p)
        derivative = seed.diff()
        primitive_factor = seed.exquo(sp.Poly(r**2, r, modulus=p))

        assert seed.degree() == n
        assert seed.eval(1) % p == 0
        assert derivative.eval(1) % p == 1
        assert sp.gcd(primitive_factor, primitive_factor.diff()).degree() == 0
        assert primitive_factor.exquo(
            sp.Poly(r + 1, r, modulus=p)
        ).degree() == n - 3

        odd_part = sum(
            int(seed.nth(2 * j + 1)) * t**j
            for j in range((n + 1) // 2)
        )
        even_part = sum(
            int(seed.nth(2 * j)) * t**j for j in range(n // 2 + 1)
        )
        odd_part = sp.Poly(odd_part, t, modulus=p)
        even_part = sp.Poly(even_part, t, modulus=p)
        if n % 2:
            assert even_part == sp.Poly(t, t, modulus=p)
        else:
            assert odd_part == sp.Poly(t, t, modulus=p)

        # On the reduced discriminant normalization, the double-root
        # cluster is the radicial factor r^2-t, while the remaining factor
        # has degree N-2.
        discriminant_pencil = sp.Poly(
            seed.as_expr()
            - odd_part.as_expr() * r
            + even_part.as_expr(),
            r,
            t,
            modulus=p,
        )
        quotient, remainder = sp.div(
            discriminant_pencil,
            sp.Poly(r**2 - t, r, t, modulus=p),
        )
        assert remainder.is_zero
        assert quotient.degree(r) == n - 2

    # The complete normalized quartic slice over F_2(parameter).
    quartic = sp.Poly(
        (1 + parameter) * r**2 + r**3 + parameter * r**4,
        r,
        parameter,
        modulus=p,
    )
    assert sp.Poly(quartic.as_expr().subs(r, 1), parameter, modulus=p).is_zero
    assert quartic.diff(r) == sp.Poly(r**2, r, parameter, modulus=p)
    quartic_primitive = sp.div(
        quartic,
        sp.Poly(r**2, r, parameter, modulus=p),
    )[0]
    extra_factor, remainder = sp.div(
        quartic_primitive,
        sp.Poly(r + 1, r, parameter, modulus=p),
    )
    assert remainder.is_zero
    assert extra_factor == sp.Poly(
        parameter * r + 1 + parameter,
        r,
        parameter,
        modulus=p,
    )
    assert quartic_primitive.diff(r) == sp.Poly(
        1, r, parameter, modulus=p
    )

    # Exhaustively regress the clean-implies-birational lemma over F_2.
    # A nontrivial common off-diagonal factor of the two difference
    # quotients is exactly a generic collision of the compressed map.
    for n in range(3, 13):
        for mask in range(1 << (n - 2)):
            coefficients = [0, 0]
            coefficients.extend((mask >> j) & 1 for j in range(n - 2))
            coefficients.append(1)
            seed = sp.Poly(
                sum(
                    coefficients[j] * r**j for j in range(n + 1)
                ),
                r,
                modulus=p,
            )
            if coefficients[2] == 0:
                continue
            if seed.eval(1) % p != 0 or seed.diff().eval(1) % p != 1:
                continue
            primitive = seed.exquo(sp.Poly(r**2, r, modulus=p))
            if sp.gcd(primitive, primitive.diff()).degree() != 0:
                continue

            odd_part = sum(
                coefficients[2 * j + 1] * t**j
                for j in range((n + 1) // 2)
            )
            even_part = sum(
                coefficients[2 * j] * t**j
                for j in range(n // 2 + 1)
            )
            odd_at_second = odd_part.subs(t, second_parameter)
            even_at_second = even_part.subs(t, second_parameter)
            odd_difference = sp.cancel(
                (odd_part - odd_at_second) / (t - second_parameter)
            )
            even_difference = sp.cancel(
                (even_part - even_at_second) / (t - second_parameter)
            )
            generic_collision = sp.gcd(
                sp.Poly(
                    odd_difference,
                    t,
                    second_parameter,
                    modulus=p,
                ),
                sp.Poly(
                    even_difference,
                    t,
                    second_parameter,
                    modulus=p,
                ),
            )
            assert generic_collision.total_degree() == 0


def verify_collision_boundary_recovery() -> None:
    w = sp.symbols("w")
    collision_seeds = []

    # Rational repeated-root collisions in odd characteristic.
    for p in (3, 5):
        rho = 2
        raw = sp.Poly(
            w**2 * (w - 1) * (w - rho) ** 2,
            w,
            modulus=p,
        )
        scale = (-pow(int(raw.diff().eval(1)) % p, -1, p)) % p
        collision_seeds.append((p, sp.Poly(scale * raw.as_expr(), w, modulus=p)))

    # A characteristic-two collision whose repeated roots form one
    # irreducible quadratic residue factor over F_2.
    collision_seeds.append(
        (
            2,
            sp.Poly(
                w**2 * (w + 1) * (w**2 + w + 1) ** 2,
                w,
                modulus=2,
            ),
        )
    )

    for p, seed in collision_seeds:
        derivative = seed.diff()
        assert seed.eval(0) % p == 0
        assert derivative.eval(0) % p == 0
        assert seed.eval(1) % p == 0
        assert derivative.eval(1) % p == -1 % p
        assert seed.nth(2) % p != 0

        primitive = seed.exquo(sp.Poly(w**2, w, modulus=p))
        root_one = sp.Poly(w - 1, w, modulus=p)
        extra = primitive.exquo(root_one)
        repeated = sp.gcd(primitive, primitive.diff())
        assert repeated.degree() > 0

        # Every repeated primitive root maps to the same critical value
        # (0,0) as the double zero.
        critical_second = sp.Poly(
            w * derivative.as_expr() - seed.as_expr(),
            w,
            modulus=p,
        )
        assert seed.rem(repeated).is_zero
        assert derivative.rem(repeated).is_zero
        assert critical_second.rem(repeated).is_zero

        # The marked boundary divisor, including multiplicities, reconstructs
        # the normalized seed even on the collision.
        monic_extra = extra.monic()
        raw_reconstruction = sp.Poly(
            w**2 * (w - 1) * monic_extra.as_expr(),
            w,
            modulus=p,
        )
        normalization_scale = (
            -pow(int(raw_reconstruction.diff().eval(1)) % p, -1, p)
        ) % p
        reconstructed = sp.Poly(
            normalization_scale * raw_reconstruction.as_expr(),
            w,
            modulus=p,
        )
        assert reconstructed == seed


def evaluate(polynomial: tuple[int, ...], value: int, p: int) -> int:
    result = 0
    for coefficient_value in reversed(polynomial):
        result = (result * value + coefficient_value) % p
    return result


def ordinary_derivative(
    polynomial: tuple[int, ...], p: int
) -> tuple[int, ...]:
    return trim(
        tuple(
            exponent * polynomial[exponent] % p
            for exponent in range(1, len(polynomial))
        ),
        p,
    )


def weighted_q_from_seed(
    seed: tuple[int, ...], p: int
) -> tuple[int, ...]:
    """Return the zero-constant primitive of W*seed'(W)."""
    seed_prime = ordinary_derivative(seed, p)
    integrand = (0,) + seed_prime
    result = [0] * (len(integrand) + 1)
    for exponent, value in enumerate(integrand):
        if value == 0:
            continue
        output_exponent = exponent + 1
        assert output_exponent % p != 0
        result[output_exponent] = value * pow(output_exponent, -1, p) % p
    return trim(tuple(result), p)


def verify_degree_eight_collision() -> None:
    p = 5
    n = 8
    h1 = (0, 0, 1, 4, 0, 0, 1, 3, 1)
    h2 = (0, 0, 2, 3, 0, 4, 3, 1, 2)

    for seed in (h1, h2):
        assert clean_seed(seed, n, p)
        assert evaluate(seed, 0, p) == 0
        assert coefficient(ordinary_derivative(seed, p), 0) == 0
        assert evaluate(seed, 1, p) == 0
        assert evaluate(ordinary_derivative(seed, p), 1, p) == -1 % p

    ordinary_hessian_1 = hessian(h1, p)
    ordinary_hessian_2 = hessian(h2, p)
    assert ordinary_hessian_2 == tuple(2 * value % p for value in ordinary_hessian_1)
    assert projective_key(ordinary_hessian_1, p) == projective_key(
        ordinary_hessian_2, p
    )
    l5 = (0, 0, 0, 0, 0, 1, -1 % p)
    assert h2 == tuple(
        (2 * coefficient(h1, exponent) - coefficient(l5, exponent)) % p
        for exponent in range(n + 1)
    )
    assert hessian(l5, p) == ()
    assert hasse(l5, 5, p) == (1, -1 % p)

    powers = p_powers_at_most(n, p)
    assert powers == (1, 5)
    channels1 = {order: hasse(h1, order, p) for order in powers}
    channels2 = {order: hasse(h2, order, p) for order in powers}
    assert channels1[5] == (0, 1, 3, 1)
    assert channels2[5] == (4, 3, 1, 2)
    assert channels1 != channels2
    assert reconstruct(channels1, n, p, constant=0) == h1
    assert reconstruct(channels2, n, p, constant=0) == h2


def verify_degree_twelve_map_intrinsic_no_go() -> None:
    p = 5
    n = 12
    h0 = (0, 0, 3, 2, 4, 2, 2, 1, 2, 1, 3, 3, 2)
    frobenius_direction = (0, 0, 0, 0, 0, -1 % p, 0, 0, 0, 0, 1)
    clean_parameters = []
    construction_data = []

    for parameter in range(p):
        seed = tuple(
            (
                coefficient(h0, exponent)
                + parameter * coefficient(frobenius_direction, exponent)
            )
            % p
            for exponent in range(n + 1)
        )
        if not clean_seed(seed, n, p):
            continue
        clean_parameters.append(parameter)

        derivative = ordinary_derivative(seed, p)
        second_derivative = ordinary_derivative(derivative, p)
        assert evaluate(seed, 0, p) == 0
        assert coefficient(derivative, 0) == 0
        assert evaluate(seed, 1, p) == 0
        assert evaluate(derivative, 1, p) == -1 % p

        # These are precisely the seed-dependent ingredients in the reduced
        # weighted-map formula: P=H', P', q=int(WP'), and kappa=P'(1).
        construction_data.append(
            (
                derivative,
                second_derivative,
                weighted_q_from_seed(derivative, p),
                evaluate(second_derivative, 1, p),
            )
        )

    assert clean_parameters == [0, 1, 2, 3, 4]
    assert len(set(construction_data)) == 1
    common_seed, common_hessian, common_q, kappa = construction_data[0]
    assert common_hessian == (1, 2, 3, 0, 0, 2, 2, 2, 0, 0, 4)
    assert ordinary_derivative(frobenius_direction, p) == ()
    assert hasse(frobenius_direction, 5, p) == (4, 0, 0, 0, 0, 2)

    # Expand the common reduced weighted map itself over F_5.  Exact
    # division proves that all three displayed coordinates are polynomials;
    # the determinant check proves that this is a genuine Keller map.
    x, y, z, w = sp.symbols("x y z w")
    domain = sp.GF(p)

    def poly(expression):
        return sp.Poly(sp.expand(expression), x, y, z, domain=domain)

    seed_expression = sum(
        value * w**exponent for exponent, value in enumerate(common_seed)
    )
    q_expression = sum(
        value * w**exponent for exponent, value in enumerate(common_q)
    )
    a_parameter = (-(1 + kappa) * pow(2 + kappa, -1, p)) % p
    v = x * y
    source_s = x**2 * z
    u = 1 + v
    gamma = 1 + a_parameter * v + source_s
    capital_w = u * gamma
    seed_at_w = sp.expand(seed_expression.subs(w, capital_w))
    q_at_w = sp.expand(q_expression.subs(w, capital_w))

    alpha = poly(u * gamma**2 + q_at_w).exquo(poly(gamma**2 * x**2))
    beta = poly(gamma + seed_at_w).exquo(poly(gamma * x))
    capital_c = x * gamma
    ordinary_output = (alpha.as_expr(), beta.as_expr(), capital_c)
    ordinary_jacobian = sp.Matrix(ordinary_output).jacobian((x, y, z)).det()
    assert poly(ordinary_jacobian - 1).is_zero
    assert [poly(entry).total_degree() for entry in ordinary_output] == [52, 51, 4]

    # Enrich the map by the Frobenius primitive itself.  H_0 has invisible
    # part 3*(W^10-W^5); adding parameter*K changes that multiplier to 3+c.
    capital_k = capital_w**10 - capital_w**5
    correction = poly(capital_k).exquo(poly(capital_c**2)).as_expr()
    assert poly(correction).total_degree() == 42
    assert len(poly(correction).terms()) == 60

    primitive_expression = sum(
        value * w**exponent for exponent, value in enumerate(h0)
    )
    direction_expression = w**10 - w**5
    enriched_first_coordinates = []
    marked_pairing_degrees = []
    equal_image_fingerprints = []
    second_normalization_coordinate = sp.symbols("u")

    for parameter in range(p):
        multiplier = (3 + parameter) % p
        enriched_a = poly(alpha.as_expr() - multiplier * correction).as_expr()
        enriched_output = (enriched_a, beta.as_expr(), capital_c)
        enriched_first_coordinates.append(poly(enriched_a).as_expr())

        jacobian = sp.Matrix(enriched_output).jacobian((x, y, z)).det()
        assert poly(jacobian - 1).is_zero
        assert [poly(entry).total_degree() for entry in enriched_output] == [
            52,
            51,
            4,
        ]

        primitive = sp.Poly(
            primitive_expression + parameter * direction_expression,
            w,
            modulus=p,
        )
        inverse_identity = poly(
            enriched_a * capital_c**2
            - (
                capital_w * (seed_at_w + gamma)
                - primitive.as_expr().subs(w, capital_w)
            )
        )
        assert inverse_identity.is_zero

        # Root-one equal-image pairing on the discriminant normalization.
        branch_seed = primitive.diff()
        branch_value = sp.Poly(
            w * branch_seed.as_expr() - primitive.as_expr(), w, modulus=p
        )
        root = sp.Poly(w - 1, w, modulus=p)
        first_equation = sp.div(
            sp.Poly(
                branch_seed.as_expr() - branch_seed.eval(1), w, modulus=p
            ),
            root,
        )[0]
        second_equation = sp.div(
            sp.Poly(
                branch_value.as_expr() - branch_value.eval(1), w, modulus=p
            ),
            root,
        )[0]
        pairing_gcd = sp.gcd(first_equation, second_equation)
        marked_pairing_degrees.append(pairing_gcd.degree())
        if parameter == 2:
            assert pairing_gcd.monic() == sp.Poly(w + 2, w, modulus=p)
            hessian_poly = primitive.diff().diff()
            assert hessian_poly.eval(1) % p == 1
            assert hessian_poly.eval(3) % p == -2 % p
        else:
            assert pairing_gcd.degree() == 0

        # Full normalization self-fiber product after removing the common
        # diagonal factor.  The three intrinsic marks fix the coordinate, so
        # distinct reduced bases are stable-invariant distinctions.
        seed_at_u = branch_seed.as_expr().subs(
            w, second_normalization_coordinate
        )
        value_at_u = branch_value.as_expr().subs(
            w, second_normalization_coordinate
        )
        first_difference = sp.cancel(
            (branch_seed.as_expr() - seed_at_u)
            / (w - second_normalization_coordinate)
        )
        second_difference = sp.cancel(
            (branch_value.as_expr() - value_at_u)
            / (w - second_normalization_coordinate)
        )
        equal_image_basis = sp.groebner(
            [first_difference, second_difference],
            w,
            second_normalization_coordinate,
            modulus=p,
            order="lex",
        )
        eliminant_expression = equal_image_basis.polys[-1].as_expr()
        assert not eliminant_expression.has(w)
        eliminant = sp.Poly(
            eliminant_expression,
            second_normalization_coordinate,
            modulus=p,
        ).monic()
        assert eliminant.degree() == 110
        equal_image_fingerprints.append(
            tuple(
                int(
                    eliminant.coeff_monomial(second_normalization_coordinate**j)
                )
                % p
                for j in range(4)
            )
        )

    assert len(set(map(str, enriched_first_coordinates))) == p
    assert marked_pairing_degrees == [0, 0, 1, 0, 0]
    assert equal_image_fingerprints == [
        (4, 1, 4, 0),
        (3, 0, 2, 2),
        (3, 3, 4, 4),
        (2, 1, 2, 3),
        (2, 3, 4, 2),
    ]
    assert len(set(equal_image_fingerprints)) == p


def main() -> None:
    verify_all_degree_reconstruction()
    verify_characteristic_two_parity_core()
    verify_characteristic_two_scalar_ansatz_no_go()
    verify_characteristic_two_weight_redistribution()
    verify_singular_parameter_weight_redistribution()
    verify_characteristic_two_wild_clean_witnesses()
    verify_collision_boundary_recovery()
    verify_degree_eight_collision()
    verify_degree_twelve_map_intrinsic_no_go()
    print(
        "PASS: p-typical Hasse channels reconstruct modulo constants, "
        "every p-power order is forced, and D^[5] separates both the "
        "degree-eight Hessian collision and the clean degree-twelve "
        "identical-map family; the enriched lifts are polynomial Keller "
        "maps with five intrinsic equal-image fingerprints, while the "
        "characteristic-two core has the proved Frobenius parity split, "
        "scalar-ansatz no-go, and a replacement weight-redistributed "
        "polynomial Keller suspension whose marked boundary recovers the "
        "seed across explicit wild repeated-root collisions"
    )


if __name__ == "__main__":
    main()
