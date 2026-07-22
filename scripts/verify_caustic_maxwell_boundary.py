#!/usr/bin/env python3
"""Exact checks for the caustic and Maxwell root-boundary classes."""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def free_orders(k: int) -> tuple[int, int, int]:
    caustic = (k - 1) * (k - 2)
    maxwell = (k - 1) * (k - 2) * (k - 3) // 2
    ll_discriminant = k * (k - 1) * (k - 2)
    return caustic, maxwell, ll_discriminant


def zero_orders(k: int) -> tuple[int, int, int]:
    caustic = k * (k - 1)
    maxwell = k * (k - 1) ** 2 // 2
    ll_discriminant = k * (k - 1) * (k + 2)
    return caustic, maxwell, ll_discriminant


def infinity_poles(n: int, k: int) -> tuple[int, int, int]:
    q = n - k
    degree = n + 2
    caustic = q * (q - 1 + 2 * k)
    ll_discriminant = q * (
        degree * (q - 1 + 2 * k) + k * (k - 1)
    )
    maxwell = q * ((n - 1) * (q - 1 + 2 * k) + k * (k - 1)) // 2
    return caustic, maxwell, ll_discriminant


def check_invariant_keel_presentation(n: int) -> None:
    """Check the orbit rank and unique relation from labelled Keel relations."""

    marked_count = n + 2
    marks = set(range(marked_count))
    divisors = [
        frozenset(part)
        for size in range(2, marked_count - 1)
        for part in combinations(range(marked_count), size)
        if 0 in part
    ]
    divisor_index = {part: index for index, part in enumerate(divisors)}

    def separation_vector(a: int, b: int, c: int, d: int) -> list[int]:
        answer = []
        for part in divisors:
            complement = marks - set(part)
            separated = (
                a in part and b in part and c in complement and d in complement
            ) or (
                a in complement and b in complement and c in part and d in part
            )
            answer.append(int(separated))
        return answer

    relations = []
    for a, b, c, d in combinations(range(marked_count), 4):
        ab_cd = separation_vector(a, b, c, d)
        ac_bd = separation_vector(a, c, b, d)
        ad_bc = separation_vector(a, d, b, c)
        relations.append([x - y for x, y in zip(ab_cd, ac_bd)])
        relations.append([x - y for x, y in zip(ab_cd, ad_bc)])
    keel = sp.Matrix(relations)

    def canonical(part: set[int]) -> frozenset[int]:
        return frozenset(part if 0 in part else marks - part)

    orbit_vectors: list[list[int]] = []
    for k in range(2, n + 1):
        vector = [0] * len(divisors)
        for block in combinations(range(2, marked_count), k):
            vector[divisor_index[canonical(set(block))]] += 1
        orbit_vectors.append(vector)
    for k in range(1, n):
        vector = [0] * len(divisors)
        for block in combinations(range(2, marked_count), k):
            vector[divisor_index[canonical({0, *block})]] += 1
        orbit_vectors.append(vector)

    orbit_matrix = sp.Matrix(orbit_vectors)
    kernel_of_relations = sp.Matrix.hstack(*keel.nullspace())
    quotient_relation_matrix = orbit_matrix * kernel_of_relations
    assert quotient_relation_matrix.rank() == 2 * n - 3

    coefficients = [
        -k * (k - 1) for k in range(2, n + 1)
    ] + [k * (n - k) for k in range(1, n)]
    relation = sp.Matrix([coefficients]).T
    assert quotient_relation_matrix.T * relation == sp.zeros(
        quotient_relation_matrix.cols, 1
    )


for root_count in range(2, 7):
    check_invariant_keel_presentation(root_count)


for n in range(2, 14):
    degree = n + 2
    caustic_degree = n * (n - 1)
    maxwell_degree = n * (n - 1) ** 2 // 2
    ll_degree = degree * n * (n - 1)
    assert 3 * caustic_degree + 2 * maxwell_degree == ll_degree

    for k in range(2, n + 1):
        c_order, m_order, d_order = free_orders(k)
        assert 3 * c_order + 2 * m_order == d_order

    for k in range(1, n):
        q = n - k
        c_zero, m_zero, d_zero = zero_orders(k)
        c_pole, m_pole, d_pole = infinity_poles(n, k)
        assert 3 * c_zero + 2 * m_zero == d_zero
        assert 3 * c_pole + 2 * m_pole == d_pole

        # Average over the k choices of normalized root on the zero side and
        # the q choices on the infinity side.
        c_boundary = (k * c_pole - q * c_zero) // n
        m_boundary = sp.Rational(k * m_pole - q * m_zero, n)
        d_boundary = (k * d_pole - q * d_zero) // n
        assert c_boundary == k * q
        assert m_boundary == sp.Rational(k * q * (n + k - 2), 2)
        assert d_boundary == k * q * (n + k + 1)
        assert 3 * c_boundary + 2 * m_boundary == d_boundary

    # Eliminating Delta_(0,n-1) with the unique invariant Keel relation gives
    # the announced genuine-basis coefficients.
    for k in range(2, n + 1):
        redundant_c_free = -(k - 1) * (k - 2)
        basis_c_free = redundant_c_free + k * (k - 1)
        assert basis_c_free == 2 * (k - 1)

        redundant_m_free = -(k - 1) * (k - 2) * (k - 3) // 2
        basis_m_free = sp.Rational(redundant_m_free) + sp.Rational(
            (2 * n - 3) * k * (k - 1), 2
        )
        expected = sp.Rational(k - 1, 2) * (
            (2 * n - 3) * k - (k - 2) * (k - 3)
        )
        assert basis_m_free == expected

    for k in range(1, n - 1):
        redundant = sp.Rational(k * (n - k) * (n + k - 2), 2)
        eliminated = redundant - sp.Rational(
            (2 * n - 3) * k * (n - k), 2
        )
        assert eliminated == -sp.Rational(k * (n - k) * (n - k - 1), 2)

    # The positive kappa presentation of C and the collision-readable one
    # differ by the invariant Keel relation.
    for k in range(2, n + 1):
        positive_free = sp.Rational(2 * (k - 1) * (n + 1 - k), n + 1)
        collision_free = -(k - 1) * (k - 2)
        relation_multiple = sp.Rational(n - 1, n + 1) * k * (k - 1)
        assert positive_free == collision_free + relation_multiple
    for k in range(1, n):
        positive_zero = sp.Rational(2 * k * (n - k), n + 1)
        collision_zero = k * (n - k)
        relation_multiple = -sp.Rational(n - 1, n + 1) * k * (n - k)
        assert positive_zero == collision_zero + relation_multiple


# Low-degree elimination checks.  For H=W^2 prod(W-a_i), remove the fixed
# critical point W=0 by using K=H'/W.  The critical-value discriminant must
# factor as caustic^3 * Maxwell^2.
w, value = sp.symbols("w value")
for n in (2, 3):
    roots = sp.symbols(f"a1:{n + 1}")
    seed = w**2
    for root in roots:
        seed *= w - root
    critical = sp.cancel(sp.diff(seed, w) / w)
    caustic = sp.factor(sp.discriminant(critical, w))
    critical_values = sp.resultant(critical, value - seed, w)
    ll_discriminant = sp.factor(sp.discriminant(critical_values, value))

    factors = sp.factor_list(ll_discriminant)[1]
    normalized_caustic = sp.Poly(caustic, *roots).monic()
    caustic_factors = [
        (factor, exponent)
        for factor, exponent in factors
        if sp.Poly(factor, *roots).monic() == normalized_caustic
    ]
    assert len(caustic_factors) == 1
    assert caustic_factors[0][1] == 3
    maxwell_factors = [
        (factor, exponent)
        for factor, exponent in factors
        if sp.Poly(factor, *roots).monic() != normalized_caustic
    ]
    assert len(maxwell_factors) == 1
    maxwell, maxwell_exponent = maxwell_factors[0]
    assert maxwell_exponent == 2

    assert sp.Poly(caustic, *roots).total_degree() == n * (n - 1)
    assert sp.Poly(maxwell, *roots).total_degree() == n * (n - 1) ** 2 // 2

    scale = sp.symbols("scale")
    scaled_c = sp.expand(caustic.subs(dict(zip(roots, (scale * r for r in roots)))))
    scaled_m = sp.expand(maxwell.subs(dict(zip(roots, (scale * r for r in roots)))))
    assert sp.expand(scaled_c - scale ** (n * (n - 1)) * caustic) == 0
    assert sp.expand(
        scaled_m - scale ** (n * (n - 1) ** 2 // 2) * maxwell
    ) == 0


# Quartic sanity check: C=9a^2-14ab+9b^2 and M=a+b.
a, b = sp.symbols("a b")
quartic = w**2 * (w - a) * (w - b)
quartic_critical = sp.cancel(sp.diff(quartic, w) / w)
quartic_caustic = sp.factor(sp.discriminant(quartic_critical, w))
assert sp.factor(quartic_caustic - (9 * a**2 - 14 * a * b + 9 * b**2)) == 0

quartic_values = sp.resultant(quartic_critical, value - quartic, w)
quartic_ll = sp.factor(sp.discriminant(quartic_values, value))
quotient = sp.factor(quartic_ll / quartic_caustic**3)
assert sp.Poly(quotient, a, b).total_degree() == 2
assert sp.factor(quotient / (a + b) ** 2).is_number


print("PASS caustic class: 2*kappa_1 and both boundary presentations")
print("PASS Maxwell class: collision valuations and invariant Keel reduction")
print("PASS LL discriminant: caustic^3 * Maxwell^2 in degrees four and five")
