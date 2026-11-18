#!/usr/bin/env python3
"""Finite audits for the maximal-order primitive norm class-error argument.

This program checks finite group saturation, mixed prime/prime-square CRT
identities, exact rational Selberg minimization, and finite random marking.
It does not prove the imported genus theorem, prime-class Mertens estimate,
or infinite convolution theorem. All examples and proof domains are fixed.
"""
from __future__ import annotations


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import gcd, prod
from pathlib import Path

F = Fraction


def subsets(items):
    return [tuple(c) for k in range(len(items) + 1) for c in combinations(items, k)]


def add(a, b, moduli):
    return tuple((x + y) % m for x, y, m in zip(a, b, moduli))


def neg(a, moduli):
    return tuple((-x) % m for x, m in zip(a, moduli))


def group_sum(items, moduli):
    result = (0,) * len(moduli)
    for item in items:
        result = add(result, item, moduli)
    return result


def repeated_order(a, moduli):
    zero = (0,) * len(moduli)
    current = a
    order = 1
    while current != zero:
        current = add(current, a, moduli)
        order += 1
        assert order <= prod(moduli)
    return order


def reject(label, mutated_statement, controls):
    """A deliberately corrupted mathematical statement must fail."""
    assert not mutated_statement, label
    controls.append(label)


def group_audit():
    examples = [(3,), (4,), (8,), (6,), (2, 2), (4, 2), (3, 3), (6, 2)]
    records = []
    controls = []
    cofactor_cases = 0
    assignment_cases = 0
    for moduli in examples:
        elements = list(product(*(range(m) for m in moduli)))
        zero = (0,) * len(moduli)
        twice_group = {add(x, x, moduli) for x in elements}
        pairs = [x for x in elements if add(x, x, moduli) != zero
                 and x < neg(x, moduli)]
        copies = []
        thresholds = []
        for g in pairs:
            twice_g = add(g, g, moduli)
            m = repeated_order(twice_g, moduli) - 1
            copies.extend([g] * m)
            thresholds.append([list(g), m])
        # Independent brute force: flip every reserved prime separately.
        orientation_sums = set()
        for signs in product((-1, 1), repeat=len(copies)):
            oriented = [g if sign == 1 else neg(g, moduli)
                        for sign, g in zip(signs, copies)]
            orientation_sums.add(group_sum(oriented, moduli))
            assignment_cases += 1
        base = group_sum(copies, moduli)
        assert orientation_sums == {add(base, t, moduli) for t in twice_group}
        for c in elements:
            actual = {add(c, x, moduli) for x in orientation_sums}
            genus_coset = {add(add(c, base, moduli), t, moduli) for t in twice_group}
            assert actual == genus_coset
            assert (zero in actual) == (zero in genus_coset)
            cofactor_cases += 1
        records.append({
            'cyclic_factors': list(moduli),
            'group_order': len(elements),
            'twice_group_order': len(twice_group),
            'nonambiguous_class_pairs': len(pairs),
            'thresholds': thresholds,
            'reserved_prime_count': len(copies),
            'independent_sign_assignments': 2 ** len(copies),
            'distinct_signed_classes': len(orientation_sums),
            'all_cofactor_classes_checked': len(elements),
            'saturation': True,
        })
    # A single C8 block needs ord(2g)-1=3 copies. Two copies miss a class.
    moduli = (8,)
    short = {(sum(signs) % 8,) for signs in product((-1, 1), repeat=2)}
    expected = {(0,), (2,), (4,), (6,)}
    reject('insufficient C8 single-block threshold', short == expected, controls)
    # Enough signs cover the correct genus coset, not all of the group.
    reject('principal conclusion without principal-genus hypothesis',
           (0,) in {(1,), (3,)}, controls)
    return {
        'examples': records,
        'cofactor_cases': cofactor_cases,
        'independent_sign_assignments': assignment_cases,
        'negative_controls': controls,
        'scope': 'Full elements and all reserved-prime signs in the eight stated finite groups.',
    }


def inert_bad(q, k, n):
    # c=n+2k, v=-k parameterizes c+2v=n. Both slopes are q-units.
    return (n + 2 * k) % q == 0 or k % q == 0


def exact_one_bad(q, k, n):
    c = n + 2 * k
    return c % q == 0 and c % (q * q) != 0


def crt_pair(x, modulus_x, y, modulus_y):
    assert gcd(modulus_x, modulus_y) == 1
    return (x + modulus_x * ((y - x) * pow(modulus_x, -1, modulus_y) % modulus_y)) % (modulus_x * modulus_y)


def root_audit():
    controls = []
    single_cases = 0
    mixed_cases = 0
    direct_classes = 0
    for q in (3, 5, 7, 11, 13):
        for valuation in range(4):
            n = q ** valuation
            inert = {k for k in range(q) if inert_bad(q, k, n)}
            marked = {k for k in range(q * q) if exact_one_bad(q, k, n)}
            assert len(inert) == 2 - (valuation > 0)
            assert len(marked) == q - 1
            assert F(len(marked), q * q) == F(q - 1, q * q)
            # Both zero and deeper-than-one valuations are allowed by the mark.
            zero_root = (-n * pow(2, -1, q * q)) % (q * q)
            assert not exact_one_bad(q, zero_root, n)
            if valuation >= 2:
                assert not exact_one_bad(q, 0, n)
            single_cases += 1
    for p in (3, 5, 7):
        for q in (11, 13):
            modulus = p * q * q
            for vp in range(4):
                for vq in range(4):
                    n = p ** vp * q ** vq
                    roots_p = {k for k in range(p) if inert_bad(p, k, n)}
                    roots_q = {k for k in range(q * q) if exact_one_bad(q, k, n)}
                    built = {crt_pair(a, p, b, q * q) for a in roots_p for b in roots_q}
                    direct = {k for k in range(modulus)
                              if inert_bad(p, k, n) and exact_one_bad(q, k, n)}
                    assert built == direct
                    assert len(direct) == (2 - (vp > 0)) * (q - 1)
                    # The union size challenges using the product for a union.
                    union = sum(inert_bad(p, k, n) or exact_one_bad(q, k, n)
                                for k in range(modulus))
                    assert union == len(roots_p) * q * q + len(roots_q) * p - len(direct)
                    mixed_cases += 1
                    direct_classes += modulus
    reject('forgot inert target-root merger',
           sum(inert_bad(3, k, 3) for k in range(3)) == 2, controls)
    reject('replaced exact-one density by divisibility density',
           F(sum(exact_one_bad(5, k, 1) for k in range(25)), 25) == F(1, 5), controls)
    reject('forbade a higher prime power at a marked prime',
           exact_one_bad(5, 0, 125), controls)
    e0 = 1 / (1 - F(2, 3))
    en = 1 / (1 - F(1, 3))
    assert en == e0 / 2
    reject('inverted the inert target correction J', en == e0 * 2, controls)
    return {
        'single_prime_cases': single_cases,
        'mixed_CRT_cases': mixed_cases,
        'direct_CRT_residues_examined': direct_classes,
        'target_valuations_checked': [0, 1, 2, 3],
        'parameterization': 'c=n+2k, v=-k; c+2v=n',
        'inert_primes_in_mixed_cases': [3, 5, 7],
        'marked_primes_in_mixed_cases': [11, 13],
        'marked_density': '(q-1)/q^2; the q^2-divisible residue is allowed',
        'negative_controls': controls,
        'scope': 'Primes have formal disjoint sieve roles; no actual field prime classes are assigned here.',
    }


def linear_solve(matrix, rhs):
    """Independent exact Gaussian elimination; positive pivots challenge convexity."""
    n = len(rhs)
    work = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    pivots = []
    for i in range(n):
        pivot = work[i][i]
        assert pivot > 0
        pivots.append(pivot)
        for j in range(i + 1, n):
            factor = work[j][i] / pivot
            for k in range(i, n + 1):
                work[j][k] -= factor * work[i][k]
    answer = [F(0)] * n
    for i in range(n - 1, -1, -1):
        answer[i] = (work[i][n] - sum(work[i][j] * answer[j]
                                    for j in range(i + 1, n))) / work[i][i]
    return answer, pivots


def selberg_and_marking_audit():
    inert = (3, 7)
    omega = (5, 11)
    z = 25
    controls = []
    selberg_cases = 0
    majorant_cases = 0
    perturbations = 0
    target_hash = hashlib.sha256()
    for n in range(2, 102):
        survivor_counts = []
        for mark in subsets(omega):
            primes = inert + mark
            nu = {q: 2 - (n % q == 0) if q in inert else q - 1 for q in primes}
            mod = {q: q if q in inert else q * q for q in primes}
            f = {q: F(nu[q], mod[q]) for q in primes}
            w = {q: f[q] / (1 - f[q]) for q in primes}
            indices = sorted((prod(t), t) for t in subsets(primes) if prod(t) <= z)
            wh = {h: prod((w[q] for q in t), start=F(1)) for h, t in indices}
            G = sum(wh.values(), F(0))
            lam = {}
            for h, t in indices:
                gh = sum((wh[k] for k, _ in indices if k <= z // h and gcd(k, h) == 1), F(0))
                lam[h] = (-1) ** len(t) * prod((1 / (1 - f[q]) for q in t), start=F(1)) * gh / G
            assert lam[1] == 1
            qmat = [[prod((f[q] for q in set(t) | set(u)), start=F(1))
                     for _, u in indices] for _, t in indices]
            # Solve the constrained rational minimization without using w or G.
            solved, pivots = linear_solve([row[1:] for row in qmat[1:]],
                                          [-row[0] for row in qmat[1:]])
            assert all(p > 0 for p in pivots)
            weights = [lam[h] for h, _ in indices]
            assert weights == [F(1)] + solved
            quad = sum(weights[i] * weights[j] * qmat[i][j]
                       for i in range(len(indices)) for j in range(len(indices)))
            assert quad == 1 / G
            endpoint = sum(abs(weights[i] * weights[j])
                           * prod(nu[q] for q in set(t) | set(u))
                           for i, (_, t) in enumerate(indices)
                           for j, (_, u) in enumerate(indices))

            def event(q, k):
                return (k % q == 0 or (n - k) % q == 0) if q in inert else (k % q == 0 and k % (q * q) != 0)

            direct_count = sum(all(not event(q, k) for q in primes) for k in range(1, n))
            direct_square = sum(sum(lam[h] for h, t in indices if all(event(q, k) for q in t)) ** 2
                                for k in range(1, n))
            assert direct_count <= direct_square
            assert abs(direct_square - (n - 1) * quad) <= endpoint
            # A nontrivial feasible perturbation must strictly increase Q.
            perturbed = weights[:]
            perturbed[1] += F(1, 7)
            qbad = sum(perturbed[i] * perturbed[j] * qmat[i][j]
                       for i in range(len(indices)) for j in range(len(indices)))
            assert qbad > quad
            perturbations += 1
            survivor_counts.append(direct_count)
            target_hash.update(f'{n}:{mark}:{quad}:{direct_count}:{direct_square}\n'.encode())
            selberg_cases += 1
        expectation = F(sum(survivor_counts), len(subsets(omega)))
        for m in (1, 2, 3):
            deficient = sum(all(k % q and (n - k) % q for q in inert)
                            and sum(k % q == 0 and k % (q * q) != 0 for q in omega) <= m - 1
                            for k in range(1, n))
            assert deficient <= 2 ** (m - 1) * expectation
            majorant_cases += 1
    # An exact independent average over the four mark subsets.
    f = {q: F(q - 1, q * q) for q in omega}
    avg = sum((prod((1 - f[q] for q in mark), start=F(1))
               for mark in subsets(omega)), F(0)) / 4
    assert avg == prod((1 - f[q] / 2 for q in omega), start=F(1))
    # Independent singleton: a=5,b=1 has one exact-one Omega divisor and
    # no inert divisor. It survives exactly half the four mark subsets.
    point_a, point_b = 5, 1
    assert all(point_a % q and point_b % q for q in inert)
    point_probability = F(sum(all(not (point_a % q == 0 and point_a % (q * q) != 0)
                                     for q in mark) for mark in subsets(omega)), 4)
    assert point_probability == F(1, 2)
    reject('omitted the bounded-deficiency probability factor',
           F(1) <= point_probability, controls)
    reject('dropped the prime-square term in Euler expectation',
           avg == prod((1 - F(1, 2 * q) for q in omega), start=F(1)), controls)
    reject('corrupted feasible Selberg weights are still optimal', qbad == quad, controls)
    return {
        'target_range': [2, 101],
        'inert_prime_roles': list(inert),
        'marked_class_prime_roles': list(omega),
        'squarefree_index_cutoff': z,
        'mark_sets_per_target': 4,
        'exact_Selberg_cases': selberg_cases,
        'independent_positive_pivot_optimization_cases': selberg_cases,
        'direct_interval_square_and_endpoint_cases': selberg_cases,
        'strict_feasible_weight_perturbations': perturbations,
        'deficiency_thresholds': [1, 2, 3],
        'probabilistic_majorant_cases': majorant_cases,
        'mark_probability': '1/2',
        'exact_inverse_Euler_expectation': str(avg),
        'exact_case_data_sha256': target_hash.hexdigest(),
        'negative_controls': controls,
    }


def distinct_prime_factors(n):
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    answer = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                answer = -answer
        a, n = n, a
        if a % 4 == n % 4 == 3:
            answer = -answer
        a %= n
    return answer if n == 1 else 0


def valuation_unit(a, p):
    assert a != 0
    e = 0
    while a % p == 0:
        a //= p
        e += 1
    return e, a


def hilbert_integer(a, d, p):
    assert a > 0 and d < 0
    e, u = valuation_unit(a, p)
    beta, v = valuation_unit(d, p)
    if p == 2:
        eps = lambda z: ((z - 1) // 2) % 2
        omega = lambda z: ((z * z - 1) // 8) % 2
        exponent = eps(u) * eps(v) + e * omega(v) + beta * omega(u)
        return -1 if exponent % 2 else 1
    answer = -1 if (e * beta * ((p - 1) // 2)) % 2 else 1
    if beta % 2:
        answer *= jacobi(u, p)
    if e % 2:
        answer *= jacobi(v, p)
    return answer


def splitting_sign(d, p):
    if d % p == 0:
        return 0
    if p == 2:
        return 1 if d % 8 == 1 else -1
    return jacobi(d, p)


def principal_form(d):
    return (1, 1, (1 - d) // 4) if d % 4 == 1 else (1, 0, -d // 4)


def norm_value(form, x, y):
    a, b, c = form
    return a * x * x + b * x * y + c * y * y


def direct_primitive_local_image(d, p, base_exp, target_exp):
    """Complete coordinates, compressed only by an exact finite Taylor lift.

    Every target vector is (x,y)+p^base_exp*(u,v). The quadratic remainder
    vanishes because 2*base_exp>=target_exp. For each base vector, the two
    linear derivatives generate precisely the gcd subgroup in the lift
    index. No Hilbert or target-residue criterion is used here.
    """
    assert 2 * base_exp >= target_exp
    base = p ** base_exp
    target = p ** target_exp
    lift = target // base
    a, b, c = principal_form(d)
    cosets = set()
    for x in range(base):
        for y in range(base):
            if x % p == y % p == 0:
                continue
            derivative_gcd = gcd(gcd(2 * a * x + b * y, b * x + 2 * c * y), lift)
            period = base * derivative_gcd
            cosets.add((norm_value((a, b, c), x, y) % period, period))
    flags = bytearray(target)
    for residue, period in cosets:
        for value in range(residue, target, period):
            flags[value] = 1
    if target <= 400:
        brute = bytearray(target)
        for x in range(target):
            for y in range(target):
                if x % p or y % p:
                    brute[norm_value((a, b, c), x, y) % target] = 1
        assert brute == flags
    return flags


def cyclic_selfsum(flags):
    modulus = len(flags)
    bits = sum(1 << n for n, present in enumerate(flags) if present)
    full_mask = (1 << modulus) - 1
    covered = 0
    for a, present in enumerate(flags):
        if present:
            covered |= ((bits << a) | (bits >> (modulus - a))) & full_mask
    return bytearray((covered >> n) & 1 for n in range(modulus))


def local_hilbert_flag(d, p, value):
    sign = splitting_sign(d, p)
    if sign == 1:
        return True
    if sign == -1:
        return value % p != 0
    if value == 0:
        return False
    e, _ = valuation_unit(value, p)
    return e <= 1 and hilbert_integer(value, d, p) == 1


def finite_global_primitive_flags(d, bound):
    from math import isqrt
    form = principal_form(d)
    _, b, _ = form
    flags = bytearray(bound + 1)
    ymax = isqrt((4 * bound) // (-d))
    for y in range(-ymax, ymax + 1):
        h = isqrt(4 * bound + d * y * y)
        lo = -((h + b * y) // 2)
        hi = (h - b * y) // 2
        for x in range(lo, hi + 1):
            if gcd(x, y) == 1:
                n = norm_value(form, x, y)
                if 0 < n <= bound:
                    flags[n] = 1
    return flags


def positive_shadow_flag(d, n, bad_primes):
    assert n > 0
    remaining = n
    for p in distinct_prime_factors(n):
        e, _ = valuation_unit(n, p)
        sign = splitting_sign(d, p)
        if sign == -1 or (sign == 0 and e > 1):
            return False
        while remaining % p == 0:
            remaining //= p
    assert remaining == 1
    return all(hilbert_integer(n, d, p) == 1 for p in bad_primes)


def histogram(a_values, b_values, a_coeff, b_coeff, modulus):
    counts = [0] * modulus
    for a in a_values:
        for b in b_values:
            counts[(a_coeff * a + b_coeff * b) % modulus] += 1
    return counts


def actual_field_audit():
    fields = (-3, -4, -7, -8, -11, -15, -20, -23, -24, -40,
              -56, -84, -120, -164, -260, -231)
    records = []
    controls = []
    frozen_state_witness = None
    prefix_label_witness = None
    missing_lift_witness = None
    split2_witness = None
    nonprincipal_ramification_witness = None
    for d in fields:
        odd_ramified = [p for p in distinct_prime_factors(-d) if p != 2]
        bad_primes = [2] + odd_ramified
        r_product = prod(odd_ramified)
        level = 8 * r_product
        prefix_support = 2 * r_product
        target_period = 16 * r_product * r_product
        prefixes = sorted(prod(s) for s in subsets([2] + odd_ramified)
                          if not (splitting_sign(d, 2) == -1 and 2 in s))
        units = [a for a in range(1, level) if gcd(a, level) == 1]
        state = lambda a: tuple(hilbert_integer(a, d, p) for p in bad_primes)
        theta = {r: state(r) for r in prefixes}
        h_units = [a for a in units if jacobi(d, a) == 1]
        assert len(h_units) * 2 == len(units)
        for a in units:
            assert prod(state(a)) == jacobi(d, a)
        for r in prefixes:
            assert prod(theta[r]) == 1
        residue_sets = {r: [a for a in units if state(a) == theta[r]] for r in prefixes}
        assert all(residue_sets[r] and set(residue_sets[r]) <= set(h_units) for r in prefixes)
        representatives = {}
        for q in h_units:
            representatives.setdefault(state(q), q)
        representatives = sorted(set(representatives.values()) | set(h_units[1:2]))

        local_records = []
        coarse_target_flags = []
        for p in bad_primes:
            base_exp, target_exp, coarse = (4, 7, 16) if p == 2 else (2, 3, p * p)
            flags = direct_primitive_local_image(d, p, base_exp, target_exp)
            expected = bytearray(local_hilbert_flag(d, p, n) for n in range(len(flags)))
            assert flags == expected
            sums = cyclic_selfsum(flags)
            for n in range(coarse):
                assert all(sums[n + k * coarse] == sums[n]
                           for k in range(len(sums) // coarse))
            coarse_target_flags.append((coarse, sums[:coarse]))
            local_records.append({'prime': p, 'coordinate_base_modulus': p ** base_exp,
                                  'full_local_modulus': p ** target_exp,
                                  'direct_vector_image_count': sum(flags),
                                  'direct_selfsum_target_count': sum(sums),
                                  'stabilized_target_modulus': coarse,
                                  'direct_image_sha256': hashlib.sha256(flags).hexdigest(),
                                  'direct_selfsum_sha256': hashlib.sha256(sums).hexdigest(),
                                  'small_modulus_full_vector_crosscheck': len(flags) <= 400})
        direct_target_mask = bytearray(all(flags[n % modulus] for modulus, flags in coarse_target_flags)
                                      for n in range(target_period))
        profile_mask = bytearray(target_period)
        profile_count = 0
        full_inflation_profiles = 0
        affine_lift_fiber_profiles = 0
        inflation_target_checks = 0
        transport_cases = 0
        bad_coefficient_profiles = 0
        minimum_gamma = None
        for r in prefixes:
            for s in prefixes:
                common = gcd(r, s)
                a_coeff, b_coeff = r // common, s // common
                assert gcd(a_coeff, b_coeff) == 1
                if a_coeff > 1 and b_coeff > 1:
                    bad_coefficient_profiles += 1
                ar, bs = residue_sets[r], residue_sets[s]
                counts = histogram(ar, bs, a_coeff, b_coeff, level)
                for t, count in enumerate(counts):
                    if count:
                        gamma = F(level * count, len(h_units) ** 2)
                        minimum_gamma = gamma if minimum_gamma is None else min(minimum_gamma, gamma)
                        for n in range(common * t, target_period, common * level):
                            profile_mask[n] = 1
                # Independent finite inflation: every coordinate has B unit lifts.
                ell = b_coeff * level
                lifted_ar = [a + k * level for a in ar for k in range(b_coeff)]
                lifted_bs = [b + k * level for b in bs for k in range(b_coeff)]
                assert all(gcd(a, ell) == 1 and state(a) == theta[r] for a in lifted_ar)
                assert all(gcd(b, ell) == 1 and state(b) == theta[s] for b in lifted_bs)
                if b_coeff <= 16:
                    inflated = histogram(lifted_ar, lifted_bs, a_coeff, b_coeff, ell)
                    assert all(inflated[t] == b_coeff * counts[t % level] for t in range(ell))
                    inflation_target_checks += ell
                    full_inflation_profiles += 1
                else:
                    # Complete affine lift-index fiber, valid for any base root:
                    # shifting its carry merely rotates this whole histogram.
                    fiber = [0] * b_coeff
                    for i in range(b_coeff):
                        for j in range(b_coeff):
                            fiber[(a_coeff * i + b_coeff * j) % b_coeff] += 1
                    assert fiber == [b_coeff] * b_coeff
                    affine_lift_fiber_profiles += 1
                    inflation_target_checks += b_coeff
                if b_coeff > 1 and any(counts) and missing_lift_witness is None:
                    t = next(t for t, count in enumerate(counts) if count)
                    missing_lift_witness = [d, r, s, t, counts[t], b_coeff * counts[t]]
                for q in representatives:
                    qstate = state(q)
                    alpha = tuple(x * y for x, y in zip(theta[r], qstate))
                    beta = tuple(x * y for x, y in zip(theta[s], qstate))
                    new_ar = [a for a in units if state(a) == alpha]
                    new_bs = [b for b in units if state(b) == beta]
                    assert {(q * a) % level for a in new_ar} == set(ar)
                    assert {(q * b) % level for b in new_bs} == set(bs)
                    directive = histogram(new_ar, bs, a_coeff * q, b_coeff, level)
                    assert directive == counts
                    common_counts = histogram(new_ar, new_bs, a_coeff, b_coeff, level)
                    inverse = pow(q, -1, level)
                    assert all(common_counts[(t * inverse) % level] == counts[t]
                               for t in range(level))
                    if frozen_state_witness is None:
                        wrong = histogram(ar, bs, a_coeff * q, b_coeff, level)
                        if wrong != counts:
                            t = next(t for t in range(level) if wrong[t] != counts[t])
                            frozen_state_witness = [d, r, s, q, t, counts[t], wrong[t]]
                    transport_cases += 1
                if prefix_label_witness is None:
                    wrong = histogram(residue_sets[a_coeff], residue_sets[b_coeff],
                                      a_coeff, b_coeff, level)
                    if wrong != counts:
                        t = next(t for t in range(level) if wrong[t] != counts[t])
                        prefix_label_witness = [d, r, s, common, t, counts[t], wrong[t]]
                profile_count += 1
        assert profile_mask == direct_target_mask
        bound = 4096
        primitive = finite_global_primitive_flags(d, bound)
        shadow = bytearray(bound + 1)
        for n in range(1, bound + 1):
            shadow[n] = positive_shadow_flag(d, n, bad_primes)
            assert not primitive[n] or shadow[n]
        hidden = [n for n in range(1, bound + 1) if shadow[n] and not primitive[n]]
        if splitting_sign(d, 2) == 1 and split2_witness is None:
            n = next(n for n in range(4, bound + 1, 4) if shadow[n])
            split2_witness = [d, n]
        if d == -20:
            assert primitive[1] and not primitive[2]
            nonprincipal_ramification_witness = [d, 1, 2]
        records.append({'discriminant': d, 'principal_form': list(principal_form(d)),
                        'level_L': level, 'target_period_16R2': target_period,
                        'prefixes': prefixes, 'genus_state_count_on_H': len({state(a) for a in h_units}),
                        'profiles': profile_count, 'both_bad_coefficients_profiles': bad_coefficient_profiles,
                        'full_inflated_pair_profiles': full_inflation_profiles,
                        'complete_affine_lift_fiber_profiles': affine_lift_fiber_profiles,
                        'inflation_target_or_fiber_checks': inflation_target_checks,
                        'state_transport_cases': transport_cases,
                        'minimum_positive_Gamma': str(minimum_gamma),
                        'primitive_local_target_residue_count': sum(direct_target_mask),
                        'profile_target_mask_sha256': hashlib.sha256(profile_mask).hexdigest(),
                        'local_coordinate_images': local_records,
                        'global_coordinate_bound': bound, 'primitive_global_values': sum(primitive),
                        'individual_Hilbert_shadow_values': sum(shadow),
                        'hidden_individual_class_failure_count_in_bound': len(hidden),
                        'first_hidden_failures': hidden[:10],
                        'global_coordinate_flags_sha256': hashlib.sha256(primitive).hexdigest(),
                        'individual_shadow_flags_sha256': hashlib.sha256(shadow).hexdigest()})
    assert frozen_state_witness is not None and prefix_label_witness is not None
    assert missing_lift_witness is not None and split2_witness is not None
    assert nonprincipal_ramification_witness is not None
    reject('froze the genus state under a directive',
           frozen_state_witness[-2] == frozen_state_witness[-1], controls)
    reject('replaced original prefix labels after cancelling their gcd',
           prefix_label_witness[-2] == prefix_label_witness[-1], controls)
    reject('omitted the B-fold inflation fiber',
           missing_lift_witness[-2] == missing_lift_witness[-1], controls)
    reject('truncated all individual split-2 valuations', split2_witness[1] % 4 != 0, controls)
    reject('assumed a ramified prime has a principal generator', nonprincipal_ramification_witness is None, controls)
    return {'fields': records, 'field_count': len(records),
            'negative_controls': controls,
            'corruption_witnesses': {'frozen_state': frozen_state_witness,
                                    'cancelled_prefix_labels': prefix_label_witness,
                                    'missing_inflation_factor': missing_lift_witness,
                                    'unbounded_individual_split2': split2_witness,
                                    'nonprincipal_ramification': nonprincipal_ramification_witness},
            'scope': ('Complete stated local residue images and finite profile masks, plus bounded global '
                      'coordinate-to-Hilbert-shadow inclusion. Higher local digits use exact finite Taylor '
                      'lifting, independently challenged by full vector enumeration at moduli<=400. '
                      'No ideal class labels or infinite prime-distribution assertion are inferred.')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Write the deterministic JSON receipt to this path.')
    args = parser.parse_args()
    report = {
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'proof_domains_fixed': True,
        'scope': ('Finite independent audits of group saturation, exact prime-square CRT, '
                  'rational Selberg optimization and finite random marking. The infinite '
                  'arithmetic/genus identification, prime-class Mertens theorem and '
                  'unweighted class-error bound are separate written/imported mathematics. '
                  'No external review or formal verification is claimed.'),
        'group_saturation': group_audit(),
        'mixed_local_roots': root_audit(),
        'Selberg_and_marking': selberg_and_marking_audit(),
        'actual_fields': actual_field_audit(),
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded, encoding='utf-8')
        print(json.dumps({'output': str(args.output),
                          'checker_sha256': report['checker_sha256'],
                          'receipt_sha256': hashlib.sha256(encoded.encode()).hexdigest(),
                          'group_examples': len(report['group_saturation']['examples']),
                          'Selberg_cases': report['Selberg_and_marking']['exact_Selberg_cases'],
                          'random_majorant_cases': report['Selberg_and_marking']['probabilistic_majorant_cases']},
                         sort_keys=True))
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
