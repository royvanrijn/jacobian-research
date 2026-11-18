#!/usr/bin/env python3
"""Finite audits for the maximal-order primitive norm class-error argument.

This program checks finite group saturation, mixed prime/prime-square CRT
identities, exact rational Selberg minimization, and finite random marking.
It does not prove the imported genus theorem, prime-class Mertens estimate,
or infinite convolution theorem. All examples and proof domains are fixed.
"""
from __future__ import annotations

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
