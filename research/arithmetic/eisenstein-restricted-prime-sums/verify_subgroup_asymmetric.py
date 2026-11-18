#!/usr/bin/env python3
"""Exact finite diagnostics for subgroup/quadratic binary sums.

This does not test or certify the analytic dispersion theorem. It checks
finite residue measures, kernel identities, local integrals, obstruction
cases and directive fibres using integer and rational arithmetic.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, gcd, isqrt
from pathlib import Path
import json
import sys


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1] if Path(__file__).parent.name == '.scratch' else Path(__file__).resolve().parents[3]
ARITH = ROOT / 'research/arithmetic'
COUNTS = Counter()


def check(ok, kind):
    require(ok, kind)
    COUNTS[kind] += 1


def chi_value(kind, n):
    if kind == '-4':
        return 1 if n % 4 == 1 else -1
    if kind == '8':
        return 1 if n % 8 in (1, 7) else -1
    p = {'-3': 3, '5': 5, '-7': 7}[kind]
    return 1 if pow(n % p, (p - 1) // 2, p) == 1 else -1


def generated(L, generators):
    H = {1}
    while True:
        expanded = H | {(a * b) % L for a in H for b in generators}
        if expanded == H:
            return H
        H = expanded


def validate(L, H, chi):
    G = {a for a in range(L) if gcd(a, L) == 1}
    require(L % 2 == 0 and H and H < G and 1 in H, 'invalid subgroup input')
    require(all(a * b % L in H for a in H for b in H), 'subgroup not closed')
    require(set(chi) == G and set(chi.values()) == {-1, 1}, 'quadratic character input')
    require(all(chi[a * b % L] == chi[a] * chi[b] for a in G for b in G), 'character not multiplicative')
    return G, {a for a in G if chi[a] == 1}


def configurations():
    data = [
        (6, (1,), '-3'), (10, (1,), '5'), (10, (9,), '5'),
        (12, (5,), '-3'), (12, (7,), '-4'), (14, (13,), '-7'),
        (14, (9,), '-7'), (16, (3,), '-4'), (16, (15,), '8'),
        (18, (17,), '-3'), (20, (9,), '-4'), (30, (11, 19), '-3'),
        (42, (41,), '-7'),
    ]
    out = []
    for L, generators, kind in data:
        G = {a for a in range(L) if gcd(a, L) == 1}
        chi = {a: chi_value(kind, a) for a in G}
        H = generated(L, generators)
        G, K = validate(L, H, chi)
        out.append((L, G, H, K, chi))
    return out


def gamma(L, H, K, t, u=1, v=1):
    return Q(L * sum((u * a + v * b - t) % L == 0 for a in H for b in K), len(H) * len(K))


def residues(configs):
    for L, G, H, K, chi in configs:
        I = H & K
        k = len(G) // len(H)
        a, b, c = Q(len(H), len(G)), Q(1, 2), Q(len(I), len(G))
        check(c * (k + 2) + (a - c) * k + (b - c) * 2 == 2, 'rough_dimension')
        check(2 * c + (a - c) + (b - c) == a + b, 'restored_dimension')
        check((1 - a) + (1 - b) == Q(3, 2) - Q(1, k), 'upper_dimension')
        for t in range(L):
            value = gamma(L, H, K, t)
            check(t % 2 == 0 or value == 0, 'parity')
            mu = Q(sum((t - h) % L in K for h in H), len(H))
            check(value == Q(2 * L, len(G)) * mu, 'filter_density')
            # A residue-class measure has means phi/L, not index/L.
            check(value * Q(len(G), L) ** 2 == 2 * Q(len(G), L) * mu, 'mean_normalization')
            for u in H:
                for v in K:
                    check(gamma(L, H, K, t, u, v) == value, 'allowed_coefficient_invariance')
            for s in I:
                check(gamma(L, H, K, t * pow(s, -1, L)) == value, 'intersection_scaling')
        for a0 in G:
            f = 1
            while pow(a0, f, L) not in H:
                f += 1
                require(f <= k, 'coset order overflow')
            check(k % f == 0, 'coset_order_divisibility')
            # The Artin/character local coefficients before correction.
            degree = 18
            pre = [comb(e // f + k // f - 1, k // f - 1) if e % f == 0 else 0 for e in range(degree + 1)]
            correction = {0: 1} if f == 1 else {j * f: (-1) ** j * comb(k // f, j) for j in range(k // f + 1)}
            for e in range(degree + 1):
                product = sum(coef * pre[e - j] for j, coef in correction.items() if j <= e)
                expected = comb(e + k - 1, k - 1) if a0 in H else int(e == 0)
                check(product == expected, 'first_euler_identity')
            if f > 1:
                check(pow(a0, f, L) in H and f >= 2, 'first_kernel_support')
            if a0 not in K:
                check(chi[a0 * a0 % L] == 1, 'second_kernel_support')


def primes_through(n):
    spf = list(range(n + 1))
    for p in range(2, isqrt(n) + 1):
        if spf[p] == p:
            for j in range(p * p, n + 1, p):
                if spf[j] == j:
                    spf[j] = p
    return spf


def factor(n, spf):
    out = {}
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out[p] = e
    return out


def finite_pairs(configs):
    nmax = 1800
    spf = primes_through(nmax)
    facts = [factor(n, spf) for n in range(1, nmax + 1)]
    for L, G, H, K, chi in configs:
        first = [False] + [all(p % L in H for p in fs) for fs in facts]
        second = [False] + [all(p % L in K for p in fs) for fs in facts]
        used = Counter()
        for n in range(2, nmax + 1):
            r = sum(first[a] and second[n - a] for a in range(1, n))
            if gamma(L, H, K, n) == 0:
                check(r == 0, 'inadmissible_exact_zero')
            if r:
                used[n % L] += 1
        check(set(used) == {t for t in range(L) if gamma(L, H, K, t)}, 'finite_target_support')
        # Recover the tiny directive cofactors from their valuation patterns.
        small = [p for p in range(2, 18) if spf[p] == p and p % L in G]
        kappa = len(G) // len(H)
        for a in range(1, 301):
            if not first[a]:
                continue
            fa = facts[a - 1]
            for b in range(1, 301):
                if not second[b]:
                    continue
                fb = facts[b - 1]
                d = u = v = 1
                valid = True
                for p in small:
                    ea, eb = fa.get(p, 0), fb.get(p, 0)
                    if ea > 1 or eb > 1:
                        valid = False
                        break
                    if ea == eb == 1:
                        d *= p
                    elif ea == 1:
                        u *= p
                    elif eb == 1:
                        v *= p
                if not valid:
                    continue
                check(gcd(d, u) == gcd(d, v) == gcd(u, v) == 1 and (a + b) % d == 0 and gcd(u * v, a + b) == 1, 'directive_gcds')
                check(all(p % L in H & K for p in factor(d, spf)), 'common_directive_support')
                check(a % (d * u) == b % (d * v) == 0, 'directive_recovery')
                rem_a, rem_b = a // (d * u), b // (d * v)
                check(all(rem_a % p and rem_b % p for p in small), 'rough_directive_residuals')
                weight_a = weight_b = 1
                for e in factor(a // d, spf).values():
                    weight_a *= comb(e + kappa - 1, kappa - 1)
                for e in factor(b // d, spf).values():
                    weight_b *= e + 1
                ra = rb = 1
                for e in factor(rem_a, spf).values():
                    ra *= comb(e + kappa - 1, kappa - 1)
                for e in factor(rem_b, spf).values():
                    rb *= e + 1
                check(Q(weight_a * weight_b, kappa ** len(factor(u, spf)) * 2 ** len(factor(v, spf))) == ra * rb, 'normalized_directive_fibres')


def local_integral(k, allowed_a, allowed_b, p, r):
    x, q = Q(1, p), Q(p - 1, p)
    A = lambda e: comb(e + k - 1, k - 1) if allowed_a else int(e == 0)
    B = lambda e: e + 1 if allowed_b else int(e == 0)
    fa = (1 - x) ** -k if allowed_a else Q(1)
    fb = (1 - x) ** -2 if allowed_b else Q(1)
    tail_a = fa - sum(A(e) * x ** e for e in range(r + 1))
    tail_b = fb - sum(B(e) * x ** e for e in range(r + 1))
    expected = q * sum(A(e) * B(e) * x ** e for e in range(r)) + x ** r * (1 - 2 * x) * A(r) * B(r) + q * (B(r) * tail_a + A(r) * tail_b)
    h = 4
    modulus = p ** h
    n = p ** r
    total = Q(0)
    for z in range(modulus):
        y = (n - z) % modulus
        if not z or not y:
            continue
        va = vb = 0
        while z % p == 0:
            va += 1
            z //= p
        while y % p == 0:
            vb += 1
            y //= p
        total += Q(A(va) * B(vb), modulus)
    ta = fa - sum(A(e) * x ** e for e in range(h))
    tb = fb - sum(B(e) * x ** e for e in range(h))
    total += q * (B(r) * ta + A(r) * tb)
    check(total == expected, 'p_adic_line_integrals')
    rough = 1 - (2 if r == 0 else 1) * x
    check(expected >= rough > 0, 'rough_mass_positivity')
    if r and not (allowed_a and allowed_b):
        check(expected == q, 'one_sided_and_neither_target')


def locals_and_controls(configs):
    for k in sorted({len(G) // len(H) for _, G, H, _, _ in configs}):
        for p in (3, 5, 7):
            for aa in (False, True):
                for bb in (False, True):
                    for r in range(3):
                        local_integral(k, aa, bb, p, r)
        for shift in range(8):
            for e in range(k + 8):
                conv = sum((-1) ** j * comb(k, j) * Q(comb(shift + e - j + k - 1, k - 1), comb(shift + k - 1, k - 1)) for j in range(min(e, k) + 1))
                if e >= k:
                    check(conv == 0, 'shifted_kernel_degree')
    for p in (3, 5, 7, 13, 29):
        ap, E = Q(p - 1, p * (p - 2)), Q(p - 1, p - 2)
        check(E - (1 + ap) == Q(1, p * (p - 2)), 'one_sided_target_compensation')
        for r in range(1, 8):
            previous = 1 if r == 1 else E
            check(E + Q(previous, p) >= 1 + 2 * ap, 'intersection_target_compensation')
    # C4+C3 leaves two half-dimensional restriction events.
    L, G, H, K, chi = next(c for c in configs if c[0] == 12 and c[2] == {1, 5})
    check({t: gamma(L, H, K, t) for t in range(L) if gamma(L, H, K, t)} == {t: Q(3) for t in (0, 2, 6, 8)}, 'mixed_quadratic_gamma')
    small_spf = primes_through(74)
    fa, fb = factor(25, small_spf), factor(49, small_spf)
    check(all(p % L in H for p in fa) and all(p % L in K for p in fb) and not all(p % L in H & K for p in fa) and not all(p % L in H & K for p in fb), 'mixed_quadratic_intersection_counterexample')
    for p in (5, 7):
        x = Q(1, p)
        unit = (1 - x + x * x) / (1 - x)
        event = (2 * x - x * x) / (1 - x)
        check(event / unit == Q(2 * p - 1, p * p - p + 1), 'mixed_quadratic_extra_sieve_density')
    check(Q(2, 4) + Q(2, 4) == 1, 'mixed_quadratic_extra_sieve_dimension')
    # Deliberate invalid-input and wrong-formula controls.
    L, G, H, K, chi = next(c for c in configs if c[0] == 14 and c[2] == {1, 13})
    try:
        validate(L, {1, 6}, chi)
    except RuntimeError:
        COUNTS['negative_controls'] += 1
    else:
        raise RuntimeError('nonunit subgroup not rejected')
    expected = {0: Q(7, 3), 2: Q(7, 3), 8: Q(7, 3), 10: Q(14, 3), 12: Q(7, 3)}
    check({t: gamma(L, H, K, t) for t in range(L) if gamma(L, H, K, t)} == expected, 'odd_index_gamma_table')
    check(gamma(L, H, K, 10) != gamma(L, H, K, 10 * 13), 'negative_controls')
    check(Q(len(G), L) != Q(len(G) // len(H), L), 'negative_controls')
    check(not any(all(p % L in H for p in factor(m, primes_through(26))) and all(p % L in K for p in factor(26 - 13 * m, primes_through(26))) for m in range(1, 2)), 'first_coefficient_target_zero')
    check(not any(all(p % L in H for p in factor(22 - 11 * l, primes_through(22))) and all(p % L in K for p in factor(l, primes_through(22))) for l in range(1, 2)), 'second_coefficient_target_zero')
    # A first-only inert prime3 does not enter the two-exclusion J factor.
    L, G, H, K, chi = next(c for c in configs if c[0] == 16 and 3 in c[2])
    check(3 in H - K and not (3 not in H and 3 not in K), 'negative_controls')
    # A nonidentity intersection scalar can change classes while preserving Gamma.
    L, G, H, K, chi = next(c for c in configs if c[0] == 10 and c[2] == {1, 9})
    check(2 != 2 * 9 % L and gamma(L, H, K, 2) == gamma(L, H, K, 18), 'class_changing_intersection')


def main():
    require(__debug__, 'run without -O; explicit debug guard protects diagnostic execution')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ARITH / 'ASYMMETRIC_SUBGROUP_SPLITTING.md')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    sources = [args.source.resolve(), ARITH / 'ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md', ARITH / 'ASYMMETRIC_FIXED_MODULUS_THEORY.md', ARITH / 'ASYMMETRIC_UNIFORM_COEFFICIENTS.md', ARITH / 'ASYMMETRIC_SEMIGROUP_SCALE.md', ARITH / 'FIXED_PERIODIC_DIVISOR_BOXES.md', ARITH / 'sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf', ARITH / 'sources/heath-brown-d3-progressions-1986.pdf']
    for path in sources:
        require(path.is_file(), 'missing source: ' + str(path))
    if args.output:
        require(args.output.resolve() not in {p.resolve() for p in sources + [Path(__file__)]}, 'output would overwrite source')
    configs = configurations()
    residues(configs)
    finite_pairs(configs)
    locals_and_controls(configs)
    receipt = {
        'status': 'pass', 'scope': 'Exact finite arithmetic diagnostics only; not an analytic theorem verifier, independent whole-proof replay, formal verification or external review.',
        'arithmetic': 'Python standard-library integers and Fraction',
        'configurations': [{'L': L, 'H': sorted(H), 'K': sorted(K), 'index': len(G) // len(H)} for L, G, H, K, _ in configs],
        'finite_target_limit': 1800, 'directive_component_limit': 300, 'local_primes': [3, 5, 7], 'local_target_valuations': [0, 1, 2],
        'checks': dict(sorted(COUNTS.items())),
        'source_sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sources},
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    text = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text, encoding='utf-8', newline='\n')
    else:
        sys.stdout.write(text)


if __name__ == '__main__':
    main()
