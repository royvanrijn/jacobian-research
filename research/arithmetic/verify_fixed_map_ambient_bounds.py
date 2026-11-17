#!/usr/bin/env python3
"""Exact controls for FIXED_MAP_AMBIENT_COUNT_BOUNDS.md; not a sieve proof."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import gcd
from pathlib import Path
import platform

import sympy as sp


def compose(g, h):
    return tuple(g[h[i]] for i in range(5))


def inverse(g):
    return tuple(g.index(i) for i in range(5))


def even(g):
    return sum(g[i] > g[j] for i in range(5) for j in range(i + 1, 5)) % 2 == 0


def cycle_type(g):
    seen = set()
    lengths = []
    for i in range(5):
        if i in seen:
            continue
        n, j = 0, i
        while j not in seen:
            seen.add(j)
            n += 1
            j = g[j]
        lengths.append(n)
    return tuple(sorted(lengths))


def mobius(n):
    sign, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
        while n % p == 0:
            n //= p
        p += 1
    return -sign if n > 1 else sign


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    checks = []

    def check(name, condition, **details):
        assert condition, name
        checks.append({'name': name, 'passed': True, **details})

    group = list(permutations(range(5)))
    block = frozenset((0, 1))
    K = {g for g in group if frozenset(g[i] for i in block) == block}
    H = {g for g in K if even(g)}
    C = {g for g in group if cycle_type(g) in {(1, 1, 1, 1, 1), (1, 1, 3), (1, 2, 2)}}
    union = {compose(compose(g, h), inverse(g)) for g in group for h in H}
    check('S5_order', len(group) == 120)
    check('quadratic_block_stabilizer', len(K) == 12)
    check('resolvent_block_stabilizer', len(H) == 6)
    check('all_H_permutations_even', all(even(h) for h in H))
    check('H_is_subgroup', all(compose(g, h) in H for g in H for h in H))
    check('two_block_H_orbits', {h[0] for h in H} == {0, 1}
          and {h[2] for h in H} == {2, 3, 4})
    check('quadratic_factor_cover_degree', len(group) // len(K) == 10)
    check('sign_refined_cover_degree', len(group) // len(H) == 20)
    check('conjugacy_union_exact', union == C and len(C) == 36)
    class_sizes = {str(t): sum(cycle_type(g) == t for g in group)
                   for t in sorted({cycle_type(g) for g in group})}
    check('allowed_class_sizes', class_sizes[str((1, 1, 1, 1, 1))] == 1
          and class_sizes[str((1, 1, 3))] == 20
          and class_sizes[str((1, 2, 2))] == 15)
    root_C = {g for g in group if any(g[i] == i for i in range(5))}
    check('fixed_point_subset', len(root_C) == 76 and C <= root_C and len(root_C) < 120)
    check('Zywina_Hasse_exponents', Fraction(36, 120) == Fraction(3, 10)
          and 3 - 1 + Fraction(3, 10) == Fraction(23, 10)
          and 3 + Fraction(3, 10) == Fraction(33, 10))
    pairs = set(combinations(range(5), 2))
    orbit = {tuple(sorted((g[0], g[1]))) for g in group}
    check('quadratic_factor_action_transitive', orbit == pairs and len(orbit) == 10)

    S, p, b, c, a, v, w = sp.symbols('S p b c a v w')
    E = (p**5*S**5 - sp.Rational(3, 2)*p**4*S**4
         + sp.Rational(3, 2)*p*S**3
         + (-sp.Rational(5, 4)*p - sp.Rational(9, 32)*b)*S**2
         + sp.Rational(9, 16)*S - sp.Rational(9, 32)*c)
    f = E/p**5
    check('generic_monic_degree_five', sp.Poly(f, S).degree() == 5
          and sp.Poly(f, S).LC() == 1)
    check('quadratic_tilt_nonzero_linear_term',
          sp.diff(f, S).subs(S, 0) == sp.Rational(9, 16)/p**5)
    q, h = S**2 + v, S**3 + a*S + w
    product_disc = sp.discriminant(q*h, S)
    factor_disc = sp.discriminant(q, S)*sp.discriminant(h, S)*sp.resultant(q, h, S)**2
    check('symbolic_factor_discriminant_identity', sp.factor(product_disc-factor_disc) == 0)
    line = (S**3-a)*(S**2+S+1)
    check('original_line_square_discriminant',
          sp.factor(sp.discriminant(line, S)-(9*a*(a-1)**2)**2) == 0)
    check('scalar_discriminant_power_is_square', (2*5-2)*5 == 40)

    plane_counts = []
    for B in range(1, 17):
        direct = sum(gcd(gcd(q0, abs(m)), abs(n)) == 1
                     for q0 in range(1, B+1)
                     for m in range(-B, B+1) for n in range(-B, B+1))
        sieved = sum(mobius(d)*(B//d)*(2*(B//d)+1)**2 for d in range(1, B+1))
        check(f'primitive_plane_count_B{B}', direct == sieved)
        if B in (1, 4, 8, 16):
            plane_counts.append({'B': B, 'rational_plane_targets': direct})

    here = Path(__file__).resolve().parent
    proof = here/'FIXED_MAP_AMBIENT_COUNT_BOUNDS.md'
    primary_paths = [here/'sources/zywina-quantitative-hit.pdf',
                     here/'sources/buggenhout-et-al-thin-projective-2506.13471v2.pdf']
    result = {
        'status': 'passed', 'python': platform.python_version(), 'sympy': sp.__version__,
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'proof_sha256': sha256(proof.read_bytes()).hexdigest(),
        'primary_sources': {p.name: sha256(p.read_bytes()).hexdigest() for p in primary_paths},
        'check_count': len(checks), 'checks': checks, 'S5_class_sizes': class_sizes,
        'plane_counts': plane_counts,
        'scope': 'Exact symbolic and complete finite permutation controls, and primitive plane '
                 'counts through B16. No proof replay of the imported monodromy, sieve, or '
                 'projective thin-set theorems; no asymptotic inferred from finite data.',
        'formal_verification': False, 'independent_whole_proof_replay': False,
        'external_review': False,
    }
    if args.write_receipt:
        (here/'fixed_map_ambient_bounds_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'check_count': len(checks),
                      'plane_counts': plane_counts}, indent=2))


if __name__ == '__main__':
    main()
