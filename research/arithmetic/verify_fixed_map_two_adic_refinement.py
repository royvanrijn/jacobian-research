#!/usr/bin/env python3
"""Exact finite controls for FIXED_MAP_TWO_ADIC_REFINEMENT.md; not a Hasse search."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


def valuation2(number):
    if number == 0:
        return None
    number = abs(number)
    answer = 0
    while number % 2 == 0:
        answer += 1
        number //= 2
    return answer


def lower_hull(points):
    hull = []
    for point in points:
        while len(hull) >= 2:
            a, b = hull[-2:]
            cross = (b[0]-a[0])*(point[1]-b[1])-(b[1]-a[1])*(point[0]-b[0])
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def edges(coefficients):
    hull = lower_hull([(i, valuation2(a)) for i, a in enumerate(coefficients) if a])
    return list(zip(hull, hull[1:]))


def coefficients(p, b, c):
    return [-9*c, 18, -(40*p+9*b), 48*p, -48*p**4, 32*p**5]


def shifted_coefficients(p, b, c):
    assert b % 4 == 0
    return [p**5-3*p**4-4*p+9-9*c-9*b//4,
            2*(p-1)*(5*p**4-7*p**3-7*p*p-7*p-9)-9*b,
            8*p*(p-1)*(5*p**3-4*p*p-4*p-4)-9*b,
            16*p*(5*p**4-6*p**3+3), 16*p**4*(5*p-3), 32*p**5]


def residue_after_scale(coeff, r, power):
    scaled = [Fraction(a*2**power, 2**(r*i)) for i, a in enumerate(coeff)]
    assert all(a.denominator == 1 for a in scaled)
    return [int(a) % 2 for a in scaled]


def has_irreducible_block(coeff, degrees):
    return any(b[0]-a[0] in degrees and
               Fraction(b[1]-a[1], b[0]-a[0]).denominator == b[0]-a[0]
               for a, b in edges(coeff))


def symbolic_controls():
    p, b, c, S, T, z, B, C = sp.symbols('p b c S T z B C')
    F = sum(a*S**i for i, a in enumerate([
        -9*c, 18, -(40*p+9*b), 48*p, -48*p**4, 32*p**5]))
    translated = sp.Poly(sp.expand(F.subs(S, T+sp.Rational(1, 2))), T)
    L = 5*p**3-4*p**2-4*p-4
    M = 5*p**4-7*p**3-7*p**2-7*p-9
    expected = [p**5-3*p**4-4*p+9-9*c-sp.Rational(9, 4)*b,
                2*(p-1)*M-9*b, 8*p*(p-1)*L-9*b,
                16*p*(5*p**4-6*p**3+3), 16*p**4*(5*p-3), 32*p**5]
    assert all(sp.expand(translated.nth(i)-value) == 0 for i, value in enumerate(expected))
    assert sp.Poly(L.subs(p, 2*z+1), z, modulus=2).as_expr() == 1
    assert sp.Poly(M.subs(p, 2*z+1), z, modulus=2).as_expr() == 1
    g3_unit = sp.Poly(sp.expand((5*p**4-6*p**3+3).subs(p, 2*z+1)-2), z)
    assert all(int(a) % 4 == 0 for a in g3_unit.all_coeffs())
    parity = sp.Poly(sp.expand(expected[0].subs({p: 2*z+1, b: 8*B, c: 2*C+1})), z, B, C)
    assert all(int(a) % 2 == 0 for a in parity.coeffs())
    return {'translated_coefficient_identities': 6, 'unit_and_parity_identities': 4}


def even_pi_controls():
    counts = Counter()
    odd_units = [-3, -1, 1, 3]
    bs = sorted({0} | {2**f*u for f in range(8) for u in odd_units})
    cs = sorted({0} | {2**h*u for h in range(7) for u in odd_units})
    for e in range(1, 9):
        for unit in odd_units:
            p = 2**e*unit
            for b in bs:
                f = valuation2(b)
                allowed = ((e == 1 and (f == 1 or f is None or f >= 3)) or
                           (e == 2 and (f == 0 or f is None or f >= 2)) or
                           (e == 3 and (f is None or f >= 1)) or e >= 4)
                for c in cs:
                    coeff = coefficients(p, b, c)
                    if not allowed:
                        assert has_irreducible_block(coeff, {3})
                        counts['forbidden_tame_cubic_blocks'] += 1
                    elif f is None or f >= 3:
                        assert residue_after_scale(coeff, 1, 0) == [c % 2, 1, 0, 0, 0, 0]
                        counts['automatic_linear_residues'] += 1
                    elif f == 0:
                        if e == 2:
                            actual = residue_after_scale(coeff, 5, 10)
                            assert actual == [0, 0, 1, 0, 0, 1]
                        else:
                            assert e >= 4
                            actual = residue_after_scale(coeff, e+4, 2*(e+4))
                            assert actual == [0, 0, 1, 1, 0, 0]
                        assert sum(actual) % 2 == 0
                        assert sum(i*actual[i] for i in range(1, 6)) % 2 == 1
                        counts['automatic_odd_b_residues'] += 1
                    elif f == 1 and e == 1:
                        assert residue_after_scale(coeff, 3, 5) == [0, 0, 1, 0, 0, 1]
                        counts['automatic_e_one_f_one'] += 1
                    else:
                        assert (f == 1 and e >= 3) or (f == 2 and e >= 2)
                        assert ((2, f), (3, 4+e)) in edges(coeff)
                        counts['automatic_adjacent_edges'] += 1
    return dict(counts), {'pi_valuations': [1, 8], 'pi_units': odd_units,
                          'b_valuations': [0, 7], 'c_valuations': [0, 6],
                          'b_and_c_include_zero': True}


def odd_even_controls():
    counts = Counter()
    for p in range(-15, 16, 2):
        for b in range(-32, 33):
            f = valuation2(b)
            for c in range(-16, 17, 2):
                coeff = coefficients(p, b, c)
                if f in (0, 1):
                    assert has_irreducible_block(coeff, {3})
                    counts['forbidden_tame_cubic_blocks'] += 1
                else:
                    assert b % 4 == 0
                    assert residue_after_scale(coeff, 1, 0) == [0, 1, (b//4) % 2, 0, 1, 1]
                    counts['automatic_simple_root_zero'] += 1
    return dict(counts)


def odd_odd_controls():
    counts = Counter()
    for p in range(-63, 64, 2):
        for b in range(-128, 129, 8):
            for c in range(-63, 64, 2):
                g = shifted_coefficients(p, b, c)
                a0 = g[0]
                case_one = p % 8 == 5 and b % 16 == 0 and a0 % 8 == 0
                case_two = p % 16 == 1 and b % 32 == 0 and a0 % 32 == 0
                assert not (case_one and case_two)
                if case_one:
                    assert valuation2(g[1]) == 3
                    if a0:
                        assert ((0, valuation2(a0)), (1, 3)) in edges(g)
                    else:
                        assert g[1] != 0
                    counts['case_one_simple_local_roots'] += 1
                elif case_two:
                    assert all(a % 32 == 0 for a in g)
                    assert all((g[i]//32) % 2 == 1 for i in (3, 4, 5))
                    counts['case_two_unresolved_integral_reductions'] += 1
                elif a0 == 0:
                    counts['excluded_rational_half_roots'] += 1
                else:
                    certificates = [bpt[0]-apt[0] for apt, bpt in edges(g)
                                    if bpt[0]-apt[0] in (3, 4, 5) and
                                    Fraction(bpt[1]-apt[1], bpt[0]-apt[0]).denominator == bpt[0]-apt[0]]
                    assert certificates, (p, b, c, g, edges(g))
                    counts[f'excluded_degree_{certificates[0]}_block'] += 1
    assert sum(counts.values()) == 135168
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    symbolic = symbolic_controls()
    even, even_domain = even_pi_controls()
    odd_even = odd_even_controls()
    odd_odd = odd_odd_controls()
    here = Path(__file__).resolve().parent
    note = here/'FIXED_MAP_TWO_ADIC_REFINEMENT.md'
    inputs = [here/'FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md', here/'sources/milne-ant300.pdf']
    report = {'status': 'passed', 'python': platform.python_version(), 'sympy': sp.__version__,
              'command': '.venv/bin/python research/arithmetic/verify_fixed_map_two_adic_refinement.py --write-receipt',
              'proof_sha256': sha256(note.read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in inputs},
              'symbolic': symbolic, 'even_pi_controls': even, 'even_pi_domain': even_domain,
              'odd_pi_even_c_controls': odd_even,
              'odd_pi_even_c_domain': {'pi': 'odd integers -15 through15', 'b': [-32, 32], 'c': 'even -16 through16'},
              'odd_pi_odd_c_controls': odd_odd,
              'odd_pi_odd_c_domain': {'pi': 'odd -63 through63', 'b': 'multiples of8 from-128 through128',
                                     'c': 'odd -63 through63'},
              'scope': ('Exact coefficient and finite residue controls, plus bounded independent Newton '
                        'certificates. They do not prove the unrestricted valuation cases, the local '
                        'resolvent argument, any CaseII local criterion, or integral Hasse existence/nonexistence.'),
              'formal_verification': False, 'independent_whole_proof_replay': False, 'external_review': False}
    if args.write_receipt:
        (here/'fixed_map_two_adic_refinement_checks.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
