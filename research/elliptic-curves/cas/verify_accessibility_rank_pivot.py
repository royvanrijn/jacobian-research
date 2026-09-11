#!/usr/bin/env python3
"""Check the finite singleton-rule theorem and equation-only rank criterion.

Reuses the sealed independent slope replay; does not repeat point arithmetic
or claim that numerical height matrices are certified enclosures.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import gzip
import json
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
SUB = ART / 'rank_accessibility_subsets_v1'
RANK = ART / '11952_rank_bound_audit_v1'
OUT = ART / 'accessibility_rank_pivot_v1'


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    raw = p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix == '.gz' else raw)


def closure():
    verified = read(SUB / 'verified.json')
    assert verified['status'] == 'PASS'
    assert verified['verifier_sha256'] == digest(ROOT / 'elliptic-curves/cas/verify_rank_accessibility_subsets.py')
    answer = []
    for cid, q, expected in [('302', 14, {(1, 10), (1, 14), (2, 7), (2, 10), (6, 12)}),
                             ('11952', 8, {(4, 2)})]:
        path = SUB / f'{cid}-original.json.gz'
        replay = next(r for r in verified['results'] if r['file'] == path.name)
        assert replay['sha256'] == digest(path)
        data = read(path)
        rules = []
        cells = 0
        for row, envelope in zip(data['target_landscapes'], replay['minimal_support_envelopes']):
            j = row['target']
            assert j == envelope['target']
            hs = [max(abs(int(r['n'])), int(r['d'])) for r in row['measurements']]
            support_min = {}
            for idx, mask in row['native_anchor_supports']:
                assert mask & (1 << j) == 0
                support_min[mask] = min(support_min.get(mask, hs[idx]), hs[idx])
            h0 = support_min[0]
            assert h0 == int(envelope['initial_H'])
            singletons = {s: h for s, h in support_min.items() if s.bit_count() == 1}
            # Every two-exceptional chart is dominated by an empty/singleton
            # support: this proves the formula for ALL subsets symbolically.
            for s, h in support_min.items():
                assert h >= min([h0] + [v for t, v in singletons.items() if t & s == t])
            for s, h in singletons.items():
                if h < h0:
                    matching = [o for o in envelope['offers'] if o['exceptional_support_mask'] == s]
                    assert len(matching) == 1 and int(matching[0]['H']) == h
                    rules.append({'source': s.bit_length(), 'target': j + 1,
                                  'initial_H': str(h0), 'H': str(h),
                                  'height_ratio': str(Fraction(h0, h)),
                                  'anchor_word_in_M17_E': matching[0]['anchor_word_in_original_M17_E']})
            assert row['subset_masks'] == [s for s in range(1 << q) if not s & (1 << j)]
            for s, idx in zip(row['subset_masks'], row['winning_anchor_ids']):
                assert hs[idx] == min([h0] + [v for t, v in singletons.items() if t & s == t])
                cells += 1
        edges = {(r['source'], r['target']) for r in rules}
        assert edges == expected
        sources = {a for a, b in edges}
        targets = {b for a, b in edges}
        assert sources.isdisjoint(targets)
        assert cells == q * 2 ** (q - 1)
        answer.append({'fibre': cid, 'rules': rules, 'cells_rechecked': cells,
                       'no_directed_path_of_length_two': True,
                       'all_two_exceptional_supports_dominated': True,
                       'improvement_function': 'max({0} union {log2(initial_H/H_ij): i in T})',
                       'monotone_submodular': True})
    return answer


def arithmetic():
    protocol = read(RANK / 'protocol.json')
    for name, value in protocol['inputs'].items():
        assert digest(ROOT / name) == value
    field = read(RANK / 'field.json')
    bnf = read(RANK / 'class_group.json')
    for name in ['field', 'class_group']:
        r = read(RANK / (name + '.json'))
        assert r['program_sha256'] == protocol['programs'][name] == digest(RANK / (name + '.gp'))
        assert r['log_sha256'] == digest(RANK / (name + '.log'))
    assert field['valid_completion'] and field['maximal_order_certified']
    log = (RANK / 'field.log').read_text()
    fields = dict(line.split('|', 1) for line in log.splitlines() if '|' in line)
    values = lambda k: json.loads(fields[k])
    short = list(map(int, read(RANK / 'equation.json')['a_invariants']))
    model = values('MINIMAL_MODEL')
    u, r, s, t = values('CHANGE_TO_MINIMAL')
    assert (u, r, s, t) == (12, 48, 0, 0)
    assert model[:3] == [0, 1, 0]
    assert 3*r == model[1]*u**2
    assert 3*r*r+short[3] == model[3]*u**4
    assert r**3+short[3]*r+short[4] == model[4]*u**6
    c, b, a, leading = values('CUBIC_COEFFICIENTS_ASCENDING')
    assert (a, b, c, leading) == (model[1], model[3], model[4], 1)
    disc_f = a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    assert 16*disc_f == int(fields['MINIMAL_DISCRIMINANT'])
    assert disc_f == int(fields['FIELD_INDEX'])**2*int(fields['FIELD_DISCRIMINANT'])
    assert int(fields['MINIMAL_DISCRIMINANT'])*u**12 == -16*(4*short[3]**3+27*short[4]**2)
    assert all((x**3+a*x*x+b*x+c) % 23 for x in range(23))
    factors = [json.loads(line.split('|', 1)[1]) for line in log.splitlines() if line.startswith('DISC_FACTOR|')]
    assert prod(p**e for p, e in factors) == 16*disc_f
    local = [json.loads(line.split('|', 1)[1]) for line in log.splitlines() if line.startswith('LOCAL|')]
    assert [(r[0], r[1]) for r in local] == [tuple(r) for r in factors]
    multiplicative_even = []
    additive = []
    n = 0
    for p, vdisc, conductor_exp, kodaira, number_above, contribution in local:
        if conductor_exp == 1:
            assert int(fields['C4']) % p != 0
            assert contribution == (vdisc % 2 == 0)
            if contribution:
                multiplicative_even.append(p)
        else:
            assert int(fields['C4']) % p == 0 and conductor_exp > 1
            assert contribution == number_above-1
            additive.append([p, number_above])
        n += contribution
    assert multiplicative_even == [3, 7, 13, 19, 83]
    assert additive == [[2, 3]]
    assert disc_f > 0 and n == 7 and field['bound_offset'] == 9
    assert not bnf['valid_completion'] and bnf['unconditional_rank_upper_bound'] is None
    assert 'bnfinit: the PARI stack overflows' in (RANK / 'class_group.log').read_text()
    return {'status': 'EXACT_LOCAL_CRITERION_GLOBAL_CLASS_BOUND_UNKNOWN',
            'minimal_model': list(map(str, model)), 'cubic_discriminant': str(disc_f),
            'field_discriminant': fields['FIELD_DISCRIMINANT'],
            'u': 2, 'multiplicative_even_primes': multiplicative_even,
            'additive_prime_counts': additive, 'n': n, 'bound': 'rank <= dim Sel2 <= g+9',
            'sufficient_independent_class_2rank_upper': 16,
            'class_2rank_upper': None, 'rank_upper': None,
            'class_probe_outcome': 'PARI_STACK_OVERFLOW_BEFORE_INVARIANTS',
            'independence': 'Python integer model/discriminant/offset and mod23 checks. Maximal order, primality and local reduction rely on the pinned PARI computations, not an independent implementation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    paths = [SUB / 'verified.json', SUB / '302-original.json.gz', SUB / '11952-original.json.gz',
             SUB / 'events.json', RANK / 'protocol.json', RANK / 'field.json', RANK / 'class_group.json',
             Path(__file__).resolve()]
    result = {'schema': 'accessibility-rank-pivot.v1', 'status': 'PASS',
              'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in paths},
              'finite_accessibility': closure(), 'equation_only_rank_audit': arithmetic(),
              'scope': 'Finite-atlas diminishing returns, not a point-recovery or rank-amplification theorem. No new upper bound.'}
    raw = (json.dumps(result, sort_keys=True, indent=2) + '\n').encode()
    path = OUT / 'verified.json'
    if args.check:
        assert path.read_bytes() == raw
    else:
        OUT.mkdir(parents=True, exist_ok=True)
        if path.exists():
            assert path.read_bytes() == raw
        else:
            path.write_bytes(raw)
    print('PASS: 115712 finite cells; singleton domination; no two-step paths; equation-only g+9 criterion. Rank upper UNKNOWN.')


if __name__ == '__main__':
    main()
