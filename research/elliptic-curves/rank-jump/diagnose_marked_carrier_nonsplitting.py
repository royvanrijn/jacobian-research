#!/usr/bin/env python3
"""Read the completed carrier bank; certify extension counts and local misses.

No enumeration, Riemann--Roch, point search or interrupted-candidate replay.
Only squarefree polynomials of degree one or two are accepted. For these,
the existing extension-key protocol reduces to a monic polynomial together
with its leading coefficient modulo rational squares. Factoring the large
rational constants is unnecessary: compare their ratios by exact isqrt.
"""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BANK = ROOT / 'artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1'
DEFAULT = ROOT / 'artifacts/generated-results/elliptic-curves/marked_carrier_closure_v1'


def read(p):
    return json.loads(p.read_text())


def rational_sqrt(q):
    q = Q(q)
    if q < 0:
        return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return Q(a, b) if a*a == q.numerator and b*b == q.denominator else None


def local_class(q, p):
    """Exact Q_p square test, including negative valuations and p=2 units."""
    assert q and p >= 2
    a, b, v = q.numerator, q.denominator, 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    modulus = 8 if p == 2 else p
    unit = a * pow(b, -1, modulus) % modulus
    unit_square = unit == 1 if p == 2 else pow(unit, (p-1)//2, p) == 1
    return {'prime': p, 'valuation': v, 'unit_modulus': modulus,
            'unit_residue': unit, 'unit_is_square': unit_square,
            'is_square': v % 2 == 0 and unit_square}


def compute():
    records = read(BANK/'carriers.json')['records']
    replay = {r['index']: r for r in read(BANK/'independent-replay.json')['nonsplit_at_control']}
    assert len(records) == len(replay) == 182
    t0 = Q(3, 17)
    primes = [p for p in range(2, 1010)
              if all(p % d for d in range(2, isqrt(p)+1))]
    groups, values, output = {}, [], []
    for r in sorted(records, key=lambda r: r['index']):
        c = list(map(Q, r['branch']))
        assert len(c) in (2, 3) and c[-1]
        assert len(c) == 2 or c[1]**2 - 4*c[0]*c[2] != 0
        monic = tuple(x/c[-1] for x in c)
        q = sum((a*t0**i for i, a in enumerate(c)), Q(0))
        assert q and rational_sqrt(q) is None
        # D(T,U)=sum c_i*T^i*U^(2-i) includes infinity for a linear d.
        homogeneous = sum((a*3**i*17**(2-i) for i, a in enumerate(c)), Q(0))
        reciprocal = homogeneous/9  # T=1 chart, U=17/3; valid also at p=17.
        assert homogeneous == 17**2*q
        local = [local_class(q, p) for p in primes]
        for p, datum in zip(primes, local):
            assert local_class(homogeneous, p)['is_square'] == datum['is_square']
            assert local_class(reciprocal, p)['is_square'] == datum['is_square']
        prior = replay[r['index']]['local_nonsplitting']
        assert not local_class(q, prior['prime'])['is_square']
        assert any(not v['is_square'] for v in local)
        bucket = groups.setdefault(monic, [])
        for g in bucket:
            root = rational_sqrt(c[-1]/g['constant'])
            if root is not None:
                g['indices'].append(r['index'])
                extension, multiplier = g['representative'], root
                break
        else:
            extension, multiplier = r['index'], Q(1)
            bucket.append({'constant': c[-1], 'indices': [r['index']],
                           'representative': r['index']})
        for g in values:
            root = rational_sqrt(q/g['value'])
            if root is not None:
                g['indices'].append(r['index'])
                value_class = g['indices'][0]
                break
        else:
            value_class = r['index']
            values.append({'value': q, 'indices': [r['index']]})
        witness = {'place': 'real', 'sign': -1} if q < 0 else next(
            {'place': str(v['prime']), **v} for v in local if not v['is_square'])
        output.append({'index': r['index'], 'shell_norm': r['divisor']['norm'],
                       'branch_ascending': list(map(str, c)),
                       'monic_branch_ascending': list(map(str, monic)),
                       'constant_representative': str(c[-1]),
                       'extension_representative_index': extension,
                       'square_multiplier_to_extension_representative': str(multiplier),
                       'value_at_control': str(q),
                       'control_squareclass_representative_index': value_class,
                       'homogeneous_value_at_3_17': str(homogeneous),
                       'reciprocal_chart_value_at_17_over_3': str(reciprocal),
                       'negative_real_sign': q < 0, 'small_witness': witness,
                       'local_classes': local})
    counts = {str(p): sum(not row['local_classes'][i]['is_square'] for row in output)
              for i, p in enumerate(primes)}
    ext = [g for bucket in groups.values() for g in bucket]
    result = {
        'schema': 'marked-carrier.nonsplitting-closure.v1',
        'inputs': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (BANK/'carriers.json', BANK/'independent-replay.json')},
        'control_parameter': str(t0), 'carrier_count': len(output),
        'degree_counts': dict(Counter(len(r['branch_ascending'])-1 for r in output)),
        'distinct_monic_branch_polynomials': len(groups),
        'distinct_quadratic_extensions_over_fixed_t': len(ext),
        'extension_multiplicities': dict(Counter(len(g['indices']) for g in ext)),
        'distinct_control_rational_squareclasses': len(values),
        'negative_real_count': sum(r['negative_real_sign'] for r in output),
        'tested_primes': primes, 'nonsplitting_counts_by_prime': counts,
        'common_obstruction_primes_in_tested_list': [p for p in primes if counts[str(p)] == len(output)],
        'boundary': 'Only the 182 completed maps. Both shells remain incomplete; the 73 split residuals and one interrupted candidate are unchanged. No class labels were tested. Absence of a common obstruction is asserted only for the real place and listed primes, not all places.',
        'records': output,
    }
    stream = io.StringIO(newline='')
    fields = ['index', 'shell_norm', 'branch_ascending', 'constant_representative',
              'monic_branch_ascending', 'extension_representative_index',
              'value_at_control', 'homogeneous_value_at_3_17', 'small_witness']
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for r in output:
        writer.writerow({k: json.dumps(r[k], separators=(',', ':')) if isinstance(r[k], (list, dict))
                         else r[k] for k in fields})
    return result, stream.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT)
    parser.add_argument('--check', action='store_true', help='Recompute only this diagnostic and compare its saved result.')
    args = parser.parse_args()
    result, table = compute()
    if args.check:
        assert read(args.output/'result.json') == json.loads(json.dumps(result))
        assert (args.output/'nonsplitting.csv').read_bytes() == table.encode()
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        (args.output/'result.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        (args.output/'nonsplitting.csv').write_bytes(table.encode())
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('records', 'inputs', 'tested_primes', 'nonsplitting_counts_by_prime')}, indent=2))


if __name__ == '__main__':
    main()
