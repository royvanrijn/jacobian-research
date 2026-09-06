#!/usr/bin/env python3
"""Fixed finite smoothness comparison for the exact norm-form preparation."""
import argparse
import json
from math import gcd, prod
from pathlib import Path
import sys
import time

import prepare_small_conductor_norm_form as forms
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT, ART = forms.ROOT, forms.ART
D = ROOT / 'artifacts/local/elliptic-curves/small-conductor-norm-pilot-v1'
OUT = ART / 'small_conductor_norm_smoothness_pilot_v1.json'


def primes_to(bound):
    sieve = bytearray(b'\x01')*(bound+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, int(bound**0.5)+1):
        if sieve[p]:
            sieve[p*p:bound+1:p] = b'\x00'*(((bound-p*p)//p)+1)
    return [p for p in range(2, bound+1) if sieve[p]]


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve norm pilot protocol')
    expected = forms.expected()
    if json.loads(json.dumps(expected)) != forms.target.original.cert.read(forms.OUT):
        raise ArithmeticError('exact norm-form gate failed')
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.small-conductor-norm-pilot.v1',
        'sources': {str(p.relative_to(ROOT)): forms.target.original.cert.hashed(p)
                    for p in [Path(__file__).resolve(), forms.OUT]},
        'box': 64, 'smooth_bound': 400000, 'seconds': 120,
        'rss_bytes': 1610612736, 'workers': 1,
        'selection': 'Every coprime(m,n) with -64<=m<=64 and1<=n<=64, in lexicographic order, on each of the original monic and new maximal-order reduced binary norm polynomials.',
        'gate': 'Only this fixed curve; exact maximal-order norm preparation reduces maximum coefficient size165to65bits. Test whether this change produces complete smooth values before any class-relation campaign.',
        'method': 'Repeated gcd with the exact primorial through400000 removes all and only supported prime factors. Unresolved remainders are retained without generic factorization.',
        'claim_boundary': 'Finite norm-smoothness yield on different explicitly mapped principal-element sets. No class relation independence, generating factor base, class-rank upper bound, calibrated runtime for full descent or improved rank follows.'})


def protocol():
    p = forms.target.original.cert.read(D / 'protocol.json')
    if any(forms.target.original.cert.hashed(ROOT / name) != h for name, h in p['sources'].items()):
        raise ArithmeticError('frozen norm pilot sources changed')
    return p


def calculate(check=False):
    p = protocol()
    prime_list = primes_to(p['smooth_bound'])
    primorial = prod(prime_list)
    data = forms.target.original.cert.read(forms.OUT)
    rows = []
    for name, key in [('original_monic', 'original_monic_cubic_descending'),
                      ('reduced_binary', 'reduced_binary_cubic_descending')]:
        c = list(map(int, data[key]))
        started = time.monotonic()
        records = []
        for m in range(-p['box'], p['box']+1):
            for n in range(1, p['box']+1):
                if gcd(m, n) != 1:
                    continue
                value = sum(v*m**(3-i)*n**i for i, v in enumerate(c))
                if not value:
                    raise ArithmeticError('irreducible cubic cannot vanish')
                remainder = abs(value)
                while True:
                    common = gcd(remainder, primorial)
                    if common == 1:
                        break
                    remainder //= common
                record = {'m': m, 'n': n, 'value': str(value), 'remainder': str(remainder)}
                if remainder == 1:
                    q = abs(value)
                    factors = []
                    for prime in prime_list:
                        e = 0
                        while q % prime == 0:
                            q //= prime
                            e += 1
                        if e:
                            factors.append([prime, e])
                        if q == 1:
                            break
                    if q != 1:
                        raise ArithmeticError('smooth value factorization incomplete')
                    record['factorization'] = factors
                records.append(record)
            if not check and (m+p['box']) % 8 == 0:
                checkpoint(D / (name + '.checkpoint.json'), {'last_m': m, 'records': records})
        elapsed = time.monotonic()-started
        path = D / (name + '.json')
        actual = {'coefficients_descending': list(map(str, c)), 'records': records}
        if check:
            old = forms.target.original.cert.read(path)
            if any(old[k] != v for k, v in actual.items()):
                raise ArithmeticError('finite norm or smoothness replay differs')
            elapsed = old['wall_seconds']
        else:
            if path.exists():
                raise FileExistsError('preserve norm arm')
            checkpoint(path, {**actual, 'wall_seconds': elapsed})
        rows.append({'arm': name, 'primitive_pairs': len(records),
                     'smooth_values': sum(r['remainder'] == '1' for r in records),
                     'wall_seconds': elapsed, 'records_sha256': forms.target.original.cert.hashed(path)})
        print(name, rows[-1], flush=True)
    result = {'schema': 'elliptic-curves.small-conductor-norm-smoothness.v1', 'status': 'PASS',
              'protocol': p, 'arms': rows, 'claim_boundary': p['claim_boundary']}
    if check:
        if forms.target.original.cert.read(OUT) != result:
            raise ArithmeticError('smoothness summary differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve norm pilot summary')
        checkpoint(OUT, result)


def launch():
    p = protocol()
    if (D / 'ledger.json').exists():
        raise FileExistsError('preserve norm pilot ledger')
    stages = []
    for stage in ['worker', 'check']:
        result = run([sys.executable, str(Path(__file__).resolve()), stage],
                     limits=Limits(p['seconds'], p['rss_bytes']), cwd=ROOT,
                     log_path=D / (stage + '.log'), checkpoint_path=D / (stage + '.supervisor.json'))
        stages.append({'stage': stage, 'supervision': result})
        passed = result['outcome'] == 'completed' and result['returncode'] == 0
        checkpoint(D / 'ledger.json', {'status': 'PASS' if passed and len(stages) == 2 else 'RUNNING' if passed else 'FAILED_OR_CENSORED', 'stages': stages})
        if not passed:
            return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'launch', 'worker', 'check'])
    args = parser.parse_args()
    if args.stage in ('worker', 'check'):
        calculate(args.stage == 'check')
    else:
        globals()[args.stage]()
