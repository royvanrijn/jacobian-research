#!/usr/bin/env sage-python
"""Exact principal-ideal relations from the fixed18 smooth norm values."""
import argparse
import json
from math import prod
from pathlib import Path
from sage.all import QQ, pari
import prepare_small_conductor_norm_form as forms
import pilot_small_conductor_norm_smoothness as pilot

ROOT, ART = forms.ROOT, forms.ART
OUTPUT = ART / 'small_conductor_norm_relations_v1.json'


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            j = row.bit_length()-1
            if j not in pivots:
                pivots[j] = row
                break
            row ^= pivots[j]
    return len(pivots)


def expected():
    data = forms.expected()
    if json.loads(json.dumps(data)) != forms.target.original.cert.read(forms.OUT):
        raise ArithmeticError('norm-form proof differs')
    p = pilot.protocol()
    summary = forms.target.original.cert.read(pilot.OUT)
    if forms.target.original.cert.read(pilot.D / 'ledger.json')['status'] != 'PASS':
        raise ArithmeticError('finite norm pilot replay required')
    arm = forms.target.original.cert.read(pilot.D / 'reduced_binary.json')
    if summary['arms'][1]['records_sha256'] != forms.target.original.cert.hashed(pilot.D / 'reduced_binary.json'):
        raise ArithmeticError('pilot records changed')
    inputs = [forms.OUT, pilot.OUT, pilot.D / 'protocol.json', pilot.D / 'reduced_binary.json', Path(__file__).resolve()]
    proof = forms.target.original.cert.read(forms.target.PROOF)
    polynomial = pari.Pol(list(map(int, data['original_monic_cubic_descending'])))
    primes = [int(q) for q, e in proof['discriminant_factorization']]
    pari.addprimes(primes)
    nf = pari.nfinit([polynomial, primes])
    if len(pari.nfcertify(nf)):
        raise ArithmeticError('maximal order certification failed')
    if str(nf.disc()) != data['field_discriminant']:
        raise ArithmeticError('field discriminant differs')
    x = pari.Mod('x', polynomial)
    a = int(data['integral_norm_generator']['fixed_a'])
    w = sum(pari(QQ(c))*x**i for i, c in enumerate(data['integral_norm_generator']['w_power_basis']))
    M = data['sl2_matrix']
    relations = []
    support = set()
    for record in arm['records']:
        if record['remainder'] != '1':
            continue
        value = int(record['value'])
        factors = record['factorization']
        if prod(q**e for q, e in factors) != abs(value) or any(q > p['smooth_bound'] or not pari.isprime(q) for q, e in factors):
            raise ArithmeticError('smooth factorization failed')
        m, n = record['m'], record['n']
        u, v = M[0][0]*m+M[0][1]*n, M[1][0]*m+M[1][1]*n
        beta = a*u+v*w
        if pari.nfeltnorm(nf, beta) != a*a*value:
            raise ArithmeticError('principal norm identity failed')
        factor = pari.idealfactor(nf, beta)
        if pari.idealhnf(nf, pari.idealfactorback(nf, factor)) != pari.idealhnf(nf, beta):
            raise ArithmeticError('principal ideal factorization product failed')
        row = []
        for i in range(factor.nrows()):
            ideal, exponent = factor[i, 0], int(factor[i, 1])
            q = int(ideal[0])
            if q > p['smooth_bound'] or exponent < 0:
                raise ArithmeticError('relation outside fixed factor support')
            support.add(q)
            row.append({'rational_prime': q, 'ideal_hnf': str(pari.idealhnf(nf, ideal)), 'exponent': exponent})
        relations.append({'m': m, 'n': n, 'value': str(value),
                          'beta_power_basis': [str(pari.polcoef(pari.lift(beta), i)) for i in range(3)],
                          'factorization': row})
    columns = []
    canonical = []
    for q in sorted(support):
        ideals = pari.idealprimedec(nf, q)
        product = pari.idealhnf(nf, 1)
        row = 0
        for ideal in ideals:
            e = int(ideal[2])
            hnf = str(pari.idealhnf(nf, ideal))
            row |= (e % 2) << len(columns)
            columns.append({'rational_prime': q, 'ideal_hnf': hnf})
            product = pari.idealmul(nf, product, pari.idealpow(nf, ideal, e))
        if pari.idealhnf(nf, product) != pari.idealhnf(nf, q):
            raise ArithmeticError('canonical rational principal relation failed')
        canonical.append(row)
    lookup = {(r['rational_prime'], r['ideal_hnf']): i for i, r in enumerate(columns)}
    actual = []
    for relation in relations:
        row = 0
        for r in relation['factorization']:
            row ^= (r['exponent'] % 2) << lookup[(r['rational_prime'], r['ideal_hnf'])]
        actual.append(row)
    baseline = rank(canonical)
    gain = rank(canonical+actual)-baseline
    return {'schema': 'elliptic-curves.small-conductor-norm-relations.v1', 'status': 'PASS',
            'sources': {str(q.relative_to(ROOT)): forms.target.original.cert.hashed(q) for q in inputs},
            'pari_version': str(pari.version()), 'field_discriminant': data['field_discriminant'],
            'relations': relations, 'columns': columns,
            'canonical_rows': list(map(str, canonical)), 'relation_rows': list(map(str, actual)),
            'canonical_rank': baseline, 'additional_relation_rank': gain,
            'factor_base_generation_certified': False, 'class_rank_upper_bound': None,
            'claim_boundary': 'Exact supported principal-ideal relations and their additional parity-row rank modulo canonical rational-prime relations, on the displayed finite columns only. The columns are not proved to generate Cl(K)/2; this is not a class-group, Selmer or curve-rank upper bound.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = expected()
    if args.check:
        if result != forms.target.original.cert.read(OUTPUT):
            raise ArithmeticError('principal-relation proof differs')
    else:
        if OUTPUT.exists():
            raise FileExistsError('preserve principal relations')
        forms.target.original.cert.write(OUTPUT, result)
    print('EXACT PRINCIPAL RELATIONS', len(result['relations']), 'ADDITIONAL MOD2 ROW RANK', result['additional_relation_rank'])
