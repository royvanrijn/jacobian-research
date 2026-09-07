#!/usr/bin/env python3
"""Exact replay of root intervals, norm enumeration and retained ideal circuits."""
import argparse
from fractions import Fraction as Q
from math import gcd, prod, floor
from pathlib import Path
import hashlib
import retrospective as r
import matched103b2_root_circuits as source
import matched103b2_class_boundary as prior
from verify_unpointed_governing_norm import Algebra

OUTPUT = r.OUT/'rank_jump_matched103b2_root_circuits_verification_v1.json'


def compute():
    from sage.all import QQ, ZZ, GF, matrix, pari, prime_range
    data = r.read(source.OUTPUT)
    for path, sha in data['bindings'].items():
        assert r.digest((r.ROOT/path).read_bytes()) == sha
    spec = r.read(source.PROTOCOL)['limits']
    pari.allocatemem(64000000, spec['pari_stack_bytes'], silent=True)
    primorial = prod(map(int, prime_range(spec['smooth_bound']+1)))
    rows = []
    for row in data['rows']:
        token = row['token']; norm = row['norm']; audit = row['audit']
        form = next(x for x in r.read(source.INPUT)['rows'] if x['token'] == token)['form']
        c = list(map(int, form['binary_coefficients_descending']))
        def value(m, n):
            return sum(ci*m**(3-i)*n**i for i, ci in enumerate(c))
        # Three disjoint sign-changing rational intervals account for all cubic roots.
        intervals = [(Q(x['lower']), Q(x['upper'])) for x in norm['roots']]
        assert len(intervals) == 3
        assert all(lo < hi and value(lo, 1)*value(hi, 1) < 0 for lo, hi in intervals)
        assert all(intervals[j][1] < intervals[j+1][0] for j in range(2))
        digest = hashlib.sha256(); seen = set(); accepted = []; best = []
        for j, ((lo, hi), root) in enumerate(zip(intervals, norm['roots'])):
            assert len(root['floors']) == spec['denominators']
            for n, m0 in enumerate(root['floors'], 1):
                assert floor(n*lo) == floor(n*hi) == m0
                for shift in (-1, 0, 1):
                    m = m0+shift
                    if gcd(m, n) != 1 or (m, n) in seen:
                        continue
                    seen.add((m, n)); v = value(m, n); rem = abs(v)
                    assert rem
                    while (g := gcd(rem, primorial)) > 1:
                        rem //= g
                    digest.update(f'{j},{m},{n},{v},{rem}\n'.encode())
                    best = sorted(best+[(rem.bit_length(), j, m, n, str(rem))])[:10]
                    if rem.bit_length() <= spec['residual_bits']:
                        accepted.append(dict(root_index=j, m=m, n=n, value=str(v), remainder=str(rem)))
        assert digest.hexdigest() == norm['digest'] and len(seen) == norm['count']
        assert len(accepted) == norm['accepted_count']
        assert accepted[:spec['accepted_relations_per_case']] == norm['accepted']
        assert [list(x) for x in best] == norm['smallest_remainders']
        records = list(norm['accepted'])
        seed = next(x for x in r.read(source.SEED)['rows'] if x['token'] == token)
        for x in seed['retained_principal_relations']:
            if not any((y['m'], y['n']) == (x['m'], x['n']) for y in records):
                records.append(dict(m=x['m'], n=x['n'], value=x['polynomial_value'],
                                    remainder=x['residual'], source='retained_box'))
        if not records:
            assert audit['status'] == 'NO_RELATIONS' and audit['relation_count'] == 0
            rows.append(dict(token=token, norm_values=len(seen), root_floors=1536,
                             principal_ideals=0, additional_strict_image_dimension=0))
            continue
        oldrow, f, pts, primes, nf = prior.setup(token)
        S = oldrow['local']['S_finite']; K = Algebra(form['cubic_ascending'])
        a = Q(form['fixed_a']); w = K.elt(form['w_ascending']); M = form['SL2_matrix']
        assert len(records) == audit['relation_count'] == len(audit['relations'])
        divisors = []
        for rec, saved in zip(records, audit['relations']):
            assert all(saved[k] == v for k, v in rec.items())
            m, n = rec['m'], rec['n']; u, v = M[0]*m+M[1]*n, M[2]*m+M[3]*n
            beta = K.add(K.elt([a*u]), tuple(v*x for x in w))
            assert list(map(str, beta)) == saved['beta_ascending']
            N = K.norm(beta)
            assert N == a*a*value(m, n) and N.denominator == 1
            factors = saved['norm_factors']
            assert prod(p**e for p, e in factors) == abs(N.numerator)
            b = pari(f.parent()(list(map(lambda x: QQ(str(x)), beta))))(pari.Mod('z', pari(f)))
            I = pari.idealhnf(nf, 1); vals = []
            for p, e in factors:
                assert ZZ(p).is_prime(proof=True) and e > 0
                for j, P in enumerate(pari.idealprimedec(nf, p)):
                    v = int(pari.idealval(nf, b, P)); assert v >= 0
                    if v:
                        vals.append([p, j, v]); I = pari.idealmul(nf, I, pari.idealpow(nf, P, v))
            assert vals == saved['ideal_valuations']
            assert pari.idealhnf(nf, b) == pari.idealhnf(nf, I)
            divisors.append({(p, j): e for p, j, e in vals})
        places = sorted(set().union(*(set(v) for v in divisors)))
        A = [[v.get(p, 0) for v in divisors] for p in places]
        assert A == audit['integer_valuation_matrix']
        assert [list(p) for p in places] == audit['valuation_columns']
        rQ, r2 = int(matrix(QQ, A).rank()), int(matrix(GF(2), A).rank())
        outside = [row for p, row in zip(places, A) if p[0] not in S]
        outside_rank = int(matrix(GF(2), outside).rank())
        assert rQ == audit['valuation_rank_Q'] and r2 == audit['valuation_rank_F2']
        assert rQ-r2 == audit['ideal_two_torsion_capacity']
        # Full outside-S rank alone excludes every nonzero retained product,
        # even after multiplying arbitrary generic point classes.
        assert outside_rank == len(records) and oldrow['local']['strict_generic_dimension'] == 0
        assert audit['strict_product_kernel'] == [] and audit['certified_independent_strict_classes'] == 0
        rows.append(dict(token=token, norm_values=len(seen), root_floors=1536,
                         principal_ideals=len(records), valuation_rank_Q=rQ,
                         valuation_rank_F2=r2, outside_S_parity_rank=outside_rank,
                         additional_strict_image_dimension=0))
    files = (Path(__file__), source.OUTPUT, Path(prior.__file__),
             Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return dict(schema='rank-jump.matched103b2-root-circuit-verification.v1', status='PASS',
                bindings={str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in files},
                rows=rows, scope='Independent rational interval and norm replay; full retained ideal identities and outside-S parity. Local signatures are unnecessary to the zero-image certificate and are not independently replayed here.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check'])
    args = parser.parse_args(); result = compute()
    if args.mode == 'build':
        r.write_new(OUTPUT, result)
    else:
        assert result == r.read(OUTPUT)
    print(result['status'], result['rows'], flush=True)
