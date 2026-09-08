#!/usr/bin/env sage-python
"""Prove partial conductor divisors and exact completions; no factor search."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from sage.all import QQ, ZZ, EllipticCurve, pari
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
WORK = ROOT/'artifacts/local/elliptic-curves/conductor-inventory-continuation-v2'
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/inventory201_conductor_bounds_v1.json'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/inventory201_conductors_v2'


def build(check=False):
    audit = cert.read(INPUT)
    rows = audit['rows']
    results = []
    for r in rows:
        directory = WORK/'factors'/r['id']
        output = OUT/(r['id']+'.json')
        if check:
            saved = cert.read(output)
            discovered = saved['discovery_factors']
        else:
            state = directory/('state.json' if (directory/'state.json').exists() else 'initial.json')
            discovered = cert.read(state)['factors']
        if ZZ.prod(ZZ(n) for n in discovered) != ZZ(r['remaining_cofactor']):
            raise ArithmeticError('discovery factor product differs')
        candidates = {str(q['prime']) for q in r['local_data']
                      if q['displayed_discriminant_valuation']}
        candidates.update(n for n in discovered if pari(ZZ(n)).ispseudoprime())
        certificates = {}
        known_conductor, remaining = ZZ(1), abs(ZZ(r['discriminant']))
        model = [QQ(a) for a in r['integral_curve']]
        E, ep = EllipticCurve(QQ,model), pari.ellinit(model)
        if E.discriminant() != ZZ(r['discriminant']) or ep.disc() != E.discriminant():
            raise ArithmeticError('model discriminant differs')
        local = []
        for text in sorted(candidates,key=int):
            p = ZZ(text)
            proof = pari(saved['prime_certificates'][text]) if check else pari(p).primecert()
            root = ZZ(proof) if proof.type() == 't_INT' else ZZ(proof[0][0])
            if root != p or not pari.primecertisvalid(proof):
                raise ArithmeticError('exact primality certificate failed')
            certificates[text] = str(proof)
            e = 0
            while remaining % p == 0:
                remaining //= p
                e += 1
            if not e:
                raise ArithmeticError('listed prime does not divide discriminant')
            sage_local = E.local_data(p,algorithm='generic',proof=True)
            pari_local = ep.elllocalred(p)
            f = int(sage_local.conductor_valuation())
            if f != int(pari_local[0]):
                raise ArithmeticError('two exact local conductor implementations disagree')
            known_conductor *= p**f
            local.append({'prime':text,'displayed_discriminant_valuation':e,
                'minimal_discriminant_valuation':int(sage_local.discriminant_valuation()),
                'conductor_exponent':f,'kodaira':str(sage_local.kodaira_symbol()),
                'pari_local_reduction':str(pari_local)})
        reconstructed = remaining * ZZ.prod(ZZ(q['prime'])**q['displayed_discriminant_valuation'] for q in local)
        if reconstructed != abs(E.discriminant()):
            raise ArithmeticError('partial discriminant reconstruction failed')
        if remaining % 2 == 0 or remaining % 3 == 0:
            raise ArithmeticError('unresolved wild prime')
        record = audit['benchmarks'][str(r['rank_lower_bound'])]
        benchmark = ZZ(record['recorded_conductor'])
        result = {'id':r['id'],'rank_lower_bound':r['rank_lower_bound'],
            'family':r['family'],'parameter':r['parameter'],'integral_curve':r['integral_curve'],
            'discriminant':r['discriminant'],'discovery_factors':discovered,
            'prime_certificates':certificates,'local_data':local,
            'unresolved_discriminant_cofactor':str(remaining),
            'conductor_divisor':str(known_conductor),
            'conductor_upper_bound':str(known_conductor*remaining),
            'exact_conductor':str(known_conductor) if remaining == 1 else None,
            'bad_primes':[q['prime'] for q in local if q['conductor_exponent']] if remaining == 1 else None,
            'status':'EXACT' if remaining == 1 else 'UNKNOWN',
            'recorded_benchmark':record,
            'proved_above_recorded_benchmark':known_conductor > benchmark,
            'proved_below_recorded_benchmark':known_conductor*remaining < benchmark,
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__),INPUT]},
            'sage_version':version,'pari_version':str(pari.version()),
            'argument':'Every listed prime has an exact primality certificate. Sage generic Tate '
                'and PARI elllocalred independently agree on its conductor exponent. Their '
                'product divides the true conductor. The unresolved cofactor is prime to2,3 '
                'and all processed primes, and bounds the remaining conductor by f_p<=v_p(Delta). '
                'Thus divisor <= N <= upper_bound. If the cofactor is1, N and the bad-prime list '
                'are exact, allowing nonminimal input models through local minimization. '
                'Public benchmarks are reported metadata, not universal records. No new rank claim.'}
        if check:
            if saved != result:
                raise ArithmeticError('certificate replay differs: '+r['id'])
        elif output.exists():
            if cert.read(output) != result:
                raise FileExistsError('preserve prior certificate: '+r['id'])
        else:
            cert.write(output,result)
        results.append(result)
    summary = {'schema':'elliptic-curves.inventory201-conductors.v2','status':'PASS',
        'curves':len(results),'exact_count':sum(r['status']=='EXACT' for r in results),
        'proved_above_own_rank_benchmark':sum(r['proved_above_recorded_benchmark'] for r in results),
        'improvement_ids':[r['id'] for r in results if r['proved_below_recorded_benchmark']],
        'priority_ids':cert.read(WORK/'factor_protocol.json')['priority_ids'],
        'benchmarks':audit['benchmarks'],'best_by_rank':{},
        'certificates':{r['id']:{'path':str((OUT/(r['id']+'.json')).relative_to(ROOT)),
            'sha256':cert.hashed(OUT/(r['id']+'.json'))} for r in results},
        'claim_boundary':'Partial factorizations remain UNKNOWN. Only a certified conductor divisor '
            'above the recorded minimum excludes improvement at that rank. No full-inventory '
            'exclusion follows merely from large upper bounds; rank means certified lower bound.'}
    for rank in range(22,28):
        cohort = [r for r in results if r['rank_lower_bound'] >= rank]
        exact = [r for r in cohort if r['status']=='EXACT']
        best = min(cohort,key=lambda r:int(r['conductor_upper_bound']))
        exact_best = min(exact,key=lambda r:int(r['exact_conductor'])) if exact else None
        summary['best_by_rank'][str(rank)] = {
            'best_upper_bound_id':best['id'],'best_upper_bound':best['conductor_upper_bound'],
            'best_exact_id':exact_best['id'] if exact_best else None,
            'best_exact_conductor':exact_best['exact_conductor'] if exact_best else None,
            'proved_above_benchmark_count':sum(ZZ(r['conductor_divisor']) > ZZ(audit['benchmarks'][str(rank)]['recorded_conductor']) for r in cohort),
            'cohort_size':len(cohort)}
    if check:
        if cert.read(OUT/'summary.json') != summary:
            raise ArithmeticError('screen summary differs')
    elif (OUT/'summary.json').exists():
        if cert.read(OUT/'summary.json') != summary:
            raise FileExistsError('preserve screen summary')
    else:
        cert.write(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('certificates','benchmarks','best_by_rank')},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    build(args.check)
