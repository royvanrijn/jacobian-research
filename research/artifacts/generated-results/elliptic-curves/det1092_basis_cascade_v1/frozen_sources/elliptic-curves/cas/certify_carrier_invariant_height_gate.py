#!/usr/bin/env python3
"""Factor-free height obstruction on the fixed twelve-word carrier sample."""
import argparse,json
from math import gcd,isqrt
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'soluble_pair_carrier_height_v1.json';OUT=ART/'soluble_pair_carrier_invariant_height_v1.json'

def expected():
    source=cert.read(INPUT);rows=[];limit=2**360-1
    for r in source['rows']:
        if 'curve' not in r:
            rows.append({'word':r['word'],'status':'NO_SPECIALIZED_MODEL_IN_INPUT'});continue
        model=tuple(map(cert.F,r['curve']));inv=cert.weierstrass_invariants(model)
        if any(a.denominator!=1 for a in model) or not inv['discriminant']:raise ArithmeticError('integral nonsingular source required')
        c4,c6=int(inv['c4']),int(inv['c6']);G=gcd(abs(c4),abs(c6))
        if not G:raise ArithmeticError('nonzero invariant gcd required')
        # Any Q-isomorphism with scale u has integral target c4/u^4,c6/u^6.
        # In lowest terms u=a/b, a^4 divides G; hence |u|^12<=G^3.
        square_numerator=c6*c6;square_denominator=G**3
        threshold=square_denominator*(1224*limit+521)**2
        square_ceiling=(square_numerator+square_denominator-1)//square_denominator
        lower_c6=isqrt(square_ceiling)
        if lower_c6**2<square_ceiling:lower_c6+=1
        lower_coefficient=max(0,(lower_c6-521+1223)//1224)
        rows.append({'word':r['word'],'compact_parameter':r['compact_parameter'],
            'recorded_model_bits':r['model_coefficient_bits'],'invariant_gcd':str(G),
            'c4':str(c4),'c6':str(c6),
            'normalized_coefficient_size_lower_bound':str(lower_coefficient),
            'normalized_coefficient_bits_lower_bound':lower_coefficient.bit_length(),
            'all_normalized_360_bit_integral_models_excluded':square_numerator>threshold,
            'status':'PROVED_NORMALIZED_HEIGHT_OBSTRUCTION' if square_numerator>threshold else 'NOT_EXCLUDED'})
    checked=[r for r in rows if 'compact_parameter' in r];excluded=[r for r in checked if r['all_normalized_360_bit_integral_models_excluded']]
    if len(checked)!=11 or len(excluded)!=10 or [r['word'] for r in checked if r not in excluded]!=[[1,1]]:
        raise ArithmeticError('fixed sample height gate differs')
    return {'schema':'elliptic-curves.soluble-pair-carrier-invariant-height.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,Path(cert.__file__))},
        'coefficient_bit_limit':360,'rows':rows,'distinct_parameters_excluded':len({r['compact_parameter'] for r in excluded}),
        'lemma':'For integral source invariants c4,c6 let G=gcd(|c4|,|c6|). An isomorphic integral equation has invariants c4/u^4,c6/u^6. Writing the rational scale u=a/b in lowest terms, a^4 divides G, hence |u|^12<=G^3 and |c6_target|^2>=c6^2/G^3. A normalized integral equation with a1,a3 in{0,1}, a2 in{-1,0,1}, and |a4|,|a6|<=M satisfies |b2|<=5, |b4|<=2M+1, |b6|<=4M+1, whence |c6_target|<=125+180(2M+1)+216(4M+1)=1224M+521. The recorded integer strict inequality therefore excludes every normalized integral equation within the360-bit gate, including a normalized global minimal model. It also excludes an integral short model within that gate, since then |c6_target|<=864M.',
        'claim_boundary':'All ten nonanchor specialized rows, representing seven distinct parameters, fail the360-bit gate even after arbitrary rational scaling and integral normalization. The anchor is not excluded; one larger parameter was censored before specialization. No factorization, explicit global minimal model, new parameter search or point search is performed. This does not exclude small models on other carrier points or arbitrary unnormalized equations with all five coefficients bounded.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('invariant height replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve invariant height proof')
        checkpoint(OUT,d)
    print('CARRIER NORMALIZED HEIGHT EXCLUSIONS',d['distinct_parameters_excluded'],[(r['word'],r.get('normalized_coefficient_bits_lower_bound')) for r in d['rows']],flush=True)
