#!/usr/bin/env sage-python
"""Exact conductor upper bound for the new27 curve, without residual factoring."""
import argparse, json, sys
from pathlib import Path
from math import gcd
from sage.all import QQ, EllipticCurve
from sage.version import version
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as minimal
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'full11952_high_rank_models_v1.json'
CAT=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json'
OUT=ART/'full11952_conductor_bound_v1.json'

def expected():
    data=cert.read(INPUT)
    row=next(r for r in data['curves'] if r['id']=='11952-0962587')
    model=tuple(map(int,row['integral_curve']))
    proof=minimal.minimality(model)
    if row['minimality']!={'status':'PROVED_GLOBAL_MINIMAL','certificate':proof}:
        raise ArithmeticError('global minimality replay differs')
    inv=cert.weierstrass_invariants(model)
    c4,c6,delta=[int(inv[k]) for k in ('c4','c6','discriminant')]
    E=EllipticCurve(QQ,model)
    if tuple(map(int,E.c_invariants()))!=(c4,c6) or int(E.discriminant())!=delta:
        raise ArithmeticError('independent invariant arithmetic differs')
    if gcd(abs(c4),abs(c6))!=75 or delta%3==0:
        raise ArithmeticError('only5 may have additive bad reduction')
    local=E.local_data(5,algorithm='generic')
    local5={'prime':5,'kodaira':str(local.kodaira_symbol()),
            'conductor_exponent':int(local.conductor_valuation()),
            'discriminant_exponent':int(local.discriminant_valuation())}
    if local5!={'prime':5,'kodaira':'IV','conductor_exponent':2,'discriminant_exponent':4}:
        raise ArithmeticError('fixed tame local5 calculation differs')
    remainder=abs(delta);factors=[]
    for p in _primes_up_to(10000):
        e=0
        while remainder%p==0:remainder//=p;e+=1
        if e:factors.append([p,e])
    reconstructed=remainder;bound=25*remainder
    for p,e in factors:
        reconstructed*=p**e
        if p!=5:
            if c4%p==0:raise ArithmeticError('multiplicative invariant gate failed')
            bound*=p
    if reconstructed!=abs(delta) or gcd(remainder,75)!=1:
        raise ArithmeticError('exact discriminant decomposition differs')
    # For every p>3, c4=0 mod p and Delta=0 mod p force c6=0 mod p.
    # The gcd therefore makes all remaining bad primes multiplicative.
    catalogue=cert.read(CAT);cohort=[r for r in catalogue['curves'] if int(r['rank_lower_bound'])>=27]
    known=[r for r in cohort if r.get('conductor') is not None]
    comparison={'rank_lower_bound_threshold':27,'catalogue_curves':len(cohort),
        'with_recorded_conductor':len(known),
        'missing_conductor_ids':sorted(r['id'] for r in cohort if r.get('conductor') is None),
        'recorded_conductor_strictly_above_bound_ids':sorted(r['id'] for r in known if int(r['conductor'])>bound),
        'recorded_conductor_at_most_bound_ids':sorted(r['id'] for r in known if int(r['conductor'])<=bound)}
    paths=[Path(__file__).resolve(),INPUT,CAT,Path(cert.__file__),Path(minimal.__file__),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']
    return {'schema':'elliptic-curves.full11952-conductor-bound.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'sage_version':version,
        'inventory_id':'new-20260906-186','source_id':row['id'],'minimal_model':list(map(str,model)),
        'minimality_certificate':proof,'c4':str(c4),'c6':str(c6),'discriminant':str(delta),
        'local5':local5,'trial_prime_bound':10000,'trial_prime_factors':factors,
        'remaining_cofactor':str(remainder),'conductor_upper_bound':str(bound),
        'exact_conductor':'UNKNOWN','catalogue_comparison':comparison,
        'argument':'The integral equation is globally minimal. The invariant gcd is75 and3 is good. At5 the exact generic Tate calculation gives typeIV and conductor exponent2. Every other bad prime is multiplicative and has conductor exponent1: at2 c4 is odd; at primes above3 the invariant identity forces common c4,c6 divisibility if c4 and Delta vanish. Thus N=25 times the product of the displayed bad primes other than5 times rad(R). Replacing rad(R) by the positive remaining cofactor R gives the stated upper bound. No primality or squarefreeness claim on R is required.',
        'claim_boundary':'This is an unconditional upper bound on the conductor of the certified rank-at-least27 curve, not an exact conductor or rank. Comparisons use the recorded conductors and rank lower bounds in the pinned593-entry catalogue without independently certifying those catalogue conductors. Missing catalogue values, other inventory curves, universal novelty and record status are not resolved.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('conductor bound replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve conductor bound certificate')
        checkpoint(OUT,d)
    print('CONDUCTOR UPPER BOUND',d['conductor_upper_bound'],d['catalogue_comparison'],flush=True)
