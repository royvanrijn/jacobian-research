#!/usr/bin/env sage-python
"""Fixed finite conductor bounds for every certified V14 inventory curve."""
import argparse,json,sys
from pathlib import Path
from sage.all import QQ,ZZ,EllipticCurve
from sage.version import version
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'new_high_rank_curve_index_v14.json'
CAT=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json'
D=ROOT/'artifacts/local/elliptic-curves/inventory187-conductor-bounds-v1'
OUT=ART/'inventory187_conductor_bounds_v1.json'

def sources():
    paths=[Path(__file__).resolve(),INPUT,CAT,Path(cert.__file__),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve conductor audit protocol')
    rows=cert.read(INPUT)['curves']
    if len(rows)!=187 or len({r['id'] for r in rows})!=187:raise ArithmeticError('fixed187 required')
    checkpoint(D/'protocol.json',{'sources':sources(),'roster':[r['id'] for r in rows],
        'trial_prime_bound':10000,'workers':1,'wall_seconds':120,'rss_bytes':1073741824,
        'scope':'Audit all187 V14 conductor upper bounds to close the coverage gap beyond the early31 and next12 audits. Exact local Tate data at2,3 and every discriminant prime through10000; leave every larger cofactor unfactored. No parameter search, point search, new rank assertion, descent or automatic factorization escalation. A bound above a recorded benchmark does not exclude a smaller true conductor.'})

def expected(write=False):
    protocol=cert.read(D/'protocol.json')
    if protocol['sources']!=sources():raise ArithmeticError('frozen conductor inputs changed')
    inputs=cert.read(INPUT)['curves'];catalogue=cert.read(CAT)['curves']
    if [r['id'] for r in inputs]!=protocol['roster']:raise ArithmeticError('roster differs')
    benchmarks={}
    for rank in sorted({r['rank_lower_bound'] for r in inputs}):
        cohort=[r for r in catalogue if r['rank_lower_bound']>=rank]
        listed=[r for r in cohort if r.get('conductor') and str(r['conductor']).isdigit()]
        best=min(listed,key=lambda r:int(r['conductor']))
        benchmarks[str(rank)]={'id':best['id'],'recorded_conductor':best['conductor'],
            'catalogue_curves':len(cohort),'missing_conductor_ids':sorted(r['id'] for r in cohort if not r.get('conductor'))}
    result={'schema':'elliptic-curves.inventory187-conductor-bounds.v1','status':'RUNNING',
        'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'sage_version':version,
        'benchmarks':benchmarks,'rows':[],
        'argument':'Displayed short equations have coefficient denominators supported on2,3 and integral discriminant. At2,3 use exact local Tate exponents even if absent from displayed Delta. For every small bad prime use exact local Tate data. At all unprocessed primes p>=5 the equation is integral, and the conductor-discriminant inequality f_p<=v_p(Delta_min)<=v_p(Delta_displayed) bounds the remaining conductor by its full positive discriminant cofactor. Thus the actual conductor divides each reported upper bound. No squarefreeness or primality assumption on the cofactor is used.',
        'claim_boundary':protocol['scope']+' Catalogue conductors are recorded comparisons, not independently verified by this audit; missing values and universal records remain unresolved.'}
    for r in inputs:
        model=tuple(map(cert.F,r['curve']));E=EllipticCurve(QQ,[QQ(str(a)) for a in model]);inv=cert.weierstrass_invariants(model)
        delta=cert.F(str(E.discriminant()))
        if delta!=inv['discriminant'] or delta.denominator!=1 or not delta:raise ArithmeticError('integral independent discriminant gate')
        for a in model:
            den=a.denominator
            for p in (2,3):
                while den%p==0:den//=p
            if den!=1:raise ArithmeticError('nonintegral away from2,3')
        remaining=abs(delta.numerator);factor=1;local=[]
        for p in _primes_up_to(protocol['trial_prime_bound']):
            e=0
            while remaining%p==0:remaining//=p;e+=1
            if not e and p not in (2,3):continue
            data=E.local_data(p,algorithm='generic',proof=True);f=int(data.conductor_valuation());factor*=p**f
            local.append({'prime':p,'displayed_discriminant_valuation':e,
                'minimal_discriminant_valuation':int(data.discriminant_valuation()),
                'conductor_exponent':f,'kodaira':str(data.kodaira_symbol())})
        reconstructed=remaining
        for q in local:reconstructed*=q['prime']**q['displayed_discriminant_valuation']
        if reconstructed!=abs(delta.numerator):raise ArithmeticError('factor reconstruction differs')
        bound=factor*remaining;benchmark=int(benchmarks[str(r['rank_lower_bound'])]['recorded_conductor'])
        row={'id':r['id'],'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':r['rank_lower_bound'],
            'curve':r['curve'],'discriminant':str(delta),'local_data':local,'remaining_cofactor':str(remaining),
            'conductor_upper_bound':str(bound),'exact_conductor':str(bound) if remaining==1 else 'UNKNOWN',
            'upper_bound_below_recorded_minimum':bound<benchmark}
        result['rows'].append(row)
        if write:checkpoint(D/'checkpoint.json',result)
    result['status']='PASS_COMPLETE_FIXED_LOCAL_AUDIT'
    result['strict_improvement_ids']=[r['id'] for r in result['rows'] if r['upper_bound_below_recorded_minimum']]
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        if a.stage=='build' and OUT.exists():raise FileExistsError('preserve conductor result')
        d=expected(write=a.stage=='build')
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('conductor audit replay differs')
        else:checkpoint(OUT,d)
        print('INVENTORY187 CONDUCTOR BOUNDS',len(d['rows']),'STRICT IMPROVEMENTS',d['strict_improvement_ids'],flush=True)
