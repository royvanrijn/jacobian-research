#!/usr/bin/env sage-python
"""Finite conductor-upper-bound screen for the31 certified new rank>=22 curves."""
from pathlib import Path
import sys
from sage.all import QQ,ZZ,EllipticCurve,prime_range,prod
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
NAMES=('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json','compact_r17_top64_interim_curves_v1.json',
       'compact_r17_largest_gain_curve_v1.json','compact_atlas_new_curves_v1.json','prospective_mw16_results_v1.json')
OUT=ROOT/'artifacts/local/elliptic-curves/new-curve-conductor-bounds-v1/result.json'

def run():
    if OUT.exists():raise FileExistsError('preserve bounded conductor audit')
    inputs=[]
    for name in NAMES:
        for r in cert.read(ART/name)['curves']:
            rank=r['rank_certificate']['rank_lower_bound']
            if rank>=22 and not r.get('icarm_matches') and not r.get('previous_matches'):
                inputs.append((name,r,rank))
    if len(inputs)!=31:raise ArithmeticError('frozen31-input roster changed')
    catalogue=cert.read(cert.DATABASE)['curves'];benchmarks={}
    for rank in sorted({rank for _,_,rank in inputs}):
        candidates=[r for r in catalogue if r['rank_lower_bound']>=rank and str(r['conductor']).isdigit()]
        best=min(candidates,key=lambda r:ZZ(r['conductor']));benchmarks[str(rank)]={'id':best['id'],'conductor':best['conductor']}
    result={'schema':'elliptic-curves.new-curve-conductor-bounds.v1','status':'RUNNING','rows':[],'benchmarks':benchmarks,
        'source_sha256':cert.hashed(Path(__file__).resolve()),'certificate_hashes':{name:cert.hashed(ART/name) for name in NAMES},
        'catalogue_sha256':cert.hashed(cert.DATABASE),
        'argument':'At primes>=5 the displayed short models are integral, so Ogg conductor-discriminant formula gives f_p<=v_p(Delta). Exact local computations replace every small-prime exponent, including2 and3 even when absent from displayed Delta. The unprocessed discriminant cofactor bounds the remaining conductor factor. This proves N divides the stated upper bound; it is not an exact conductor unless the cofactor is1.',
        'claim_boundary':'Finite comparison to listed catalogue conductors only. Failure of an upper bound to beat the benchmark does not exclude a smaller actual conductor.'}
    checkpoint(OUT,result)
    for name,r,rank in inputs:
        model=[QQ(q) for q in r['curve']]
        for a in model:
            d=ZZ(a.denominator())
            for p in (2,3):d//=p**d.valuation(p)
            if d!=1:raise ArithmeticError('coefficient denominator outside2,3')
        E=EllipticCurve(model);D=E.discriminant()
        if D.denominator()!=1:raise ArithmeticError('displayed discriminant not integral')
        remaining=abs(ZZ(D));local=[];known=ZZ(1)
        for p in prime_range(2,10001):
            e=int(remaining.valuation(p))
            if not e and p not in (2,3):continue
            remaining//=p**e;data=E.local_data(p,proof=True);f=int(data.conductor_valuation());known*=p**f
            local.append({'prime':str(p),'displayed_discriminant_valuation':e,'conductor_valuation':f,'kodaira_symbol':str(data.kodaira_symbol())})
        upper=known*remaining;benchmark=ZZ(benchmarks[str(rank)]['conductor'])
        row={'certificate':name,'family':r.get('family','published-R17'),'parameter':r['parameter'],'curve':r['curve'],'rank_lower_bound':rank,
             'displayed_discriminant':str(D),'local_data':local,'unprocessed_cofactor':str(remaining),'conductor_upper_bound':str(upper),
             'exact_conductor':remaining==1,'upper_bound_below_listed_benchmark':bool(upper<benchmark),'upper_to_benchmark_ratio':str(QQ(upper)/benchmark)}
        result['rows'].append(row);checkpoint(OUT,result)
        print('CONDUCTOR BOUND',row['family'],row['parameter'],'rank',rank,'upper digits',len(str(upper)),'benchmark digits',len(str(benchmark)),'strict improvement',upper<benchmark,flush=True)
    result['status']='COMPLETE_BOUNDED_LOCAL_AUDIT';checkpoint(OUT,result)

if __name__=='__main__':run()
