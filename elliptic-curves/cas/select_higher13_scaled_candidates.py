#!/usr/bin/env python3
"""Fixed S1 selection conditioned on an exact smaller integral model at13."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'higher13-scaled-selection-v1';PARENT=LOCAL/'higher32768-r17-extended-v1';SCALING=ART/'higher_displayed_reduction_scalings_v1.json';EXCLUDED=[LOCAL/n/'protocol.json' for n in ('higher24-r17-pari-v1','product24-r17-pari-v1','productfirst24-r17-pari-v1')]
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),PARENT/'protocol.json',PARENT/'result.json',SCALING,*EXCLUDED)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve model-size conditional protocol')
    gate=ART/'product_first_portable_replay_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('complete product-first experiment required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher13-scaled-selection.v1','sources':sources(),'gate_sha256':cert.hashed(gate),'prime':13,'scale_exponent':1,'selection_prime_bound':32749,'validation_interval':[32771,65521],'per_family':4,'gate':'The exact scaling audit exhibits766 smaller integral models in the union of retained pools. The original higher24 cohort included eight of them, including the curve that later certified26; the paired and product-first cohorts included one and zero. This is not evidence of causation. Test this measured model-size difference directly on the already saved original S1 pool, retaining its score and the same point policy rather than changing scores and model-size conditioning together.','selection':'From the original6144 S1-extended addresses, retain only exact13-scaled models in the audited certificate. Exclude the union of all addresses in the three completed higher-height point protocols, uniformly and without reading their ranks or points. Within each family select four by the unchanged S1 score through32749, good count descending, denominator and signed numerator. Validation primes never break ties. No public curve or record parameter, rank label, new trace, new parameter scan or result-dependent refill.','geometry':'A=13^4*a and B=13^6*b gives the explicit Q-isomorphism x=13^2*X,y=13^3*Y. The displayed coefficients shrink exactly, without sacrificing the inherited17-point subgroup. This says nothing about global minimality, conductor, rank incidence or exceptional-point visibility.','future_point_scope':'After exact selection replay, at most24 generic17-only attempts with the original43/49 maximum-parity policy, height125000 and ten seconds per chart, under a separate frozen point protocol.','limits':{'maximum_input_rows':6144,'maximum_selected_rows':24,'wall_seconds':120,'rss_bytes':536870912,'workers':1},'boundaries':'A bounded conditional model-size experiment, not a demonstrated rank predictor or optimality claim. It cannot establish rank upper bounds, exact ranks, saturation or universal novelty.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen conditional selector changed')
    return p
def expected():
    p=protocol();parent=cert.read(PARENT/'result.json');sc=cert.read(SCALING)
    if parent['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or sc['status']!='PASS_EXACT_SCALING_AUDIT':raise ArithmeticError('complete saved scores and exact scaling required')
    scaled={(r['family'],r['parameter']):r for r in sc['rows'] if r['prime']==13 and r['scale_exponent']==1};excluded={(r['family'],r['parameter']) for path in EXCLUDED for r in cert.read(path)['rows']};eligible=[]
    for r in parent['rows']:
        key=r['family'],r['parameter']
        if key not in scaled or key in excluded:continue
        q=scaled[key];A,B=map(int,r['model'][3:])
        if A!=13**4*int(q['scaled_A']) or B!=13**6*int(q['scaled_B']):raise ArithmeticError('exact smaller model differs')
        eligible.append({**r,'local_scaling':{'prime':13,'exponent':1,'scaled_model':['0','0','0',q['scaled_A'],q['scaled_B']]}})
    selection={}
    for f in sorted({r['family'] for r in parent['rows']}):
        pool=sorted((r for r in eligible if r['family']==f),key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        if len(pool)<4:raise ArithmeticError('insufficient fixed conditional candidates')
        selection[f]=[r['retained_index'] for r in pool[:4]]
    return {'schema':'elliptic-curves.higher13-scaled-selection-result.v1','status':'COMPLETE_FROZEN_CONDITIONAL_SELECTION','protocol_hash':digest(p),'input_rows':len(parent['rows']),'excluded_addresses':len(excluded),'eligible_rows':len(eligible),'rows':eligible,'selection':selection,'claim_boundary':p['boundaries']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        out=D/'result.json';d=expected()
        if a.stage=='replay':
            if cert.read(out)!=d:raise ArithmeticError('conditional selector replay differs')
        else:
            if out.exists():raise FileExistsError('preserve conditional result')
            checkpoint(out,d)
        print('EXACT13-CONDITIONAL SELECTION',d['eligible_rows'],'eligible;24 selected',d['selection'],flush=True)
