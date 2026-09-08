#!/usr/bin/env python3
"""Retrospective coefficient/point-exposure comparison, with no rank prediction."""
import argparse
from pathlib import Path
from collections import Counter
from math import gcd
from fractions import Fraction
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as minimal
import certify_higher24_r17_results as cohort
from research_runtime.store import checkpoint
ROOT=cohort.ROOT;ART=cohort.ART;OUT=ART/'higher24_visibility_cost_v1.json'
def size(model):
    model=tuple(map(cert.F,model));inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
    short=(cert.F(0),)*3+(-inv['c4']/48,-inv['c6']/864)
    try:normal=minimal.integral(short)
    except ArithmeticError:normal=None
    result={'j_projective_height_bits':max(abs(j.numerator),j.denominator).bit_length(),'j_invariant':str(j),'normalized_integral_curve':list(map(str,normal)) if normal else None,'normalized_coefficient_bits':max(abs(int(x)).bit_length() for x in normal) if normal else None}
    if normal is None:result['minimality_status']='NO_INTEGRAL_NORMALIZATION_AT_CURRENT_SCALE';return result
    result['invariant_gcd']=gcd(abs(int(inv['c4'])),abs(int(inv['c6'])))
    try:result['minimality']=minimal.minimality(normal);result['minimality_status']='PROVED_GLOBAL_MINIMAL'
    except ArithmeticError as e:result['minimality_status']='UNRESOLVED_WITHIN_EXISTING_CHEAP_GATE';result['minimality_gate_message']=str(e)
    return result

def statistics(rows,key):
    values=sorted(r[key] for r in rows if r.get(key) is not None)
    return {'count':len(values),'minimum':values[0] if values else None,'median':str(Fraction(values[(len(values)-1)//2]+values[len(values)//2],2)) if values else None,'maximum':values[-1] if values else None}
def expected():
    if cert.read(cohort.D/'verification-ledger.json')['status']!='PASS':raise ArithmeticError('all higher24 point replays required before retrospective comparison')
    paths=[Path(__file__).resolve(),Path(minimal.__file__).resolve(),cohort.OUT,ART/'retention24_r17_results_v1.json',ART/'public_compact_parameter_heights_v1.json',cohort.D/'verification-ledger.json'];groups={}
    for name,path in [('prior_h4096',ART/'retention24_r17_results_v1.json'),('higher_h32768_slices',cohort.OUT)]:
        d=cert.read(path);rows=[]
        for r in d['curves']:
            rawpath=ROOT/r['discovery_witness']['path'];cloudpath=ROOT/r['complete_cloud_certificate']['path'];raw=cert.read(rawpath);cloud=cert.read(cloudpath);paths += [rawpath,cloudpath]
            if cert.hashed(rawpath)!=r['discovery_witness']['sha256'] or cert.hashed(cloudpath)!=r['complete_cloud_certificate']['sha256'] or cloud['input_sha256']!=cert.hashed(rawpath):raise ArithmeticError('point proof binding differs')
            H={c['search']['height_bound'] for c in raw['charts']};T={c['search']['timeout_seconds'] for c in raw['charts']}
            if H!={125000} or T!={10}:raise ArithmeticError('paired exposure budgets differ')
            t=cert.F(r['parameter']);rows.append({'id':r['id'],'family':r['family'],'parameter':r['parameter'],'parameter_height':max(abs(t.numerator),t.denominator),'rank_lower_bound':r['rank_lower_bound'],'retained_points':len(cloud['points']),'charts':len(raw['charts']),'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in raw['charts']),'gp_cpu_ms':sum(c['search']['search_cpu_ms'] or 0 for c in raw['charts']),**size(r['curve'])})
        if len(rows)!=24:raise ArithmeticError('fixed24 cohorts required')
        groups[name]={'rows':rows,'rank_lower_bound_counts':{str(k):v for k,v in sorted(Counter(r['rank_lower_bound'] for r in rows).items())},'statistics':{k:statistics(rows,k) for k in ('parameter_height','j_projective_height_bits','normalized_coefficient_bits','retained_points','gp_cpu_ms')}}
    public=[r for r in cert.read(ART/'public_compact_parameter_heights_v1.json')['rows'] if r['status']=='PASS_EXACT_COMPACT_Q_ISOMORPHISM' and r['current_reported_rank_lower_bound']>=28]
    rows=[{'catalogue_id':r['id'],'reported_rank_lower_bound':r['current_reported_rank_lower_bound'],**size(r['public_curve'])} for r in public]
    if len(rows)!=8:raise ArithmeticError('fixed eight mapped reported>=28 public successes required')
    groups['public_successes_reported_ge28']={'rows':rows,'statistics':{k:statistics(rows,k) for k in ('j_projective_height_bits','normalized_coefficient_bits')}}
    return {'schema':'elliptic-curves.higher24-visibility-cost.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'groups':groups,'claim_boundary':'Retrospective descriptive comparison of two fixed24 candidate cohorts at identical125000/10-second per-chart exposure and eight mapped public successes. Exact j height and normalized coefficient bit sizes; global minimality only where the existing bounded gcd proof closes it. Point ranks are lower bounds from separately replayed proofs, public ranks remain catalogue metadata. Different selected populations and coefficient sizes confound causality: no density, calibrated rank predictor, rank upper bound or proof that misses result from poor visibility. No selection/refill or new point search occurs here.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('coefficient/visibility comparison differs')
    else:
        if OUT.exists():raise FileExistsError('preserve descriptive audit')
        checkpoint(OUT,d)
    for name,r in d['groups'].items():print('HIGHER VISIBILITY COST',name,r.get('rank_lower_bound_counts'),r['statistics'],flush=True)
