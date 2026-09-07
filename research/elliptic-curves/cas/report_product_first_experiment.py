#!/usr/bin/env python3
"""Exact aggregate of product-first retention, cached traces and completed point exposure."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'productfirst24-r17-pari-v1';OUT=ART/'product_first_experiment_v1.json'
def expected():
    paths={Path(__file__).resolve()}
    def read(p):paths.add(p);return cert.read(p)
    population=read(LOCAL/'higher32768-product-first-v1/result.json');extended=read(LOCAL/'higher32768-product-first-extended-v1/result.json');old=read(LOCAL/'higher-r17-product-score-v1/result.json');p=read(D/'protocol.json');v=read(D/'verification-ledger.json');o=read(D/'odd-cloud-audit/ledger.json');proof=read(ART/'productfirst24_r17_results_v1.json');sc=read(ART/'higher_displayed_reduction_scalings_v1.json')
    for folder,name in [(LOCAL/'higher32768-product-first-v1','replay'),(LOCAL/'higher32768-product-first-extended-v1','replay'),(D/'post-batch','certify'),(D/'post-batch','proof-replay'),(D/'post-batch','geometry')]:
        s=read(folder/(name+'.supervisor.json'))
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required exact replay failed or censored')
    if population['status']!='COMPLETE_FIXED_PRODUCT_FIRST_POPULATION' or extended['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or len(population['rows'])!=6144 or len(extended['rows'])!=6144 or v['status']!='PASS' or o['status']!='PASS':raise ArithmeticError('complete population and point audits required')
    keyed={(r['family'],r['parameter']):r for r in old['rows']};overlap=0
    for r in extended['rows']:
        previous=keyed.get((r['family'],r['parameter']))
        if previous:
            if (r['combined_selection_units'],r['validation_units'],r['combined_good'])!=(previous['product_selection_units'],previous['product_validation_units'],previous['combined_good']):raise ArithmeticError('independent product-score implementations disagree')
            overlap+=1
        if bool(previous)!=bool(r['cached_trace']):raise ArithmeticError('cache membership differs')
    if overlap!=1806 or population['previous_pool_overlap']!=overlap:raise ArithmeticError('retention overlap differs')
    vk={r['id']:r for r in v['rows']};ok={r['id']:r for r in o['rows']};qk={r['id']:r for r in proof['curves']};rows=[]
    if len(qk)!=24 or set(vk)!=set(qk) or set(ok)!=set(qk) or proof['within_batch_isomorphic_pairs']:raise ArithmeticError('fixed distinct24 proof roster differs')
    for row in p['rows']:
        name=row['id'];raw=read(D/name/'result.json');q=qk[name];odd=ok[name];cloud=read(ROOT/vk[name]['cloud_certificate']);l=read(ROOT/odd['output']);bound=q['rank_lower_bound']
        if (row['family'],row['parameter']) in keyed or q['icarm_matches'] or q['previous_matches']:raise ArithmeticError('new-retention or catalogue exclusions differ')
        if cert.hashed(D/name/'result.json')!=vk[name]['input_sha256'] or cert.hashed(ROOT/vk[name]['cloud_certificate'])!=vk[name]['cloud_sha256'] or cert.hashed(ROOT/odd['output'])!=odd['output_sha256']:raise ArithmeticError('bound point proof changed')
        if raw['status']!='COMPLETE_DECLARED_POINT_ATTEMPT' or raw['rank_lower_bound']!=bound or cloud['rank_lower_bound']!=bound or odd['mod2_lower_bound']!=bound or {a['modulus']:a['finite_column_rank'] for a in l['audits']}!={3:bound,5:bound}:raise ArithmeticError('full-cloud bounds differ')
        if len(raw['charts'])!=len(p['generic_masks'][row['family']]) or any(c['search']['status']!='bounded_search_complete' or c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('uniform completed coverage differs')
        rows.append({'id':name,'family':row['family'],'parameter':row['parameter'],'rank_lower_bound':bound,'completed_boxes':len(raw['charts']),'retained_points':len(cloud['points'])})
    if sum(r['completed_boxes'] for r in rows)!=1080 or Counter(r['rank_lower_bound'] for r in rows)!=Counter({17:23,18:1}):raise ArithmeticError('completed cohort result differs')
    if sc['status']!='PASS_EXACT_SCALING_AUDIT' or sc['address_models']!=10482 or sc['status_counts']!={'SCALED_DISPLAY_STILL_SINGULAR':766} or sc['prime_counts']!={'13':766}:raise ArithmeticError('displayed-reduction audit differs')
    scaled={(r['family'],r['parameter']) for r in sc['rows']};cohorts={}
    for name in ('higher24_r17_results_v1.json','product22_r17_results_v1.json','productfirst24_r17_results_v1.json'):
        cohort=read(ART/name)['curves'];cohorts[name]={'curves':len(cohort),'models_with_13_scaling':sum((r['family'],r['parameter']) in scaled for r in cohort)}
    return {'schema':'elliptic-curves.product-first-experiment.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'primitive_addresses':sum(r['summary'][3] for r in population['shards']),'small_exhaustive_scores':sum(r['summary'][3] for r in population['small_checks']),'retained_addresses':6144,'overlap_and_cross_implementation_agreement':overlap,'newly_retained_addresses':4338,'fresh_extension_traces':4338*5978,'reused_extension_traces':1806*5978,'distinct_point_curves':24,'completed_boxes':1080,'rows':rows,'catalogue_snapshot_curves':proof['catalogue']['curve_count'],'previous_address_equations':len(proof['previous_equations']),'catalogue_or_previous_matches':0,'local_scaling_audit':{'retained_union':10482,'scaled_models':766,'prime':13,'recovered_good_displays':0,'cohort_counts':cohorts},'claim_boundary':'Fixed product-first retention and longer product scoring expose24 previously discarded addresses, without a high-rank inventory addition. These lower bounds do not establish true-rank distributions, optimal selectors, saturation or rank upper bounds. The13-scaling cohort contrast is descriptive and does not establish a cause or a calibrated predictor. The whole122368792-address population was not audited for local model scaling. No new near-record or universal novelty.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('product-first aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve product-first report')
        checkpoint(OUT,d)
    print('PRODUCT-FIRST REPORT PASS:4338 newly retained;1080 completed boxes;23x17 and1x18',flush=True)
