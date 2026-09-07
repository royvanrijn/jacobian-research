#!/usr/bin/env python3
"""Bind two frozen local-feature cohorts, the new25 proof and its completed follow-up."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUT=ART/'local_feature_experiments_v1.json'

def expected():
    paths={Path(__file__).resolve()}
    def read(p):paths.add(p);return cert.read(p)
    def terminal(p):
        s=read(p)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required terminal replay failed or censored')
    selection=read(LOCAL/'higher13-scaled-selection-v1/result.json')
    nodes=read(LOCAL/'higher-split-node-score-v1/result.json')
    for name in ('higher13-scaled-selection-v1','higher-split-node-score-v1'):
        terminal(LOCAL/name/'replay.supervisor.json')
    if selection['status']!='COMPLETE_FROZEN_CONDITIONAL_SELECTION' or (selection['input_rows'],selection['excluded_addresses'],selection['eligible_rows'])!=(6144,70,454):raise ArithmeticError('conditional selection differs')
    if nodes['status']!='COMPLETE_FIXED_SPLIT_NODE_SCORING' or len(nodes['rows'])!=10482 or len(nodes['prospective_candidates'])!=24 or sum(r['excluded_prior_address'] for r in nodes['rows'])!=94:raise ArithmeticError('fixed split-node scoring differs')
    cohorts=[]
    for folder,prefix,histogram,previous in [('scaled13_24-r17-pari-v1','scaled13_24_r17',{17:13,18:4,19:4,20:2,25:1},472),('splitnode24-r17-pari-v1','splitnode24_r17',{17:23,18:1},496)]:
        d=LOCAL/folder;p=read(d/'protocol.json');v=read(d/'verification-ledger.json');o=read(d/'odd-cloud-audit/ledger.json');proof=read(ART/(prefix+'_results_v1.json'))
        for label in ('certify','proof-replay','geometry'):terminal(d/'post-batch'/(label+'.supervisor.json'))
        if v['status']!='PASS' or o['status']!='PASS' or proof['within_batch_isomorphic_pairs'] or len(proof['previous_equations'])!=previous:raise ArithmeticError('complete distinct cohort required')
        vk={r['id']:r for r in v['rows']};ok={r['id']:r for r in o['rows']};qk={r['id']:r for r in proof['curves']};rows=[]
        if len(qk)!=24 or set(vk)!=set(qk) or set(ok)!=set(qk) or {r['id'] for r in p['rows']}!=set(qk):raise ArithmeticError('fixed24 roster differs')
        for row in p['rows']:
            name=row['id'];raw=read(d/name/'result.json');q=qk[name];odd=ok[name];cloud=read(ROOT/vk[name]['cloud_certificate']);ell=read(ROOT/odd['output']);bound=q['rank_lower_bound']
            if q['icarm_matches'] or q['previous_matches'] or (q['family'],q['parameter'])!=(row['family'],row['parameter']):raise ArithmeticError('curve identity or catalogue exclusion differs')
            if cert.hashed(d/name/'result.json')!=vk[name]['input_sha256'] or cert.hashed(ROOT/vk[name]['cloud_certificate'])!=vk[name]['cloud_sha256'] or cert.hashed(ROOT/odd['output'])!=odd['output_sha256']:raise ArithmeticError('bound point proof changed')
            if raw['status']!='COMPLETE_DECLARED_POINT_ATTEMPT' or any(b!=bound for b in (raw['rank_lower_bound'],cloud['rank_lower_bound'],odd['mod2_lower_bound'])) or {a['modulus']:a['finite_column_rank'] for a in ell['audits']}!={3:bound,5:bound}:raise ArithmeticError('full-cloud bounds differ')
            if len(raw['charts'])!=len(p['generic_masks'][row['family']]) or any(c['search']['status']!='bounded_search_complete' or c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('uniform completed coverage differs')
            rows.append({'id':name,'family':row['family'],'parameter':row['parameter'],'rank_lower_bound':bound,'completed_boxes':len(raw['charts']),'retained_points':len(cloud['points'])})
        if sum(r['completed_boxes'] for r in rows)!=1080 or Counter(r['rank_lower_bound'] for r in rows)!=Counter(histogram):raise ArithmeticError('fixed cohort result differs')
        cohorts.append({'cohort':prefix,'rows':rows,'completed_boxes':1080,'rank_lower_bound_counts':dict(sorted(histogram.items())),'catalogue_snapshot_curves':proof['catalogue']['curve_count'],'previous_address_equations':previous,'catalogue_or_previous_matches':0})
    coverage=read(ART/'scaled13_rank25_11952_300_adaptive_coverage_v1.json');follow=read(ART/'scaled13_25_followup_results_v1.json');minimal=read(ART/'scaled13_24_rank25_minimal_proof_v1.json')
    for label in ('replay','cloud-audit','geometry','certify','certify-check'):terminal(LOCAL/'scaled13-rank25-11952-300-adaptive-v1'/(label+'.supervisor.json'))
    if (coverage['initial_charts'],coverage['adaptive_declared_charts'],coverage['adaptive_attempted_charts'],coverage['adaptive_completed_boxes'],coverage['retained_point_count'],coverage['mod2_lower_bound'])!=(49,301,301,301,3260,25) or coverage['odd_modulus_lower_bounds']!={'3':25,'5':25}:raise ArithmeticError('completed follow-up differs')
    if len(follow['curves'])!=1 or len(minimal['curves'])!=1 or minimal['status']!='PASS':raise ArithmeticError('single new25 proof required')
    m=minimal['curves'][0];f=follow['curves'][0]
    if any((r['id'],r['parameter'],r['rank_lower_bound'])!=('11952-300','102/1525',25) or r['icarm_matches'] or r['previous_matches'] for r in (m,f)):raise ArithmeticError('new25 identity differs')
    inventory=read(ART/'new_high_rank_curve_index_v12.json');replay=read(ART/'new_high_rank_curve_index_v12_memory_replay_v1.json');incidence=read(ART/'inventory101_incidence_v1.json');new=next(r for r in inventory['curves'] if r['id']=='new-20260906-101')
    if len(inventory['curves'])!=101 or replay['status']!='PASS' or replay['curves_checked']!=101 or new['source_certificate']!='scaled13_25_followup_results_v1.json' or new['rank_lower_bound']!=25 or incidence['status']!='PASS' or incidence['pairs_checked']!=1212:raise ArithmeticError('new101 inventory/incidence binding differs')
    return {'schema':'elliptic-curves.local-feature-experiments.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'conditional_selection':{'input_rows':6144,'excluded_addresses':70,'eligible_rows':454,'selected_rows':24},'split_node_score':{'address_models':10482,'excluded_addresses':94,'selected_rows':24,'bonus_constant':[7,5]},'cohorts':cohorts,'new_curve':{'id':new['id'],'family':m['family'],'parameter':m['parameter'],'rank_lower_bound':25,'minimal_curve':m['minimal_curve'],'minimal_model_proof':'scaled13_24_rank25_minimal_proof_v1.json','sage_export':'new_scaled13_rank25_curve_11952.sage'},'followup':{'completed_boxes':301,'initial_boxes':49,'all_retained_points':3260,'mod2_lower_bound':25,'odd_modulus_lower_bounds':coverage['odd_modulus_lower_bounds']},'inventory_curves':101,'inventory_rank_lower_bound_counts':dict(sorted(Counter(r['rank_lower_bound'] for r in inventory['curves']).items())),'claim_boundary':'One additional catalogue-unmatched curve with25 exactly independent rational points. The conditional13-scaling experiment and the fixed split-node bonus experiment measure two distinct bounded policies; they do not prove selector superiority or a causal model-size effect. All2160 initial boxes and301 adaptive boxes completed, with no rank upper bound or saturation conclusion. Twelve-presentation incidence is not universal novelty. The strongest new inventory lower bound remains27, and a new near-record remains open.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=__import__('json').loads(__import__('json').dumps(d)):raise ArithmeticError('local-feature aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve local-feature report')
        checkpoint(OUT,d)
    print('LOCAL FEATURES PASS: one new25;101 curves;2160 initial and301 adaptive boxes',flush=True)
