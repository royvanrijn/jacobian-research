#!/usr/bin/env python3
"""Bind completed paired selection, search, quotient and visibility evidence."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2]; ART=ROOT/'artifacts/generated-results/elliptic-curves'; LOCAL=ROOT/'artifacts/local/elliptic-curves'; D=LOCAL/'product24-r17-pari-v1'; OUT=ART/'product22_comparison_v1.json'
def expected():
    paths={Path(__file__).resolve()}
    def read(p):
        paths.add(p); return cert.read(p)
    selection=read(LOCAL/'higher-r17-product-score-v1/result.json'); p=read(D/'protocol.json'); v=read(D/'verification-ledger.json'); odd=read(D/'odd-cloud-audit/ledger.json'); proof=read(ART/'product22_r17_results_v1.json')
    if selection['status']!='COMPLETE_FIXED_PAIRED_SCORE_COMPARISON' or len(selection['rows'])!=6144 or v['status']!='PASS' or odd['status']!='PASS': raise ArithmeticError('complete selection and quotient replays required')
    for name in ('certify','proof-replay','geometry'):
        s=read(D/'post-batch'/(name+'.supervisor.json'))
        if s['outcome']!='completed' or s['returncode']!=0: raise ArithmeticError('post-batch exact replay failed')
    keyed={r['id']:r for r in proof['curves']}; vk={r['id']:r for r in v['rows']}; ok={r['id']:r for r in odd['rows']}; rows=[]
    if len(keyed)!=22 or set(keyed)!=set(vk) or set(keyed)!=set(ok) or proof['within_batch_isomorphic_pairs']: raise ArithmeticError('fixed distinct22 roster differs')
    chosen={f"{r['family']}-{r['retained_index']:03}":r for r in selection['prospective_candidates']}
    if set(chosen)!=set(keyed): raise ArithmeticError('paired score selection differs')
    for row in p['rows']:
        name=row['id']; raw=read(D/name/'result.json'); q=keyed[name]; o=ok[name]; c=read(ROOT/vk[name]['cloud_certificate']); l=read(ROOT/o['output'])
        if row['parameter']!=chosen[name]['parameter'] or row['arms']!=chosen[name]['arms'] or q['icarm_matches'] or q['previous_matches']: raise ArithmeticError('address, arms or catalogue novelty differs')
        if cert.hashed(D/name/'result.json')!=vk[name]['input_sha256'] or cert.hashed(ROOT/vk[name]['cloud_certificate'])!=vk[name]['cloud_sha256'] or cert.hashed(ROOT/o['output'])!=o['output_sha256']: raise ArithmeticError('bound evidence changed')
        if raw['status']!='COMPLETE_DECLARED_POINT_ATTEMPT' or vk[name]['status']!='PASS' or o['status']!='PASS': raise ArithmeticError('complete point attempt required')
        if any(r['search']['status']!='bounded_search_complete' or r['search']['height_bound']!=125000 or r['search']['timeout_seconds']!=10 for r in raw['charts']): raise ArithmeticError('uniform completed coverage differs')
        bound=q['rank_lower_bound']
        if raw['rank_lower_bound']!=bound or o['mod2_lower_bound']!=bound or any(a['finite_column_rank']!=bound for a in l['audits']): raise ArithmeticError('all quotient lower bounds differ')
        rows.append({'id':name,'family':row['family'],'parameter':row['parameter'],'arms':row['arms'],'completed_boxes':len(raw['charts']),'rank_lower_bound':bound,'retained_points':len(c['points']),'odd_lower_bounds':o['odd_lower_bounds']})
    if sum(r['completed_boxes'] for r in rows)!=988 or Counter(r['rank_lower_bound'] for r in rows)!=Counter({17:21,18:1}): raise ArithmeticError('fixed completed result differs')
    visibility=read(ART/'generic_point_box_visibility_v1.json'); control=read(ART/'native28_generic_visibility_v1.json')
    if visibility['status']!='COMPLETE_EXACT_VISIBILITY_AUDIT' or visibility['discrepancy_count'] or control['status']!='PASS_EXACT_OBSERVATIONS' or control['discrepancies'] or control['generic_sections_with_a_visible_sign']!=0 or control['independently_certified_recovered_lower_bound']!=28: raise ArithmeticError('exact visibility diagnostic differs')
    arms={a:{'curves':sum(a in r['arms'] for r in rows),'lower_bound_histogram':dict(sorted(Counter(str(r['rank_lower_bound']) for r in rows if a in r['arms']).items()))} for a in ('product','original_s1')}
    if any(a['curves']!=12 for a in arms.values()): raise ArithmeticError('paired12-versus12 roster differs')
    return {'schema':'elliptic-curves.product22-comparison.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in sorted(paths)},'scored_addresses':6144,'distinct_curves':22,'completed_boxes':988,'height':125000,'seconds_per_chart':10,'rows':rows,'arms':arms,'shared_curves':[r['id'] for r in rows if len(r['arms'])==2],'catalogue_comparison':{'snapshot_curves':proof['catalogue']['curve_count'],'previous_address_equations':len(proof['previous_equations']),'matches':0},'generic_visibility':{'observations':visibility['observations_checked'],'cohorts':visibility['cohort_summary'],'omitted_in_box_points':0,'positive_control':{'observations':control['point_chart_observations'],'visible_generic_representatives':0,'certified_recovered_lower_bound':28}},'claim_boundary':'The saved6144 pool was already truncated by short-prime S1. Each predeclared score arm has eleven bounds17 and one18, with the same18 curve in both arms. These are bounded detection outcomes, not true-rank comparisons, calibrated probabilities, saturation or upper bounds. The known28 control has zero visible original generic representatives: their direct visibility is not a valid diagnostic for exceptional-direction failure. No new near-record curve, universal novelty or proof of optimal selection.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d: raise ArithmeticError('paired report differs')
    else:
        if OUT.exists(): raise FileExistsError('preserve paired report')
        checkpoint(OUT,d)
    print('PAIRED22 REPORT PASS:988 completed boxes; both arms 11x17 +1x18; visibility countercheck retained',flush=True)
