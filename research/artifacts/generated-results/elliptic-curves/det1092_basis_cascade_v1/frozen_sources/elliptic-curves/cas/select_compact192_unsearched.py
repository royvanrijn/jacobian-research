#!/usr/bin/env python3
"""Fixed deeper point exposure in the previously scored compact H4096 population."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'compact192-unsearched-selection-v1'
POOLS=[LOCAL/n/'result.json' for n in ('r17-retained-extended-primes-v1','r17-discarded-shards-extended-v1','r17-retention512-extended-v1')];PRIOR=ART/'skew8_r17_results_v1.json';ENDPOINTS=ART/'compact_atlas_endpoints_v2.json';GATE=ART/'retention24_r17_results_v1.json'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),Path(cert.__file__).resolve(),*POOLS,PRIOR,ENDPOINTS,GATE]}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve compact192 selection protocol')
    gate=cert.read(GATE)
    if sum(r['rank_lower_bound']>=27 for r in gate['curves'])<2:raise ArithmeticError('productive same-population control required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.compact192-unsearched-selection.v1','sources':sources(),'per_family':32,'maximum_input_rows':6144,'selection_prime_bound':32749,'validation_interval':[32771,65521],'gate':'The same H4096 retention512 population has6144 exactly replayed extended-prime rows. Its preceding24-curve retention cohort produced two new27-point curves. Many high-scoring compact addresses remain unsearched; this finite deeper exposure tests them instead of changing score, extending parameter height or rescanning the population. The previous outcome motivates the budget and is not a rank-density estimate or independent validation.','selection':'Pool all768 original,768 discarded and4608 wider-retained extended scores. Exclude rational isomorphs of all528 prior measured address-equations and all21 separately audited nonsingular endpoints, using equations only. In lexicographic family order choose32 each by unchanged combined S1 through32749, selection-band good count, denominator and signed numerator, excluding exact rational isomorphs of earlier selections. Canonical family order chooses representatives; it is not an independent sampling assertion. Saved validation scores and all measured ranks/points are excluded from ordering. No public catalogue is newly loaded and no failed or repeated point attempt is refilled.','future_point_scope':'After exact selection replay, a separate protocol may freeze192 generic17-only point attempts, all43/49 exact maximum parity classes, at most8640 initial boxes, height125000 and ten seconds per chart, stop at28 pending replay. No adaptive wave or automatic escalation is included.','limits':{'wall_seconds':120,'rss_bytes':1073741824,'maximum_workers':1},'boundaries':'Deeper point exposure within an already finite score-retained population, not a new rank predictor, statistical density law, full parameter coverage, point absence, exact rank, rank upper bound or universal novelty claim.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen compact192 selection inputs differ')
    return p

def j(model):
    v=cert.weierstrass_invariants(tuple(map(cert.F,model)));return v['c4']**3/v['discriminant']

def expected():
    p=protocol();pool=[];seen=set()
    for path in POOLS:
        for r in cert.read(path)['rows']:
            key=r['family'],str(cert.F(r['parameter']))
            if key in seen or max(abs(r['numerator']),r['denominator'])>4096:raise ArithmeticError('unique compact4096 population required')
            seen.add(key);pool.append({**r,'source_score_file':str(path.relative_to(ROOT)),'original_retained_index':r['retained_index']})
    if len(pool)!=6144:raise ArithmeticError('complete6144 saved compact scores required')
    prior=cert.read(PRIOR);old=prior['previous_equations']+[{'address':PRIOR.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in prior['curves']]
    endpoints=[r for r in cert.read(ENDPOINTS)['rows'] if r['status']=='CERTIFIED_SPECIALIZED_SUBGROUP']
    if len(old)!=528 or len(endpoints)!=21:raise ArithmeticError('fixed measured-equation exclusions differ')
    excluded={}
    for r in old+[{'address':'endpoint:'+r['family']+':'+r['endpoint'],'curve':r['curve']} for r in endpoints]:excluded.setdefault(j(r['curve']),[]).append(r)
    rows=[];skips=[];selection={};chosen={}
    for family in sorted({r['family'] for r in pool}):
        ordered=sorted((r for r in pool if r['family']==family),key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']));accepted=[]
        for position,r in enumerate(ordered):
            inv=j(r['model']);matches=[q['address'] for q in excluded.get(inv,[]) if cert.isomorphic(r['model'],q['curve'])];aliases=[q['id'] for q in chosen.get(inv,[]) if cert.isomorphic(r['model'],q['model'])]
            if matches or aliases:skips.append({'family':family,'parameter':r['parameter'],'score_order':position+1,'previous_matches':matches,'selected_aliases':aliases});continue
            q={**r,'id':family+f'-{len(accepted)+1:03}','retained_index':len(accepted)+1,'score_order':position+1};accepted.append(q);rows.append(q);chosen.setdefault(inv,[]).append(q)
            if len(accepted)==p['per_family']:break
        if len(accepted)!=32:raise ArithmeticError('insufficient fixed unsearched population')
        selection[family]=[r['retained_index'] for r in accepted]
    if len(rows)!=192:raise ArithmeticError('fixed192 cohort differs')
    return {'schema':'elliptic-curves.compact192-unsearched-selection-result.v1','status':'COMPLETE_FROZEN_UNSEARCHED_SELECTION','protocol_hash':digest(p),'input_rows':6144,'previous_address_equations':528,'endpoint_exclusions':21,'rows':rows,'selection':selection,'skipped_before_roster_completion':skips,'claim_boundary':p['boundaries']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected();out=D/'result.json'
        if a.stage=='replay':
            if cert.read(out)!=json.loads(json.dumps(d)):raise ArithmeticError('compact192 selection replay differs')
        else:
            if out.exists():raise FileExistsError('preserve compact192 selected cohort')
            checkpoint(out,d)
        print('EXACT192 UNSearched COMPACT COHORT',len(d['skipped_before_roster_completion']),'prior/alias skips',flush=True)
