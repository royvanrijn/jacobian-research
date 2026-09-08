#!/usr/bin/env python3
"""Freeze sixty highest unsearched near-finalists from already verified late scores."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'retained-nearcutoff-mw16-v1'
SCORES=LOCAL/'corrected-mw16-higher-scores-v1'
ORIGINAL=ART/'corrected_mw16_higher_selection_v1.json'
MATCHED=ART/'retained_mw16_score_strata_selection_v1.json'
OUT=ART/'retained_mw16_nearcutoff_selection_v1.json'


def completion_gate():
    portable=cert.read(ART/'strata60_mw16_point_portable_replay_v1.json')
    audit=cert.read(ART/'retained_mw16_score_strata_accounting_replay_v1.json')
    if portable['status']!='PASS' or portable['logical_stages']!=240 or audit['status']!='PASS':raise ArithmeticError('completed independently checked comparison required')


def sources():
    paths=[Path(__file__).resolve(),Path(cert.__file__),SCORES/'protocol.json',SCORES/'result.json',SCORES/'controller/ledger.json',ORIGINAL,MATCHED]
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve next retained-pool design')
    completion_gate()
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-nearcutoff-mw16-protocol.v1','sources':sources(),
        'blocks':10,'curves_per_block':6,'late_rank_window':[7,32],'maximum_curves':60,
        'point_exposure':{'generic_rank':16,'charts_per_curve':43,'height':125000,'seconds_per_chart':10,
                          'seconds_per_curve':600,'workers':2,'rank_stop':None,'maximum_boxes':2580},
        'motivation':'The completed matched comparison favours the strong initial-score tail for certified-direction yield and finds a new25 curve already ninth by the existing extended score. This motivates a separate same-size retained near-finalist discovery trial; it does not establish late-score optimality or a record-rate claim.',
        'selection':'Use saved corrected scores through65521 only. Per family/band rank the original1024 cases by descending late units, good count, denominator and signed numerator; among ranks7..32 choose the first six rational-isomorphism classes not in the two prior frozen prospective60 allocations or this new roster. No point outcome, public catalogue, target or withheld-prime value enters the selection function. No pool or rank-window expansion on failure.',
        'execution_order':'Interleave the six selected positions across the ten lexically ordered family/band blocks.',
        'scope':'A new fixed60 discovery trial from existing verified scalar candidates, not a continuation or rule change of either completed experiment. No parameter enumeration or new prime traces. All maps and independent16 baselines precede every point search. Fixed43 generic charts per curve, identical budgets, no rank stop, adaptive wave, retry, refill or automatic following campaign. Point gains and measured completion are certified independently; novelty only after terminal proofs.'})


def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen nearcutoff inputs changed')
    return p


def key(model):
    v=cert.weierstrass_invariants(tuple(map(cert.F,model)))
    if not v['discriminant']:raise ArithmeticError('singular retained equation')
    return v['c4']**3/v['discriminant']


def expected():
    p=protocol();completion_gate();original=cert.read(ORIGINAL);matched=cert.read(MATCHED)
    old=cert.read(SCORES/'protocol.json');scores=cert.read(SCORES/'result.json')
    if cert.read(SCORES/'controller/ledger.json')['status']!='PASS' or scores['status']!='COMPLETE_FROZEN10240' or len(old['rows'])!=10240 or len(scores['rows'])!=10240:raise ArithmeticError('complete original scalar trial required')
    if original['protocol_sha256']!=cert.hashed(SCORES/'protocol.json') or original['scores_sha256']!=cert.hashed(SCORES/'result.json'):raise ArithmeticError('scalar proof binding differs')
    if original['status']!='PASS_FROZEN60_SELECTION' or matched['status']!='PASS_FROZEN60_MATCHED_SELECTION':raise ArithmeticError('prior prospective allocations required')
    prior={}
    for row in original['selected']+matched['selected']:
        prior.setdefault(key(row['model']),[]).append((row['id'],row['model']))
    pools={}
    for r,s in zip(old['rows'],scores['rows']):
        if r['id']!=s['id']:raise ArithmeticError('scalar roster differs')
        pools.setdefault((r['family'],r['band']),[]).append({**r,**s})
    if len(pools)!=10:raise ArithmeticError('ten family/band blocks required')
    blocks=[];chosen_by_block=[];failures=[];new_seen={}
    for block,pool in sorted(pools.items()):
        if len(pool)!=1024:raise ArithmeticError('equal1024-row blocks required')
        pool.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
        chosen=[];covered=[];excluded=[]
        for i,row in enumerate(pool,1):
            j=key(row['model']);matches=[name for name,m in prior.get(j,[]) if cert.isomorphic(m,row['model'])]
            if matches:covered.append({'late_rank':i,'prior_ids':matches})
            if not p['late_rank_window'][0]<=i<=p['late_rank_window'][1] or len(chosen)==6:continue
            if matches:excluded.append({'late_rank':i,'reason':'PRIOR_PROSPECTIVE_ALLOCATION','ids':matches});continue
            if any(cert.isomorphic(row['model'],m) for m in new_seen.get(j,[])):
                excluded.append({'late_rank':i,'reason':'WITHIN_NEW_ROSTER_Q_ISOMORPHISM'});continue
            new={**row,'original_scalar_id':row['id'],'id':f'nearcut-b{block[1]}-{block[0]}-{i:04}',
                'late_rank':i,'retained_rank':row['retained_index']+1,'arm':'nearcutoff','block':f'{block[0]}-b{block[1]}',
                'sign':1 if row['numerator']>0 else -1,'j_invariant':str(j),
                'j_height':str(max(abs(j.numerator),j.denominator)),
                'parameter_height':max(abs(row['numerator']),row['denominator'])}
            chosen.append(new);new_seen.setdefault(j,[]).append(row['model'])
        blocks.append({'family':block[0],'band':block[1],'previously_allocated':covered,'unallocated_scalar_cases':1024-len(covered),
                       'selected_late_ranks':[r['late_rank'] for r in chosen],'excluded_in_window':excluded})
        if len(chosen)!=6:failures.append({'family':block[0],'band':block[1],'selected_count':len(chosen),'reason':'FIXED_WINDOW_INCOMPLETE_NO_POINT_SEARCH'})
        chosen_by_block.append(chosen)
    rows=[pool[i] for i in range(6) for pool in chosen_by_block if i<len(pool)]
    return {'schema':'elliptic-curves.retained-nearcutoff-mw16-selection.v1',
        'status':'PASS_FROZEN60_RETAINED_SELECTION' if not failures and len(rows)==60 else 'SELECTION_INCOMPLETE_NO_POINT_SEARCH',
        'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'selected':rows,'blocks':blocks,
        'selection_failures':failures,'prior_prospective_allocations':120,
        'original_scalar_cases':10240,'covered_scalar_cases':sum(len(r['previously_allocated']) for r in blocks),
        'claim_boundary':p['scope']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','select','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=d:raise ArithmeticError('retained nearcutoff selection replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve next retained allocation')
            checkpoint(OUT,d)
        print('RETAINED NEARCUTOFF',d['status'],'COVERED',d['covered_scalar_cases'],'OF10240',flush=True)
