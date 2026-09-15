#!/usr/bin/env sage -python
"""Order regression and a retained-cloud Sage-alarm fault injection; no search."""
from fractions import Fraction as F
import time

from cancellation_basis_early import CAS, ROOT, OUT, RAW, PRIOR, read, sha, need, new_write


def check():
    from cysignals.signals import AlarmInterrupt
    from cancellation_basis_early_order import visit_order
    import cancellation_basis_early_epoch as worker
    import pari_pointed_backend as backend
    start = time.process_time()
    rows = [{'representative':w} for w in ([1,0,0,0],[0,1,1,0],[1,0,0,1],[1,1,1,1])]
    need(visit_order(rows,2,3,True)==[2,3,1,0], 'latest/earlier/original order differs')
    need(visit_order(rows,2,2,True)==[1,2,3,0], 'two-generator block omitted')
    need(visit_order([{'representative':[1,0]},{'representative':[0,1]}],2,2,False)==[0,1], 'initial order changed')
    rejected = 0
    try: visit_order(rows,2,4,True)
    except ValueError: rejected += 1
    need(rejected==1, 'unchanged basis called a refresh')
    case=next(c for c in read(PRIOR/'inputs.json') if c['id']=='6f5afe31891f7dc1938b')
    origin=ROOT/'artifacts/local/elliptic-curves/cancellation-basis-amplification-v1/arms'/case['id']/'basis_refresh'
    result=read(origin/'result.json'); prefix=[]
    for item in result['calls']:
        need(sha(origin/item['file'])==item['sha256'], 'retained prefix changed')
        call=read(origin/item['file']);prefix.append(call)
        if call['gain']:break
    need(prefix[-1]['gain']>0 and all(c['epoch']==0 for c in prefix), 'missing first full cloud')
    plan=read(OUT/'design.json');plan['gp_sha256']=sha(__import__('pathlib').Path('/usr/bin/gp'))
    new_write(RAW/'alarm-regression-recipe.json', {'case':case,'plan':plan,'prefix_sha256':[r['sha256'] for r in result['calls'][:len(prefix)]],
        'injection':'Replay the literal completed first-cloud transcript; raise Sage AlarmInterrupt on the next rebuild. No PARI point search is called.',
        'source_sha256':sha(CAS/'check_cancellation_basis_early.py')})
    saved_execute=backend.execute;saved_rebuild=worker.rebuild;used=0;injected=0
    def replay_execute(search,mapping,height,seconds,expected_gp_hash):
        nonlocal used
        call=prefix[used];used+=1
        need(mapping==call['mapping'] and height==call['search']['height_bound'], 'retained prefix no longer matches')
        need(expected_gp_hash==call['search']['gp_binary_sha256'], 'GP binding differs')
        points=backend.replay(search,mapping,call['search'])
        return call['search'],points
    def interrupted_rebuild(packet,case,plan,folder):
        nonlocal injected
        if len(packet['points'])>18:
            injected+=1;raise AlarmInterrupt('retained-cloud bank-boundary fault injection')
        return saved_rebuild(packet,case,plan,folder)
    backend.execute=replay_execute;worker.rebuild=interrupted_rebuild
    try:
        worker.run(case,'basis_refresh',plan,RAW/'alarm-regression','DEVELOPMENT_ALARM_REPLAY:'+sha(RAW/'alarm-regression-recipe.json'))
    finally:
        backend.execute=saved_execute;worker.rebuild=saved_rebuild
    out=read(RAW/'alarm-regression/result.json');proof=read(RAW/'alarm-regression/independent-verification.json')
    need(injected==1 and used==len(prefix), 'fault was not reached through the complete cloud')
    need(out['status']=='BANK_UNKNOWN' and out['unknowns']==1 and out['new_directions']==prefix[-1]['gain'], 'certified prefix lost')
    need(proof['status']=='PASS' and proof['rank']['rank']==out['rank_lower_bound'], 'failed bank prevented final rank proof')
    row={'status':'PASS','order_cases':4,'alarm_injections':injected,'retained_calls_replayed':used,
        'rank_lower_bound_preserved':out['rank_lower_bound'],'cpu_seconds':time.process_time()-start,
        'result_sha256':sha(RAW/'alarm-regression/result.json'),'verification_sha256':sha(RAW/'alarm-regression/independent-verification.json'),
        'source_sha256':sha(CAS/'check_cancellation_basis_early.py'),'point_search_calls':0,
        'boundary':'Retained point transcripts are exact replay inputs. The injected Sage exception keeps the whole independently certified cloud, marks the bank UNKNOWN and finalizes without a search retry.'}
    new_write(OUT/'order-regressions.json',row);print(row,flush=True)


if __name__=='__main__':check()
