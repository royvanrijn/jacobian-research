"""Scheduling, proof promotion, deduplication, budget and recovery regressions."""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sqlite3
import sys

import pytest

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import high_rank_foundry_policy as policy
import high_rank_foundry_intake as intake
import high_rank_foundry_job as worker
import run_high_rank_foundry as control


def curve(rank=25,stale=0,calls=100):
    return {'id':'c','rank':rank,'stale_calls':stale,'total_calls':calls,'batches':1,
            'last_batch_gain':3,'gaining_batches':3,'bank_index':0,'history':[]}


def result(rank=26,calls=100,at=97):
    return {'rank_lower_bound':rank,'calls':calls,'packet_path':'packet.json','packet_sha256':'abc',
            'gain_timeline':[{'call':at,'before':rank-1,'after':rank}]}


def test_stalled_high_rank_cannot_dominate_late_lower_cascade():
    stalled=curve(28,600,6000);stalled['last_batch_gain']=0
    assert policy.utility(curve(24,3))>policy.utility(stalled)
    c=policy.apply_result(curve(27,290),{'rank_lower_bound':27,'calls':100,'packet_path':'p','packet_sha256':'h','gain_timeline':[]})
    assert c['state']=='COOLED' and c['rank']==27


def test_gain_recency_resets_stale_budget_and_does_not_stop_at_28():
    out=policy.apply_result(curve(27,290,500),result(28,100,97))
    assert out['rank']==28 and out['state']=='READY' and out['stale_calls']==3
    assert out['history'][-1]['call']==597
    assert policy.apply_result(curve(31),result(32))['state']=='TARGET_CERTIFIED'


def test_invalid_rank_or_unattributed_gain_never_promoted():
    with pytest.raises(ArithmeticError):policy.apply_result(curve(27),result(26))
    bad=result();bad['gain_timeline']=[]
    with pytest.raises(ArithmeticError):policy.apply_result(curve(),bad)
    with pytest.raises(ArithmeticError):policy.apply_result(curve(),result(at=101))


def test_breadth_is_reserved_by_expensive_time_not_current_rank():
    assert policy.choose_lane(59,41,policy.DEFAULTS)=='fresh'
    assert policy.choose_lane(61,39,policy.DEFAULTS)=='exploit'
    assert policy.choose_lane(100,0,policy.DEFAULTS,has_exploit=False)=='fresh'


def test_batch_work_grows_in_steps_without_a_daily_ceiling():
    config={**policy.DEFAULTS,'daily_point_calls':None,'daily_worker_seconds':None}
    assert policy.effective_batch_calls(config,0)==100
    assert policy.effective_batch_calls(config,99)==100
    assert policy.effective_batch_calls(config,100)==125
    assert policy.effective_batch_calls(config,308)==175
    assert policy.effective_batch_calls(config,2000)==300


def test_parent_refresh_escapes_maximum_span_and_traverses_disjoint_samples():
    maximum=list(range(1,44))
    blocks=[policy.parent_masks('103b2',maximum,i)[0] for i in range(3)]
    assert set.union(*map(set,blocks))==set(maximum)
    assert not (set(blocks[0])&set(blocks[1]))
    a,_=policy.parent_masks('103b2',maximum,3)
    b,_=policy.parent_masks('103b2',maximum,4)
    assert not set(a)&set(b) and not set(a)&set(maximum)
    assert any(m >= 65536 for m in a)
    assert policy.next_bank(curve(stale=100))==1


def test_height_shells_cover_every_reduced_signed_parameter_once():
    found=[(p,q) for h in range(1,9) for p,q in intake.shell(h)]
    assert len(found)==len(set(F(p,q) for p,q in found))
    expected={F(p,q) for q in range(1,9) for p in range(-8,9)}
    assert {F(p,q) for p,q in found}==expected


def test_exact_dedup_distinguishes_quadratic_twists():
    model=['0','0','0','-1','1']
    scaled=['0','0','0','-16','64']
    twist=['0','0','0','-4','8']
    assert intake.jkey(model)==intake.jkey(scaled)==intake.jkey(twist)
    assert intake.matches(model,[scaled])
    assert not intake.matches(model,[twist])


def test_score_independent_order_ignores_score_and_rank_fields():
    pool=[{'family':f,'parameter':str(i),'numerator':i,'denominator':1,
           'combined_selection_units':i,'combined_good':1} for f in policy.FAMILIES for i in range(1,7)]
    before=intake.streams(pool)
    for r in pool:r.update(combined_selection_units=-r['combined_selection_units'],rank=999)
    after=intake.streams(pool)
    for f in policy.FAMILIES:
        assert [r['parameter'] for r in before[f]['all']]==[r['parameter'] for r in after[f]['all']]


def test_conductor_gate_never_claims_a_record_from_discriminant():
    model=['0','0','0','-1','1']
    catalogue={'curves':[{'rank_lower_bound':20,'conductor':'100000'},{'rank_lower_bound':30,'conductor':None}]}
    gate=intake.conductor_gate(model,20,catalogue)
    assert gate['eligible'] and gate['missing_catalogue_conductors']==1
    assert 'exact_conductor' not in gate
    assert intake.conductor_gate(model,31,catalogue)['reason']=='NO_REPORTED_BENCHMARK'


def test_cached_calls_and_cloud_directions_are_not_double_counted():
    terminal={'stages':[{'epoch':0,'charts':80,'before':24,'after':25},
                        {'epoch':1,'charts':27,'before':25,'after':26}]}
    cloud={'gains':[{'chart':'epoch-01/chart-0026.json','after':27}]}
    events=worker.segment_events(terminal,cloud,inherited=100,used=12)
    assert [(e['call'],e['after']) for e in events]==[(19,26),(19,27)]
    with pytest.raises(ArithmeticError):
        worker.segment_events(terminal,{'gains':[{'chart':'epoch-00/chart-0000.json','after':28}]},100,0)


def test_lock_and_process_identity(tmp_path):
    fd=control.exclusive(tmp_path/'lock')
    try:
        with pytest.raises(RuntimeError):control.exclusive(tmp_path/'lock')
    finally:os.close(fd)
    p=control.process_info(os.getpid())
    assert control.same_process(os.getpid(),p['start_token'])
    assert not control.same_process(os.getpid(),'wrong')


def test_ledger_reopens_without_duplicate_dispatch(tmp_path):
    db=control.connect(tmp_path)
    with db:
        control.putcurve(db,{'id':'c','j_key':'1','rank':27})
        control.putjob(db,{'id':1,'cid':'c','state':'PENDING'})
        control.putmeta(db,'intake_cursor',{'ordinal':37})
    db.close();db=control.connect(tmp_path)
    assert control.meta(db,'intake_cursor')=={'ordinal':37}
    assert control.jobs(db,('PENDING',))==[{'id':1,'cid':'c','state':'PENDING'}]
    db.close()


def test_failed_phase_is_preserved_before_retry(tmp_path,monkeypatch):
    worker.save(tmp_path/'request.json',{'config':{**policy.DEFAULTS,'sage':'sage'}})
    directory=tmp_path/'phases/generic-00-0'
    worker.save(directory/'supervisor.json',{'outcome':'parent_interrupted','returncode':None})
    (directory/'worker.log').write_text('original failure')
    def supervise(*a,**kw):
        record={'outcome':'completed','returncode':0}
        worker.save(kw['checkpoint_path'],record);return record
    monkeypatch.setattr(worker,'supervise',supervise)
    assert worker.step(tmp_path,'generic')['outcome']=='completed'
    assert len(list((tmp_path/'phases').glob('generic-00-0-attempt-*/worker.log')))==1


def test_corrupted_export_is_rejected_before_ledger_promotion(tmp_path):
    frozen=tmp_path/'runtime';job=frozen/'job';job.mkdir(parents=True)
    control.save(tmp_path/'config.json',{'frozen_root':str(frozen),'export':str(tmp_path/'exports')})
    control.save(job/'request.json',{})
    control.save(job/'point.json',{'point':'altered'})
    r={'request_sha256':control.sha(job/'request.json'),'evidence':{'point.json':'bad digest'}}
    with pytest.raises(ValueError,match='evidence changed'):
        control.export_result(tmp_path,{'path':'job','id':1,'cid':'c'},r)
    assert not (tmp_path/'exports').exists()


def test_seed_cloud_and_restarts_share_one_amplification_budget(tmp_path,monkeypatch):
    monkeypatch.setattr(worker,'ROOT',tmp_path)
    row={'id':'fresh','rank':17,'family':'f','parameter':'1','model':['0','0','0','-1','1']}
    request={'candidate':row,'config':policy.DEFAULTS,'kind':'fresh','allowance':100,
             'bank_index':0,'catalogue':'catalogue.json'}
    worker.save(tmp_path/'catalogue.json',{})
    observed=[]
    def packet(rank):
        return {'curve':row['model'],'points':[[str(i),'1'] for i in range(rank)],
                'rank_lower_bound':rank,'generic_rank':17,'proof':{},'gains':[]}
    def step(job,phase,index=0,**kw):
        if phase=='generic':
            worker.save(job/'generic-seed-gate.json',{'status':'PASS_GENERIC_SEED_GATE'})
        elif phase=='seed-cloud':
            p=packet(20);p['reconciliation_gains']=[{'seed_call':1,'after':19},{'seed_call':1,'after':20}]
            worker.save(job/'cases/fresh/seed-reconciled.json',p)
            worker.save(job/'cases/fresh/seed-search/terminal.json',{'charts':1,'point_timeouts':0,'map_attempts':[],'rank_lower_bound':18})
        elif phase=='search':
            allowance=worker.read(job/f'segment-{index:02d}.json')['allowance'];observed.append(allowance)
            before,after,calls,stop=((20,21,7,'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION') if index==0
                                     else (23,24,93,'CHART_BUDGET_EXHAUSTED'))
            t=packet(after);t.update(charts=calls,stop_reason=stop,stages=[{'epoch':0,'before':before,'after':after,'charts':calls,
                                      'censored':0,'censored_maps':0,'odd_ranks':{'3':23 if index==0 else 24}}])
            worker.save(job/f'search-{index:02d}/terminal.json',t);worker.save(job/f'search-{index:02d}/protocol.json',{})
        elif phase=='cloud':
            p=packet(23 if index==0 else 24)
            if index==0:p['gains']=[{'chart':'epoch-00/chart-0006.json','after':22},{'chart':'epoch-00/chart-0006.json','after':23}]
            worker.save(job/f'reconciled-{index:02d}.json',p)
    monkeypatch.setattr(worker,'step',step)
    monkeypatch.setattr(worker,'novelty',lambda *_:{})
    monkeypatch.setattr(worker,'conductor_gate',lambda *_:{})
    import high_rank_foundry_certificate
    monkeypatch.setattr(high_rank_foundry_certificate,'compact_packet',lambda p:p)
    r=worker.search_job(tmp_path,request)
    assert observed==[100,93] and r['calls']==101 and r['rank_lower_bound']==24
    assert [e['call'] for e in r['gain_timeline']]==[1,1,1,8,8,8,101]


def test_pending_intent_materializes_request_after_controller_crash(tmp_path,monkeypatch):
    import io
    frozen=tmp_path/'runtime';frozen.mkdir()
    config={**policy.DEFAULTS,'frozen_root':str(frozen),'sage':'sage'}
    control.save(tmp_path/'config.json',config)
    db=control.connect(tmp_path)
    j={'id':1,'cid':'c','state':'PENDING','path':'job','request':{'frozen':'request'}}
    class Process:
        pid=999999999
    monkeypatch.setattr(control.subprocess,'Popen',lambda *a,**kw:Process())
    control.spawn_job(tmp_path,db,j)
    assert control.read(frozen/'job/request.json')==j['request']
    assert control.jobs(db,('RUNNING',))[0]['id']==1
    db.close()


def test_existing_worker_lease_is_adopted_without_duplicate_process(tmp_path,monkeypatch):
    frozen=tmp_path/'runtime';job=frozen/'job';job.mkdir(parents=True)
    control.save(tmp_path/'config.json',{'frozen_root':str(frozen)})
    token=control.process_info(os.getpid())['start_token']
    control.save(job/'lease.json',{'pid':os.getpid(),'token':token,'started_at':1})
    db=control.connect(tmp_path)
    j={'id':1,'cid':'c','state':'PENDING','path':'job','request':{}}
    monkeypatch.setattr(control.subprocess,'Popen',lambda *a,**kw:pytest.fail('duplicate worker'))
    control.spawn_job(tmp_path,db,j)
    assert j['pid']==os.getpid() and j['token']==token and j['state']=='RUNNING'
    db.close()


def test_compaction_preserves_independence_and_requires_exact_replay():
    from high_rank_foundry_certificate import compact_packet
    p={'rank_lower_bound':3,'proof':{'signatures':[
        {'prime':3,'rows':[[1,0,0]]},{'prime':5,'rows':[[1,0,0]]},
        {'prime':7,'rows':[[0,1,0]]},{'prime':11,'rows':[[1,1,1]]},
        {'prime':13,'rows':[[0,0,0]]}]}}
    c=compact_packet(p)
    assert [s['prime'] for s in c['proof']['signatures']]==[3,7,11]
    assert len(p['proof']['signatures'])==5
    assert c['status'].endswith('PENDING_REPLAY')
    p['rank_lower_bound']=4
    with pytest.raises(ArithmeticError):compact_packet(p)


def test_excluded_candidates_do_not_consume_fresh_family_or_height_quota():
    pool=[]
    for f in policy.FAMILIES:
        for i in (1,2,3):
            pool.append({'family':f,'parameter':str(i),'numerator':i,'denominator':1,
                         'combined_selection_units':4-i,'combined_good':1,
                         'model':['0','0','0',str(-i),str(i)]})
    stream=intake.Intake({'pool':pool,'excluded_models':[],'reserved_addresses':[['074d9','1'],['074d9','2']]})
    row,skipped=stream.next({})
    assert (row['family'],row['parameter'],row['lane'])==('074d9','3','low')
    assert len(skipped)==2 and stream.cursor['ordinal']==1


def test_complete_point_return_survives_crash_before_chart_publication(tmp_path):
    from types import SimpleNamespace
    from high_rank_foundry_point_calls import install
    invoked=[]
    def execute(*args):
        invoked.append(1);return {'exact':'retained point transcript'},(('1','2'),)
    backend=SimpleNamespace(execute=execute,replay=lambda *args:(('1','2'),),sources=lambda:{'source':'fixed'})
    install(tmp_path,backend)
    search=SimpleNamespace(state=SimpleNamespace(key='unchanged-state'))
    first=backend.execute(search,{},125000,10,'gp')
    second=backend.execute(search,{},125000,10,'gp')
    assert first==second and len(invoked)==1
    assert len(list(tmp_path.glob('point-invocations/*/attempt-*.json')))==1
    assert len(list(tmp_path.glob('point-invocations/*/result.json')))==1
    backend.execute(search,{'different':'map'},125000,10,'gp')
    assert len(invoked)==2


def test_interrupted_point_invocation_retains_intent_and_is_not_a_miss(tmp_path):
    from types import SimpleNamespace
    from high_rank_foundry_point_calls import install
    def fail(*args):raise RuntimeError('worker interrupted')
    backend=SimpleNamespace(execute=fail,replay=lambda *args:(),sources=lambda:{})
    install(tmp_path,backend)
    with pytest.raises(RuntimeError):
        backend.execute(SimpleNamespace(state=SimpleNamespace(key='s')),{},125000,10,'gp')
    assert len(list(tmp_path.glob('point-invocations/*/attempt-*.json')))==1
    assert not list(tmp_path.glob('point-invocations/*/result.json'))


def test_due_revival_runs_even_when_all_previous_curves_are_cooled(tmp_path,monkeypatch):
    from types import SimpleNamespace
    control.save(tmp_path/'config.json',{**policy.DEFAULTS,'frozen_root':str(tmp_path),'catalogue':'catalogue.json'})
    c={**curve(rank=17,stale=200),'state':'COOLED','j_key':'1','revivals':0}
    db=control.connect(tmp_path)
    with db:
        control.putcurve(db,c)
        control.putmeta(db,'fresh_completed',40)
        control.putjob(db,{'id':1,'cid':'old','state':'DONE','kind':'fresh','charged_seconds':100,'created_at':0})
    dispatched=[]
    monkeypatch.setattr(control,'spawn_job',lambda folder,db,j:dispatched.append(j))
    def no_fresh(*args):pytest.fail('a due revival was starved by an empty ready queue')
    control.dispatch(tmp_path,db,SimpleNamespace(cursor={},next=no_fresh))
    req=dispatched[0]['request']
    assert req['revival'] and req['bank_index']>=3 and req['candidate']['revivals']==1
    assert req['kind']=='exploit' and control.meta(db,'last_revival')==40
    db.close()


def test_raw_point_journal_is_bound_into_export_evidence(tmp_path,monkeypatch):
    worker.save(tmp_path/'request.json',{'kind':'fresh'})
    journal=tmp_path/'point-invocations/hash/result.json'
    worker.save(journal,{'complete':'raw return'})
    monkeypatch.setattr(worker,'search_job',lambda *_:{'status':'engineering-control','calls':0})
    worker.run(tmp_path)
    r=worker.read(tmp_path/'result.json')
    assert r['evidence']['point-invocations/hash/result.json']==worker.sha(journal)
    assert 'lease.json' not in r['evidence'] and 'result.json' not in r['evidence']
