"""The real -1 specializations are input misses, not controller failures."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import high_rank_foundry_seed_gate as gate
import high_rank_foundry_job as worker
import run_high_rank_foundry as control
from high_rank_foundry_policy import DEFAULTS


def candidate(family='074d9',parameter='-1'):
    family_row=next(f for f in gate.read(gate.atlas.ATLAS)['families'] if f['family']==family)
    model,_=gate.atlas.specialize(family_row,parameter)
    return control.fresh_record({'id':'control-'+family,'family':family,'parameter':parameter,
                                 'model':list(map(str,model)),'j_key':'control','lane':'raw'})


@pytest.mark.parametrize('family',['074d9','07ca9','08234','08f72','103b2','11952'])
def test_actual_small_specializations_stay_unknown(family):
    result=gate.assess(candidate(family))
    assert result['status']=='UNRESOLVED_GENERIC_SEED'
    assert result['admitted_columns']<17 and result['point_searches']==0
    assert 'rank_lower_bound' not in result and 'rank_upper_bound' not in result


def test_mismatched_model_is_still_an_engineering_failure():
    row=candidate();row['model'][4]=str(int(row['model'][4])+1)
    with pytest.raises(ValueError,match='selected equation differs'):gate.assess(row)


def test_independent_known_seed_still_passes():
    r=gate.assess(candidate('07ca9','-1508/909'))
    assert r['status']=='PASS_GENERIC_SEED_GATE' and r['admitted_columns']==17
    assert r['no_rational_2_torsion_prime'] is not None


def test_seed_miss_stops_before_point_search_and_does_not_break_controller(tmp_path,monkeypatch):
    monkeypatch.setattr(worker,'ROOT',tmp_path)
    config={**DEFAULTS,'frozen_root':str(tmp_path),'export':str(tmp_path/'exports'),'sage':'sage'}
    control.save(tmp_path/'config.json',config)
    control.save(tmp_path/'manifest.json',{'scope':'isolated test'})
    row=candidate();gate_result=gate.assess(row)
    observed=[]
    def step(job,phase,**kwargs):
        observed.append(phase)
        assert phase=='generic'
        worker.save(job/'generic-seed-gate.json',gate_result)
    monkeypatch.setattr(worker,'step',step)
    monkeypatch.setattr(control,'guard',lambda _:tmp_path)
    db=control.connect(tmp_path)
    with db:
        control.putmeta(db,'status','RUNNING');control.putmeta(db,'consecutive_failures',2)
        control.putcurve(db,row)
    # Six sequential misses must not reach the engineering halt threshold;
    # they also must not erase the existing engineering failure streak.
    for n in range(1,7):
        job=tmp_path/f'job-{n}'
        req={'kind':'fresh','candidate':row,'config':config,'allowance':100}
        worker.save(job/'request.json',req);worker.run(job)
        j={'id':n,'cid':row['id'],'state':'RUNNING','path':job.name,'created_at':0,
           'initial_rank':17,'reserved_calls':200,'kind':'fresh'}
        control.consume(tmp_path,db,j)
        assert j['state']=='SKIPPED' and j['charged_calls']==0
    assert observed==['generic']*6 and control.meta(db,'status')=='RUNNING'
    assert control.meta(db,'consecutive_failures')==2
    c=control.curves(db)[0]
    assert c['rank'] is None and c['state']=='UNCERTIFIED_SEED' and not c['packet']
    assert not list(tmp_path.rglob('point-invocations'))
    db.close()


def test_four_backend_failures_still_halt(tmp_path,monkeypatch):
    config={**DEFAULTS,'frozen_root':str(tmp_path),'export':str(tmp_path/'exports')}
    control.save(tmp_path/'config.json',config)
    monkeypatch.setattr(control,'guard',lambda _:tmp_path)
    db=control.connect(tmp_path)
    row=candidate()
    with db:control.putcurve(db,row)
    for n in range(1,5):
        job=tmp_path/f'job-{n}';control.save(job/'request.json',{})
        control.save(job/'result.json',{'status':'UNRESOLVED_JOB_FAILURE','calls':None,'wall_seconds':1,
            'request_sha256':control.sha(job/'request.json'),'evidence':{}})
        j={'id':n,'cid':row['id'],'path':job.name,'created_at':0,'reserved_calls':200}
        control.consume(tmp_path,db,j)
    assert control.meta(db,'status')=='HALTED_REPEATED_FAILURE'
    db.close()


def test_migration_preserves_cursor_budget_and_cooling(tmp_path,monkeypatch):
    import high_rank_foundry_migration as migration
    import high_rank_foundry_arithmetic as arithmetic
    old=tmp_path/'parent';oldroot=old/'runtime';oldroot.mkdir(parents=True)
    new=tmp_path/'new';(new/'foundry-inputs').mkdir(parents=True)
    monkeypatch.setattr(control,'guard',lambda _:oldroot)
    control.save(old/'config.json',{'frozen_root':str(oldroot),'intake':'intake.json','catalogue':'catalogue.json'})
    control.save(old/'manifest.json',{'original':'manifest'})
    for name in ('intake','catalogue'):control.save(oldroot/(name+'.json'),{'original':name})
    row=candidate();row.update(rank=26,state='COOLED',head='search-old',stale_calls=350,total_calls=501)
    control.save(oldroot/'packet.json',{'exact':'original bytes'})
    row.update(packet='packet.json',packet_sha256=control.sha(oldroot/'packet.json'))
    db=control.connect(old)
    with db:
        control.putcurve(db,row);control.putmeta(db,'status','HALTED_REPEATED_FAILURE')
        control.putmeta(db,'intake_cursor',{'ordinal':34,'shells':{'074d9':[1,1]}})
        control.putmeta(db,'fresh_completed',30)
        control.putjob(db,{'id':71,'cid':row['id'],'state':'DONE','path':'job-old','kind':'fresh',
                           'charged_calls':101,'charged_seconds':123,'created_at':1})
    db.close();migration.freeze(old,new)
    monkeypatch.setattr(arithmetic,'verify_packet',lambda p,c:None)
    target=tmp_path/'target';target.mkdir();db=control.connect(target)
    assert migration.restore(new,db)==1
    after=control.curves(db)[0]
    assert after['rank']==26 and after['state']=='COOLED' and after['stale_calls']==350
    assert after['total_calls']==501 and after['head'] is None and after['bank_index']==1
    assert control.sha(new/after['packet'])==row['packet_sha256']
    assert control.meta(db,'intake_cursor')['ordinal']==34
    j=control.jobs(db)[0]
    assert j['id']==71 and j['charged_calls']==101 and j['charged_seconds']==123
    assert j['path']==str(oldroot/'job-old')
    db.close()
