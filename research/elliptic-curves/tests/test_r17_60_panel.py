import importlib.util
import json
from pathlib import Path
import sys

import pytest

CAS = Path(__file__).resolve().parents[1] / 'cas'
sys.path.insert(0, str(CAS))
import select_r17_60_panel as select
import run_r17_60_panel as run


def population():
    rows = []
    for f, family in enumerate(select.FAMILIES):
        for high in (False, True):
            for i in (1, 2, 3):
                k = 20*f+5*high+i
                rows.append({'family': family, 'parameter': str(2000+i if high else i),
                             'numerator': 2000+i if high else i, 'denominator': 1,
                             'model': ['0','0','0',str(-k),str(k)],
                             'combined_selection_units': 100-i, 'combined_good': 10})
    return {'pool': rows, 'excluded_models': [], 'reserved_addresses': []}


def test_stratified_roster_is_outcome_independent():
    d = population()
    before = select.select(d, 1)['rows']
    for row in d['pool']:
        row['rank'] = 99 if row['numerator'] % 2 == 0 else 17
        row['validation_score'] = -row['combined_selection_units']
    after = select.select(d, 1)['rows']
    assert [(r['family'],r['parameter']) for r in before] == [(r['family'],r['parameter']) for r in after]
    assert len(before) == 12
    assert sum(r['stratum']=='low' for r in before) == 6


def test_reservations_and_isomorphic_duplicates_are_excluded():
    d = population()
    first = d['pool'][0]
    d['reserved_addresses'] = [[first['family'],first['parameter']]]
    d['pool'][2]['model'] = d['pool'][1]['model']
    rows = select.select(d, 1)['rows']
    assert next(r for r in rows if r['family']==first['family'] and r['stratum']=='low')['parameter']=='2'
    d['excluded_models'] = [d['pool'][1]['model']]
    with pytest.raises(ValueError, match='not enough fresh'):
        select.select(d, 1)


def test_gain_metrics_include_reconciliation_without_extra_calls():
    t = {'charts':7, 'stages':[{'epoch':0,'before':20,'after':21,'charts':3},
                              {'epoch':1,'before':21,'after':22,'charts':4}]}
    rc = {'gains':[{'chart':'epoch-01/chart-0003.json','after':23}]}
    timeline, calls = run.call_metrics(2,18,[(t,rc)], [{'seed_call':2,'after':19},{'seed_call':2,'after':20}])
    assert calls==7
    assert [g['call'] for g in timeline if g['phase'].startswith('complement')]==[3,7,7]
    assert timeline[-1]['after']==23


def test_cumulative_budget_survives_reconciliation_restart(tmp_path, monkeypatch):
    monkeypatch.setattr(run, 'D', tmp_path)
    monkeypatch.setattr(run, 'ROOT', tmp_path)
    run.atomic(tmp_path/'roster.json', {'rows':[{'id':'c','family':'f','parameter':'1','stratum':'low'}]})
    run.atomic(tmp_path/'protocol.json', {})
    observed = []
    folder = tmp_path/'cases/c'
    def step(case, phase, index=0, remaining=100):
        observed.append((phase,index,remaining))
        if phase=='seed-prepare':
            run.atomic(folder/'seed-search/protocol.json', {})
        elif phase=='seed-search':
            run.atomic(folder/'seed-search/terminal.json', {'charts':2,'rank_lower_bound':18,'status':'FIRST_M18_CERTIFIED','point_timeouts':0})
        elif phase=='seed-replay':
            run.atomic(folder/'seed-search/verified.json', {})
        elif phase=='seed-cloud':
            run.atomic(folder/'seed-reconciled.json', {'rank_lower_bound':18,'reconciliation_gains':[]})
        elif phase=='complement-prepare':
            run.atomic(folder/'complement-preparation/prepared.json', {})
        elif phase=='search':
            charts = 5 if index==0 else 95
            stage={'epoch':0,'before':18 if index==0 else 20,'after':19 if index==0 else 21,
                   'charts':charts,'censored':0,'censored_maps':0}
            run.atomic(folder/f'complement-{index:02d}/terminal.json', {'charts':charts,
                       'rank_lower_bound':stage['after'],'stages':[stage],
                       'stop_reason':'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION' if index==0 else 'CHART_BUDGET_EXHAUSTED'})
        elif phase=='replay':
            run.atomic(folder/f'complement-{index:02d}/verified.json', {})
        elif phase=='reconcile':
            run.atomic(folder/f'reconciled-{index:02d}.json', {'rank_lower_bound':20 if index==0 else 21,
                       'gains':[{'chart':'epoch-00/chart-0004.json','after':20}] if index==0 else []})
        elif phase=='rebuild':
            run.atomic(folder/f'preparation-{index+1:02d}/prepared.json', {})
    monkeypatch.setattr(run,'step',step)
    result=run.case_run('c')
    assert [(i,remaining) for phase,i,remaining in observed if phase=='search']==[(0,100),(1,95)]
    assert result['complement_calls']==100 and result['total_calls']==102
    assert result['rank_lower_bound']==21
    assert result['last_complement_gain_call']==100
    assert sum(phase=='reconcile' for phase,_,_ in observed)==4
    before=len(observed)
    assert run.case_run('c')==result
    assert len(observed)==before


def test_exclusive_case_lock(tmp_path):
    fd=run.lock(tmp_path/'lock')
    try:
        with pytest.raises(RuntimeError):run.lock(tmp_path/'lock')
    finally:run.os.close(fd)


def test_ranking_keeps_resource_stop_unresolved(tmp_path, monkeypatch):
    monkeypatch.setattr(run,'D',tmp_path)
    run.atomic(tmp_path/'roster.json',{'rows':[{'id':'c'}]})
    run.atomic(tmp_path/'protocol.json',{})
    run.atomic(tmp_path/'cases/c/result.json',{'id':'c','status':'RESOURCE_STOP_UNRESOLVED','reason':'cap'})
    r=run.report()
    assert r['completed']==1 and r['independently_verified']==0
    assert r['status']=='PARTIAL_R17_60_PANEL'
