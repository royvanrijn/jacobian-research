#!/usr/bin/env python3
"""Independent read-only accounting audit; no report-generator calculations."""
import argparse
import json
import math
from collections import Counter
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'nearcut60-mw16-pari-v2'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
REPORT=ART/'nearcut60v2_mw16_experiment_v1.json'
OUT=ART/'nearcut60v2_mw16_accounting_replay_v1.json'


def same(actual,expected):
    if isinstance(expected,float):
        if not isinstance(actual,(int,float)) or not math.isfinite(actual) or not math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-12):
            raise ArithmeticError('accounting numeric mismatch')
    elif actual!=expected:raise ArithmeticError('accounting value mismatch')


def totals(rows):
    known=[r['gain'] for r in rows if r['gain'] is not None]
    boxes=[r['boxes'] for r in rows if r['boxes'] is not None]
    spent=sum(r['discovery'] for r in rows)
    inclusive=spent+sum(r['verification'] for r in rows)
    gains_complete=len(known)==len(rows); exposure_complete=len(boxes)==len(rows)
    return {'allocated_curves':len(rows),'certified_gain_sum_known':sum(known),
        'unresolved_gain_curves':len(rows)-len(known),
        'certified_gain_distribution':dict(Counter(map(str,known))),
        'completed_boxes_known':sum(boxes),'unresolved_exposure_curves':len(rows)-len(boxes),
        'allocated_boxes':43*len(rows),'discovery_worker_seconds':spent,
        'including_verification_seconds':inclusive,
        'gain_per_discovery_second':sum(known)/spent if gains_complete and spent else None,
        'gain_per_inclusive_second':sum(known)/inclusive if gains_complete and inclusive else None,
        'gain_per_completed_box':sum(known)/sum(boxes) if gains_complete and exposure_complete and sum(boxes) else None,
        'completion_fraction':sum(boxes)/(43*len(rows)) if exposure_complete and rows else None,
        'censored_or_failed_workers':sum(r['worker_status']!='COMPLETE_DECLARED_POINT_ATTEMPT' for r in rows)}


def audit():
    bindings={}
    def read(path):
        raw=path.read_bytes();bindings[str(path.relative_to(ROOT))]=sha256(raw).hexdigest();return json.loads(raw)
    def bind(path,h):
        raw=path.read_bytes();actual=sha256(raw).hexdigest()
        if actual!=h:raise ArithmeticError('changed accounting evidence: '+str(path))
        bindings[str(path.relative_to(ROOT))]=actual
    def elapsed(s):
        if s['outcome']=='running' or s['wall_seconds']<0:raise ArithmeticError('nonterminal or negative cost')
        bind(Path(s['log']),s['log_sha256']);return s['wall_seconds']
    report=read(REPORT)
    for name,h in report['sources'].items():bind(ROOT/name,h)
    p=read(D/'protocol.json')
    selected=read(ART/'retained_mw16_nearcutoff_selection_v1.json')['selected']
    same(p['rows'],selected)
    same(len(selected),60)
    same(Counter(r['arm'] for r in selected),Counter(nearcutoff=60))
    same(p['stop_rank'],None);same(p['height'],125000);same(p['seconds_per_chart'],10)
    same(p['worker_wall_seconds'],600);same(p['maximum_workers'],2)
    ledgers={}
    for name in ('map','baseline','verification'):
        path=D/(name+'-ledger.json')
        ledgers[name]=read(path) if path.exists() else {'rows':[]}
    path=D/'ledger.json';ledgers['point']=read(path) if path.exists() else {'rows':[]}
    for d in ledgers.values():
        if d.get('status')=='RUNNING':raise ArithmeticError('active computation cannot be an outcome')
    indices={name:{r['id']:r for r in d['rows']} for name,d in ledgers.items()}
    same([r['id'] for r in report['rows']],[r['id'] for r in selected])
    independent=[]
    for choice,reported in zip(selected,report['rows']):
        ident=choice['id'];folder=D/ident
        for key in ('id','arm','block','family','band','sign','parameter','retained_rank','late_rank','j_height','parameter_height'):
            same(reported[key],choice[key])
        m=indices['map'][ident];map_time=elapsed(m['supervision']);base=indices['baseline'].get(ident)
        base_time=elapsed(base['supervision']) if base else 0
        point=indices['point'].get(ident);point_time=elapsed(point['point_supervision']) if point else 0
        worker_status=point['status'] if point else 'NOT_RUN_PRESEARCH_GATE'
        verify=indices['verification'].get(ident);verify_time=base_time;gain=lower=None
        attempts=returned=completed=timeouts=failures=inflight=unattempted=cpu=missing_cpu=None
        exposure_valid=False
        if verify:
            same(read(folder/'verification.json'),verify)
            for name,h in verify['inputs'].items():bind(folder/name,h)
            verify_time+=sum(elapsed(s['supervision']) for s in verify['steps'])
            exposure_valid=any(s['name']=='geometry' and s['status']=='PASS' for s in verify['steps'])
            if exposure_valid:
                exposure=read(folder/'exposure.json')['charts'];attempts=len(exposure)
                returned=completed=timeouts=failures=cpu=missing_cpu=0
                for i,c in enumerate(exposure):
                    same(c['index'],i)
                    if c['status']=='STARTED':missing_cpu+=1;continue
                    same(c['status'],'RETURNED');returned+=1
                    search=c['search'];state=search['status']
                    if state=='bounded_search_complete':
                        same(search['returncode'],0)
                        if '***' in search['stderr']:raise ArithmeticError('failed PARI call counted complete')
                        completed+=1
                    elif state=='bounded_search_timeout':timeouts+=1
                    elif state=='backend_failure':failures+=1
                    else:raise ArithmeticError('unknown chart outcome')
                    if search['search_cpu_ms'] is None:missing_cpu+=1
                    else:cpu+=search['search_cpu_ms']
                inflight=attempts-returned;unattempted=43-attempts
                if not 0<=inflight<=1 or unattempted<0:raise ArithmeticError('invalid exposure prefix')
            if verify['status']=='PASS':
                same([s['name'] for s in verify['steps']],['history','geometry','cloud-build','cloud-check','odd-build','odd-check','provenance'])
                for step in verify['steps']:
                    same(step['status'],'PASS');same(step['supervision']['outcome'],'completed');same(step['supervision']['returncode'],0)
                cloud=read(ROOT/verify['cloud_path']);odd=read(ROOT/verify['odd_path'])
                lower=max([cloud['rank_lower_bound']]+[a['finite_column_rank'] for a in odd['audits']]);gain=lower-16
                if gain<0:raise ArithmeticError('lost certified baseline')
                same(verify['certified_gain'],gain)
        fields={'map_seconds':map_time,'point_worker_seconds':point_time,'discovery_worker_seconds':map_time+point_time,
            'verification_seconds':verify_time,'worker_status':worker_status,'certified_gain':gain,'rank_lower_bound':lower,
            'allocated_boxes':43,'attempted_boxes':attempts,'returned_charts':returned,'completed_boxes':completed,
            'chart_timeouts':timeouts,'backend_failures':failures,'in_flight_at_termination':inflight,
            'unattempted_boxes':unattempted,'known_search_cpu_ms':cpu,'charts_without_cpu_time':missing_cpu,
            'exposure_verified':exposure_valid,'verification_status':verify['status'] if verify else 'NOT_RUN'}
        for key,value in fields.items():same(reported[key],value)
        independent.append({'id':ident,'arm':choice['arm'],'family':choice['family'],'block':choice['block'],
            'gain':gain,'boxes':completed,'discovery':map_time+point_time,'verification':verify_time,'worker_status':worker_status})
    def compare_summary(actual,rows):
        computed=totals(rows)
        for key,value in computed.items():same(actual[key],value)
        return computed
    for arm in ('nearcutoff',):
        compare_summary(report['arms'][arm],[r for r in independent if r['arm']==arm])
    for family in sorted({r['family'] for r in independent}):
        for arm in ('nearcutoff',):
            compare_summary(report['families'][family][arm],[r for r in independent if r['family']==family and r['arm']==arm])
    aggregate=compare_summary(report['totals'],independent)
    bindings[str(Path(__file__).resolve().relative_to(ROOT))]=sha256(Path(__file__).read_bytes()).hexdigest()
    return {'schema':'elliptic-curves.nearcut60v2-independent-accounting.v1','status':'PASS','sources':bindings,
        'allocated_curves':60,'allocated_boxes':2580,'independent_rows':independent,'totals':aggregate,
        'scope':'Independent read-only reconstruction of all60 allocated curve costs, certified gain bindings, terminal chart completion, missing CPU timing, arm/family rates and whole-cohort totals. Exact rank and geometry proofs are separately checked by their recorded workers; this audit checks accounting, not point search completeness or a rank-density theorem. No report-generator functions imported and no new search.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();data=audit()
    if args.check:
        if json.loads(OUT.read_text())!=data:raise ArithmeticError('independent accounting replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve independent accounting audit')
        OUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('INDEPENDENT STRATA ACCOUNTING',data['status'],data['totals']['certified_gain_sum_known'],flush=True)
