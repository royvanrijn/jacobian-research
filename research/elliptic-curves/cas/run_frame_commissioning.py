#!/usr/bin/env python3
"""Bounded arm-prefix commissioning over an unchanged frozen search runtime."""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import fcntl
import gzip
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time

from run_class1_prospective_search import read, write, sha, selection, allowance


def panel_rows(rows, ranked, controls):
    used = {False: 0, True: 0}
    limits = {False: ranked, True: controls}
    result = []
    for row in selection(rows):
        arm = row['control']
        if used[arm] < limits[arm]:
            result.append(row)
            used[arm] += 1
    if used != limits:
        raise ValueError('Insufficient rows for fixed arm quotas')
    return result


def decision(terminals):
    # A partially failed campaign cannot certify a complete null panel.
    gains = [r for r in terminals if (r.get('rank_lower_bound') or 0) >= 20]
    if gains:
        return 'CERTIFIED_GE20_ESCALATED_REVIEW_BEFORE_NEXT_FRAME'
    if any(r['status'] != 'COMPLETE_BOUNDED' for r in terminals):
        return 'UNKNOWN_INCOMPLETE_COMMISSIONING'
    return 'NO_EVIDENCE_CURRENT_SEARCH_PRODUCTIVE'


def prepare(folder, ranked, controls):
    plan = read(folder/'plan.json')
    window = Path(plan['root'])/'ordinary-search/window-000'
    rows = panel_rows(read(window/'scores.json')['rows'], ranked, controls)
    allowed = {r['index'] for r in rows}
    dispatched = {int(p.parent.name) for p in (window/'cases').glob('*/selection.json')}
    assert dispatched <= allowed, 'Existing dispatch exceeds requested arm prefix'
    amendment = {
        'schema': 'frame.commissioning-cap.v1',
        'ranked_cap': ranked, 'control_cap': controls,
        'indices': [r['index'] for r in rows],
        'plan_sha256': sha(folder/'plan.json'),
        'scores_sha256': sha(window/'scores.json'),
        'queue_sha256': sha(window/'queue.json'),
        'controller_sha256': sha(Path(__file__)),
        'scoring_change': False,
        'selection': 'Existing frozen ranked prefix and score-independent control prefix; retain within-arm order.',
        'ge20': 'Finish existing adaptive escalation before any next-frame realization.',
        'no_ge20': 'NO_EVIDENCE_CURRENT_SEARCH_PRODUCTIVE; bounded protocol outcome, never a rank upper bound.',
        'incomplete': 'UNKNOWN_INCOMPLETE_COMMISSIONING; no automatic frame release.'}
    write(folder/'commissioning-cap.json', amendment, True)


def freeze_baseline(folder, window, rows, terminals, outcome):
    out = folder/'matched-baseline'
    out.mkdir(exist_ok=True)
    # Keep every score in the full fixed window, not just successful/selected rows.
    for source, name in [(window/'scores.json', 'scores.json.gz'),
                         (window/'queue.json', 'queue.json.gz')]:
        dest = out/name
        data = gzip.compress(source.read_bytes(), mtime=0)
        if dest.exists():
            assert dest.read_bytes() == data
        else:
            dest.write_bytes(data)
    packets = []
    for row in rows:
        case = window/'cases'/str(row['index'])
        for p in sorted(case.rglob('*.json')):
            packets.append({'path': str(p.relative_to(window)), 'sha256': sha(p)})
    write(out/'results.json', {'outcome': outcome, 'parameters': terminals,
          'selected_rows': rows, 'artifacts': packets,
          'interpretation': 'Operational evidence under a bounded fixed protocol; no rank upper bound or parent exclusion theorem.'}, True)
    write(out/'manifest.json', {p.name: sha(p) for p in sorted(out.iterdir()) if p.name != 'manifest.json'}, True)
    return out


def run(folder):
    amendment = read(folder/'commissioning-cap.json')
    assert sha(Path(__file__)) == amendment['controller_sha256']
    plan = read(folder/'plan.json')
    rt = Path(plan['root'])
    window = rt/'ordinary-search/window-000'
    for p, key in [(folder/'plan.json','plan_sha256'), (window/'scores.json','scores_sha256'),
                   (window/'queue.json','queue_sha256')]:
        assert sha(p) == amendment[key]
    # The old controller drains its current batches, then releases this same lock.
    lock = (folder/'controller.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX)
    if (folder/'commissioning-result.json').exists():
        return
    previous = read(folder/'STATUS.json') if (folder/'STATUS.json').exists() else {}
    if previous.get('status') not in ('STOPPED_AFTER_DRAIN', 'COMMISSIONING', 'DRAINING', None):
        raise RuntimeError('Unexpected predecessor state: '+str(previous))
    frozen = rt/'elliptic-curves/cas/run_class1_prospective_search.py'
    sys.path.insert(0, str(frozen.parent))
    spec = importlib.util.spec_from_file_location('frozen_search', frozen)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.guard(folder)
    write(folder/'commissioning-predecessor.json', previous, True)
    (folder/'STOP').unlink(missing_ok=True)
    write(folder/'controller.json', {'pid': os.getpid(), 'started_at': time.time(),
          'mode': 'FIXED_ARM_COMMISSIONING', 'strict_class_required': False})
    rows = panel_rows(read(window/'scores.json')['rows'], amendment['ranked_cap'], amendment['control_cap'])
    assert [r['index'] for r in rows] == amendment['indices']
    pending = [r for r in rows if not (window/'cases'/str(r['index'])/'terminal.json').exists()]
    # Resume any certified >=20 trajectory before dispatching fresh fibres.
    def priority(row):
        path = window/'cases'/str(row['index'])/'progress.json'
        return -(read(path).get('rank_lower_bound',17) if path.exists() else 17)
    pending.sort(key=priority)
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        active = {}
        while pending or active:
            while pending and len(active)<plan['workers'] and not (folder/'STOP').exists():
                row = pending.pop(0)
                active[pool.submit(module.fibre,folder,window,row)] = row['index']
            if not active:
                break
            done,_ = wait(active,timeout=10,return_when=FIRST_COMPLETED)
            for future in done:
                index = active.pop(future)
                result = future.result()
                print('FIBRE',index,result['status'],'rank',result['rank_lower_bound'],flush=True)
            write(folder/'STATUS.json', {'status':'DRAINING' if (folder/'STOP').exists() else 'COMMISSIONING',
                'active_indices':list(active.values()),'remaining':len(pending),
                'completed':sum((window/'cases'/str(r['index'])/'terminal.json').exists() for r in rows),
                'ranked_cap':amendment['ranked_cap'],'control_cap':amendment['control_cap'],
                'strict_class_required':False,'controller_cpu_seconds':time.process_time(),'updated_at':time.time()})
    if (folder/'STOP').exists():
        write(folder/'STATUS.json', {'status':'STOPPED_AFTER_DRAIN','updated_at':time.time()})
        return
    terminals = [read(window/'cases'/str(r['index'])/'terminal.json') for r in rows]
    outcome = decision(terminals)
    # module.fibre does not return COMPLETE_BOUNDED until its adaptive cap is met.
    for r in terminals:
        if r['status']=='COMPLETE_BOUNDED' and (r.get('rank_lower_bound') or 0)>=20:
            assert r['calls']>=allowance(r['rank_lower_bound']) or r['rank_lower_bound']>=32
    baseline = freeze_baseline(folder,window,rows,terminals,outcome)
    result = {'status':outcome,'ranked':amendment['ranked_cap'],'controls':amendment['control_cap'],
        'best_certified_lower_bound':max((r.get('rank_lower_bound') or 0) for r in terminals),
        'calls':sum(r['calls'] for r in terminals),
        'search_cpu_seconds':sum(r['search_cpu_seconds'] for r in terminals),
        'baseline':str(baseline),'baseline_manifest_sha256':sha(baseline/'manifest.json'),
        'next_frame_released':outcome=='NO_EVIDENCE_CURRENT_SEARCH_PRODUCTIVE',
        'boundary':'Search productivity classification only; no mathematical rank upper bound.',
        'updated_at':time.time()}
    write(folder/'commissioning-result.json',result,True)
    write(folder/'STATUS.json',result)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['prepare','run'])
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--ranked',type=int,default=64)
    p.add_argument('--controls',type=int,default=16)
    a=p.parse_args();a.folder=a.folder.resolve()
    if a.mode=='prepare':prepare(a.folder,a.ranked,a.controls)
    else:run(a.folder)
