#!/usr/bin/env python3
"""Successor adapter: preserve already compatible endpoint model coordinates.

The first audit stopped before target evaluation because `short` also removes
an integral scaling. This version checks the supplied short equation exactly;
it changes no anchors, representatives, limits, or finite comparison rule.
"""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import time

import cancellation_basis_accessibility as base
from cancellation_basis_accessibility import ROOT, LOCAL, PRIOR, CAS, F, read, sha, need, new_write, loaded_inputs
from finite_cancellation_corpus import write

FAILED = base.OUT
base.OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v2'
base.RAW = LOCAL/'cancellation-basis-accessibility-v2'
OUT, RAW = base.OUT, base.RAW


def same_model(model, points):
    curve = tuple(map(F, model))
    need(len(curve) == 5 and not any(curve[:3]), 'audit requires the supplied compatible short model')
    values = [tuple(map(F, p)) for p in points]
    need(all(y*y == x*x*x+curve[3]*x+curve[4] for x, y in values), 'endpoint point membership failed')
    return list(map(str, curve)), [list(map(str, p)) for p in values]


# Explicit adapter of the endpoint normalization stage; all remaining stages
# are the byte-pinned v1 implementation. The original module/file is unedited.
base.short = same_model


def freeze():
    plan = read(FAILED/'protocol.json'); failed = read(FAILED/'supervision.json')
    need(failed['status'] == 'STOPPED_RETAINED_FAILURE' and len(failed['records']) == 1, 'unexpected predecessor state')
    record = failed['records'][0]
    need(record['status'] == 'WORKER_FAILURE' and 'scan.json' not in record['files'], 'an outcome preceded the correction')
    packets = {r['case']: r['packet'] for r in read(PRIOR/'endpoint-packets.json')['rows']}
    for row in plan['cases']:
        data = loaded_inputs(row)
        need(list(map(str, map(F, packets[row['id']]['curve']))) == data['new_bank']['seed']['curve'], 'not already compatible')
    plan['source_sha256'].update({str(Path(__file__).relative_to(ROOT)): sha(Path(__file__)),
        'elliptic-curves/cas/verify_cancellation_basis_accessibility_v2.py': sha(CAS/'verify_cancellation_basis_accessibility_v2.py')})
    plan['oracle'] = plan['oracle'].replace('normalize every retained certified endpoint point onto the search model',
        'require the retained certified endpoint equation to equal the search model and preserve its exact coordinates')
    plan['predecessor'] = {'protocol_sha256': sha(FAILED/'protocol.json'), 'supervision_sha256': sha(FAILED/'supervision.json'),
        'charged_cpu_seconds': record['charged_cpu_seconds'],
        'boundary': 'One retained model-normalization failure after map preparation and before target evaluation. No point searches or outcome-driven roster change.'}
    new_write(OUT/'protocol.json', plan)
    print(json.dumps({'status': 'FROZEN_CORRECTED_MODEL_ADAPTER', 'cases': len(plan['cases']),
                      'protocol_sha256': sha(OUT/'protocol.json')}), flush=True)


def worker(case_id):
    plan = base.guard(); row = next(r for r in plan['cases'] if r['id'] == case_id)
    resource.setrlimit(resource.RLIMIT_CPU, (plan['hard_cpu_seconds_per_case'], plan['hard_cpu_seconds_per_case']))
    folder = RAW/case_id; new_write(folder/'start.json', {'protocol_sha256': sha(OUT/'protocol.json')})
    maps = base.prepare_maps(row, plan, folder)
    base.evaluate(row, plan, folder, maps)
    from verify_cancellation_basis_accessibility_v2 import verify
    new_write(folder/'verification.json', verify(row, plan, folder))
    print(json.dumps({'case': case_id, 'status': 'PASS', 'cpu_seconds': base.cpu()}), flush=True)


def run():
    plan = base.guard(); records = []
    for row in plan['cases']:
        folder = RAW/row['id']; need(not folder.exists(), 'no audit retry or refill'); folder.mkdir(parents=True)
        t = time.monotonic(); before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with (folder/'worker.log').open('w') as log:
            try:
                child = subprocess.run(['sage', '-python', str(Path(__file__)), 'worker', '--case', row['id']],
                    stdout=log, stderr=subprocess.STDOUT, timeout=plan['wall_seconds_per_case'])
                code = child.returncode; status = 'COMPLETE' if code == 0 else 'WORKER_FAILURE'
            except subprocess.TimeoutExpired:
                code = None; status = 'WALL_LIMIT'
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        item = {'case': row['id'], 'status': status, 'returncode': code,
            'charged_cpu_seconds': after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'wall_seconds': time.monotonic()-t,
            'files': {p.name: sha(p) for p in folder.iterdir() if p.is_file()}}
        new_write(folder/'receipt.json', item); records.append(item)
        write(OUT/'supervision.json', {'status': 'RUNNING', 'records': records}); print(json.dumps(item), flush=True)
        if status != 'COMPLETE':
            write(OUT/'supervision.json', {'status': 'STOPPED_RETAINED_FAILURE', 'records': records}); return
    write(OUT/'supervision.json', {'status': 'COMPLETE', 'records': records})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['freeze', 'worker', 'run', 'report'])
    parser.add_argument('--case'); args = parser.parse_args()
    if args.command == 'worker': worker(args.case)
    elif args.command == 'report': base.report()
    else: globals()[args.command]()
