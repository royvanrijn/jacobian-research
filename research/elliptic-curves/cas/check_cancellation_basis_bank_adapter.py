#!/usr/bin/env sage -python
"""Meter one two-rank live-producer commissioning of the JSON-only adapter."""
import copy
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

import cancellation_basis_early as early
from cancellation_basis_epoch import read, sha, need, native_rank
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, canonical
from verify_cancellation_basis_bank_v2 import rebuild, verify_selection

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2'
RAW = LOCAL/'cancellation-basis-bank-certificate-v2'
PRE = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1'
PRE_RAW = LOCAL/'cancellation-basis-bank-certificate-v1'


def worker():
    resource.setrlimit(resource.RLIMIT_CPU,(25,30)); started = time.process_time()
    plan,cases = early.guard(); case = next(c for c in cases if c['id'] == 'fcf6ce8e6a518eb715d2')
    rows = []
    for epoch in (1,2):
        old = early.RAW/'arms'/case['id']/f'basis_refresh/epoch-{epoch:02d}'
        packet = read(old/'seed.json'); rank = native_rank(packet)
        folder = RAW/'adapter-check'/f'epoch-{epoch:02d}'
        bank,receipt = rebuild(packet,case,plan,folder)
        need(canonical(bank) == canonical(read(old/'bank.json')) and sha(folder/'bank.json') == sha(old/'bank.json'),
             'live producer serialization changes bank bytes')
        selection = read(folder/'producer/selection.json'); wrong = copy.deepcopy(selection)
        wrong['centres'].reverse()
        try: verify_selection(packet,case,plan,folder/'producer',wrong)
        except ArithmeticError as error:
            need(str(error) == 'selection object differs from bound file','unexpected mismatch rejection')
        else: raise ArithmeticError('JSON adapter accepts unbound data')
        rows.append({'epoch':epoch,'rank':rank,'bank_sha256':sha(folder/'bank.json'),
                     'verification_sha256':sha(folder/'verification.json'),'receipt':receipt,
                     'unbound_selection':'REJECTED','live_bank_python_equality':bank == read(old/'bank.json')})
    new_write(RAW/'adapter-check/result.json',{'status':'PASS_SERIALIZATION_ADAPTER','rows':rows,
        'internal_cpu_seconds':time.process_time()-started,'point_search_calls':0})


def run():
    old = read(PRE/'supervision.json'); need(old['status'] == 'STOPPED_FAILURE' and len(old['records']) == 1,'unexpected predecessor')
    previous = PRE_RAW/'arms/fcf6ce8e6a518eb715d2-e01/rational_rebuild'
    retained = early.RAW/'arms/fcf6ce8e6a518eb715d2/basis_refresh/epoch-01'
    need(read(previous/'failure.json')['error'] == 'ArithmeticError: rebuilt exported bank differs from retained bank','unexpected predecessor failure')
    need(sha(previous/'bank/bank.json') == sha(retained/'bank.json'),'predecessor bank bytes differ')
    bindings = {str(p.relative_to(ROOT)):sha(p) for p in
        [PRE/'protocol.json',PRE/'supervision.json',previous/'failure.json',previous/'bank/bank.json',retained/'bank.json']}
    new_write(OUT/'failed-predecessor.json',{'status':'RETAINED_SERIALIZATION_COMPARISON_FAILURE',
        'protocol_sha256':sha(PRE/'protocol.json'),'charged_child_cpu_seconds':old['charged_child_cpu_seconds'],
        'workers':1,'candidate_workers':0,'equivalent_exported_bank_bytes':True,'input_sha256':bindings,
        'classification':'Live Python CVP minimum tuples compare unequal to stored JSON arrays. The baseline completed independent reconstruction and wrote byte-identical bank JSON; the harness rejected only its in-memory-versus-JSON equality. The predecessor stays failed and is never resumed. No arithmetic source, bank policy, inputs or point search are changed by the successor.'})
    for name in ('development.json','corruption-checks.json','binding-checks.json'):
        target = OUT/name; need(not target.exists(),'preserve commissioning copy'); shutil.copyfile(PRE/name,target)
    recipe = {'status':'FROZEN_LIVE_PRODUCER_ADAPTER_CHECK','epochs':[1,2],
        'source_sha256':source_closure([Path(__file__),CAS/'cancellation_basis_bank_cost_v2.py']),
        'cpu_soft_seconds':25,'cpu_hard_seconds':30,'wall_seconds':90,'point_search_calls':0,
        'failed_predecessor_sha256':sha(OUT/'failed-predecessor.json')}
    new_write(RAW/'adapter-check/recipe.json',recipe)
    for name in recipe['source_sha256']:
        target = RAW/'adapter-check/sources'/name; target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,target)
    before = resource.getrusage(resource.RUSAGE_CHILDREN); wall = time.monotonic()
    with (RAW/'adapter-check/worker.log').open('w') as log:
        process = subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],stdout=log,stderr=subprocess.STDOUT,timeout=90)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = read(RAW/'adapter-check/result.json') if (RAW/'adapter-check/result.json').exists() else {'status':'RETAINED_ADAPTER_FAILURE'}
    new_write(OUT/'adapter-check.json',{**result,'returncode':process.returncode,
        'full_child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        'wall_seconds':time.monotonic()-wall,'recipe_sha256':sha(RAW/'adapter-check/recipe.json'),
        'worker_log_sha256':sha(RAW/'adapter-check/worker.log'),'adapter_sha256':sha(CAS/'verify_cancellation_basis_bank_v2.py')})
    need(process.returncode == 0 and result['status'] == 'PASS_SERIALIZATION_ADAPTER','adapter commissioning failed')
    print({k:v for k,v in read(OUT/'adapter-check.json').items() if k != 'rows'},flush=True)


if __name__ == '__main__':
    worker() if len(sys.argv)>1 and sys.argv[1] == 'worker' else run()
