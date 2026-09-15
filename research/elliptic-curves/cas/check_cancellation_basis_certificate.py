#!/usr/bin/env sage -python
"""Frozen transcript-only integration of the certified-bank epoch engine."""
import argparse
import copy
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time
import traceback

import cancellation_basis_early as early
from cancellation_basis_bank_cost_v2 import guard as cost_guard, OUT as COST
from cancellation_basis_epoch import read, sha, need
from cancellation_basis_certificate_receipts import CHECKERS
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from cancellation_scheduler_cpu import atomic
from finite_cancellation_corpus import ROOT, LOCAL, canonical, digest, write

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_certificate_integration_v1'
RAW = LOCAL/'cancellation-basis-certificate-integration-v1'


def seal():
    plan,cases = early.guard(); cost_guard()
    need(read(COST/'completion.json')['engineering_gate'] is True,'bank cost prerequisite failed')
    cases = {c['id']:c for c in cases}; rows = []
    descriptions = [('single_refresh','fcf6ce8e6a518eb715d2','basis_refresh',False),
                    ('single_fixed','fcf6ce8e6a518eb715d2','fixed_bank',False),
                    ('double_refresh','208dde8d631796bb6548','basis_refresh',False),
                    ('bank_alarm','fcf6ce8e6a518eb715d2','basis_refresh',True)]
    bindings = {}
    for label,ident,arm,alarm in descriptions:
        origin = early.RAW/'arms'/ident/arm; result = read(origin/'result.json')
        need(result['status'] == 'CERTIFIED_TARGET','integration source does not reach target')
        calls = result['calls']
        if alarm:
            stop = next(i for i,c in enumerate(calls) if c['gain'])+1; calls = calls[:stop]
        inputs = []
        for item in calls:
            p = origin/item['file']; need(sha(p) == item['sha256'],'retained transcript changed')
            bindings[str(p.relative_to(ROOT))] = sha(p); inputs.append({'file':str(p.relative_to(ROOT)),'sha256':sha(p)})
        bindings[str((origin/'result.json').relative_to(ROOT))] = sha(origin/'result.json')
        rows.append({'id':label,'case':cases[ident],'arm':arm,'alarm':alarm,'calls':inputs,
            'origin':str(origin.relative_to(ROOT)),
            'expected_directions':result['first_cloud_gain'] if alarm else result['new_directions'],
            'expected_later':0 if alarm else result['later_cloud_directions'],
            'expected_epochs':[18] if alarm else [e['rank'] for e in result['epochs']]})
    configuration = {key:plan[key] for key in ('arms','bank_policy','maximum_extension_dimension','maximum_centres',
        'prime_bound','bank_wall_seconds','target_directions','maximum_calls','height','point_wall_seconds','gp_sha256')}
    plan = {**configuration,'status':'FROZEN_TRANSCRIPT_ONLY_CERTIFICATE_INTEGRATION',
        'hard_process_cpu_seconds':60,'working_cpu_seconds':40,
        'bank_certificate_sources':{'elliptic-curves/cas/'+n:sha(CAS/n) for n in CHECKERS}}
    sources = source_closure([Path(__file__),CAS/'cancellation_basis_certificate_epoch.py',
        CAS/'verify_cancellation_basis_certificate.py']+list((ROOT/'elliptic-curves/ecsearch').glob('*.py')))
    protocol = {'status':'FROZEN_FOUR_INTEGRATION_WORKERS','jobs':rows,'plan':plan,'source_sha256':sources,
        'retained_input_sha256':bindings,'bank_cost_protocol_sha256':sha(COST/'protocol.json'),
        'maximum_workers':4,'worker_hard_cpu_seconds':60,'worker_wall_seconds':120,'campaign_child_cpu_seconds':240,
        'point_search_calls':0,
        'injection':'Only the point executor is replaced by literal retained transcripts, verified against each requested map, box and subgroup before returning. The alarm case raises Sage AlarmInterrupt after the first enlarged producer has written its selection, before its independent certificate. All other bank construction, map preparation, cloud reconciliation, native rank and final replay paths execute normally.',
        'gate':'All three positive jobs reproduce their full source call counts, gain counts and search-bank rank sequences, including M18 to M20 and an uncapped four-direction target overshoot. The real Sage alarm preserves the first complete certified cloud, leaves its unfinished bank unused, and finishes independent replay. Updated outer hashes must not conceal an obsolete seed, an omitted cloud point or a changed minimum packet.',
        'negative_checks':['obsolete_seed','omitted_cloud_point','changed_minimum_packet'],
        'boundary':'Software and arithmetic integration on known retained transcripts; no new point search or search-performance observation. Complete minimum exhaustion occurs before bank use, and final replay rechecks its bound direct witnesses. All prior protocols, failure records and timing results remain unchanged.'}
    new_write(OUT/'protocol.json',protocol)
    for name in sources:
        p = RAW/'sources'/name; p.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,p)
    print(json.dumps({'status':protocol['status'],'sources':len(sources),'jobs':len(rows),'protocol_sha256':sha(OUT/'protocol.json')}),flush=True)


def guard():
    protocol = read(OUT/'protocol.json')
    for name,h in protocol['source_sha256'].items(): need(sha(ROOT/name) == h,'integration source changed: '+name)
    for name,h in protocol['retained_input_sha256'].items(): need(sha(ROOT/name) == h,'retained transcript changed: '+name)
    need(sha(COST/'protocol.json') == protocol['bank_cost_protocol_sha256'],'bank cost protocol changed')
    need(sha(Path('/usr/bin/gp')) == protocol['plan']['gp_sha256'],'GP executable changed')
    return protocol


def worker(label):
    import pari_pointed_backend as backend
    import verify_cancellation_basis_bank_v2 as bank_module
    import cancellation_basis_certificate_epoch as engine
    from cysignals.signals import AlarmInterrupt
    protocol = guard(); job = next(r for r in protocol['jobs'] if r['id'] == label)
    dest = RAW/'jobs'/label; need(not dest.exists(),'preserve started integration')
    used = 0; injected = 0; old_execute = backend.execute; old_verify = bank_module.verify_selection

    def tape(search,mapping,height,seconds,expected_gp_hash):
        nonlocal used
        need(used < len(job['calls']),'engine requested an unretained query')
        row = read(ROOT/job['calls'][used]['file'])
        need(mapping == row['mapping'] and height == row['search']['height_bound'] and
             expected_gp_hash == row['search']['gp_binary_sha256'],'retained query binding differs')
        points = backend.replay(search,mapping,row['search']); used += 1
        return row['search'],points

    def verifying(packet,case,plan,folder,selection):
        nonlocal injected
        if job['alarm'] and len(packet['points']) > 18:
            need((folder/'selection.json').exists(),'alarm precedes retained failed producer')
            injected += 1; raise AlarmInterrupt('injected after enlarged producer, before independent certification')
        return old_verify(packet,case,plan,folder,selection)

    backend.execute = tape; bank_module.verify_selection = verifying
    try:
        engine.run(job['case'],job['arm'],protocol['plan'],dest,sha(OUT/'protocol.json'))
    except BaseException as error:
        new_write(RAW/'failures'/f'{label}.json',{'status':'RETAINED_INTEGRATION_FAILURE','job':label,
            'retained_calls_used':used,'alarm_injections':injected,'error':type(error).__name__+': '+str(error),
            'traceback':traceback.format_exc(),'protocol_sha256':sha(OUT/'protocol.json')})
        raise
    finally:
        backend.execute = old_execute; bank_module.verify_selection = old_verify
    result = read(dest/'result.json'); proof = read(dest/'independent-verification.json')
    need(used == len(job['calls']) and result['new_directions'] == job['expected_directions'] and
         result['later_cloud_directions'] == job['expected_later'] and
         [e['rank'] for e in result['epochs']] == job['expected_epochs'],'integration trace differs')
    need(proof['status'] == 'PASS' and result['unknowns'] == int(job['alarm']) and injected == int(job['alarm']), 'integration proof or fault differs')
    need(result['status'] == ('BANK_UNKNOWN' if job['alarm'] else 'CERTIFIED_TARGET'),'integration terminal state differs')
    if job['alarm']:
        need((dest/'epoch-01/producer/selection.json').exists() and not (dest/'epoch-01/verification.json').exists(),
             'failed enlarged bank not preserved as unverified')
    need(not list(dest.glob('epoch-*/reference')),'a reference landscape was fabricated')
    new_write(dest/'integration.json',{'status':'PASS_TRANSCRIPT_INTEGRATION','job':label,'retained_calls':used,
        'new_directions':result['new_directions'],'later_directions':result['later_cloud_directions'],
        'bank_ranks':[e['rank'] for e in result['epochs']],'alarm_injections':injected,
        'tail_witness_calls':proof['tail_witness_calls'],'result_sha256':sha(dest/'result.json'),
        'verification_sha256':sha(dest/'independent-verification.json'),'point_search_calls':0})


def run():
    protocol = guard(); need(not (OUT/'supervision.json').exists(),'preserve integration supervision')
    records = []
    for job in protocol['jobs']:
        label = job['id']; log_path = RAW/f'{label}.log'; need(not log_path.exists(),'preserve started integration log')
        before = resource.getrusage(resource.RUSAGE_CHILDREN); wall = time.monotonic()
        command = ['sage','-python',str(Path(__file__).resolve()),'worker','--job',label]
        with log_path.open('w') as log:
            child = subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try: code = child.wait(timeout=120)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGTERM)
                try: child.wait(timeout=2)
                except subprocess.TimeoutExpired: os.killpg(child.pid,signal.SIGKILL); child.wait()
                code = child.returncode
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        row = {'job':label,'returncode':code,'status':'COMPLETE' if code == 0 else 'RETAINED_FAILURE',
            'full_child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'wall_seconds':time.monotonic()-wall,'log_sha256':sha(log_path)}
        path = RAW/'jobs'/label/'integration.json'
        if path.exists(): row['integration_sha256'] = sha(path)
        records.append(row); atomic(OUT/'supervision.json',{'status':'RUNNING' if code == 0 else 'STOPPED_FAILURE',
            'records':records,'protocol_sha256':sha(OUT/'protocol.json'),'point_search_calls':0})
        print(json.dumps(row),flush=True)
        if code != 0: return
    need(sum(r['full_child_cpu_seconds'] for r in records) <= 240,'integration budget exceeded')
    atomic(OUT/'supervision.json',{'status':'COMPLETE','records':records,'protocol_sha256':sha(OUT/'protocol.json'),'point_search_calls':0})


def tamper():
    from verify_cancellation_basis_certificate import verify_arm
    protocol = guard(); need(read(OUT/'supervision.json')['status'] == 'COMPLETE','positive integrations incomplete')
    rows = []; start = time.process_time()
    for label in protocol['negative_checks']:
        source_job = 'double_refresh' if label == 'omitted_cloud_point' else 'single_refresh'
        job = next(j for j in protocol['jobs'] if j['id'] == source_job)
        dest = RAW/'negative'/label; shutil.copytree(RAW/'jobs'/source_job,dest)
        result = read(dest/'result.json'); events = read(dest/'events.json')
        if label == 'omitted_cloud_point':
            event = next(e for e in events if e['kind'] == 'call' and e['gain'])
            call = read(dest/event['file']); rank = read(dest/call['rank_file']); rank['packet']['points'].pop()
            write(dest/call['rank_file'],rank); call['rank_sha256'] = sha(dest/call['rank_file']); write(dest/event['file'],call)
            event['sha256'] = sha(dest/event['file'])
            next(c for c in result['calls'] if c['file'] == event['file'])['sha256'] = event['sha256']
            expected = 'cloud certificate drops or changes points'
        else:
            folder = dest/'epoch-01'; receipt = read(folder/'verification.json')
            if label == 'obsolete_seed':
                old = read(dest/'epoch-00/seed.json'); bank = read(folder/'bank.json')
                bank['seed'] = old; bank['rank'] = 18; write(folder/'seed.json',old); write(folder/'bank.json',bank)
                receipt['seed_sha256'] = sha(folder/'seed.json'); receipt['bank_sha256'] = sha(folder/'bank.json')
                expected = 'bank does not use complete current subgroup'
            else:
                proof = read(folder/'selection-proof.json'); proof['cvp_certificates'][0]['norm'] += 1
                write(folder/'selection-proof.json',proof); receipt['selection_proof_sha256'] = sha(folder/'selection-proof.json')
                expected = 'minimum witness packet differs'
            write(folder/'verification.json',receipt)
            event = next(e for e in events if e['kind'] == 'bank_ready' and e['epoch'] == 1)
            event['verification_sha256'] = sha(folder/'verification.json'); event['bank_sha256'] = sha(folder/'bank.json')
            epoch = next(e for e in result['epochs'] if e['epoch'] == 1)
            epoch['verification_sha256'] = event['verification_sha256']; epoch['bank_sha256'] = event['bank_sha256']
        write(dest/'events.json',events); result['events_sha256'] = sha(dest/'events.json'); write(dest/'result.json',result)
        try: verify_arm(job['case'],job['arm'],dest,protocol['plan'],result)
        except ArithmeticError as error:
            need(expected in str(error),'unexpected negative-check rejection: '+str(error))
            rows.append({'test':label,'status':'REJECTED_AS_REQUIRED','error':str(error),'result_sha256':sha(dest/'result.json')})
        else: raise ArithmeticError('corrupted integration accepted: '+label)
    new_write(OUT/'negative-checks.json',{'status':'PASS_UPDATED_HASH_CORRUPTION_CHECKS','rows':rows,
        'cpu_seconds':time.process_time()-start,'protocol_sha256':sha(OUT/'protocol.json'),'point_search_calls':0})
    records = read(OUT/'supervision.json')['records']; results = [read(RAW/'jobs'/j['id']/'integration.json') for j in protocol['jobs']]
    new_write(OUT/'summary.json',{'status':'PASS_COMPLETE_CERTIFICATE_ENGINE_INTEGRATION','jobs':results,
        'full_child_cpu_seconds':sum(r['full_child_cpu_seconds'] for r in records),
        'negative_check_component_cpu_seconds':time.process_time()-start,'negative_checks_sha256':sha(OUT/'negative-checks.json'),
        'protocol_sha256':sha(OUT/'protocol.json'),'point_search_calls':0,'boundary':protocol['boundary']})
    print(json.dumps(read(OUT/'summary.json'),indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command',choices=['seal','run','worker','tamper']); parser.add_argument('--job')
    args = parser.parse_args(); worker(args.job) if args.command == 'worker' else globals()[args.command]()
