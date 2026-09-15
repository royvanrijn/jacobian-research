#!/usr/bin/env sage -python
"""Serialization-corrected successor; the failed first protocol is preserved."""
import argparse
from collections import Counter, defaultdict
import fcntl
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
from cancellation_basis_epoch import read, sha, need, native_rank
from cancellation_scheduler_cpu import atomic
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, write, canonical, digest

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2'
RAW = LOCAL/'cancellation-basis-bank-certificate-v2'
ARMS = ['rational_rebuild', 'certificate_rebuild']


def seal():
    import gmpy2
    from sage.all import pari
    from sage.env import SAGE_VERSION
    import sys
    parent_plan, cases = early.guard()
    need(read(early.OUT/'supervision.json')['status'] == 'COMPLETE', 'retained experiment incomplete')
    need(read(OUT/'development.json')['status'] == 'PASS_RETAINED_BANK_COMMISSIONING', 'commissioning missing')
    need(read(OUT/'corruption-checks.json')['status'] == 'PASS_POSITIVE_AND_CORRUPTION_CHECKS', 'corruption checks missing')
    need(read(OUT/'binding-checks.json')['status'] == 'PASS_BOUND_SELECTION_SUCCESSOR', 'bound-selection successor missing')
    need(read(OUT/'binding-checks.json')['checker_sha256'] == sha(CAS/'verify_cancellation_basis_bank.py'), 'bound checker changed')
    need(read(OUT/'adapter-check.json')['status'] == 'PASS_SERIALIZATION_ADAPTER', 'live producer adapter uncommissioned')
    need(read(OUT/'adapter-check.json')['adapter_sha256'] == sha(CAS/'verify_cancellation_basis_bank_v2.py'), 'adapter changed')
    transitions = []; bindings = dict(read(OUT/'failed-predecessor.json')['input_sha256'])
    for case in cases:
        folder = early.RAW/'arms'/case['id']/'basis_refresh'
        result = read(folder/'result.json')
        need(read(folder/'independent-verification.json')['status'] == 'PASS', 'retained arm unverified')
        for name in ('result.json', 'independent-verification.json'):
            bindings[str((folder/name).relative_to(ROOT))] = sha(folder/name)
        for epoch in result['epochs']:
            if epoch['epoch'] == 0:
                continue
            old = folder/f'epoch-{epoch["epoch"]:02d}'
            receipt = read(old/'verification.json'); bank = read(old/'bank.json')
            need(receipt['status'] == 'PASS_INDEPENDENT_BANK' and sha(old/'verification.json') == epoch['verification_sha256'] and
                 sha(old/'bank.json') == epoch['bank_sha256'] == receipt['bank_sha256'], 'retained transition changed')
            packet = read(old/'seed.json'); need(packet == bank['seed'], 'retained seed binding differs')
            ident = case['id']+f'-e{epoch["epoch"]:02d}'
            transitions.append({'id':ident,'case':case,'epoch':epoch['epoch'],'rank':len(packet['points']),
                                'packet':packet,'retained':str(old.relative_to(ROOT))})
            for p in [old/'seed.json',old/'bank.json',old/'verification.json']+sorted((old/'producer').iterdir()):
                if p.is_file(): bindings[str(p.relative_to(ROOT))] = sha(p)
    need(len(transitions) == 34 and Counter(t['rank'] for t in transitions) == {19:17,20:17}, 'retained transition census differs')
    new_write(OUT/'inputs.json', transitions)
    env = {'sage':SAGE_VERSION,'pari':str(pari('version()')),'python':sys.version,'gmpy2':gmpy2.version(),
           'threads':{'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}}
    new_write(OUT/'environment.json', env)
    names = source_closure([Path(__file__),CAS/'verify_cancellation_basis_bank.py',CAS/'verify_cancellation_basis_bank_v2.py',
                            CAS/'verify_parity_minimum.py',CAS/'check_cancellation_basis_bank.py',
                            CAS/'check_cancellation_basis_bank_bindings.py',CAS/'check_cancellation_basis_bank_adapter.py'])
    plan = {'status':'FROZEN_ALL_RETAINED_TRANSITIONS_COST_GATE','arms':ARMS,'transitions':34,
        'primary_protocol_sha256':sha(early.OUT/'protocol.json'),
        'failed_predecessor':read(OUT/'failed-predecessor.json'),
        'adapter_commissioning_sha256':sha(OUT/'adapter-check.json'),
        'selection':'Every successful positive-cloud bank transition of the completed early-exposure candidate, in its frozen case order and increasing epoch. Seventeen M19 and seventeen M20 banks. No missing targets enter construction. This is a development engineering population, not a new search-validation population.',
        'hypothesis':'An independent exact fixed-radius certificate check can verify the identical complete selection more cheaply than repeating reference discovery.',
        'comparison':'Both isolated arms freshly verify the entire supplied subgroup using native Sage finite group arithmetic, then rebuild the same retained bank from the same seed and parent. The baseline executes the unmodified full producer/reference/native rebuild. The candidate executes the unmodified producer followed by the new independent fixed-radius, full-selection, native-word and prescribed-map checker. Both compare the complete rebuilt selection and exported bank against the same retained output only after rebuilding. No reference output is used to construct either bank.',
        'bank_plan':{k:parent_plan[k] for k in ('bank_policy','maximum_extension_dimension','maximum_centres')},
        'process_cpu_soft_seconds':30,'process_cpu_hard_seconds':35,'process_wall_seconds':90,
        'campaign_child_cpu_seconds':600,'maximum_processes':68,'point_search_calls':0,
        'execution':'Sequential isolated workers; alternate the two arm orders by transition ordinal. Charge complete child CPU including imports, source/input checks, native rank, construction, verification, comparisons and serialization. Report supervisor CPU separately. Reserve both full hard limits before a new pair. Retain every started job; no failed or unreceipted worker retry, input substitution, budget extension, or outcome-dependent roster change.',
        'engineering_gate':'All 68 workers must finish and certify identical complete selections and exported banks. Aggregate full-child CPU saving must be at least25 percent; the central97.5-percent paired whole-curve family-stratified bootstrap interval for baseline/candidate CPU must have lower endpoint above1, with10000 draws seed20260914. Both M19 and M20 strata must have baseline/candidate CPU ratio above1. An incomplete pair or any certificate failure blocks the engineering gate.',
        'bootstrap':{'draws':10000,'seed':20260914,'quantiles':[0.0125,0.9875]},
        'boundary':'Implementation-cost evidence on retained known-control transitions only. Rounded numerical heights remain a scoring metric. No point search, later-gain rate, new-rank discovery, whole-coset exclusion or external state-of-the-art claim. This gate does not alter the failed charged search comparisons. A future search integration needs its own honest certificate receipt handling and newly sealed disjoint comparison.',
        'source_sha256':names,'retained_input_sha256':bindings,
        'input_sha256':{n:sha(OUT/n) for n in ('inputs.json','environment.json','development.json','corruption-checks.json','binding-checks.json','adapter-check.json','failed-predecessor.json')}}
    new_write(OUT/'protocol.json', plan)
    for name in names:
        target = RAW/'sources'/name
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/name,target)
    print(json.dumps({'status':plan['status'],'protocol_sha256':sha(OUT/'protocol.json'),
                      'sources':len(names),'transitions':len(transitions),'rank_counts':dict(Counter(t['rank'] for t in transitions))}), flush=True)


def guard():
    plan = read(OUT/'protocol.json')
    for name,h in plan['source_sha256'].items(): need(sha(ROOT/name) == h, 'sealed source changed: '+name)
    for name,h in plan['input_sha256'].items(): need(sha(OUT/name) == h, 'sealed input changed: '+name)
    for name,h in plan['retained_input_sha256'].items(): need(sha(ROOT/name) == h, 'retained input changed: '+name)
    need(sha(early.OUT/'protocol.json') == plan['primary_protocol_sha256'], 'primary protocol changed')
    return plan, read(OUT/'inputs.json')


def worker(ident, arm):
    # The supervisor's full-child meter starts before this Python process.
    resource.setrlimit(resource.RLIMIT_CPU, (30,35))
    started = time.process_time(); plan, inputs = guard()
    need(arm in plan['arms'], 'unsealed arm')
    item = next(x for x in inputs if x['id'] == ident); dest = RAW/'arms'/ident/arm
    need(not (dest/'result.json').exists() and not (dest/'failure.json').exists(), 'preserve started worker')
    try:
        from run_complement_seed_v3 import normalized_selection
        if arm == 'rational_rebuild':
            from cancellation_basis_epoch import rebuild
        else:
            from verify_cancellation_basis_bank_v2 import rebuild
        tick = time.process_time(); rank = native_rank(item['packet'])
        need(rank['rank'] == item['rank'], 'subgroup rank differs')
        write(dest/'native-rank.json', {'status':'PASS','certificate':rank,'seed_sha256':digest(canonical(item['packet']))})
        rank_cpu = time.process_time()-tick
        bank, receipt = rebuild(item['packet'],item['case'],plan['bank_plan'],dest/'bank')
        retained = ROOT/item['retained']; tick = time.process_time()
        selection = read(dest/'bank/producer/selection.json')
        need(canonical(bank) == canonical(read(retained/'bank.json')) and sha(dest/'bank/bank.json') == sha(retained/'bank.json'), 'rebuilt exported bank differs from retained bank')
        need(normalized_selection(selection) == normalized_selection(read(retained/'producer/selection.json')),
             'complete rebuilt selection differs from retained selection')
        result = {'status':'PASS_IDENTICAL_CERTIFIED_BANK','id':ident,'arm':arm,'rank':item['rank'],
            'case':item['case']['id'],'family':item['case']['family'],'epoch':item['epoch'],
            'native_rank_sha256':sha(dest/'native-rank.json'),'verification_sha256':sha(dest/'bank/verification.json'),
            'bank_sha256':sha(dest/'bank/bank.json'),'retained_bank_sha256':sha(retained/'bank.json'),
            'complete_centres':len(selection['centres']),'exported_centres':len(bank['centres']),
            'cvp_packets':sum(len(a['refined']) for a in selection['anchors']),
            'components':{'native_rank_cpu_seconds':rank_cpu,'comparison_cpu_seconds':time.process_time()-tick,
                          'bank_receipt':receipt},'protocol_sha256':sha(OUT/'protocol.json'),
            'internal_cpu_seconds':time.process_time()-started,'point_search_calls':0}
        new_write(dest/'result.json', result)
    except BaseException as error:
        new_write(dest/'failure.json', {'status':'RETAINED_WORKER_FAILURE','id':ident,'arm':arm,
            'error':type(error).__name__+': '+str(error),'traceback':traceback.format_exc(),
            'internal_cpu_seconds':time.process_time()-started,'protocol_sha256':sha(OUT/'protocol.json')})
        raise


def run():
    plan, inputs = guard(); started = time.process_time(); ph = sha(OUT/'protocol.json')
    RAW.mkdir(parents=True, exist_ok=True)
    with (RAW/'supervisor.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        previous = read(OUT/'supervision.json') if (OUT/'supervision.json').exists() else {}
        if previous.get('status') in ('COMPLETE','STOPPED_FAILURE','STOPPED_BUDGET','STOPPED_UNRECEIPTED'):
            print(json.dumps({'status':previous['status'],'action':'RETAINED_TERMINAL_NO_RETRY'}), flush=True); return
        records = []; old_supervisor = previous.get('supervisor_cpu_seconds',0)

        def checkpoint(status):
            atomic(OUT/'supervision.json', {'status':status,'records':records,'protocol_sha256':ph,
                'charged_child_cpu_seconds':sum(r['charged_cpu_seconds'] for r in records),
                'supervisor_cpu_seconds':old_supervisor+time.process_time()-started})

        def child_cpu():
            r = resource.getrusage(resource.RUSAGE_CHILDREN); return r.ru_utime+r.ru_stime

        for ordinal,item in enumerate(inputs):
            order = plan['arms'][::1 if ordinal % 2 == 0 else -1]
            pending = [arm for arm in order if not (RAW/'arms'/item['id']/arm/'supervisor.json').exists()]
            charged = sum(r['charged_cpu_seconds'] for r in records)
            charged += sum(read(RAW/'arms'/item['id']/arm/'supervisor.json')['charged_cpu_seconds']
                           for arm in order if arm not in pending)
            if charged+len(pending)*plan['process_cpu_hard_seconds'] > plan['campaign_child_cpu_seconds']:
                checkpoint('STOPPED_BUDGET'); return
            for arm in order:
                dest = RAW/'arms'/item['id']/arm; dest.mkdir(parents=True,exist_ok=True)
                if (dest/'supervisor.json').exists():
                    receipt = read(dest/'supervisor.json'); need(receipt['protocol_sha256'] == ph,'resume protocol differs')
                    records.append(receipt)
                    if receipt['status'] != 'COMPLETE': checkpoint('STOPPED_FAILURE'); return
                    continue
                if (dest/'start.json').exists() or (dest/'result.json').exists():
                    checkpoint('STOPPED_UNRECEIPTED'); return
                command = ['sage','-python',str(Path(__file__).resolve()),'worker','--id',item['id'],'--arm',arm]
                new_write(dest/'start.json',{'command':command,'protocol_sha256':ph,'ordinal':ordinal})
                before = child_cpu(); wall = time.monotonic()
                env = {**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
                with (dest/'worker.log').open('w') as log:
                    process = subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
                    try:
                        rc = process.wait(timeout=plan['process_wall_seconds']); status = 'COMPLETE' if rc == 0 else 'WORKER_FAILURE'
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGTERM)
                        try: process.wait(timeout=2)
                        except subprocess.TimeoutExpired: os.killpg(process.pid,signal.SIGKILL); process.wait()
                        rc = process.returncode; status = 'OUTER_WALL_TIMEOUT'
                receipt = {'id':item['id'],'case':item['case']['id'],'family':item['case']['family'],'rank':item['rank'],
                    'epoch':item['epoch'],'arm':arm,'status':status,'returncode':rc,'charged_cpu_seconds':child_cpu()-before,
                    'wall_seconds':time.monotonic()-wall,'protocol_sha256':ph,'worker_log_sha256':sha(dest/'worker.log')}
                if (dest/'result.json').exists(): receipt['result_sha256'] = sha(dest/'result.json')
                if (dest/'failure.json').exists(): receipt['failure_sha256'] = sha(dest/'failure.json')
                new_write(dest/'supervisor.json',receipt); records.append(receipt); checkpoint('RUNNING')
                print(json.dumps({'workers':len(records),'total':68,'id':item['id'],'arm':arm,
                                  'status':status,'cpu_seconds':receipt['charged_cpu_seconds']}), flush=True)
                if status != 'COMPLETE': checkpoint('STOPPED_FAILURE'); return
        checkpoint('COMPLETE')


def report():
    import numpy as np
    plan, inputs = guard(); supervision = read(OUT/'supervision.json'); ph = sha(OUT/'protocol.json')
    need(supervision['status'] == 'COMPLETE' and len(supervision['records']) == 68, 'comparison incomplete; no engineering gate')
    rows = []; totals = defaultdict(float); groups = defaultdict(dict)
    expected = {(t['id'],a) for t in inputs for a in plan['arms']}; seen = set()
    for receipt in supervision['records']:
        key = (receipt['id'],receipt['arm']); need(key in expected and key not in seen, 'worker census differs'); seen.add(key)
        dest = RAW/'arms'/receipt['id']/receipt['arm']; result = read(dest/'result.json')
        need(receipt['protocol_sha256'] == ph == result['protocol_sha256'] and
             receipt['result_sha256'] == sha(dest/'result.json') and receipt['status'] == 'COMPLETE' and
             result['status'] == 'PASS_IDENTICAL_CERTIFIED_BANK', 'worker verification binding differs')
        need(all(result[k] == receipt[k] for k in ('id','arm','rank','case','family','epoch')) and
             receipt['worker_log_sha256'] == sha(dest/'worker.log') and receipt['charged_cpu_seconds'] > 0,
             'worker identity or meter binding differs')
        need(result['bank_sha256'] == result['retained_bank_sha256'] == sha(dest/'bank/bank.json') and
             result['verification_sha256'] == sha(dest/'bank/verification.json') and
             result['native_rank_sha256'] == sha(dest/'native-rank.json'), 'bank receipt differs')
        rows.append({**receipt,'components':result['components'],'cvp_packets':result['cvp_packets']})
        totals[receipt['arm']] += receipt['charged_cpu_seconds']
        vector = groups[receipt['family']].setdefault(receipt['case'],[0.,0.])
        vector[plan['arms'].index(receipt['arm'])] += receipt['charged_cpu_seconds']
    need(seen == expected,'missing worker'); base,candidate = [totals[a] for a in plan['arms']]
    rng = np.random.default_rng(plan['bootstrap']['seed']); draws = []
    arrays = [np.array(list(cases.values())) for _,cases in sorted(groups.items())]
    for _ in range(plan['bootstrap']['draws']):
        samples = np.zeros(2)
        for values in arrays: samples += values[rng.integers(len(values),size=len(values))].sum(axis=0)
        draws.append(samples[0]/samples[1])
    interval = list(map(float,np.quantile(draws,plan['bootstrap']['quantiles'])))
    strata = {}
    for rank in (19,20):
        sums = {a:sum(r['charged_cpu_seconds'] for r in rows if r['arm']==a and r['rank']==rank) for a in plan['arms']}
        strata[str(rank)] = {'cpu_seconds':sums,'baseline_over_candidate':sums[plan['arms'][0]]/sums[plan['arms'][1]]}
    saving = 1-candidate/base
    gate = saving >= .25 and interval[0] > 1 and all(s['baseline_over_candidate'] > 1 for s in strata.values())
    summary = {'status':'COMPLETE_IDENTICAL_BANK_COST_COMPARISON','protocol_sha256':ph,'rows':rows,
        'workers':68,'transitions':34,'whole_curve_groups':sum(map(len,groups.values())),
        'family_curve_counts':{f:len(v) for f,v in groups.items()},'full_child_cpu_seconds':dict(totals),
        'baseline_over_candidate_cpu':base/candidate,'candidate_cpu_saving_fraction':saving,
        'paired_family_whole_curve_bootstrap_central_97_5_percent':interval,'rank_strata':strata,
        'engineering_gate':bool(gate),'supervisor_cpu_seconds':supervision['supervisor_cpu_seconds'],
        'development_component_cpu_seconds':read(OUT/'development.json')['cpu_seconds'],
        'corruption_check_component_cpu_seconds':read(OUT/'corruption-checks.json')['cpu_seconds'],
        'binding_check_component_cpu_seconds':read(OUT/'binding-checks.json')['cpu_seconds'],
        'adapter_full_child_cpu_seconds':read(OUT/'adapter-check.json')['full_child_cpu_seconds'],
        'failed_predecessor_full_child_cpu_seconds':plan['failed_predecessor']['charged_child_cpu_seconds'],
        'saving_charging_recorded_development_and_failed_predecessor_to_candidate':1-(candidate+read(OUT/'development.json')['cpu_seconds']+read(OUT/'corruption-checks.json')['cpu_seconds']+read(OUT/'binding-checks.json')['cpu_seconds']+read(OUT/'adapter-check.json')['full_child_cpu_seconds']+plan['failed_predecessor']['charged_child_cpu_seconds'])/base,
        'point_search_calls':0,'later_gain_performance':'UNTESTED_BY_THIS_COST_GATE','boundary':plan['boundary']}
    new_write(OUT/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('rows','boundary')},indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command',choices=['seal','run','worker','report'])
    parser.add_argument('--id'); parser.add_argument('--arm'); args = parser.parse_args()
    worker(args.id,args.arm) if args.command == 'worker' else globals()[args.command]()
