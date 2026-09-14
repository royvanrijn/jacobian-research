#!/usr/bin/env sage -python
"""Independent successor round: all declared jobs, unchanged per-arm CPU.

Reuses the source-frozen V3 arithmetic and worker. Configuration is explicit,
versioned and hashed; no V3 input, result, fit or script is overwritten.
"""
import argparse
import fcntl
import gzip
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

import cancellation_scheduler_prepare as preparation
from finite_cancellation_corpus import ROOT, LOCAL, OUT as CORPUS, canonical, digest, write

OUT=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4'
RAW=LOCAL/'cancellation-scheduler-v4'
PREVIOUS=preparation.OUT
CAS=Path(__file__).resolve().parent
ENTRY=Path(__file__).resolve()


def prepare():
    previous=json.loads((PREVIOUS/'summary.json').read_text())
    if previous['status']!='COMPLETE_FROZEN_ADAPTIVE_COMPARISON':raise ArithmeticError('previous round incomplete')
    totals=previous['scopes']['all']['totals']
    candidate=max(('adaptive_local','adaptive_uniform'),key=lambda a:(totals[a]['gains']/totals[a]['cpu_seconds'],a))
    baseline=totals['factor_free']
    relative=(totals[candidate]['gains']/totals[candidate]['cpu_seconds'])/(baseline['gains']/baseline['cpu_seconds'])
    if relative<1.25:raise ArithmeticError('no material development-rate reason for successor')
    receipts=json.loads((PREVIOUS/'supervision.json').read_text())['records']
    cap_stops=[]
    for receipt in receipts:
        if receipt['arm']!=candidate:continue
        path=preparation.RAW/'arms'/receipt['case']/candidate/'result.json'
        row=json.loads(path.read_text())
        if row['status']=='CALL_CAP' and row['search_phase_cpu_seconds']<40:
            cap_stops.append({'case':receipt['case'],'search_cpu_seconds':row['search_phase_cpu_seconds'],
                              'result_sha256':digest(path.read_bytes())})
    if not cap_stops:raise ArithmeticError('call-cap failure mechanism not present')
    design=json.loads((PREVIOUS/'design.json').read_text())
    design.update(status='SUCCESSOR_DESIGN_BEFORE_NEW_OUTCOMES',arms=['factor_free',candidate],
        maximum_calls=432,
        rationale='V3 is preserved as a failed primary comparison. Its best observed recovery/CPU adaptive arm becomes the single development-selected candidate. Permit all 48 anchors times 3 models times 3 heights, instead of terminating at 144 calls below the CPU budget. Keep exactly the same 40-CPU-second search allowance and H125000 ceiling. Do not rerun any previous CPU curve.',
        previous_summary_sha256=digest((PREVIOUS/'summary.json').read_bytes()),
        candidate=candidate,development_relative_rate=relative,development_cap_stops=cap_stops,
        selection='24 previously untested shallow j groups, four per each of six families, from the retained eligibility audit in SHA256(cancellation-scheduler-v4/ + j) order; lexical first bank per j. Exclude all V1, V2 and V3 CPU curves and Curve302. Refit the identical mixture after excluding these 24 groups.',
        primary_gate='One preselected candidate versus factor-free V3: at least as many independently certified directions, at least 10 percent higher recoveries per complete CPU, and central 97.5 percent paired family bootstrap interval wholly above one. 10000 draws, seed 20260914. Success is a finite policy validation, not cancellation-specific attribution.',
        fresh_gate='A passing new withheld comparison permits a separately frozen high-rank transfer gate on retained M27/M29/M30 control banks. Fresh-fibre point search still requires that transfer gate plus generic-only new inputs, full cold-construction costs and a separate finite paired budget.',
        boundary='Development-selected successor evaluated on disjoint whole-j CPU controls. No retroactive promotion of V3, no larger CPU budget, and no universal non-predictiveness or rank32 claim.')
    preparation.new_write(OUT/'design.json',design)
    old_roster=json.loads((PREVIOUS/'roster.json').read_text());excluded=set(old_roster['excluded_training_j_groups'])
    candidates=json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes()))
    banks={}
    for row in sorted(candidates,key=lambda r:r['source_landscape']):
        if row['stratum']=='shallow' and row['j_group'] not in excluded:
            banks.setdefault(row['j_group'],row)
    chosen=[]
    for family in preparation.FAMILIES:
        rows=sorted((r for r in banks.values() if r['family']==family),
                    key=lambda r:digest(('cancellation-scheduler-v4/'+r['j_group']).encode()))
        chosen.extend(dict(r) for r in rows[:4])
    if len(chosen)!=24:raise ArithmeticError('missing successor controls; no replacement')
    corpus={c['id']:c for c in json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))}
    inputs=[];oracle=[]
    for row in chosen:
        for key,hkey in [('source_seed','seed_sha256'),('source_landscape','landscape_sha256')]:
            if digest((ROOT/row[key]).read_bytes())!=row[hkey]:raise ArithmeticError('retained input changed')
        seed=json.loads((ROOT/row['source_seed']).read_text());landscape=json.loads((ROOT/row['source_landscape']).read_text())
        if seed['points']!=landscape['basis']:raise ArithmeticError('basis differs')
        ident=digest(canonical(['v4',row['j_group'],row['source_landscape']]))[:20];row['id']=ident
        inputs.append({'id':ident,'j_group':row['j_group'],'family':row['family'],'stratum':row['stratum'],
            'seed':{k:seed[k] for k in ('curve','points','proof')},
            'centres':[c['representative'] for c in landscape['centres'][:48]]})
        endpoint=corpus[row['corpus_id']]
        oracle.append({'id':ident,'curve':endpoint['curve'],'points':endpoint['generic_points']+endpoint['targets'],
            'source':endpoint['source'],'proof_sha256':endpoint['proof_sha256']})
    preparation.new_write(OUT/'inputs.json',inputs);preparation.new_write(OUT/'oracle.json',oracle)
    preparation.new_write(OUT/'roster.json',{'chosen':chosen,'prior_cpu_j_groups':sorted(excluded),
        'excluded_training_j_groups':sorted(excluded|{c['j_group'] for c in inputs}),
        'inputs_sha256':digest((OUT/'inputs.json').read_bytes()),'oracle_sha256':digest((OUT/'oracle.json').read_bytes()),
        'previous_roster_sha256':digest((PREVIOUS/'roster.json').read_bytes())})
    print(json.dumps({'status':'FROZEN_SUCCESSOR_ROSTER','cases':len(inputs),'candidate':candidate}),flush=True)


def train():
    preparation.OUT=OUT;preparation.train()


def preflight():
    preparation.OUT=OUT;preparation.preflight()


def training_replay():
    import verify_cancellation_scheduler_training as verifier
    verifier.OUT=OUT;verifier.main()


def seal():
    plan=json.loads((OUT/'design.json').read_text());fit=json.loads((OUT/'fit.json').read_text())
    pre=json.loads((OUT/'preflight.json').read_text());replay=json.loads((OUT/'training-replay.json').read_text())
    if pre['status']!='PASS' or replay['status']!='PASS' or replay['fit_sha256']!=digest((OUT/'fit.json').read_bytes()):
        raise ArithmeticError('preflight/replay prerequisite')
    inputs=json.loads((OUT/'inputs.json').read_text())
    if set(fit['training_j_groups'])&{c['j_group'] for c in inputs}:raise ArithmeticError('fit leakage')
    old=json.loads((PREVIOUS/'protocol.json').read_text())
    sources=dict(old['source_sha256'])
    sources.update(json.loads((PREVIOUS/'dependency-audit.json').read_text())['extra_source_sha256'])
    for name,h in sources.items():
        if digest((ROOT/name).read_bytes())!=h:raise ArithmeticError('V3 arithmetic source changed')
    sources[str(Path(__file__).relative_to(ROOT))]=digest(Path(__file__).read_bytes())
    plan.update(status='FROZEN_SUCCESSOR_BEFORE_CPU',cases=len(inputs),source_sha256=sources,
        gp_sha256=old['gp_sha256'],fit={k:fit[k] for k in old['fit']},
        input_sha256={n:digest((OUT/n).read_bytes()) for n in ('design.json','inputs.json','oracle.json','roster.json','fit.json','preflight.json','training-replay.json')},
        worker_configuration={'OUT':str(OUT.relative_to(ROOT)),'RAW':str(RAW.relative_to(ROOT)),
            'implementation':'Unmodified V3 cancellation_scheduler_cpu.run; explicit versioned input/output module constants set by this sealed entry point.'})
    preparation.new_write(OUT/'protocol.json',plan);print('SEALED',digest((OUT/'protocol.json').read_bytes()),flush=True)


def worker(case,arm):
    import cancellation_scheduler_cpu as cpu_worker
    cpu_worker.OUT=OUT;cpu_worker.RAW=RAW
    cpu_worker.run(case,arm)


def run():
    RAW.mkdir(parents=True,exist_ok=True)
    with (RAW/'supervisor.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
        for name,h in plan['source_sha256'].items():
            if digest((ROOT/name).read_bytes())!=h:raise ArithmeticError('source changed')
        inputs=json.loads((OUT/'inputs.json').read_text());records=[]
        if digest((OUT/'inputs.json').read_bytes())!=plan['input_sha256']['inputs.json']:raise ArithmeticError('input changed')
        def cpu():
            r=resource.getrusage(resource.RUSAGE_CHILDREN);return r.ru_utime+r.ru_stime
        for i,case in enumerate(inputs):
            for arm in plan['arms'][::1 if i%2==0 else -1]:
                dest=RAW/'arms'/case['id']/arm;dest.mkdir(parents=True,exist_ok=True)
                if (dest/'supervisor.json').exists():
                    receipt=json.loads((dest/'supervisor.json').read_text())
                    if receipt['protocol_sha256']!=ph:raise ArithmeticError('resume protocol differs')
                    records.append(receipt);continue
                if (dest/'start.json').exists() or (dest/'result.json').exists():raise ArithmeticError('unreceipted arm retained')
                command=['sage','-python',str(ENTRY), 'worker','--case',case['id'],'--arm',arm]
                write(dest/'start.json',{'command':command,'protocol_sha256':ph})
                before=cpu();wall=time.monotonic()
                with (dest/'worker.log').open('w') as log:
                    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                    try:
                        rc=process.wait(timeout=plan['arm_wall_seconds']);status='COMPLETE' if rc==0 else 'WORKER_FAILURE'
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGTERM)
                        try:process.wait(timeout=2)
                        except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                        rc=process.returncode;status='OUTER_WALL_TIMEOUT'
                receipt={'case':case['id'],'arm':arm,'family':case['family'],'stratum':case['stratum'],
                    'status':status,'returncode':rc,'charged_cpu_seconds':cpu()-before,
                    'wall_seconds':time.monotonic()-wall,'protocol_sha256':ph}
                result={}
                if (dest/'result.json').exists():
                    receipt['result_sha256']=digest((dest/'result.json').read_bytes());result=json.loads((dest/'result.json').read_text())
                write(dest/'supervisor.json',receipt);records.append(receipt)
                write(OUT/'supervision.json',{'status':'RUNNING','protocol_sha256':ph,'records':records})
                print(json.dumps({'arms':len(records),'total':2*len(inputs),'case':case['id'],'arm':arm,
                    'status':status,'gain':result.get('success'),'calls':len(result.get('calls',[])),
                    'cpu':receipt['charged_cpu_seconds']}),flush=True)
                if status!='COMPLETE':
                    write(OUT/'supervision.json',{'status':'STOPPED_INFRASTRUCTURE_UNKNOWN','records':records,'protocol_sha256':ph});return
        write(OUT/'supervision.json',{'status':'COMPLETE','records':records,'protocol_sha256':ph})


def report():
    import report_cancellation_scheduler as reporting
    plan=json.loads((OUT/'protocol.json').read_text());reporting.ARMS=plan['arms']
    supervision=json.loads((OUT/'supervision.json').read_text())
    if supervision['status']!='COMPLETE' or len(supervision['records'])!=2*plan['cases']:raise ArithmeticError('incomplete successor')
    rows=[];unknowns=0
    for receipt in supervision['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm'];result=json.loads((dest/'result.json').read_text())
        verification=json.loads((dest/'independent-verification.json').read_text())
        if receipt['status']!='COMPLETE' or receipt['result_sha256']!=digest((dest/'result.json').read_bytes()):raise ArithmeticError('bad receipt')
        if verification['status']!='PASS' or result['independent_verification_sha256']!=digest((dest/'independent-verification.json').read_bytes()):raise ArithmeticError('bad proof receipt')
        unknowns+=verification['preparation_unknowns']
        rows.append({'case':receipt['case'],'arm':receipt['arm'],'family':receipt['family'],'stratum':receipt['stratum'],
            'success':result['success'],'cpu_seconds':receipt['charged_cpu_seconds'],'calls':len(result['calls']),
            'literal_hit':False,'rank_lower_bound':result['rank_lower_bound'],'initial_rank':result['initial_rank'],
            'status':result['status']})
    comparison=reporting.compare(rows,plan['candidate'],'factor_free');totals=reporting.summarize(rows)
    for aggregate in totals.values():aggregate.pop('literal_hit_cases')
    for row in rows:row.pop('literal_hit')
    offline={n:json.loads((OUT/n).read_text())['cpu_seconds'] for n in ('fit.json','preflight.json','training-replay.json')}
    a,b=totals[plan['candidate']],totals['factor_free']
    ratio=(a['gains']/(a['cpu_seconds']+sum(offline.values())))/(b['gains']/b['cpu_seconds']) if b['gains'] else None
    write(OUT/'summary.json',{'status':'COMPLETE_INDEPENDENT_SUCCESSOR','candidate':plan['candidate'],'rows':rows,
        'totals':totals,'comparison':comparison,'offline_cpu_seconds':offline,
        'conservative_campaign_rate_ratio':ratio,'high_rank_transfer_gate':comparison['gate_passed'] is True and unknowns==0,
        'fresh_fibres_run':0,'rank32':'UNKNOWN','protocol_sha256':digest((OUT/'protocol.json').read_bytes()),
        'literal_oracle_hits':'NOT_ANALYSED',
        'boundary':'24 new whole-j controls; inherited V3 outcomes used only for development choice and rationale. Exact same CPU and height budgets, all possible declared jobs permitted. Literal oracle-hit analysis not computed here. No rank record or fresh-fibre result.'})
    print(json.dumps({'totals':totals,'comparison':comparison,'conservative_campaign_rate_ratio':ratio},indent=2),flush=True)


def pack():
    import pack_cancellation_scheduler as packing
    packing.OUT=OUT;packing.RAW=RAW;packing.pack()


def replay():
    import verify_cancellation_scheduler as verifier
    verifier.OUT=OUT;verifier.RAW=RAW;verifier.main()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','train','preflight','training_replay','seal','run','worker','report','pack','replay'])
    parser.add_argument('--case');parser.add_argument('--arm');args=parser.parse_args()
    worker(args.case,args.arm) if args.command=='worker' else globals()[args.command]()
