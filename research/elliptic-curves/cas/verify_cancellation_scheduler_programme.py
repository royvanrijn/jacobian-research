#!/usr/bin/env sage -python
"""Audit all four scheduler stages; optional arithmetic replay makes no searches.

The default checks retained independent receipts and complete byte bundles.
It does not describe that integrity audit as a new arithmetic proof.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

from finite_cancellation_corpus import ROOT, LOCAL, digest, write

OUT=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_programme_v1'
STAGES=(('cancellation_scheduler_v3','cancellation-scheduler-v3',123),
        ('cancellation_scheduler_v4','cancellation-scheduler-v4',48),
        ('cancellation_scheduler_transfer_v1','cancellation-scheduler-transfer-v1',6),
        ('cancellation_scheduler_fresh_v1','cancellation-scheduler-fresh-v1',8))


def read(path):return json.loads(path.read_text())
def sha(path):return digest(path.read_bytes())
def need(condition,message):
    if not condition:raise ArithmeticError(message)


def audit():
    import pack_cancellation_scheduler as bundle
    rows=[];bindings={};groups=[];summaries=[];start=time.process_time()
    for name,rawname,expected in STAGES:
        folder=OUT.parent/name;raw=LOCAL/rawname
        plan=read(folder/'protocol.json');ph=sha(folder/'protocol.json')
        supervisor=read(folder/'supervision.json');summary=read(folder/'summary.json')
        need(supervisor['status']=='COMPLETE' and len(supervisor['records'])==expected,'incomplete stage')
        need(summary['protocol_sha256']==ph==supervisor['protocol_sha256'],'protocol link differs')
        for path,h in plan['source_sha256'].items():need(sha(ROOT/path)==h,'source/input changed: '+path)
        for path,h in plan['input_sha256'].items():need(sha(folder/path)==h,'frozen input differs')
        if name.endswith(('v3','v4')):
            inputs=read(folder/'inputs.json');test={c['j_group'] for c in inputs}
            fit=read(folder/'fit.json');need(not test&set(fit['training_j_groups']),'fit/test overlap')
            need(all(not test&g for g in groups),'repeated CPU holdout')
            groups.append(test)
        totals={a:{'gains':0,'cpu_seconds':0.,'calls':0,'cloud_directions':0} for a in plan['arms']}
        statuses=Counter();certified=0;unknowns=0;seen=set()
        for receipt in supervisor['records']:
            key=receipt['case'],receipt['arm'];need(key not in seen,'duplicate arm');seen.add(key)
            dest=raw/'arms'/key[0]/key[1];result=read(dest/'result.json')
            need(receipt==read(dest/'supervisor.json'),'supervisor receipt differs')
            need(receipt['status']=='COMPLETE' and receipt['result_sha256']==sha(dest/'result.json'),'result hash differs')
            total=totals[key[1]];total['cpu_seconds']+=receipt['charged_cpu_seconds']
            total['calls']+=len(result['calls']);total['gains']+=bool(result['success'])
            if result['rank_lower_bound'] is None:
                need(not result['success'] and not result['calls'],'invalid unresolved result');unknowns+=1;continue
            verification=read(dest/'independent-verification.json')
            need(verification['status']=='PASS' and result['independent_verification_sha256']==sha(dest/'independent-verification.json'),'independent proof receipt differs')
            need(verification['rank']['rank']==result['rank_lower_bound']==result['initial_rank']+int(result['success']),'rank accounting differs')
            certified+=result['success'];unknowns+=verification['preparation_unknowns']
            for call in result['calls']:
                need(call['sha256']==sha(dest/call['file']),'point call bytes differ');statuses[call['status']]+=1
            if (dest/'cloud-verification.json').exists():
                cloud=read(dest/'cloud-verification.json')
                need(cloud['certificate_sha256']==sha(dest/'cloud-certificate.json'),'cloud hash differs')
                need(cloud['result_sha256']==sha(dest/'result.json') and cloud['point_search_calls_added']==0,'cloud/search boundary differs')
                total['cloud_directions']+=cloud['rank']['rank']-result['initial_rank']
        reported_totals=summary.get('totals',summary.get('scopes',{}).get('all',{}).get('totals'))
        if reported_totals is not None:
            for arm,total in totals.items():
                for key in ('gains','cpu_seconds','calls'):
                    need(abs(reported_totals[arm][key]-total[key])<1e-9,'aggregate differs')
        if 'cpu_totals' in summary:
            for arm,total in totals.items():need(abs(summary['cpu_totals'][arm]-total['cpu_seconds'])<1e-9,'transfer aggregate differs')
        bundle.OUT=folder;bundle.RAW=raw;bundle.verify()
        for path in (folder/'protocol.json',folder/'summary.json',folder/'replay-manifest.json'):
            bindings[str(path.relative_to(ROOT))]=sha(path)
        rows.append({'stage':name,'arms':expected,'certified_first_directions':certified,
            'unknowns':unknowns,'point_statuses':dict(statuses),'totals':totals})
        summaries.append(summary)
    need(summaries[0]['promote_local_policy'] is False,'original failed gate was changed')
    need(summaries[1]['high_rank_transfer_gate'] is True,'V4 gate not passed')
    need(summaries[2]['fresh_fibre_gate'] is True,'transfer gate not passed')
    fresh=summaries[3]
    need(all(r.get('absent_retained_corpus') for r in fresh['rows']),'freshness not established')
    need(len({r['j_group'] for r in fresh['rows']})==4,'fresh fibres collide')
    selected=read(OUT.parent/'cancellation_scheduler_v4/protocol.json')
    write(OUT/'policy.json',{'schema':'cancellation-next-direction-policy-v1',
        'status':'VALIDATED_FINITE_CANDIDATE','arm':selected['candidate'],
        'implementation':'elliptic-curves/cas/cancellation_scheduler.py',
        'worker':'elliptic-curves/cas/cancellation_scheduler_cpu.py',
        'fit':selected['fit'],
        'profile':{k:selected[k] for k in ('heights','maximum_centres','maximum_calls','search_cpu_seconds','point_wall_seconds','gp_sha256')},
        'input_contract':'Certified short-model subgroup and compatible ordered centre words; pre-freeze recipe, source hashes, isolated limits and independent certificate endpoint. The worker accepts versioned OUT/RAW module configuration as in the retained stage entry points.',
        'selected_protocol_sha256':sha(OUT.parent/'cancellation_scheduler_v4/protocol.json'),
        'boundary':'Optimizes first independent direction, not final rank or total cloud yield. Four fresh fibres give only 4.25 percent cold saving and fewer cloud directions than V3. Keep a paired baseline and new held-out validation for extensions; rebuild banks after basis changes. No automatic search or global default change.'})
    write(OUT/'audit.json',{'status':'PASS_RETAINED_RECEIPTS_AND_BUNDLES','stages':rows,
        'bindings_sha256':bindings,'checker_sha256':sha(Path(__file__)),
        'policy_sha256':sha(OUT/'policy.json'),
        'cpu_seconds':time.process_time()-start,'rank32':'UNKNOWN',
        'boundary':'Integrity and aggregate audit of already independently replayed arithmetic. Default execution makes no point calls or new arithmetic proof. Complete search CPU includes independent proofs inside each original arm. Optional --arithmetic reconstructs those proofs and cold banks without searching. Timing reproducibility means replayable protocol and receipts, not bit-identical runtimes.'})
    print(json.dumps({'status':'PASS_RETAINED_RECEIPTS_AND_BUNDLES','stages':rows},indent=2),flush=True)


def arithmetic():
    import verify_cancellation_scheduler as verifier
    import cancellation_scheduler_fresh as fresh
    for name,rawname,_ in STAGES[:-1]:
        verifier.OUT=OUT.parent/name;verifier.RAW=LOCAL/rawname;verifier.main()
    fresh.replay()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--arithmetic',action='store_true')
    if parser.parse_args().arithmetic:arithmetic()
    audit()
