#!/usr/bin/env python3
"""Detached conditional class1 closure -> exact next frame -> small panel.

No model calls; every realization stage has its own time/RSS cap and receipt.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

from run_class1_prospective_search import read,write,sha

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'


def prepare(folder, predecessor):
    from run_euclidean_seed_foundry import SUFFIXES
    folder.mkdir(parents=True,exist_ok=False)
    rt=folder/'runtime/research'
    for directory in (ROOT/'elliptic-curves/cas',ROOT/'elliptic-curves/ecsearch',ROOT/'elkies-k3/scripts'):
        for p in directory.rglob('*'):
            if p.is_file() and p.suffix in SUFFIXES and '__pycache__' not in p.parts:
                dest=rt/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in (ROOT/'elliptic-curves').glob('*.py'):
        dest=rt/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    names=['curve302_recovered_mw17_parent_v1.json','det1092_pruned_rootless_j2_census_v1.json',
        'det1092_frame_realization_priority_v1.json','det1092_pruned_anchor_packets_v1.zip',
        'curve302_parent_degree2_multisection_orbits_v1.tsv']
    for name in names:
        dest=rt/'artifacts/generated-results/elliptic-curves'/name
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ART/name,dest)
    ranked=read(ART/'det1092_frame_realization_priority_v1.json')['ranked_new_types']
    next_class=next(r['class_index'] for r in ranked if r['class_index']!=1)
    assert next_class==3
    plan={'schema':'x1092.frame-transition.v1','predecessor':str(predecessor),
        'predecessor_cap_sha256':sha(predecessor/'commissioning-cap.json'),
        'runtime_root':str(rt),'publish_root':str(ROOT),'target_class':next_class,
        'release_requires':'NO_EVIDENCE_CURRENT_SEARCH_PRODUCTIVE',
        'ranked_panel':32,'control_panel':8,'address_offset':65536,'address_count':65536,
        'score_prime_bound':997,'workers':2,'sage':shutil.which('sage'),
        'stages':[['discover',120],['marking',120],['trace',600],['equation',600],
                  ['section_plan',120],['sections',600],['normalize',120],['compact',120],['verify',180]],
        'stage_rss_bytes':3*1024**3,
        'failure':'Retain UNKNOWN and all checkpoints; do not launch a parameter panel after failed/incomplete realization.',
        'source_files':{str(p.relative_to(rt)):sha(p) for p in sorted(rt.rglob('*')) if p.is_file()},
        'scoring':'Same frozen scoring function, rational addresses, control ordering and point exposure as class1; only the generic parent changes.',
        'arithmetic_strict_gate':'Optional and UNKNOWN; not a prerequisite for ordinary search.',
        'next_panel_auto_extension':False}
    write(folder/'plan.json',plan,True)
    print('PREPARED',folder,flush=True)


def run(folder):
    plan=read(folder/'plan.json');rt=Path(plan['runtime_root'])
    frozen=rt/'elliptic-curves/cas'/Path(__file__).name
    if Path(__file__).resolve()!=frozen.resolve():
        os.execv(sys.executable,[sys.executable,str(frozen),'run','--folder',str(folder)])
    lock=(folder/'transition.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    for name,digest in plan['source_files'].items():assert sha(rt/name)==digest,name
    if (folder/'result.json').exists():return
    predecessor=Path(plan['predecessor'])
    assert sha(predecessor/'commissioning-cap.json')==plan['predecessor_cap_sha256']
    write(folder/'controller.json',{'pid':os.getpid(),'started_at':time.time()})
    while not (predecessor/'commissioning-result.json').exists():
        write(folder/'STATUS.json',{'status':'WAITING_FOR_CAPPED_CLASS1','updated_at':time.time()})
        if (folder/'STOP').exists():return
        time.sleep(10)
    prior=read(predecessor/'commissioning-result.json')
    if prior['status']!=plan['release_requires']:
        result={'status':'HOLD_PREDECESSOR_GAIN_OR_UNKNOWN','predecessor':prior,'next_frame_started':False}
        write(folder/'result.json',result,True);write(folder/'STATUS.json',result);return
    # Portable matched baseline: full frozen score distribution + exact outcome references.
    publish=Path(plan['publish_root'])/'artifacts/generated-results/elliptic-curves'
    baseline=publish/'class1_commissioning_64_16_baseline_v1'
    baseline.mkdir(exist_ok=True)
    for p in (predecessor/'matched-baseline').iterdir():
        dest=baseline/p.name
        if dest.exists():assert sha(dest)==sha(p)
        else:shutil.copyfile(p,dest)
    write(baseline/'commissioning-result.json',prior,True)
    from research_runtime.supervisor import run as supervised,Limits
    costs=[];prefix='x1092_class'+str(plan['target_class'])+'_realization'
    output=rt/'artifacts/generated-results/elliptic-curves'
    for stage,seconds in plan['stages']:
        job=folder/'stages'/stage;job.mkdir(parents=True,exist_ok=True)
        if (job/'cost.json').exists():
            cost=read(job/'cost.json')
        else:
            if (job/'started.json').exists():
                write(folder/'STATUS.json',{'status':'UNKNOWN_INTERRUPTED_STAGE','stage':stage});return
            write(job/'started.json',{'time':time.time()},True)
            script='verify_x1092_next_frame.sage' if stage=='verify' else 'realize_x1092_next_frame.sage'
            command=[plan['sage'],'-python',str(rt/'elkies-k3/scripts'/script)]
            if stage!='verify':command.append(stage)
            command+=['--class-index',str(plan['target_class'])]
            before=resource.getrusage(resource.RUSAGE_CHILDREN)
            write(folder/'STATUS.json',{'status':'REALIZING_FRAME','class_index':plan['target_class'],'stage':stage,'updated_at':time.time()})
            result=supervised(command,limits=Limits(seconds,plan['stage_rss_bytes']),
                log_path=job/'worker.log',checkpoint_path=job/'supervisor.json',cwd=rt,
                env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            cost={'stage':stage,'outcome':result['outcome'],'returncode':result['returncode'],
                'child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                'supervisor':result}
            write(job/'cost.json',cost,True)
        costs.append(cost)
        if cost['outcome']!='completed' or cost['returncode']!=0:
            result={'status':'UNKNOWN_BOUNDED_REALIZATION','stage':stage,'costs':costs,'panel_started':False}
            write(folder/'result.json',result,True);write(folder/'STATUS.json',result);return
        if stage in ('discover','section_plan'):
            name='discovery' if stage=='discover' else stage
            if not read(output/(prefix+'_'+name+'_v1.json'))['status'].startswith('PASS_'):
                result={'status':'UNKNOWN_BOUNDED_REALIZATION','stage':stage,'costs':costs,'panel_started':False}
                write(folder/'result.json',result,True);write(folder/'STATUS.json',result);return
    certificate={'status':'PASS_SEPARATE_EXACT_REPLAY','class_index':plan['target_class'],
        'generic_rank':17,'nonisometric_to_classes':[1,6],
        'files':{p.name:sha(p) for p in sorted(output.glob(prefix+'_*_v1.json'))},
        'producer_sha256':sha(rt/'elkies-k3/scripts/realize_x1092_next_frame.sage'),
        'checker_sha256':sha(rt/'elkies-k3/scripts/verify_x1092_next_frame.sage'),
        'reused_compiler_sha256':sha(rt/'elkies-k3/scripts/realize_x1092_class1.sage'),
        'stages':costs,'arithmetic_strict_classes':'UNKNOWN',
        'scope':'One exact rational J2 realization, no J1 classification, new K3, or specialized-rank claim.'}
    write(output/(prefix+'_manifest_v1.json'),certificate,True)
    for p in output.glob(prefix+'_*_v1.json'):
        dest=publish/p.name
        if dest.exists():assert sha(dest)==sha(p)
        else:shutil.copyfile(p,dest)
    parent=output/(prefix+'_compact_parent_v1.json')
    panel=folder/'commissioning-panel'
    command=[sys.executable,str(rt/'elliptic-curves/cas/run_class1_prospective_search.py')]
    if not (panel/'plan.json').exists():
        subprocess.run(command+['prepare','--folder',str(panel),'--parent',str(parent),'--parent-sha256',sha(parent),
            '--workers',str(plan['workers']),'--window',str(plan['address_count']),
            '--offset',str(plan['address_offset']),'--prime-bound',str(plan['score_prime_bound'])],check=True)
    panel_root=Path(read(panel/'plan.json')['root'])
    if not (panel_root/'ordinary-search/window-000/scores.json').exists():
        subprocess.run(command+['score','--folder',str(panel)],check=True)
    capped=[sys.executable,str(rt/'elliptic-curves/cas/run_frame_commissioning.py')]
    if not (panel/'commissioning-cap.json').exists():
        subprocess.run(capped+['prepare','--folder',str(panel),'--ranked',str(plan['ranked_panel']),
            '--controls',str(plan['control_panel'])],check=True)
    write(folder/'STATUS.json',{'status':'COMMISSIONING_NEXT_FRAME','class_index':plan['target_class'],'panel':str(panel),'updated_at':time.time()})
    subprocess.run(capped+['run','--folder',str(panel)],check=True)
    if not (panel/'commissioning-result.json').exists():return
    result={'status':'NEXT_FRAME_SMALL_PANEL_FINISHED','class_index':plan['target_class'],
        'panel':read(panel/'commissioning-result.json'),'realization_child_cpu_seconds':sum(c['child_cpu_seconds'] for c in costs),
        'no_automatic_panel_extension':True}
    write(folder/'result.json',result,True);write(folder/'STATUS.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','run'])
    p.add_argument('--folder',type=Path,required=True);p.add_argument('--predecessor',type=Path)
    a=p.parse_args();a.folder=a.folder.resolve()
    if a.mode=='prepare':prepare(a.folder,a.predecessor.resolve())
    else:run(a.folder)
