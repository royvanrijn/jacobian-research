#!/usr/bin/env python3
"""Sequential, checkpointed complete-CPU comparison; no silent arm reruns."""
import fcntl
from itertools import permutations
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

from cancellation_scheduler_prepare import OUT, RAW, CAS
from finite_cancellation_corpus import ROOT, digest, write


def cpu():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime+r.ru_stime


def main():
    RAW.mkdir(parents=True,exist_ok=True)
    with (RAW/'supervisor.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        run()


def run():
    plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
    for name,h in plan['source_sha256'].items():
        if digest((ROOT/name).read_bytes())!=h:raise ArithmeticError(f'source changed: {name}')
    inputs=json.loads((OUT/'inputs.json').read_text())
    if digest((OUT/'inputs.json').read_bytes())!=plan['input_sha256']['inputs.json']:
        raise ArithmeticError('inputs changed')
    orders=list(permutations(plan['arms']));records=[]
    for i,case in enumerate(inputs):
        for arm in orders[i%len(orders)]:
            dest=RAW/'arms'/case['id']/arm;dest.mkdir(parents=True,exist_ok=True)
            if (dest/'supervisor.json').exists():
                receipt=json.loads((dest/'supervisor.json').read_text())
                if receipt['protocol_sha256']!=ph:raise ArithmeticError('resume protocol differs')
                records.append(receipt);continue
            if (dest/'start.json').exists() or (dest/'result.json').exists():
                raise ArithmeticError('Interrupted/unreceipted arm retained. Reconcile it; do not rerun or delete it.')
            command=['sage','-python',str(CAS/'cancellation_scheduler_cpu.py'),'--case',case['id'],'--arm',arm]
            write(dest/'start.json',{'protocol_sha256':ph,'command':command})
            before=cpu();wall=time.monotonic()
            with (dest/'worker.log').open('w') as log:
                process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:
                    rc=process.wait(timeout=plan['arm_wall_seconds'])
                    status='COMPLETE' if rc==0 else 'WORKER_FAILURE'
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=2)
                    except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                    rc=process.returncode;status='OUTER_WALL_TIMEOUT'
            receipt={'case':case['id'],'family':case['family'],'stratum':case['stratum'],'arm':arm,
                'status':status,'returncode':rc,'charged_cpu_seconds':cpu()-before,
                'wall_seconds':time.monotonic()-wall,'protocol_sha256':ph,
                'includes_independent_rank_certification':status=='COMPLETE'}
            result={}
            if (dest/'result.json').exists():
                receipt['result_sha256']=digest((dest/'result.json').read_bytes())
                result=json.loads((dest/'result.json').read_text())
            write(dest/'supervisor.json',receipt);records.append(receipt)
            write(OUT/'supervision.json',{'status':'RUNNING','protocol_sha256':ph,'records':records})
            print(json.dumps({'arms_complete':len(records),'total':len(inputs)*len(plan['arms']),
                'case':case['id'],'arm':arm,'status':status,'gain':result.get('success'),
                'calls':len(result.get('calls',[])),'cpu':receipt['charged_cpu_seconds']}),flush=True)
            if status!='COMPLETE':
                write(OUT/'supervision.json',{'status':'STOPPED_INFRASTRUCTURE_UNKNOWN','protocol_sha256':ph,'records':records})
                return
    write(OUT/'supervision.json',{'status':'COMPLETE','protocol_sha256':ph,'records':records})


if __name__=='__main__':main()
