#!/usr/bin/env python3
"""Sequential complete-cost supervisor for the frozen three-policy validation."""
from itertools import permutations
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

from finite_cancellation_corpus import LOCAL,digest,write
from finite_cancellation_validation_audit import OUT

RAW=LOCAL/'finite-cancellation-validation-v2'


def cpu():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime+r.ru_stime


def main():
    plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
    inputs=json.loads((OUT/'inputs.json').read_text())
    assert digest((OUT/'inputs.json').read_bytes())==plan['input_sha256']['inputs.json']
    assert digest(Path(__file__).read_bytes())==plan['source_sha256'][Path(__file__).name]
    orders=list(permutations(plan['arms']));records=[]
    for i,case in enumerate(inputs):
        for arm in orders[i%len(orders)]:
            dest=RAW/'arms'/case['id']/arm;dest.mkdir(parents=True,exist_ok=True)
            if (dest/'supervisor.json').exists():
                old=json.loads((dest/'supervisor.json').read_text());assert old['protocol_sha256']==ph
                records.append(old);continue
            assert not (dest/'result.json').exists(),'Unreceipted earlier results must remain separate.'
            before=cpu();wall=time.monotonic()
            command=['sage','-python',str(Path(__file__).with_name('finite_cancellation_validation_cpu.py')),'--case',case['id'],'--arm',arm]
            with (dest/'worker.log').open('w') as log:
                p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:
                    rc=p.wait(timeout=plan['arm_wall_seconds']);status='COMPLETE' if rc==0 else 'WORKER_FAILURE'
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid,signal.SIGTERM)
                    try:p.wait(timeout=2)
                    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                    rc=p.returncode;status='OUTER_WALL_TIMEOUT'
            record={'case':case['id'],'family':case['family'],'stratum':case['stratum'],'arm':arm,
                'status':status,'returncode':rc,'charged_cpu_seconds':cpu()-before,'wall_seconds':time.monotonic()-wall,
                'command':command,'protocol_sha256':ph,'supervisor_source_sha256':digest(Path(__file__).read_bytes())}
            if (dest/'result.json').exists():record['result_sha256']=digest((dest/'result.json').read_bytes())
            write(dest/'supervisor.json',record);records.append(record)
            write(OUT/'supervision.json',{'status':'RUNNING','protocol_sha256':ph,'records':records})
            result=json.loads((dest/'result.json').read_text()) if (dest/'result.json').exists() else {}
            print(json.dumps({'completed_arms':len(records),'total_arms':3*len(inputs),
                'case':case['id'],'stratum':case['stratum'],'arm':arm,'status':status,
                'success':result.get('success'),'calls':result.get('point_search_calls'),
                'cpu_seconds':record['charged_cpu_seconds']}),flush=True)
    write(OUT/'supervision.json',{'status':'COMPLETE','protocol_sha256':ph,'records':records})


if __name__=='__main__':main()
