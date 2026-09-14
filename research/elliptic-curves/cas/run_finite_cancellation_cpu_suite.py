#!/usr/bin/env python3
"""Sequential outer supervisor; charge interpreter and descendant CPU as well."""
import json
from pathlib import Path
import resource
import subprocess
import time

from finite_cancellation_corpus import OUT, digest, write

def cpu():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime+r.ru_stime

def main():
    folder=OUT/'cpu';plan=json.loads((folder/'protocol.json').read_text());inputs=json.loads((folder/'inputs.json').read_text())
    assert digest((folder/'inputs.json').read_bytes())==plan['inputs_sha256']
    records=[]
    for i,case in enumerate(inputs):
        arms=plan['arms'] if i%2==0 else plan['arms'][::-1]
        for arm in arms:
            dest=folder/'arms'/case['id']/arm
            if (dest/'supervisor.json').exists():
                record=json.loads((dest/'supervisor.json').read_text());assert record['protocol_sha256']==digest((folder/'protocol.json').read_bytes())
                records.append(record);continue
            assert not (dest/'result.json').exists(),'unreceipted earlier timing must remain a separate pilot'
            dest.mkdir(parents=True,exist_ok=True);start=cpu();wall=time.monotonic()
            command=['sage','-python',str(Path(__file__).with_name('finite_cancellation_cpu.py')),'run','--case',case['id'],'--arm',arm]
            with (dest/'worker.log').open('w') as log:
                try:
                    p=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=90)
                    status='COMPLETE' if p.returncode==0 else 'WORKER_FAILURE';returncode=p.returncode
                except subprocess.TimeoutExpired:status='OUTER_WALL_TIMEOUT';returncode=None
            record={'case':case['id'],'family':case['family'],'arm':arm,'status':status,'returncode':returncode,
                'charged_cpu_seconds':cpu()-start,'wall_seconds':time.monotonic()-wall,'command':command,
                'protocol_sha256':digest((folder/'protocol.json').read_bytes()),'supervisor_source_sha256':digest(Path(__file__).read_bytes())}
            if (dest/'result.json').exists():record['result_sha256']=digest((dest/'result.json').read_bytes())
            write(dest/'supervisor.json',record);records.append(record)
            write(folder/'supervision.json',{'status':'RUNNING','records':records})
            print(json.dumps({k:v for k,v in record.items() if k not in ['command','supervisor_source_sha256']}),flush=True)
    write(folder/'supervision.json',{'status':'COMPLETE','records':records})

if __name__=='__main__':main()
