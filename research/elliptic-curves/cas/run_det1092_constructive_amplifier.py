#!/usr/bin/env python3
"""One source-bound, finite V3 attempt from the separately constructed M18."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import det1092_funnel as f
from v3_warm_support import atomic,bindings,read,require,sha
from research_runtime.supervisor import run,Limits


def worker(folder):
    plan=read(folder/'amplifier-launch.json')
    bindings(f.ROOT,plan['bindings'])
    atomic(folder/'amplifier-status.json',dict(status='RUNNING',pid=os.getpid()))
    command=[shutil.which('sage'),'-python',str(f.CAS/'det1092_funnel_worker.py'),
             'amplify','--directory',str(folder),'--case','conic-small-01']
    result=run(command,limits=Limits(plan['wall_seconds'],plan['rss_bytes']),
               log_path=folder/'amplifier.log',checkpoint_path=folder/'amplifier.supervisor.json',
               cwd=f.ROOT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','PYTHONUNBUFFERED':'1'})
    queue=folder/'amplifiers/conic-small-01/queue.json'
    passed=result['outcome']=='completed' and result['returncode']==0 and queue.exists()
    atomic(folder/'amplifier-status.json',dict(status='COMPLETE_REPLAYED' if passed else 'ERROR_OR_CENSORED',
            supervision=result,queue=read(queue) if passed else None))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['launch','worker']);parser.add_argument('--directory',type=Path,required=True)
    a=parser.parse_args();folder=a.directory.resolve()
    require(folder.is_relative_to(f.ROOT/'artifacts/local/elliptic-curves'),'local evidence directory required')
    if a.action=='worker':
        worker(folder);return
    require(not (folder/'amplifier-launch.json').exists(),'attempt already launched; preserve its limits and evidence')
    require(read(folder/'constructive-result.json')['status']=='CERTIFIED_SMALL_CONIC_M18','missing constructed seed')
    p=read(folder/'protocol.json');bindings(f.ROOT,p['sources']);bindings(f.ROOT,p['inputs'])
    paths=[Path(__file__).resolve(),f.CAS/'det1092_funnel_worker.py',folder/'protocol.json',
           folder/'seeds/conic-small-01/m18.json',folder/'constructive-result.json']
    plan=dict(case='conic-small-01',wall_seconds=14400,rss_bytes=3*1024**3,
              bindings={str(p.relative_to(f.ROOT)):sha(p) for p in paths},
              scope='One unchanged V3 cascade and independent replay from the smaller formula-derived M18. Separate constructive lane; no population refill.')
    atomic(folder/'amplifier-launch.json',plan,immutable=True)
    with (folder/'amplifier-controller.log').open('ab',buffering=0) as log:
        proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker','--directory',str(folder)],
                              cwd=f.ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    print('LAUNCHED_CONSTRUCTIVE_V3',proc.pid,flush=True)


if __name__=='__main__':
    main()
