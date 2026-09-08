#!/usr/bin/env python3
"""Finish the two leading unpublished rank25 conductor candidates with GNFS."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT/'artifacts/local/elliptic-curves/conductor-record-screen-v1'
SOURCE = ROOT/'artifacts/local/elliptic-curves/cado-nfs-source'
IDS = ['new-20260906-54','new-20260906-43']


def run(item):
    index, identifier = item
    directory = WORK/'factors'/identifier
    state_path = directory/'state.json'
    frozen = directory/'state_ecm.json'
    if not frozen.exists():
        shutil.copyfile(state_path,frozen)
    state = json.loads(frozen.read_text())
    composite = max(map(int,state['factors']))
    directory = directory/'nfs'
    directory.mkdir(exist_ok=True)
    cpus = sorted(os.sched_getaffinity(0))[index*8:(index+1)*8]
    argv = ['taskset','-c',','.join(map(str,cpus)),
        str(SOURCE/'cado-nfs.venv/bin/python3'),str(SOURCE/'build/local/cado-nfs.py'),
        str(composite),'tasks.workdir='+str(directory/'work'),'tasks.threads=8',
        'slaves.nrclients=4','tasks.sieve.las.threads=2','server.address=127.0.0.1',
        'server.port='+str(8654+index)]
    protocol = {'id':identifier,'input':str(composite),'argv':argv,'wall_seconds':1800,
        'cpus':cpus,'source_commit':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),
        'previous_state_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'motivation':'These are the two smallest unpublished rank25 conductor bounds, '
            'and the only two of13 priorities unresolved after the fixed ECM schedule. '
            'Resolve repeated factors before excluding a small exact conductor.',
        'failure':'Retain all relations/checkpoints. Incomplete factors remain UNKNOWN.'}
    path = directory/'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        path.write_text(json.dumps(protocol,indent=2)+'\n')
    start = time.monotonic()
    with (directory/'cado.stdout').open('w') as out,(directory/'cado.stderr').open('w') as err:
        p = subprocess.run(['timeout','--signal=INT','--kill-after=30','1800']+argv,stdout=out,stderr=err)
    (directory/'result.json').write_text(json.dumps({'returncode':p.returncode,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
    if p.returncode:
        print(identifier,'NFS INCOMPLETE',p.returncode,flush=True)
        return
    factors = list(map(int,(directory/'cado.stdout').read_text().split()))
    assert math.prod(factors) == composite
    prior = list(map(int,state['factors']))
    prior.remove(composite)
    prior.extend(factors)
    assert math.prod(prior) == math.prod(map(int,state['factors']))
    state.update(factors=list(map(str,sorted(prior))),status='NFS_FACTORS_RETURNED_PRIMALITY_PENDING')
    state_path.write_text(json.dumps(state,indent=2)+'\n')
    print(identifier,'NFS FACTORS',[len(str(p)) for p in factors],flush=True)


if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,enumerate(IDS)))
