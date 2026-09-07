#!/usr/bin/env python3
"""Deterministic, bounded ECM stages with per-curve/factor checkpoints."""
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT/'artifacts/local/elliptic-curves/submitted627-630-conductors-v1'
ECM = '/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/ecm'
SCHEDULE = [2000]*25 + [11000]*50 + [50000]*75 + [250000]*75 + [1000000]*50


def prp(n):
    p = subprocess.run(['gp','-fq'],input=f'print(ispseudoprime({n}));\n',
                       text=True,capture_output=True,timeout=5)
    if p.returncode or p.stdout.strip() not in ('0','1'):
        raise ArithmeticError('probable-prime screen failed')
    return p.stdout.strip() == '1'


def run(row):
    identifier = row['icarm_id']
    out = WORK/f'ecm_{identifier}'
    out.mkdir(exist_ok=True)
    path = out/'state.json'
    prior = json.loads((WORK/f'pari_partial_{identifier}.json').read_text())
    if path.exists():
        state = json.loads(path.read_text())
    else:
        factors = [int(f['factor']) for f in prior['factors'] for _ in range(f['exponent'])]
        if not factors:
            factors = [int(row['remaining_cofactor'])]
        state = {'icarm_id':identifier,'factors':list(map(str,factors)),
                 'completed_attempts':0,'elapsed_seconds':0,'status':'RUNNING'}
    if state['status'] in ('PROBABLE_PRIME_FACTORIZATION','BUDGET_EXHAUSTED'):
        return
    factors = list(map(int,state['factors']))
    def checkpoint():
        if math.prod(factors) != int(row['remaining_cofactor']):
            raise ArithmeticError('ECM factor product mismatch')
        state['factors'] = list(map(str,sorted(factors)))
        tmp = path.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(path)
    for attempt in range(state['completed_attempts'],len(SCHEDULE)):
        composite = next((n for n in sorted(factors) if not prp(n)),None)
        if composite is None:
            state['status']='PROBABLE_PRIME_FACTORIZATION';checkpoint()
            print(identifier,'FACTORIZATION COMPLETE (primality proofs pending)',state['factors'],flush=True)
            return
        if state['elapsed_seconds'] >= 600:
            break
        b1 = SCHEDULE[attempt]
        sigma = 1000000 + identifier*1000 + attempt
        argv = [ECM,'-sigma',str(sigma),'-maxmem','256',str(b1)]
        start=time.monotonic()
        try:
            p=subprocess.run(argv,input=str(composite)+'\n',text=True,capture_output=True,
                             timeout=min(60,600-state['elapsed_seconds']))
            stdout,stderr=p.stdout,p.stderr
            status='RETURNED';returncode=p.returncode
        except subprocess.TimeoutExpired as e:
            stdout,stderr=e.stdout or b'',e.stderr or b''
            stdout=stdout.decode() if isinstance(stdout,bytes) else stdout
            stderr=stderr.decode() if isinstance(stderr,bytes) else stderr
            status='TIMEOUT';returncode=None
        elapsed=time.monotonic()-start;state['elapsed_seconds']+=elapsed
        (out/f'{attempt:04d}.stdout').write_text(stdout)
        (out/f'{attempt:04d}.stderr').write_text(stderr)
        record={'argv':argv,'input':str(composite),'wall_seconds':elapsed,'status':status,'returncode':returncode}
        (out/f'{attempt:04d}.json').write_text(json.dumps(record,indent=2)+'\n')
        for found in re.findall(r'Factor found in step \d+:\s*(\d+)',stdout):
            q=math.gcd(int(found),composite)
            if 1<q<composite:
                factors.remove(composite);factors.extend([q,composite//q])
                print(identifier,'SPLIT',len(str(composite)),'digits ->',len(str(q)),len(str(composite//q)),flush=True)
                break
        state['completed_attempts']=attempt+1;checkpoint()
        if (attempt+1)%25==0:
            print(identifier,'ECM',attempt+1,'attempts;',round(state['elapsed_seconds'],1),'seconds',flush=True)
    state['status']='BUDGET_EXHAUSTED';checkpoint()
    print(identifier,'ECM BUDGET EXHAUSTED',state['factors'],flush=True)


if __name__=='__main__':
    (WORK/'ecm_schedule.json').write_text(json.dumps({'B1':SCHEDULE,'sigma':'1000000+icarm_id*1000+attempt',
        'max_workers':2,'seconds_per_curve':600,'per_attempt_seconds':60},indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,json.loads((WORK/'protocol.json').read_text())['roster']))
