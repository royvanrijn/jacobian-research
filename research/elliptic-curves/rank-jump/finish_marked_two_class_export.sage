#!/usr/bin/env sage-python
"""Repair Sage-integer JSON export within the ORIGINAL candidate CPU budgets.

The original arithmetic implementation and all failed receipts are preserved.
Only export failures can continue, with their CPU subtracted from 30 seconds.
The interrupted child with no usage receipt is charged its entire allowance.
"""
import hashlib
import json
import os
from pathlib import Path
import resource
import runpy
import signal
import time
import traceback
from collections import Counter

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/local/elliptic-curves/marked-two-class-v1'


def encode(value):
    from sage.rings.integer import Integer
    if isinstance(value,Integer):return int(value)
    raise TypeError(f'Unsupported certificate value: {type(value)}')


def write(path,data):
    with path.with_suffix(path.suffix+'.tmp').open('w') as stream:
        json.dump(data,stream,indent=2,sort_keys=True,default=encode);stream.write('\n')
    path.with_suffix(path.suffix+'.tmp').replace(path)


def main():
    start=time.process_time()
    oldsource=OUT/'source-snapshots/elliptic-curves/rank-jump/compile_marked_two_class.sage'
    binding=json.loads((OUT/'compilation-protocol.json').read_text())
    assert hashlib.sha256(oldsource.read_bytes()).hexdigest()==next(iter(binding['sources'].values()))
    code=runpy.run_path(str(oldsource))
    # The frozen source computes ROOT from its path. Rebind its filesystem
    # constants only; setup and compile_one arithmetic remain byte identical.
    for fn in [code['setup'],code['compile_one']]:
        fn.__globals__['ROOT']=ROOT;fn.__globals__['OUT']=OUT
    old=json.loads((OUT/'compilation-progress.json').read_text())['receipts']
    assert len(old)==171 and all(r['status']=='COMPILATION_ERROR' for r in old)
    charges={r['index']:r['cpu_seconds'] for r in old}
    for r in old:
        packet=json.loads((OUT/'compiled'/f"{r['index']:03d}.json").read_text())
        assert packet['exception']=="TypeError('Object of type Integer is not JSON serializable')"
    interrupted=171;charges[interrupted]=30.0
    directory=OUT/'compiled-export-v2';directory.mkdir(exist_ok=False)
    write(OUT/'export-repair-protocol.json',{'status':'FROZEN_EXPORT_REPAIR','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'original_arithmetic_sha256':hashlib.sha256(oldsource.read_bytes()).hexdigest(),
        'fixed_candidate_sha256':hashlib.sha256((OUT/'frozen-candidates.json').read_bytes()).hexdigest(),
        'charged_prior_cpu_seconds':charges,'interrupted_child_charged_full_budget':interrupted,
        'change':'JSON serialization of Sage Integer only; exact arithmetic reused. Export failures continue within the original total30 CPU seconds per candidate. The interrupted child is not repeated.',
        'maximum_workers':1,'new_candidate_count':0})
    snapshot=OUT/'source-snapshots'/Path(__file__).relative_to(ROOT);snapshot.write_bytes(Path(__file__).read_bytes())
    ctx=code['setup']();receipts=[]
    candidates=json.loads((OUT/'frozen-candidates.json').read_text())['rows']
    for index,row in enumerate(candidates):
        path=directory/f'{index:03d}.json';prior=charges.get(index,0.0)
        remaining=30-prior
        if remaining<=0:
            write(path,{'index':index,'divisor':row,'status':'UNKNOWN_INTERRUPTED_BUDGET_CHARGED'})
            receipts.append({'index':index,'status':'UNKNOWN_INTERRUPTED_BUDGET_CHARGED','prior_charged_cpu_seconds':prior,'new_cpu_seconds':0,'total_charged_cpu_seconds':prior})
            continue
        wall=time.monotonic();pid=os.fork()
        if pid==0:
            signal.signal(signal.SIGPROF,signal.SIG_DFL)
            signal.setitimer(signal.ITIMER_PROF,max(0.001,remaining-0.05))
            resource.setrlimit(resource.RLIMIT_CPU,(30,31))
            try:
                write(path,code['compile_one'](index,row,ctx));os._exit(0)
            except BaseException as exc:
                write(path,{'index':index,'divisor':row,'status':'COMPILATION_ERROR','exception':repr(exc),'traceback':traceback.format_exc()});os._exit(1)
        _,status,usage=os.wait4(pid,0)
        if not path.exists():write(path,{'index':index,'divisor':row,'status':'UNKNOWN_CPU_LIMIT','wait_status':status})
        result=json.loads(path.read_text());cpu=usage.ru_utime+usage.ru_stime
        receipt={'index':index,'status':result['status'],'prior_charged_cpu_seconds':prior,'new_cpu_seconds':cpu,
            'total_charged_cpu_seconds':prior+cpu,'wall_seconds':time.monotonic()-wall,'max_rss_kib':usage.ru_maxrss,
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        receipts.append(receipt);write(OUT/'export-repair-progress.json',{'receipts':receipts})
        print(json.dumps(receipt),flush=True)
        if result['status']=='COMPILATION_ERROR':
            print('Stopped on an implementation error; unattempted candidates remain UNKNOWN.',flush=True);break
    write(OUT/'compilation-v2.json',{'status':'FROZEN_PREFIX_ATTEMPTED' if len(receipts)==256 else 'IMPLEMENTATION_STOP',
        'receipts':receipts,'histogram':dict(Counter(r['status'] for r in receipts)),
        'total_charged_child_cpu_seconds':sum(r['total_charged_cpu_seconds'] for r in receipts),
        'parent_cpu_seconds':time.process_time()-start,'shells_complete':False,
        'boundary':'Export-repair arithmetic uses only the frozen candidate prefix, with prior failed work charged. No control splitting or class labels influenced the bank.'})


if __name__=='__main__':main()
