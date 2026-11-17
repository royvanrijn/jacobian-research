#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source=(p/'unmixed_modular_depth.sing').read_text()
source=source.replace('print("MODULAR_RESOLUTION_START");','ideal V=minbase(Up);print("MINIMAL_GENERATORS");print(size(V));print("MINIMIZED_RESOLUTION_START");').replace('mres(Up,0)','mres(V,0)')
fpath=p/'unmixed_modular_minimized_depth.sing';fpath.write_text(source)
log=p/'unmixed_modular_minimized_depth_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Good-prime resolution after exact minimal-generator pruning, preserving the leading-coefficient specialization check; a depth bound requires the written lifting argument.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'unmixed_modular_minimized_depth_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
