#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source=(p/'unmixed_depth_audit.sing').read_text().split('print("UNMIXED_DEPTH_START")')[0]
source+='ideal L=std(lead(U));print("MONOMIAL_DEPTH_START");\nresolution RU=mres(L,0);\nprint("RESOLUTION_LENGTH");print(size(RU));\nprint("BETTI");print(betti(RU));\nprint("END");quit;\n'
fpath=p/'unmixed_monomial_depth.sing';fpath.write_text(source)
log=p/'unmixed_monomial_depth_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Minimal resolution of the leading monomial ideal of the retained isolated-primary intersection; only a proved depth bound may transfer by Groebner degeneration.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'unmixed_monomial_depth_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
