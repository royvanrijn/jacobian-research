#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source=(p/'unmixed_depth_audit.sing').read_text().split('print("UNMIXED_DEPTH_START")')[0]
source+='\nideal LU=lead(U);\nring rp=32003,(a,b,c,d,e,f),dp;\nmap specialize=r,a,b,c,d,e,f;\nideal Up=specialize(U);ideal LUp=specialize(LU);\nprint("LEADING_SPECIALIZATION_CHECK");\nfor(int i=1;i<=size(Up);i++){\n if(leadmonom(Up[i])!=leadmonom(LUp[i])){print("BAD_LEADING_COEFFICIENT");quit;}\n}\nprint("LEADING_COEFFICIENTS_ALL_UNITS");print(size(Up));\nattrib(Up,"isSB",1);\nprint("MODULAR_RESOLUTION_START");\nresolution RU=mres(Up,0);\nprint("RESOLUTION_LENGTH");print(size(RU));\nprint("BETTI");print(betti(RU));\nprint("END");quit;\n'
fpath=p/'unmixed_modular_depth.sing';fpath.write_text(source)
log=p/'unmixed_modular_depth_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact good-prime resolution of the isolated-primary intersection, with leading-coefficient specialization checks; transfer to characteristic zero requires the written integral flatness and Koszul-rank argument.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'unmixed_modular_depth_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
