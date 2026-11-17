#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
text=(p/'three_factor_minimal_primary_stdout.txt').read_text()
basis=text.split('UNMIXED_BASIS\n')[1].split('FULL_IN_UNMIXED')[0].strip()
source='ring r=0,(a,b,c,d,e,f),dp;\noption(redSB);\nideal U='+basis+';\nattrib(U,"isSB",1);\n'
source+=r'''print("UNMIXED_DEPTH_START");
resolution RU=mres(U,0);
print("RESOLUTION_LENGTH");print(size(RU));
print("BETTI");print(betti(RU));
print("END");quit;
'''
fpath=p/'unmixed_depth_audit.sing';fpath.write_text(source)
log=p/'unmixed_depth_audit_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Minimal free resolution/depth of the retained exact isolated-primary intersection U, solely to test whether the quaternary generic-section torsion lower bound is exact.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'unmixed_depth_audit_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
