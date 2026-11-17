#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='ring r=0,(x,y,z,w),dp;\noption(redSB);\nideal I=w*x + y*z,w*y^3 + x*z^3 - 3*y^2*z^2,y^3*z^3;\nLIB "primdec.lib";\nideal G=std(I);print("STANDARD_BASIS");print(G);print("HILBERT_FIRST");print(hilb(G,1));\nprint("PRIMARY_START");list D=primdecGTZ(I);\nprint("PRIMARY_COUNT");print(size(D));\nideal Inter=1;\nfor(int j=1;j<=size(D);j++){print("COMPONENT");print(j);print("PRIMARY_IDEAL");print(std(D[j][1]));print("ASSOCIATED_PRIME");print(std(D[j][2]));Inter=intersect(Inter,D[j][1]);}\nideal GI=std(Inter);print("INTERSECTION_FORWARD");print(size(reduce(G,GI)));print("INTERSECTION_REVERSE");print(size(reduce(GI,G)));\nprint("END");quit;\n'
stem='binary_cubic_coefficients'
fpath=p/(stem+'.sing');fpath.write_text(source)
log=p/(stem+'_stdout.txt')
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact primary decomposition of the binary coefficient moment ideal; full intersection must pass before promotion.','wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and result.rstrip().endswith('END') and '?' not in result,'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3500:]}
(p/(stem+'_run.json')).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
