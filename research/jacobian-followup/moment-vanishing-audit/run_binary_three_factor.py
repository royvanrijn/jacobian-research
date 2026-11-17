#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='ring r=0,(a,b,c,d,e,f),dp;\noption(redSB);\nideal I=a^2*b*e*f^2 + a^2*c*e^2*f + a*b^2*d*f^2 + 4*a*b*c*d*e*f + a*c^2*d*e^2 + b^2*c*d^2*f + b*c^2*d^2*e,-3*a^4*b^2*e^2*f^4 - 5*a^4*b*c*e^3*f^3 - 3*a^4*c^2*e^4*f^2 - 5*a^3*b^3*d*e*f^4 - 18*a^3*b^2*c*d*e^2*f^3 - 18*a^3*b*c^2*d*e^3*f^2 - 5*a^3*c^3*d*e^4*f - 3*a^2*b^4*d^2*f^4 - 18*a^2*b^3*c*d^2*e*f^3 - 33*a^2*b^2*c^2*d^2*e^2*f^2 - 18*a^2*b*c^3*d^2*e^3*f - 3*a^2*c^4*d^2*e^4 - 5*a*b^4*c*d^3*f^3 - 18*a*b^3*c^2*d^3*e*f^2 - 18*a*b^2*c^3*d^3*e^2*f - 5*a*b*c^4*d^3*e^3 - 3*b^4*c^2*d^4*f^2 - 5*b^3*c^3*d^4*e*f - 3*b^2*c^4*d^4*e^2,a^6*b^3*e^3*f^6 + 3*a^6*b^2*c*e^4*f^5 + 3*a^6*b*c^2*e^5*f^4 + a^6*c^3*e^6*f^3 + 3*a^5*b^4*d*e^2*f^6 + 15*a^5*b^3*c*d*e^3*f^5 + 24*a^5*b^2*c^2*d*e^4*f^4 + 15*a^5*b*c^3*d*e^5*f^3 + 3*a^5*c^4*d*e^6*f^2 + 3*a^4*b^5*d^2*e*f^6 + 24*a^4*b^4*c*d^2*e^2*f^5 + 60*a^4*b^3*c^2*d^2*e^3*f^4 + 60*a^4*b^2*c^3*d^2*e^4*f^3 + 24*a^4*b*c^4*d^2*e^5*f^2 + 3*a^4*c^5*d^2*e^6*f + a^3*b^6*d^3*f^6 + 15*a^3*b^5*c*d^3*e*f^5 + 60*a^3*b^4*c^2*d^3*e^2*f^4 + 93*a^3*b^3*c^3*d^3*e^3*f^3 + 60*a^3*b^2*c^4*d^3*e^4*f^2 + 15*a^3*b*c^5*d^3*e^5*f + a^3*c^6*d^3*e^6 + 3*a^2*b^6*c*d^4*f^5 + 24*a^2*b^5*c^2*d^4*e*f^4 + 60*a^2*b^4*c^3*d^4*e^2*f^3 + 60*a^2*b^3*c^4*d^4*e^3*f^2 + 24*a^2*b^2*c^5*d^4*e^4*f + 3*a^2*b*c^6*d^4*e^5 + 3*a*b^6*c^2*d^5*f^4 + 15*a*b^5*c^3*d^5*e*f^3 + 24*a*b^4*c^4*d^5*e^2*f^2 + 15*a*b^3*c^5*d^5*e^3*f + 3*a*b^2*c^6*d^5*e^4 + b^6*c^3*d^6*f^3 + 3*b^5*c^4*d^6*e*f^2 + 3*b^4*c^5*d^6*e^2*f + b^3*c^6*d^6*e^3;\nLIB "primdec.lib";\nideal G=std(I);print("STANDARD_BASIS");print(G);print("HILBERT_FIRST");print(hilb(G,1));\nprint("PRIMARY_START");list D=primdecGTZ(I);\nprint("PRIMARY_COUNT");print(size(D));\nideal Inter=1;\nfor(int j=1;j<=size(D);j++){print("COMPONENT");print(j);print("PRIMARY_IDEAL");print(std(D[j][1]));print("ASSOCIATED_PRIME");print(std(D[j][2]));Inter=intersect(Inter,D[j][1]);}\nideal GI=std(Inter);print("INTERSECTION_FORWARD");print(size(reduce(G,GI)));print("INTERSECTION_REVERSE");print(size(reduce(GI,G)));\nprint("END");quit;\n'
stem='binary_three_factor'
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
