#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='ring r=0,(x,y,z,w),dp;\noption(redSB);\nideal I=w*x + y*z,w*y^3 + x*z^3 - 3*y^2*z^2,y^3*z^3;\nLIB "elim.lib";\nideal G=std(I);\nideal Q1=w,\nz;\nideal Q2=y,\nx;\nideal Q3=yz+xw,\nw3,\nz4-y2w2-3xzw2,\nxz3+y3w-3x2w2,\nz3w2,\ny4w-x2z2w-3x2yw2,\ny9;\nprint("AXIS_X_RADICAL_FORWARD");print(size(reduce(Q3,std(ideal(y,z,w)))));\nprint("AXIS_X_RADICAL_REVERSE");print(size(reduce(ideal(y12,z12,w12),std(Q3))));\nprint("AXIS_X_REGULAR");print(size(reduce(quotient(Q3,ideal(x)),std(Q3))));\nideal Q4=yz+xw,\nx3,\nxz3+y3w-3x2w2,\ny4-x2z2-3x2yw,\nx2y3,\nz9;\nprint("AXIS_W_RADICAL_FORWARD");print(size(reduce(Q4,std(ideal(x,y,z)))));\nprint("AXIS_W_RADICAL_REVERSE");print(size(reduce(ideal(x12,y12,z12),std(Q4))));\nprint("AXIS_W_REGULAR");print(size(reduce(quotient(Q4,ideal(w)),std(Q4))));\nideal Q5=yz+xw,\nw4,\nxz3+y3w-3x2w2,\nx4,\ny4w-x2z2w-3x2yw2,\nx3w3,\nx3y3w,\nz9,\ny9;\nprint("ORIGIN_RADICAL_REVERSE");print(size(reduce(ideal(x12,y12,z12,w12),std(Q5))));\nfor(int omit=1;omit<=5;omit++){\n ideal Rest=1;\n if(omit!=1){Rest=intersect(Rest,Q1);}if(omit!=2){Rest=intersect(Rest,Q2);}if(omit!=3){Rest=intersect(Rest,Q3);}if(omit!=4){Rest=intersect(Rest,Q4);}if(omit!=5){Rest=intersect(Rest,Q5);}\n print("IRREDUNDANT_COUNT");print(omit);print(size(reduce(std(Rest),G)));\n}\nideal OriginSat=std(intersect(intersect(sat(G,ideal(x)),sat(G,ideal(y))),intersect(sat(G,ideal(z)),sat(G,ideal(w)))));\nprint("FULL_HILBERT");print(hilb(G,1));print("ORIGIN_SAT_HILBERT");print(hilb(OriginSat,1));\nideal J=x*z,x*w,y*z,y*w;ideal Pow=J;\nfor(int nn=1;nn<=12;nn++){\n print("RADICAL_POWER");print(nn);print(size(Pow));print(size(reduce(Pow,G)));\n if(size(reduce(Pow,G))==0){break;}\n Pow=std(Pow*J);\n}\nprint("END");quit;\n'
stem='binary_cubic_primary_audit'
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
rec={'scope':'Validate the five primary ideals by radical-power and parameter-regularity checks, irredundancy, intrinsic origin torsion and complete nilradical-power tests.','wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and result.rstrip().endswith('END') and '?' not in result,'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3500:]}
(p/(stem+'_run.json')).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
