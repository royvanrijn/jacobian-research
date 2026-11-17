#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='LIB "elim.lib";\nring rLQ_NONISO_D1ZERO_H1=0,(d2,d3),dp;\noption(redSB);\nideal All=-7*d2^2 - 20*d2 - 40*d3,6*d2^3 + 63*d2^2 + 126*d2*d3 + 63*d3^2 - 28*d3,11*d2^4 + 16*d2^3 + 32*d2^2*d3 + 176*d2^2 - 88*d2*d3^2 + 1024*d2*d3 + 1056*d3^2,-130*d2^5 - 2585*d2^4 - 5170*d2^3*d3 - 5720*d2^3 - 2145*d2^2*d3^2 - 25080*d2^2*d3 - 34320*d2*d3^2 + 2288*d2*d3 - 14300*d3^3 + 5720*d3^2;\nfor(int kk=2;kk<=4;kk++){ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\nprint("LQ_NONISO_D1ZERO_H1_GATE");print(2*(kk+1));ideal G=std(I);print("LQ_NONISO_D1ZERO_H1_BASIS_SIZE");print(size(G));print("LQ_NONISO_D1ZERO_H1_BASIS");print(G);\nfor(int pp=1;pp<=4;pp++){poly remainderPoly=reduce(d3^pp,G);print("LQ_NONISO_D1ZERO_H1_POWER");print(pp);print(size(remainderPoly));if(remainderPoly==0){break;}}\nideal Sat=sat(G,ideal(d3));print("LQ_NONISO_D1ZERO_H1_SATURATED");print(Sat);if(reduce(1,Sat)==0){break;}\n}\nring rLQ_NONISO_D1H1ZERO_D2=0,(d3),dp;\noption(redSB);\nideal All=-1,21*d3^2 + 2,1 - 8*d3^2,-33*d3^2 - 2;\nfor(int kk=2;kk<=4;kk++){ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\nprint("LQ_NONISO_D1H1ZERO_D2_GATE");print(2*(kk+1));ideal G=std(I);print("LQ_NONISO_D1H1ZERO_D2_BASIS_SIZE");print(size(G));print("LQ_NONISO_D1H1ZERO_D2_BASIS");print(G);\nfor(int pp=1;pp<=4;pp++){poly remainderPoly=reduce(d3^pp,G);print("LQ_NONISO_D1H1ZERO_D2_POWER");print(pp);print(size(remainderPoly));if(remainderPoly==0){break;}}\nideal Sat=sat(G,ideal(d3));print("LQ_NONISO_D1H1ZERO_D2_SATURATED");print(Sat);if(reduce(1,Sat)==0){break;}\n}\nring rLQ_ISO_H1=0,(d1,d2,d3),dp;\noption(redSB);\nideal I=-7*d1^2 - 20*d2,3*d1^3 - 14*d3,11*d1^4 + 104*d1^2*d2 - 32*d1*d3 + 176*d2^2;\nideal G=std(I);print("LQ_ISO_H1_BASIS");print(G);print("LQ_ISO_H1_D3_POWER");print(reduce(d3^4,G));\nprint("END");quit;\n'
fpath=p/'linear_quadratic_other_invariant_charts.sing';fpath.write_text(source)
log=p/'linear_quadratic_other_invariant_charts_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(30,30))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=30)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact invariant D/H coefficient charts for a linear form times a nondegenerate quadratic in three Gaussian variables; finite exclusions require the written Gaussian generating-function identity and exhaustive rescaling cover.',
'wall_limit_seconds':30,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'linear_quadratic_other_invariant_charts_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
