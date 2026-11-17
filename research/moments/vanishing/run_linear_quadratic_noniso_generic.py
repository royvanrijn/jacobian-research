#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='LIB "elim.lib";\nring rLQ_NONISO_D1=0,(d2,d3,h1),dp;\noption(redSB);\nideal All=-112*d2^2 - 320*d2*h1^2 + 896*d2*h1 - 360*d2 - 640*d3*h1 + 1120*d3 - 112*h1^2 + 480*h1 - 495,96*d2^3 + 1008*d2^2*h1^2 - 1728*d2^2*h1 + 144*d2^2 + 2016*d2*d3*h1 - 2016*d2*d3 + 792*d2*h1^2 - 2736*d2*h1 + 2466*d2 + 1008*d3^2 - 448*d3*h1^3 + 2016*d3*h1^2 - 2520*d3*h1 + 672*d3 + 96*h1^3 - 549*h1^2 + 1080*h1 - 723,2816*d2^4 + 4096*d2^3*h1^2 - 90112*d2^3*h1 + 81152*d2^3 + 8192*d2^2*d3*h1 - 101376*d2^2*d3 + 45056*d2^2*h1^4 - 262144*d2^2*h1^3 + 451584*d2^2*h1^2 - 198656*d2^2*h1 - 20960*d2^2 - 22528*d2*d3^2 + 262144*d2*d3*h1^3 - 1155072*d2*d3*h1^2 + 1511424*d2*d3*h1 - 478720*d2*d3 + 26624*d2*h1^4 - 212992*d2*h1^3 + 590592*d2*h1^2 - 657920*d2*h1 + 235728*d2 + 270336*d3^2*h1^2 - 901120*d3^2*h1 + 760320*d3^2 - 8192*d3*h1^4 + 131072*d3*h1^3 - 583680*d3*h1^2 + 1026560*d3*h1 - 634304*d3 + 2816*h1^4 - 28672*h1^3 + 102592*h1^2 - 155520*h1 + 85227,-33280*d2^5 - 661760*d2^4*h1^2 + 1664000*d2^4*h1 - 629760*d2^4 - 1323520*d2^3*d3*h1 + 1830400*d2^3*d3 - 1464320*d2^3*h1^4 + 6983680*d2^3*h1^3 - 10343680*d2^3*h1^2 + 5351680*d2^3*h1 - 1410560*d2^3 - 549120*d2^2*d3^2 - 6420480*d2^2*d3*h1^3 + 23485440*d2^2*d3*h1^2 - 23591040*d2^2*d3*h1 + 4942080*d2^2*d3 - 1436160*d2^2*h1^4 + 9297920*d2^2*h1^3 - 20986080*d2^2*h1^2 + 18694080*d2^2*h1 - 4675760*d2^2 - 8785920*d2*d3^2*h1^2 + 24161280*d2*d3^2*h1 - 14551680*d2*d3^2 + 585728*d2*d3*h1^5 - 4055040*d2*d3*h1^4 + 6589440*d2*d3*h1^3 + 7912960*d2*d3*h1^2 - 28168800*d2*d3*h1 + 18271968*d2*d3 - 112640*d2*h1^5 + 713600*d2*h1^4 - 1125120*d2*h1^3 - 1323280*d2*h1^2 + 5051440*d2*h1 - 3599550*d2 - 3660800*d3^3*h1 + 6223360*d3^3 + 1464320*d3^2*h1^4 - 8785920*d3^2*h1^3 + 17571840*d3^2*h1^2 - 11897600*d3^2*h1 + 377520*d3^2 + 168960*d3*h1^5 - 1745920*d3*h1^4 + 6920320*d3*h1^3 - 13273920*d3*h1^2 + 12423400*d3*h1 - 4576000*d3 - 33280*h1^5 + 354080*h1^4 - 1502080*h1^3 + 3181375*h1^2 - 3363600*h1 + 1418895;\nfor(int kk=2;kk<=4;kk++){ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\nprint("LQ_NONISO_D1_GATE");print(2*(kk+1));ideal G=std(I);print("LQ_NONISO_D1_BASIS_SIZE");print(size(G));print("LQ_NONISO_D1_BASIS");print(G);\nfor(int pp=1;pp<=4;pp++){poly remainderPoly=reduce(d3^pp,G);print("LQ_NONISO_D1_POWER");print(pp);print(size(remainderPoly));if(remainderPoly==0){break;}}\nideal Sat=sat(G,ideal(d3));print("LQ_NONISO_D1_SATURATED");print(Sat);if(reduce(1,Sat)==0){break;}\n}\nprint("END");quit;\n'
fpath=p/'linear_quadratic_noniso_generic.sing';fpath.write_text(source)
log=p/'linear_quadratic_noniso_generic_stdout.txt'
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
(p/'linear_quadratic_noniso_generic_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
