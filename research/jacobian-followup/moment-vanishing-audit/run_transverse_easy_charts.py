#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='LIB "elim.lib";\nring rTRANSVERSE_A0D1E1=0,(b,f),dp;\noption(redSB);\nideal All=-5*b^2 - 20*b*f - 12*f^2 + 20,-f^3 - 7*f,11*b^4 + 88*b^3*f + 224*b^2*f^2 - 88*b^2 + 192*b*f^3 - 352*b*f + 48*f^4 - 144*f^2 + 176,55*b^2*f^3 + 429*b^2*f + 220*b*f^4 + 1716*b*f^2 + 120*f^5 + 704*f^3 - 1716*f;\npoly delta=f;\nfor(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n ideal G=std(I);print("TRANSVERSE_A0D1E1_GATE");print(2*(kk+1));\n print("TRANSVERSE_A0D1E1_BASIS");print(G);\n poly remainderPoly=reduce(delta^12,G);print("TRANSVERSE_A0D1E1_DET_REMAINDER_SIZE");print(size(remainderPoly));\n ideal Sat=sat(G,ideal(delta));print("TRANSVERSE_A0D1E1_SATURATED");print(Sat);\n if(reduce(1,Sat)==0){break;}\n}\nring rTRANSVERSE_DECOUPLED=0,(f),dp;\noption(redSB);\nideal All=-7*f^4 - 38*f^3 - 108*f^2 - 90*f + 195,4*f^6 + 38*f^5 + 199*f^4 + 552*f^3 + 918*f^2 + 810*f + 1319,11*f^8 + 92*f^7 + 320*f^6 + 164*f^5 - 3130*f^4 - 13372*f^3 - 27864*f^2 - 21060*f + 34887,-260*f^10 - 3634*f^9 - 26613*f^8 - 120210*f^7 - 356970*f^6 - 664950*f^5 - 610292*f^4 + 201866*f^3 + 1093230*f^2 + 1295280*f + 4617705;\npoly delta=-f^2 - 2*f - 5;\nfor(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n ideal G=std(I);print("TRANSVERSE_DECOUPLED_GATE");print(2*(kk+1));\n print("TRANSVERSE_DECOUPLED_BASIS");print(G);\n poly remainderPoly=reduce(delta^12,G);print("TRANSVERSE_DECOUPLED_DET_REMAINDER_SIZE");print(size(remainderPoly));\n ideal Sat=sat(G,ideal(delta));print("TRANSVERSE_DECOUPLED_SATURATED");print(Sat);\n if(reduce(1,Sat)==0){break;}\n}\nprint("END");quit;\n'
fpath=p/'transverse_easy_charts.sing';fpath.write_text(source)
log=p/'transverse_easy_charts_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact line-conic covariance charts, beginning with the first-six pure moments and extending only if the determinant-saturated ideal is nonunit; interpretation requires the complete written chart cover.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'transverse_easy_charts_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
