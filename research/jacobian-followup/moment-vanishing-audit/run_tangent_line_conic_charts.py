#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='LIB "elim.lib";\nring rTANGENT_A1B1D0=0,(e,f),dp;\noption(redSB);\nideal All=-28*e^4 - 21*e^3 - 9*e^2 + 7*e*f^2 + e + 3*f^2 + 1,396*e^6 + 858*e^5 + 831*e^4 - 66*e^3*f^2 + 384*e^3 - 36*e^2*f^2 + 94*e^2 - 18*e*f^2 - 2*e + 9*f^4 + 8*f^2 - 1,728*e^8 + 442*e^7 - 325*e^6 - 494*e^5*f^2 - 858*e^5 - 520*e^4*f^2 - 645*e^4 - 260*e^3*f^2 - 258*e^3 + 65*e^2*f^4 - 16*e^2*f^2 - 47*e^2 + 52*e*f^4 + 34*e*f^2 + 2*e + 11*f^4 + 8*f^2 + 1,-12920*e^10 - 32300*e^9 - 36720*e^8 + 6460*e^7*f^2 - 19941*e^7 + 12852*e^6*f^2 - 1939*e^6 + 13107*e^5*f^2 + 5307*e^5 - 1020*e^4*f^4 + 6645*e^4*f^2 + 3845*e^4 - 935*e^3*f^4 + 1770*e^3*f^2 + 1305*e^3 - 465*e^2*f^4 + 30*e^2*f^2 + 183*e^2 + 85*e*f^6 + 15*e*f^4 - 41*e*f^2 - 3*e + 35*f^6 + 45*f^4 + 9*f^2 - 1;\npoly delta=-5*e^2 - 2*e - f^2 - 1;\nfor(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n ideal G=std(I);print("TANGENT_A1B1D0_GATE");print(2*(kk+1));\n print("TANGENT_A1B1D0_BASIS");print(G);\n poly remainderPoly=reduce(delta^12,G);print("TANGENT_A1B1D0_DET_REMAINDER_SIZE");print(size(remainderPoly));\n ideal Sat=sat(G,ideal(delta));print("TANGENT_A1B1D0_SATURATED");print(Sat);\n if(reduce(1,Sat)==0){break;}\n}\nring rTANGENT_A1B0D0=0,(e,f),dp;\noption(redSB);\nideal All=-4*e^4 + e*f^2,132*e^6 - 22*e^3*f^2 + 3*f^4,56*e^8 - 38*e^5*f^2 + 5*e^2*f^4,-152*e^10 + 76*e^7*f^2 - 12*e^4*f^4 + e*f^6;\npoly delta=-f^2;\nfor(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n ideal G=std(I);print("TANGENT_A1B0D0_GATE");print(2*(kk+1));\n print("TANGENT_A1B0D0_BASIS");print(G);\n poly remainderPoly=reduce(delta^12,G);print("TANGENT_A1B0D0_DET_REMAINDER_SIZE");print(size(remainderPoly));\n ideal Sat=sat(G,ideal(delta));print("TANGENT_A1B0D0_SATURATED");print(Sat);\n if(reduce(1,Sat)==0){break;}\n}\nring rTANGENT_A0B0D1=0,(c,f),dp;\noption(redSB);\nideal All=f,c,f^2,c*f;\npoly delta=-c;\nfor(int kk=2;kk<=4;kk++){\n ideal I=0;for(int jj=1;jj<=kk;jj++){I=I,All[jj];}\n ideal G=std(I);print("TANGENT_A0B0D1_GATE");print(2*(kk+1));\n print("TANGENT_A0B0D1_BASIS");print(G);\n poly remainderPoly=reduce(delta^12,G);print("TANGENT_A0B0D1_DET_REMAINDER_SIZE");print(size(remainderPoly));\n ideal Sat=sat(G,ideal(delta));print("TANGENT_A0B0D1_SATURATED");print(Sat);\n if(reduce(1,Sat)==0){break;}\n}\nprint("END");quit;\n'
fpath=p/'tangent_line_conic_charts.sing';fpath.write_text(source)
log=p/'tangent_line_conic_charts_stdout.txt'
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
(p/'tangent_line_conic_charts_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
