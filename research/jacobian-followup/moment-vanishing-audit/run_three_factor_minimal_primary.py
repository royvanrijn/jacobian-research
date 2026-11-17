#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent
start=time.monotonic()
source=(p/'three_factor_colon_geometry_v3.sing').read_text().split('ideal G8=std(I);')[0]
source+='\nLIB "elim.lib";\nideal Iall=I,P;ideal G=std(Iall);\nprint("SAT_START");\nlist LR=sat_with_exp(G,ideal(b*c));ideal QR=std(LR[1]);\nprint("ROW_EXPONENT");print(LR[2]);print("ROW_BASIS");print(QR);\nlist LS=sat_with_exp(G,ideal(a*d*e));ideal QS=std(LS[1]);\nprint("STAR_EXPONENT");print(LS[2]);print("STAR_BASIS");print(QS);\nmap ab=r,b,a,c,d,f,e;map ac=r,c,b,a,f,e,d;\nideal QB=ab(QR);ideal QC=ac(QR);ideal SB=ab(QS);ideal SC=ac(QS);\nprint("UNMIXED_START");\nideal U=std(intersect(intersect(intersect(intersect(intersect(QR,QB),QC),QS),SB),SC));\nprint("UNMIXED_BASIS");print(U);\nprint("FULL_IN_UNMIXED");print(size(reduce(G,U)));\nprint("UNMIXED_IN_FULL");print(size(reduce(U,G)));\nprint("FULL_HILBERT");print(hilb(G,1));\nprint("UNMIXED_HILBERT");print(hilb(U,1));\nprint("END");quit;\n'
fpath=p/'three_factor_minimal_primary.sing';fpath.write_text(source)
log=p/'three_factor_minimal_primary_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Two characteristic-zero isolated-primary saturations and intersection of their six symmetric copies; primaryness proved separately by normalized-family freeness.',
'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'three_factor_minimal_primary_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
