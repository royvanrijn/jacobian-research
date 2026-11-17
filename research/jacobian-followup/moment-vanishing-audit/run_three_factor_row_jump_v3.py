#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent
start=time.monotonic()
import sympy as S
from verify_nodal_product_reduction import moments
a,b,c,d,e,f,z=S.symbols('a b c d e f z')
raw=moments([[a,d,e],[d,b,f],[e,f,c]],10)
ar=-2*(d*d+e*e+4*z*d*e)/(1+2*z*z)
row=[S.fraction(S.cancel(v.subs({b:1,c:1,f:z}).subs(a,ar)))[0] for v in raw[1::2]]
source='ring r=0,(d,e,z),dp;\noption(redSB);\nideal I='+','.join(str(v).replace('**','^') for v in row)+';\n'
res=S.factor(S.resultant(row[1].subs(e,1),row[2].subs(e,1),d))
factors=S.factor_list(res,z)[1]
D=next(q for q,k in factors if S.degree(q,z)==18)
assert S.degree(D,z)==18 and S.gcd(D,1+2*z*z)==1
source+='poly D='+str(D).replace('**','^')+';\n'
source+='\nprint("ROW_BASE");ideal H=std(I);print(size(H));print(hilb(H,1));\nprint("SPECIAL_RESULTANT_FIBER");ideal F=std(I+ideal(D));print(vdim(F));print(hilb(F,1));\nprint("ORDINARY_FIBER");ideal F0=std(I+ideal(z));print(vdim(F0));\nprint("END");quit;\n'
fpath=p/'three_factor_row_jump_v3.sing';fpath.write_text(source)
log=p/'three_factor_row_jump_v3_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(35,35))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=35)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Row slice exact generic and resultant-fiber Hilbert data; no full primary claim.',
'wall_limit_seconds':35,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'three_factor_row_jump_v3_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
