#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent
start=time.monotonic()
import sympy as S
from verify_nodal_product_reduction import moments
a,b,c,d,e,f=S.symbols('a b c d e f');raw=moments([[a,d,e],[d,b,f],[e,f,c]],10)
source=''
for label,params,variables in [('PLANE_AD','a,d','b,c,e,f'),('AXIS_A','a','b,c,d,e,f'),('AXIS_D','d','a,b,c,e,f'),('ORIGIN','','a,b,c,d,e,f')]:
 field='(0,'+params+')' if params else '0'
 source+='ring r'+label+'='+field+',('+variables+'),dp;\noption(redSB);\n'
 source+='ideal I='+','.join(str(v).replace('**','^') for v in raw[1::2])+';\n'
 source+='ideal G=std(I);ideal m='+variables+';\n'
 source+='print("'+label+'_START");ideal K=std(quotient(G,m));\n'
 source+='print("'+label+'_SOCLE_GENERATORS");print(size(reduce(K,G)));\n'
 source+='print("'+label+'_FORWARD");print(size(reduce(K*m,G)));print("'+label+'_REVERSE");print(size(reduce(G,K)));\n'
 source+='print("'+label+'_WITNESS");for(int jj'+label+'=1;jj'+label+'<=size(K);jj'+label+'++){if(reduce(K[jj'+label+'],G)!=0){print(reduce(K[jj'+label+'],G));break;}}\n'
 source+='print("'+label+'_I_HILBERT");print(hilb(G,1));print("'+label+'_K_HILBERT");print(hilb(K,1));\n'
source+='print("END");quit;\n'
fpath=p/'three_factor_associated_primes_v2.sing';fpath.write_text(source)
log=p/'three_factor_associated_primes_v2_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Four symmetry-reduced localized socle tests for candidate embedded coordinate primes; coordinate-support theorem proved separately.',
'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'three_factor_associated_primes_v2_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
