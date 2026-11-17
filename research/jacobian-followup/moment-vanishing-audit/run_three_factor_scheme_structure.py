#!/usr/bin/env python3
"""Bounded exact follow-up: colon via syzygies, radical, depth and nilpotence."""
from pathlib import Path
import subprocess,resource,time,json,hashlib
import sympy as S
from verify_nodal_product_reduction import moments

p=Path(__file__).resolve().parent
start=time.monotonic()
a,b,c,d,e,f=S.symbols('a b c d e f')
mm=moments([[a,d,e],[d,b,f],[e,f,c]],10)
def text(v):return str(v).replace('**','^')
source='ring r=0,(a,b,c,d,e,f),dp;\noption(redSB);\nideal I='+','.join(text(mm[k-1]) for k in [2,4,6,8])+';\npoly P='+text(mm[9])+';\n'
source+='''
ideal G8=std(I);
ideal Iall=I,P;
ideal G10=std(Iall);
print("GLOBAL_BASES");print(size(G8));print(size(G10));
ideal Q=std(quotient(G8,ideal(P)));
module Z=syz(Iall);
ideal Q2;
for(int j=1;j<=ncols(Z);j++){Q2[j]=Z[5,j];}
Q2=std(Q2);
print("SYZYGY_COLON_MUTUAL_REDUCTION");print(reduce(Q,Q2));print(reduce(Q2,Q));
ideal Rad=d,e,f,a*b,a*c,b*c;
print("COLON_RADICAL_FORWARD");print(reduce(Q,std(Rad)));
ideal radicalWitness=d6,e6,f6,5a2b2-12abd2,5a2c2-12ace2,5b2c2-12bcf2;
print("COLON_RADICAL_REVERSE_WITNESSES");print(reduce(radicalWitness,Q));
print("COLON_LINEAR_FORM_REGULAR");print(reduce(quotient(Q,ideal(a+b+c)),Q));
print("COLON_HILBERT");print(hilb(Q,2));
resolution R10=mres(Iall,0);
print("I10_RESOLUTION_LENGTH");print(size(R10));print(betti(R10));
ideal J=a*f,b*e,c*d,a*b*c,d*e*f;
ideal K=J;
for(int k=1;k<=19;k++){
 if(k>=8){print("NIL_POWER");print(k);print(size(K));
 print("I8_REMAINDER_SIZE");print(size(reduce(K,G8)));
 print("I10_REMAINDER_SIZE");print(size(reduce(K,G10)));}
 if(k<19){K=std(K*J);}
}
print("END");quit;
'''
fpath=p/'three_factor_scheme_structure.sing';fpath.write_text(source)
log=p/'three_factor_scheme_structure_stdout.txt'
def limit():
    resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
    resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
    proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
    print('LIVE_SINGULAR_PID',proc.pid,flush=True)
    timeout=False
    try:code=proc.wait(timeout=50)
    except subprocess.TimeoutExpired:
        timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
complete=code==0 and 'END' in result and '?' not in result
rec={'scope':'Exact universal six-Gram-variable three-factor scheme structure; completed stages require audit. Resource limits never imply a mathematical result.',
     'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':complete,
     'classification':'EXACT_OUTPUT_PENDING_AUDIT' if complete else 'UNKNOWN_OR_PARTIAL_OUTPUT_PENDING_AUDIT',
     'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),
     'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-7000:]}
(p/'three_factor_scheme_structure_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
