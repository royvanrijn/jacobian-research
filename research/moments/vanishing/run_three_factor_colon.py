#!/usr/bin/env python3
"""One checkpointed characteristic-zero colon computation, 50 seconds/768MiB."""
from pathlib import Path
import json,hashlib,subprocess,resource,time
import sympy as S
from verify_nodal_product_reduction import moments
p=Path(__file__).resolve().parent;start=time.monotonic();a,b,c,d,e,f=S.symbols('a b c d e f')
mm=moments([[a,d,e],[d,b,f],[e,f,c]],10)
def singular(v):return str(v).replace('**','^')
source='ring r=0,(a,b,c,d,e,f),dp;\noption(redSB);\nideal I='+','.join(singular(mm[k-1]) for k in [2,4,6,8])+';\npoly P='+singular(mm[9])+';\nprint("BEGIN_STD");\nideal G=std(I);\nprint("STD_DONE");\nprint(size(G));\nideal Q=quotient(G,ideal(P));\nQ=std(Q);\nprint("COLON_DONE");\nprint(size(Q));\nprint(Q);\nprint("MULTIPLIER_CHECK");\nprint(reduce(Q*P,G));\nprint("END");\nquit;\n'
input_path=p/'three_factor_colon.sing';input_path.write_text(source)
log=p/'three_factor_colon_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(input_path)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_SINGULAR_PID',proc.pid,flush=True)
 timed_out=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:
  timed_out=True;proc.kill();code=proc.wait()
text=log.read_text(errors='replace')
complete=code==0 and 'COLON_DONE' in text and 'END' in text and '?' not in text
receipt={'scope':'I8:M10 in the universal six-variable characteristic-zero Gram ring; a terminal resource limit is UNKNOWN, not a unit or theorem.',
         'command':['/usr/bin/Singular','-q',str(input_path)],'wall_limit_seconds':50,'memory_limit_MiB':768,
         'process_pid':proc.pid,'terminal':True,'timed_out':timed_out,'exit_code':code,'completed':complete,
         'classification':'EXACT_OUTPUT_PENDING_AUDIT' if complete else 'UNKNOWN',
         'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),
         'elapsed_seconds':time.monotonic()-start,'last_stdout':text[-1500:]}
(p/'three_factor_colon_exact_run.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
