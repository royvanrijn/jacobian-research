#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source=(p/'quaternary_modular_linear_saturation.sing').read_text().split('print("LINEAR_SATURATION_START")')[0]
source+='ideal L=std(lead(G));print("REGULARITY_START");\nprint("GENERATOR_COUNTS");print(size(G));print(size(L));\nprint("LEADING_HILBERT");print(hilb(L,1));\nresolution RL=mres(L,0);intmat B=betti(RL);\nprint("BETTI");print(B);\nint regBound=-1;\nfor(int bi=1;bi<=nrows(B);bi++){\n for(int bj=1;bj<=ncols(B);bj++){if(B[bi,bj]!=0){regBound=bi-1;}}\n}\nprint("REGULARITY_BOUND");print(regBound);\nprint("RESOLUTION_LENGTH");print(size(RL));\nprint("END");quit;\n'
fpath=p/'quaternary_modular_regularity.sing';fpath.write_text(source)
log=p/'quaternary_modular_regularity_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Monomial leading-ideal regularity of the retained finite-field section basis. Exact Hilbert matching and the written integral Koszul-rank argument are required for a characteristic-zero regularity bound.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'quaternary_modular_regularity_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
