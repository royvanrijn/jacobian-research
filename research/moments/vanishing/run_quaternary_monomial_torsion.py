#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source=(p/'quaternary_modular_regularity.sing').read_text().split('ideal L=std(lead(G));')[0]
source+='ideal L=std(lead(G));print("MONOMIAL_SATURATION_START");\nideal LS=sat(L,maxideal(1));\nprint("LEADING_HILBERT");print(hilb(L,1));\nprint("SATURATION_HILBERT");print(hilb(LS,1));\nprint("LEADING_SATURATION_BASIS");print(LS);\nprint("END");quit;\n'
fpath=p/'quaternary_monomial_torsion.sing';fpath.write_text(source)
log=p/'quaternary_monomial_torsion_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Origin saturation of the monomial leading ideal of the finite-field section basis, for a rigorously lifted upper bound on characteristic-zero origin torsion.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'quaternary_monomial_torsion_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
