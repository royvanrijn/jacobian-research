#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
text=(p/'quaternary_modular_section_stdout.txt').read_text()
basis=text.split('SECTION_BASIS\n')[1].split('SATURATION_START')[0].strip()
source='ring r=32003,(a,b,d,e,f),dp;\noption(redSB);\nLIB "elim.lib";\nideal G='+basis+';\nattrib(G,"isSB",1);\n'
source+='print("LINEAR_SATURATION_START");\nideal S=sat(G,ideal(a+2b+3d+5e+7f));\nprint("SECTION_HILBERT");print(hilb(G,1));\nprint("SATURATION_HILBERT");print(hilb(S,1));\nprint("SATURATION_BASIS_SIZE");print(size(S));\nprint("SATURATION_BASIS");print(S);\nprint("END");quit;\n'
fpath=p/'quaternary_modular_linear_saturation.sing';fpath.write_text(source)
log=p/'quaternary_modular_linear_saturation_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Principal linear saturation of the retained exact finite-field section standard basis; finite Hilbert difference is required before identifying origin torsion.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'quaternary_modular_linear_saturation_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
