#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent
start=time.monotonic()
source=(p/'three_factor_colon_geometry_v3.sing').read_text().split('ideal G8=std(I);')[0]
source+='\nLIB "primdec.lib";\nideal Iall=I,P;\nprint("PRIMARY_START");\nlist D=primdecGTZ(Iall);\nprint("PRIMARY_COUNT");print(size(D));\nfor(int j=1;j<=size(D);j++){\n print("COMPONENT");print(j);print("PRIMARY_IDEAL");print(std(D[j][1]));print("ASSOCIATED_PRIME");print(std(D[j][2]));\n}\nprint("END");quit;\n'
fpath=p/'three_factor_full_primary.sing';fpath.write_text(source)
log=p/'three_factor_full_primary_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Full characteristic-zero I10 primary decomposition; exploratory output requires exact intersection and primary certification before theorem promotion.',
'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'three_factor_full_primary_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
