#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent
start=time.monotonic()
analysis=json.loads((p/'origin_linear_map_analysis.json').read_text())
c=json.loads((p/'three_factor_colon_geometry_certificate.json').read_text())
num=c['first_Hilbert_numerators']['I10'][:]
while len(num)<37:num.append(0)
from math import comb
for ds,mul in analysis['kernel_hilbert'].items():
 for j in range(6):num[int(ds)+1+j]+=int(mul)*comb(5,j)*(-1)**j
while num[-1]==0:num.pop()
# Singular's first Hilbert numerator intvec has a final zero sentinel.
source=(p/'ternary_quaternary_null_section.sing').read_text().split('print("SECTION_STANDARD_START")')[0]
source+='intvec hi='+','.join(map(str,num+[0]))+';\n'
source+=r'''print("HILBERT_DRIVEN_START");
ideal G=std(H,hi,intvec(1,1,1,1,1));
print("SECTION_BASIS_SIZE");print(size(G));
print("SECTION_HILBERT");print(hilb(G,1));
print("SECTION_BASIS");print(G);
print("GENERATOR_INCLUSION");print(size(reduce(H,G)));
print("END");quit;
'''
fpath=p/'quaternary_hilbert_section.sing';fpath.write_text(source)
log=p/'quaternary_hilbert_section_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=50)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Hilbert-driven exact rank-one-functional section standard basis, with independently derived kernel Hilbert correction on the certified77-dimensional origin torsion.',
'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'Hilbert_numerator':num,'stdout_tail':result[-500:]}
(p/'quaternary_hilbert_section_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
