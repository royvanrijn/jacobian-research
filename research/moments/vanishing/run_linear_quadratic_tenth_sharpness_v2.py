#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,resource,subprocess,time
p=Path(__file__).resolve().parent;start=time.monotonic()
source='LIB "elim.lib";\nring r=0,(d2,d3,h1),dp;\noption(redSB);\nideal I=-112*d2^2 - 320*d2*h1^2 + 896*d2*h1 - 360*d2 - 640*d3*h1 + 1120*d3 - 112*h1^2 + 480*h1 - 495,96*d2^3 + 1008*d2^2*h1^2 - 1728*d2^2*h1 + 144*d2^2 + 2016*d2*d3*h1 - 2016*d2*d3 + 792*d2*h1^2 - 2736*d2*h1 + 2466*d2 + 1008*d3^2 - 448*d3*h1^3 + 2016*d3*h1^2 - 2520*d3*h1 + 672*d3 + 96*h1^3 - 549*h1^2 + 1080*h1 - 723,2816*d2^4 + 4096*d2^3*h1^2 - 90112*d2^3*h1 + 81152*d2^3 + 8192*d2^2*d3*h1 - 101376*d2^2*d3 + 45056*d2^2*h1^4 - 262144*d2^2*h1^3 + 451584*d2^2*h1^2 - 198656*d2^2*h1 - 20960*d2^2 - 22528*d2*d3^2 + 262144*d2*d3*h1^3 - 1155072*d2*d3*h1^2 + 1511424*d2*d3*h1 - 478720*d2*d3 + 26624*d2*h1^4 - 212992*d2*h1^3 + 590592*d2*h1^2 - 657920*d2*h1 + 235728*d2 + 270336*d3^2*h1^2 - 901120*d3^2*h1 + 760320*d3^2 - 8192*d3*h1^4 + 131072*d3*h1^3 - 583680*d3*h1^2 + 1026560*d3*h1 - 634304*d3 + 2816*h1^4 - 28672*h1^3 + 102592*h1^2 - 155520*h1 + 85227;\npoly discPoly=d2^2-4*d3-4*d2^3-27*d3^2+18*d2*d3;\npoly goodPoly=d3*discPoly*(4*d2+4*h1-5);\nideal G=std(I);print("EIGHTH_PREFIX_DIMENSION");print(dim(G));\nprint("EIGHTH_PREFIX_LENGTH");print(vdim(G));\nideal Good=sat(G,ideal(goodPoly));\nprint("GOOD_WITNESS_LENGTH");print(vdim(Good));print("GOOD_WITNESS_BASIS");print(Good);\nideal Bad=std(Good+ideal(goodPoly));print("BAD_LOCUS_IDEAL");print(Bad);\nprint("END");quit;\n'
fpath=p/'linear_quadratic_tenth_sharpness_v2.sing';fpath.write_text(source)
log=p/'linear_quadratic_tenth_sharpness_v2_stdout.txt'
def limit():
 resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
 resource.setrlimit(resource.RLIMIT_CPU,(20,20))
with log.open('w') as out:
 proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
 print('LIVE_PID',proc.pid,flush=True);timeout=False
 try:code=proc.wait(timeout=20)
 except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact algebraic existence certificate for first-eight-zero/first-nine-zero transverse line-conic cubics: full invariant prefix through8 with nonzero determinant, distinct quadratic eigenvalues and nonzero H2. Realization uses the written partial-fraction construction.',
'wall_limit_seconds':20,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,'completed':code==0 and 'END' in result and '?' not in result,
'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-3000:]}
(p/'linear_quadratic_tenth_sharpness_v2_run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
