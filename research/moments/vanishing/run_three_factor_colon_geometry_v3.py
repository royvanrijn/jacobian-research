#!/usr/bin/env python3
"""Bounded colon geometry and nilpotence run; no raw-generator syzygy stage."""
from pathlib import Path
import subprocess,resource,time,json,hashlib
import sympy as S
from verify_nodal_product_reduction import moments

p=Path(__file__).resolve().parent;start=time.monotonic()
a,b,c,d,e,f=S.symbols('a b c d e f');mm=moments([[a,d,e],[d,b,f],[e,f,c]],10)
def singular(v):return str(v).replace('**','^')
source='ring r=0,(a,b,c,d,e,f),dp;\noption(redSB);\nideal I='+','.join(singular(mm[k-1]) for k in [2,4,6,8])+';\npoly P='+singular(mm[9])+';\n'
source+='''
ideal G8=std(I);
ideal Iall=I,P;
ideal G10=std(Iall);
print("GLOBAL_BASES");print(size(G8));print(size(G10));
ideal J=a*f,b*e,c*d,a*b*c,d*e*f;
ideal K=J;
for(int k=1;k<=14;k++){
 if(k>=12){print("NIL_POWER");print(k);print(size(K));
 print("I8_REMAINDER_SIZE");print(size(reduce(K,G8)));
 print("I10_REMAINDER_SIZE");print(size(reduce(K,G10)));}
 if(k<14){K=std(K*J);}
}
print("NIL_DONE");
ideal Q=std(quotient(G8,ideal(P)));
ideal Rad=d,e,f,a*b,a*c,b*c;
print("COLON_GENERATORS");print(Q);
print("COLON_FORWARD_MULTIPLIER");print(size(reduce(Q*P,G8)));
print("HILBERT_FIRST_I8");print(hilb(G8,1));
print("HILBERT_FIRST_I10");print(hilb(G10,1));
print("HILBERT_FIRST_COLON");print(hilb(Q,1));
print("COLON_RADICAL_FORWARD");print(size(reduce(Q,std(Rad))));
ideal reverseWitness=d6,e6,f6,5a2b2-12abd2,5a2c2-12ace2,5b2c2-12bcf2;
print("COLON_RADICAL_REVERSE");print(size(reduce(reverseWitness,Q)));
print("COLON_LINEAR_REGULAR");print(size(reduce(quotient(Q,ideal(a+b+c)),Q)));
print("COLON_HILBERT");print(hilb(Q,2));
LIB "elim.lib";
list Lsat=sat_with_exp(Q,ideal(a+b+c));
ideal Qsat=std(Lsat[1]);
int exponent=Lsat[2];
print("SATURATION_POWER");print(exponent);
print("SATURATION_HILBERT");print(hilb(Qsat,2));
list La=sat(Q,ideal(a));ideal Qa=std(La[1]);
list Lb=sat(Q,ideal(b));ideal Qb=std(Lb[1]);
list Lc=sat(Q,ideal(c));ideal Qc=std(Lc[1]);
ideal AxisPart=std(intersect(intersect(Qa,Qb),Qc));
print("SATURATION_INTERSECTION_CHECK");print(size(reduce(AxisPart,Qsat)));print(size(reduce(Qsat,AxisPart)));
ideal OriginPart=std(Q+ideal((a+b+c)^exponent));
ideal AllPart=std(intersect(Qsat,OriginPart));
print("PRIMARY_DECOMPOSITION_CHECK");print(size(reduce(AllPart,Q)));print(size(reduce(Q,AllPart)));
print("AXIS_PRIMARY_BASIS");print(Qa);
print("AXIS_RADICAL_FORWARD");print(size(reduce(Qa,std(ideal(b,c,d,e,f)))));
print("AXIS_RADICAL_REVERSE");print(size(reduce(ideal(b6,c6,d6,e6,f6),Qa)));
print("AXIS_PARAMETER_REGULAR");print(size(reduce(quotient(Qa,ideal(a)),Qa)));
print("SATURATION_EXPONENT_CHECK");print(size(reduce(Qsat*((a+b+c)^exponent),Q)));
print("ORIGIN_COMPONENT_LENGTH");print(vdim(OriginPart));
print("AXIS_GENERIC_LENGTH");print(vdim(std(Q+ideal(a-1))));
ideal N8=std(J^13);ideal N10=std(J^12);
print("I8_LOWER_WITNESS");
for(int ii=1;ii<=size(N8);ii++){if(reduce(N8[ii],G8)!=0){print(N8[ii]);print(reduce(N8[ii],G8));break;}}
print("I10_LOWER_WITNESS");
for(int jj=1;jj<=size(N10);jj++){if(reduce(N10[jj],G10)!=0){print(N10[jj]);print(reduce(N10[jj],G10));break;}}
print("END");quit;
'''
fpath=p/'three_factor_colon_geometry_v3.sing';fpath.write_text(source)
log=p/'three_factor_colon_geometry_v3_stdout.txt'
def limit():
    resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
    resource.setrlimit(resource.RLIMIT_CPU,(50,50))
with log.open('w') as out:
    proc=subprocess.Popen(['/usr/bin/Singular','-q',str(fpath)],stdout=out,stderr=subprocess.STDOUT,preexec_fn=limit)
    print('LIVE_SINGULAR_PID',proc.pid,flush=True);timeout=False
    try:code=proc.wait(timeout=50)
    except subprocess.TimeoutExpired:timeout=True;proc.kill();code=proc.wait()
result=log.read_text(errors='replace')
rec={'scope':'Exact radical-power and colon geometry calculation; raw syzygy timeout remains UNKNOWN; invalid v1/v2 inputs preserved; v3 uses documented sat_with_exp and Hilbert colon completeness.',
     'wall_limit_seconds':50,'memory_limit_MiB':768,'terminal':True,'timed_out':timeout,'exit_code':code,
     'completed':code==0 and 'END' in result and '?' not in result,
     'input_sha256':hashlib.sha256(fpath.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),
     'elapsed_seconds':time.monotonic()-start,'stdout_tail':result[-7000:]}
(p/'three_factor_colon_geometry_v3_run.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
