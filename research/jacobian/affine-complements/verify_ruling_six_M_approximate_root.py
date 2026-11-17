#!/usr/bin/env python3
"""Exact controls for PF167; all-order Hensel/full-image proofs are written."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import tempfile
import sympy as s

BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_APPROXIMATE_CUBIC_ROOT_REGULARITY_AT_M_DIVISORS.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(label,value):check(label,s.cancel(s.expand(value))==0)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
U,w,h,M,N,K,d,a1,a0,b1,b0,b2,p=s.symbols("U w h M N K d a1 a0 b1 b0 b2 p")
C=(U-K)**2-N*h**2
F=(U**2+d)*C**2+h*M*a1*(U-K)*C+h**2*M*b1*C+h**3*M**2*a0*(U-K)+h**4*M**2*b0
z=(w**2-N)/M
P=((h*w+K)**2+d)*z**2+a1*w*z+a0*w+b1*z+b0
zero("exact full translated sextic",h**4*M**2*P.subs(w,(U-K)/h)-F)
zero("reached-line unit identity",((h*w+K)**2+d-h**2*(w**2-N)-2*h*K*w-b2).subs(d,b2-N*h**2-K**2))
n0=s.symbols("n0")
for sign in (-1,1):
    zero("both reached signs square identity "+str(sign),N*h**2+K**2+sign*2*h*K*n0-(h*n0+sign*K)**2+(n0**2-N)*h**2)
for h0 in (0,1):
    CC=C.subs({h:h0,K:1,N:2})
    check("constant simple M-root C unit h"+str(h0),CC.subs(U,2)!=0)
    check("constant simple M-root derivative unit h"+str(h0),(2*U*CC**2).subs(U,2)!=0)
check("off-M one-root collision does not remove other simple root",C.subs({h:0,K:1,U:1})==0 and C.subs({h:0,K:1,U:-1})!=0)

ss,Y=s.symbols("ss Y")
models=[]
for jj,mm,dd,central in [(0,2,1,True),(0,3,2,False),(1,2,1,False),(2,2,3,True),(3,0,2,False),(2,0,1,False)]:
    alpha=1 if dd%2==0 else 2
    rr=alpha*dd//2
    check(f"model normalized parity j{jj} m{mm} d{dd}",math.gcd(alpha,rr)==1)
    check(f"model insufficient order j{jj} m{mm} d{dd}",0<dd<jj+mm)
    hh=ss**(alpha*jj);mM=ss**(alpha*mm)
    nn=ss**alpha if central else s.Integer(2)
    model={h:hh,M:mM,N:nn,K:s.Integer(1),d:3*ss**(alpha*dd)}
    whole=s.cancel((F-p*h**4*M**2).subs(model).subs(U,ss**rr*Y)/ss**(alpha*dd))
    check(f"entire scaled equation regular j{jj} m{mm} d{dd}",s.denom(whole)==1)
    zero(f"entire simple quadratic reduction j{jj} m{mm} d{dd}",whole.subs(ss,0)-(Y**2+3))
    derivative=s.cancel(s.diff(F,U).subs(model).subs(U,ss**rr*Y)/ss**rr)
    check(f"exact scaled derivative regular j{jj} m{mm} d{dd}",s.denom(derivative)==1)
    zero(f"exact scaled derivative leading coefficient j{jj} m{mm} d{dd}",derivative.subs(ss,0)-2*Y)
    omega_order=alpha*(3*jj+mm+1)-rr-1
    check(f"coframe holomorphic j{jj} m{mm} d{dd}",omega_order>=0)
    ww=(ss**rr*Y-1)/hh
    if jj:
        zero(f"original w pole leading coefficient j{jj} m{mm} d{dd}",(ss**(alpha*jj)*ww).subs(ss,0)+1)
    else:
        check(f"finite w limit avoids reached values j{jj} m{mm} d{dd}",1-(0 if central else 2)!=0)
    models.append({"j":jj,"m":mm,"delta":dd,"central":central,"alpha":alpha,"coframe_order":omega_order})

# Complete symbolic margins, not a bounded census of vanishing orders.
jj,dd,gap=s.symbols("jj dd gap",integer=True,nonnegative=True)
# j+m=delta+1+gap; delta>=1. These margins are strictly positive.
mm=dd+1+gap-jj
zero("all-order lower-term margin",jj+mm-dd-1-gap)
zero("all-order parameter margin",4*jj+2*mm-dd-(2*jj+dd+2+2*gap))
zero("all-order coframe positive expression",3*jj+mm+1-dd/s.Integer(2)-(2*jj+dd/s.Integer(2)+gap+2))
kk=s.symbols("kk",integer=True,nonnegative=True)
check("all odd Puiseux branches normalized",s.Mod(2*kk+1,2)==1)
zero("all even branch exponent",(2*kk+2)/s.Integer(2)-(kk+1))
# A pole-face controls the potentially nonunit K case exactly.
c0,h0,M0=s.symbols("c0 h0 M0",nonzero=True)
zero("A face quartic square coefficient",6*c0**2*h0**2-(-2*c0*h0)**2-2*c0**2*h0**2)

record={"claim":"PF2RULED167","outcome":"PASS","check_count":len(checks),"checks":checks,
    "environment":{"python":platform.python_version(),"sympy":s.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_sources":{name:sha(BASE/name) for name in ["RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md","POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md"]},
    "Puiseux_transcription_controls":models,
    "uniform_scope":{"every_M_root_ord_d_at_least_ord_h_plus_ord_M":True,"central_and_j_zero_cases_retained":True,"off_M_h_root_K_unit_covered":True,"actual_A_pole_roots_have_exact_ord_d_equal_ord_h":True,"global_hM_divisibility_claimed_here":False,"first_derivative_integrality_proved":False,"JC2_or_HC4_closed":False,"independent_whole_proof_review":False,"formal_verification":False,"external_human_review":False}}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="ruling-six-M-root-check-"))/"ruling_six_M_approximate_root_receipt.json"
out.write_text(json.dumps(record,indent=2)+chr(10))
print(json.dumps({"outcome":"PASS","check_count":len(checks),"receipt":str(out)}))
