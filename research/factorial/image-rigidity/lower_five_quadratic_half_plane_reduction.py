#!/usr/bin/env python3
"""Exact algebra controls; finite stability samples are not a uniform proof."""
from pathlib import Path
import os,json,hashlib,math,platform
import sympy as S
from flint import fmpz, fmpq
HERE=Path(__file__).resolve().parent
INPUTS={
"CUBIC_ADJOINT_STABILITY_SFC.tex":"ea2510e511545b9dd7a49a0796d75c45575687e60770c3359d38c5bc99cb4125",
"LOWER_FIVE_QUADRATIC_RESIDUE_REDUCTION.tex":"343cd52be51b99fe8a47d5ff2cf493abe1e5503e47da466cf44480808981e5aa"}
for name,digest in INPUTS.items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
checks=0
n,u,z,y=S.symbols("n u z y",real=True)
M=7*n+1;s=2*(n+1)
A=S.prod(M+j for j in range(2,7));b=3*s*(M+4)**2
k=(M+4)*(M+5)*(M+6);c=9*(n+1)
U=A+b*u+s*u**2;W=s*u*(k+c*u)
lam=M*(b*b-2*A*s)-s*c*A+s*k*b
expected=2*(n+1)**2*(7*n+3)*(7*n+5)*(343*n**3+2303*n**2+2926*n+1002)
assert S.expand(lam-expected)==0
checks+=1
lhs=S.expand(M*U.subs(u,S.I*y)*U.subs(u,-S.I*y)
 +(W.subs(u,S.I*y)*U.subs(u,-S.I*y)+W.subs(u,-S.I*y)*U.subs(u,S.I*y))/2)
assert S.expand(lhs-(M*A*A+expected*y*y+s*s*(M+c)*y**4))==0
checks+=1
# Generic even/odd raising identity, without substituting a finite B.
b0=S.Function("b0")(u);b1=S.Function("b1")(u)
bb=b0.subs(u,z*z)+z*b1.subs(u,z*z)
cc=S.expand((z+M+1)*bb-z*S.diff(bb,z))
odd=S.simplify((cc-cc.subs(z,-z))/(2*z))
target=(b0+M*b1-2*u*S.diff(b1,u)).subs(u,z*z)
assert S.simplify(odd-target)==0
checks+=1
def routh(coeffs):
    arr=[fmpq(int(x)) for x in reversed(coeffs)]
    r=arr[::2]; t=arr[1::2]
    width=max(len(r),len(t)); r += [fmpq(0)]*(width-len(r));t += [fmpq(0)]*(width-len(t))
    if r[0]<=0:return False
    for j in range(len(coeffs)-1):
        if t[0]<=0:return False
        q=[(t[0]*r[i+1]-r[0]*t[i+1])/t[0] for i in range(width-1)]+[fmpq(0)]
        r,t=t,q
    return True
def bc(nn):
    gg=5*nn+1
    B=S.Integer(1)
    for j in range(nn):
        B=S.expand((z+2*j+gg+1)*B-z*S.diff(B,z))
    C=S.expand((z+2*nn+gg+1)*B-z*S.diff(B,z))
    return S.Poly(B,z),S.Poly(C,z)
samples=[]
for nn in [1,2,3,5,7,10,11,13,20,30]:
    B,C=bc(nn)
    B1=sum(B.nth(2*j+1)*u**j for j in range((nn-1)//2+1))
    C1=sum(C.nth(2*j+1)*u**j for j in range(nn//2+1))
    hh=S.Poly(S.expand(U.subs(n,nn)*C1+W.subs(n,nn)*B1),u)
    assert routh([hh.nth(j) for j in range(hh.degree()+1)])
    samples.append({"n":nn,"H1_hurwitz":True,"degree":hh.degree()})
    checks+=1
out={"passed":True,"uniform_H1_half_plane_exclusion_proved":True,
"all_index_support_five_seven_sfc":False,
"symbolic_imaginary_axis_identity":True,"linear_factor":str(expected),
"checks":checks,"finite_H1_Routh_controls":samples,
"source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
 for name in list(INPUTS)+["LOWER_FIVE_QUADRATIC_HALF_PLANE_REDUCTION.tex",Path(__file__).name]},
"python":platform.python_version(),"sympy":S.__version__,
"command":"python research/factorial/image-rigidity/lower_five_quadratic_half_plane_reduction.py",
"independent_replay":False,"formal_verification":False,"external_review":False}
(HERE/"lower_five_quadratic_half_plane_reduction_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"passed":True,"checks":checks,"finite_Routh_indices":[x["n"] for x in samples]}))
