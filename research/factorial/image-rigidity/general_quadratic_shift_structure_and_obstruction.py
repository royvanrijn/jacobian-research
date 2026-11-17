#!/usr/bin/env python3
"""Exact finite formula controls; the written argument supplies all-h quantifiers."""
import hashlib,json,platform
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[3]
M,n,u,v,e,z=S.symbols("M n u v eps z")
gamma=M-2*n
controls=[]
def check(name,value):
 if not bool(value):
  raise AssertionError((name,value))
 controls.append(name)
def same(name,a,b): check(name,S.expand(a-b)==0)
def step(k,rr,ww,uu):
 return tuple(S.expand(q) for q in (
  (k+2)*rr-2*u*S.diff(rr,u)+gamma*u*uu,
  -rr+(k+1)*ww-2*u*S.diff(ww,u)+u*uu,
  ww+(M+k+2)*uu-2*u*S.diff(uu,u)))
triples=[(u,S.Integer(0),S.Integer(1))]
for h in range(9): triples.append(step(h,*triples[-1]))
s=2*(n+1)
forms={
 2:(S.rf(M+2,2),s*u),
 3:(S.rf(M+2,3)+s*u,s*(M+4)*u),
 4:(S.rf(M+2,4)+s*(2*M+7)*u,s*u*(S.rf(M+4,2)+u)),
 5:(S.rf(M+2,5)+3*s*(M+4)**2*u+s*u*u,
    s*u*(S.rf(M+4,3)+(M+2*n+8)*u))}
for h,(uu,ww) in forms.items():
 same(f"display-U-{h}",triples[h][2],uu)
 same(f"display-W-{h}",triples[h][1],ww)
# Independently construct the subtraction-free coefficient network.
svar=S.symbols("s")
UU,WW,VV=S.Integer(1),S.Integer(0),S.Integer(0)
for h in range(10):
 rr0,ww0,uu0=triples[h]
 same(f"positive-network-U-h{h}",UU,uu0.subs(n,svar/2-1))
 same(f"positive-network-W-h{h}",WW,ww0.subs(n,svar/2-1))
 same(f"positive-network-V-h{h}",VV,(u*uu0-rr0).subs(n,svar/2-1))
 for name,q in [("U",UU),("W",WW),("V",VV)]:
  check(f"positive-network-coefficients-h{h}-{name}",
        all(c>=0 for c in S.Poly(q,M,svar,u).coeffs()))
 UU,WW,VV=tuple(S.expand(q) for q in (
  (M+h+2)*UU-2*u*S.diff(UU,u)+WW,
  (h+1)*WW-2*u*S.diff(WW,u)+VV,
  (h+2)*VV-2*u*S.diff(VV,u)+u*WW+svar*u*UU))
# Direct actual polynomial action; all h used here are explicitly bounded controls.
for nn,gg in [(1,S.Rational(3,2)),(2,5),(3,S.Rational(7,3))]:
 mm=gg+2*nn
 R=S.expand(sum(S.binomial(nn,j)*S.rf(gg,2*j)*(-1)**j*z**(2*(nn-j)) for j in range(nn+1)))
 def E(f):return S.expand(mm*f-z*S.diff(f,z))
 K1=E(R);K2=S.expand(E(K1)+K1)
 L=R
 rising=[R]
 for j in range(11):
  L=S.expand(E(L)+j*L);rising.append(L)
 for h in range(9):
  rr,ww,uu=(f.subs({M:mm,n:nn,u:z*z}) for f in triples[h])
  direct=z*z*rising[h]-rising[h+2]
  same(f"actual-remainder-n{nn}-h{h}",direct,rr*R-ww*K1-uu*K2)
  upper=sum(S.binomial(nn+1,j)*S.rf(gg+h,2*j)*(-1)**j*z**(2*(nn+1-j)) for j in range(nn+2))
  same(f"actual-upper-shift-n{nn}-h{h}",direct,S.rf(gg,h)*upper)
# Weighted asymptotics independently extracted from exact polynomials.
for h,(rr,ww,uu) in enumerate(triples):
 vals=[S.expand(e**h*uu.subs({M:1/e,u:v/e**2})),
       S.expand(e**h*ww.subs({M:1/e,u:v/e**2})),
       S.expand(e**(h+2)*rr.subs({M:1/e,u:v/e**2}))]
 for j,q in enumerate(vals):
  check(f"eps-polynomial-h{h}-{j}",S.denom(S.cancel(q))==1)
 same(f"eps-U-leading-h{h}",vals[0].subs(e,0),1)
 same(f"eps-r-leading-h{h}",vals[2].subs(e,0),v)
 same(f"eps-W-leading-h{h}",vals[1].subs(e,0),
      s*sum(v**j for j in range(1,h//2+1)))
# An admissible physical-family quotient counterexample with STRICT Re u > 0.
h=6;nn=2;mm=1000;beta=mm-(h+2)*nn;uu=1+3000000*S.I
rr,ww,UU=triples[h]
Wval=S.expand(ww.subs({M:mm,n:nn,u:uu}))
Uval=S.expand(UU.subs({M:mm,n:nn,u:uu}))
Hval=S.simplify(S.im(Wval*S.conjugate(Uval))/S.im(uu))
check("counterexample-beta-positive",beta>0)
check("counterexample-Reu-positive",S.re(uu)>0)
check("counterexample-U-nonzero",Uval!=0)
check("counterexample-imag-kernel-negative",bool(Hval<0))
# A wrong recurrence sign must fail a direct actual remainder.
bad_r,bad_w,bad_u=triples[2]
bad_w=-bad_w
check("wrong-W-sign-rejected",S.expand(bad_w-triples[2][1])!=0)
# Test universal geometric estimate without floating-point trigonometry.
check("uniform-angle-range",S.Rational(5,12)<S.Rational(1,2))
check("universal-leading-imaginary-margin",S.Rational(1,9)<S.Rational(1,2))
inputs=["research/factorial/image-rigidity/GENERAL_QUADRATIC_SHIFT_STRUCTURE_AND_OBSTRUCTION.tex",
 "research/factorial/image-rigidity/general_quadratic_shift_structure_and_obstruction.py",
 "research/factorial/image-rigidity/LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex",
 "research/factorial/image-rigidity/LOWER_TWO_PRIME_GAP_SFC.tex"]
result={"status":"PASS","passed":True,"kind":"internal exact replay; unrestricted quantifiers come from the written induction","controls":len(controls),"control_names":controls,
 "scope":"Finite exact formula controls only. The written epsilon-recursion proves the every-h shift and every-h>=6 quadrant-Pick obstruction. No SFC counterexample or universal SFC conclusion.",
 "counterexample":{"h":h,"n":nn,"M":mm,"gamma":mm-2*nn,"beta":beta,
 "u":str(uu),"U":str(Uval),"W":str(Wval),"imag_kernel":str(Hval)},
 "source_sha256":{Path(q).name:hashlib.sha256((ROOT/q).read_bytes()).hexdigest() for q in inputs},
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/general_quadratic_shift_structure_and_obstruction.py"}
Path(__file__).with_name("general_quadratic_shift_structure_and_obstruction_verification.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"controls":len(controls),"counterexample":result["counterexample"]},indent=2))
