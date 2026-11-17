#!/usr/bin/env python3
"""Exact symbolic finite Laurent controls for the exceptional AB correction."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,platform,tempfile
HERE=Path(__file__).resolve().parent
# Variables (t,rho,a,b,tau); only t has possibly negative exponents.
ZERO=(0,0,0,0,0)
ONE={ZERO:Q(1)}
def clean(p):return {e:c for e,c in p.items() if c}
def add(*polys):
 out={}
 for p in polys:
  for e,c in p.items():out[e]=out.get(e,Q(0))+c
 return clean(out)
def scale(p,q):return clean({e:c*Q(q) for e,c in p.items()})
def mul(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   ef=tuple(a+b for a,b in zip(e,f));out[ef]=out.get(ef,Q(0))+c*d
 return clean(out)
def power(p,n):
 out=ONE
 for _ in range(n):out=mul(out,p)
 return out
def mon(n,variable=0):
 e=[0]*5;e[variable]=n;return {tuple(e):Q(1)}
def coeff(p,n):return {(0,)+e[1:]:c for e,c in p.items() if e[0]==n}
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 rho,a,b,tau=[mon(1,i) for i in range(1,5)]
 controls=[]
 for p,q0,q1 in [(2,2,2),(4,6,16)]:
  for eps0 in [Q(1),Q(2)]:
   for N in [3,5,9]:
    D=2*N+1;pN=p*D//2
    q0N=q0*D//2;q1N=q1*D//2
    Lprev=(p-1)*q0N;LI=(p-1)*(q0N+q1N)
    r=LI-pN;gamma=-pN+N+1;qI=r+N+1
    assert p*q1N-qI==N
    if p==2:
     # W=z*eps0*(1+rho/eps0^2*z^-1)^(1/2).
     # Omit its trace term -rho^2/(8eps0^3)*t^-D.
     Y=add(scale(mon(D),eps0),scale(rho,Q(1,2)/eps0),
           mul(a,mon(gamma)),mul(b,mon(gamma-1)),mul(tau,mon(-pN)))
     baseline=scale(power(rho,3),Q(1,2)/eps0**2)
    else:
     # W=z^3*eps0*(1+rho/eps0^4*z^-4)^(1/4);
     # its t^-pN coefficient is zero.
     Y=add(scale(mon(3*D),eps0),scale(mul(rho,mon(-D)),Q(1,4)/eps0**3),
           mul(a,mon(gamma)),mul(b,mon(gamma-1)),mul(tau,mon(-pN)))
     baseline={}
    P=mon(pN)
    F1=add(power(Y,p),scale(power(P,q0),-eps0**p))
    F2=add(power(F1,p),scale(mul(power(rho,p),power(P,q1)),-1))
    lam=scale(power(rho,p-1),p*p*eps0**(p-1))
    assert max(e[0] for e in F1)==q1N and coeff(F1,q1N)==rho
    assert max(e[0] for e in F2)==qI
    A=coeff(F2,qI);B=coeff(F2,qI-1);C=coeff(F2,r)
    assert A==mul(lam,a) and B==mul(lam,b)
    assert coeff(Y,-pN)==tau
    trace=mul(lam,tau)
    # Cross-multiplied identity avoids any numerical rho specialization:
    # C-baseline-lambda*tau = (p-1)/(p*rho^p)*A*B.
    assert scale(mul(power(rho,p),add(C,scale(baseline,-1),scale(trace,-1))),p)==scale(mul(A,B),p-1)
    assert B and add(C,scale(baseline,-1),scale(trace,-1))
    controls.append({"p":p,"epsilon0":str(eps0),"N":N,"D":D,
       "final_order":qI,"distinguished_exponent":r,"last_gap":N,
       "indeterminate_parameters":["rho!=0","a!=0","b","tau"],
       "baseline_nonzero":bool(baseline),"nonzero_B_retained":True,
       "AB_and_trace_identity":True})
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 receipt={"status":"PASS","claim":"Exceptional final Laurent coefficient has exact AB correction and optional trace correction",
  "environment":{"python":platform.python_version(),"arithmetic":"Independent Fraction Laurent polynomials in t,rho,a,b,tau"},
  "formula":"C_N-lambda_I*tau_N-(p-1)/(p*epsilon_last^p)*A_N*B_N=-lambda_I*[z^(-p/2)]W",
  "controls":controls,
  "scope":"All-N conditional theorem is written; finite scalar Laurent controls are not operator/DC-pair realizations.",
  "hashes":{"EXCEPTIONAL_LAURENT_CORRECTION.md":sha(HERE/"EXCEPTIONAL_LAURENT_CORRECTION.md"),
   "check_exceptional_laurent_correction.py":sha(Path(__file__)),
   "CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md":sha(HERE/"CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md"),
   "sources/zheglov-2410.06959v5.pdf":sha(HERE/"sources/zheglov-2410.06959v5.pdf")},
  "formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="exceptional-laurent-correction-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"exceptional_laurent_correction_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: 12 symbolic nonzero-B Laurent controls; AB and trace corrections, including nonzero baseline")
 print("Receipt:",dest)
if __name__=="__main__":main()
