#!/usr/bin/env python3
"""Exact PBW transport and conditional final-coefficient comparison controls."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,hashlib,json,platform,tempfile
HERE=Path(__file__).resolve().parent
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
 out={tuple(0 for _ in next(iter(p))):Q(1)}
 for _ in range(n):out=mul(out,p)
 return out
def pbw(p,q):
 out={}
 for (i,j),c in p.items():
  for (k,l),d in q.items():
   for r in range(min(j,k)+1):
    e=(i+k-r,j+l-r)
    out[e]=out.get(e,Q(0))+c*d*(comb(j,r)*factorial(k)//factorial(k-r))
 return clean(out)
def ppow(p,n):
 out={(0,0):Q(1)}
 for _ in range(n):out=pbw(out,p)
 return out
def shear(poly,N):
 out={}
 for (i,j),c in poly.items():out=add(out,scale(pbw(ppow({(1,0):Q(1),(0,N):Q(1)},i),{(0,j):Q(1)}),c))
 return out
def homogeneous(poly,order):
 return {i:c for (i,j),c in poly.items() if j-i==order}
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 controls=[]
 for a,b in [(2,1),(3,2),(8,4),(23,12)]:
  for N in [2,3,5]:
   rho=Q(-2);A=Q(7,5)
   P={(a,b):Q(1),(a-1,b):a*rho,(a-1,b-1):A,(a-2,b):Q(11)}
   image=shear(P,N);p=a*N+b
   assert max(j-i for i,j in image)==p
   assert all(not homogeneous(image,p-j) for j in range(1,N))
   assert homogeneous(image,p-N)=={0:a*rho}
   assert homogeneous(image,p-N-1)=={1:Q(a),0:A+N*comb(a,2)}
   controls.append({"a":a,"b":b,"N":N,"gap":N,
     "head_minus_N":str(a*rho),"raw_adjacent_constant":str(A+N*comb(a,2))})
 # Seven-variable polynomial identity (N,a,b,c,A,B,k0), using Fraction only.
 zero=(0,)*7;one={zero:Q(1)}
 def v(i):
  e=[0]*7;e[i]=1;return {tuple(e):Q(1)}
 N,a,b,c,A,B,k0=[v(i) for i in range(7)]
 minus=lambda x,y:add(x,scale(y,-1))
 p=add(mul(a,N),b)
 # Terminal e=(b(c-1)+a)/a; q=(a*c*N+b(c-1)+a)/a.
 qnum=add(mul(mul(a,c),N),mul(b,minus(c,one)),a)
 Delta=minus(a,b)
 # T=2*a*p*(conditional normalized final Step7 coefficient).
 T=add(mul(mul(a,minus(b,a)),add(p,scale(one,-1),scale(mul(p,k0),-2))),
       scale(power(qnum,2),-1),mul(mul(a,p),qnum),scale(mul(A,qnum),-2),
       scale(mul(mul(mul(N,a),minus(a,one)),qnum),-1),
       scale(mul(mul(a,p),B),2),mul(mul(mul(mul(a,N),c),minus(c,one)),p))
 Z=add(scale(mul(A,c),-2),scale(mul(B,a),2),
       scale(mul(power(a,2),k0),2),scale(power(a,2),-1),
       mul(mul(a,b),c),scale(mul(mul(a,b),k0),-2),mul(a,b),
       scale(mul(a,c),-2),a,scale(mul(b,power(c,2)),-1),
       scale(mul(b,c),2),scale(b,-1))
 residual=mul(Delta,add(mul(a,b),scale(A,-2)))
 assert not add(T,scale(mul(p,Z),-1),scale(residual,-1))
 assert all(e[0]<=1 for e in T)
 # Omitting PBW N-binomial terms leaves the spurious quadratic numerator.
 oldT=add(mul(mul(a,minus(b,a)),add(p,scale(one,-1),scale(mul(p,k0),-2))),
         scale(power(qnum,2),-1),mul(mul(a,p),qnum),scale(mul(A,qnum),-2),
         scale(mul(mul(a,p),B),2))
 oldN2={(0,)+e[1:]:q for e,q in oldT.items() if e[0]==2}
 assert oldN2==mul(mul(power(a,2),c),minus(a,c))
 # Rational numerical Schur controls; unequal centroids retained.
 comparisons=[]
 aa,bb,cc,ee=8,4,23,12
 assert aa*ee-bb*cc==aa-bb
 for Adj in [Q(aa*bb,2),Q(aa*bb,2)+1]:
  values=[]
  for NN in [3,5,9]:
   pp=aa*NN+bb;qq=cc*NN+ee;rawA=Adj+NN*comb(aa,2);rawB=Q(2,3)+NN*comb(cc,2)
   gam2=Q(-aa,2*pp);gam1=-rawA/pp+Q(aa,2)
   assert aa+2*pp*gam2==0
   assert rawA+pp*pp*gam2+pp*gam1==0
   assert cc+2*qq*gam2==Q(bb-aa,pp)
   val=(bb-aa)*(Q(1,2)-Q(1,2*pp)-Q(5,7))+gam2*qq*qq+gam1*qq+rawB
   const=Q(2,3)-Q(cc,aa)*Adj+(bb-aa)*(Q(1,2)-Q(5,7))+Q(bb*cc,2)-cc+Q(1,2)-Q(bb*(cc-1)**2,2*aa)
   assert val==const+Q((aa-bb)*(aa*bb-2*Adj),2*aa*pp)
   values.append(val)
  assert (len(set(values))==1)==(Adj==Q(aa*bb,2))
  comparisons.append({"A_source":str(Adj),"values":[str(q) for q in values]})
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 receipt={"status":"PASS","claim":"Exact tame-shear first components, head-N versus head-1 distinction, and conditional corrected Step7 comparison",
  "environment":{"python":platform.python_version(),"arithmetic":"Independent Fraction normal-PBW multiplication and seven-variable polynomial identity"},
  "PBW_controls":controls,"symbolic_final_identity":True,
  "residual":"(a-b)*(a*b-2*A_source)/(2*a*(a*N+b))",
  "vanishing_condition":"A_source=a*b/2 for a>b>0",
  "centroid_scope":"Unequal edge centroids affect head-N; monomial-face gap gives head-1 zero for N>1",
  "Schur_comparisons":comparisons,
  "hashes":{"PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md":sha(HERE/"PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md"),
   "check_pbw_transport_final_coefficient.py":sha(Path(__file__)),
   "EXCEPTIONAL_LAURENT_CORRECTION.md":sha(HERE/"EXCEPTIONAL_LAURENT_CORRECTION.md"),
   "sources/zheglov-2410.06959v5.pdf":sha(HERE/"sources/zheglov-2410.06959v5.pdf")},
  "scope":"First-component and Schur algebra plus conditional use of source body/tail/averaging formula; no actual-pair construction, full theorem refutation or DC1 conclusion.",
  "formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="pbw-transport-final-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"pbw_transport_final_coefficient_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: PBW N-binomial terms, exact residual identity and zero/nonzero subprincipal controls")
 print("Receipt:",dest)
if __name__=="__main__":main()
