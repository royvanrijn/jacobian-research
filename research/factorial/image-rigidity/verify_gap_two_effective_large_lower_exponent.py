#!/usr/bin/env python3
"""Exact numerical/algebraic budgets for the complete gap-two endpoint proof.

This guard supports a written proof importing separately audited analytic
contours. It does not certify those contours by sampling moments or roots.
"""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: exact guards must execute.')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse,json,math
HERE=Path(__file__).resolve().parent
ROWS=[]
def ck(name,left,right,op='<'):
 left,right=Q(left),Q(right)
 assert(left<right if op=='<'else left==right),(name,left,op,right)
 ROWS.append(dict(name=name,left=str(left),relation=op,right=str(right)))
def add(a,b):
 out=dict(a)
 for k,v in b.items():out[k]=out.get(k,Q(0))+v
 return{k:v for k,v in out.items()if v}
def scale(a,c):return{k:v*c for k,v in a.items()if v*c}
def mul(a,b):
 out={}
 for i,u in a.items():
  for j,v in b.items():
   k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,Q(0))+u*v
 return{k:v for k,v in out.items()if v}
def var(n,i):return{tuple(int(i==j)for j in range(n)):Q(1)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'gap_two_effective_large_lower_exponent_verification.json');args=ap.parse_args()
 m0=10**16;amax=4*10**16;delta=Q(1,10**8);r=Q(41,100);R=Q(17,40)
 ck('actual threshold delta',4,amax*delta**2,'=')
 ck('soft boundary exact',Q(4),2*2,'=')
 ck('m cube-root lower',200000**3,m0)
 ck('m two-thirds lower',200000**2,4*10**10,'=')
 ck('true enlarged h',Q(47,50)/(Q(9,20)-R),38)
 ck('true annular G',R*38/((R-r)*(R*R-Q(1,8))),20000)
 ck('true annular G derivative',R*38/((R-r)**2*(R*R-Q(1,8))),1300000)
 ck('compact critical-reference comparison',450*36000,2*10**7)
 ck('large critical-reference comparison',Q(16,3),2*10**7)
 ck('right circle connector loss',Q(1,40000),Q(1,20000)-Q(41,2500000))
 ck('left circle connector loss',Q(1,10000),Q(1,7000)-Q(41,2500000))
 ck('compact circle connector loss',Q(1,40000),Q(1,100)-Q(41,2500000))
 ck('original circle endpoint loss',Q(1,40000),Q(66,100000)-Q(41,2500000))
 ck('left circle real derivative',4*Q(8,25)**2-r*r-Q(1,8),Q(233,2000),'=')
 # pi<22/7 only used for a simple upper bound on connector lengths.
 connector=Q(62,1000)+3*Q(22,7)*r
 ck('all added connector lengths',connector,10)
 rho=Q(3,8)
 ck('weighted negative tail elementary coefficient',365/(3*rho)+365**2/(9*rho**4),10**6)
 ck('weighted negative compact coefficient',4*10**6*450*36000,10**15)
 ck('weighted negative large coefficient',4*10**6*Q(16,3),10**15)
 ck('positive weighted tail coefficient',2*10**10,10**15)
 ck('old compact weighted tail coefficient',2*10**6,10**15)
 ck('positive and compact tail exponential',Q(1,100000),Q(1,40000))
 ck('negative tail stronger exponential',Q(1,40000),Q(1,14000))
 # Every listed term has already been normalized by 2pi*m^-1/3.
 ck('canonical A/B tail rounded coefficient',Q(3*10**15,6),5*10**14,'=')
 ck('A/B connector coefficient',Q(3*10*2*10**7,6),10**9)
 ck('G derivative connector coefficient',Q(1300000*10*2*10**7,6),5*10**13)
 ck('two physical G endpoint coefficient',Q(2*20000*2*10**7,6),2*10**11)
 total=60000+5*10**14+10**9+5*10**13+2*10**11
 ck('aggregate other coefficient',total,10**16)
 ck('actual original nonlocal exponent',Q(1,40000),Q(33,50000))
 ck('sixth Taylor exponential budget',Q(720*10**46,m0**4),1)
 ck('local compact complete remainder',Q(1300000*2*10**7,6)+1,10**14)
 ck('local large complete remainder',Q(1300000*300,6)+1,10**9)
 Cc,Cl=10**14,10**9
 compact_lambda=(Q(52,100)/200000+Q(Cc,m0))/(Q(98,100)-Q(Cc,m0))
 large_K=(Q(52,100)+Q(Cl,4*10**10))/(Q(98,100)-Q(Cl,m0))
 compact_K=(Q(52,100)+Q(Cc,4*10**10))/(Q(98,100)-Q(Cc,m0))
 ck('compact actual lambda',compact_lambda,Q(11,1000))
 ck('large rescaled K',large_K,Q(56,100))
 ck('uniform rescaled K',compact_K,2600)
 ck('compact Robin product',4*Q(11,1000)**2,Q(1,4))
 ck('large Robin product',Q(56,100)**2/8,Q(1,4))
 ck('compact Robin cubic',Q(11,1000)**3,Q(1,16))
 ck('large Robin cubic',(Q(56,100)/200000)**3,Q(1,16))
 ck('Robin path bootstrap',Q(1,4)+2*Q(1,16),Q(1,2))
 ck('Robin path reciprocal',Q(8,5),2)
 ck('radial geometric lower',Q(1,2),Q(3,4)-Q(15,4)*delta)
 ck('endpoint window real margin',Q(1,10**7)+3*delta*(1+3*delta),Q(1,10**5))
 ck('endpoint window imaginary margin',Q(21,10)*(1+3*delta**2),Q(21001,10000))
 ck('actual adjacent-order drift versus two bands',6*2600,10**8)
 # Universal identities: all indeterminates are formal, not sampled.
 D,V,Z=var(3,0),var(3,1),var(3,2);one={(0,0,0):Q(1)}
 # (1+D*V/2)^2+(-1+D*Z)=D*(Z+V+D*V^2/4).
 lhs=add(mul(add(one,scale(mul(D,V),Q(1,2))),add(one,scale(mul(D,V),Q(1,2)))),add(scale(one,-1),mul(D,Z)))
 rhs=mul(D,add(add(Z,V),scale(mul(D,mul(V,V)),Q(1,4))))
 assert lhs==rhs
 # Exact radial order-shift after clearing n^2*delta.
 N,n,Z,D=[var(4,j)for j in range(4)]
 lhs=add(mul(mul(n,n),mul(D,Z)),mul(add(mul(N,N),scale(mul(n,n),-1)),add(mul(D,Z),{(0,0,0,0):Q(-1)})))
 rhs=add(mul(mul(N,N),mul(D,Z)),add(mul(n,n),scale(mul(N,N),-1)))
 assert lhs==rhs
 # Cubic radial tail lower affine tangent remainder is q^2*(3rho+q).
 q,rho=var(2,0),var(2,1)
 cube=mul(add(q,rho),mul(add(q,rho),add(q,rho)))
 tangent=add(mul(rho,mul(rho,rho)),scale(mul(mul(rho,rho),q),3))
 assert add(cube,scale(tangent,-1))==mul(mul(q,q),add(scale(rho,3),q))
 # Actual scalar normalization: literal small integer coefficients are only
 # controls for the universal integral expansion, not an all-order inference.
 controls=[]
 for a,m in[(4,1),(4,2),(6,3),(8,4)]:
  actual=[Q(math.comb(m,k)*math.factorial(a*m+2*k),math.factorial(a*m))for k in range(m+1)]
  rising=[]
  for k in range(m+1):
   product=math.prod(range(a*m+1,a*m+2*k+1));rising.append(Q(math.comb(m,k)*product))
  assert actual==rising;controls.append(dict(a=a,m=m,coefficients=[str(v)for v in actual]))
 # Three deliberately false alternatives must fail.
 assert compact_K**2/8>Q(1,4) # The compact C cannot be used globally.
 assert compact_lambda>Q(1,1000) # A falsely tiny compact Robin number fails.
 wrong=[math.comb(2,k)*math.prod(range(4+1,4+2*k+1))for k in range(3)]
 correct=[math.comb(2,k)*math.prod(range(8+1,8+2*k+1))for k in range(3)]
 assert wrong!=correct and wrong[1]==60 and correct[1]==180
 proof=HERE/'GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md'
 result=dict(schema=1,status='PASS_EXACT_BUDGETS_ONLY',scope='Complete-period aggregation and Robin arithmetic supporting the written theorem and separately proved contour inputs. No all-order conclusion follows from these finite guards alone.',proof_sha256=sha256(proof.read_bytes()).hexdigest(),checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),rational_guards=ROWS,universal_identities=3,actual_moment_normalization_controls=controls,negative_controls=3)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='PASS',rational_guards=len(ROWS),universal_identities=3,normalization_controls=len(controls),negative_controls=3)))
if __name__=='__main__':main()
