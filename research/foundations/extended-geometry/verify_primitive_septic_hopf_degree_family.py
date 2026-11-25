#!/usr/bin/env python3
"""Finite exact checks; all-order conclusions are proved in the companion note."""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as Q
from math import factorial,comb
import argparse,hashlib,json,platform

def add(P,Qp):
 out=dict(P)
 for e,c in Qp.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def mul(P,Qp):
 out=defaultdict(int)
 for a,c in P.items():
  for b,d in Qp.items():out[tuple(x+y for x,y in zip(a,b))]+=c*d
 return {e:c for e,c in out.items() if c}
def scale(P,c):return {e:c*v for e,v in P.items() if c*v}
def power(P,n):
 out={(0,)*len(next(iter(P))):1}
 for _ in range(n):out=mul(out,P)
 return out
R={(1,0,1,0):1,(0,1,0,1):1}
Z={(1,0,0,1):1}
W={(0,1,1,0):2}
T={(1,0,1,0):1,(0,1,0,1):-1}
L=add(R,Z)
D=add(scale(mul(W,power(R,2)),2),scale(mul(add(scale(R,2),Z),power(T,2)),-1))
E=add(scale(mul(W,R),2),scale(power(T,2),-1))
F4=mul(L,D);ww=mul(T,L);gg=mul(L,E);shift=add(power(R,2),ww);F7=mul(F4,gg)

def original_operator(P):
 # Direct two-variable derivatives for balanced input: all surviving outputs are constants.
 assert all(e[0]+e[1]==e[2]+e[3] for e in P), "balanced input required"
 return sum(c*factorial(e[2])*factorial(e[3]) for e,c in P.items() if e[0]==e[2] and e[1]==e[3])
def signed_moment(m,a):
 h=a*m
 if h%2:return Q(0)
 return sum((Q((-1)**j*comb(m,j),h+2*j+1) for j in range(m+1)),Q(0))
def shifted_moment(m,a):return Q(2**((a+2)*m)*factorial((a+1)*m)*factorial(m),factorial((a+2)*m+1))
def evaluate(P,z):
 result=0
 for e,c in P.items():
  v=c
  for i in range(len(e)):v*=z[i]**e[i]
  result+=v
 return result

def specialization(P):
 out=defaultdict(int)
 for e,c in P.items():
  if e[3]==0:out[e[1]]+=c
 return {e:c for e,c in out.items() if c}

def main():
 ap=argparse.ArgumentParser(description='Exact finite checks of the primitive septic and joint Hopf degree family.')
 ap.add_argument('--output',type=Path)
 args=ap.parse_args()
 assert mul(R,F7)==add(power(F4,2),mul(F4,power(ww,2)))
 checks=[]
 for variant in ['signed','shifted']:
  for a,b,mmax in [(0,1,5),(0,0,4),(1,0,4),(2,0,3),(1,1,3),(0,2,3)]:
   F=mul(mul(F4,power(ww if variant=='signed' else shift,a)),power(gg,b));degree=4+2*a+3*b
   assert {(e[0]+e[1],e[2]+e[3]) for e in F}=={(degree,degree)}
   expected_radial=-4*(-2)**a*4**b
   assert evaluate(F,[1,1,1,-1])==expected_radial
   # In the specified ORIGINAL specialization, W2=1, B2=2, G3=4u-1.
   spec=specialization(F)
   univar=lambda xs:{(k,):v for k,v in xs.items()}
   want=mul(univar({0:-2,1:4}),power(univar({0:-1,1:4}),b));want=scale(want,2**a if variant=='shifted' else 1)
   assert {(k,):v for k,v in spec.items()}==want
   P={(0,0,0,0):1}
   for m in range(1,mmax+1):
    P=mul(P,F);pure=original_operator(P);marked=original_operator(mul(Z,P))
    angular=signed_moment(m,a) if variant=='signed' else shifted_moment(m,a)
    expected=Q(factorial(degree*m+2))*angular
    assert pure==0,(variant,a,b,m,pure)
    assert marked==expected,(variant,a,b,m,marked,expected)
    checks.append({'variant':variant,'a':a,'b':b,'degree':degree,'m':m,'pure':str(pure),'Z_marked':str(marked),'expected':str(expected),'monomials':len(P)})
 # Separate literal angular multiplication, not using ORIGINAL polynomials or the all-order endpoint law.
 U={(0,0):1,(1,0):1};x={(1,0):1};t={(0,1):1};y={(-1,0):Q(1,2),(-1,2):Q(-1,2)}
 e=mul(U,add(scale(y,2),scale(mul(add({(0,0):2},x),power(t,2)),-1)))
 g=mul(U,add(scale(y,2),scale(power(t,2),-1)));f=mul(e,g)
 integral=lambda P:sum((Q(c,h+1) for (z,h),c in P.items() if z==0 and h%2==0),Q(0))
 P={(0,0):1};angular_checks=[]
 for m in range(1,9):
  P=mul(P,f);pure=integral(P);marked=integral(mul(x,P));assert pure==0;assert marked==shifted_moment(m,0)
  angular_checks.append({'m':m,'pure':str(pure),'x_marked':str(marked)})
 assert evaluate(F7,[1,1,1,-1])==-16
 assert specialization(F7)=={2:16,1:-12,0:2}
 # One literal corruption must destroy the first pure moment.
 assert original_operator(add(F7,power(R,7)))==factorial(8)!=0
 record={'schema':1,'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Fresh literal ORIGINAL-variable contractions and separate literal Laurent controls; written endpoint/binomial proof supplies every unrestricted assertion. No finite-prefix extrapolation.','original_variable_order':['xi1','xi2','z1','z2'],'F7_terms':[{'powers':list(e),'coefficient':str(c)} for e,c in sorted(F7.items())],'original_checks':checks,'angular_checks':angular_checks,'radial_zero_value':-16,'specialization_lowfirst':[2,-12,16],'corrupted_F7_plus_R7_first_pure_moment':str(factorial(8)),'canonical_writes':False,'formal_verification':False,'independent_replay':False,'external_review':False,'actual_DAG_cases':0}
 record['python_version']=platform.python_version()
 record.pop('canonical_writes',None)
 record.pop('actual_DAG_cases',None)
 if args.output:args.output.write_text(json.dumps(record,indent=2)+'\n')
 print('PASS: 44 original-variable contractions, eight Laurent controls, radial and simple-root certificates, and one rejected corruption.')

if __name__=='__main__':main()
