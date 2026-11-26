#!/usr/bin/env python3
"""Supplementary finite native/deformed polygon controls; not a universal-prefix proof."""
if not __debug__:raise SystemExit('Optimized execution is refused; use ordinary Python.')
import argparse
from pathlib import Path
from math import comb,factorial
from fractions import Fraction
import json,hashlib
def vp(n,p):
 if n==0:return None
 n=abs(n);a=0
 while n%p==0:n//=p;a+=1
 return a
def hull(points):
 out=[]
 for pt in points:
  while len(out)>1:
   a,b=out[-2];c,d=out[-1];e,f=pt
   if (d-b)*(e-c)>=(f-d)*(c-a):out.pop()
   else:break
  out.append(pt)
 return out
def positive(h):return [(a,b) for a,b in zip(h,h[1:]) if b[1]>a[1]]
fac=[factorial(j) for j in range(121)]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',type=Path,help='Write a fresh supplementary receipt; never overwrite an existing file.')
args=ap.parse_args()
cases=0;coefficients=0;generic=0;divisible=0;nonzero_constant_exceptions=[]
for m in range(1,61):
 for C in [-7,-2,-1,1,2,7]:
  actual=[comb(m,i)*sum(comb(m-i,l)*fac[2*i+l]*C**(m-i-l) for l in range(m-i+1)) for i in range(m+1)]
  for p in [3,5,7,11,13]:
   if C%p==0:continue
   native=[vp(comb(m,i)*fac[2*i],p) for i in range(m+1)]
   HV=hull(list(enumerate(native)));AV=hull([(i,vp(c,p)) for i,c in enumerate(actual) if c])
   r=min(m%p,(p-1)//2)
   assert positive(HV)==positive(AV),(m,C,p,HV,AV)
   assert all(vp(actual[i],p) is None or vp(actual[i],p)>=h for i,h in enumerate(native))
   assert vp(actual[r],p)==0
   assert next(i for i,c in enumerate(actual) if c)<=r
   for (i,bi),(j,bj) in positive(HV):
    assert (i==m or (m-i)%p==0 or (2*i+1)%p==0)
    assert vp(actual[i],p)==bi and vp(actual[j],p)==bj
   if m%p==0:assert AV==HV and vp(actual[0],p)==0;divisible+=1
   else:generic+=1
   cases+=1;coefficients+=m+1
# Independent small exact witnesses for the required qualifications.
m=15;p=3;H=hull([(i,vp(comb(m,i)*fac[2*i],p)) for i in range(m+1)])
assert H==[(0,0),(3,2),(12,10),(13,11),(15,14)] and 13%p!=m%p and (2*13+1)%p==0
F1=[comb(1,i)*sum(comb(1-i,l)*fac[2*i+l]*(-1)**(1-i-l) for l in range(2-i)) for i in range(2)]
assert F1==[0,2] # F_1(A;-1)=2A, so its constant coefficient really vanishes.
assert hull([(1,0)])!=hull([(0,0),(1,0)])
H6=hull([(i,vp(comb(6,i)*fac[2*i],3)) for i in range(7)])
den=[Fraction(b[1]-a[1],b[0]-a[0]).denominator for a,b in zip(H6,H6[1:])]
assert den==[3,1]
out={'status':'PASS','scope':'Supplementary finite exact native and full deformed coefficient controls only; universal result is the written monotone-height and vertex-unit proof.','m_range':[1,60],'C_values':[-7,-2,-1,1,2,7],'p_values':[3,5,7,11,13],'p_unit_cases':cases,'p_divides_m_cases':divisible,'generic_initial_block_cases':generic,'coefficient_values_checked':coefficients,'qualification_controls':[{'false_claim':'all positive vertices i=m modp','witness':{'m':15,'p':3,'hull':H}},{'false_claim':'genericpnotdividingm has full equality including initialblock','witness':'m1,Cminus1,p3,F1=2A'}, {'false_claim':'arbitrarycarryingnativefaces all have pdivisible degree','witness':{'m':6,'p':3,'hull':H6,'face_denominators':den}}],'actual_engine_cases':0,'canonical_mutations':False,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
encoded=json.dumps(out,indent=2)+'\n'
if args.output:
 with args.output.open('x',encoding='utf-8') as f:f.write(encoded)
print(encoded,end='')