from pathlib import Path
import sys
from functools import lru_cache
from math import comb, factorial
from fractions import Fraction as F
from collections import Counter
import sympy as S
root=Path('/home/royvanrijn/src/jacobian-research')
sys.path.insert(0,str(root/'research/jacobian-followup/moment-vanishing-audit'))
from verify_exceptional_sl2 import group
x=S.symbols('x')
@lru_cache(None)
def derivative(l,q):
 return S.Poly(S.diff(S.legendre(l,x),x,q),x,domain=S.QQ)
def expand(k,j):
 beta=S.Poly(sum((-1)**r*comb(2*k,r)*(comb(k+2*r,j+k) if 0<=j+k<=k+2*r else 0)*x**(2*r) for r in range(2*k+1)),x,domain=S.QQ)
 if j<0:
  beta,remainder=beta.div(S.Poly((1-x*x)**(-2*j),x))
  assert remainder.is_zero
 q=2*abs(j); result={}
 while not beta.is_zero:
  l=beta.degree()+q
  c=beta.LC()/derivative(l,q).LC()
  result[l]=c
  beta-=c*derivative(l,q)
 return result
def bvector(k):
 coeff={j:expand(k,j) for j in range(-k,k+1)}
 b={}
 for l in range(0,4*k+1,2):
  b[l]=(coeff[0].get(l,0)**2+2*sum(coeff[j].get(l,0)*coeff[-j].get(l,0)*F(factorial(l+2*j),factorial(l-2*j)) for j in range(1,min(k,l//2)+1)))/(2*l+1)
 assert sum(b.values())==0,(k,b)
 return b
def wvector(vals,d,L):
 def add(a,b):return(a[0]+b[0],a[1]+b[1])
 def mul(a,b):return(a[0]*b[0]+d*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
 def scale(a,c):return(a[0]*c,a[1]*c)
 tau=Counter(add(add(mul(a,a),mul(b,b)),scale(add(mul(c,c),mul(e,e)),-1)) for a,b,c,e in vals)
 out=[(F(0),F(0))]*(L+1)
 for t,count in tau.items():
  old=(F(0),F(0)); now=(F(1),F(0))
  for l in range(L+1):
   out[l]=add(out[l],scale(now,F(count,len(vals))))
   new=scale(add(scale(mul(t,now),2*l+1),scale(old,-l)),F(1,l+1))
   old,now=now,new
 return out
O=group(2,[[(0,F(1,2)),(0,F(1,2)),(0,0),(0,0)]])[0]
I=group(5,[[(0,0),(F(1,2),0),(F(1,4),-F(1,4)),(-F(1,4),-F(1,4))]])[0]
ws={'2O':wvector(O,2,24),'2I':wvector(I,5,24)}
for k in range(1,7):
 b=bvector(k)
 print('k',k,'b',b,flush=True)
 for name,w in ws.items():
  result=tuple(sum(F(bv)*w[l][i] for l,bv in b.items()) for i in (0,1))
  print(name,[str(v) for v in result],flush=True)
