"""Fresh triple residual and linear pole-chart source, exact characteristic zero."""
from pathlib import Path
from fractions import Fraction
from math import lcm,gcd
from functools import reduce
import json,hashlib,time,sympy as S
from flint import fmpz_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();s,a,b=S.symbols('s a b');ctx=fmpz_mpoly_ctx.get(['d','a','s','b']);d,an,sn,bn=ctx.gens()
def native(e):return ctx.from_dict({(0,ex[1],ex[0],ex[2]):int(c) for ex,c in S.Poly(e,s,a,b).terms()})
def coefficient(F,k,j):return ctx.from_dict({ex[:k]+(0,)+ex[k+1:]:int(c) for ex,c in F.to_dict().items() if ex[k]==j})
G2=a*a+a+s*b+s**3;G3=2*s*G2+a+s**3+s*s*b;G4=3*s*s*G2+b+s**4;X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian([s,a,b])));gn=list(map(native,[G2,G3,G4,X1]));X=[gn[3],gn[0]+sn*gn[3],gn[1]+sn**2*gn[3],gn[2]+sn**3*gn[3]];t=sn+d;Y1,Y2,Y3,Y4=X;bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4;at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt;R=at*at+at+t*bt+t**3+Y1*t-Y2
src=json.loads((HERE.parent/'gao-f7-quadruple/quadruple-input.json').read_text());f2=native(S.sympify(src['equations'][0]['polynomial']));c2=coefficient(R,0,2);ex=next(iter(f2.to_dict()));ratio=Fraction(int(c2.to_dict()[ex]),int(f2.to_dict()[ex]));assert ratio.denominator*c2==ratio.numerator*f2
Q=ctx.from_dict({(ex[0]-3,)+ex[1:]:int(c) for ex,c in R.to_dict().items() if ex[0]>=3});assert R==d**3*Q+d*d*c2;assert coefficient(Q,0,9)==1 and Q.degrees()[3]==2
A,B,C=[coefficient(f2,3,j) for j in [2,1,0]];D=sn*(3*sn**3-2);assert A==4*D*D
# On D=0 the quadratic term vanishes. The linear coefficient reduces to
# 2*(2a+1)*(6s^3-1), a unit times(2a+1) at every carrier value.
U=B*B*coefficient(Q,3,0)-B*C*coefficient(Q,3,1)+C*C*coefficient(Q,3,2)
def pole_reduce(F,cubic):
 out={}
 for (i,j,k,l),v in F.to_dict().items():
  assert l==0
  if cubic:
   power,rem=divmod(int(k),3);value=Fraction(int(v))*Fraction(2,3)**power;key=(int(i),int(j),rem)
  else:
   if k:continue
   value=Fraction(int(v));key=(int(i),int(j),0)
  out[key]=out.get(key,Fraction(0))+value
 out={k:v for k,v in out.items() if v};den=reduce(lcm,(v.denominator for v in out.values()),1);ints={ex:int(v*den) for ex,v in out.items()};content=reduce(gcd,(abs(v) for v in ints.values()),0);ints={ex:v//content for ex,v in ints.items()};return {'denominator_clearing':str(den),'integer_content_removed':str(content),'terms':[[list(ex),str(v)] for ex,v in sorted(ints.items())],'degrees':[max(ex[k] for ex in ints) for k in range(3)]}
charts=[{'name':'s-zero','base_modulus_coefficients':[0,1],**pole_reduce(U,False)},{'name':'s-cubic','base_modulus_coefficients':[-2,0,0,3],**pole_reduce(U,True)}]
rec={'status':'UNVERIFIED_SOURCE_CANDIDATE','variables':['d','a','s','b'],'f2_ratio_in_Taylor_coefficient':str(ratio),'Q9_terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(Q.to_dict().items())],'pole_charts':charts,'source_sha256':hashlib.sha256((HERE.parent/'gao-f7-quadruple/quadruple-input.json').read_bytes()).hexdigest(),'scope':'On D0=s*(3s^3-2)=0 and B!=0, chart polynomial is a nonzero rational scalar times B^2*Q9(d;s,a,-C/B). No exclusion follows from this generation.'}
(HERE/'pole-source-candidate.json').write_text(json.dumps(rec,separators=(',',':'))+'\n')
print(json.dumps({'Q9_degrees':list(map(int,Q.degrees())),'Q9_terms':len(Q.to_dict()),'pole_charts':[{k:v for k,v in row.items() if k!='terms'} for row in charts],'seconds':time.monotonic()-START}))
