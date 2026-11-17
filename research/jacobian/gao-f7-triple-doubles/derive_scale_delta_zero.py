"""Exact r!=0 scale quotient and a complete denominator split of delta=0.
Weights(r,u,v,z,w)=(1,1,2,1,3); all equation layers differ by3.
Set U=u/r,V=v/r^2,Z=z/r,W=w/r^3,x=r^-3. The pure equation
is P5(U,V,Z)+32*x*H(U,V,Z)=0, H=3+18U-10UZ-12V+Z^2.
"""
from pathlib import Path
from fractions import Fraction
from math import gcd,lcm
from functools import reduce
import json,hashlib,time
from flint import fmpz_mpoly_ctx,fmpq_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();data=json.loads((HERE/'delta-zero-fivevar-input.json').read_text());ctx=fmpz_mpoly_ctx.get(['u','v','z','w','x']);u,v,z,w,x=ctx.gens();normalized=[];degrees=[]
for row in data['equations']:
 wt=[sum(int(e)*b for e,b in zip(ex,[1,1,2,1,3])) for ex,c in row['terms']];D=max(wt);assert all((D-j)%3==0 for j in wt);terms={}
 for (ex,c),weight in zip(row['terms'],wt):
  key=(ex[1],ex[2],ex[3],ex[4],(D-weight)//3);terms[key]=terms.get(key,0)+int(c)
 normalized.append(ctx.from_dict(terms));degrees.append(D)
assert degrees==[12,11,10,9,8,7,5]
def xc(F,k):return ctx.from_dict({ex[:4]+(0,):int(c) for ex,c in F.to_dict().items() if ex[4]==k})
E7=normalized[-1];P5=xc(E7,0);M=xc(E7,1);H=3+18*u-10*u*z-12*v+z*z;assert E7.degrees()[4]==1 and M==32*H
# H=0: V has a constant nonzero coefficient, so this chart loses no cases.
qc=fmpq_mpoly_ctx.get(['u','z','w','x']);uq,zq,wq,xq=qc.gens();V0=(3+18*uq-10*uq*zq+zq*zq)/12

def primitive(terms):
 terms={ex:Fraction(str(c)) for ex,c in terms.items() if c};den=reduce(lcm,(c.denominator for c in terms.values()),1);ints={ex:int(c*den) for ex,c in terms.items()};content=reduce(gcd,(abs(c) for c in ints.values()),0);return {ex:c//content for ex,c in ints.items()},{'denominator_clearing':str(den),'integer_content_removed':str(content)}
for branch in ['H-zero','H-nonzero']:
 equations=[]
 for row,F in zip(data['equations'],normalized if branch=='H-zero' else normalized[:-1]):
  if branch=='H-zero':
   out=qc.constant(0)
   for ex,c in F.to_dict().items():out+=qc.from_dict({(ex[0],ex[2],ex[3],ex[4]):int(c)})*V0**ex[1]
   td,norm=primitive(out.to_dict());variables=['u','z','w','x']
  else:
   degree=int(F.degrees()[4]);out=sum((xc(F,j)*(-P5)**j*M**(degree-j) for j in range(degree+1)),ctx.constant(0));assert out.degrees()[4]<=0;td,norm=primitive({ex[:4]:int(c) for ex,c in out.to_dict().items()});variables=['u','v','z','w']
  four=fmpz_mpoly_ctx.get(variables);out=four.from_dict(td);equations.append({'source_coefficient':row['old_coefficient'],'normalization':norm,'degrees':list(map(int,out.degrees())),'total_degree':int(max(sum(ex) for ex in td)),'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(td.items())],'polynomial':str(out)})
 target=qc.gen(3) if branch=='H-zero' else None
 rec={'status':'UNVERIFIED_WORKING_INPUT','branch':branch,'variables':variables,'source_sha256':hashlib.sha256((HERE/'delta-zero-fivevar-input.json').read_bytes()).hexdigest(),'weights':[1,1,2,1,3],'weighted_degrees':degrees,'scale_recovery':'u_old=r*u,v_old=r^2*v,z_old=r*z,w_old=r^3*w,x=r^-3; r!=0 covered by GAOF7TRP1 pole exclusion','H':'3+18*u-10*u*z-12*v+z^2','P5_terms':[[list(map(int,ex[:3])),str(c)] for ex,c in sorted(P5.to_dict().items())],'equations':equations,'nonzero_conditions':['x'] if branch=='H-zero' else ['H','P5'],'target_hint':'x^k' if branch=='H-zero' else '(H*P5)^k','boundary':'x=0 is an introduced infinity component, not an original source solution. A literal carrier-power identity is sufficient only after the stated branch and pole coverage.'}
 (HERE/(branch+'-scale-input.json')).write_text(json.dumps(rec,separators=(',',':'))+'\n');print(json.dumps({'branch':branch,'variables':variables,'equations':len(equations),'terms':[len(e['terms']) for e in equations],'degrees':[e['total_degree'] for e in equations],'seconds':time.monotonic()-START}),flush=True)
