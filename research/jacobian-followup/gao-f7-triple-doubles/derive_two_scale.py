"""Correct two-scale r!=0 reduction, retaining all grading and pole cases.
The source has weight layers differing by3 and7. Use x=r^-3,y=r^-7,
y^3=x^7; r=x^2/y. Every normalized coefficient is linear in y, and
coefficient5 has a constant y coefficient, so y can be eliminated globally.
"""
from pathlib import Path
from fractions import Fraction
from math import gcd,lcm
from functools import reduce
import json,hashlib,time
from flint import fmpz_mpoly_ctx,fmpq_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();data=json.loads((HERE/'delta-zero-fivevar-input.json').read_text());ctx=fmpz_mpoly_ctx.get(['u','v','z','w','x','y']);u,v,z,w,x,y=ctx.gens();scaled={};degrees=[]
for row in data['equations']:
 wt=[sum(int(e)*b for e,b in zip(ex,[1,1,2,1,3])) for ex,c in row['terms']];D=max(wt);terms={}
 for (ex,c),weight in zip(row['terms'],wt):
  loss=D-weight
  if loss%3==0:xp,yp=loss//3,0
  else:assert loss>=7 and (loss-7)%3==0;xp,yp=(loss-7)//3,1
  key=(ex[1],ex[2],ex[3],ex[4],xp,yp);terms[key]=terms.get(key,0)+int(c)
 scaled[row['old_coefficient']]=ctx.from_dict(terms);degrees.append(D)
assert degrees==[12,11,10,9,8,7,5]
def yc(F,j):return ctx.from_dict({ex[:5]+(0,):int(c) for ex,c in F.to_dict().items() if ex[5]==j})
K=yc(scaled[5],1);assert len(K.to_dict())==1 and next(iter(K.to_dict()))==(0,0,0,0,0,0);ki=int(next(iter(K.to_dict().values())));T=yc(scaled[5],0);assert ki!=0
small=fmpz_mpoly_ctx.get(['u','v','z','w','x']);us,vs,zs,ws,xs=small.gens()
def primitive(F):
 td={ex:Fraction(str(c)) for ex,c in F.to_dict().items() if c};den=reduce(lcm,(c.denominator for c in td.values()),1);out={ex:int(c*den) for ex,c in td.items()};g=reduce(gcd,(abs(c) for c in out.values()),0);return small.from_dict({ex:int(c//g) for ex,c in out.items()}),{'denominator_clearing':str(den),'integer_content_removed':str(g)}
remaining={};ynorm={}
for index,F in scaled.items():
 if index==5:continue
 assert F.degrees()[5]<=1
 N=ki*yc(F,0)-yc(F,1)*T;NN=small.from_dict({ex[:5]:int(c) for ex,c in N.to_dict().items()});remaining[index],ynorm[index]=primitive(NN)
Tsmall=small.from_dict({ex[:5]:int(c) for ex,c in T.to_dict().items()});toric=Tsmall**3+ki**3*xs**7;toric,toricnorm=primitive(toric)
def xc(F,j):return small.from_dict({ex[:4]+(0,):int(c) for ex,c in F.to_dict().items() if ex[4]==j})
E7=remaining[7];P5,M=xc(E7,0),xc(E7,1);H=3+18*us-10*us*zs-12*vs+zs*zs;ratio,rem=divmod(M,H);assert not rem and ratio and all(sum(ex)==0 for ex in ratio.to_dict());assert E7.degrees()[4]==1
normalized_source={'status':'UNVERIFIED_NORMALIZED_WORKING_INPUT','variables':['u','v','z','w','x'],'source_sha256':hashlib.sha256((HERE/'delta-zero-fivevar-input.json').read_bytes()).hexdigest(),'weights':[1,2,1,3,3],'y_coefficient_constant':ki,'y_constant_terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(Tsmall.to_dict().items())],'equations':[{'source_coefficient':i,'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(F.to_dict().items())],'polynomial':str(F)} for i,F in remaining.items()],'toric_relation_terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(toric.to_dict().items())],'nonzero_conditions':['x'],'boundary':'The toric relation is retained as data but omitted from the necessary coefficient ideal for a stronger carrier-power target. Literal x^k membership would exclude the whole delta-zero case since r!=0 by GAOF7TRP1.'}
(HERE/'two-scale-normalized-input.json').write_text(json.dumps(normalized_source,separators=(',',':'))+'\n')
qc=fmpq_mpoly_ctx.get(['u','z','w','x']);uq,zq,wq,xq=qc.gens();V0=(3+18*uq-10*uq*zq+zq*zq)/12
for branch in ['pole','generic']:
 equations=[]
 for index,F in list(remaining.items())+[(99,toric)]:
  if branch=='generic' and index==7:continue
  if branch=='pole':
   out=qc.constant(0)
   for ex,c in F.to_dict().items():out+=qc.from_dict({(ex[0],ex[2],ex[3],ex[4]):int(c)})*V0**ex[1]
   td={ex:Fraction(str(c)) for ex,c in out.to_dict().items()};den=reduce(lcm,(c.denominator for c in td.values()),1);integers={ex:int(c*den) for ex,c in td.items()};g=reduce(gcd,(abs(c) for c in integers.values()),0);ints={ex:c//g for ex,c in integers.items()};norm={'denominator_clearing':str(den),'integer_content_removed':str(g)};variables=['u','z','w','x']
  else:
   degree=int(F.degrees()[4]);out=sum((xc(F,j)*(-P5)**j*M**(degree-j) for j in range(degree+1)),small.constant(0));td={ex[:4]:int(c) for ex,c in out.to_dict().items()};g=reduce(gcd,(abs(c) for c in td.values()),0);ints={ex:c//g for ex,c in td.items()};norm={'integer_content_removed':str(g),'x_degree_cleared':degree};variables=['u','v','z','w']
  four=fmpz_mpoly_ctx.get(variables);poly=four.from_dict(ints);equations.append({'source_coefficient':index,'normalization':norm,'degrees':list(map(int,poly.degrees())),'total_degree':int(max(sum(ex) for ex in ints)),'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(ints.items())],'polynomial':str(poly)})
 result={'status':'UNVERIFIED_WORKING_INPUT','branch':branch,'variables':variables,'source_sha256':hashlib.sha256((HERE/'delta-zero-fivevar-input.json').read_bytes()).hexdigest(),'weights':[1,1,2,1,3],'weighted_degrees':degrees,'scale_change':'u_old=r*u,v_old=r^2*v,z_old=r*z,w_old=r^3*w,x=r^-3,y=r^-7,y^3=x^7,r=x^2/y','y_coefficient_constant':ki,'y_constant_terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(Tsmall.to_dict().items())],'y_substitution_normalizations':ynorm,'toric_relation_normalization':toricnorm,'H':'3+18*u-10*u*z-12*v+z^2','E7_M_terms':[[list(map(int,ex[:3])),str(c)] for ex,c in sorted(M.to_dict().items())],'E7_P_terms':[[list(map(int,ex[:3])),str(c)] for ex,c in sorted(P5.to_dict().items())],'equations':equations,'nonzero_conditions':['x'] if branch=='pole' else ['H','P5'],'boundary':'Pole H0 uses constant-coefficient V recovery; generic H!=0 uses x=-P5/M and requires P5!=0. Both retain the toric scale relation. No global exclusion without exact replay.'}
 (HERE/('two-scale-'+branch+'-input.json')).write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({'branch':branch,'y_coefficient':ki,'equations':len(equations),'terms':[len(e['terms']) for e in equations],'degrees':[e['total_degree'] for e in equations],'seconds':time.monotonic()-START}),flush=True)
(HERE/'one-scale-failure.json').write_text(json.dumps({'outcome':'rejected_by_exact_weight_check','script':'derive_scale_delta_zero.py','reason':'The tb(t) terms introduce a weight loss7, so r^-3 alone does not parametrize the normalized equations. No one-scale input was written. The corrected model retains x=r^-3,y=r^-7 and y^3=x^7.'},indent=2)+'\n')
