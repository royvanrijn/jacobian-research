"""Fresh source transport and two complete five-variable coefficient branches.
Untrusted working input, no global exclusion. The old six-variable polynomials
are verified from the eliminant before any substitution is used.
"""
from pathlib import Path
from fractions import Fraction
from math import gcd,lcm
from functools import reduce
import hashlib,json,time,sympy as S
from flint import fmpq_mpoly_ctx,fmpz_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();oldpath=HERE.parent/'gao-f7-three-point/linear_cube_cubic_square-input.json';old=json.loads(oldpath.read_text());C=fmpq_mpoly_ctx.get(['t','r','u','v','w','z','h','k']);t,r,u,v,w,z,h,k=C.gens()
def co(F,j,n):return C.from_dict({ex[:j]+(0,)+ex[j+1:]:c for ex,c in F.to_dict().items() if ex[j]==n})
P=(t-r)**3*(t**3+u*t*t+v*t+w)**2*(t**3+z*t*t+h*t+k);Y1=-co(P,0,11)/4;Y2=(co(P,0,10)-4*Y1*Y1)/6;pc9=co(P,0,9);assert co(pc9,7,1)==1;kv=-12*Y1*Y2-2-co(pc9,7,0);P=co(P,7,0)+co(P,7,1)*kv;Y4=(6*Y1+9*Y2*Y2-co(P,0,8))/2;Y3=(co(P,0,6)-14*Y1*Y2+6*Y2*Y4-2)/2;bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4;at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt;R=at*at+at+t*bt+t**3+Y1*t-Y2;difference=R-P
for j in [12,11,10,9,8,6]:assert co(difference,0,j)==0
rr,uu,vv,ww,zz,hh=S.symbols('r u v w z h');old_native=[];ratios=[]
for row in old['equations']:
 F=C.from_dict({(0,)+ex+(0,):int(c) for ex,c in S.Poly(S.sympify(row['polynomial']),rr,uu,vv,ww,zz,hh).terms()});expected=co(difference,0,row['coefficient']);ex=next(iter(F.to_dict()));ratio=expected.to_dict()[ex]/F.to_dict()[ex];assert ratio and expected==F*ratio;old_native.append(F);ratios.append(str(ratio))
Z=fmpz_mpoly_ctx.get(['r','u','v','z','delta','w']);zc=Z.gens();Q=fmpq_mpoly_ctx.get(['r','u','v','z','delta','w']);rn,un,vn,zn,dn,wn=Q.gens();H=vn+((zn-un)*(3*rn+zn)-dn)/2

def primitive(F):
 td={ex:Fraction(str(c)) for ex,c in F.to_dict().items()};den=reduce(lcm,(c.denominator for c in td.values()),1);co={ex:int(c*den) for ex,c in td.items()};g=reduce(gcd,(abs(c) for c in co.values()),0);return Z.from_dict({ex:c//g for ex,c in co.items()}),{'denominator_clearing':str(den),'integer_content_removed':str(g)}
transformed=[];normalizations=[]
for F in old_native:
 out=Q.constant(0)
 for ex,c in F.to_dict().items():
  assert ex[0]==ex[7]==0;out+=Q.from_dict({(ex[1],ex[2],ex[3],ex[5],0,ex[4]):c})*H**ex[6]
 f,n=primitive(out);transformed.append(f);normalizations.append(n)
last=transformed[-1]
def cw(F,j):return Z.from_dict({ex[:5]+(0,):int(c) for ex,c in F.to_dict().items() if ex[5]==j})
K,B=cw(last,1),cw(last,0);assert last.degrees()[5]==1;kd=K.to_dict();assert len(kd)==1 and next(iter(kd))==(0,0,0,0,1,0);kconst=int(next(iter(kd.values())))
branches=[]
for branch in ['delta-zero','delta-nonzero']:
 equations=[]
 for index,F in enumerate(transformed if branch=='delta-zero' else transformed[:-1]):
  if branch=='delta-zero':G=Z.from_dict({ex:int(c) for ex,c in F.to_dict().items() if ex[4]==0});td={(ex[0],ex[1],ex[2],ex[3],ex[5]):int(c) for ex,c in G.to_dict().items()};power=0
  else:
   degree=int(F.degrees()[5]);G=sum((cw(F,j)*(-B)**j*K**(degree-j) for j in range(degree+1)),Z.constant(0));power=int(min(ex[4] for ex in G.to_dict()));td={ex[:4]+(ex[4]-power,):int(c) for ex,c in G.to_dict().items()};assert all(ex[5]==0 for ex in G.to_dict())
  content=reduce(gcd,(abs(c) for c in td.values()),0);td={ex:c//content for ex,c in td.items()};variables=['r','u','v','z','w'] if branch=='delta-zero' else ['r','u','v','z','delta'];z5=fmpz_mpoly_ctx.get(variables);g=z5.from_dict(td);equations.append({'old_coefficient':old['equations'][index]['coefficient'],'removed_delta_power':power,'integer_content_removed':str(content),'degrees':list(map(int,g.degrees())),'total_degree':int(max(sum(ex) for ex in td)),'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(td.items())],'polynomial':str(g)})
  print(json.dumps({'branch':branch,'coefficient':old['equations'][index]['coefficient'],'terms':len(td),'total_degree':equations[-1]['total_degree'],'seconds':time.monotonic()-START}),flush=True)
 record={'status':'UNVERIFIED_WORKING_INPUT','variables':variables,'branch':branch,'coordinate':'delta=-2*h-3*r*u+3*r*z-u*z+2*v+z^2','h_recovery':'v+((z-u)*(3*r+z)-delta)/2','w_recovery_nonzero_branch':{'coefficient_constant':kconst,'constant_terms':[[list(map(int,ex)),str(c)] for ex,c in B.to_dict().items()]},'source_sha256':hashlib.sha256(oldpath.read_bytes()).hexdigest(),'fresh_source_ratios':ratios,'h_substitution_normalizations':normalizations,'equations':equations,'boundary':'Necessary coefficient equations for R=(t-r)^3*T3^2*S3. No distinctness/squarefreeness assumed; delta-zero and delta-nonzero are both retained. Any exclusion requires literal rational identities and independent replay. Pole chart r*(3*r^3-2)=0 is handled separately by GAOF7TRP1.'}
 (HERE/(branch+'-fivevar-input.json')).write_text(json.dumps(record,separators=(',',':'))+'\n');branches.append({'branch':branch,'equations':len(equations),'total_terms':sum(len(e['terms']) for e in equations)})
print(json.dumps({'branches':branches,'seconds':time.monotonic()-START}))
