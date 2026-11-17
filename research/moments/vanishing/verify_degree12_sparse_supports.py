#!/usr/bin/env python3
"""Exact Haar and polynomial-identity certificates for all <=4-weight V12 supports.

Default: replay frozen identities. --freeze: construct them by tracked Buchberger.
Standard library only; never infers a tail from a finite zero prefix.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from math import factorial
from collections import Counter
import hashlib,json,sys

import verify_haar_shear_witnesses as S
from verify_degree12_four_weights import groups

ROOT=Path(__file__).resolve().parent
FROZEN=ROOT/'degree12_sparse_certificates.json'

def add(a,b):
 out=a.copy()
 for e,c in b.items():
  out[e]=out.get(e,Q(0))+c
  if not out[e]:del out[e]
 return out

def shift(a,e,c=Q(1)):
 return {tuple(x+y for x,y in zip(k,e)):v*c for k,v in a.items() if v*c}

def mul(a,b):
 out={}
 for e,c in b.items():out=add(out,shift(a,e,c))
 return out

def key(e):return (sum(e),)+tuple(-x for x in reversed(e))
def lead(a):return max(a,key=key)
def divides(a,b):return all(x<=y for x,y in zip(a,b))
def subexp(a,b):return tuple(x-y for x,y in zip(a,b))
def encode(a):return [[list(e),str(c)] for e,c in sorted(a.items())]
def decode(a):return {tuple(e):Q(c) for e,c in a}

def bezout(polys):
 """Return a tracked identity sum multipliers[i]*polys[i]=1 over Q."""
 n=len(polys);nv=len(next(iter(polys[0])));zero=(0,)*nv
 basis=[];reps=[];pairs=[]
 def insert(f,r):
  # Polynomial division, updating its expression in the input generators.
  pending=f.copy();remainder={}
  while pending:
   e=lead(pending);c=pending[e]
   for g,v in zip(basis,reps):
    eg=lead(g)
    if divides(eg,e):
     d=subexp(e,eg);q=-c/g[eg]
     pending=add(pending,shift(g,d,q))
     r=[add(x,shift(y,d,q)) for x,y in zip(r,v)]
     break
   else:
    remainder[e]=c;del pending[e]
  if not remainder:return None
  c=Q(1)/remainder[lead(remainder)]
  remainder=shift(remainder,zero,c);r=[shift(v,zero,c) for v in r]
  if remainder=={zero:Q(1)}:return r
  at=len(basis);pairs.extend((i,at) for i in range(at))
  basis.append(remainder);reps.append(r)
  return None
 for i,f in enumerate(polys):
  r=[{} for _ in range(n)];r[i]={zero:Q(1)}
  found=insert(f,r)
  if found is not None:return found
 while pairs:
  i,j=pairs.pop(0);a,b=basis[i],basis[j];ea,eb=lead(a),lead(b)
  # Relatively prime leading terms have zero S-remainder (product criterion).
  if all(min(x,y)==0 for x,y in zip(ea,eb)):continue
  lcm=tuple(max(x,y) for x,y in zip(ea,eb))
  da,db=subexp(lcm,ea),subexp(lcm,eb)
  ca,cb=1/a[ea],-1/b[eb]
  f=add(shift(a,da,ca),shift(b,db,cb))
  r=[add(shift(x,da,ca),shift(y,db,cb)) for x,y in zip(reps[i],reps[j])]
  found=insert(f,r)
  if found is not None:return found
 raise AssertionError('No unit-ideal certificate: do not claim exclusion')

class Model:
 def __init__(self,name):
  coeffs={2:Q(1),6:Q(-2),10:Q(1)} if name=='2O' else {1:Q(1),6:Q(11),11:Q(-1)}
  field,gens,values=groups(name)
  S.invariant_form(field,gens,coeffs,12)
  assert S.invariant_dimensions(field,values,12)[12]==1
  self.basis=[{(12-i,i,0,0):c for i,c in coeffs.items()}]
  for j in range(1,13):self.basis.append(S.scale(S.lower(self.basis[-1]),Q(1,j)))
  assert self.basis[-1]==S.compact_star(self.basis[0])
  assert not S.lower(self.basis[-1])
  for j,p in enumerate(self.basis):
   assert all(sum(e)==12 and e[0]+e[1]-e[2]-e[3]==12-2*j for e in p)
   # The matrix Weyl element (x,y,z,w)->(z,w,-x,-y).
   reflected={ (e[2],e[3],e[0],e[1]):c*(-1)**(e[2]+e[3]) for e,c in p.items() }
   assert reflected==S.scale(self.basis[12-j],(-1)**j)
  self.powers={(i,0):{(0,0,0,0):Q(1)} for i in range(13)}
  self.integrals={}
 def power(self,i,m):
  if (i,m) not in self.powers:self.powers[i,m]=S.mul(self.power(i,m-1),self.basis[i])
  return self.powers[i,m]
 def integral(self,support,exponents):
  k=tuple((i,m) for i,m in zip(support,exponents) if m)
  if k not in self.integrals:
   p={(0,0,0,0):Q(1)}
   for i,m in k:p=S.mul(p,self.power(i,m))
   self.integrals[k]=S.haar(p)
  return self.integrals[k]
 def moment(self,support,m):
  r=len(support)-2;weights=[6-i for i in support];out={}
  def exponents(n,total):
   if not n:yield ();return
   for b in range(total+1):
    for rest in exponents(n-1,total-b):yield (b,)+rest
  for mid in exponents(r,m):
   num=-weights[-1]*m-sum((w-weights[-1])*b for w,b in zip(weights[1:-1],mid))
   den=weights[0]-weights[-1]
   if num%den:continue
   a=num//den;d=m-a-sum(mid)
   if a<0 or d<0:continue
   powers=(a,)+mid+(d,);factor=Q(factorial(m))
   for b in powers:factor/=factorial(b)
   c=factor*self.integral(support,powers)
   if c:out[mid]=c
  if not out:return None
  valuation=tuple(min(e[i] for e in out) for i in range(r))
  out={subexp(e,valuation):c for e,c in out.items()}
  scalar=out[lead(out)];out={e:c/scalar for e,c in out.items()}
  return valuation,scalar,out

def supports(size):
 allsets={q for q in combinations(range(13),size) if q[0]<6<q[-1]}
 reflect=lambda q:tuple(12-i for i in reversed(q))
 reps=sorted(q for q in allsets if q<=reflect(q))
 assert set(reps)|{reflect(q) for q in reps}==allsets
 return reps,len(allsets)

def main():
 freeze='--freeze' in sys.argv
 frozen={} if freeze else json.loads(FROZEN.read_text())['cases']
 cases={};summary={}
 for name in ['2O','2I']:
  model=Model(name);rows=[];counts={};orders={}
  for size in [2,3,4]:
   reps,total=supports(size);counts[str(size)]={'crossing_supports':total,'weyl_representatives':len(reps)}
   ordercount=Counter()
   for support in reps:
    r=size-2;nv=r+1;zero=(0,)*nv;unit={zero:Q(1)}
    saturation={tuple([1]*nv):Q(1),zero:Q(-1)}
    if not r:saturation={(1,):Q(1),(0,):Q(-1)}
    row={'support':list(support),'moments':[]}
    source=None if freeze else next(x for x in frozen[name] if x['support']==list(support))
    polys=[];bound=11 if source is None else source['moments'][-1]['m']
    assert 2<=bound<=11
    for m in range(2,bound+1):
     result=model.moment(support,m)
     if result is None:continue
     val,scalar,poly=result
     item={'m':m,'valuation':list(val),'scalar':str(scalar),'polynomial':encode(poly)}
     row['moments'].append(item)
     polys.append({e+(0,):c for e,c in poly.items()})
     if freeze:
      try:cert=bezout(polys+[saturation])
      except AssertionError:continue
      row['multipliers']=[encode(x) for x in cert];break
    else:
     if freeze:raise AssertionError(('unresolved',name,support))
    if not freeze:
     assert row['moments']==source['moments']
     row['multipliers']=source['multipliers']
    cert=[decode(x) for x in row['multipliers']]
    assert len(cert)==len(polys)+1
    identity={}
    for a,b in zip(cert,polys+[saturation]):identity=add(identity,mul(a,b))
    assert identity==unit
    ordercount[row['moments'][-1]['m']]+=1;rows.append(row)
   orders[str(size)]=dict(sorted(ordercount.items()))
  hzero=S.haar(model.power(6,2))
  assert hzero==(Q(32,13) if name=='2O' else Q(275,13))
  controls={}
  for support,exps,label in [((0,11),(5,6),'prefix11'),((5,12),(6,1),'marked_sharp')]:
   value=model.integral(support,exps);assert value
   if label=='prefix11':
    for m in range(1,11):assert model.moment(support,m) is None
    value*=Q(factorial(11),factorial(5)*factorial(6))
   controls[label]=str(value)
  cases[name]=rows
  summary[name]={'support_counts':counts,'certificate_orders':orders,'certificates':len(rows),
                 'zero_weight_square':str(hzero),'controls':controls,
                 'cached_balanced_integrals':len(model.integrals)}
 if freeze:FROZEN.write_text(json.dumps({'cases':cases},indent=2)+'\n')
 else:assert cases==frozen
 receipt={'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'certificate_sha256':hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
          'cases':summary,'all_exact_identities_pass':True,
          'full_V12_decided':False,'full_exceptional_MZ_decided':False,
          'independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'degree12_sparse_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
