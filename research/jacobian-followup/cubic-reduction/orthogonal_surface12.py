#!/usr/bin/env python3
"""Exact two-parameter completion of the previously tested orthogonal line."""
from pathlib import Path
import json,hashlib,time
import sympy as s
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def representation(a,b,c,d,sign=1):
 L=s.eye(12)
 L[10,:]=s.zeros(1,12);L[10,10]=a*a;L[10,1]=2*a*b;L[10,6]=b*b/2
 L[1,:]=s.zeros(1,12);L[1,10]=a*c;L[1,1]=a*d+b*c;L[1,6]=b*d/2
 L[6,:]=s.zeros(1,12);L[6,10]=2*c*c;L[6,1]=4*c*d;L[6,6]=d*d
 for i in (1,6,10):L[i,:]=sign*L[i,:]
 L[9,:]=s.eye(12)[9,:]+s.eye(12)[1,:]-L[1,:]
 return L

def run(name,params,abcd,sign=1):
 out=HERE/'orthogonal-surface12'/name;out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());X=s.symbols('x1:13')
 F=s.Matrix([sum(s.Rational(v['coefficient'])*s.prod(x**e for x,e in zip(X,v['exponent'])) for v in row) for row in raw['coefficients']])
 a,b,c,d=abcd;L=representation(a,b,c,d,sign);inverse=representation(d,-b,-c,a,sign)
 assert (L*inverse).applyfunc(s.expand)==s.eye(12)
 new=L*s.Matrix(X)
 assert s.expand(2*new[1]**2-new[10]*new[6]-(2*X[1]**2-X[10]*X[6]))==0
 assert s.expand(new[1]+new[9]-X[1]-X[9])==0
 G=(inverse*F.subs(dict(zip(X,new)),simultaneous=True)).applyfunc(s.expand)
 factors=[];index={};groups={};coefficients=[]
 for ri,f in enumerate(G):
  row=[]
  for ex,co in s.Poly(f,*X).terms():
   row.append({'exponent':list(ex),'coefficient':str(co)})
   if sum(ex)<2:continue
   _,fac=s.factor_list(co,*params);labels=[]
   for f0,power in fac:
    f0=s.Poly(f0,*params).monic().as_expr();key=str(f0)
    if key not in index:index[key]=len(factors);factors.append(key)
    labels.append(index[key])
   labels=tuple(sorted(set(labels)));groups[labels]=groups.get(labels,0)+1
  coefficients.append(row)
 payload={'status':'symbolic-model-complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'parameters':list(map(str,params)),'sign':sign,'matrix':[[str(x) for x in row] for row in L.tolist()],'inverse':[[str(x) for x in row] for row in inverse.tolist()],'factors':factors,'factor_groups':[{'factors':list(k),'multiplicity':v} for k,v in sorted(groups.items())],'coefficients':coefficients,'generic_nonlinear_count':sum(groups.values()),'constant_occurrences':groups.get((),0),'seconds':round(time.monotonic()-start,3),'scope':'Complete SL2 congruence subgroup preserving2*x2^2-x11*x7 andx2+x10, modulo right diagonal source scaling; chart decomposition requires separate proof and exact stratum coverage. Not full linear orbit.'}
 (out/'model.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
 print(name,'GENERIC',sum(groups.values()),'CONSTANT',groups.get((),0),'FACTORS',list(enumerate(factors)),flush=True)
 print('GROUPS',payload['factor_groups'],flush=True)
if __name__=='__main__':
 t,u=s.symbols('t u');run('open',(t,u),(1+t*u,t,u,s.Integer(1)))
 a=s.symbols('a');run('boundary',(a,),(a,s.Integer(1),s.Integer(-1),s.Integer(0)))
 run('open-minus',(t,u),(1+t*u,t,u,s.Integer(1)),-1)
 run('boundary-minus',(a,),(a,s.Integer(1),s.Integer(-1),s.Integer(0)),-1)
