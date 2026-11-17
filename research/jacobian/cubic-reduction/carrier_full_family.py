#!/usr/bin/env python3
"""Explicit eight-parameter cubic carrier family; completeness not assumed."""
from pathlib import Path
import hashlib,json,time
import sympy as S
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ip=HERE/'carrier-reparam-full/input.json';d=json.loads(ip.read_text());Z=S.symbols(' '.join(d['variables']));p,q,s,t,u,v,w,z=S.symbols('p q s t u v w z');values=[S.Integer(0)]*135
 active={1:1/p,2:1/p,3:1,4:p*q,5:1/(p*q),6:1-3*u/(p**2*q),7:p,8:p,9:p,10:q,11:p*q+3*t,12:q,13:s/(2*p),17:s,18:t/q,19:-3*u/(p*q),22:t,26:(u-3*v/(p*q))/p,27:w/p,28:v/p,30:z/q,34:u/q,35:v/q,38:p*q*s+3*z,39:-p*q*z,42:-p*q*t,45:p*q*w,79:u-3*v/(p*q),80:w,83:v,96:z,117:u,118:v}
 for j,value in active.items():values[j-1]=S.sympify(value)
 sub=dict(zip(Z,values));start=time.monotonic()
 for terms in d['high_equations']:
  f=sum(S.Rational(c)*S.prod(sub[zz]**k for zz,k in zip(Z,e) if k) for e,c in terms);assert S.cancel(f)==0,f
 print('HIGH SUBSTITUTION PASS',len(d['high_equations']),flush=True);rows=[]
 for a in d['low_coefficients']:
  f=S.factor(S.sympify(a['coefficient'],locals=dict(zip(map(str,Z),Z))).subs(sub,simultaneous=True))
  if f:rows.append(dict(row=a['row'],exponent=a['exponent'],coefficient=str(f)))
 payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(ip),variables=list(map(str,[p,q,s,t,u,v,w,z])),nonzero=['p','q'],parameterization={str(zz):str(value) for zz,value in zip(Z,values)},low_coefficients=rows,all_high_coefficients_zero=True,seconds=round(time.monotonic()-start,3),scope='Explicit eight-parameter sufficient cubic family in the complete135-atom triangular support. No completeness or global minimum assumed from the primary basis. Any improved witness needs independent exact composition/determinant/collision verification.');(ip.parent/'family.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('GENERIC LOW',len(rows),flush=True)
 for f in sorted({str(S.factor(S.sympify(a['coefficient']).subs({p:1,q:1}))) for a in rows}):print(f)
if __name__=='__main__':main()
