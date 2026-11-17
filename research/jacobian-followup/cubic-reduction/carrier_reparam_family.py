#!/usr/bin/env python3
"""Save the explicit complete family suggested by the47-parameter basis."""
from pathlib import Path
import hashlib,json
import sympy as S
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ip=HERE/'carrier-reparam-wide/input.json';d=json.loads(ip.read_text());Z=S.symbols(' '.join(d['variables']));p,q,s,t=S.symbols('p q s t');values=[S.Integer(0)]*47
 for j,v in {1:1/p,2:1/p,3:1,4:p*q,5:1/(p*q),6:1,7:p,8:p,9:p,10:q,11:p*q+3*t,12:q,13:s/(2*p),17:s,18:t/q,22:t,38:p*q*s,42:-p*q*t}.items():values[j-1]=S.sympify(v)
 sub=dict(zip(Z,values));rows=[]
 for a in d['low_coefficients']:
  f=S.factor(S.sympify(a['coefficient'],locals=dict(zip(map(str,Z),Z))).subs(sub,simultaneous=True))
  if f:rows.append(dict(row=a['row'],exponent=a['exponent'],coefficient=str(f)))
 payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(ip),variables=['p','q','s','t'],nonzero=['p','q'],parameterization={str(v):str(x) for v,x in zip(Z,values)},low_coefficients=rows,scope='Explicit proposed complete cubic family. Completeness requires independent43-membership replay and written elimination; sufficiency and exact support formula require independent substitution of every high/low coefficient.')
 (ip.parent/'family.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('FAMILY',len(rows),'generic nonlinear occurrences')
if __name__=='__main__':main()
