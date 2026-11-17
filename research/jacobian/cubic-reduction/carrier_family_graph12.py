#!/usr/bin/env python3
"""Changed raw graph slices in the full cubic carrier family, p=q=1."""
from pathlib import Path
import hashlib,json,time
import sympy as S
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 source=HERE/'carrier-reparam-wide/family.json';d=json.loads(source.read_text());X=S.symbols('x1:13');s,t=S.symbols('s t');p,q=S.symbols('p q');G=list(X)
 for a in d['low_coefficients']:G[a['row']]+=S.sympify(a['coefficient']).subs({p:1,q:1})*S.prod(x**k for x,k in zip(X,a['exponent']))
 out=HERE/'carrier-family-graphs12';out.mkdir(exist_ok=True);cases=[]
 for j in range(3,12):
  start=time.monotonic();graph=S.expand(X[j]-G[j]);assert not graph.has(X[j]);retained={i:S.expand(f.subs(X[j],graph)) for i,f in enumerate(G) if i!=j};high=[c for f in retained.values() for e,c in S.Poly(f,*X).terms() if sum(e)>3 and c];eqs=list(dict.fromkeys(high));B=S.groebner(eqs,s,t);result=dict(pivot=j+1,graph=str(graph),parameters=['s','t'],high_equations=list(map(str,eqs)),basis=list(map(str,B.polys)),raw_cubic_status='unit-ideal-no-raw-cubic' if list(B)==[1] else 'proper-ideal-needs-extraction',source_sha256=sha(Path(__file__)),input_sha256=sha(source),seconds=round(time.monotonic()-start,3));(out/('pivot-%02d.json'%(j+1))).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');cases.append(dict(pivot=j+1,status=result['raw_cubic_status'],high_equations=len(eqs)));print(cases[-1],flush=True)
 (out/'checkpoint.json').write_text(json.dumps(dict(status='complete',source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,scope='Changed raw one-carrier graph restrictions of the proved cubic family atp=q=1; primary coefficient basis only, no independent replay or target repairs. No global dimension bound.'),indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
