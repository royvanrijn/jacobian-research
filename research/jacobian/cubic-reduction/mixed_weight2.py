#!/usr/bin/env python3
"""Complete primary weighted models(2,-1,-r),4<=r<=13,degree<=7."""
from pathlib import Path
from itertools import product
import sympy as s
import json,hashlib,time
HERE=Path(__file__).resolve().parent;OUT=HERE/'mixed-weight2';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*vs).terms() if c]
def main():
 start=time.monotonic();a,b,c,u,v=s.symbols('a b c u v');cp={'status':'running','source_sha256':sha(Path(__file__)),'degree_cap':7,'r_values':list(range(4,14)),'models':[],'preflight':'KCUB32 r1,KCUB33 r2,KCUB17 r3 read as separate closed scopes; complete registry/source/archive scans find no r>=4 closure. KCUB34 has different positive weight1. No previous UNKNOWN/timeouts replayed.'};path=OUT/'checkpoint.json'
 for r in range(4,14):
  basis=[[e for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and 2*e[0]-e[1]-r*e[2]==t] for t in [2,-1,-r]];names=s.symbols('c1:'+str(sum(map(len,basis))+1));pos=0;F=[a,b,c]
  for j,row in enumerate(basis):
   for e in row:F[j]+=names[pos]*a**e[0]*b**e[1]*c**e[2];pos+=1
  quot=[s.expand(z.subs({a:1,b:u,c:v})) for z in F];f,g,h=quot;bracket=lambda p,q:s.diff(p,u)*s.diff(q,v)-s.diff(p,v)*s.diff(q,u);D=s.expand(f*bracket(g,h)+s.Rational(1,2)*g*bracket(f,h)+s.Rational(r,2)*h*bracket(g,f));J=s.expand(s.det(s.Matrix(F).jacobian([a,b,c])));assert s.expand(J.subs(a,1)-D.subs({u:b,v:c}))==0
  # Complete source weights uniquely recover every a exponent; parity makes the half-cover descend.
  model={'status':'complete-primary-model','source_sha256':cp['source_sha256'],'r':r,'weights':[2,-1,-r],'basis':[[list(e) for e in row] for row in basis],'parameter_names':list(map(str,names)),'source_outputs':[enc(z,[a,b,c]+list(names)) for z in F],'double_cover_outputs':[enc(z,[u,v]+list(names)) for z in quot],'quotient_jacobian':enc(D,[u,v]+list(names)),'native_generators':[enc(z,list(names)) for z in s.Poly(D-1,u,v).coeffs()],'slots':list(map(len,basis)),'classification_status':'unresolved-primary-model-only'}
  p=OUT/('r'+str(r)+'.json');p.write_text(json.dumps(model,indent=2)+'\n');cp['models'].append({'r':r,'input':p.name,'input_sha256':sha(p),'slots':list(map(len,basis)),'native_generators':len(model['native_generators'])});path.write_text(json.dumps(cp,indent=2)+'\n');print('MODEL',r,model['slots'],len(model['native_generators']),flush=True)
 cp['status']='complete';cp['elapsed_seconds']=time.monotonic()-start;path.write_text(json.dumps(cp,indent=2)+'\n')
if __name__=='__main__':main()
