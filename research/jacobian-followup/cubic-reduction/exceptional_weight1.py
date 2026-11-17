#!/usr/bin/env python3
"""Complete primary degree-seven models for weights(1,-1,-r),r=4..7."""
from pathlib import Path
import sympy as sp
import json,hashlib,time
from itertools import product
HERE=Path(__file__).resolve().parent;OUT=HERE/'exceptional-weight1';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{"exponent":list(e),"coefficient":str(c)} for e,c in sp.Poly(sp.expand(p),*vs).terms() if c]
def main():
 started=time.monotonic();a,b,c,u,w=sp.symbols('a b c u w');checkpoint={'status':'running','source_sha256':sha(Path(__file__)),'degree_cap':7,'weights':[4,5,6,7],'models':[],'prior_scope_preflight':'KCUB15 full112,KCUB18 full113,KCUB33 m>=2 exception; complete current registry and source/archive text scan found no114..117 closure. Original UNKNOWN/timeouts are not reused.'};path=OUT/'checkpoint.json'
 for r in [4,5,6,7]:
  basis=[[e for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and e[0]-e[1]-r*e[2]==target] for target in [1,-1,-r]];names=[sp.Symbol('c'+str(j+1)) for j in range(sum(map(len,basis)))];offset=0;F=[]
  for j,es in enumerate(basis):
   F.append([a,b,c][j]+sum(names[offset+k]*a**e[0]*b**e[1]*c**e[2] for k,e in enumerate(es)));offset+=len(es)
  quotient=[sp.expand(F[0].subs({a:1,b:u,c:w})),sp.expand(F[1].subs({a:1,b:u,c:w})),sp.expand(F[2].subs({a:1,b:u,c:w}))];f,g,h=quotient
  # a=1 loses no coefficient: each negative monomial fixes its a exponent uniquely.
  bracket=lambda p,q:sp.diff(p,u)*sp.diff(q,w)-sp.diff(p,w)*sp.diff(q,u)
  D=sp.expand(g*bracket(f,h)+r*h*bracket(g,f)+f*bracket(g,h));J=sp.expand(sp.det(sp.Matrix(F).jacobian([a,b,c])));assert sp.expand(J-D.subs({u:a*b,w:a**r*c}))==0
  model={'status':'complete-primary','source_sha256':checkpoint['source_sha256'],'r':r,'weights':[1,-1,-r],'basis':[[list(e) for e in es] for es in basis],'parameter_names':list(map(str,names)),'source_outputs':[enc(p,[a,b,c]+names) for p in F],'quotient_outputs':[enc(p,[u,w]+names) for p in quotient],'quotient_jacobian':enc(D,[u,w]+names),'native_generators':[enc(v,names) for v in sp.Poly(D-1,u,w).coeffs()],'nonlinear_slots':sum(map(len,basis)),'per_output_slots':list(map(len,basis)),'classification_status':'unresolved-primary-model-only'}
  p=OUT/('r'+str(r)+'.json');p.write_text(json.dumps(model,indent=2)+'\n');checkpoint['models'].append({'r':r,'input':p.name,'input_sha256':sha(p),'nonlinear_slots':model['nonlinear_slots'],'generators':len(model['native_generators'])});path.write_text(json.dumps(checkpoint,indent=2)+'\n');print('MODEL',r,list(map(len,basis)),len(model['native_generators']),flush=True)
 checkpoint['status']='complete';checkpoint['elapsed_seconds']=time.monotonic()-started;path.write_text(json.dumps(checkpoint,indent=2)+'\n')
if __name__=='__main__':main()
