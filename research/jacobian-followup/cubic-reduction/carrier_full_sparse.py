#!/usr/bin/env python3
"""Exact rational support query on the six-parameter linear coefficient family."""
from pathlib import Path
from collections import Counter
import hashlib,json,time
import sympy as S
import z3
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 source=HERE/'carrier-reparam-full/family.json';d=json.loads(source.read_text());p,q=S.symbols('p q');P=S.symbols('s t u v w z');lookup=dict(zip(map(str,P),P));groups=Counter();units=0
 for a in d['low_coefficients']:
  f=S.Poly(S.sympify(a['coefficient'],locals=lookup).subs({p:1,q:1}),*P,domain=S.QQ);assert f.total_degree()<=1
  if f.total_degree()==0:units+=1
  else:groups[str(f.monic().as_expr())]+=1
 A={str(v):z3.Real(str(v)) for v in P};solver=z3.Solver();solver.set(timeout=20000);constraints=[];clauses=[]
 for text,weight in groups.items():
  f=S.Poly(S.sympify(text,locals=lookup),*P);expr=z3.RealVal(0)
  for e,c in f.terms():
   term=z3.RealVal(str(c))
   for k,n in zip(P,e):
    if n:term=term*A[str(k)]
   expr+=term
  constraints.append((text,weight,expr));clauses.append((expr!=0,weight))
 out=HERE/'carrier-reparam-full';tests=[]
 for bound in [35,34,33,32,31,30]:
  solver.push();solver.add(z3.PbLe(clauses,bound-units));started=time.monotonic();answer=solver.check();row=dict(bound=bound,status=str(answer),seconds=round(time.monotonic()-started,3));
  if answer==z3.sat:
   m=solver.model();values={k:str(m.eval(v,model_completion=True)) for k,v in A.items()};total=units+sum(weight for text,weight,expr in constraints if z3.is_true(m.eval(expr!=0)));row['parameters']=values;row['nonlinear_total']=total;(out/('sparse-candidate-'+str(bound)+'.json')).write_text(json.dumps(dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters={'p':'1','q':'1'}|values,nonlinear_total=total,scope='Exact primary affine coefficient model; requires independent full polynomial composition/collision verification before a new witness.'),indent=2,sort_keys=True)+'\n')
  elif answer==z3.unknown:row['reason']=solver.reason_unknown()
  solver.pop();tests.append(row);print(row,flush=True)
  if answer!=z3.sat:break
 receipt=dict(status='complete' if tests[-1]['status']!='unknown' else 'unknown',source_sha256=sha(Path(__file__)),input_sha256=sha(source),constant_nonlinear_occurrences=units,hyperplanes=[dict(polynomial=text,multiplicity=weight) for text,weight in groups.items()],tests=tests,scope='Support queries only on p=q1 explicit six-parameter sufficient family. Solver UNSAT is not an independent rational certificate or a global minimum. SAT requires full witness replay.');(out/'sparse-checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('UNITS',units,'HYPERPLANES',len(groups),flush=True)
if __name__=='__main__':main()
