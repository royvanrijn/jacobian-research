#!/usr/bin/env python3
"""Exact real-algebraic sparse search in the enlarged19-atom weighted tail model."""
from pathlib import Path
import hashlib,json,time
import sympy as s
import z3
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'weighted-tail7/nra';out.mkdir(exist_ok=True);source=HERE/'weighted-tail7/generic/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));Z={str(v):z3.Real(str(v)) for v in V}
 def convert(poly):
  p=s.Poly(poly,*V);terms=[]
  for e,c in p.terms():
   term=z3.RealVal(str(c))
   for name,k in zip(d['variables'],e):
    if k:term*=Z[name]**k
   terms.append(term)
  return z3.Sum(terms)
 eqs=[sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in f) for f in d['equations']]
 original=json.loads((HERE/'weighted-tail7/input.json').read_text());byname={str(v):v for v in V};coeff=[s.sympify(d['substitution'][a['coefficient']],locals=byname) for a in original['tail_atoms']]
 solver=z3.SolverFor('QF_NRA');solver.set(timeout=40000)
 for f in eqs:solver.add(convert(f)==0)
 active=[z3.Bool('active_%d'%i) for i in range(17)]
 for a,f in zip(active,coeff):solver.add(z3.Implies(z3.Not(a),convert(f)==0))
 solver.add(z3.PbLe([(a,1) for a in active],10))
 receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),nonlinear_bound=12,timeout_ms=40000,scope='Real-algebraic sparse search in the enlarged19-atom fixed-first-output model. An exact candidate still needs full coefficient, determinant and collision replay; UNSAT is not a complex-field exclusion.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');started=time.monotonic();answer=solver.check();receipt.update(status=str(answer),seconds=round(time.monotonic()-started,3),solver_version=z3.get_version_string())
 if answer==z3.sat:
  model=solver.model();values={name:model.eval(value,model_completion=True) for name,value in Z.items()};receipt['values']={name:str(value) for name,value in values.items()};receipt['value_expressions']={name:value.sexpr() for name,value in values.items()};receipt['all_rational']=all(z3.is_rational_value(v) for v in values.values())
  if receipt['all_rational']:
   sub={v:s.Rational(receipt['values'][str(v)]) for v in V};receipt['original_coefficients']={a['coefficient']:str(f.subs(sub)) for a,f in zip(original['tail_atoms'],coeff)};receipt['original_auxiliary']={name:str(s.sympify(expr,locals=byname).subs(sub)) for name,expr in d['substitution'].items() if name in ['beta','sat']}
 elif answer==z3.unknown:receipt['reason']=solver.reason_unknown()
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
