#!/usr/bin/env python3
"""Bounded degree7, raw<=12 weighted Keller search using an automatic collision.

For one positive and two ordered negative weights, every nonlinear first
component factors as a*f with f(0,b,c)=1. Its second zero-fibre component
maps into the same plane as a=0, whose tail map is triangular and onto.
Thus determinant1 plus a nonlinear first component already forces a
collision over an algebraic closure; collision variables are unnecessary.
Actual maps/points must still be extracted and independently replayed.
"""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
ALLOWED=[(1,-1,-3),(2,-1,-3),(1,-2,-3)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 weights=tuple(map(int,sys.argv[1].split(',')));mode=sys.argv[2];assert weights in ALLOWED and mode in ['prepare','nra','groebner'];tag='_'.join(map(str,weights));out=HERE/'weighted-sparse7'/tag;out.mkdir(parents=True,exist_ok=True)
 if mode=='prepare':
  X=sp.symbols('a b c');supports=[]
  for row,target in enumerate(weights):
   for a in range(8):
    for b in range(8):
     for c in range(8):
      if 2<=a+b+c<=7 and sum(w*e for w,e in zip(weights,[a,b,c]))==target:supports.append((row,(a,b,c)))
  C=sp.symbols('c1:'+str(len(supports)+1));F=list(X)
  for c,(row,e) in zip(C,supports):F[row]+=c*sp.prod(x**k for x,k in zip(X,e))
  det=sp.expand(sp.det(sp.Matrix(F).jacobian(X))-1);eqs=[c for e,c in sp.Poly(det,*X).terms() if c];linear=[f for f in eqs if sp.Poly(f,*C).total_degree()==1];M,b=sp.linear_eq_to_matrix(linear,C);R,pivots=M.row_join(b).rref();assert all(i<len(C) for i in pivots);V=tuple(C[i] for i in range(len(C)) if i not in pivots);sub={C[p]:sp.expand(R[k,len(C)]-sum(R[k,j]*C[j] for j in range(len(C)) if j not in pivots)) for k,p in enumerate(pivots)};reduced=[];seen=set()
  for f in eqs:
   f=sp.expand(f.subs(sub,simultaneous=True))
   if f:
    f=sp.Poly(f,*V).monic().as_expr()
    if str(f) not in seen:seen.add(str(f));reduced.append(f)
  payload=dict(source_sha256=sha(Path(__file__)),weights=list(weights),max_degree=7,raw_nonlinear_bound=12,variables=list(map(str,V)),linear_rank=len(pivots),substitution={str(c):str(sub.get(c,c)) for c in C},atoms=[dict(row=row,exponent=list(e),coefficient=str(c)) for c,(row,e) in zip(C,supports)],equations=[enc(f,V) for f in reduced],scope='Complete degree<=7 weighted determinant coefficient model with identity linear part. Search raw nonlinear occurrences<=12 and a nonlinear first component, which implies a collision over an algebraic closure by the written triangular zero-fibre argument. No point or counterexample is claimed until independent exact determinant/support/collision replay.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(tag,'ATOMS',len(C),'FREE',len(V),'EQS',len(reduced),'LINEAR_RANK',len(pivots),flush=True);return
 source=out/'input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],V));eqs=[dec(p,V) for p in d['equations']];started=time.monotonic();receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),weights=list(weights),max_degree=7,raw_nonlinear_bound=12,mode=mode,scope='A proper coefficient ideal or solver SAT is not a verified map; extraction and exact replay required. Timeouts remain UNKNOWN. No global minimum claim.')
 if mode=='groebner':
  script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n')
  try:
   r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;receipt['status']='primary-coefficient-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
  except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');receipt['status']='timeout'
  (out/'model.log').write_text(text);receipt['scope']='Unrestricted determinant coefficient basis; no nonlinear-first or raw<=12 restrictions imposed in this exploratory basis. Any point requires actual map replay. Proper ideal dimension does not establish a counterexample.'
 else:
  import z3
  Z={str(v):z3.Real(str(v)) for v in V}
  def cv(f):
   result=[]
   for e,c in sp.Poly(f,*V).terms():
    term=z3.RealVal(str(c))
    for v,k in zip(V,e):
     if k:term*=Z[str(v)]**k
    result.append(term)
   return z3.Sum(result)
  solver=z3.SolverFor('QF_NRA');solver.set(timeout=30000)
  for f in eqs:solver.add(cv(f)==0)
  active=[];first=[]
  for i,a in enumerate(d['atoms']):
   f=sp.sympify(d['substitution'][a['coefficient']],locals=by);atom=z3.Bool('active%d'%i);active.append(atom);solver.add(z3.Implies(z3.Not(atom),cv(f)==0))
   if a['row']==0:first.append(cv(f)!=0)
  solver.add(z3.PbLe([(v,1) for v in active],12),z3.Or(first));answer=solver.check();receipt['status']=str(answer);receipt['solver_version']=z3.get_version_string()
  if answer==z3.sat:
   model=solver.model();values={str(v):model.eval(Z[str(v)],model_completion=True) for v in V};receipt['values']={name:str(v) for name,v in values.items()};receipt['sexprs']={name:v.sexpr() for name,v in values.items()};receipt['all_rational']=all(z3.is_rational_value(v) for v in values.values())
   if receipt['all_rational']:
    point={v:sp.Rational(receipt['values'][str(v)]) for v in V};receipt['coefficients']={name:str(sp.sympify(text,locals=by).subs(point)) for name,text in d['substitution'].items()}
  elif answer==z3.unknown:receipt['reason']=solver.reason_unknown()
 receipt['seconds']=round(time.monotonic()-started,3);(out/(mode+'-checkpoint.json')).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
