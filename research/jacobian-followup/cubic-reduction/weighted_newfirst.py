#!/usr/bin/env python3
"""Prepare and search new first-output terms in the negative-weight degree6/7 families."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def decode(terms,V):return s.Add(*(s.Rational(c)*s.prod(x**a for x,a in zip(V,e)) for e,c in terms))
def prepare(degree,out):
 u,w=s.symbols('u w');supports=[]
 for row,target in enumerate([1,-1,-2]):
  atoms=[]
  for b in range(degree+1):
   for c in range(degree+1):
    a=target+b+2*c
    if a<0 or not 2<=a+b+c<=degree:continue
    if row==0 and (b,c) in [(1,0),(0,1)]:continue
    atoms.append((b,c,(a,b,c)))
  supports.append(sorted(atoms,reverse=True))
 groups=[s.symbols(prefix+'1:'+str(len(atoms)+1)) for prefix,atoms in zip(['f','g','h'],supports)];C=tuple(x for group in groups for x in group)
 f=1-3*u-w+sum(a*u**b*w**c for a,(b,c,e) in zip(groups[0],supports[0]));g=u+sum(a*u**b*w**c for a,(b,c,e) in zip(groups[1],supports[1]));h=w+sum(a*u**b*w**c for a,(b,c,e) in zip(groups[2],supports[2]))
 J=s.expand(g*(s.diff(f,u)*s.diff(h,w)-s.diff(f,w)*s.diff(h,u))+2*h*(s.diff(g,u)*s.diff(f,w)-s.diff(g,w)*s.diff(f,u))+f*(s.diff(g,u)*s.diff(h,w)-s.diff(g,w)*s.diff(h,u))-1)
 jac=[s.expand(c) for e,c in s.Poly(J,u,w).terms() if c];eqs=jac+[s.expand(f.subs({u:-1,w:4})),s.expand(g.subs({u:-1,w:4})+1)];linear=[p for p in eqs if p and s.Poly(p,*C).total_degree()==1]
 M,b=s.linear_eq_to_matrix(linear,C);R,pivots=M.row_join(b).rref();assert all(i<len(C) for i in pivots);free=[i for i in range(len(C)) if i not in pivots];V=tuple(C[i] for i in free);sub={C[p]:s.expand(R[k,len(C)]-sum(R[k,j]*C[j] for j in free)) for k,p in enumerate(pivots)}
 reduced=[];seen=set()
 for p in eqs:
  p=s.expand(p.subs(sub,simultaneous=True))
  if not p:continue
  p=s.Poly(p,*V).monic().as_expr()
  if str(p) not in seen:seen.add(str(p));reduced.append(p)
 atoms=[dict(row=0,exponent=[2,1,0],coefficient='-3'),dict(row=0,exponent=[3,0,1],coefficient='-1')]
 atoms += [dict(row=i,exponent=list(e),coefficient=str(a),quotient_exponent=[b,c]) for i,(coefs,mon) in enumerate(zip(groups,supports)) for a,(b,c,e) in zip(coefs,mon)]
 payload=dict(source_sha256=sha(Path(__file__)),degree=degree,variables=list(map(str,V)),atoms=atoms,substitution={str(c):str(sub.get(c,c)) for c in C},equations=[encode(p,V) for p in reduced],linear_rank=len(pivots),scope='All weight-compatible atoms through the declared degree, with coefficients of a^2b and a^3c fixed at-3,-1 and collision between(0,1,c0) and(-1,1,4). The third collision coordinate is recovered linearly.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('PREPARED',degree,len(atoms),len(V),len(reduced),flush=True)
def main():
 degree=int(sys.argv[1]);mode=sys.argv[2];assert degree in [6,7] and mode in ['prepare','nra','groebner'];out=HERE/('weighted-newfirst%d'%degree);out.mkdir(exist_ok=True)
 if mode=='prepare':prepare(degree,out);return
 source=out/'input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));eqs=[decode(p,V) for p in d['equations']];started=time.monotonic()
 receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),degree=degree,mode=mode,scope='New-first-output weighted coefficient search. Exact candidate determinant/collision verification required; solver timeouts remain UNKNOWN.')
 if mode=='groebner':
  script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n')
  try:
   r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;receipt['status']='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
  except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');receipt['status']='timeout'
  (out/'model.log').write_text(text)
 else:
  import z3
  Z={str(v):z3.Real(str(v)) for v in V}
  def convert(f):
   terms=[]
   for e,c in s.Poly(f,*V).terms():
    t=z3.RealVal(str(c))
    for name,k in zip(d['variables'],e):
     if k:t*=Z[name]**k
    terms.append(t)
   return z3.Sum(terms)
  solver=z3.SolverFor('QF_NRA');solver.set(timeout=30000)
  for f in eqs:solver.add(convert(f)==0)
  if degree==7:
   byname={str(v):v for v in V};active=[]
   for i,atom in enumerate(d['atoms']):
    expr=s.sympify(d['substitution'].get(atom['coefficient'],atom['coefficient']),locals=byname);a=z3.Bool('active%d'%i);active.append(a);solver.add(z3.Implies(z3.Not(a),convert(expr)==0))
   solver.add(z3.PbLe([(a,1) for a in active],12));receipt['nonlinear_bound']=12
  answer=solver.check();receipt['status']=str(answer);receipt['solver_version']=z3.get_version_string()
  if answer==z3.sat:
   model=solver.model();values={name:model.eval(v,model_completion=True) for name,v in Z.items()};receipt['values']={name:str(v) for name,v in values.items()};receipt['expressions']={name:v.sexpr() for name,v in values.items()};receipt['all_rational']=all(z3.is_rational_value(v) for v in values.values())
   if receipt['all_rational']:
    point={v:s.Rational(receipt['values'][str(v)]) for v in V};byname={str(v):v for v in V};receipt['coefficients']={name:str(s.sympify(f,locals=byname).subs(point)) for name,f in d['substitution'].items()}
  elif answer==z3.unknown:receipt['reason']=solver.reason_unknown()
 receipt['seconds']=round(time.monotonic()-started,3);(out/(mode+'-checkpoint.json')).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
