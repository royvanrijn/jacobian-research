#!/usr/bin/env python3
"""Joint O(2*b^2-U*S) conjugacy preserving e=b+T, with exact line roots."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 source=HERE/'cubic12-integer-witness.json';d=json.loads(source.read_text());X=s.symbols('x1:13');t=s.Symbol('t');F=s.Matrix([sum(s.Rational(a['coefficient'])*s.prod(x**k for x,k in zip(X,a['exponent'])) for a in row) for row in d['coefficients']]);L=s.eye(12);L[1,6]=t;L[9,6]=-t;L[10,1]=4*t;L[10,6]=2*t**2;inv=L.subs(t,-t);assert L*inv==s.eye(12);new=L*s.Matrix(X);assert s.expand(2*new[1]**2-new[10]*new[6])==2*X[1]**2-X[10]*X[6] and s.expand(new[1]+new[9])==X[1]+X[9];G=(inv*F.subs(dict(zip(X,new)),simultaneous=True)).applyfunc(s.expand);rows=[s.Poly(f,*X) for f in G];coeff=[c for f in rows for e,c in f.terms() if sum(e)>=2];roots={s.Rational(0)}
 for c in coeff:
  p=s.Poly(c,t)
  if p.degree()>0:roots.update(p.ground_roots())
 scores=[dict(t=str(r),nonlinear_total=sum(c.subs(t,r)!=0 for c in coeff)) for r in sorted(roots)];out=HERE/'orthogonal-line12';out.mkdir(exist_ok=True);payload=dict(status='complete',source_sha256=sha(Path(__file__)),input_sha256=sha(source),source_matrix=list(map(str,L)),inverse_matrix=list(map(str,inv)),variables=list(map(str,X)),generic_support=len(coeff),scores=scores,best=min(x['nonlinear_total'] for x in scores),coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in f.terms()] for f in rows],scope='All rational support-cancellation values on this joint polynomial-matrix linear conjugacy line. No arbitrary linear orbit or nonlinear carrier completion optimum.');(out/'checkpoint.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('GENERIC',len(coeff),'SCORES',scores,flush=True)
 if payload['best']<36:
  r=next(s.Rational(a['t']) for a in scores if a['nonlinear_total']==payload['best']);candidate=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameter=str(r),coefficients=[[dict(exponent=list(e),coefficient=str(c.subs(t,r))) for e,c in f.terms() if c.subs(t,r)] for f in rows],scope='Requires independent conjugacy/determinant/collision replay before acceptance.');(out/'candidate.json').write_text(json.dumps(candidate,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
