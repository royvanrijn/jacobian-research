#!/usr/bin/env python3
"""Exact rational translation-line search with tangent normalization of the cubic."""
from pathlib import Path
import hashlib,json,time
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/"translation-lines12";out.mkdir(exist_ok=True);source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text())
 X=s.symbols("x1:13");t=s.Symbol("t");F=s.Matrix([sum(s.Rational(v["coefficient"])*s.prod(x**e for x,e in zip(X,v["exponent"])) for v in row) for row in d["coefficients"]])
 J=F.jacobian(X);best=36;logs=[];started=time.monotonic();answer=None
 for i in range(12):
  centre=s.zeros(12,1);centre[i]=t
  point=dict(zip(X,centre));Jp=J.subs(point,simultaneous=True);inverse=Jp.inv()
  assert s.simplify(Jp.det())==1
  shifted=(inverse*(F.subs(dict(zip(X,s.Matrix(X)+centre)),simultaneous=True)-F.subs(point,simultaneous=True))).applyfunc(s.expand)
  assert shifted.jacobian(X).subs(dict.fromkeys(X,0))==s.eye(12)
  polys=[s.Poly(f,*X) for f in shifted];coefficients=[c for p in polys for e,c in p.terms() if sum(e)>=2]
  roots={s.Rational(-2),s.Rational(-1),s.Rational(-1,2),s.Rational(1,2),s.Rational(1),s.Rational(2)}
  for c in coefficients:
   poly=s.Poly(c,t)
   if poly.degree()>0:roots.update(poly.ground_roots())
  local=36
  for r in sorted(roots):
   if r==0:continue
   score=sum(c.subs(t,r)!=0 for c in coefficients);local=min(local,score)
   if score<best:
    best=score;new=[s.Poly(f.subs(t,r),*X) for f in shifted]
    answer=dict(status="translation-candidate",source_sha256=sha(Path(__file__)),input_sha256=sha(source),
      dimension=12,degree=3,nonlinear_total=score,coordinate=i,parameter=str(r),
      coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in p.terms()] for p in new],
      output_matrix=[[str(c.subs(t,r)) for c in row] for row in inverse.tolist()],
      centre=[str(c.subs(t,r)) for c in centre],
      scope="Exact primary translation and tangent normalization candidate, requiring independent determinant and collision verification.")
    (out/"candidate.json").write_text(json.dumps(answer,indent=2,sort_keys=True)+"\n");print("BEST",score,"coordinate",i,"parameter",r,flush=True)
  logs.append(dict(coordinate=i,roots=len(roots),line_minimum=local,best=best,seconds=round(time.monotonic()-started,3)))
  receipt=dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),best=best,lines=logs,
    scope="Rational coordinate translations of the36-term map, followed by inverse derivative target normalization. No arbitrary translations or linear frames.")
  (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(logs[-1]),flush=True)
 receipt["status"]="candidate" if answer else "completed-lines-no-improvement"
 (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
