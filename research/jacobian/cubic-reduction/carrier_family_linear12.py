#!/usr/bin/env python3
"""Greedy rational linear conjugacy from changed nonlinear carrier frames."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
import sympy as S
from verify import variable,add
from optimize12 import family,at,serialize
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def score(F):return sum(sum(e)>=2 for f in F for e in f)
def main():
 source=HERE/'carrier-reparam-wide/family.json';d=json.loads(source.read_text());ss,tt=S.symbols('s t');p,q=S.symbols('p q');u=S.Symbol('u');out=HERE/'carrier-family-linear12';out.mkdir(exist_ok=True);cache={};cases=[];globalbest=36;started=time.monotonic();tested=0;stop=False
 seeds=[(Q(-3,2),Q(0)),(Q(0),Q(-1,3)),(Q(-3,2),Q(-1,3)),(Q(1),Q(0)),(Q(0),Q(1)),(Q(1),Q(1))]
 for sv,tv in seeds:
  F=[variable(i,12) for i in range(12)]
  for a in d['low_coefficients']:
   v=Q(S.sympify(a['coefficient']).subs({p:1,q:1,ss:S.Rational(sv),tt:S.Rational(tv)}));
   if v:F[a['row']]=add(F[a['row']],{tuple(a['exponent']):v})
  initial=score(F);moves=[];steps=[]
  for step in range(6):
   baseline=score(F);best=(baseline,None,None,None);lines=0
   for i in range(12):
    for j in range(12):
     if i==j:continue
     if time.monotonic()-started>60:stop=True;break
     fam=family(F,i,j);roots=set()
     for row in fam:
      for e,poly in row.items():
       if sum(e)<2 or max(poly)==0:continue
       leading=poly[max(poly)];signature=tuple(poly.get(a,Q(0))/leading for a in range(max(poly)+1))
       if signature not in cache:
        sp=S.Poly(sum(S.Rational(c)*u**a for a,c in enumerate(signature)),u);cache[signature]={Q(r) for r in sp.ground_roots()}
       roots|=cache[signature]
     for root in roots-{Q(0)}:
      tested+=1;sc=sum(sum(e)>=2 and bool(sum((c*root**a for a,c in poly.items()),Q(0))) for row in fam for e,poly in row.items())
      if sc<best[0]:best=(sc,i,j,root)
     lines+=1
    if stop:break
   steps.append(dict(step=step,baseline=baseline,best=best[0],lines=lines))
   if best[1] is None:break
   sc,i,j,root=best;F=at(family(F,i,j),root);moves.append(dict(i=i,j=j,t=str(root)))
   if sc<globalbest:
    globalbest=sc;payload=dict(status='candidate-needs-independent-verification',source_sha256=sha(Path(__file__)),input_sha256=sha(source),seed=dict(p='1',q='1',s=str(sv),t=str(tv)),moves=moves,coefficients=serialize(F),nonlinear_total=sc,dimension=12,degree=3,scope='Composed nonlinear carrier change plus rational linear conjugacy candidate. Requires full independent composition, determinant/collision replay.');(out/'candidate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('BEST',sc,payload['seed'],moves,flush=True)
   if stop:break
  case=dict(seed=dict(s=str(sv),t=str(tv)),initial=initial,final=score(F),moves=moves,steps=steps);cases.append(case);receipt=dict(status='timeout' if stop else 'running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,global_best=globalbest,tested=tested,seconds=round(time.monotonic()-started,3),scope='Six fixed changed nonlinear carrier frames, each up to6 greedy rational elementary-linear conjugacy passes, with60-second global budget. All rational support-cancellation roots on completed lines; no simultaneous frames, arbitrary carrier deformations or global support optimum.');(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('SEED',sv,tv,'COUNTS',initial,score(F),flush=True)
  if stop:break
 if not stop:receipt['status']='complete'
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
