#!/usr/bin/env python3
"""Bounded nonmonotone elementary-conjugation beam from the36-term cubic."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,random,time
import sympy as sp
from verify import decode
from optimize12 import family,at,serialize
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(F):return tuple(tuple(sorted(f.items())) for f in F)
def support(F):return sum(sum(e)>=2 for f in F for e in f)
def main():
 out=HERE/"linear-beam12";out.mkdir(exist_ok=True);source=HERE/"cubic12-integer-witness.json"
 data=json.loads(source.read_text());F=[decode(r,12) for r in data["coefficients"]];assert support(F)==36
 rng=random.Random(261119);cache={};t=sp.Symbol("t");seen={key(F)};beam=[(F,[])];best=36;log=[];tested=0;started=time.monotonic();stopped=False
 for depth in range(1,9):
  candidates=[]
  for current,moves in beam:
   for i in range(12):
    for j in range(12):
     if i==j:continue
     if time.monotonic()-started>90:stopped=True;break
     fam=family(current,i,j);roots={Q(-2),Q(-1),Q(-1,2),Q(1,2),Q(1),Q(2)}
     for row in fam:
      for e,p in row.items():
       if sum(e)<2 or max(p)==0:continue
       leading=p[max(p)];signature=tuple(p.get(a,Q(0))/leading for a in range(max(p)+1))
       if signature not in cache:
        poly=sp.Poly(sum(sp.Rational(c.numerator,c.denominator)*t**a for a,c in enumerate(signature)),t)
        cache[signature]={Q(int(r.p),int(r.q)) for r in poly.ground_roots()}
       roots|=cache[signature]
     for root in sorted(roots-{Q(0)}):
      tested+=1;powers=[root**a for a in range(5)];score=0
      for row in fam:
       for e,p in row.items():
        if sum(e)>=2 and sum((c*powers[a] for a,c in p.items()),Q(0)):score+=1
       if score>best+4:break
      if score>best+4:continue
      new=at(fam,root);signature=key(new)
      if signature in seen:continue
      seen.add(signature)
      height=max(max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for row in new for c in row.values())
      if height>20:continue
      history=moves+[dict(i=i,j=j,t=str(root))];candidates.append((score,height,rng.random(),new,history))
      if score<best:
       best=score;payload=dict(status="exact-conjugation-candidate",source_sha256=sha(Path(__file__)),input_sha256=sha(source),
          dimension=12,degree=3,nonlinear_total=best,coefficients=serialize(new),moves=history,
          scope="Bounded beam candidate. No improved Keller witness accepted until independent exact conjugacy, determinant and collision replay.")
       (out/"candidate.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n");print("BEST",best,history,flush=True)
    if stopped:break
   if stopped:break
  row=dict(depth=depth,parents=len(beam),children=len(candidates),best=best,tested=tested,seconds=round(time.monotonic()-started,3));log.append(row)
  receipt=dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),best=best,log=log,tested=tested,
    beam_width=12,max_depth=8,timeout_seconds=90,seed=261119,max_support_over_best=4)
  (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(row),flush=True)
  if stopped or not candidates:break
  candidates.sort(key=lambda c:(c[0],c[1],c[2]));selected=candidates[:4]
  medium=[c for c in candidates[4:] if c[0]<=best+2];rng.shuffle(medium);selected+=medium[:4]
  rest=[c for c in candidates if all(c is not s for s in selected)];rng.shuffle(rest);selected+=rest[:12-len(selected)]
  beam=[(c[3],c[4]) for c in selected]
 receipt.update(status="candidate" if (out/"candidate.json").exists() else "bounded-no-improvement",seconds=round(time.monotonic()-started,3),
    distinct_maps=len(seen),coefficient_polynomial_cache=len(cache),
    scope="Finite seeded nonmonotone beam with rational cancellation roots and six generic values per elementary line, bounded height and support. Not exhaustive over affine orbits, translations or all conjugations.")
 (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt),flush=True)
if __name__=="__main__":main()
