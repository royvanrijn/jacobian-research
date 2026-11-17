#!/usr/bin/env python3
"""Exact rational elementary-conjugation descent for the cubic witness."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json
import sympy as sp
from verify import add,decode
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def family(F,i,j):
    # Rows indexed by (source exponent, parameter exponent).
    rows=[]
    for f in F:
        out={}
        for e,c in f.items():
            for a in range(e[i]+1):
                ee=list(e);ee[i]-=a;ee[j]+=a;key=(tuple(ee),a)
                out[key]=out.get(key,Q(0))+c*comb(e[i],a)
        rows.append({key:c for key,c in out.items() if c})
    result=dict(rows[i])
    for (e,a),c in rows[j].items():result[(e,a+1)]=result.get((e,a+1),Q(0))-c
    rows[i]={key:c for key,c in result.items() if c}
    grouped=[]
    for row in rows:
        by_e={}
        for (e,a),c in row.items():by_e.setdefault(e,{})[a]=c
        grouped.append(by_e)
    return grouped
def at(rows,t):
    return [{e:sum((c*t**a for a,c in p.items()),Q(0)) for e,p in row.items()
             if sum((c*t**a for a,c in p.items()),Q(0))} for row in rows]
def serialize(F):return [[dict(exponent=list(e),coefficient=str(c)) for e,c in sorted(f.items())] for f in F]
def main():
    path=HERE/"de-bondt12-witness.json";d=json.loads(path.read_text());n=d["dimension"]
    F=[decode(row,n) for row in d["coefficients"]];moves=[];log=[];cache={};t=sp.Symbol("t")
    initial=sum(sum(e)>=2 for f in F for e in f)
    for step in range(12):
        baseline=sum(sum(e)>=2 for f in F for e in f);best=(baseline,None,None,None);checked=0
        for i in range(n):
            for j in range(n):
                if i==j:continue
                fam=family(F,i,j);candidates={Q(0)}
                for row in fam:
                    for e,p in row.items():
                        if sum(e)<2 or max(p)==0:continue
                        leading=p[max(p)];key=tuple(p.get(a,Q(0))/leading for a in range(max(p)+1))
                        if key not in cache:
                            poly=sp.Poly(sum(sp.Rational(c.numerator,c.denominator)*t**a for a,c in enumerate(key)),t)
                            cache[key]={Q(int(r.p),int(r.q)) for r in poly.ground_roots()}
                        candidates|=cache[key]
                for root in candidates:
                    if not root:continue
                    score=sum(sum(e)>=2 and bool(sum((c*root**a for a,c in p.items()),Q(0))) for row in fam for e,p in row.items())
                    if score<best[0]:best=(score,i,j,root)
                checked+=1
        row=dict(step=step,baseline=baseline,lines=checked,best=best[0],root=None if best[3] is None else str(best[3]));log.append(row);print(json.dumps(row),flush=True)
        if best[1] is None:break
        score,i,j,root=best;F=at(family(F,i,j),root);moves.append(dict(i=i,j=j,t=str(root)))
        witness=dict(status="rational-descent-candidate",source_sha256=sha(Path(__file__)),input_sha256=sha(path),
          dimension=n,degree=max(sum(e) for f in F for e in f),nonlinear_total=score,
          coefficients=serialize(F),moves=moves,log=log,
          scope="Greedy rational elementary conjugations only; no global optimum claim.")
        (HERE/"de-bondt12-sparse-candidate.json").write_text(json.dumps(witness,indent=2,sort_keys=True)+"\n")
    receipt=dict(status="complete",source_sha256=sha(Path(__file__)),input_sha256=sha(path),
      initial_support=initial,final_support=sum(sum(e)>=2 for f in F for e in f),moves=moves,log=log,
      distinct_coefficient_polynomials=len(cache),boundary="All rational cancellation values on each tested elementary line; greedy descent with at most12 iterations, not full affine-orbit coverage.")
    (HERE/"linear-descent12.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

