"""Exact bounded literal graph deletion plus one target shear stage."""
from __future__ import annotations
from itertools import combinations_with_replacement
import json
import time
import sympy as s
from flint import fmpq_mat
from verify import X, ROOT, cubic12


def high(f):
    return {m:c for m,c in s.Poly(s.expand(f),*X[:12]).terms() if sum(m)>3}


def rank_pair(cols, target):
    mons=sorted(set(target).union(*(c.keys() for c in cols)))
    rows=[[str(c.get(m,0)) for c in cols]+[str(target.get(m,0))] for m in mons]
    matrix=fmpq_mat(rows)
    rr,aug=matrix.rref()
    pivots=[]
    for i in range(aug):
        pivots.append(next(j for j in range(len(cols)+1) if rr[i,j]))
    rank=aug-(len(cols) in pivots)
    solution=None
    if rank==aug:
        solution=[s.Rational(0)]*len(cols)
        for i,p in enumerate(pivots):solution[p]=s.Rational(str(rr[i,len(cols)]))
    return rank,aug,solution,len(mons)


def main():
    _,k,_=cubic12(2)
    results=[]
    for j in range(3,12):
        start=time.monotonic()
        graph=X[j]-k[j]
        assert not graph.has(X[j])
        kept={i:s.expand(f.subs(X[j],graph)) for i,f in enumerate(k) if i!=j}
        bad=[i for i,f in kept.items() if high(f)]
        for bound in [2,3]:
            allok=True
            for i in bad:
                allowed=[a for a in kept if a!=i]
                indices=[ix for d in range(1,bound+1) for ix in combinations_with_replacement(allowed,d)]
                cols=[high(s.prod(kept[a] for a in ix)) for ix in indices]
                r,aug,solution,nrows=rank_pair(cols,high(kept[i]))
                row={'pivot':j+1,'component':i+1,'target_degree':bound,
                     'columns':len(cols),'rows':nrows,'rank':r,'augmented_rank':aug}
                if solution is not None:
                    row['solution']=[{'indices':[a+1 for a in ix],'coefficient':str(co)} for ix,co in zip(indices,solution) if co]
                    expr=kept[i]-sum(co*s.prod(kept[a] for a in ix) for ix,co in zip(indices,solution) if co)
                    assert not high(expr)
                else:allok=False
                results.append(row)
                print(row if solution is None else {**row,'solution':'found'},flush=True)
                (ROOT/'graph_audit.json').write_text(json.dumps(results,indent=2)+'\n')
                if solution is None:break
            if allok:
                print('CANDIDATE',j+1,bound,flush=True)
                break
        print('pivot_seconds',j+1,round(time.monotonic()-start,3),flush=True)


if __name__=='__main__':main()
