"""Constructive exact polarization, with a bounded union-size search."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path
from itertools import product
import json
import random
import sympy as s
from verify import X, ROOT, cubic12, homogeneous18, grouped18


def canonical(v):
    v = tuple(int(a) for a in v)
    g = reduce(gcd,v,0)
    if not g:
        return v,0
    g *= next(1 if a>0 else -1 for a in v if a)
    return tuple(a//g for a in v),g**3


def linear(v):
    return s.Add(*(a*x for a,x in zip(v,X)))


def sumforms(*forms):
    return tuple(map(sum,zip(*forms)))


def scaled(v,c):
    return tuple(c*a for a in v)


def decomposition(factors, mode=0):
    """Return a rational cube identity for a degree-three product."""
    out=defaultdict(Fraction)
    def add(v,c):
        v,scale=canonical(v)
        if scale:
            out[v]+=Fraction(c)*scale
    if factors[0]==factors[2]:
        add(factors[0],1)
    elif len(set(factors))==2:
        a=next(v for v in factors if factors.count(v)==2)
        b=next(v for v in factors if factors.count(v)==1)
        add(sumforms(a,b),Fraction(1,6))
        add(sumforms(a,scaled(b,-1)),Fraction(-1,6))
        add(b,Fraction(-1,3))
    elif mode==0:
        a,b,c=factors
        add(sumforms(a,b,c),Fraction(1,6))
        for u,v in [(a,b),(a,c),(b,c)]:
            add(sumforms(u,v),Fraction(-1,6))
        for v in factors:
            add(v,Fraction(1,6))
    elif mode==1:
        a,b,c=factors
        for sb,sc in product([1,-1],repeat=2):
            add(sumforms(a,scaled(b,sb),scaled(c,sc)),Fraction(sb*sc,24))
    else:
        # Seven-power formula on signed arguments; same product up to sign.
        a,b,c=factors
        sb,sc=[(1,-1),(-1,1),(-1,-1)][mode-2]
        out=defaultdict(Fraction,{v:co*sb*sc for v,co in decomposition([a,scaled(b,sb),scaled(c,sc)],0).items()})
    return {v:c for v,c in out.items() if c}


def terms(grouped=True, scale11=2, mask=511):
    _,k,_=cubic12(scale11)
    h,_=homogeneous18(k)
    if grouped:
        h,ys,forms=grouped18(h,scale11,mask)
    else:
        ys,forms=(),()
    variables=X+tuple(ys)
    vs=[tuple(int(s.expand(f).coeff(x)) for x in X) for f in list(X)+list(forms)]
    ts=[]
    for row,expr in enumerate(h+[ ]):
        for mon,co in s.Poly(s.expand(expr),*variables).terms():
            if not co:
                continue
            factors=[v for v,e in zip(vs,mon) for _ in range(e)]
            assert len(factors)==3
            ts.append((row,Fraction(co),factors))
    # Add t^3 as in the author's rank-completing trick.
    ts.append((17,Fraction(1),[vs[17]]*3))
    return ts


def build(ts,choices):
    rows=[defaultdict(Fraction) for _ in range(18)]
    for (row,co,factors),choice in zip(ts,choices):
        for v,c in decomposition(factors,choice).items():
            rows[row][v]+=co*c
    forms=sorted({v for row in rows for v,c in row.items() if c})
    return forms, [[row.get(v,Fraction()) for v in forms] for row in rows]


def search(ts,restarts=30,seed=261024889):
    rng=random.Random(seed)
    options=[[decomposition(f,m) for m in range(5)] if len(set(f))==3 else [decomposition(f)] for _,_,f in ts]
    best=None
    for restart in range(restarts):
        choices=[0 if restart==0 else rng.randrange(len(op)) for op in options]
        used=defaultdict(int)
        for op,ch in zip(options,choices):
            for v in op[ch]: used[v]+=1
        for sweep in range(12):
            changed=False
            order=list(range(len(ts))); rng.shuffle(order)
            for i in order:
                for v in options[i][choices[i]]: used[v]-=1
                costs=[sum(used[v]==0 for v in op) for op in options[i]]
                current=choices[i]
                candidates=[j for j,c in enumerate(costs) if c==min(costs)]
                choice=current if current in candidates else rng.choice(candidates)
                choices[i]=choice
                changed |= choice!=current
                for v in options[i][choice]: used[v]+=1
            if not changed: break
        forms,B=build(ts,choices)
        if best is None or len(forms)<best[0]:
            best=(len(forms),choices,forms,B)
            print('best',best[0],'restart',restart,flush=True)
    return best


def main():
    raw=terms(False)
    grouped=terms(True)
    baseline=len(build(raw,[0]*len(raw))[0])
    grouped_baseline=len(build(grouped,[0]*len(grouped))[0])
    print('baseline',baseline,'grouped',grouped_baseline,flush=True)
    assert baseline==104 and grouped_baseline==90
    best=search(grouped)
    n,choices,forms,B=best
    _,k,_=cubic12(2); h,_=homogeneous18(k); h[-1]=X[-1]**3
    polys=[s.Poly(linear(v)**3,*X) for v in forms]
    for i in range(18):
        reconstructed={}
        for co,p in zip(B[i],polys):
            for mon,c in p.terms():
                reconstructed[mon]=reconstructed.get(mon,0)+co*c
        assert s.Poly.from_dict({m:c for m,c in reconstructed.items() if c},X)==s.Poly(h[i],*X)
    assert s.Matrix(B).rank()==18
    payload={'baseline':baseline,'grouped_baseline':grouped_baseline,'size':n,
             'forms':forms,'B':[[str(v) for v in row] for row in B],
             'choices':choices,'seed':261024889,'restarts':30}
    (ROOT/'polarization.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('EXACT IDENTITY VERIFIED',n,flush=True)


if __name__=='__main__':main()
