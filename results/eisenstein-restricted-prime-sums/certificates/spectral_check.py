"""Independent finite-group checks of the discriminant-3 Hecke normalization.

No theorem is inferred from the floating-point spectra; these checks detect
normalization, stabilizer, and orbit-counting mistakes in the proposed proof.
"""
from itertools import product
from collections import Counter
import numpy as np


def check(p):
    def canonical(g):
        g=tuple(x%p for x in g)
        u=pow(next(x for x in g if x),-1,p)
        return tuple(x*u%p for x in g)
    def rawmul(x,y):
        a,b,c,d=x;e,f,g,h=y
        return ((a*e+b*g)%p,(a*f+b*h)%p,(c*e+d*g)%p,(c*f+d*h)%p)
    def mul(x,y):return canonical(rawmul(x,y))
    def add(x,y):return tuple((a+b)%p for a,b in zip(x,y))
    def det(x):return (x[0]*x[3]-x[1]*x[2])%p
    I=(1,0,0,1);W=(0,p-1,1,p-1)
    a,b=next((a,b) for a,b in product(range(p),repeat=2) if (a*a+a*b+b*b+1)%p==0)
    J=(a,b,(a+b)%p,-a%p)
    assert rawmul(J,J)==tuple(-x%p for x in I)
    assert rawmul(J,W)==rawmul(add(tuple(-x%p for x in I),tuple(-x%p for x in W)),J)
    W2=rawmul(W,W)
    gamma={canonical(x) for x in [I,W,W2,J,rawmul(W,J),rawmul(W2,J)]}
    group={canonical(g) for g in product(range(p),repeat=4) if det(g)}
    assert len(gamma)==6 and len(group)==p*(p*p-1)
    def coset(g):return min(mul(h,g) for h in gamma)
    cosets=sorted({coset(g) for g in group})
    index={g:i for i,g in enumerate(cosets)}
    alphas=[canonical(add(I,rawmul(w,J))) for w in [I,W,W2]]
    matrix=np.zeros((len(cosets),len(cosets)),dtype=np.int64)
    for i,g in enumerate(cosets):
        for alpha in alphas:matrix[i,index[coset(mul(alpha,g))]]+=1
    assert np.array_equal(matrix,matrix.T)
    assert np.all(matrix.sum(axis=1)==3)
    char=np.array([1 if pow(det(g),(p-1)//2,p)==1 else -1 for g in cosets])
    leg2=1 if pow(2,(p-1)//2,p)==1 else -1
    assert np.array_equal(matrix@char,3*leg2*char)
    torus={canonical(add(tuple(a*x%p for x in I),tuple(b*x%p for x in W)))
           for a,b in product(range(p),repeat=2) if (a,b)!=(0,0)}
    normalizer=torus|{mul(x,J) for x in torus}
    orbit={coset(g) for g in normalizer}
    assert len(normalizer)==2*(p+1) and len(orbit)==(p+1)//3
    f=np.array([int(g in orbit) for g in cosets])
    assert f@char==0
    Tf=matrix@f
    assert np.all(Tf[f==1]==Tf[index[coset(I)]])
    assert f@Tf==len(orbit)*Tf[index[coset(I)]]
    vals=np.linalg.eigvalsh(matrix)
    nontriv=vals[np.abs(vals)<3-1e-8]
    assert np.max(np.abs(nontriv))<=2*np.sqrt(2)+1e-8
    print({'p':p,'size':len(cosets),'orbit':len(orbit),
           'T2_spectrum':dict(sorted(Counter(np.round(vals,8)).items())),
           'largest_nontrivial_abs':float(np.max(np.abs(nontriv))),
           'all_checks':'passed'})


for p in [5,11]:check(p)
