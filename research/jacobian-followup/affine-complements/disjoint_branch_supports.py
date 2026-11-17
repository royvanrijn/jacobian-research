"""Small exact controls for PF121; the proof is unbounded and written."""
from itertools import combinations


def subset_signatures(m,k,i,j):
    assert 2<=k<=m-2 and i!=j
    family=[frozenset((i,)+z) for z in combinations([x for x in range(m) if x not in (i,j)],k-1)]
    return family,[tuple(x in A for A in family) for x in range(m)]


def polar(x,y,a):
    mask=(1<<a)-1
    return (((x&mask)&(y>>a)).bit_count()+((y&mask)&(x>>a)).bit_count())%2


def quadratic(x,a,coeff=0):
    mask=(1<<a)-1
    return (((x&mask)&(x>>a)).bit_count()+(x&coeff).bit_count())%2


def arf(coeff,a):
    mask=(1<<a)-1
    return ((coeff&mask)&(coeff>>a)).bit_count()%2


def transvection(x,v,a):
    return x^(v if polar(x,v,a) else 0)


def minus_forms(a):
    return [c for c in range(1<<(2*a)) if arf(c,a)==1]


def form_permutation(a,v):
    forms=minus_forms(a);where={q:i for i,q in enumerate(forms)}
    # Pullback by the involution: q(t_v(x))=q(x)+(q(v)+1)B(x,v).
    dual=(v>>a)|((v&((1<<a)-1))<<a)
    return [where[c^(dual if quadratic(v,a,c)==0 else 0)] for c in forms]


def vector_permutation(a,negative,epsilon,v):
    coeff=(1|(1<<a)) if negative else 0
    assert quadratic(v,a,coeff)==1
    points=[x for x in range(1<<(2*a)) if x and quadratic(x,a,coeff)==epsilon]
    where={x:i for i,x in enumerate(points)}
    return points,[where[transvection(x,v,a)] for x in points]


def support(p):
    return {i for i,j in enumerate(p) if i!=j}
