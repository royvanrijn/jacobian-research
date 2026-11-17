"""Exact finite controls for the central prime-cover Keller obstruction."""
from itertools import product

def sl2_prime(p):
    elems=[a for a in product(range(p),repeat=4) if (a[0]*a[3]-a[1]*a[2])%p==1]
    I=(1,0,0,1);Z=(p-1,0,0,p-1)
    def mul(a,b):
        return ((a[0]*b[0]+a[1]*b[2])%p,(a[0]*b[1]+a[1]*b[3])%p,
                (a[2]*b[0]+a[3]*b[2])%p,(a[2]*b[1]+a[3]*b[3])%p)
    neg=lambda a:tuple((-v)%p for v in a)
    canon=lambda a:min(a,neg(a))
    quotient=sorted({canon(a) for a in elems})
    qmul=lambda a,b:canon(mul(a,b))
    def act(a):
        out=[]
        for x in range(p):
            den=(a[2]*x+a[3])%p
            out.append(p if den==0 else ((a[0]*x+a[1])*pow(den,-1,p))%p)
        out.append(p if a[2]==0 else (a[0]*pow(a[2],-1,p))%p)
        return tuple(out)
    return dict(p=p,elements=elems,I=I,Z=Z,mul=mul,neg=neg,canon=canon,
      quotient=quotient,qI=canon(I),qmul=qmul,act=act)

def generated(gens,I,mul):
    out={I};todo=[I]
    while todo:
        a=todo.pop()
        for b in gens:
            c=mul(a,b)
            if c not in out:out.add(c);todo.append(c)
    return frozenset(out)

def gf2_rank(rows):
    basis={}
    for row in rows:
        while row:
            pivot=row.bit_length()-1
            if pivot not in basis:
                basis[pivot]=row;break
            row^=basis[pivot]
    return len(basis)

def h2_mod2(elems,I,mul):
    order=[I]+sorted(set(elems)-{I});lookup={a:i for i,a in enumerate(order)}
    n=len(order);table=[[lookup[mul(a,b)] for b in order] for a in order]
    pairs={(a,b):j for j,(a,b) in enumerate(product(range(1,n),repeat=2))}
    def cocycle_rows():
        for a,b,c in product(range(1,n),repeat=3):
            row=0
            for pair in [(b,c),(a,table[b][c]),(table[a][b],c),(a,b)]:
                if 0 not in pair:row^=1<<pairs[pair]
            if row:yield row
    def coboundary_rows():
        for a,b in pairs:
            row=0
            for c in [a,b,table[a][b]]:
                if c:row^=1<<(c-1)
            yield row
    r2=gf2_rank(cocycle_rows());r1=gf2_rank(coboundary_rows())
    return dict(order=n,normalized_C2_dimension=len(pairs),cocycle_constraint_rank=r2,
      coboundary_rank=r1,H2_dimension=len(pairs)-r2-r1)

def permutation_inverse(a):
    out=[0]*len(a)
    for i,j in enumerate(a):out[j]=i
    return tuple(out)

def permutation_subgroups(group):
    n=len(next(iter(group)));I=tuple(range(n));mul=lambda a,b:tuple(a[i] for i in b)
    found={frozenset([I])};todo=list(found)
    while todo:
        H=todo.pop()
        for a in group-set(H):
            K=generated(set(H)|{a},I,mul)
            if K not in found:found.add(K);todo.append(K)
    return found

def permutation_order(a):
    from math import lcm
    unseen=set(range(len(a)));lengths=[]
    while unseen:
        i=next(iter(unseen));j=i;n=0
        while j in unseen:
            unseen.remove(j);n+=1;j=a[j]
        lengths.append(n)
    return lcm(*lengths)
