#!/usr/bin/env python3
"""Exact H3 root lifts over Z[tau], tau^2=tau+1, and residue fields.

No floating point algebraic numbers, approximate groups, or prime cutoff
is used in a completed all-characteristic certificate.
"""

from itertools import permutations

ZERO, ONE, TAU = (0,0), (1,0), (0,1)


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def neg(a):
    return -a[0],-a[1]


def mul(a,b):
    return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]


IDENTITY = tuple(ONE if i==j else ZERO for i in range(3) for j in range(3))


def matmul(a,b):
    return tuple(add(add(mul(a[3*i],b[j]),mul(a[3*i+1],b[3+j])),mul(a[3*i+2],b[6+j]))
                 for i in range(3) for j in range(3))


def matvec(a,v):
    return tuple(add(add(mul(a[3*i],v[0]),mul(a[3*i+1],v[1])),mul(a[3*i+2],v[2])) for i in range(3))


def permutation_product(a,b):
    return tuple(a[b[i]] for i in range(5))


def permutation_inverse(a):
    return tuple(a.index(i) for i in range(5))


def permutation_order(a):
    b,n=a,1
    while b != tuple(range(5)):
        b=permutation_product(b,a)
        n+=1
    return n


def h3_model():
    identity=tuple(range(5))
    double=[p for p in permutations(range(5)) if permutation_order(p)==2 and sum(p[i]!=i for i in range(5))==4]
    triple=next((a,b,c) for a in double for b in double for c in double
                if permutation_order(permutation_product(a,b))==5
                and permutation_order(permutation_product(b,c))==3
                and permutation_order(permutation_product(a,c))==2)
    cartan=[[(2,0),neg(TAU),ZERO],[neg(TAU),(2,0),(-1,0)],[ZERO,(-1,0),(2,0)]]
    simple=[]
    for i in range(3):
        matrix=list(IDENTITY)
        for j in range(3):
            matrix[3*i+j]=add(matrix[3*i+j],neg(cartan[i][j]))
        simple.append(tuple(matrix))
    model={(identity,0):IDENTITY}
    queue=[(identity,0)]
    while queue:
        key=queue.pop()
        for p,s in zip(triple,simple):
            target=(permutation_product(key[0],p),key[1]^1)
            matrix=matmul(model[key],s)
            if target in model:
                assert model[target]==matrix
            else:
                model[target]=matrix
                queue.append(target)
    assert len(model)==len(set(model.values()))==120
    roots={}
    for key,g in model.items():
        ginv=model[(permutation_inverse(key[0]),key[1])]
        assert matmul(g,ginv)==IDENTITY
        for i,p in enumerate(triple):
            reflection=(permutation_product(permutation_product(key[0],p),permutation_inverse(key[0])),1)
            root=tuple(g[3*j+i] for j in range(3))
            assert matmul(matmul(g,simple[i]),ginv)==model[reflection]
            assert matvec(model[reflection],root)==tuple(neg(v) for v in root)
            roots.setdefault(reflection,root)
    assert len(roots)==15
    return model,roots,{"simple_A5_components":triple,"simple_matrices":simple,"group_order":120,"reflections":15}


def root_matrix(matrices,roots,words):
    n=len(matrices)
    output=[]
    for word in words:
        block=[[ZERO]*n for _ in range(3)]
        prefix=IDENTITY
        for letter in word:
            i=abs(letter)-1
            # These affine root reflections are involutions, including in char2.
            image=matvec(prefix,roots[i])
            for row in range(3):
                block[row][i]=add(block[row][i],image[row])
            prefix=matmul(prefix,matrices[i])
        assert prefix==IDENTITY
        output.extend(block)
    return output


def integral_expansion(matrix):
    output=[]
    for row in matrix:
        first,second=[],[]
        for a,b in row:
            first.extend([a,b])
            second.extend([b,a+b])
        output.extend([first,second])
    return output


def factor_exact(n,trial_limit=1000000):
    n=abs(n)
    factors=[]
    d=2
    while d*d<=n:
        if d>trial_limit:
            raise ValueError("factorization budget exhausted; all-field audit INCOMPLETE")
        if n%d==0:
            factors.append(d)
            while n%d==0:
                n//=d
        d+=1
    if n>1:
        factors.append(n)
    return factors


class ResidueField:
    def __init__(self,prime,root):
        self.p=prime
        self.root=root

    def reduce(self,a):
        if self.root is None:
            return a[0]%self.p,a[1]%self.p
        return (a[0]+self.root*a[1])%self.p,0

    def add(self,a,b):
        return self.reduce(add(a,b))

    def neg(self,a):
        return self.reduce(neg(a))

    def mul(self,a,b):
        return self.reduce(mul(a,b))

    def inv(self,a):
        if self.root is not None:
            assert a[1]==0
            return pow(a[0],-1,self.p),0
        norm=(a[0]*a[0]+a[0]*a[1]-a[1]*a[1])%self.p
        scale=pow(norm,-1,self.p)
        return self.reduce(((a[0]+a[1])*scale,-a[1]*scale))


def residue_fields(prime,root_search_limit=1000000):
    if prime>root_search_limit:
        raise ValueError("residue-root budget exhausted; all-field audit INCOMPLETE")
    roots=[i for i in range(prime) if (i*i-i-1)%prime==0]
    if prime==5:
        assert roots==[3]
    else:
        assert len(roots) in [0,2]
    return [ResidueField(prime,r) for r in roots] if roots else [ResidueField(prime,None)]


def kernel(matrix,columns,field):
    rows=[[field.reduce(v) for v in row] for row in matrix]
    pivots=[]
    rank=0
    for col in range(columns):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]!=ZERO),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        scale=field.inv(rows[rank][col])
        rows[rank]=[field.mul(v,scale) for v in rows[rank]]
        for i in range(len(rows)):
            if i!=rank and rows[i][col]!=ZERO:
                scale=rows[i][col]
                rows[i]=[field.add(v,field.neg(field.mul(scale,w))) for v,w in zip(rows[i],rows[rank])]
        pivots.append(col)
        rank+=1
    basis=[]
    for free in range(columns):
        if free not in pivots:
            vector=[ZERO]*columns
            vector[free]=ONE
            for i,pivot in enumerate(pivots):
                vector[pivot]=field.neg(rows[i][free])
            basis.append(vector)
    return basis
