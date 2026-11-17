"""Marked-root polynomial family and exact colored-sheet group generators.

Construction: Harish, as reproduced in Fonollosa arXiv:2609.33491v1.
The arithmetic color condition and block calculations are derived in
MARKED_WREATH_ACTIONS.md. No polynomial factorization is inferred from blocks.
"""
from itertools import combinations
from math import gcd
import sympy as sp

def polynomial_map(n):
    if n<3:raise ValueError("n must be at least3")
    x,y=sp.symbols("x y")
    zs=sp.symbols("z3:"+str(n+1));T=sp.Symbol("T");p=1+x*y
    z=dict(zip(range(3,n+1),zs))
    phi=sp.Poly(sp.expand((1-y*T)**(n-1)*(p*T-x)+(p*T-x)**2*
      (n*y*(y*T+sp.Rational(1,2))+z[3]*(p*T+x/2)+
       sum(z[j]*T**(j-2) for j in range(4,n+1)))),T)
    outputs=[phi.nth(j) for j in range(3,n+1)]+[phi.nth(2),-phi.nth(0)]
    return (x,y)+zs,T,phi,outputs

def parity(sigma):
    return sum(sigma[i]>sigma[j] for i in range(len(sigma)) for j in range(i+1,len(sigma)))%2

def cyclotomic_twist(n,u):
    # c^2=(-1)^(n(n-1)/2). For even n, c=1 or i; i is in mu_m when needed.
    return ((u-1)//2)%2 if n%4==2 else 0

def sheet_action(n,v,sigma=None,u=1):
    m=n-2;sigma=tuple(range(n)) if sigma is None else sigma
    return tuple(sigma[i]*m+(u*j+v[sigma[i]])%m for i in range(n) for j in range(m))

def colored_allowed(n,v,sigma,u=1):
    m=n-2
    return gcd(u,m)==1 and (m%2==1 or sum(v)%2==(parity(sigma)+cyclotomic_twist(n,u))%2)

def geometric_generators(n):
    m=n-2;zero=[0]*n;out=[]
    if m%2:
        v=zero.copy();v[0]=1;out.append(sheet_action(n,v))
    else:
        for i in range(n-1):
            v=zero.copy();v[i]=1;v[i+1]=-1;out.append(sheet_action(n,v))
        v=zero.copy();v[0]=2;out.append(sheet_action(n,v))
    for i in range(n-1):
        sigma=list(range(n));sigma[i],sigma[i+1]=sigma[i+1],sigma[i]
        v=zero.copy()
        if m%2==0:v[i]=1
        out.append(sheet_action(n,v,sigma))
    return out

def arithmetic_generators(n):
    m=n-2;out=geometric_generators(n)
    for u in range(1,max(m,2)):
        if gcd(u,m)!=1:continue
        sigma=list(range(n))
        if m%2==0 and cyclotomic_twist(n,u):sigma[0],sigma[1]=sigma[1],sigma[0]
        out.append(sheet_action(n,[0]*n,sigma,u))
    return out

def generated_block(generators,a,b):
    """Smallest invariant equivalence relation identifying a and b."""
    size=len(generators[0]);parent=list(range(size))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):
        i,j=find(i),find(j)
        if i==j:return False
        parent[j]=i;return True
    union(a,b)
    while True:
        changed=False;classes={}
        for i in range(size):classes.setdefault(find(i),[]).append(i)
        for indices in classes.values():
            base=indices[0]
            for g in generators:
                for i in indices[1:]:changed=union(g[base],g[i]) or changed
        if not changed:break
    return frozenset(i for i in range(size) if find(i)==find(a))

def permutation_cycles(g):
    todo=set(range(len(g)));out=[]
    while todo:
        first=next(iter(todo));i=first;cycle=[]
        while i in todo:
            todo.remove(i);cycle.append(i);i=g[i]
        out.append(tuple(cycle))
    return out
