"""Exact finite controls for the primitive Weyl resultant and multiplicity-26 face."""

from fractions import Fraction as Fr
from itertools import permutations
from math import gcd

def add(out, p, c=1):
    for k,v in p.items():
        out[k]=out.get(k,0)+v*c
        if not out[k]:
            del out[k]

def mul(p,q):
    z={}
    for (i,j),v in p.items():
        add(z,{(i+k,j+l):w for (k,l),w in q.items()},v)
    return z

def determinant(M):
    n=len(M)
    out={}
    for sigma in permutations(range(n)):
        ans={(0,0):1}
        for i in range(n):
            ans=mul(ans,M[i][sigma[i]])
            if not ans:
                break
        add(out,ans,(-1)**sum(sigma[i]>sigma[j] for i in range(n) for j in range(i+1,n)))
    return out

def power(F,j):
    H={(0,0):1}
    for _ in range(j):
        H=mul(H,F)
    return H

num=0
for a in range(2,6):
    for b in range(2,12):
        if gcd(a,b)!=1:
            continue
        for seed in (1,2):
            gamma=seed+1
            G={r: ((r*r+7*seed)%11)-5 for r in range(1-a,b)}
            G[b]=gamma
            M=[[{} for _ in range(a)] for _ in range(a)]
            for i in range(a):
                M[i][i][0,1]=1
                for r,c in G.items():
                    p,q=divmod(i+r,a)
                    if c:
                        add(M[q][i],{(p,0):-c})
            H=determinant(M)
            assert max(a*i+b*j for i,j in H)==a*b
            assert {ij:c for ij,c in H.items() if a*ij[0]+b*ij[1]==a*b}=={(0,a):1,(b,0):-gamma**a}
            G0={(r,0):c for r,c in G.items() if c}
            out, deriv={},{}
            for (i,j),c in H.items():
                add(out,mul({(a*i,0):1},power(G0,j)),c)
                if j:
                    add(deriv,mul({(a*i,0):1},power(G0,j-1)),c*j)
            assert not out
            deg=(a-1)*b
            assert max(i for i,j in deriv)==deg
            assert deriv[deg,0]==a*gamma**(a-1)
            num+=1
assert num==50

# Exact high-root ODE arithmetic at three nonzero choices of the simple root parameter.
def dx(p):
    return {(i-1,j):i*c for (i,j),c in p.items() if i}
for beta in (Fr(1),Fr(-2),Fr(3,7)):
    r={(6,0):1,(5,0):beta}
    G=mul(mul({(1,0):1},{(1,0):1,(0,0):beta}),{(1,0):1,(0,0):-beta/2})
    G={i:-2*c/(15*beta**2) for i,c in G.items()}
    W=mul(r,G)
    lhs={}
    add(lhs,mul(W,dx(r)),15)
    add(lhs,mul(dx(W),r),-10)
    assert lhs==mul(r,r)
    assert max(i for i,j in W)==9
# Arithmetic of the boundary partitions after the e=1 exclusion.
partitions=[(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]
assert [xs for xs in partitions if gcd(6,26,*xs)==2]==[(4,),(2,2)]
print(f'{num} exact Laurent resultants, three high-root face identities, and the boundary partitions verified.')
