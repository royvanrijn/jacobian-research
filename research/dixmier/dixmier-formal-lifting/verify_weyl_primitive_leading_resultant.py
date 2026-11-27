"""Exact finite controls for the primitive Weyl resultant and multiplicity-26 face."""

from fractions import Fraction as Fr
from itertools import permutations
from math import comb, factorial, gcd

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
    # The new (P,W) face has a nonzero first bracket and zero second bracket.
    p=power(r,2)
    chain={}
    add(chain,mul(W,dx(p)),15)
    add(chain,mul(dx(W),p),-20)
    assert chain=={ij:2*c for ij,c in power(r,3).items()}
    second={}
    add(second,mul(chain,dx(p)),30)
    add(second,mul(dx(chain),p),-20)
    assert not second
    # The first binomial resultant has nonzero differential degree 34.
    A=-Fr(2)/(15*beta**2)
    resultant=power(W,4)
    add(resultant,power(r,6),-A**4)
    inner=mul(power({(1,0):1,(0,0):beta},2),
              power({(1,0):1,(0,0):-beta/2},4))
    add(inner,{(6,0):1},-1)
    assert max(i for i,j in inner)==4
    assert inner[4,0]==-3*beta**2/2
    factored=mul({(24,0):A**4},
                 mul(power({(1,0):1,(0,0):beta},6),inner))
    assert resultant==factored
    assert min(i for i,j in resultant)==24
    assert max(i for i,j in resultant)==34
    assert resultant[34,0]==-3*A**4*beta**2/2
# Arithmetic of the boundary partitions after the e=1 exclusion.
partitions=[(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]
assert [xs for xs in partitions if gcd(6,26,*xs)==2]==[(4,),(2,2)]
assert all(gcd(3,*(m//2 for m in (26,*xs)))==1 for xs in [(4,),(2,2)])

# Universal normal ordering in Q[kappa,c]<P,Q>/(QP-PQ-kappa).
# Keys are the exponents of P, Q, kappa and c, in that order.
def weyl_mul(p,q):
    out={}
    for (i,j,k,l),v in p.items():
        for (a,b,c,d),w in q.items():
            for t in range(min(j,a)+1):
                coefficient=v*w*comb(j,t)*factorial(a)//factorial(a-t)
                add(out,{(i+a-t,j+b-t,k+c+t,l+d):coefficient})
    return out

def weyl_power(p,n):
    out={(0,0,0,0):1}
    for _ in range(n):
        out=weyl_mul(out,p)
    return out

def commutator(p,q):
    out=weyl_mul(p,q)
    add(out,weyl_mul(q,p),-1)
    return out

P={(1,0,0,0):1}
Q={(0,1,0,0):1}
K={(0,0,1,0):1}
C0={(0,0,0,1):1}
assert commutator(Q,P)==K
W=weyl_power(Q,2)
add(W,weyl_power(P,3),-1)
add(W,weyl_mul(C0,P),2)
C=commutator(W,P)
assert C=={ij:2*c for ij,c in weyl_mul(K,Q).items()}
assert commutator(C,P)=={(0,0,2,0):2}
squared_rhs=dict(W)
add(squared_rhs,weyl_power(P,3))
add(squared_rhs,weyl_mul(C0,P),-2)
squared_rhs=weyl_mul(weyl_power(K,2),squared_rhs)
assert weyl_power(C,2)=={ij:4*c for ij,c in squared_rhs.items()}
remaining={(2,0,0,0):3,(0,0,0,1):-2}
remaining=weyl_mul(weyl_power(K,2),remaining)
assert commutator(W,C)=={ij:2*c for ij,c in remaining.items()}
# A corrupted defining identity must fail even though its length-two chain survives.
bad_W=dict(W)
add(bad_W,weyl_mul(C0,P),-2)
assert commutator(commutator(bad_W,P),P)=={(0,0,2,0):2}
bad_rhs=dict(bad_W)
add(bad_rhs,weyl_power(P,3))
add(bad_rhs,weyl_mul(C0,P),-2)
bad_rhs=weyl_mul(weyl_power(K,2),bad_rhs)
assert weyl_power(C,2)!={ij:4*c for ij,c in bad_rhs.items()}
print(f'{num} exact Laurent resultants; three high-root faces, their length-two '
      'brackets and order-34 resultants; boundary partitions; universal Weyl '
      'chain identities and a corrupted-identity control verified.')
