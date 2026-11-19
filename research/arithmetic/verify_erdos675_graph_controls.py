#!/usr/bin/env python3
"""Exact finite controls for attempted graph extensions; no recurrence proof."""
if not __debug__:
    raise RuntimeError('Exact controls require assertions.')
from fractions import Fraction as Q
from math import isqrt,lcm
import json

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def vp(n,p):
    e=0
    while n%p==0:n//=p;e+=1
    return e
def norm(n):
    for p in range(2,isqrt(n)+1):
        if n%p:continue
        e=0
        while n%p==0:n//=p;e+=1
        if p%4==3 and e%2:return False
    return n%4!=3

assert norm(1+16) and not norm(1+16*17) and prime(17) and 17%16==1
d,q,p=43201,115201,576001
assert all(prime(x) and x%14400==1 for x in (d,q,p))
assert p==5*q-4 and p==14400*40+1
n=d*p;nn=d*5*q;R=d*q
assert nn-n==4*d
assert n==24883819201 and nn==24883992005 and R==4976798401
assert 5*q-3==576002==2*7*41143 and vp(576002,7)==1
assert 5*d-3==216002==59**2+461**2
assert not norm(5*q-3) and norm(5*d-3)
assert (Q(int(norm(5*q-3)))+Q(int(norm(5*d-3))))/2==Q(1,2)

rho=Q(1,100); k=4
c={(i,j):Q(i+j+1,3) for i in range(k) for j in range(i+1,k)}
law={():1-k*rho+sum(c.values())*rho*rho}
for i in range(k):law[(i,)]=rho-sum(v for ij,v in c.items() if i in ij)*rho*rho
for ij,v in c.items():law[ij]=v*rho*rho
assert all(v>=0 for v in law.values()) and sum(law.values())==1
for i in range(k):assert sum(v for s,v in law.items() if i in s)==rho
for ij,v in c.items():assert sum(w for s,w in law.items() if all(i in s for i in ij))==v*rho*rho
assert all(len(s)<=2 for s in law)

hadamard=0
for dim in (4,16):
    H=[[1]]
    while len(H)<dim:
        H=[r+r for r in H]+[r+[-x for x in r] for r in H]
    for i in range(dim):
        for j in range(dim):assert sum(H[i][a]*H[j][a] for a in range(dim))==(dim if i==j else 0)
    root=isqrt(dim);assert root*root==dim
    B=[[Q(H[i][j] if H[i][j]>0 else 0,root) for j in range(dim)] for i in range(dim)]
    out=[sum(r,Q(0))/root for r in B]
    assert sum(x*x for x in out)==Q(dim+3,4)>1
    hadamard+=1

atoms=0
for N in range(1,21):
    A=[a for a in range(1,N+1) if norm(a)];K=len(A);M=4*lcm(*range(1,N+1))**2
    for p in (3,7,11,19,23,31,43):
        if p<=N:continue
        mod=p*p
        bad=sum(any((1+(M//a)*u)%p==0 and (1+(M//a)*u)%mod!=0 for a in A) for u in range(mod))
        assert bad==K*(p-1),(N,p,bad,K)
        atoms+=1
print(json.dumps(dict(status='PASS',prime_bank_controls=3,
    raw_permutation_control='0 -> 1/2',sparse_probability_marginals=k,
    sparse_probability_pairs=len(c),hadamard_positive_schur_controls=hadamard,
    square_atom_root_counts=atoms,all_N_existence_tested=False),indent=2))
