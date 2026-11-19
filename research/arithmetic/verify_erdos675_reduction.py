#!/usr/bin/env python3
"""Finite exact diagnostics for the Erdős675 reduction; not an all-N proof."""
if not __debug__:
    raise RuntimeError('Exact checks require assertions; disable optimization.')
from fractions import Fraction as Q
from math import isqrt, lcm
import json

def primes_up_to(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]

def vp(n,p):
    assert n>0
    e=0
    while n%p==0:n//=p;e+=1
    return e

def norm_criterion(n):
    assert n>0
    for p in range(2,isqrt(n)+1):
        if n%p:continue
        e=0
        while n%p==0:n//=p;e+=1
        if p%4==3 and e%2:return False
    return n%4!=3

limit=4096
image={x*x+y*y for x in range(isqrt(limit)+1) for y in range(isqrt(limit)+1)
       if 0<x*x+y*y<=limit}
for n in range(1,limit+1):assert (n in image)==norm_criterion(n)

states=negative=dyadic=divisibility=coefficient_states=0
for N in range(1,101):
    inert=[p for p in primes_up_to(N) if p%4==3]
    M=4*lcm(*range(1,N+1))**2
    for a in range(1,N+1):
        assert M%a==0
        divisibility+=1
        for p in primes_up_to(N):
            assert (M//a)%p==0
            coefficient_states+=1
        witness=next((p for p in inert if vp(a,p)%2),None)
        assert (witness is None)==(a in image)
        for u in range(1,17):
            b=M*u+a
            for p in inert:
                assert vp(b,p)==vp(a,p)
                states+=1
            e=vp(a,2)
            assert vp(b,2)==e and ((b>>e)-(a>>e))%4==0
            dyadic+=1
            if witness:
                assert vp(b,witness)%2==1
                negative+=1

enclosures=0
for N in range(1,21):
    A=[a for a in range(1,N+1) if a in image]
    k=len(A)
    # Inert p>N is automatically coprime to all position differences.
    for p in (q for q in (3,7,11,19,23) if q>N):
        for depth in (2,3):
            modulus=p**depth
            # Scale M to1 modulo p^depth by the invertible u-variable map.
            # Undetermined zero residues are discarded for a lower bound.
            good=[]
            for u in range(modulus):
                vals=[(u+a)%modulus for a in A]
                good.append(all(v and vp(v,p)%2==0 for v in vals))
            lower=Q(sum(good),modulus)
            upper=lower+Q(k,modulus)
            predicted=1-Q(k,p+1)
            assert lower<=predicted<=upper,(N,p,depth,lower,predicted,upper)
            enclosures+=1

L=9*59**2
unit=(1,0)
period=0
while True:
    x,y=unit
    unit=((2*x+3*y)%L,(x+2*y)%L)
    period+=1
    if unit==(1,0):break
    assert period<10*L
assert period>0

def mul(a,b):return (a[0]*b[0]+3*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def power(a,n):
    value=(1,0)
    while n:
        if n%2:value=mul(value,a)
        a=mul(a,a);n//=2
    return value

pell_checks=0
for j in range(3):
    z,k=mul((5,3),power((2,1),period*j))
    assert z*z-3*k*k==-2
    assert (z%L,k%L)==(5,3)
    m=2*k*k-1;t=m*m-1
    assert t>0
    representations={1:(m,0),2:(m,1),4:(m-1,2*k),5:(m,2),8:(m-3,2*z)}
    for a,(x,y) in representations.items():
        assert t+a==x*x+y*y
        pell_checks+=1
    assert vp(t+3,3)==vp(t+6,3)==vp(t+7,59)==1
    pell_checks+=3

assert all(((288+a) in image)==(a in image) for a in range(1,9))
assert tuple(value%294 for value in power((2,1),168))==(1,0)
assert tuple(value%441 for value in power((2,1),504))==(1,0)
modular_pell=0
for j in range(1,4):
    z,k=mul((1,1),power((2,1),504*j))
    assert z*z-3*k*k==-2
    m=2*k*k-1;t=m*m-1
    assert t>0 and t%441==0
    for a,(x,y) in {1:(m,0),2:(m,1),4:(m-1,2*k),5:(m,2),8:(m-3,2*z)}.items():
        assert t+a==x*x+y*y
        modular_pell+=1
    assert vp(t+3,3)==vp(t+6,3)==vp(t+7,7)==1
    modular_pell+=3
preferred=0
for j in range(3):
    z,k=mul((1,1),power((2,1),28+168*j))
    assert z*z-3*k*k==-2
    assert k%2==1 and k%3==0 and (k*k-1)%49==0
    m=2*k*k-1;t=m*m-1
    assert t>0 and t%14112==0
    for a,(x,y) in {1:(m,0),2:(m,1),4:(m-1,2*k),5:(m,2),8:(m-3,2*z)}.items():
        assert t+a==x*x+y*y
        preferred+=1
    assert vp(t+3,3)==vp(t+6,3)==vp(t+7,7)==1
    preferred+=3
print(json.dumps(dict(status='PASS',two_square_classifications=limit,
                     small_inert_state_checks=states,negative_witness_checks=negative,
                     dyadic_state_checks=dyadic,lcm_affine_divisibility_checks=divisibility,
                     small_prime_unit_coefficient_checks=coefficient_states,
                     rational_local_density_enclosures=enclosures,
                     pell_identity_and_valuation_checks=pell_checks,pell_return_period=period,
                     modulus_14112_pell_checks=preferred,
                     modulus_441_pell_checks=modular_pell,
                     first_explicit_translation=288,arithmetic='exact integers and fractions',
                     all_N_recurrence_tested=False),indent=2))
