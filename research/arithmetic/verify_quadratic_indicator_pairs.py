#!/usr/bin/env python3
"""Finite exact tests for five class-number-one ordinary norm indicators.

This checks integer classification through LIMIT and p-adic residue enclosures.
It does not prove or numerically test any correlation asymptotic.
"""
if not __debug__:
    raise RuntimeError("Run without Python optimization: the exact checks require assertions.")

from fractions import Fraction as Q
from math import isqrt
import json

DS=(1,2,3,4,7)
DISCR={1:-4,2:-8,3:-3,4:-4,7:-7}
LIMIT=65536
SHIFTS=range(1,129)

def vp(n,p):
    assert n>0
    e=0
    while n%p==0:
        e+=1;n//=p
    return e,n

def chi(d,p):
    D=DISCR[d]
    if D%p==0:return 0
    if p==2:return 1 if D%8 in (1,7) else -1
    return 1 if pow(D%p,(p-1)//2,p)==1 else -1

def local_profile(d,p,n):
    e,u=vp(n,p)
    if p==2:
        if d in (1,4):return u%4==1 and (d!=4 or e!=1)
        if d==2:return u%8 in (1,3)
        if d==3:return e%2==0
        if d==7:return e!=1
    if d==3 and p==3:return u%3==1
    if d==7 and p==7:return pow(u%7,3,7)==1
    return chi(d,p)!=-1 or e%2==0

def certify(d,p,k,r):
    # None means the residue does not specify enough valuation/unit bits.
    if r==0:
        if p==2 and d==7 and k>=2:return True
        if chi(d,p)==1 and not (p==2 and d==7):return True
        return None
    e,u=vp(r,p)
    if p==2 and d==4 and e==1:return False
    conductor=2 if p==2 and d in (1,4) else 3 if p==2 and d==2 else 1 if (d,p) in ((3,3),(7,7)) else 0
    if e+conductor>k:return None
    return local_profile(d,p,r)

def norm_image(d,p,k):
    modulus=p**k
    squares={x*x%modulus for x in range(modulus)}
    bits=sum(1<<x for x in squares)
    mask=(1<<modulus)-1
    image=0
    for a in squares:
        shift=d*a%modulus
        image |= ((bits<<shift)|(bits>>(modulus-shift))) & mask if shift else bits
    return bytearray((image>>r)&1 for r in range(modulus))

def delta(d,p,h):
    a,_=vp(h,p)
    if p==2 and d==1:return Q(1,4) if a==0 else Q(1,2)-Q(3,2**(a+2))
    if p==2 and d==2:return Q(1,4) if a<=1 else Q(1,2)-Q(3,2**(a+1))
    if p==2 and d==4:
        if a==0:return Q(1,8)
        if a==1:return Q(0)
        if a==2:return Q(5,16)
        return Q(3,8)-Q(3,2**(a+2))
    if p==2 and d==7:return Q(1,2) if a<=1 else Q(3,4)
    if (d,p) in ((3,3),(7,7)):return Q(1,2)-Q(p+1,4*p**(a+1))
    return Q(1) if chi(d,p)==1 else Q(p,p+1)*(1-Q(1,p**(a+1)))

def single(d,p):
    if p==2 and d==4:return Q(3,8)
    if p==2 and d==7:return Q(3,4)
    if (p==2 and d in (1,2)) or (d,p) in ((3,3),(7,7)):return Q(1,2)
    return Q(1) if chi(d,p)==1 else Q(p,p+1)

def spf_table(limit):
    spf=list(range(limit+1))
    for p in range(2,isqrt(limit)+1):
        if spf[p]==p:
            for n in range(p*p,limit+1,p):
                if spf[n]==n:spf[n]=p
    return spf

def value_tests():
    spf=spf_table(LIMIT)
    counts={};tested=0
    for d in DS:
        actual=bytearray(LIMIT+1)
        for y in range(isqrt(LIMIT//d)+1):
            dy=d*y*y
            for x in range(isqrt(LIMIT-dy)+1):actual[x*x+dy]=1
        for n in range(1,LIMIT+1):
            m=n;pred=True;e2=0
            while m>1:
                p=spf[m];e,m=vp(m,p)
                if p==2:e2=e
                if chi(d,p)==-1 and e%2:pred=False
            if d in (4,7) and e2==1:pred=False
            assert bool(actual[n])==pred,(d,n,actual[n],pred)
            tested+=1
        counts[d]=sum(actual[1:])
    return tested,counts

def residue_tests():
    tests=singles=stable=0;moduli={}
    cases={1:(2,3,5),2:(2,5,7),3:(2,3,5),4:(2,3),7:(2,3,7)}
    for d,primes in cases.items():
        for p in primes:
            maxa=max(vp(h,p)[0] for h in SHIFTS)
            extra=3 if p==2 and d==2 else 2 if p==2 and d in (1,4) else 1
            depths=(maxa+extra+1,maxa+extra+2)
            moduli[f'{d}:{p}']=[p**k for k in depths]
            for k in depths:
                modulus=p**k
                image=norm_image(d,p,k)
                cert=[certify(d,p,k,r) for r in range(modulus)]
                for r,c in enumerate(cert):
                    if c is not None:
                        assert image[r]==c,(d,p,k,r,image[r],c)
                        stable+=1
                lower=Q(sum(c is True for c in cert),modulus)
                upper=Q(sum(image),modulus)
                assert lower<=single(d,p)<=upper,(d,p,k,'single',lower,single(d,p),upper)
                singles+=1
                for h in SHIFTS:
                    lo=Q(sum(cert[r] is True and cert[(r+h)%modulus] is True for r in range(modulus)),modulus)
                    hi=Q(sum(image[r] and image[(r+h)%modulus] for r in range(modulus)),modulus)
                    predicted=delta(d,p,h)
                    assert lo<=predicted<=hi,(d,p,k,h,lo,predicted,hi)
                    tests+=1
    return tests,singles,stable,moduli

def normalization_tests():
    c2={1:Q(4),2:Q(8),3:Q(4,3),4:Q(4),7:Q(28,9)}
    expected={1:Q(1,4),2:Q(1,2),3:Q(1,27),4:Q(1,16),7:Q(1,28)}
    bad={1:(2,),2:(2,),3:(3,),4:(2,),7:(2,7)}
    for d in DS:
        mass=Q(1)
        for p in bad[d]:mass*=delta(d,p,1)
        assert c2[d]*mass*mass==expected[d],(d,c2[d],mass,expected[d])
    return 5

def scope_guards():
    # Two is inert for chi_-3, but split with an order restriction for d7.
    assert chi(3,2)==-1 and not local_profile(3,2,2)
    assert chi(7,2)==1 and not local_profile(7,2,2)
    assert local_profile(7,2,4) and local_profile(7,2,8)
    # Abstract inert-parity support for chi_-20 contains3, while the
    # principal form x^2+5y^2 does not. General Euler and form scopes differ.
    assert pow((-20)%3,(3-1)//2,3)==1
    assert not any(x*x+5*y*y==3 for y in range(isqrt(3//5)+1)
                   for x in range(isqrt(3)+1))
    # Move inert2 for d3 from the good-prime product to the finite bad set.
    # The new raw prefactor is multiplied by3=(2+1)/(2-1).
    for h in SHIFTS:
        a,_=vp(h,2)
        good_correction=(1-Q(1,2**(a+1)))/(1-Q(1,2))
        assert good_correction==3*delta(3,2,h),h
    return dict(character_and_order_guards=4,
                nonprincipal_form_scope_guard=1,
                bad_set_enlargement_checks=len(SHIFTS))

if __name__=='__main__':
    values,counts=value_tests()
    pairs,singles,stable,moduli=residue_tests()
    norms=normalization_tests()
    guards=scope_guards()
    print(json.dumps(dict(status='PASS',integer_limit=LIMIT,integer_classification_checks=values,
          value_counts=counts,residue_pair_enclosure_checks=pairs,singleton_enclosure_checks=singles,
          certified_residue_memberships=stable,shift_one_normalization_checks=norms,
          moduli=moduli,scope_guards=guards,arithmetic='exact rational and integer',asymptotic_tested=False),indent=2))
