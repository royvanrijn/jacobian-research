#!/usr/bin/env python3
"""Exact rational checks of constants used in the odd-norm mixed argument.

These finite calculations support the stated monotone local-factor proof.
They do not prove the odd-filter spectral or 3-adic identities themselves.
"""
from fractions import Fraction
from pathlib import Path
import hashlib,json


def require(condition,message):
    if not condition: raise RuntimeError(message)


def tau_factor_fifth(p,e):
    return Fraction((e+1)*(p-1)*p**e,p**(e+1)-1)**5/p**e


def local_maxima(primes):
    rows=[];product=Fraction(1)
    for p in primes:
        end=0
        while (end+2)**5>p*(end+1)**5: end+=1
        e=max(range(end+1),key=lambda e:tau_factor_fifth(p,e))
        value=tau_factor_fifth(p,e);product*=value
        rows.append({'prime':p,'last_exponent_required':end,'maximizing_exponent':e,
                     'maximum_fifth_power':str(value)})
    return rows,product


def main():
    end=1
    while (end+4)**5>2*(end+3)**5: end+=1
    e=max(range(1,end+1),key=lambda e:Fraction((e+3)**5,3**5*2**e))
    two=Fraction((e+3)**5,3**5*2**e)
    rows,odd=local_maxima([5,11,17,23,29])
    combined=two*odd
    require(e==4,'Wrong2-adic maximum')
    require(combined<Fraction(9,4)**5,'Spectral constant too small')
    lattice=4*Fraction(15,8)
    require(lattice==Fraction(15,2),'Wrong lattice constant')
    small_leading=Fraction(9,4)*Fraction(4,15)
    small_error=Fraction(9,4)*Fraction(17,15)
    tail_leading=Fraction(15,2)*Fraction(1,15)
    tail_error=Fraction(15,2)*Fraction(2,5)
    require(small_leading+tail_leading==Fraction(11,10),'Wrong X^-1/5 coefficient')
    require(small_error+tail_error==Fraction(111,20),'Wrong X^-6/5 coefficient')
    require(Fraction(1,3)+Fraction(3,20)==Fraction(29,60),'Wrong main term')
    require(Fraction(5,2)**5<100,'Invalid fifth-root control')
    union=Fraction(29,60)+Fraction(11,10)*Fraction(2,5)+Fraction(111,20)*Fraction(1,250)+Fraction(2,3*100**2)
    require(union==Fraction(591,625)<1,'Union threshold does not close')
    receipt={'scope':'Exact local maxima and coefficient arithmetic; not a spectral or 3-adic proof',
             'two_adic_factor':'(e+3)/(3*2^(e/5)), e>=1',
             'two_adic_monotone_tail_start':end,'two_adic_maximizing_exponent':e,
             'two_adic_maximum_fifth_power':str(two),'odd_inert_local_maxima':rows,
             'odd_product_fifth_power':str(odd),'combined_fifth_power':str(combined),
             'spectral_bound_constant':'9/4','spectral_bound_fifth_power':str(Fraction(9,4)**5),
             'spectral_strict_comparison_verified':True,
             'lattice_constant':str(lattice),'union_main_constant':'29/60',
             'union_X_minus_one_fifth_constant':str(small_leading+tail_leading),
             'union_X_minus_six_fifths_constant':str(small_error+tail_error),
             'union_X_minus_two_constant':'2/3','cutoff':100000000,
             'strict_union_bound_at_cutoff':str(union),
             'large_prime_tail_reason':'For p>=32, (e+1)^5<=32^e<=p^e. Smaller odd inert local ratios decrease after the recorded cutoff.',
             'two_adic_tail_reason':'Successive ratios (e+4)/(e+3)/2^(1/5) decrease, so the recorded first ratio<=1 proves the global maximum.',
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(receipt,indent=2))


if __name__=='__main__': main()
