"""Necessary arithmetic filters for the proved cyclic-axis torus-cusp theorem."""
from fractions import Fraction
from math import gcd


def ceil_fifth(value):
    return (value+4)//5


def local_gate(a,b,e,q,ell,delta,alpha,*,cyclic_axis_pair_certified=False,
               sole_critical_branch_certified=False,torus_cusp_certified=False):
    data=(a,b,e,q,ell,delta,alpha)
    if any(type(x) is not int or x<1 for x in data):
        raise ValueError('positive exact integer data required')
    if not 2<=a<b or gcd(a,b)!=1 or e<2 or ell<e*q:
        raise ValueError('coprime torus cusp and valid ramified local rank required')
    w=e*q;eta=ell-w
    source_hypotheses=(cyclic_axis_pair_certified and sole_critical_branch_certified
                       and torus_cusp_certified)
    divisibility=(delta%q==0 and w%delta==0 and eta%delta==0 and ell%delta==0)
    cone_allowed=(a%alpha==0 or b%alpha==0)
    chi=Fraction(ell)*(1-Fraction(1,a)-Fraction(1,b))
    balance=(chi==eta-Fraction(delta,alpha))
    universal=(w<=5*eta-2*delta)
    return dict(mass=w,excess=eta,regular_fibre_order=delta,
                base_degree=Fraction(ell,delta),unramified_boundary_count=Fraction(eta,delta),
                source_hypotheses_certified=source_hypotheses,
                necessary_arithmetic_passes=divisibility and cone_allowed and balance and universal,
                hypotheses_contradicted=source_hypotheses and not (divisibility and cone_allowed and balance and universal),
                exact_euler_balance=balance,cone_order_allowed=cone_allowed,
                fibre_order_multiple_of_local_residue=(delta%q==0),
                universal_excess_bound=universal,realization_claimed=False,
                arbitrary_finite_cover_excluded=False)


def universal_index_cap(degree):
    if type(degree) is not int or degree<7:
        raise ValueError('this current universal cap is stated for degree>=7')
    return degree-4 if degree<9 else (5*degree-17)//6


def high_index_gate(degree,index,prime_count=2):
    if any(type(x) is not int or x<1 for x in (degree,index,prime_count)):
        raise ValueError('positive integer degree/index/prime count required')
    high=2*index>degree
    excess=ceil_fifth(index+2)
    return dict(strict_high_index=high,minimum_own_cusp_excess=excess,
                necessary_budget_passes=(not high or (prime_count>=2 and index+excess+prime_count<=degree-1)),
                realization_claimed=False,source_smoothness_assumed=False)


def image_mass_gate(degree,primes,total_prime_count=None):
    if type(degree) is not int or degree<2 or not primes:
        raise ValueError('degree and a nonempty same-image prime list required')
    if any(len(p)!=2 or any(type(v) is not int or v<1 for v in p) for p in primes):
        raise ValueError('positive integer (index,residue) pairs required')
    t=len(primes);r=t+1 if total_prime_count is None else total_prime_count
    if type(r) is not int or r<t:
        raise ValueError('total prime count must retain the known primes')
    w=sum(e*f for e,f in primes);h=sum(ceil_fifth(f*(e+2)) for e,f in primes)
    high=2*w>degree
    return dict(mass=w,strict_high_mass=high,minimum_cusp_excess=h,
                necessary_budget_passes=(not high or (r>=t+1 and all(e>=2 for e,f in primes) and w+h+r-t<=degree-2)),
                unknown_primes_deleted=False,realization_claimed=False)
