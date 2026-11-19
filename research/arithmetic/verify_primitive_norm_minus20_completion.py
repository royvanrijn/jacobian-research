#!/usr/bin/env python3
"""Exact finite core/seed certificates supporting the full -20 sumset proof."""
if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")
import argparse
from fractions import Fraction as F
from math import isqrt,prod
from bisect import bisect_right
import json
import hashlib
from pathlib import Path


def primes_through(bound):
    flags=bytearray(b"\x01")*(bound+1)
    flags[:2]=b"\x00\x00"
    for p in range(2,isqrt(bound)+1):
        if flags[p]:
            flags[p*p:bound+1:p]=b"\x00"*((bound-p*p)//p+1)
    return [p for p in range(2,bound+1) if flags[p]]


def lattice_counts(bound):
    counts=[0]*(bound+1)
    even=[0]*(bound+1)
    odd=[0]*(bound+1)
    for y in range(-isqrt(bound//5),isqrt(bound//5)+1):
        for x in range(-isqrt(bound-5*y*y),isqrt(bound-5*y*y)+1):
            n=x*x+5*y*y
            counts[n]+=1
            if x%2 and y%2:odd[n]+=1
            if x%2==0 and y%2==0:even[n]+=1
    return counts,even,odd


def cube_floor(n):
    lo,hi=0,1
    while hi**3 <= n:hi*=2
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**3<=n:lo=mid
        else:hi=mid
    return lo


def maxsplit_envelope(maximum):
    split=[p for p in primes_through(1000) if p%20 in (1,3,7,9)]
    cumulative=1
    depth=0
    while cumulative<=maximum:
        cumulative*=split[depth]
        depth+=1
    # No feasible initial segment reaches the final listed prime.
    split=split[:depth]
    states=[]
    def visit(i,previous,cost,divisors):
        states.append((cost,divisors))
        if i==len(split):return
        power=1
        for a in range(1,previous+1):
            power*=split[i]
            if cost*power>maximum:break
            visit(i+1,a,cost*power,divisors*(a+1))
    visit(0,64,1,1)
    assert 3**65>maximum
    states.sort()
    costs=[];values=[];best=0
    for cost,d in states:
        costs.append(cost);best=max(best,d);values.append(best)
    def M(X):return values[bisect_right(costs,X)-1]
    return M, {'states':len(states),'split_primes_used':split,
               'max_bound':maximum,'first_excluded_prime_prefix_product':cumulative}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,
                        default=Path(__file__).with_name('primitive_norm_minus20_completion_checks.json'))
    args=parser.parse_args()
    r,_,oo=lattice_counts(35)
    cum=0
    for X in range(1,36):
        cum+=oo[X]
        assert 3*cum<=2*X
    assert [i for i,c in enumerate(oo) if c] == [6,14,30]
    assert all(oo[i]==4 for i in (6,14,30))
    assert F(99,280)+F(29,120)+F(1,36)<F(2,3)
    small=primes_through(63)
    inert=lambda p:p%20 in (11,13,17,19)
    sixth=prod(max(F((a+1)**6*p**(5*a),sum(p**j for j in range(a+1))**6)
                      for a in range(7)) for p in small if inert(p))
    sixth*=F(2**6*7**5,8**6)
    assert sixth<F(61,20)**6
    R=65536;D=2**64
    inert_primes=[p for p in primes_through(R) if inert(p)]
    tails=[0]*(len(inert_primes)+1)
    for i in range(len(inert_primes)-1,-1,-1):
        p=inert_primes[i]
        tails[i]=tails[i+1]+(D+p*p-1)//(p*p)
    B=sum(((2*D+p*(p-1)-1)//(p*(p-1)) for p in inert_primes),0)
    B=F(B,D)+F(5,2)*(F(1,5*R)+F(4,R*R))
    assert B<F(53,1000)
    M,env=maxsplit_envelope(2**52)
    records=[]
    for is_five in (False,True):
        C=F(427,50) if is_five else F(61,20)
        baseline=F(0) if is_five else F(1,5)
        fixederror=F(2) if is_five else F(24,5)
        for e in range(3):
            L=[3,3,2][e] if is_five else [6,6,4][e]
            for k in range(18,50):
                low=2**k;high=2**(k+1)-1
                err=C/cube_floor(low)
                mult=M(2**e*high)
                scale=F(L*mult,3)/err
                P=min(R-1,isqrt(scale.numerator//scale.denominator))
                j=bisect_right(inert_primes,P)
                tail=F(tails[j],D)+F(1,5*R)+F(4,R*R)
                lhs=B+baseline+err*(fixederror+3*j)+L*mult*tail
                records.append({'five_divisible':is_five,'e':e,'k':k,
                                'M':mult,'error_upper':str(err),'P':P,'prime_count':j,
                                'tail_upper':str(tail),'lhs':str(lhs),'lhs_decimal':float(lhs),
                                'passed':lhs<1})
    summary=[]
    for is_five in (False,True):
        for e in range(3):
            row=[r for r in records if r['five_divisible']==is_five and r['e']==e]
            first=next((r['k'] for i,r in enumerate(row) if all(v['passed'] for v in row[i:])),None)
            summary.append({'five_divisible':is_five,'e':e,'first_all_subsequent_k':first,
                            'first_margin':next((1-r['lhs_decimal'] for r in row if r['k']==first),None)})
    # Complete the infinite end without an M-envelope beyond its finite domain.
    assert 35**12 < 2**62
    infinite_core_unit = F(21,80)+F(1263,100)*F(2,35)+F(12447,50*2**16)
    infinite_core_five = F(1,16)+F(2631,250)*F(2,35)+F(4408,25*2**16)
    assert infinite_core_unit < 1 and infinite_core_five < 1
    starts={0:28,1:28,2:27}
    assert all(r['passed'] for r in records if r['k']>=starts[r['e']])
    # The exceptional cores3,4 propagate only through atoms of the genus-zero
    # split multiplier semigroup: a principal split prime or a pair of
    # nonprincipal split primes. Their seed targets have E<=8/sqrt(N).
    seeds=[]
    for e,L in [(0,6),(2,4)]:
        for k in range(24,30):
            low=2**k;high=2**(k+1)-1
            err=F(8,isqrt(low));mult=M(high)
            scale=F(L*mult,3)/err
            P=min(R-1,isqrt(scale.numerator//scale.denominator))
            j=bisect_right(inert_primes,P)
            tail=F(tails[j],D)+F(1,5*R)+F(4,R*R)
            lhs=B+F(1,5)+err*(F(24,5)+3*j)+L*mult*tail
            assert lhs<1
            seeds.append({'e':e,'k':k,'M':mult,'P':P,'lhs':str(lhs),'lhs_decimal':float(lhs)})
    infinite_seeds=F(21,80)+F(66,5*2**5)+F(1272,5*2**15)
    assert infinite_seeds<1
    claimed=[r for r in records if r['k']>=starts[r['e']]]
    assert len(claimed)==134
    result={'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'odd_odd_lattice':'S_oo(X)<=2X/3; finite table1..35 and row estimate above36',
            'core_E_sixth_product':str(sixth),'core_E_constant':'61/20',
            'five_core_E_constant':'427/50','B_upper':str(B),
            'split_envelope':env,'summary':summary,'core_interval_count':len(claimed),
            'intervals':claimed,
            'infinite_core_unit_upper_at_2_50':str(infinite_core_unit),
            'infinite_core_five_upper_at_2_50':str(infinite_core_five),
            'atomic_seed_intervals':seeds,'atomic_seed_infinite_upper_at_2_30':str(infinite_seeds),
            'scope':'Exact finite Euler, lattice, prime-tail and M-envelope support;134core and12atomic-seed interval inequalities. Written multiplier descent, local spectral application and infinite inequalities are separate proof arguments. No finite sumset scan is performed here.',
            'completion_dependency':'The complete2^29coordinate/genus/packed-sumset receipt is a separate indispensable finite computation.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('intervals','atomic_seed_intervals')},indent=2))


if __name__=='__main__':main()
