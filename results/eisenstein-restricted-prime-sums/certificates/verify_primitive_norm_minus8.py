#!/usr/bin/env python3
"""Exact finite certificates supporting the primitive -8 norm-sumset proof.
Standard library only. The written spectral argument and finite reductions
supply the connection between these exact computations and the theorem.
"""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
from bisect import bisect_right
from fractions import Fraction as F
from math import gcd, isqrt, prod
from pathlib import Path
import hashlib
import json


def prime_flags(limit):
    flags=bytearray(b"\x01")*(limit+1)
    flags[:2]=b"\x00\x00"
    for p in range(2,isqrt(limit)+1):
        if flags[p]:
            flags[p*p:limit+1:p]=b"\x00"*((limit-p*p)//p+1)
    return flags


def primes_below(limit):
    return [p for p in range(2,limit)
            if all(p%d for d in range(2,isqrt(p)+1))]


def cubert_floor(n):
    lo=0
    hi=1<<((n.bit_length()+2)//3)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**3<=n:lo=mid
        else:hi=mid
    assert lo**3<=n<(lo+1)**3
    return lo


def sharp_bounds():
    inert=[p for p in primes_below(64) if p%8 in (5,7)]
    product=F(1)
    factors=[]
    assert 5**5<4**6
    for p in inert:
        factor,exponent=max((F((a+1)**6*p**(5*a),((p**(a+1)-1)//(p-1))**6),a)
                            for a in range(4))
        if factor>1:
            product*=factor
            factors.append([p,exponent,factor.numerator,factor.denominator])
    assert product<F(18,5)**6
    b_small=sum((F(2,p*(p-1)) for p in inert),F(0))
    b_upper=b_small+F(5,2)*(F(1,69**2)+F(1,8*69)+F(1,71**2)+F(1,8*71))
    assert b_upper<F(3,16)
    # Complete finite point-count ranges; Q>=x*x+y*y bounds the square box.
    bound=80
    counts=[0]*(bound+1)
    odd=[0]*(bound+1)
    even_primitive=[0]*(bound+1)
    radius=isqrt(bound)
    for x in range(-radius,radius+1):
        for y in range(-radius,radius+1):
            q=x*x+2*y*y
            if 0<q<=bound:
                counts[q]+=1
                if x%2:odd[q]+=1
                if x%2==0 and y%2:even_primitive[q]+=1
    all_cumulative=[];odd_cumulative=[];even_cumulative=[]
    a=o=e=0
    for n in range(1,bound+1):
        a+=counts[n];o+=odd[n];e+=even_primitive[n]
        all_cumulative.append(a);odd_cumulative.append(o);even_cumulative.append(e)
        assert 3*a<=8*n
        assert o<=2*n and e<=n
        assert e==sum(odd[1:n//2+1])
    assert all_cumulative[2]==8
    assert 7**2<2*5**2 and 10**2>2*7**2
    assert F(110,49)+F(24,7*9)==F(386,147)<F(8,3)
    assert F(55,49)+F(17,7*5)+F(1,25)<2
    ordinary=[]
    for name,predicate,constant in [('split',lambda p:p%8 in (1,3),6),
                                     ('odd_inert',lambda p:p%8 in (5,7),7)]:
        value=F(1)
        for p in primes_below(64):
            if predicate(p):
                value*=max(F((a+1)**6,p**a) for a in range(7))
        assert value<constant**6
        ordinary.append([name,constant,value.numerator,value.denominator])
    assert 8**6<3*7**6
    core_tail=F(3,16)+F(19,64)+F(69,2**21)
    assert core_tail<F(1,2)
    assert F(4096,297)<2**4
    assert 3**4*2**11>20**4
    prime_tail=F(3,16)+F(3,10)+F(7,2048)
    assert prime_tail<F(1,2)
    return dict(E_core_bound='(18/5)*m^(-1/3)',E_product_sixth_numerator=product.numerator,
                E_product_sixth_denominator=product.denominator,E_factor_columns=['prime','exponent','numerator','denominator'],
                E_maximizing_factors=factors,E_prime_cutoff=64,E_exponent_cutoff=3,
                B_upper_numerator=b_upper.numerator,B_upper_denominator=b_upper.denominator,B_simple_upper='3/16',
                all_point_bound='8X/3',all_point_finite_bound=80,all_point_counts=all_cumulative,
                all_point_analytic_from=81,all_point_analytic_coefficient='386/147',
                odd_point_bound='2X',odd_point_finite_bound=24,odd_point_counts=odd_cumulative[:24],
                odd_point_analytic_from=25,even_primitive_point_bound='X',
                ordinary_euler_columns=['prime_family','constant','product_numerator','product_denominator'],
                ordinary_euler_certificates=ordinary,inert_core_tail_onset=1<<60,
                core_tail_bound_numerator=core_tail.numerator,core_tail_bound_denominator=core_tail.denominator,
                prime_boundary_onset=1<<22,prime_tail_bound_numerator=prime_tail.numerator,
                prime_tail_bound_denominator=prime_tail.denominator)


def divisor_frontier(split_primes,upper):
    records=[]
    assert prod(split_primes[:13])>upper
    assert 3**40>upper
    def dfs(index,cap,n,value):
        records.append((n,value))
        p=split_primes[index]
        for exponent in range(1,cap+1):
            n*=p
            if n>upper:return
            dfs(index+1,exponent,n,value*(exponent+1))
    dfs(0,upper.bit_length(),1,1)
    records.sort()
    frontier=[];best=0
    for n,value in records:
        if value>best:
            frontier.append((n,value));best=value
    return frontier,len(records)


def interval_certificate():
    limit=1<<16
    flags=prime_flags(limit)
    primes=[p for p in range(2,limit+1) if flags[p]]
    assert [p for p in primes if p<=1000]==primes_below(1001)
    inert=[p for p in primes if p%8 in (5,7)]
    split=[p for p in primes if p%8 in (1,3)]
    scale=1<<64
    suffix=[0]*(len(inert)+1)
    for i in range(len(inert)-1,-1,-1):
        p=inert[i]
        suffix[i]=suffix[i+1]+(scale+p*p-1)//(p*p)
    tail=F(1,4*limit)+F(2,limit*limit)
    frontier,states=divisor_frontier(split,1<<63)
    keys=[n for n,value in frontier]
    running=1
    for n in range(1,10001):
        rem=n;value=1;p=2
        while p*p<=rem:
            exponent=0
            while rem%p==0:
                rem//=p;exponent+=1
            if exponent and p%8 in (1,3):value*=exponent+1
            p+=1
        if rem>1 and rem%8 in (1,3):value*=2
        running=max(running,value)
        assert running==frontier[bisect_right(keys,n)-1][1]
    cases=[]
    specs=[('odd_P+P',27,1,F(8,3),F(1,2),F(1,2)),
           ('C+C_at_2m',26,2,F(2),F(1,2),F(0)),
           ('C+C_at_4m',20,4,F(4,3),F(2,3),F(0))]
    for label,start,multiplier,L,baseline,dyadic_error in specs:
        rows=[];worst_gap=F(1);worst_k=None;worst_left=None
        for k in range(start,60):
            lower=1<<k;upper=(1<<(k+1))-1
            cube=cubert_floor(lower)
            E=F(18,5*cube)
            M=frontier[bisect_right(keys,multiplier*upper)-1][1]
            lm=L*M;ratio=lm/E
            cut=min(limit-1,isqrt(ratio.numerator//ratio.denominator))
            assert 1<=cut<limit
            count=bisect_right(inert,cut)
            T=F(suffix[count],scale)+tail
            left=F(3,16)+E*(count+dyadic_error)+lm*T
            gap=baseline-left
            assert gap>F(1,256),(label,k,left)
            rows.append([k,M,cube,cut,count,T.numerator,T.denominator,left.numerator,left.denominator])
            if gap<worst_gap:worst_gap=gap;worst_k=k;worst_left=left
        cases.append(dict(label=label,start_k=start,target_multiplier=multiplier,
                          L_numerator=L.numerator,L_denominator=L.denominator,
                          baseline_numerator=baseline.numerator,baseline_denominator=baseline.denominator,
                          dyadic_error_numerator=dyadic_error.numerator,dyadic_error_denominator=dyadic_error.denominator,
                          bins=rows,worst_k=worst_k,worst_numerator=worst_left.numerator,
                          worst_denominator=worst_left.denominator,minimum_margin_numerator=worst_gap.numerator,
                          minimum_margin_denominator=worst_gap.denominator))
    return dict(prime_limit=limit,prime_flags_sha256=hashlib.sha256(flags).hexdigest(),fixed_point_scale=scale,
                inert_prime_count=len(inert),frontier_upper=1<<63,frontier_states=states,frontier=frontier,
                frontier_regression_bound=10000,
                bin_columns=['k','M(multiplier*(2^(k+1)-1))','floor_cuberoot(2^k)','P','pi_I(P)',
                             'tail_numerator','tail_denominator','criterion_numerator','criterion_denominator'],
                cases=cases,uniform_margin='1/256')



def core_membership_flags(bound):
    """Exact C membership through bound, from complete forbidden-prime removal."""
    # Complete Eratosthenes sieve. Every forbidden prime at most bound
    # is removed; a larger prime cannot divide a positive value <=bound.
    prime = bytearray(b'\x01') * (bound + 1)
    prime[:2] = b'\x00\x00'
    for p in range(2, isqrt(bound) + 1):
        if prime[p]:
            prime[p*p::p] = b'\x00' * ((bound-p*p)//p + 1)
    core = bytearray(b'\x01') * (bound + 1)
    core[::2] = b'\x00' * (bound//2 + 1)
    for p in range(3, bound + 1, 2):
        if prime[p] and p % 8 in (5, 7):
            core[p::p] = b'\x00' * (bound//p)
    del prime
    return core


def minus8_finite_coverage(bound=1 << 27, short_bound=5000,
                          regression_bound=10000):
    from hashlib import sha256
    from math import gcd, isqrt

    assert bound >= regression_bound >= 121
    assert 1 <= short_bound <= bound

    core = core_membership_flags(bound)
    assert core[0] == 0 and core[1] == 1
    assert all(core[n] == 0 for n in (2,5,7,10,13,14,25,49))
    assert all(core[n] == 1 for n in (1,3,9,11,17,27,81,121))

    # Independent membership regression: direct primitive coordinates,
    # not another prime-factor sieve. Absolute values suffice for Q_8.
    direct = set()
    for y in range(isqrt(regression_bound//2) + 1):
        for x in range(isqrt(regression_bound-2*y*y) + 1):
            n = x*x + 2*y*y
            if n > 0 and n % 2 and gcd(x,y) == 1:
                direct.add(n)

    def prefix_agrees(flags):
        return all(bool(flags[n]) == (n in direct)
                   for n in range(regression_bound + 1))

    prefix = bytearray(core[:regression_bound + 1])
    assert prefix_agrees(prefix)
    corrupt_forbidden = bytearray(prefix)
    corrupt_forbidden[5] = 1
    assert not prefix_agrees(corrupt_forbidden)
    corrupt_missing = bytearray(prefix)
    corrupt_missing[1] = 0
    assert not prefix_agrees(corrupt_missing)

    size = bound//8 + 1
    packed = bytearray(size)
    packed_double = bytearray(size)
    for n in range(1, bound + 1, 2):
        if core[n]:
            packed[n >> 3] |= 1 << (n & 7)
            if 2*n <= bound:
                packed_double[(2*n) >> 3] |= 1 << ((2*n) & 7)
    bits = int.from_bytes(packed, 'little')
    double_bits = int.from_bytes(packed_double, 'little')
    shorts = [n for n in range(1, short_bound + 1, 2) if core[n]]
    double_shorts = [2*n for n in range(1, short_bound//2 + 1, 2)
                     if core[n]]

    # All shifts have an actual positive core summand. The second odd
    # coverage family also permits the doubled summand to be short.
    even_cover = 0
    odd_cover = 0
    for a in shorts:
        even_cover |= bits << a
        odd_cover |= double_bits << a
    for twice_b in double_shorts:
        odd_cover |= bits << twice_b
    limit_mask = (1 << (bound + 1)) - 1
    even_cover &= limit_mask
    odd_cover &= limit_mask

    # 0x54 has precisely residues 2,4,6 modulo eight. 0xaa contains
    # odd residues; delete the inadmissible minimum one separately.
    expected_even_bytes = bytearray(b'\x54') * size
    expected_odd_bytes = bytearray(b'\xaa') * size
    expected_odd_bytes[0] &= 0xfd
    final_mask = (1 << ((bound % 8) + 1)) - 1
    expected_even_bytes[-1] &= final_mask
    expected_odd_bytes[-1] &= final_mask
    expected_even = int.from_bytes(expected_even_bytes, 'little')
    expected_odd = int.from_bytes(expected_odd_bytes, 'little')
    missing_even = expected_even & ~even_cover
    missing_odd = expected_odd & ~odd_cover
    forbidden_even = even_cover & ~expected_even
    forbidden_odd = odd_cover & ~expected_odd

    def first_positions(integer, count=20):
        out = []
        while integer and len(out) < count:
            low = integer & -integer
            out.append(low.bit_length()-1)
            integer ^= low
        return out

    assert not missing_even, ('C+C missing', missing_even.bit_count(),
                              first_positions(missing_even))
    assert not missing_odd, ('C+2C missing', missing_odd.bit_count(),
                             first_positions(missing_odd))
    assert not forbidden_even, ('C+C forbidden', first_positions(forbidden_even))
    assert not forbidden_odd, ('C+2C forbidden', first_positions(forbidden_odd))

    # Negative certificate checks: removing the least eligible target,
    # or injecting an excluded target, must be detected by these masks.
    assert expected_even & ~(even_cover & ~(1 << 2)) == 1 << 2
    assert expected_odd & ~(odd_cover & ~(1 << 3)) == 1 << 3
    assert (even_cover | 1) & ~expected_even
    assert (odd_cover | (1 << 1)) & ~expected_odd

    even_bytes = even_cover.to_bytes(size, 'little')
    odd_bytes = odd_cover.to_bytes(size, 'little')
    return {
        'scope': 'Exact finite strong-core coverage; unrestricted claims need the separate analytic cutoff and descent proof.',
        'bound': bound,
        'short_bound': short_bound,
        'short_core_summands': len(shorts),
        'short_doubled_core_summands': len(double_shorts),
        'C_plus_C_targets': expected_even.bit_count(),
        'C_plus_2C_targets': expected_odd.bit_count(),
        'C_plus_C_missing': 0,
        'C_plus_2C_missing': 0,
        'C_plus_C_forbidden_coverage': 0,
        'C_plus_2C_forbidden_coverage': 0,
        'independent_coordinate_membership_regression_bound': regression_bound,
        'corruption_checks': {
            'forbidden_core_prime_injection_detected': True,
            'missing_core_unit_detected': True,
            'missing_least_even_target_detected': True,
            'missing_least_odd_target_detected': True,
            'forbidden_even_target_injection_detected': True,
            'forbidden_odd_target_injection_detected': True,
        },
        'core_flags_sha256': sha256(core).hexdigest(),
        'packed_core_sha256': sha256(packed).hexdigest(),
        'packed_doubled_core_sha256': sha256(packed_double).hexdigest(),
        'even_coverage_sha256': sha256(even_bytes).hexdigest(),
        'odd_coverage_sha256': sha256(odd_bytes).hexdigest(),
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus8_checks.json'))
    args=parser.parse_args()
    result=dict(sharp_bounds=sharp_bounds(),interval_certificate=interval_certificate(),
                finite_coverage=minus8_finite_coverage(),
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Exact arithmetic certificates. The written theta/spectral sieve, dyadic reductions and split-prime descent establish the unrestricted theorem.')
    out=json.dumps(result,indent=2)+'\n'
    args.output.write_text(out)
    print(out,end='')


if __name__=='__main__':main()
