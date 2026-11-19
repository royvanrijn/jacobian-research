#!/usr/bin/env python3
"""Exact finite certificates for the -11 theta identity and divisor bounds.
The checker proves finite integer/rational statements; the accompanying
proof gives the Sturm and Euler completeness arguments. Standard library.
"""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
import hashlib
import json
from bisect import bisect_right
from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path


def primes_below(n):
    return [p for p in range(2,n)
            if all(p%d for d in range(2,isqrt(p)+1))]


def sigma(n):
    return sum(d for d in range(1,n+1) if n%d==0) if n>0 else 0


def binary_counts(bound):
    counts=[0]*(bound+1)
    for y in range(-isqrt(4*bound//11),isqrt(4*bound//11)+1):
        s=isqrt(4*bound-11*y*y)
        for x in range(-((s+y)//2),(s-y)//2+1):
            q=x*x+x*y+3*y*y
            assert 0<=q<=bound
            counts[q]+=1
    return counts


def eta11(bound):
    # eta(tau)^2 eta(11*tau)^2 = q times this product.
    f=[1]+[0]*bound
    for step in (1,11):
        for j in range(step,bound+1,step):
            for repeat in range(2):
                for n in range(bound,j-1,-1):
                    f[n]-=f[n-j]
    return [0]+f[:bound]


def euler_certificate(name,predicate,ceiling):
    product=Fraction(1)
    factors=[]
    for p in primes_below(64):
        if not predicate(p):
            continue
        a=max(range(11),key=lambda a:Fraction((a+1)**6,p**a))
        f=Fraction((a+1)**6,p**a)
        if f>1:
            product*=f
            factors.append(dict(prime=p,exponent=a,
                                factor_numerator=f.numerator,
                                factor_denominator=f.denominator))
    assert product<ceiling**6
    return dict(name=name,constant=ceiling,power=6,
                rational_product_numerator=product.numerator,
                rational_product_denominator=product.denominator,
                maximizing_factors=factors)


def prime_flags(limit):
    flags=bytearray(b"\x01")*(limit+1)
    flags[:2]=b"\x00\x00"
    for p in range(2,isqrt(limit)+1):
        if flags[p]:
            flags[p*p:limit+1:p]=b"\x00"*((limit-p*p)//p+1)
    return flags


def sharp_core_bounds():
    product=Fraction(1)
    factors=[]
    inert=[p for p in primes_below(64) if p>2 and p%11 not in (0,1,3,4,5,9)]
    assert 4**6<7*3**6
    for p in inert:
        factor,exponent=max((Fraction((a+1)**6*p**(5*a),((p**(a+1)-1)//(p-1))**6),a)
                            for a in range(3))
        if factor>1:
            product*=factor
            factors.append([p,exponent,factor.numerator,factor.denominator])
    assert product<Fraction(11,4)**6
    counts=binary_counts(48)
    cumulative=[]
    total=0
    for n in range(1,49):
        total+=counts[n]
        cumulative.append(total)
        assert 5*total<=12*n
    assert cumulative[4]==12
    assert 33**2<100*11
    assert Fraction(440,231)+Fraction(106,33*7)==Fraction(26,11)<Fraction(12,5)
    return dict(E_core_bound='(11/4)*n^(-1/3)',
                E_product_sixth_numerator=product.numerator,
                E_product_sixth_denominator=product.denominator,
                inert_primes_below_64=inert,
                E_factor_columns=['prime','exponent','numerator','denominator'],
                E_maximizing_factors=factors,E_exponent_cutoff=2,
                point_count_bound='12X/5',point_count_finite_bound=48,
                point_count_cumulative_nonzero=cumulative,
                point_count_analytic_tail_from=49,
                point_count_analytic_tail_coefficient='26/11',
                large_prime_sieve_coefficient='8/3')


def cubert_floor(n):
    lo=0
    hi=1<<((n.bit_length()+2)//3)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**3<=n:
            lo=mid
        else:
            hi=mid
    assert lo**3<=n<(lo+1)**3
    return lo


def divisor_frontier(split_primes,upper):
    records=[]
    def dfs(index,cap,n,value):
        records.append((n,value))
        if index>=len(split_primes):
            return
        p=split_primes[index]
        for a in range(1,cap+1):
            n*=p
            if n>upper:
                return
            dfs(index+1,a,n,value*(a+1))
    # Sorted exponent patterns on a prime prefix suffice by rearrangement.
    prefix=1
    for p in split_primes:
        prefix*=p
        if prefix>upper:
            break
    assert prefix>upper
    dfs(0,upper.bit_length(),1,1)
    records.sort()
    frontier=[]
    best=0
    for n,value in records:
        if value>best:
            frontier.append((n,value))
            best=value
    return frontier,len(records)


def core_interval_certificate():
    limit=1<<16
    flags=prime_flags(limit)
    prime_list=[p for p in range(2,limit+1) if flags[p]]
    assert [p for p in prime_list if p<=1000]==primes_below(1001)
    inert=[p for p in prime_list if p>2 and p%11 not in (0,1,3,4,5,9)]
    split=[p for p in prime_list if p%11 in (1,3,4,5,9)]
    assert prod(split[:12])>1<<61
    assert 3**41>1<<61
    scale=1<<64
    suffix=[0]*(len(inert)+1)
    for i in range(len(inert)-1,-1,-1):
        p=inert[i]
        suffix[i]=suffix[i+1]+(scale+p*p-1)//(p*p)
    tail=Fraction(5,22*limit)+Fraction(5,limit*limit)
    frontier,states=divisor_frontier(split,1<<61)
    keys=[n for n,value in frontier]
    # Independent trial-factor check of the divisor envelope on a small range.
    running=1
    for n in range(1,10001):
        rem=n
        value=1
        p=2
        while p*p<=rem:
            exponent=0
            while rem%p==0:
                rem//=p
                exponent+=1
            if exponent and p%11 in (1,3,4,5,9):
                value*=exponent+1
            p+=1
        if rem>1 and rem%11 in (1,3,4,5,9):
            value*=2
        running=max(running,value)
        assert running==frontier[bisect_right(keys,n)-1][1]
    rows=[]
    worst=Fraction(0)
    worst_k=None
    for k in range(25,60):
        lower=1<<k
        upper=(1<<(k+1))-1
        cube=cubert_floor(lower)
        e_bound=Fraction(11,4*cube)
        m_bound=frontier[bisect_right(keys,upper)-1][1]
        lm=Fraction(8,3)*m_bound
        ratio=lm/(Fraction(5,2)*e_bound)
        cut=isqrt(ratio.numerator//ratio.denominator)
        assert 1<=cut<limit
        count=bisect_right(inert,cut)
        tail_bound=Fraction(suffix[count],scale)+tail
        left=Fraction(9,100)+lm*tail_bound+e_bound*(3+Fraction(5,2)*count)
        assert left<Fraction(2,3)-Fraction(1,1024)
        rows.append([k,m_bound,cube,cut,count,tail_bound.numerator,tail_bound.denominator,
                     left.numerator,left.denominator])
        if left>worst:
            worst=left
            worst_k=k
    return dict(prime_limit=limit,fixed_point_scale=scale,
                prime_flags_sha256=hashlib.sha256(flags).hexdigest(),
                inert_prime_count=len(inert),frontier_upper=1<<61,
                frontier_states=states,frontier=frontier,
                frontier_regression_bound=10000,
                bin_columns=['k','M(2^(k+1)-1)','floor_cuberoot(2^k)','P','pi_I(P)',
                             'tail_numerator','tail_denominator','criterion_numerator','criterion_denominator'],
                bins=rows,worst_bin=worst_k,worst_numerator=worst.numerator,
                worst_denominator=worst.denominator,uniform_margin='1/1024')


def verify_coverage_mask(covered,expected):
    if covered!=expected:
        raise AssertionError('Coverage differs from the exact target mask.')


def core_membership_flags(bound):
    """Exact C membership through bound, from complete forbidden-prime removal."""
    flags=prime_flags(bound)
    core=bytearray(b"\x01")*(bound+1)
    core[0]=0
    for p in range(2,bound+1):
        if flags[p] and p%11 not in (1,3,4,5,9):
            core[p:bound+1:p]=b"\x00"*(bound//p)
    return core


def finite_core_coverage():
    bound=1<<25
    short_bound=5000
    core=core_membership_flags(bound)
    assert core[0]==0 and core[1]==1 and core[2]==0 and core[11]==0
    # Independent direct primitive norm enumeration, then exclude 11 factors.
    regression_bound=10000
    direct=bytearray(regression_bound+1)
    radius=isqrt(2*regression_bound)
    for x in range(-radius,radius+1):
        for y in range(-radius,radius+1):
            q=x*x+x*y+3*y*y
            if 1<=q<=regression_bound and q%11 and gcd(x,y)==1:
                direct[q]=1
    assert core[:regression_bound+1]==direct
    packed=bytearray((bound+8)//8)
    for n in range(1,bound+1,2):
        if core[n]:
            packed[n>>3]|=1<<(n&7)
    bits=int.from_bytes(packed,'little')
    shorts=[n for n in range(1,short_bound+1) if core[n]]
    covered=0
    for a in shorts:
        covered|=bits<<a
    covered&=(1<<(bound+1))-1
    target_bytes=bytearray(b"\x55")*len(packed)
    target_bytes[0]&=254
    target_bytes[-1]&=(1<<((bound%8)+1))-1
    for n in range(22,bound+1,22):
        target_bytes[n>>3]&=255^(1<<(n&7))
    target=int.from_bytes(target_bytes,'little')
    verify_coverage_mask(covered,target)
    target_count=target.bit_count()
    assert target_count==bound//2-bound//22==15252015
    assert len(shorts)==774
    negative_cases=[('missing_target_2',covered^(1<<2)),
                    ('forbidden_odd_3',covered|(1<<3)),
                    ('forbidden_multiple_11_22',covered|(1<<22))]
    rejected=[]
    for name,corrupt in negative_cases:
        try:
            verify_coverage_mask(corrupt,target)
        except AssertionError:
            rejected.append(name)
        else:
            raise AssertionError('Corrupted coverage accepted: '+name)
    covered_bytes=covered.to_bytes(len(packed),'little')
    return dict(bound=bound,short_bound=short_bound,short_count=len(shorts),
                target_count=target_count,core_value_count=bits.bit_count(),
                missing_targets=0,forbidden_targets_covered=0,
                exact_target_mask_equal=True,
                membership_regression_bound=regression_bound,
                membership_regression_method='Independent Q/gcd square-box enumeration, excluding 11 factors.',
                negative_cases_rejected=rejected,
                core_flags_sha256=hashlib.sha256(core).hexdigest(),
                coverage_sha256=hashlib.sha256(covered_bytes).hexdigest(),
                target_mask_sha256=hashlib.sha256(target_bytes).hexdigest())



def unit_target_scope_checks(theta):
    """Challenge M11.6 on its stated 11-unit domain, with exact squares."""
    checked=0
    for n in range(2,len(theta),2):
        if n%11==0:
            continue
        m=n
        t=0
        while m%2==0:
            m//=2
            t+=1
        S=sigma(n)-(sigma(n//4) if n%4==0 else 0)
        G=theta[n]-(theta[n//4] if n%4==0 else 0)
        tau_m=sum(m%d==0 for d in range(1,m+1))
        delta_square=tau_m**2*(t+3)**2*n//4
        # G >= (12S-8Delta)/5; square only when its deficit is positive.
        deficit=12*S-5*G
        assert deficit<=0 or deficit**2<=64*delta_square, n
        checked+=1
    # This input is deliberately outside the stripped-target domain.
    n=22
    G=theta[n]
    S=sigma(n)
    delta_square=352  # Delta(22)=4*sqrt(22).
    deficit=12*S-5*G
    assert G==4 and S==36 and deficit==412
    assert deficit**2>64*delta_square
    return dict(domain='even n with11 not dividingn',
                checked_through=len(theta)-1,unit_targets_checked=checked,
                excluded_ramified_target=dict(n=n,G=G,S=S,Delta='4*sqrt(22)',
                    claimed_unrestricted_lower_bound='432/5-(32/5)*sqrt(22)',
                    deficit_squared=deficit**2,
                    squared_permitted_deficit=64*delta_square,
                    unrestricted_reading_rejected=True))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus11_checks.json'))
    args=parser.parse_args()
    bound=256
    r=binary_counts(bound)
    # Independent enclosing square-box enumeration: Q(x,y)>=(x*x+y*y)/2.
    brute=[0]*(bound+1)
    radius=isqrt(2*bound)
    for x in range(-radius,radius+1):
        for y in range(-radius,radius+1):
            q=x*x+x*y+3*y*y
            if q<=bound:
                brute[q]+=1
    assert r==brute
    matrix=[[-1,-6],[2,1]]
    gram=[[2,1],[1,6]]
    matrix_square=[[sum(matrix[i][k]*matrix[k][j] for k in range(2))
                    for j in range(2)] for i in range(2)]
    transformed_gram=[[sum(matrix[k][i]*gram[k][l]*matrix[l][j]
                          for k in range(2) for l in range(2))
                       for j in range(2)] for i in range(2)]
    assert matrix_square==[[-11,0],[0,-11]]
    assert transformed_gram==[[11*z for z in row] for row in gram]
    for x in range(-16,17):
        for y in range(-16,17):
            u,v=-x-6*y,2*x+y
            assert u*u+u*v+3*v*v==11*(x*x+x*y+3*y*y)
            assert (-u-6*v,2*u+v)==(-11*x,-11*y)
            if (x*x+x*y+3*y*y)%11==0:
                assert (x+6*y)%11==0 and (-2*x-y)%11==0
    theta=[sum(r[j]*r[n-j] for j in range(n+1))
           for n in range(bound+1)]
    a=eta11(bound)
    assert theta[:3]==[1,4,4]
    assert a[:3]==[0,1,-2]
    assert a[11]==1
    for n in range(1,bound+1):
        e=sigma(n)-11*(sigma(n//11) if n%11==0 else 0)
        assert 5*theta[n]==12*e+8*a[n], n
    # Finite challenges to the U_11 eigenrelation and ramified norm identity.
    for n in range(1,bound//11+1):
        assert a[11*n]==a[n], n
        assert theta[11*n]==theta[n], n
    # Exponent 10 and prime 64 suffice for the unbounded Euler proof.
    assert 12**6<2*11**6
    assert all((a+1)<=2**a for a in range(1,101))
    square_residues={1,3,4,5,9}
    split=lambda p:p%11 in square_residues
    inert=lambda p:p!=11 and p%11 not in square_residues
    certs=[euler_certificate('split_-11',split,7),
           euler_certificate('inert_-11_including_2',inert,15),
           euler_certificate('all_primes',lambda p:True,139),
           euler_certificate('primes_except_11',lambda p:p!=11,103),
           euler_certificate('odd_primes_except_11',lambda p:p not in (2,11),29)]
    inert_odd=[p for p in primes_below(65) if p>2 and inert(p)]
    b_small=sum((Fraction(2,p*(p-1)) for p in inert_odd),Fraction(0))
    starts=(65,73,79,83,85)
    b_upper=b_small+Fraction(7,3)*sum(
        (Fraction(1,a*a)+Fraction(1,22*a) for a in starts),Fraction(0))
    assert b_upper<Fraction(9,100)
    assert 2**37>71**6
    coarse_bound=Fraction(9,100)+Fraction(5375,132*71)+Fraction(8381,6*2**24)
    assert coarse_bound<Fraction(2,3)
    assert 9**6>2*8**6 and 10**6<2*9**6
    peak=Fraction(9**6,2**6)
    assert all(Fraction((t+3)**6,2**t)<=peak for t in range(1,7))
    assert 47**12<2**67
    p_coefficient=Fraction(4,5)
    global_first=29*Fraction(5,2)*Fraction(5,22)*p_coefficient+Fraction(20,3)*7*Fraction(5,22)/p_coefficient
    global_second=29*(3+Fraction(5,2)*Fraction(3,2))+Fraction(20,3)*7*5/p_coefficient**2
    assert (global_first,global_second)==(Fraction(1745,66),Fraction(1681,3))
    threshold_bound=Fraction(9,100)+Fraction(1745,66*47)+Fraction(1681,3*2**22)
    assert threshold_bound<Fraction(2,3)
    assert Fraction(2,3)-threshold_bound>Fraction(3,1000)
    core_bound=Fraction(9,100)+Fraction(1075,66*32)+Fraction(1805,6*2**20)
    assert core_bound<Fraction(2,3)
    global_margin=Fraction(2,3)-threshold_bound
    assert global_margin==Fraction(758567697,54211379200)
    assert Fraction(18,5)*Fraction(3,1000)/196>Fraction(1,20000)
    assert 11**2<5**3
    assert 16>11
    assert Fraction(29,2**22)<1 and Fraction(10,2**20)<1
    assert [r for r in range(1,22,2) if r%11 and r%11 not in square_residues]==[7,13,17,19,21]
    residues=[7,13,17,19,21]
    offsets=[]
    for remainder in range(22):
        count=sum(r<=remainder for r in residues)
        offset=Fraction(count)-Fraction(5*remainder,22)
        offsets.append(offset)
        assert offset<=Fraction(3,2)
    assert sum(22-r for r in residues)==33
    sharp=sharp_core_bounds()
    interval_certificate=core_interval_certificate()
    finite_coverage=finite_core_coverage()
    result=dict(checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                unit_target_scope=unit_target_scope_checks(theta),
                sharp_core_bounds=sharp,
                core_interval_certificate=interval_certificate,
                finite_core_coverage=finite_coverage,
                theta_identity=dict(sturm_bound=2,theta_coefficients_through_sturm=theta[:3],
                                    eta_coefficients_through_sturm=a[:3],eta_a11=a[11],
                                    secondary_formula_check_bound=bound,
                                    theta_level=11,theta_gram_determinant=121,
                                    independent_box_radius=radius,
                                    independent_box_counts_agree=True),
                ramified_map=dict(matrix=matrix,gram=gram,
                                  matrix_square=matrix_square,
                                  transformed_gram=transformed_gram,
                                  secondary_grid_min=-16,secondary_grid_max=16,
                                  secondary_grid_points=33**2),
                euler_certificates=certs,
                inert_density_certificate=dict(primes_through_64=inert_odd,tail_starts=list(starts),
                                               upper_numerator=b_upper.numerator,
                                               upper_denominator=b_upper.denominator,
                                               simple_upper='9/100',
                                               odd_inert_residues_mod_22=residues,
                                               period_counts_checked=list(range(22)),
                                               maximum_count_offset_numerator=max(offsets).numerator,
                                               maximum_count_offset_denominator=max(offsets).denominator,
                                               used_count_offset='3/2'),
                dyadic_satake_certificate=dict(peak_integer=6,peak_value="9/2",
                                               increasing_before_peak=(9**6>2*8**6),
                                               decreasing_after_peak=(10**6<2*9**6),
                                               E_definition="(2/3)*Delta/S",
                                               E_bound="29*n^(-1/3)",
                                               large_P_constant="4/5",
                                               remainder_coefficients=["1745/66","1681/3"]),
                coarse_secondary_threshold=dict(reduced_n_onset=2**74,original_N_onset=2**78,
                                                 numerator=coarse_bound.numerator,denominator=coarse_bound.denominator,
                                                 scope="Secondary check using the coarser triangle bound."),
                conditional_threshold=dict(reduced_n_onset=2**67,
                                           original_N_onset=2**71,
                                           criterion='B+(20/3)*M*T+E*(3+5*pi/2)<2/3',
                                           numerator=threshold_bound.numerator,
                                           denominator=threshold_bound.denominator,
                                           margin_numerator=global_margin.numerator,margin_denominator=global_margin.denominator,
                                           core_reduced_n_onset=2**60,
                                           core_criterion_bound_numerator=core_bound.numerator,
                                           core_criterion_bound_denominator=core_bound.denominator,
                                           reduced_distinct_value_lower_constant='1/20000',
                                           original_distinct_value_lower_constant='1/100000',
                                           scope='Arithmetic verification conditional on the stated sieve criterion; the core cutoff is a least-counterexample reduction, not a global onset.'),
                scope='Exact finite Sturm/Euler, lattice, divisor-frontier, dyadic interval and C+C coverage certificates; the written spectral/local argument and descent connect them to the unrestricted theorem.')
    out=json.dumps(result,indent=2)+'\n'
    args.output.write_text(out)
    print(out,end='')


if __name__=='__main__':
    main()
