#!/usr/bin/env python3
"""Exact bounded primitive x^2+14y^2 self-sums; no unrestricted exception claim."""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
from hashlib import sha256
import json
from math import gcd,isqrt
from pathlib import Path


def coordinate_flags(bound,a,b,c):
    assert b*b-4*a*c==-56
    flags=bytearray(bound+1)
    if b==0:
        for y in range(isqrt(bound//c)+1):
            for x in range(isqrt((bound-c*y*y)//a)+1):
                if gcd(x,y)==1:
                    flags[a*x*x+c*y*y]=1
    else:
        for y in range(isqrt(4*a*bound//56)+1):
            root=isqrt(4*a*bound-56*y*y)
            for x in range(-((root+b*y)//(2*a)),(root-b*y)//(2*a)+1):
                if gcd(x,y)==1:
                    flags[a*x*x+b*x*y+c*y*y]=1
    flags[0]=0
    return flags


def prime_flags(bound):
    flags=bytearray(b'\1')*(bound+1)
    flags[:2]=b'\0\0'
    for p in range(2,isqrt(bound)+1):
        if flags[p]:
            flags[p*p:bound+1:p]=b'\0'*((bound-p*p)//p+1)
    return flags


def chi8(n):
    return 0 if n%2==0 else (1 if n%8 in (1,7) else -1)


def chi56(p):
    assert p not in (2,7)
    return chi8(p)*(1 if p%7 in (1,2,4) else -1)


def ideal_flags(bound,q0,q2,qodd):
    primes=prime_flags(bound)
    allowed=bytearray(b'\1')*(bound+1)
    active=bytearray(bound+1)
    fixed=bytearray(bound+1)
    counts=dict(principal_split=0,order_two_split=0,odd_class_split=0,inert=0,ramified=0)
    for p in range(2,bound+1):
        if not primes[p]:
            continue
        if p in (2,7):
            counts['ramified']+=1
            assert q2[p] and not q0[p] and not qodd[p]
            allowed[p*p:bound+1:p*p]=b'\0'*(bound//(p*p))
            if p==2:
                fixed[2:bound+1:2]=b'\1'*(bound//2)
            else:
                for n in range(p,bound+1,p):
                    fixed[n]^=1
        elif chi56(p)==-1:
            counts['inert']+=1
            assert not q0[p] and not q2[p] and not qodd[p]
            allowed[p:bound+1:p]=b'\0'*(bound//p)
        else:
            assert sum(bool(q[p]) for q in (q0,q2,qodd))==1
            if qodd[p]:
                counts['odd_class_split']+=1
                assert chi8(p)==-1
                power,exponent=p,1
                while power<=bound:
                    delta=1 if exponent%2 else -1
                    for n in range(power,bound+1,power):
                        active[n]+=delta
                    # For an even exponent, fixed class2 occurs exactly
                    # at2mod4. Toggle the cumulative bit at2,3mod4.
                    if exponent%4 in (2,3):
                        for n in range(power,bound+1,power):
                            fixed[n]^=1
                    power*=p
                    exponent+=1
            elif q2[p]:
                counts['order_two_split']+=1
                assert chi8(p)==1
                power=p
                while power<=bound:
                    for n in range(power,bound+1,power):
                        fixed[n]^=1
                    power*=p
            else:
                counts['principal_split']+=1
                assert chi8(p)==1
    allowed[0]=0
    genus=bytearray(ok and n%16 in (1,2,7,9,14,15) for n,ok in enumerate(allowed))
    class_genus=bytearray(ok and r%2==0 for ok,r in zip(allowed,active))
    require_equal(genus,class_genus)
    bad_genus=bytearray(genus);bad_genus[3]=1
    expect_rejection(lambda:require_equal(bad_genus,class_genus))
    result=bytearray(ok and r%2==0 and (r>=2 or b==0)
                     for ok,r,b in zip(allowed,active,fixed))
    return result,counts,allowed.count(1),genus


def trial_member(n,q2,qodd):
    if n<1:
        return False
    r=b=0
    p=2
    while p*p<=n:
        exponent=0
        while n%p==0:
            exponent+=1
            n//=p
        if exponent:
            if p in (2,7):
                if exponent>=2:
                    return False
                b+=exponent
            elif chi56(p)==-1:
                return False
            elif q2[p]:
                b+=exponent
            elif qodd[p]:
                if exponent%2:
                    r+=1
                else:
                    b+=exponent//2
        p=3 if p==2 else p+2
    if n>1:
        if n in (2,7):
            b+=1
        elif chi56(n)==-1:
            return False
        elif q2[n]:
            b+=1
        elif qodd[n]:
            r+=1
    return r%2==0 and (r>=2 or b%2==0)


def bit_positions(value):
    result=[]
    while value:
        low=value&-value
        result.append(low.bit_length()-1)
        value^=low
    return result


def require_equal(actual,expected):
    if actual!=expected:
        raise AssertionError('exact mask/flag comparison failed')


def expect_rejection(action):
    try:
        action()
    except AssertionError:
        return
    raise AssertionError('corruption accepted')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound',type=int,default=1<<24)
    parser.add_argument('--short-bound',type=int,default=25000)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus56_sumset_checks.json'))
    args=parser.parse_args()
    bound=args.bound
    if not bound>=10000 or not 1<=args.short_bound<=bound:
        parser.error('requires --bound >= 10000 and 1 <= --short-bound <= --bound')
    flags=coordinate_flags(bound,1,0,14)
    q2=coordinate_flags(bound,2,0,7)
    qodd=coordinate_flags(bound,3,2,5)
    ideal,prime_counts,support_count,genus=ideal_flags(bound,flags,q2,qodd)
    local_count=genus.count(1)
    require_equal(flags,ideal)
    for n in range(1,10001):
        assert bool(flags[n])==trial_member(n,q2,qodd),n
    small=1000
    square=bytearray(small+1)
    for x in range(-isqrt(small)-1,isqrt(small)+2):
        for y in range(-isqrt(small//14)-1,isqrt(small//14)+2):
            n=x*x+14*y*y
            if 1<=n<=small and gcd(x,y)==1:
                square[n]=1
    require_equal(square,flags[:small+1])
    assert not flags[2] and not flags[7] and not flags[9] and flags[14] and flags[15] and flags[18] and not flags[450]
    packed=bytearray((bound+8)//8)
    for n in range(1,bound+1):
        if flags[n]:
            packed[n>>3]|=1<<(n&7)
    bits=int.from_bytes(packed,'little')
    shorts=[n for n in range(1,args.short_bound+1) if flags[n]]
    covered=0
    for a in shorts:
        covered|=bits<<a
    covered&=(1<<(bound+1))-1
    local_bytes=bytearray(b'\xff')*len(packed)
    local_bytes[0]&=~3
    for n in range(49,bound+1,49):
        local_bytes[n>>3]&=~(1<<(n&7))
    local_bytes[-1]&=(1<<((bound%8)+1))-1
    local_mask=int.from_bytes(local_bytes,'little')
    assert local_mask.bit_count()==bound-1-bound//49
    assert not covered&~local_mask
    residual=bit_positions(local_mask&~covered)
    exceptions=[]
    witnesses=[]
    for n in residual:
        a=next((a for a in range(1,n//2+1) if flags[a] and flags[n-a]),None)
        if a is None:
            exceptions.append(n)
        else:
            assert ideal[a] and ideal[n-a]
            witnesses.append([n,a,n-a])
            covered|=1<<n
    assert len(residual)==len(exceptions)+len(witnesses)
    expected=local_mask
    for n in exceptions:
        expected&=~(1<<n)
    require_equal(covered,expected)
    assert 3 in exceptions and 2 not in exceptions
    expect_rejection(lambda:require_equal(covered&~(1<<2),expected))
    expect_rejection(lambda:require_equal(covered|(1<<3),expected))
    expect_rejection(lambda:require_equal(covered|(1<<49),expected))
    bad=bytearray(flags);bad[9]=1
    expect_rejection(lambda:require_equal(bad,ideal))
    bad=bytearray(flags);bad[450]=1
    expect_rejection(lambda:require_equal(bad,ideal))
    receipt=dict(
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),bound=bound,
        short_bound=args.short_bound,short_count=len(shorts),primitive_value_count=bits.bit_count(),
        individual_support_candidate_count=support_count,individual_local_candidate_count=local_count,
        individual_class_obstruction_count=local_count-bits.bit_count(),
        complete_coordinate_class_membership_equality_bound=bound,
        complete_genus_parity_vs_dyadic_equality_bound=bound,
        independent_trial_factor_regression_bound=10000,independent_square_box_regression_bound=small,
        prime_class_counts=prime_counts,
        local_additive_condition='2<=N<=bound and49 does not divide N',
        locally_excluded_49_multiple_count=bound//49,locally_admissible_target_count=local_mask.bit_count(),
        short_sweep_uncovered_local_targets=residual,every_short_sweep_residual_resolved=True,
        exhaustive_nonrepresentation_residual_checks=True,residual_witness_count=len(witnesses),residual_witnesses=witnesses,
        locally_admissible_global_exceptions=exceptions,exception_count=len(exceptions),largest_exception=max(exceptions,default=None),
        represented_target_count=covered.bit_count(),forbidden_targets_covered=0,exact_target_mask_equal=True,
        negative_controls=['deleted represented target2','injected exception3','injected local exclusion49',
                           'injected wrong-genus individual norm3',
                           'injected hidden class obstruction individual norm9','injected restoration-guard norm450'],
        primitive_flags_sha256=sha256(flags).hexdigest(),ideal_flags_sha256=sha256(ideal).hexdigest(),
        individual_local_genus_flags_sha256=sha256(genus).hexdigest(),
        additive_local_mask_sha256=sha256(local_bytes).hexdigest(),
        coverage_sha256=sha256(covered.to_bytes(len(packed),'little')).hexdigest(),
        scope='Exact bounded self-sum classification only for2<=N<=bound. The49-multiple local exclusions are separated from finite global exceptions. No unrestricted exception classification or cutoff is asserted.')
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('checker_sha256','bound','short_bound','short_count','primitive_value_count',
        'locally_excluded_49_multiple_count','locally_admissible_target_count','represented_target_count',
        'residual_witness_count','exception_count','largest_exception','locally_admissible_global_exceptions','scope')},indent=2))


if __name__=='__main__':
    main()
