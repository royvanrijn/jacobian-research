#!/usr/bin/env python3
"""Exact bounded primitive Q_-23 self-sum audit; no unrestricted exception claim."""


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

QR23={a*a%23 for a in range(1,23)}


def coordinate_flags(bound,a,b,c):
    assert b*b-4*a*c == -23
    flags=bytearray(bound+1)
    for y in range(isqrt(4*a*bound//23)+1):
        h=isqrt(4*a*bound-23*y*y)
        lo=-((h+b*y)//(2*a))
        hi=(h-b*y)//(2*a)
        for x in range(lo,hi+1):
            if gcd(x,y)==1:
                n=a*x*x+b*x*y+c*y*y
                if n:
                    assert n<=bound
                    flags[n]=1
    return flags


def prime_flags(bound):
    flags=bytearray(b'\1')*(bound+1)
    flags[:2]=b'\0\0'
    for p in range(2,isqrt(bound)+1):
        if flags[p]:
            flags[p*p:bound+1:p]=b'\0'*((bound-p*p)//p+1)
    return flags


def ideal_flags(bound,nonprincipal_form,principal_form):
    primes=prime_flags(bound)
    allowed=bytearray(b'\1')*(bound+1)
    active=bytearray(bound+1)
    counts=dict(principal_split=0,nonprincipal_split=0,inert=0,ramified=0)
    for p in range(2,bound+1):
        if not primes[p]:
            continue
        if p==23:
            counts['ramified']+=1
            assert principal_form[p] and not nonprincipal_form[p]
        elif p%23 not in QR23:
            counts['inert']+=1
            assert not principal_form[p] and not nonprincipal_form[p]
            allowed[p:bound+1:p]=b'\0'*(bound//p)
        else:
            assert bool(principal_form[p]) != bool(nonprincipal_form[p])
            if nonprincipal_form[p]:
                counts['nonprincipal_split']+=1
                # The cumulative updates are1 at valuations1,2mod3,
                # and0 at valuations0mod3. No factorization loop is used.
                power,exponent=p,1
                while power<=bound:
                    delta=1 if exponent%3==1 else (-1 if exponent%3==0 else 0)
                    if delta:
                        for n in range(power,bound+1,power):
                            active[n]+=delta
                    power*=p
                    exponent+=1
            else:
                counts['principal_split']+=1
    allowed[0:bound+1:529]=b'\0'*(bound//529+1)
    local_candidates=allowed.count(1)
    class_failures=sum(bool(ok) and r==1 for ok,r in zip(allowed,active))
    result=bytearray(bool(ok) and r!=1 for ok,r in zip(allowed,active))
    return result,counts,local_candidates,class_failures


def trial_member(n,q1):
    if n<1:
        return False
    r=0
    p=2
    while p*p<=n:
        exponent=0
        while n%p==0:
            exponent+=1
            n//=p
        if exponent:
            if p==23:
                if exponent>=2:
                    return False
            elif p%23 not in QR23:
                return False
            elif q1[p] and exponent%3:
                r+=1
        p=3 if p==2 else p+2
    if n>1:
        if n!=23 and n%23 not in QR23:
            return False
        if n!=23 and q1[n]:
            r+=1
    return r!=1


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
    parser.add_argument('--short-bound',type=int,default=5000)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus23_sumset_checks.json'))
    args=parser.parse_args()
    bound=args.bound
    if not bound>=10000 or not 1<=args.short_bound<=bound:
        parser.error('requires --bound >= 10000 and 1 <= --short-bound <= --bound')
    flags=coordinate_flags(bound,1,1,6)
    q1=coordinate_flags(bound,2,1,3)
    ideal,prime_counts,local_norm_count,class_failures=ideal_flags(bound,q1,flags)
    require_equal(flags,ideal)
    for n in range(1,10001):
        assert bool(flags[n])==trial_member(n,q1),n
    # Independent square-box enum challenges the ellipse signs/endpoints.
    small=1000
    square=bytearray(small+1)
    ycap=isqrt(4*small//23)+1
    xcap=isqrt(small)+ycap+1
    for x in range(-xcap,xcap+1):
        for y in range(-ycap,ycap+1):
            n=x*x+x*y+6*y*y
            if 1<=n<=small and gcd(x,y)==1:
                square[n]=1
    require_equal(square,flags[:small+1])
    assert not flags[2] and not flags[4] and flags[8] and not flags[351]
    packed=bytearray((bound+8)//8)
    for n in range(1,bound+1):
        if flags[n]:
            packed[n>>3]|=1<<(n&7)
    bits=int.from_bytes(packed,'little')
    shorts=[n for n in range(1,args.short_bound+1) if flags[n]]
    covered=0
    for a in shorts:
        covered|=bits<<a
    cap=(1<<(bound+1))-1
    covered&=cap
    local_bytes=bytearray(b'\xff')*len(packed)
    local_bytes[0]&=~3
    for n in range(529,bound+1,529):
        local_bytes[n>>3]&=~(1<<(n&7))
    local_bytes[-1]&=(1<<((bound%8)+1))-1
    local_mask=int.from_bytes(local_bytes,'little')
    assert local_mask.bit_count()==bound-1-bound//529
    assert not covered&~local_mask
    residual=bit_positions(local_mask&~covered)
    exceptions=[]
    residual_witnesses=[]
    for n in residual:
        witness=next((a for a in range(1,n//2+1) if flags[a] and flags[n-a]),None)
        if witness is None:
            exceptions.append(n)
        else:
            assert 1<=witness<=n//2 and ideal[witness] and ideal[n-witness]
            residual_witnesses.append([n,witness,n-witness])
            covered|=1<<n
    assert len(residual)==len(exceptions)+len(residual_witnesses)
    expected=local_mask
    for n in exceptions:
        expected&=~(1<<n)
    require_equal(covered,expected)
    assert 3 in exceptions and 2 not in exceptions
    expect_rejection(lambda:require_equal(covered&~(1<<2),expected))
    expect_rejection(lambda:require_equal(covered|(1<<3),expected))
    expect_rejection(lambda:require_equal(covered|(1<<529),expected))
    bad=bytearray(flags)
    bad[2]=1
    expect_rejection(lambda:require_equal(bad,ideal))
    bad=bytearray(flags)
    bad[351]=1
    expect_rejection(lambda:require_equal(bad,ideal))
    receipt=dict(
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        bound=bound,short_bound=args.short_bound,short_count=len(shorts),
        primitive_value_count=bits.bit_count(),
        local_individual_norm_candidate_count=local_norm_count,
        individual_class_obstruction_count=class_failures,
        complete_coordinate_ideal_membership_equality_bound=bound,
        independent_trial_factor_regression_bound=10000,
        independent_square_box_regression_bound=small,
        prime_class_counts=prime_counts,
        local_additive_condition='2<=N<=bound and529 does not divide N',
        locally_excluded_529_multiple_count=bound//529,
        locally_admissible_target_count=local_mask.bit_count(),
        short_sweep_uncovered_local_targets=residual,
        every_short_sweep_residual_resolved=True,
        exhaustive_nonrepresentation_residual_checks=True,
        residual_witness_count=len(residual_witnesses),
        residual_witnesses=residual_witnesses,
        locally_admissible_global_exceptions=exceptions,
        exception_count=len(exceptions),largest_exception=max(exceptions,default=None),
        represented_target_count=covered.bit_count(),
        forbidden_targets_covered=0,exact_target_mask_equal=True,
        negative_controls=['deleted represented target2','injected exception3',
                           'injected locally forbidden target529',
                           'injected nonprincipal individual norm2',
                           'injected exponent-relaxation norm351'],
        primitive_flags_sha256=sha256(flags).hexdigest(),
        ideal_flags_sha256=sha256(ideal).hexdigest(),
        coverage_sha256=sha256(covered.to_bytes(len(packed),'little')).hexdigest(),
        scope='Exact bounded self-sum classification only for2<=N<=bound. Complete local exclusions are separated from the listed finite global exceptions. No unrestricted exception classification is asserted.')
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in (
        'checker_sha256','bound','short_bound','short_count','primitive_value_count',
        'locally_excluded_529_multiple_count','locally_admissible_target_count',
        'represented_target_count','residual_witness_count','exception_count',
        'largest_exception','locally_admissible_global_exceptions','scope')},indent=2))


if __name__=='__main__':
    main()
