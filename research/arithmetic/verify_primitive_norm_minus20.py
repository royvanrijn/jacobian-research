#!/usr/bin/env python3
"""Exact bounded x^2+5y^2 primitive norm sumset scan; no extrapolation.
Independent gcd-coordinate values and splitting/genus values are compared.
"""
import argparse
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json


def coordinate_flags(bound):
    flags=bytearray(bound+1)
    for y in range(isqrt(bound//5)+1):
        for x in range(isqrt(bound-5*y*y)+1):
            if gcd(x,y)==1:
                flags[x*x+5*y*y]=1
    flags[0]=0
    return flags


def prime_flags(bound):
    flags=bytearray(b"\x01")*(bound+1)
    flags[:2]=b"\x00\x00"
    for p in range(2,isqrt(bound)+1):
        if flags[p]:flags[p*p:bound+1:p]=b"\x00"*((bound-p*p)//p+1)
    return flags


def genus_flags(bound):
    primes=prime_flags(bound)
    flags=bytearray(b"\x01")*(bound+1)
    for p in range(3,bound+1,2):
        if primes[p] and p%20 in (11,13,17,19):
            flags[p:bound+1:p]=b"\x00"*(bound//p)
    # These residues encode the complete 2-adic primitive condition.
    for residue in (0,2,3,4,7):
        flags[residue:bound+1:8]=b"\x00"*len(range(residue,bound+1,8))
    flags[0:bound+1:25]=b"\x00"*(bound//25+1)
    return flags


def class_character_member(n):
    if n<1:return False
    parity=0
    exponent2=0
    while n%2==0:n//=2;exponent2+=1
    if exponent2>1:return False
    parity=exponent2
    exponent5=0
    while n%5==0:n//=5;exponent5+=1
    if exponent5>1:return False
    p=3
    while p*p<=n:
        exponent=0
        while n%p==0:n//=p;exponent+=1
        if exponent:
            if p%20 in (11,13,17,19):return False
            if p%20 in (3,7):parity+=exponent
        p+=2
    if n>1:
        if n%20 in (11,13,17,19):return False
        if n%20 in (3,7):parity+=1
    return parity%2==0


def bit_positions(integer):
    result=[]
    while integer:
        low=integer&-integer
        result.append(low.bit_length()-1)
        integer^=low
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound',type=int,default=1<<24)
    parser.add_argument('--short-bound',type=int,default=5000)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus20_checks.json'))
    args=parser.parse_args()
    B=args.bound
    assert B>=10000 and 1<=args.short_bound<=B
    flags=coordinate_flags(B)
    genus=genus_flags(B)
    assert flags==genus
    for n in range(1,10001):assert bool(flags[n])==class_character_member(n),n
    assert [n for n in range(1,31) if flags[n]]==[1,5,6,9,14,21,29,30]
    packed=bytearray((B+8)//8)
    for n in range(1,B+1):
        if flags[n]:packed[n>>3]|=1<<(n&7)
    bits=int.from_bytes(packed,'little')
    shorts=[n for n in range(1,args.short_bound+1) if flags[n]]
    covered=0
    for a in shorts:covered|=bits<<a
    covered&=(1<<(B+1))-1
    local_bytes=bytearray(b"\xdc")*len(packed)
    local_bytes[-1]&=(1<<((B%8)+1))-1
    local_mask=int.from_bytes(local_bytes,'little')
    assert not (covered&~local_mask)
    residual=bit_positions(local_mask&~covered)
    exceptions=[]
    for n in residual:
        if any(flags[a] and flags[n-a] for a in range(1,n//2+1)):
            covered|=1<<n
        else:exceptions.append(n)
    exact_mask=local_mask
    for n in exceptions:exact_mask&=~(1<<n)
    assert covered==exact_mask
    assert exceptions==[3,4]
    # Reject both a missing true value and false insertion of an exception.
    assert (covered&~(1<<2))!=exact_mask
    assert (covered|(1<<3))!=exact_mask
    assert (covered|(1<<5))!=exact_mask
    bad_flags=bytearray(flags);bad_flags[3]=1
    assert bad_flags!=genus
    receipt=dict(checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 bound=B,short_bound=args.short_bound,short_count=len(shorts),
                 primitive_value_count=bits.bit_count(),
                 coordinate_genus_flag_equality_bound=B,
                 independent_trial_factor_class_character_regression_bound=10000,
                 short_sweep_uncovered_local_targets=residual,
                 exhaustive_residual_scan=True,
                 locally_admissible_global_exceptions=exceptions,
                 local_condition='N mod 8 in {2,3,4,6,7}',
                 locally_admissible_target_count=local_mask.bit_count(),
                 represented_target_count=covered.bit_count(),
                 forbidden_targets_covered=0,
                 exact_target_mask_equal=True,
                 negative_checks=['missing represented target 2','injected exception 3',
                                  'injected locally forbidden target 5','injected nonprincipal norm 3'],
                 primitive_flags_sha256=hashlib.sha256(flags).hexdigest(),
                 coverage_sha256=hashlib.sha256(covered.to_bytes(len(packed),'little')).hexdigest(),
                 scope='Exact only for 2<=N<=bound; does not assert an unrestricted sumset identity. The class-character norm criterion is the classical PL4 specialization.')
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
