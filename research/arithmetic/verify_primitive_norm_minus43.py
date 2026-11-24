#!/usr/bin/env python3
"""Complete finite coverage for the discriminant -43 primitive norm sumset.

The written proof imports the separately pinned Heegner analytic reduction.
Default mode covers every eligible even43-unit through2**32. --pilot covers
four named segments only and can never report a complete classification.
"""
if not __debug__:
    raise RuntimeError('Active Python assertions are required; do not use -O or -OO.')

import argparse
import hashlib
import json
import resource
import sys
import time
from math import gcd, isqrt
from pathlib import Path

D = 43
QR = frozenset(i*i%D for i in range(1,D))
BOUND = 1 << 32
SPAN = 1 << 25
SHORT_BOUND = 49999
CORE_EXCEPTIONS = [4,6,8,10,16,20,38,50,56,74]
PRIMITIVE_EXCEPTIONS = [4,6,8,10,16,20,38,50,172,258,344,430,688,860,1634,2150,2408,3182]
ASCII_BITS = bytes([48,49]+list(range(2,256)))
ANALYTIC_PINS = {
    'PRIMITIVE_NORM_HEEGNER_REFINED_BOUNDS.md':'3a63897a0b8d5438e7912a66bb1f9c22c163143007f340d5cfc5a8e28e4f2384',
    'verify_primitive_norm_heegner_refined_bounds.py':'70c2b4c0a2964464b0590ccfa0c939341677c28d9bdb780d0f3b386f44790668',
    'primitive_norm_heegner_refined_bounds_checks.json':'3cfd95315edf8f74b95a443f55b92a3009eb197d0d84ea3074e926eff39b86d2',
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def primes(limit):
    a = bytearray(b'\x01')*(limit+1)
    a[:2] = b'\0\0'
    for p in range(2,isqrt(limit)+1):
        if a[p]:
            a[p*p::p] = b'\0'*((limit-p*p)//p+1)
    return [p for p,flag in enumerate(a) if flag]


def core_flags(lo,hi,odd_inert_primes):
    assert 1 <= lo <= hi and lo%2 == hi%2 == 1
    length = (hi-lo)//2+1
    cycle = bytes(int((lo+2*i)%D in QR) for i in range(D))
    a = bytearray(cycle*(length//D)+cycle[:length%D])
    for p in odd_inert_primes:
        if p*p > hi:
            break
        first = (-lo*((p+1)//2))%p
        if first < length:
            a[first::p] = b'\0'*((length-1-first)//p+1)
    return a


def trial_core(n,all_primes):
    assert n >= 1
    rem = n
    for p in all_primes:
        if p*p > rem:
            break
        if rem%p == 0:
            if p%D not in QR:
                return False
            while rem%p == 0:
                rem //= p
    return rem == 1 or rem%D in QR


def pack_flags(a):
    return int(a[::-1].translate(ASCII_BITS),2) if a else 0


def required_mask(lo,hi):
    assert 2 <= lo <= hi and lo%2 == hi%2 == 0
    length = (hi-lo)//2+1
    cycle = bytes(49 if (lo+2*i)%D else 48 for i in range(D))
    mask = int((cycle*(length//D)+cycle[:length%D])[::-1],2)
    for e in CORE_EXCEPTIONS:
        if lo <= e <= hi:
            mask &= ~(1 << ((e-lo)//2))
    return mask


def demand_cover(covered,required):
    assert required & ~covered == 0,'A required finite target has no certified positive pair'


def equal(a,b):
    assert a == b


def expect_rejected(name,action):
    try:
        action()
    except (AssertionError,ValueError):
        return name
    raise AssertionError('Negative control accepted: '+name)


def analytic_import():
    base = Path(__file__).resolve().parent
    for name,pin in ANALYTIC_PINS.items():
        assert sha((base/name).read_bytes()) == pin,('Changed analytic premise',name)
    receipt = json.loads((base/'primitive_norm_heegner_refined_bounds_checks.json').read_text())
    assert receipt['status'] == 'PASS'
    record = next(a for a in receipt['results'] if a['D'] == D)
    assert record['analytic_onset'] == BOUND
    assert record['large_finite_coverage_performed'] is False
    assert record['small_exact_controls']['core_exceptions_through_bound'] == CORE_EXCEPTIONS
    assert record['small_exact_controls']['candidate_primitive_exceptions'] == PRIMITIVE_EXCEPTIONS
    return {
        'statement':'HeegnerCoreAndQuotientOnsets H1 and HeegnerFiniteClassificationReduction H2 at D43',
        'pins':ANALYTIC_PINS,
        'scope':'Pinned mathematical/analytic premise. This scanner verifies source identity, not a new proof of the spectral theorem or replay of its analytic checker.',
    }


def controls(all_primes,inert):
    bound = 1 << 20
    direct = bytearray(bound//2)
    for y in range(-isqrt(4*bound//D),isqrt(4*bound//D)+1):
        s = isqrt(4*bound-D*y*y)
        for x in range(-((s+y)//2),(s-y)//2+1):
            q = x*x+x*y+11*y*y
            if q and gcd(x,y) == 1:
                assert q%2 == 1
                direct[(q-1)//2] = 1
    c = core_flags(1,bound-1,inert)
    expected = bytearray(c)
    ramified_start = (D-1)//2
    slots = len(expected[ramified_start::D])
    expected[ramified_start::D] = c[:slots]
    assert direct == expected
    direct_core = bytearray(direct)
    direct_core[ramified_start::D] = b'\0'*slots
    assert direct_core == c
    for i in range(5000):
        assert bool(c[i]) == trial_core(2*i+1,all_primes)
    short = [2*i+1 for i,flag in enumerate(c[:(SHORT_BOUND+1)//2]) if flag]
    assert short[0] == 1 and short[-1] <= SHORT_BOUND
    small = 4096
    cbits,pbits = pack_flags(c[:small//2]),pack_flags(direct[:small//2])
    csums = psums = 0
    for i in range(small//2):
        if c[i]:
            csums |= cbits << i
        if direct[i]:
            psums |= pbits << i
    for n in range(2,small+1,2):
        assert bool(csums & (1 << (n//2-1))) == (n%D != 0 and n not in CORE_EXCEPTIONS)
        assert bool(psums & (1 << (n//2-1))) == (n%(D*D) != 0 and n not in PRIMITIVE_EXCEPTIONS)
    assert required_mask(2,small) == sum(1 << (n//2-1) for n in range(2,small+1,2)
                                       if n%D and n not in CORE_EXCEPTIONS)
    assert sorted([e for e in CORE_EXCEPTIONS if e not in [56,74]]+[D*e for e in CORE_EXCEPTIONS]) == PRIMITIVE_EXCEPTIONS
    assert 13 in short and 31 in short and 56 == 13+D and 74 == 31+D
    # Probe exact membership and both bit-shift directions at low/high windows.
    membership_cases = shift_cases = 0
    for lo,hi in [(1,2047),(BOUND-SPAN+1,BOUND-SPAN+2047),(BOUND-2047,BOUND-1)]:
        flags = core_flags(lo,hi,inert)
        for i in list(range(16))+list(range(len(flags)-16,len(flags))):
            assert bool(flags[i]) == trial_core(lo+2*i,all_primes)
            membership_cases += 1
        bits = pack_flags(flags)
        target_lo = lo+1 if lo == 1 else lo+257
        for a in short[:24]:
            shift = (target_lo-a-lo)//2
            shifted = bits >> shift if shift >= 0 else bits << (-shift)
            for i in range(12):
                n = target_lo+2*i-a
                expected = lo <= n <= hi and trial_core(n,all_primes)
                assert bool(shifted & (1 << i)) == expected
                shift_cases += 1
    bad_square = bytearray(c)
    bad_square[4] = 1
    bad_three = bytearray(c)
    bad_three[1] = 1
    negatives = [
        expect_rejected('character-only mask admits inert square9',lambda: equal(bad_square,direct_core)),
        expect_rejected('false summand3 would create excluded target4',lambda: equal(bad_three,direct_core)),
        expect_rejected('omitted required target2',lambda: demand_cover(csums & ~1,required_mask(2,small))),
        expect_rejected('excluded target4 wrongly required',lambda: demand_cover(csums,required_mask(2,small)|2)),
        expect_rejected('nonbinary membership byte',lambda: pack_flags(bytearray([0,2,1]))),
        expect_rejected('primitive exceptions omit3182',lambda: equal(PRIMITIVE_EXCEPTIONS,PRIMITIVE_EXCEPTIONS[:-1])),
    ]
    return {
        'coordinate_bound':bound,'all_odd_flags_compared':bound//2,
        'primitive_coordinate_flags_sha256':sha(direct),'core_flags_sha256':sha(c),
        'trial_factor_initial_flags':5000,'small_complete_sumsets_through':small,
        'additional_low_high_membership_cases':membership_cases,'bit_shift_cases':shift_cases,
        'core_exceptions':CORE_EXCEPTIONS,'primitive_exceptions':PRIMITIVE_EXCEPTIONS,
        'negative_controls':negatives,
    },short


def segment(index,inert,short):
    lo = 2+index*SPAN
    hi = min(lo+SPAN-2,BOUND)
    assert lo <= hi <= BOUND and lo%2 == hi%2 == 0
    b0 = max(1,lo-short[-1])
    flags = core_flags(b0,hi-1,inert)
    membership_hash = sha(flags)
    bits = pack_flags(flags)
    del flags
    required = required_mask(lo,hi)
    covered = 0
    for used,a in enumerate(short,1):
        shift = (lo-a-b0)//2
        assert 2*shift == lo-a-b0
        covered |= bits >> shift if shift >= 0 else bits << (-shift)
        if used%16 == 0 or used == len(short):
            if required & ~covered == 0:
                break
    demand_cover(covered,required)
    count = (hi-lo)//2+1
    covered &= (1 << count)-1
    size = (count+7)//8
    return {
        'index':index,'lower':lo,'upper':hi,
        'summand_flags_lower':b0,'summand_flags_upper':hi-1,
        'required_targets':required.bit_count(),'short_summands_used':used,
        'membership_sha256':membership_hash,
        'target_mask_sha256':sha(required.to_bytes(size,'little')),
        'coverage_sha256':sha(covered.to_bytes(size,'little')),'missing_targets':0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--pilot',action='store_true',help='Only segments0,1,126,127; never full coverage.')
    args = parser.parse_args()
    start = time.perf_counter()
    imported = analytic_import()
    ps = primes(isqrt(BOUND))
    inert = [p for p in ps if p > 2 and p != D and p%D not in QR]
    control_result,short = controls(ps,inert)
    count = BOUND//SPAN
    indices = [0,1,count-2,count-1] if args.pilot else list(range(count))
    segments = []
    scan_start = time.perf_counter()
    for i in indices:
        row = segment(i,inert,short)
        segments.append(row)
        if args.pilot or (i+1)%16 == 0:
            print(f'segment{i+1}/{count} through{row["upper"]}; used{row["short_summands_used"]}; scan_elapsed={time.perf_counter()-scan_start:.2f}s',file=sys.stderr,flush=True)
    if not args.pilot:
        assert len(segments) == count == 128
        previous = 0
        for i,row in enumerate(segments):
            assert row['index'] == i and row['lower'] == previous+2
            previous = row['upper']
        assert previous == BOUND
        assert sum(row['required_targets'] for row in segments) == BOUND//2-BOUND//(2*D)-len(CORE_EXCEPTIONS)
    elapsed = time.perf_counter()-scan_start
    result = {
        'result':'PILOT_PASS' if args.pilot else 'PASS',
        'checker_sha256':sha(Path(__file__).read_bytes()),
        'arithmetic':'Exact integers; standard library only; assertions required.',
        'analytic_import':imported,'local_and_negative_controls':control_result,
        'finite_coverage':{
            'complete':not args.pilot,'lower':2,'upper':BOUND,
            'target_span_per_segment':SPAN,'segment_count':len(segments),
            'eligible_targets_verified':sum(row['required_targets'] for row in segments),
            'excluded_core_targets':CORE_EXCEPTIONS,'missing_targets':0,
            'short_summand_bound':SHORT_BOUND,'short_summand_count':len(short),'short_summands':short,
            'maximum_summands_used':max(row['short_summands_used'] for row in segments),
            'membership_method':'For odd b<=H: character(b)=+1 and no inert prime<=sqrt(H) divides b.',
            'segments':segments,
        },
        'scope':'Four named pilot segments only; no classification claim.' if args.pilot else
                'Every eligible finite target through2^32 is covered. Global classification additionally uses the pinned written analytic and descent theorems. No formal verification or external review.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'result':result['result'],'segments':len(segments),'scan_seconds':round(elapsed,3),
                      'total_seconds':round(time.perf_counter()-start,3),
                      'estimated_full_scan_seconds':round(elapsed*count/len(segments),3),
                      'peak_rss_mib':round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,2),
                      'maximum_short_summands_used':result['finite_coverage']['maximum_summands_used'],
                      'short_summands_available':len(short),'output':str(args.output)},sort_keys=True))


if __name__ == '__main__':
    main()
