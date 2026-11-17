#!/usr/bin/env python3
"""Exact finite mixed certificate: n=a+3b, a,b supported on primes 1 mod3.

For n=6k+4, write a=6i+1 and b=6j+1.  The witness index j is valid exactly
when i=k-3j is nonnegative and both i,j are in the exact membership sieve.
No unrestricted spectral theorem is checked by this program.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import numpy as np


def digest(array):
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def membership(last):
    bound=math.isqrt(6*last+1)
    primes=np.ones(bound+1,dtype=np.bool_)
    primes[:2]=False
    for p in range(2,math.isqrt(bound)+1):
        if primes[p]: primes[p*p::p]=False
    good=np.ones(last+1,dtype=np.bool_)
    for p in np.flatnonzero(primes):
        p=int(p)
        if p%6==5: good[(5*p-1)//6::p]=False
    return good


def trial_allowed(a):
    if a<1: return False
    d=2
    while d*d<=a:
        if a%d==0:
            if d%3!=1: return False
            while a%d==0: a//=d
        d+=1
    return a==1 or a%3==1


def generate(good):
    missing=np.arange(len(good),dtype=np.int64)
    witness=np.full(len(good),-1,dtype="<i4")
    histogram=[]
    for j in np.flatnonzero(good):
        j=int(j)
        eligible=missing[missing>=3*j]
        hits=eligible[good[eligible-3*j]]
        if len(hits):
            witness[hits]=j
            histogram.append({"b":6*j+1,"new_targets":int(len(hits))})
            missing=missing[witness[missing]<0]
        if not len(missing): break
        if 3*j>missing[-1]: break
    if len(missing):
        raise RuntimeError(f"Uncovered target {6*int(missing[0])+4}")
    return witness,histogram


def verify(good,witness,chunk=1000000):
    if witness.dtype.kind not in "iu" or witness.ndim!=1 or len(witness)!=len(good):
        raise ValueError("Malformed witness array")
    for start in range(0,len(witness),chunk):
        k=np.arange(start,min(start+chunk,len(witness)),dtype=np.int64)
        j=witness[start:start+len(k)].astype(np.int64)
        if np.any(j<0): raise ValueError("Uncovered mixed target")
        i=k-3*j
        if np.any(i<0) or np.any(j>=len(good)):
            raise ValueError("Nonpositive summand or out-of-range witness")
        if not np.all(good[i]&good[j]): raise ValueError("Forbidden summand")
        a=6*i+1;b=6*j+1
        if not np.all(a+3*b==6*k+4): raise ValueError("Incorrect mixed sum")
    for j in np.unique(witness):
        if not trial_allowed(6*int(j)+1):
            raise ValueError("Small b summand has a forbidden prime factor")


def controls():
    good=membership(1666)
    mismatch=[k for k in range(len(good)) if bool(good[k])!=trial_allowed(6*k+1)]
    if mismatch: raise RuntimeError("Trial-factor membership control failed")
    small=membership(20);witness,histogram=generate(small);verify(small,witness)
    result={}
    for name,at,value in [("missing_target",0,-1),("negative_complement",0,1),
                          ("forbidden_b",12,4)]:
        corrupt=witness.copy();corrupt[at]=value
        try: verify(small,corrupt)
        except ValueError as error: result[name]=str(error)
        else: raise RuntimeError("Corrupt witness accepted: "+name)
    corrupt_good=small.copy();corrupt_good[4]=True
    corrupt=witness.copy();corrupt[12]=4
    try: verify(corrupt_good,corrupt)
    except ValueError as error: result["corrupt_membership_and_witness"]=str(error)
    else: raise RuntimeError("Corrupted membership and witness accepted")
    return {"membership_indices_checked":len(good),"largest_membership_value":6*(len(good)-1)+1,
            "independent_trial_factorization_matches":True,"negative_controls":result}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--limit",type=int,default=99999999)
    parser.add_argument("--output",type=Path,default=Path("mixed_certificate"))
    args=parser.parse_args()
    if args.limit<4: parser.error("limit must be at least4")
    start=time.monotonic();last=(args.limit-4)//6
    good=membership(last);witness,histogram=generate(good);verify(good,witness)
    args.output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output/"witnesses.npz",mixed=witness)
    receipt={"scope":"Exact finite mixed statement only; n=6k+4<=requested_limit",
             "requested_limit":args.limit,"largest_target":6*last+4,
             "number_of_targets":last+1,"allowed_summands_count":int(good.sum()),
             "unrestricted_exceptions":[],"all_witnesses_rechecked":True,
             "largest_least_b":int(6*witness.max()+1),"least_b_histogram":histogram,
             "witness_encoding":"For n=6k+4, signed little-endian int32 j means a=6(k-3j)+1,b=6j+1 and n=a+3b",
             "sha256_membership_uint8":digest(good),
             "sha256_mixed_least_b_index_int32_le":digest(witness),
             "controls":controls(),"python":sys.version,"numpy":np.__version__,
             "platform":platform.platform(),
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "elapsed_seconds":round(time.monotonic()-start,3)}
    (args.output/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!="least_b_histogram"},indent=2))


if __name__=="__main__": main()
