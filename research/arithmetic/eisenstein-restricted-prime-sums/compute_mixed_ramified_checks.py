#!/usr/bin/env python3
"""Finite exact verification of the mixed problem's ramified-prime identity.

For N=1 mod3 define B3(N) as the number of representations
Q(x,y)+Q(z,w)=N having one component norm divisible by9, including norm0.
The conjectural exact identity checked here is B3(N)=4sigma(N)+8a(N),
where sum a(N)q^N = eta(3z)^2 eta(9z)^2.

This script establishes the listed finite checks only. An unrestricted
identity needs a modular-form/theta-series proof.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=3000)
    parser.add_argument("--output", type=Path, default=Path("mixed_ramified_checks.json"))
    args = parser.parse_args()
    started = time.monotonic()
    limit = args.limit
    sigma = [0]*(limit+1)
    tau = [0]*(limit+1)
    r = [0]*(limit+1)
    r[0] = 1
    for d in range(1,limit+1):
        chi = 0 if d%3 == 0 else (1 if d%3 == 1 else -1)
        for n in range(d,limit+1,d):
            sigma[n] += d
            tau[n] += 1
            r[n] += 6*chi
    # q times these products gives eta(3z)^2 eta(9z)^2.
    product = [1]+[0]*(limit-1)
    for spacing in (3,9):
        for m in range(spacing,limit,spacing):
            for j in range(limit-1,m-1,-1):
                product[j] -= 2*product[j-m]
                if j >= 2*m:
                    product[j] += product[j-2*m]
    a = [0]+product
    errors = []
    rows = []
    exact_one_third_count = 0
    for n in range(1,limit+1,3):
        total = sum(r[t]*r[n-t] for t in range(n+1))
        bad = 2*sum(r[u]*r[n-9*u] for u in range(n//9+1))
        expected = 4*sigma[n]+8*a[n]
        if total != 12*sigma[n] or bad != expected:
            errors.append({"n":n,"total":total,"bad":bad,"expected_bad":expected,"a_n":a[n]})
        if a[n]*a[n] > tau[n]*tau[n]*n:
            errors.append({"n":n,"a_n":a[n],"error":"Ramanujan bound violation"})
        if bad*3 == total:
            exact_one_third_count += 1
        if n <= 120:
            rows.append({"n":n,"total":total,"bad3":bad,"sigma":sigma[n],"a_n":a[n]})
    # Direct coordinate enumeration independently checks the binary theta
    # coefficients and then all ramified counts in the small range.
    direct_limit = min(120,limit)
    direct_r = [0]*(direct_limit+1)
    bound = math.isqrt(2*direct_limit)+2
    for x in range(-bound,bound+1):
        for y in range(-bound,bound+1):
            q=x*x-x*y+y*y
            if q<=direct_limit:
                direct_r[q]+=1
    direct_errors=[]
    for n in range(1,direct_limit+1,3):
        bad=sum(direct_r[t]*direct_r[n-t] for t in range(n+1) if t%9==0 or (n-t)%9==0)
        if bad != 4*sigma[n]+8*a[n]:
            direct_errors.append({"n":n,"direct_bad":bad})
    receipt={
        "scope":"Finite exact checks only; no unrestricted modular identity established by this computation",
        "limit":limit,
        "ramified_bad_definition":"Q(x,y)+Q(z,w)=N, with9 dividing one component norm; zero included",
        "candidate_identity":"B3(N)=4sigma(N)+8a_N forN=1mod3",
        "modular_form":"sum a_N q^N = eta(3z)^2 eta(9z)^2 = q product(1-q^(3m))^2(1-q^(9m))^2",
        "candidate_normalized_bound":"B3(N)/R(N)<=1/3+(2/3)tau(N)sqrt(N)/sigma(N)",
        "checked_targets":len(range(1,limit+1,3)),
        "exact_one_third_count":exact_one_third_count,
        "identity_and_bound_errors":errors,
        "direct_coordinate_count_limit":direct_limit,
        "direct_coordinate_count_errors":direct_errors,
        "small_values":rows,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds":round(time.monotonic()-started,3),
    }
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!='small_values'},indent=2))
    if errors or direct_errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
