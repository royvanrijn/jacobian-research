#!/usr/bin/env python3
"""Exact controls for UNRESTRICTED_SFC_MECHANISMS.tex."""
from math import comb, factorial
from pathlib import Path
import hashlib, json

def moments(nmax):
    return [sum(comb(n,k)*factorial(2*k) for k in range(n+1))
            for n in range(nmax+1)]

def main():
    ms=moments(32)
    for n in range(2,33):
        assert ms[n]-2*n*(2*n-1)*ms[n-1]+4*n*(n-1)*ms[n-2]==1
        assert ms[n]>=factorial(2*n)
    for n in range(3,33):
        assert (ms[n]-(4*n*n-2*n+1)*ms[n-1]
                +(n-1)*(8*n-6)*ms[n-2]
                -4*(n-1)*(n-2)*ms[n-3])==0
    assert ms[2]-(2*2*(2*2-1)+1)*ms[1]+4*2*(2-1)*ms[0]!=1
    print(json.dumps({
        "status":"PASS_exact_controls_only",
        "inhomogeneous_controls":31,
        "homogeneous_controls":30,
        "first_moments":ms[:5],
        "corrupt_coefficient_rejected":True,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope":"All-index recurrence and minimal-order-three proof are written; finite controls do not establish minimality."
    },indent=2))

if __name__=="__main__":
    main()
