#!/usr/bin/env python3
"""Exact endpoint certificates for the elementary prime wheel estimates.

Checking the listed endpoints proves the discrepancy bound for all real t:
the discrepancy is linear and decreasing between jumps, and30-periodic.
The count and reciprocal-square tail inequalities then follow algebraically.
"""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def main():
    residues=(11,17,23,29)
    endpoints=[]
    for t in (0,11,17,23,29,30):
        before=sum(r<t for r in residues)
        after=sum(r<=t for r in residues)
        endpoints.append({"t":t,"left_limit":str(F(before)-F(2*t,15)),"value":str(F(after)-F(2*t,15))})
    values=[F(row[k]) for row in endpoints for k in ("left_limit","value")]
    minimum,maximum=min(values),max(values)
    assert minimum==F(-22,15) and maximum==F(2,15)
    assert maximum-minimum==F(8,5)
    main=F(1,10)+F(1,55)+F(1,136)+F(1,48)
    assert main==F(6569,44880)<F(3,20)
    mixed_leading=F(9,4)*F(4,15)+4*F(15,8)*F(1,15)
    mixed_next=F(9,4)*F(17,15)+4*F(15,8)*F(2,5)
    assert mixed_leading==F(11,10) and mixed_next==F(111,20)
    original_leading=F(5,2)*F(4,15)+6*F(15,8)*F(1,15)
    original_next=F(5,2)*F(17,15)+6*F(15,8)*F(2,5)
    assert original_leading==F(17,12) and original_next==F(22,3)
    receipt={
        "scope":"Exact finite endpoint and rational certificates for an elementary all-real-variable proof",
        "wheel_modulus":30,"wheel_residues":list(residues),"endpoints":endpoints,
        "minimum_discrepancy_or_left_limit":str(minimum),"maximum_discrepancy":str(maximum),
        "discrepancy_width":str(maximum-minimum),
        "global_argument":"A(t+30)=A(t)+4; A(t)-2t/15 is decreasing linear between its four jumps. Thus the endpoint extrema control all real t>=0.",
        "prime_coverage":"Every odd prime p=2mod3 other than5 lies in11,17,23,29mod30. Prime2 is handled separately or filtered; prime3 is not inert and is handled separately in the mixed problem.",
        "count_bound":"#{odd p=2mod3:p<=P}<=2P/15+17/15 forP>=5 (the extra1 counts prime5)",
        "tail_argument":"ForP>=5 all primesp>P lie in the wheel; A(t)-A(P)<=2(t-P)/15+8/5. Integrating2(A(t)-A(P))/t^3 over[P,infinity) gives2/(15P)+8/(5P^2).",
        "main_sum_argument":"All odd inert primes are6k-1,k>=1. Keepk=1,2,3; boundk>=4 by integral from3 toinfinity of2/(6x-2)^2, which equals1/48.",
        "main_sum_rational_upper_bound":str(main),"main_sum_claimed_upper_bound":"3/20",
        "application_cutoff":"P=2*N^(1/4)>=200 whenN>=100000000; all wheel estimates therefore apply without rounding or omitted primes",
        "mixed_small_and_large_prime_coefficients":[str(mixed_leading),str(mixed_next)],
        "original_small_and_large_prime_coefficients":[str(original_leading),str(original_next)],
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("mixed_wheel_checks_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
