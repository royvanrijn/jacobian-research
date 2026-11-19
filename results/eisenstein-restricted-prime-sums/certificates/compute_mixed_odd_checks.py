#!/usr/bin/env python3
"""Exact elementary constants and finite odd-norm/ramified checks.

The divisor-constant inequalities include a global exponent-tail argument.
The theta and CM identities below are checked only in the stated range.
"""

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


def main():
    local_two = [(e, Fraction((e+3)**5, 3**5 * 2**e)) for e in range(1,5)]
    max_two_e, product = max(local_two, key=lambda row: row[1])
    assert max_two_e == 4
    assert 8**5 < 2*7**5  # The next-term ratio decreases for all later e.
    rows=[{"prime":2,"maximizing_exponent":max_two_e,"local_fifth_power":str(product)}]
    for p in (5,11,17,23,29):
        end=0
        while (end+2)**5 > p*(end+1)**5:
            end+=1
        values=[]
        for e in range(end+1):
            factor=Fraction((e+1)*(p-1)*p**e,p**(e+1)-1)
            values.append(factor**5 / p**e)
        e=max(range(len(values)),key=values.__getitem__)
        product*=values[e]
        rows.append({"prime":p,"maximizing_exponent":e,"local_fifth_power":str(values[e]),"last_exponent_required":end})
    bound=Fraction(9,4)
    assert product<bound**5
    union=Fraction(29,60)+Fraction(11,10)*Fraction(2,5)+Fraction(111,20)*Fraction(1,250)+Fraction(2,3*10000)
    assert union==Fraction(591,625)<1

    limit=3000
    r=[0]*(limit+1);r[0]=1
    sigma=[0]*(limit+1)
    for d in range(1,limit+1):
        chi=0 if d%3==0 else (1 if d%3==1 else -1)
        for n in range(d,limit+1,d):
            r[n]+=6*chi;sigma[n]+=d
    eta=[1]+[0]*(limit-1)
    for spacing in (3,9):
        for d in range(spacing,limit,spacing):
            for j in range(limit-1,d-1,-1):
                eta[j]-=2*eta[j-d]
                if j>=2*d:eta[j]+=eta[j-2*d]
    cm=[0]+eta
    theta=[0]*(limit+1)
    bad=[0]*(limit+1)
    for n in range(1,limit+1,3):
        theta[n]=sum(r[q]*r[n-q] for q in range(n+1))
        bad[n]=2*sum(r[u]*r[n-9*u] for u in range(n//9+1))

    errors=[];lambda_rows=[];odd_checks=0;reduced_checks=0
    for n in range(4,limit+1,6):
        odd_checks+=1
        quarter=n//4 if n%4==0 else None
        total_odd=theta[n]-(theta[quarter] if quarter else 0)
        bad_odd=bad[n]-(bad[quarter] if quarter else 0)
        direct_total=sum(r[q]*r[n-q] for q in range(1,n,2))
        direct_bad=sum(r[q]*r[n-q] for q in range(1,n,2) if q%9==0 or (n-q)%9==0)
        e=0;m=n
        while m%2==0:e+=1;m//=2
        if total_odd!=direct_total or bad_odd!=direct_bad or total_odd!=18*2**e*sigma[m]:
            errors.append({"n":n,"error":"odd subtraction identity","total":total_odd,"direct_total":direct_total,"bad":bad_odd,"direct_bad":direct_bad})
        if bad_odd != 4*(sigma[n]-(sigma[quarter] if quarter else 0))+8*(cm[n]-(cm[quarter] if quarter else 0)):
            errors.append({"n":n,"error":"CM odd identity"})
        t=n;p=2;factorization={}
        while p*p<=t:
            while t%p==0:
                factorization[p]=factorization.get(p,0)+1;t//=p
            p+=1
        if t>1:factorization[t]=factorization.get(t,0)+1
        if any(p%3!=2 for p in factorization):continue
        reduced_checks+=1
        if any(exponent%2 for exponent in factorization.values()):
            lam=Fraction(0)
        else:
            sign=(-1)**sum(exponent//2 for exponent in factorization.values())
            cm_expected=sign*math.isqrt(n)
            if cm[n]!=cm_expected:errors.append({"n":n,"error":"reduced CM coefficient"})
            lam=Fraction(cm_expected,2**e*sigma[m])
        ratio=Fraction(bad_odd,total_odd)
        if ratio!=Fraction(1,3)+Fraction(2,3)*lam:
            errors.append({"n":n,"error":"reduced lambda formula","ratio":str(ratio),"lambda":str(lam)})
        if n<=120:lambda_rows.append({"n":n,"R_odd":total_odd,"B3_odd":bad_odd,"ratio":str(ratio),"lambda":str(lam)})

    direct_limit=120
    dr=[0]*(direct_limit+1)
    coordinate_bound=math.isqrt(2*direct_limit)+2
    for x in range(-coordinate_bound,coordinate_bound+1):
        for y in range(-coordinate_bound,coordinate_bound+1):
            q=x*x-x*y+y*y
            if q<=direct_limit:dr[q]+=1
    coordinate_checks=0
    for n in range(4,direct_limit+1,6):
        coordinate_checks+=1
        dt=sum(dr[q]*dr[n-q] for q in range(1,n,2))
        db=sum(dr[q]*dr[n-q] for q in range(1,n,2) if q%9==0 or (n-q)%9==0)
        quarter=n//4 if n%4==0 else None
        if dt!=theta[n]-(theta[quarter] if quarter else 0) or db!=bad[n]-(bad[quarter] if quarter else 0):
            errors.append({"n":n,"error":"coordinate enumeration"})

    receipt={
        "elementary_constant_scope":"Global exponent maxima, justified by monotone next-term ratios and p>=32 factors<=1",
        "E_bound":"E(N)<=(9/4)N^(-3/10), N=2^e*m,e>=1,m odd with prime factors2mod3",
        "local_maxima":rows,
        "product_fifth_power":str(product),
        "claimed_bound_fifth_power":str(bound**5),
        "strict_constant_inequality_verified":True,
        "assumed_union_bound":"29/60+(11/10)X^(-1/5)+(111/20)X^(-6/5)+2/(3X^2)",
        "X_range":"X>=100",
        "rational_upper_bound":str(union),
        "threshold":100000000,
        "finite_identity_scope":"Identities in this section checked finitely only; independent theta/CM proof required for all N",
        "odd_count_identity":"R_odd(N)=R(N)-R(N/4)=18*2^e*sigma(m)",
        "ramified_odd_identity":"B3_odd(N)=B3(N)-B3(N/4)",
        "reduced_lambda_formula":"B3_odd/R_odd=1/3+(2/3)lambda, lambda=0 for nonsquare reducedN; for square reducedN, lambda=(-1)^Omega(sqrtN)*sqrtN/(2^e*sigma(m))",
        "limit":limit,"odd_count_checks":odd_checks,"reduced_lambda_checks":reduced_checks,
        "coordinate_checks_through120":coordinate_checks,"errors":errors,"small_lambda_values":lambda_rows,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output=Path(__file__).with_name("mixed_odd_checks_receipt.json")
    output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!='small_lambda_values'},indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
