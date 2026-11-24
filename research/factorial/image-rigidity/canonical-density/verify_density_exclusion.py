#!/usr/bin/env python3
"""Exact arithmetic replay of the two canonical-density nonreachability proofs."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path

BASE = Path(__file__).resolve().parent
PROOF_NAME = 'canonical_density_exclusion.md'
PROOF_SHA256 = '3a4a48ca11658e7bb87985f4b06a09548b064e50e2e498efe37a17c63b8a409e'

def require(condition, message):
    if not condition:
        raise ArithmeticError(message)

def rising(a, j):
    value = 1
    for i in range(j):
        value *= a + i
    return value

def beta_moment(a, b, j):
    return Q(rising(a, j), rising(a+b, j))

def beta_leading(a, b):
    return Q((-1)**(b-1)*factorial(a+b-1), factorial(a-1)*factorial(b-1))

def central(a, b, j):
    mean = Q(a, a+b)
    return sum(Q(comb(j, i)) * (-mean)**(j-i) * beta_moment(a, b, i)
               for i in range(j+1))

def first_pivot(n, r, u):
    a,b=n*r+1,n*u+1
    c,d=(n+1)*r+1,(n+1)*u+1
    alpha=(n+1)*beta_moment(c,d,1)-n*beta_moment(a,b,1)
    return (comb(n+1,2)*beta_moment(c,d,2)
            -comb(n,2)*beta_moment(a,b,2)-alpha*n*beta_moment(a,b,1))

def current_candidates(k, D, p, q):
    result = []
    # Every later density is divisible by t^ell(1-t)^ell.
    for ell in range(1, min(p, q)+1):
        n = k+ell
        value = D-2*ell
        if value < 0 or value % (n+1):
            continue
        s = value//(n+1)
        for r in range(s+1):
            u = s-r
            if n*r+ell <= p and n*u+ell <= q:
                result.append(dict(n=n, ell=ell, r=r, u=u, s=s))
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=BASE/'density_exclusion_verification.json')
    args=parser.parse_args()
    require(sha256((BASE/PROOF_NAME).read_bytes()).hexdigest()==PROOF_SHA256, 'proof authentication')
    eps = Q(1, 10000)
    weight = Q(701, 1344)
    first = weight*beta_leading(111,143)
    second = (1-weight)*beta_leading(132,122)
    require(first>0 and second<0 and first+second!=0, 'D1 exact previous degree')
    ratio = -first/second
    require(ratio==Q(1434214893826,8212939136313), 'D1 top-term ratio')
    d1mean = (1-eps)*(weight*beta_moment(111,143,1)
             +(1-weight)*beta_moment(132,122,1))+eps*beta_moment(111,122,1)
    d1difference = d1mean-beta_moment(111,122,1)
    require(d1difference==Q(9999,59648000), 'D1 first moment mismatch')
    d1candidates=current_candidates(10,254,111,122)
    require(d1candidates==[dict(n=11,ell=1,r=10,u=11,s=21)],
            'D1 complete current-density candidates')
    # Its sole candidate has previous degree n*s, not the observed 252.
    require(d1candidates[0]['n']*d1candidates[0]['s']!=252,
            'D1 previous degree exclusion')
    d1pivot=first_pivot(11,10,11)
    require(d1pivot==Q(135575957,22855940445) and d1pivot>0,
            'D1 forced first multiplier')
    require((-1)**(12*11)>0 and beta_leading(122,134)<0,
            'D1 current leading-sign exclusion')

    Z=sum(Q(110**j,factorial(j))*central(56,67,2*j) for j in range(6))
    shift=sum(Q(110**j,factorial(j))*central(56,67,2*j+1)
              for j in range(6))/Z
    require(Z==Q(4051587308496779144355663972757,
                 3046760870601420524860394238360), 'D2 normalization')
    require(shift==Q(19243030530352707943816806380284,
                     22425535752529672564008600089209995),
            'D2 first moment mismatch')
    require(shift>0 and Z>0, 'D2 signs')
    d2candidates=current_candidates(10,134,61,73)
    require(d2candidates==[dict(n=11,ell=1,r=5,u=6,s=11)],
            'D2 complete current-density candidates')
    require(d2candidates[0]['n']*d2candidates[0]['s']!=131,
            'D2 previous degree exclusion')
    d2pivot=first_pivot(11,5,6)
    require(d2pivot==Q(588764,52371555) and d2pivot>0,
            'D2 forced first multiplier')
    require((-1)**(12*6)>0 and beta_leading(62,74)<0,
            'D2 current leading-sign exclusion')
    require(254%10!=0 and 134%10!=0, 'Current densities cannot be initial')
    # The degree-difference proof provides a separate integer contradiction.
    require(134-131==3 and 11*((134-131)-2)!=131,
            'D2 independent degree-jump contradiction')
    for k,D,p,q in [(10,254,111,122),(10,134,61,73)]:
        reflected=current_candidates(k,D,q,p)
        expected=[dict(x, r=x['u'],u=x['r']) for x in current_candidates(k,D,p,q)]
        require(reflected==expected, 'Reflection exclusion')

    receipt={
      'status':'PASS',
      'proof_sha256':sha256((BASE/'canonical_density_exclusion.md').read_bytes()).hexdigest(),
      'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'represented_degree':10,
      'D1':{'previous_degree':252,'current_degree':254,'current_endpoint_orders':[111,122],
            'top_term_ratio':str(ratio),'previous_mean_delta':str(d1difference),
            'complete_current_candidates':d1candidates,'previous_degree_match':False,
            'forced_actual_first_pivot':str(d1pivot),'current_leading_sign_match':False},
      'D2':{'previous_degree':131,'current_degree':134,'current_endpoint_orders':[61,73],
            'normalizer':str(Z),'previous_mean_delta':str(shift),
            'complete_current_candidates':d2candidates,'previous_degree_match':False,
            'forced_actual_first_pivot':str(d2pivot),'current_leading_sign_match':False},
      'scope':'Exact exclusion as canonical density pairs at k=10 for all integer profiles and depths; not a general inverse theorem for truncated moments or an SFC closure.'}
    out=args.output
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','receipt':str(out)}))

if __name__=='__main__':
    main()
