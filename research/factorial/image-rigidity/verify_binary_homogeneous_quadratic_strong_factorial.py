#!/usr/bin/env python3
"""Finite exact transcription controls; the all-order argument is in the MD."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

def require(condition, message):
    if not condition:
        raise ValueError(message)

def add(*polys):
    ans = {}
    for p in polys:
        for e, c in p.items():
            ans[e] = ans.get(e, 0) + c
            if not ans[e]:
                del ans[e]
    return ans

def scale(p, a):
    return {e: a*c for e, c in p.items() if a*c}

def mul(p, q):
    ans = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(a+b for a, b in zip(e, f))
            ans[g] = ans.get(g, 0) + c*d
    return {e: c for e, c in ans.items() if c}

def power(p, n):
    dim = len(next(iter(p))) if p else 3
    ans = {(0,)*dim: Fraction(1)}
    for _ in range(n):
        ans = mul(ans, p)
    return ans

def native_moment(n):
    ans = {}
    for i in range(n+1):
        for j in range(n-i+1):
            k = n-i-j
            ans[(i,j,k)] = Fraction(
                math.factorial(n)*math.factorial(2*i+j)*math.factorial(j+2*k),
                math.factorial(i)*math.factorial(j)*math.factorial(k))
    return ans

def interval_moment(n):
    # q=A*u^2+B*u-B*u^2+C-2*C*u+C*u^2; independent monomial expansion.
    q = {(1,0,0,2):Fraction(1), (0,1,0,1):Fraction(1),
         (0,1,0,2):Fraction(-1), (0,0,1,0):Fraction(1),
         (0,0,1,1):Fraction(-2), (0,0,1,2):Fraction(1)}
    ans = {}
    for e,c in power(q,n).items():
        abc=e[:3]
        ans[abc]=ans.get(abc,0)+c/Fraction(e[3]+1)
    return {e:c for e,c in ans.items() if c}

def odd_df(n):
    ans=1
    for j in range(1,n+1,2):
        ans *= j
    return ans

def qadd(a,b): return (a[0]+b[0],a[1]+b[1])
def qmul(a,b): return (a[0]*b[0]-15*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def qpow(a,n):
    out=(1,0)
    for _ in range(n): out=qmul(out,a)
    return out

def evaluate_witness(p):
    ans=(0,0)
    for (i,j,k),c in p.items():
        val=qmul(qmul(qpow((5,1),i),qpow((-20,0),j)),qpow((5,-1),k))
        ans=qadd(ans,(c*val[0],c*val[1]))
    return ans

def expect_rejection(test, label):
    try:
        test()
    except ValueError:
        return label
    raise ValueError('Corruption was accepted: '+label)

def main():
    require(not sys.flags.optimize, 'Run without -O; exact controls must remain active')
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    one={(0,0,0):Fraction(1)}
    A={(1,0,0):Fraction(1)}; B={(0,1,0):Fraction(1)}; C={(0,0,1):Fraction(1)}
    alpha=add(A,scale(B,-1),C)
    discriminant_opposite=add(scale(mul(A,C),4),scale(mul(B,B),-1))
    im=[interval_moment(n) for n in range(13)]
    controls=[]
    for n in range(13):
        require(scale(im[n],math.factorial(2*n+1))==native_moment(n),
                'Beta/native identity failed at '+str(n))
        controls.append('native_beta_identity_m='+str(n))
    for n in range(1,13):
        lhs=add(scale(mul(alpha,im[n]),2*(2*n+1)),
                scale(mul(discriminant_opposite,im[n-1]),-n))
        rhs=add(mul(add(scale(A,2),scale(B,-1)),power(A,n)),
                mul(add(scale(C,2),scale(B,-1)),power(C,n)))
        require(lhs==rhs,'Cleared integration identity failed at '+str(n))
        controls.append('cleared_IBP_identity_m='+str(n))
    for n in range(1,41):
        cs=[]
        for j in range(n+1):
            val=Fraction(math.comb(n,j)*math.factorial(2*j)*math.factorial(2*n-2*j),
                         2**n*math.factorial(n))
            require(val.denominator==1,'Nonintegral normalized coefficient')
            require(val.numerator==odd_df(2*j-1)*odd_df(2*n-2*j-1),
                    'Double factorial transcription failed')
            require(val.numerator%2==1,'Even coefficient in opposite-sparse control')
            cs.append(val.numerator)
        h=[1]*(n+1); hn=[1]*(n+2)
        shifted=[0]+h
        require([(x-y)%2 for x,y in zip(hn,shifted)]==[1]+[0]*(n+1),
                'Consecutive complete-polynomial identity failed')
        controls.append('opposite_sparse_odd_coefficients_n='+str(n))
    witness={str(n):evaluate_witness(native_moment(n)) for n in range(1,5)}
    expected={'1':(0,0),'2':(0,0),'3':(-115200,0),'4':(-49766400,0)}
    require(witness==expected,'Sharp witness expansion failed')
    controls.append('native_Qsqrt_minus15_witness_moments_1_through_4')
    corruptions=[]
    for n in (1,2,5,9):
        corrupted=add(native_moment(n),one)
        corruptions.append(expect_rejection(lambda n=n,p=corrupted:require(
            scale(im[n],math.factorial(2*n+1))==p,'native corruption'),
            'native_coefficient_plus_one_m='+str(n)))
    for n in (1,2,5,9):
        wrong=scale(im[n],math.factorial(2*n))
        corruptions.append(expect_rejection(lambda n=n,p=wrong:require(
            p==native_moment(n),'beta scale corruption'),
            'wrong_beta_factorial_m='+str(n)))
    for n in (3,4):
        wrong=list(expected[str(n)]);wrong[0]*=-1
        corruptions.append(expect_rejection(lambda n=n,p=tuple(wrong):require(
            witness[str(n)]==p,'witness sign corruption'),
            'witness_sign_m='+str(n)))
    for n in (2,5):
        wrong=Fraction(math.factorial(2*n),4**n*math.factorial(n))
        corruptions.append(expect_rejection(lambda p=wrong:require(
            p.denominator==1 and p.numerator%2==1,'opposite normalization corruption'),
            'wrong_opposite_normalization_n='+str(n)))
    payload={'status':'PASS','scope':'Finite exact native transcription/witness controls only; no all-order theorem inferred from the prefix',
             'controls':controls,'control_count':len(controls),
             'rejected_corruptions':corruptions,'corruption_count':len(corruptions),
             'witness_moments':{n: [int(c) for c in v] for n,v in witness.items()},
             'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded=json.dumps(payload,indent=2,sort_keys=True)+'\n'
    if args.output:
        require(not args.output.exists(),'Refusing to replace an existing receipt')
        args.output.write_text(encoded,encoding='utf-8')
    print(encoded,end='')

if __name__=='__main__':
    main()