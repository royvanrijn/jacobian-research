#!/usr/bin/env python3
"""Exact degree-case arithmetic for the separate two-height prime corollary."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from pathlib import Path
from math import isqrt
import argparse,hashlib,json
def prime(n):return n>=2 and all(n%q for q in range(2,isqrt(n)+1))
def add(a,b):
    out=dict(a)
    for m,v in b.items():out[m]=out.get(m,0)+v
    return {m:v for m,v in out.items()if v}
def mul(a,b):
    out={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(x+y for x,y in zip(m,n));out[k]=out.get(k,0)+v*w
    return {m:v for m,v in out.items()if v}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();here=Path(__file__).resolve().parent
    primes=[p for p in range(2,13)if prime(p)]
    assert primes==[2,3,5,7,11]
    odd_candidates=[(a,d)for d in primes for a in range(1,d+1,2)if 2*a+4<=d]
    assert odd_candidates==[(1,7),(1,11),(3,11)]
    a3=[]
    for h in range(8):
        for p in range(1,4):
            for j in range(5):
                if h+2*p>7 or h+2*j>9:continue
                if max(h+2*j+2,h+2*p+4)!=11:continue
                compensating=[(palpha,halpha)for palpha in range(1,p+1)
                    for halpha in range(3*palpha+1,h+1)]
                if compensating:a3.append((h,p,j))
    assert a3==[(5,1,0),(5,1,1),(5,1,2)]
    auxiliary=[(3+4*j,5+2*j,max(3+4*j,7+2*j))for j in range(3)]
    assert auxiliary==[(3,5,7),(7,7,9),(11,9,11)]
    equal=[(g,1+6//g)for g in range(1,7)if 6%g==0]
    assert equal==[(1,7),(2,4),(3,3),(6,2)]
    residual=[]
    for g,c in equal:
        if c==2:continue # Explicit low-degree source, not a numerical inference.
        for H in range(4,c+1,2):
            for n in range(c,2*(c-H)+1):
                if 6*c>=5*n+3*H:residual.append((g,c,H,n))
    assert residual==[]
    # Formal polynomial identity for the even-start degree factorization.
    # Coefficients are for1,a,p,a*p.
    one={(0,0):1};A={(1,0):1};P={(0,1):1};Aplus2=add(A,{(0,0):2})
    left=add(add(one,mul(Aplus2,P)),add(A,one)) # 1+(a+2)p+a+1
    right=mul(Aplus2,add(P,one)) # (a+2)(p+1)
    degree_left=[left.get(m,0)for m in((0,0),(1,0),(0,1),(1,1))]
    degree_right=[right.get(m,0)for m in((0,0),(1,0),(0,1),(1,1))]
    assert degree_left==degree_right
    # Existing upper-bound parameters use no primality conjecture.
    a=1;n=2;degree=(2*a+3)*n+a+2
    assert degree==13 and prime(degree)
    proof=here/'TWO_HEIGHT_MINIMUM_UNSAFE_PRIME_THIRTEEN.md'
    ledger_path=here/'two_height_minimum_unsafe_prime_thirteen_source_ledger.json'
    ledger=json.loads(ledger_path.read_text());pins=[]
    for resource in ledger['retained_resources']:
        source=here/resource['retained_path'];got=hashlib.sha256(source.read_bytes()).hexdigest()
        assert got==resource['sha256'];pins.append({'path':resource['retained_path'],'sha256':got,'exact_use':resource['exact_use']})
    out={'passed':True,'scope':'Exhaustive arithmetic of primesbelow13, necessary odd-start candidates, the a3compensating-root degree box and exact auxiliary/equal-degree cases. All-order source theorems, even-start infinity argument and marker safety remain the separate written proof.',
         'primes_below13':primes,'odd_start_candidates':odd_candidates,
         'a3_degree_boxes':a3,'a1_auxiliary_degrees':auxiliary,
         'equal_degree7_quotients':equal,'reflection_residuals':residual,
         'even_degree_polynomial_coefficients':degree_left,
         'existing_F13_parameters':{'a':a,'n':n,'degree':degree},
         'proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest()}
    out.update({'schema_version':1,'status':'passed','retained_inputs':pins,
        'packet_sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest()for q in[proof,Path(__file__).resolve(),ledger_path]},
        'independent_replay':False,'external_review':False,'formal_verification':False})
    path=args.output or here/'two_height_minimum_unsafe_prime_thirteen_verification.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'prime_cases':len(primes),'odd_start_candidates':len(odd_candidates),
                     'a3_degree_boxes':len(a3),'output':str(path)}))
if __name__=='__main__':main()
