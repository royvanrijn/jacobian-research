#!/usr/bin/env python3
"""Exact arithmetic controls for the written degree-seventeen exclusion.

No moment-prefix inference: coefficient tests replay a forced rational
identity, after the separate mathematical normal-form reduction.
"""
if not __debug__:
    raise RuntimeError('Assertions must be enabled; -O is unsupported')
from pathlib import Path
from fractions import Fraction
import argparse, hashlib, json

HERE=Path(__file__).resolve().parent
ROOT=HERE
PROOF=HERE/'ODD_CUBIC_DEGREE_SEVENTEEN_SAFETY.md'
LEDGER=HERE/'odd_cubic_degree_seventeen_safety_source_ledger.json'

def add(a,b,scale=1):
    out=dict(a)
    for m,v in b.items(): out[m]=out.get(m,0)+scale*v
    return {m:v for m,v in out.items() if v}
def mul(a,b):
    out={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(x+y for x,y in zip(m,n));out[k]=out.get(k,0)+v*w
    return {m:v for m,v in out.items() if v}
def power(a,n):
    out={(0,)*len(next(iter(a))):1}
    for _ in range(n):out=mul(out,a)
    return out
def term(m,v=1):return {m:v}
def scale(a,v):return {m:c*v for m,c in a.items() if c*v}
def coefficient(a,i):return {m[1:]:v for m,v in a.items() if m[0]==i}

def integer_branches():
    pole=[(g,c,r) for g in range(1,15) for c in range(2,20)
          for r in range(1,20) if g*(3*c+2*r-3)==14]
    assert pole==[(1,3,4),(1,5,1),(2,2,2)]
    moving=[(g,c,3*c-4) for g in range(1,15) for c in range(2,20)
            if g*(3*c-5)==14]
    assert moving==[(2,4,8),(14,2,2)]
    fixed=[(g,1+14//g) for g in range(1,15) if 14%g==0]
    assert fixed==[(1,15),(2,8),(7,3),(14,2)]
    # The only not-immediately-low-degree/parity fixed D-dominant case
    # has cbar8, s7, degubar>=8, Hbar>=4. Clear the denominator in Q2.
    assert 6*8 < 5*8+3*4
    candidates=[]
    for g in range(1,17):
        if 16%g:continue
        c=1+16//g
        if c<3:continue
        for e_original in range(1,16,2):
            if (e_original-1)%g:continue
            e=1+(e_original-1)//g
            if e>c-2 or (c-e)%2:continue
            s=(3*c-e-2)//2
            for h in range(4,c+1,2):
                if 5*e<3*c+6*h:continue
                # Exact rank bound and reflection even-degree requirement.
                for n in range(s+1,2*(c-h)+1):
                    if n%2 or 6*c<5*n+3*h:continue
                    candidates.append((g,c,e,h,s,n))
    assert candidates==[(1,17,15,4,17,18)]
    pair_solutions=[]
    # Equality in Q2 is m=2k/3,n=k+1; there are two positive pairs.
    for k in range(1,9):
        ell=9-k
        if k<=ell and k%3==ell%3==0:
            pair_solutions.append(((2*k//3,k+1,k),(2*ell//3,ell+1,ell)))
    assert pair_solutions==[((2,4,3),(4,7,6))]
    return {'D_dominant_poles':pole,'D_dominant_moving_zeros':moving,
            'D_dominant_fixed_zeros':fixed,'C_dominant_residual':candidates,
            'reflection_pairs':pair_solutions}

def reflection_identity():
    # The ordered variables are z,r,s. Compute the common-denominator
    # numerator from the product rule, not an asserted coefficient list.
    z=term((1,0,0));r=term((0,1,0));s=term((0,0,1))
    zr=add(power(z,2),power(r,2),-1)
    zs=add(power(z,2),power(s,2),-1)
    derivative_num=add(add(mul(zr,zs),mul(power(z,2),zs),-6),
                       mul(power(z,2),zr),-12)
    expected={(4,0,0):-17,(2,2,0):11,(2,0,2):5,(0,2,2):1}
    assert derivative_num==expected
    reciprocal_num=mul(power(add(z,r,-1),2),power(add(z,s,-1),3))
    even={m:v for m,v in reciprocal_num.items() if m[0]%2==0}
    expected_even={(4,1,0):-2,(4,0,1):-3,(2,2,1):-3,
                   (2,1,2):-6,(2,0,3):-1,(0,2,3):-1}
    assert even==expected_even
    # Constants fix K=-A*s. Residual=(-s)*derivative_num-even.
    residual=add(mul(scale(s,-1),derivative_num),even,-1)
    assert coefficient(residual,0)=={}
    assert coefficient(residual,4)=={(1,0):2,(0,1):20}
    assert coefficient(residual,2)=={(2,1):-8,(1,2):6,(0,3):-4}
    # Dividing by 2 and -2s, then setting t=r/s yields the forced
    # equations L=t+10 and q=4t²-3t+2. Replay a literal Bezout identity.
    L={(1,):1,(0,):10};q={(2,):4,(1,):-3,(0,):2}
    certificate=add(q,mul({(1,):4,(0,):-43},L),-1)
    assert certificate=={(0,):432}
    # A changed constant must fail the SAME exact identity.
    corrupted=dict(q);corrupted[(0,)]+=1
    assert add(corrupted,mul({(1,):4,(0,):-43},L),-1)!={(0,):432}
    return {'product_rule_derivative_terms':len(derivative_num),
            'even_reciprocal_terms':len(even),'literal_nonzero_constant':432,
            'corrupt_coefficient_rejected':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    out={'passed':True,'scope':'Arithmetic controls for the separate universal proof. Exhaustive degree17 quotient integer branches and the entire forced residual rational derivative identity. No finite moment-prefix assertion, no classification beyond the stated odd-cubic class.',
         'integer_branches':integer_branches(),'exact_reflection':reflection_identity()}
    paths=[PROOF,Path(__file__).resolve(),LEDGER]
    if LEDGER.exists():
        ledger=json.loads(LEDGER.read_text())
        for resource in ledger['retained_resources']:
            path=(ROOT/resource['retained_path']) if resource.get('location')!='scratch' else HERE/resource['retained_path']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==resource['sha256']
            paths.append(path)
    out['schema_version']=1;out['status']='passed';out['independent_replay']=False;out['external_review']=False;out['formal_verification']=False
    out['pins']=[{'path':str(p.relative_to(HERE)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
    path=args.output or HERE/'odd_cubic_degree_seventeen_safety_verification.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'C_dominant_residual_cells':len(out['integer_branches']['C_dominant_residual']),
                      'literal_nonzero_constant':432,'output':str(path)}))
if __name__=='__main__':main()
