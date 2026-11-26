#!/usr/bin/env python3
"""Supplementary exact transcriptions and finite criterion frontier; not a universal-prefix proof."""
if not __debug__:
    raise SystemExit('Optimized execution is refused; use ordinary Python 3.')
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path

def check(ok, message):
    if not ok:
        raise ValueError(message)

def add(*polys):
    out={}
    for poly in polys:
        for e,c in poly.items():
            out[e]=out.get(e,0)+c
    return {e:c for e,c in out.items() if c}

def scale(poly,c):
    return {e:c*v for e,v in poly.items() if c*v}

def mul(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            g=tuple(x+y for x,y in zip(e,f))
            out[g]=out.get(g,0)+c*d
    return {e:c for e,c in out.items() if c}

def power(p,n):
    out={(0,0,0):F(1)}
    for _ in range(n):
        out=mul(out,p)
    return out

def native_abc(j):
    return {(a,b,j-a-b):F(factorial(j)*factorial(2*a+b),
             factorial(a)*factorial(b)*factorial(j-a-b))
             for a in range(j+1) for b in range(j-a+1)}

def native_fixed(j,C):
    # First choose x^2 factors, then x factors, and apply Lambda(x^r)=r!.
    return [F(comb(j,i)*sum(comb(j-i,t)*factorial(2*i+t)*C**(j-i-t)
              for t in range(j-i+1))) for i in range(j+1)]

def vp(x,p):
    x=F(x)
    if not x:
        return None
    def val(n):
        n=abs(n);r=0
        while n%p==0:
            n//=p;r+=1
        return r
    return val(x.numerator)-val(x.denominator)

def residue(x,p):
    x=F(x)
    check(x.denominator%p!=0,'Nonintegral residue')
    return x.numerator*pow(x.denominator,-1,p)%p

def lower_hull(coefficients,p):
    hull=[]
    for i,c in enumerate(coefficients):
        if not c:
            continue
        pt=(i,vp(c,p))
        while len(hull)>=2:
            a,b=hull[-2];u,v=hull[-1];x,y=pt
            if (v-b)*(x-u)>=(y-v)*(u-a):
                hull.pop()
            else:
                break
        hull.append(pt)
    return hull

def factor(n):
    out={};p=2
    while p*p<=n:
        if n%p==0:
            t=0
            while n%p==0:
                n//=p;t+=1
            out[p]=t
        p+=1
    if n>1:
        out[n]=1
    return out

def digits(n,p):
    out=[]
    while n:
        out.append(n%p);n//=p
    return out or [0]

def gate(q,p,t):
    z=1;L=0
    while z*p<=2*q:
        z*=p;L+=1
    B=t*p**t-(p**t-1)//(p-1)
    return {'p':p,'t':t,'L':L,'B':B,'digits_low_first':digits(q,p),
            'admissible':p%2==1 and all(x<=(p-1)//2 for x in digits(q,p)) and L<B}

def remaining_data(q):
    bits=[i for i in range((q-1).bit_length()) if (q-1)>>i&1]
    if q%2==0 or len(bits)<4:
        return {'q':q,'remaining':False}
    D=2**bits[3];J=1;gs=[]
    for p,t in factor(q).items():
        g=gate(q,p,t);gs.append(g)
        if g['admissible']:
            J*=p**t
    return {'q':q,'D':D,'J':J,'prime_gates':gs,'remaining':D*J<=q-2}

def rejected(fn,label):
    try:
        fn()
    except ValueError:
        return label
    raise ValueError('Corruption accepted: '+label)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,help='Write a fresh receipt; existing files are never overwritten.')
    args=ap.parse_args()
    controls=[];negative=[]
    A={(1,0,0):F(1)};B={(0,1,0):F(1)};C={(0,0,1):F(1)}
    disc=add(mul(B,B),scale(mul(A,C),-4))
    for j in range(2,11):
        core=add(scale(mul(A,native_abc(j-1)),2*j*(2*j-1)),
                 scale(mul(disc,native_abc(j-2)),j*(j-1)))
        plus=add(core,power(C,j),scale(mul(B,power(C,j-1)),j))
        check(plus==native_abc(j),'Symbolic native PLUS recurrence failed')
        controls.append('native_abc_recurrence_j='+str(j))
        if j in (2,3,5,8):
            wrong=add(core,power(C,j),scale(mul(B,power(C,j-1)),-j))
            negative.append(rejected(lambda w=wrong,j=j:check(w==native_abc(j),'wrong forcing sign'),
                                     'forcing_minus_j='+str(j)))
    ps=[[1],[1,1]]
    for j in range(2,9):
        row=[0]+ps[-1]
        for i,c in enumerate(ps[-2]):
            row[i]^=c
        ps.append(row)
    for c0 in (-4,0,2):
        for j in range(9):
            row=[v/F(2**i*factorial(j)) for i,v in enumerate(native_fixed(j,c0))]
            check([residue(v,2) for v in row]==ps[j],'Even-constant reduction failed')
            controls.append('even_constant_C='+str(c0)+'_j='+str(j))
    height=[];coefficient_controls=0
    for q,p,t in ((3,3,1),(9,3,2),(55,5,1),(63,7,1),(87,29,1),(135,5,1)):
        check(gate(q,p,t)['admissible'],'Expected admissible example failed')
        f=native_fixed(q,-q)
        hurwitz=[f[q-j]*F(factorial(j),factorial(q)**2) for j in range(q+1)]
        ordinary=[a/factorial(j) for j,a in enumerate(hurwitz)]
        for j,a in enumerate(hurwitz):
            direct=sum((F(comb(j,k)*(-q)**k*factorial(2*q-j-k),
                          factorial(q)*factorial(q-j)) for k in range(j+1)),F(0))
            check(a==direct,'Native Hurwitz transcription differs')
            check(not a or vp(a,p)>=0,'Hurwitz coefficient nonintegral in admissible example')
            coefficient_controls+=1
        prefix=[0];x=0;ds=digits(q,p)
        for s in range(len(ds)-1,-1,-1):
            if ds[s]:
                x+=ds[s]*p**s;prefix.append(x)
        check(all(vp(hurwitz[j],p)==0 for j in prefix),'Prefix unit failed')
        expected=[(j,-vp(factorial(j),p)) for j in prefix]
        observed=lower_hull(ordinary,p)
        check(observed==expected,'Full native prefix polygon differs')
        height.append({'q':q,**gate(q,p,t),'computed_polygon':observed})
        controls.append('native_height_prefix_q='+str(q)+'_p='+str(p))
        wrong=ordinary[:];wrong[0]=F(0)
        negative.append(rejected(lambda w=wrong,p=p,e=expected:check(lower_hull(w,p)==e,'deleted endpoint'),
                                 'deleted_height_endpoint_q='+str(q)))
    binary=[]
    for m in (2,6,14,30,54,62,90):
        bits=[r for r in range(m.bit_length()) if m>>r&1]
        prefix=[(0,0)];x=0
        for h,r in enumerate(bits,1):
            x+=2**r;prefix.append((x,2*x-h))
        coeff=native_fixed(m,-7)
        check(lower_hull(coeff,2)==prefix,'Binary coefficient polygon differs')
        for i,c in enumerate(coeff):
            bound=2*i+(m-i).bit_count()-m.bit_count()
            check(not c or vp(c,2)>=bound+(i%2),'Binary coefficient bound differs')
            if i%2==0:
                check(vp(c,2)==bound,'Even-index equality differs')
        binary.append({'m':m,'computed_polygon':prefix})
        controls.append('native_binary_polygon_m='+str(m))
    fourth=F(-1)+F(2-3)+F(14-2,16)
    check(fourth==F(-5,4),'Fourth-cluster boundary changed')
    negative.append(rejected(lambda:check(fourth>F(-1),'false fourth-cluster exclusion'),
                             'overstrong_fourth_cluster_closure'))
    frontier=[q for q in range(3,92) if remaining_data(q)['remaining']]
    check(frontier==[91],'Displayed R-set frontier differs')
    odd_composites=[q for q in range(55,91,2)
                    if len(factor(q))>1 or next(iter(factor(q).values()))>1]
    check(odd_composites==[55,57,63,65,69,75,77,81,85,87],'Composite list differs')
    examples=[remaining_data(q) for q in (55,63,87,91)]
    check([(r['D'],r['J']) for r in examples]==[(32,55),(16,7),(64,29),(64,1)],'Product degree arithmetic differs')
    controls.append('exact_displayed_R_set_first_q91_through_91')
    for wrong,label in (([87],'wrong_first_q87'),([], 'mutated_displayed_R_empty')):
        negative.append(rejected(lambda w=wrong:check(frontier==w,'incorrect frontier'),label))
    check(all(k+vp(factorial(21-k),5)-4>=0 for k in range(22)),
          'Counterfamily finite head inequality differs')
    controls.append('counterfamily_22_small_head_inequalities_only')
    out={'status':'PASS','scope':'Supplementary finite exact native identities and explicitly finite frontier arithmetic only. Universal moment/degree assertions are proved in the written source, not inferred from a tested prefix. No enormous counterfamily coefficient is expanded.',
         'controls':controls,'control_count':len(controls),'native_hurwitz_coefficients_crosschecked':coefficient_controls,
         'height_examples':height,'binary_examples':binary,
         'first_remaining_under_displayed_criterion_q':91,'first_start_under_displayed_criterion':89,
         'criterion_examples':examples,'odd_composites_55_leq_q_lt91':odd_composites,
         'rejected_corruptions':negative,'corruption_count':len(negative),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.output:
        with args.output.open('x',encoding='utf-8') as f:
            f.write(encoded)
    print(encoded,end='')

if __name__=='__main__':
    main()