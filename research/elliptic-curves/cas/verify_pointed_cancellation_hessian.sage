#!/usr/bin/env sage -python
"""Covariant cancellation corollary; no model selection or point search.

Use the established pointed-quartic Hessian covering map to reconstruct N,D,
then verify exact derivative valuations and the fixed gcd squareclass on all
retained V1 signed target/model observations. This does not change V2 policy.
"""
from collections import Counter
from fractions import Fraction as F
import gzip
import json
from math import gcd,isqrt
from pathlib import Path
import time

from sage.all import PolynomialRing,QQ
from finite_cancellation_corpus import OUT as V1,canonical,digest,write
from finite_cancellation_validation_audit import OUT


def hessian(q):
    e,d,c,b,a=q
    return [24*c*e-9*d*d,72*b*e-12*c*d,144*a*e+18*b*d-12*c*c,
            72*a*d-12*b*c,24*a*c-9*b*b]


def val(x,p):
    if not x:return 10**9
    k=0
    while x%p==0:x//=p;k+=1
    return k


def universal():
    R=PolynomialRing(QQ,names=('a','b','A','u','v'));a,b,A,u,v=R.gens();B=b*b-a**3-A*a
    d=u**4-6*a*u*u*v*v-8*b*u*v**3-(3*a*a+4*A)*v**4
    n=a*u**4+4*b*u**3*v+(6*a*a+4*A)*u*u*v*v+4*a*b*u*v**3+(a**3+4*B)*v**4
    H=d.derivative(u,2)*d.derivative(v,2)-d.derivative(u).derivative(v)**2
    assert H==-144*n
    S=PolynomialRing(QQ,names=('e','d','c','b','a','t'));e,d,c,b,a,t=S.gens()
    q=e+d*t+c*t*t+b*t**3+a*t**4
    hh=sum(z*t**i for i,z in enumerate(hessian([e,d,c,b,a])))
    assert hh==12*q*q.derivative(t,2)-9*q.derivative(t)**2
    return {'raw_Hessian_equals_minus_144_N':True,'affine_Hessian_equals_12_q_qsecond_minus_9_qprime_squared':True}


def main():
    start=time.process_time();identities=universal();count=Counter();mus=Counter();hs=Counter();hashes={}
    primes=(2,3,5,7,11,13,17,19,23,29,31)
    for path in sorted((V1/'cases').glob('*.json.gz')):
        raw=path.read_bytes();hashes[path.name]=digest(raw);packet=json.loads(gzip.decompress(raw))
        assert packet['status']=='PASS_EXACT_ACCESSIBILITY'
        prep=packet['prepared'];assert all(F(x).denominator==1 for x in prep['curve'])
        data=[]
        for model in prep['models']:
            q=[int(v) for v in model['q']];H=hessian(q)
            a,b,c,d=map(F,model['mapping']['matrix']);rho=F(model['mapping']['square_ratio'])
            mu=rho/(144*(a*d-b*c)**2);alpha,beta=mu.numerator,mu.denominator
            assert alpha>0 and isqrt(alpha)**2==alpha and isqrt(beta)**2==beta
            pair=[-alpha*z for z in H]+[beta*z for z in q];h=gcd(*pair)
            assert [z//h for z in pair]==list(map(int,model['N']+model['D']))
            data.append((q,alpha,beta,h));mus[str(mu)]+=1
            hs['square' if isqrt(h)**2==h else 'nonsquare']+=1
            count['models']+=1
        assert len({F(a,b) for q,a,b,h in data})==1
        for row in packet['observations']:
            q,alpha,beta,h=data[row['model_index']];m,n=map(int,row['coordinate']);g=int(row['g'])
            qq=sum(q[i]*m**i*n**(4-i) for i in range(5))
            assert qq>0 and isqrt(qq)**2==qq
            hv=sum(c*m**i*n**(4-i) for i,c in enumerate(hessian(q)))
            assert gcd(-alpha*hv,beta*qq)//h==g
            assert isqrt(h*g)**2==h*g
            dm=sum(i*q[i]*m**(i-1)*n**(4-i) for i in range(1,5))
            dn=sum((4-i)*q[i]*m**i*n**(3-i) for i in range(4))
            for p in primes:
                if (6*alpha*beta)%p==0:
                    count['exceptional_prime_locations']+=1;continue
                derivative=dm if n%p else dn
                expected=min(val(qq,p),2*val(derivative,p))-val(h,p)
                assert expected==val(g,p)
                count['derivative_valuation_checks']+=1
            count['squareclass_checks']+=1
        count['anchors']+=1
        if count['anchors']%1000==0:print(json.dumps({'anchors':count['anchors'],'cpu_seconds':time.process_time()-start}),flush=True)
    result={'status':'PASS_POINTED_HESSIAN_CANCELLATION_COROLLARY', 'universal_identities':identities,
        'counts':dict(count),'mu_counts':dict(mus),'normalization_content_squareclasses':dict(hs),
        'input_cases_sha256':hashes,'source_sha256':digest(Path(__file__).read_bytes()),'cpu_seconds':time.process_time()-start,
        'boundary':'Established Hessian covering identity, applied to cancellation. For mu=alpha/beta and h=content(-alpha Hess(q),beta q), the primitive pair is exactly (-alpha Hess(q),beta q)/h. Reconstructed the full gcd at every target, including all exceptional and unprocessed primes. On integral short curves h*g is a square for every rational point in the patch. At p not dividing 6 alpha beta, v_p(g)=min(v_p(q),2v_p(qprime))-v_p(h), using the reciprocal chart at infinity. Fixed exceptional primes remain outside the derivative shortcut. The initial certificate is preserved separately; this successor adds full-gcd reconstruction. No C/L optimization, target-based policy tuning, new points or additional CPU arms.'}
    write(OUT/'hessian-cancellation.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='input_cases_sha256'},indent=2))


if __name__=='__main__':main()
