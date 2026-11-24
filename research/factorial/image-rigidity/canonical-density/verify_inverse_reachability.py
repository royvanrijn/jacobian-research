#!/usr/bin/env python3
"""Exact controls for finite canonical-density reachability and its inverse."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as S

t,z=S.symbols('t z')
HERE=Path(__file__).resolve().parent
PROOF_NAME='inverse_canonical_reachability.md'
PROOF_SHA256='6345b503c1754630fe8e44cdb162cca672f05933c4b0094234d082cae91ea2ba'

def require(value,message):
    if not bool(value): raise ArithmeticError(message)

def P(expr): return S.Poly(expr,t,domain=S.QQ)

def mass(f,j=0):
    return sum(c/S.Rational(p[0]+j+1) for p,c in f.terms())

def primitive(f):
    return P(sum(c*t**(p[0]+1)/(p[0]+1) for p,c in f.terms()))

def beta(r,u):
    return P(t**r*(1-t)**u*S.factorial(r+u+1)
             /(S.factorial(r)*S.factorial(u)))

def forward(g,f,k):
    require(mass(g)==mass(f)==1,'input normalization')
    a=(k+1)*mass(g,1)-k*mass(f,1)
    K=primitive(primitive(g-f))-primitive(f*P(t-a)).mul_ground(S.Rational(1,k+1))
    b=k*(k+1)*mass(K)
    require(b!=0,'nondefective selected step')
    return K.mul_ground(1/mass(K)),a,b

def endpoints(f):
    a=min(p[0] for p,c in f.terms())
    reversed_f=P(f.as_expr().subs(t,1-t))
    b=min(p[0] for p,c in reversed_f.terms())
    return a,b

def degree(f): return int(f.degree())

def fits(f,n,r,u,ell):
    s=r+u
    expected=n*s if ell==0 else (n+1)*s+2*ell
    if degree(f)!=expected or mass(f)!=1:return False
    factor=P(t**(n*r+ell)*(1-t)**(n*u+ell))
    q,rem=S.div(f,factor)
    return rem.is_zero and degree(q)<=(0 if ell==0 else s)

def backward(g,f,k):
    D=degree(f);d=degree(g);m=k+1
    require(D>=4 and d==D-2,'backward degree pattern')
    b=f.LC();b1=f.coeff_monomial(t**(D-1))
    c=g.LC();c1=g.coeff_monomial(t**(d-1))
    bb=-S.Rational(m*(m+D),D*(D-1))*c/b
    aa=((m+d+1)*c1+bb/m*(D-1)*d*b1)/(d*c)
    h=(g.mul_ground(m+2)+P(t-aa)*g.diff()+f.diff().diff().mul_ground(bb/m)).mul_ground(S.Rational(1,m+1))
    require(bb!=0 and degree(h)<=d-2,'backward top-coefficient cancellation')
    return h,aa,bb

def candidates(g,f,k):
    D=degree(f)
    if D<2:return []
    a0,a1=endpoints(f);delta=D-a0-a1
    cap=(D-2)//(k+2)
    result=[]
    for s in range(delta,cap+1):
        numerator=D-(k+1)*s
        if numerator%(s+2):continue
        ell=numerator//(s+2)
        if ell<1:continue
        n=k+ell
        for r in range(s+1):
            u=s-r
            if fits(f,n,r,u,ell) and fits(g,n,r,u,ell-1):
                result.append((n,r,u,ell))
    return result

def membership(g,f,k):
    require(k>=1 and mass(g)==mass(f)==1,'membership input')
    found=[]
    # Initial pairs are handled separately, including the zero profile.
    D=degree(f)
    if D%k==0:
        s=D//k
        for r in range(s+1):
            u=s-r
            if f==beta(k*r,k*u) and g==beta((k+1)*r,(k+1)*u):
                found.append((k,r,u,0))
    for n,r,u,ell in candidates(g,f,k):
        gg,ff=g,f;kk=k;ok=True
        for j in range(ell,1,-1):
            h,aa,bb=backward(gg,ff,kk)
            if not fits(h,n,r,u,j-2):ok=False;break
            reconstructed,actual_a,actual_b=forward(h,gg,kk+1)
            require(reconstructed==ff and aa==actual_a and bb==actual_b,
                    'inverse operation reconstructs forward kernel')
            ff,gg=gg,h;kk+=1
        if not ok:continue
        lower=beta(n*r,n*u);upper=beta((n+1)*r,(n+1)*u)
        try:
            initial,a,b=forward(upper,lower,n)
        except ArithmeticError as error:
            if str(error)=='nondefective selected step':
                continue
            raise
        if gg==lower and ff==initial:
            found.append((n,r,u,ell))
    return found

def moment(f,k):
    return S.Poly(sum((-1)**j*S.binomial(k,j)*mass(f,j)*z**(k-j)
                      for j in range(k+1)),z,domain=S.QQ)

def actual_control(n,r,u,target_ell):
    g=beta((n+1)*r,(n+1)*u);f=beta(n*r,n*u)
    rows=[]
    for ell in range(target_ell):
        new,a,b=forward(g,f,n-ell)
        require(fits(new,n,r,u,ell+1),'actual exact degree and endpoint factors')
        g,f=f,new
        rows.append({'ell':ell+1,'degree':degree(f),'pivot':str(b)})
    k=n-target_ell
    found=membership(g,f,k)
    require((n,r,u,target_ell) in found,'actual pair recognized')
    return {'n':n,'r':r,'u':u,'ell':target_ell,'k':k,'steps':rows,
            'complete_candidates':[list(x) for x in candidates(g,f,k)],
            'accepted':[list(x) for x in found]}

def invisible_control():
    g=P(1);f=P(6*t*(1-t));k=2
    delta=P(S.diff(t**4*(1-t)**4,t,3))
    displayed=P(24*t*(1-t)*(1-9*t+21*t**2-14*t**3))
    require(delta==displayed,'invisible derivative identity')
    altered=f+delta.mul_ground(S.Rational(1,1000))
    require(all(mass(delta,j)==0 for j in range(k+1)),
            'truncated invisible moments')
    require(moment(altered,k)==moment(f,k),'identical represented quadratic')
    quotient,rem=S.div(delta,P(t*(1-t)))
    require(rem.is_zero and sum(abs(c) for c in quotient.all_coeffs())==1080,
            'positive-density coefficient bound')
    lower_bound=6-S.Rational(1080,1000)
    require(lower_bound==S.Rational(123,25)>0,'strict density positivity')
    require(degree(altered)==5 and endpoints(altered)==(1,1),'altered degree fingerprint')
    original=membership(g,f,k);fake=membership(g,altered,k)
    require(original==[(3,0,0,1)] and not fake,'noninjective membership separation')
    return {'original_canonical_profile':[3,0,0,1],
            'modified_degree':5,'modified_endpoint_orders':[1,1],
            'invisible_moments':[str(mass(delta,j)) for j in range(3)],
            'positive_modifier_lower_bound':str(lower_bound),
            'same_moment_polynomial':str(moment(f,k).as_expr()),
            'modified_pair_canonical_reachable':False,
            'modified_pair_moment_polynomial_reachable':True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=HERE/'inverse_reachability_verification.json')
    args=ap.parse_args()
    proof_digest=sha256((HERE/PROOF_NAME).read_bytes()).hexdigest()
    require(proof_digest==PROOF_SHA256,'proof authentication')
    controls=[actual_control(*x) for x in [(4,2,3,3),(5,1,0,4),(4,0,0,3),(3,0,2,2),(3,1,1,1)]]
    report={'status':'PASS','proof_sha256':proof_digest,'sympy_version':S.__version__,
            'actual_controls':controls,'noninjectivity_control':invisible_control(),
            'scope':'Exact finite canonical-density membership and inverse controls; no SFC closure.'}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','actual_controls':len(controls),'output':str(args.output)}))

if __name__=='__main__':main()
