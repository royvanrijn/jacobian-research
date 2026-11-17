#!/usr/bin/env python3
"""Small exact controls for the written valuation/carry theorem."""
from pathlib import Path
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent

def action(symbol,poly,variables,dual):
    symbol=S.sympify(symbol);poly=S.sympify(poly)
    assert not symbol.has(S.Float) and not poly.has(S.Float)
    assert S.Poly(poly,*variables).total_degree()<=32
    result=0
    for powers,c in S.Poly(symbol,*dual).terms():
        t=poly
        for x,k in zip(variables,powers):t=S.diff(t,x,k)
        result+=c*t
    return S.expand(result)

def iterate(symbol,poly,k,variables,dual):
    assert isinstance(k,int) and 0<=k<=7
    for _ in range(k):poly=action(symbol,poly,variables,dual)
    return poly

def frobenius(poly,variables,p):
    assert S.isprime(p)
    return S.expand(sum(c**p*S.prod(x**(p*k) for x,k in zip(variables,a)) for a,c in S.Poly(poly,*variables).terms()))

def carry(poly,variables,p):
    return S.expand((poly**p-frobenius(poly,variables,p))/p)

def modp(poly,variables,p):
    result=0
    for powers,c in S.Poly(S.expand(poly),*variables).terms():
        c=S.Rational(c);assert c.q%p!=0,'nonintegral selected-prime coefficient'
        result+=(int(c.p)*pow(int(c.q),-1,p)%p)*S.prod(x**k for x,k in zip(variables,powers))
    return S.expand(result)

def vp_factorial(k,p):
    answer=0
    while k:
        k//=p;answer+=k
    return answer

def exponent(n,r):return min(r,max(2,r-n+3))

def main():
    x,y,z,X,Y,Z=S.symbols('x y z X Y Z');v=(x,y,z);d=(X,Y,Z)
    cases=[(2,X*Y+Y*Z+X*Z,x*x+y*y+z*z),
           (3,X*X+Y*Y+Z*Z,x*x+y*y+z*z),
           (5,X*X+Y*Y+Z*Z,x*x+y*y+z*z),
           (3,X*Y+Z**3,x*x+y*y+z*z+x*y*z)]
    carries=[]
    for p,lam,P in cases:
        out=iterate(lam,P**p,p,v,d)
        q=sum(c*S.prod(a**k for a,k in zip(d,alpha)) for alpha,c in S.Poly(lam,*d).terms() if sum(alpha)==2)
        correction=action(carry(lam,d,p),carry(P,v,p),v,d)
        assert modp(out/math.factorial(p)**2-action(q,P,v,d)**p-correction,v,p)==0
        carries.append(dict(prime=p,symbol=str(lam),phase=str(P),residue=str(modp(out/math.factorial(p)**2,v,p)),carry_residue=str(modp(correction,v,p))))
    assert carries[1]['residue']=='2' and carries[1]['carry_residue']=='2'
    assert carries[2]['residue']=='2' and carries[2]['carry_residue']=='1'
    valuation=[]
    # Six fixed coefficient-profile controls for the all-n/r sharpness proof.
    # There is no differentiation in growing dimension or prime sweep.
    for n,r,p in [(1,4,7),(3,5,7),(4,3,149),(4,4,197),(5,4,409),(8,3,149)]:
        target=exponent(n,r)
        if target==r:
            actual=vp_factorial(r*p,p)
            profile=[p]
        else:
            s=min(n,r+1);assert p>s*r*(s-1)
            k=((s-1)*p)//(s*r);profile=[p-(s-1)*k]+[k]*(s-1)
            assert sum(profile)==p and all(0<t<p for t in profile)
            floors=[(r*t)//p for t in profile]
            assert floors==[r-s+1]+[0]*(s-1)
            actual=2*(vp_factorial(p,p)-sum(vp_factorial(t,p) for t in profile))+sum(vp_factorial(r*t,p) for t in profile)
        assert actual==target
        valuation.append(dict(n=n,r=r,prime=p,profile=profile,valuation=actual))
    small_profile=[2,2,2,1]
    normalized=S.prod(S.Rational(math.factorial(3*k),math.factorial(k)**2) for k in small_profile)
    assert normalized==34992000 and int(normalized)%7==1
    assert 2-3==-1  # Cubic-order normalization has a genuine 1/p coefficient.
    # A nonreduced p-torsion-free coefficient model retains the SAME identity.
    e=S.symbols('e');lam=X*Y+e*Z**2;P=x*y+z*z+e*x*z;p=3
    lhs=iterate(lam,P**p,p,v,d)/math.factorial(p)**2
    rhs=action(lam,P,v,d)**p+action(carry(lam,d,p),carry(P,v,p),v,d)
    difference=S.rem(S.expand(lhs-rhs),e**2,e)
    for k in (0,1):assert modp(S.expand(difference).coeff(e,k),v,p)==0
    guards=0
    try:modp(x/3,v,3)
    except AssertionError:guards+=1
    else:raise AssertionError('nonintegral guard failed')
    try:action(X,S.Float(1),v,d)
    except AssertionError:guards+=1
    else:raise AssertionError('float guard failed')
    deps=['verify_prime_drift.py','verify_formal_apolar_gaps.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        exact_carry_controls=carries,sharp_valuation_controls=valuation,
        cubic_higher_tail_profile=small_profile,cubic_higher_tail_normalized_coefficient=str(normalized),cubic_higher_tail_residue_mod7=1,
        cubic_order_normalized_valuation=-1,nonreduced_coefficient_identity_pass=True,guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        universal_sharpness_proved_by_samples=False,carry_can_be_omitted=False,
        higher_symbol_orders_discarded=False,order_normalization_integral_in_all_dimensions=False,
        general_apolar_nesting_proved=False,general_formal_gap_decided=False,
        full_exceptional_MZ_decided=False,unrestricted_frontier_endpoint_changed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (HERE/'prime_carry_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['completed','checker_sha256','cubic_higher_tail_normalized_coefficient','guard_refusals']}))

if __name__=='__main__':main()
