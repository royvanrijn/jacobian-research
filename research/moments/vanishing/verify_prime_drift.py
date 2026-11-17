#!/usr/bin/env python3
"""Bounded exact controls for the all-order prime congruence and drift proof."""
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
        term=poly
        for x,k in zip(variables,powers):term=S.diff(term,x,k)
        result+=c*term
    return S.expand(result)

def iterate(symbol,poly,k,variables,dual):
    assert isinstance(k,int) and 0<=k<=7
    for _ in range(k):poly=action(symbol,poly,variables,dual)
    return poly

def modp(poly,variables,p):
    result=0
    for powers,c in S.Poly(S.expand(poly),*variables).terms():
        c=S.Rational(c);assert c.q%p!=0,'nonintegral coefficient at selected prime'
        reduced=(int(c.p)%p)*pow(int(c.q)%p,-1,p)%p
        result+=reduced*S.prod(x**k for x,k in zip(variables,powers))
    return S.expand(result)

def main():
    x,y,z,X,Y,Z=S.symbols('x y z X Y Z');v=(x,y);d=(X,Y)
    controls=[]
    cases=[(2,x+y,X+Y+X*Y),(3,x*x+y,X+Y+X*X),
           (5,x*y+x,X+2*Y+X*Y),(3,x*y+y*y,S.Rational(1,2)*X+Y+X*Y*Y)]
    for p,P,lam in cases:
        f=iterate(lam,P**p,p,v,d)/math.factorial(p)
        b=S.Poly(lam,*d).coeff_monomial(X)*X+S.Poly(lam,*d).coeff_monomial(Y)*Y
        rhs=action(b,P,v,d)**p
        assert modp(f-rhs,v,p)==0
        controls.append(dict(prime=p,phase=str(P),symbol=str(lam),normalized_reduction=str(modp(f,v,p))))
    u=S.symbols('u');hasse=[];P=x*x+x*y+y
    for p in (2,3,5):
        lhs=S.expand(P.subs(x,x+u)**p).coeff(u,p)
        rhs=S.diff(P,x)**p
        assert modp(lhs-rhs,v,p)==0
        hasse.append(dict(prime=p,equal=True))
    bad_symbol=X+X**2/2;bad_phase=x*x+x
    bad_normalized=iterate(bad_symbol,bad_phase**2,2,v,d)/2
    actual=modp(bad_normalized,v,2);expected=modp(S.diff(bad_phase,x)**2,v,2)
    assert actual==0 and expected==1
    e=S.symbols('e');nv=(x,y,z);nd=(X,Y,Z)
    phase=z+e*(x-y*y/2);lam=X+Y**2
    reduce_e=lambda f:S.expand(S.rem(S.Poly(S.expand(f),e).as_expr(),e**2,e))
    nil=[]
    for j in (1,2,3,4):
        assert reduce_e(action(lam,phase**j,nv,nd))==0
        assert reduce_e(iterate(lam,phase**j,j,nv,nd))==0
        nil.append(j)
    assert S.diff(phase,x)==e
    assert phase.subs(e,0)==z
    # Implicit recursion: h^2+2h-u=0. It introduces only powers of 2
    # into denominators; the general all-degree statement is written.
    h=0;coefficients=[]
    for k in range(1,7):
        known=S.expand(h*h+2*h-u).coeff(u,k)
        coefficient=-known/2;h+=coefficient*u**k
        assert S.denom(coefficient).is_Integer
        denominator=int(S.denom(coefficient));assert denominator>0 and denominator&(denominator-1)==0
        coefficients.append(str(coefficient))
    assert all(S.expand(h*h+2*h-u).coeff(u,k)==0 for k in range(7))
    # A small pure mixed-order pair with passive drift; no all-degree
    # or optimal-cutoff inference is made from these few evaluations.
    passive=[]
    for m in (1,2,3,4):
        assert iterate(X+Y**2,(y*z)**m,m,nv,nd)==0
        passive.append(m)
    boundary=iterate(X+Y**2,x*(y*z)**2,2,nv,nd)
    assert boundary==4*z**2
    guards=0
    for symbol,poly in [(X,S.Float(1)),(X,x**33)]:
        try:action(symbol,poly,v,d)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    try:modp(x/3,v,3)
    except AssertionError:guards+=1
    else:raise AssertionError('nonintegral reduction guard failed')
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['verify_phase_power_quotients.py','verify_formal_binary_gvc.py','verify_algebraic_formal_symbols.py']},
        prime_congruence_controls=controls,hasse_frobenius_controls=hasse,
        uncontrolled_denominator_countercontrol=dict(symbol=str(bad_symbol),phase=str(bad_phase),prime=2,actual=str(actual),false_prediction=str(expected)),
        dual_number_power_controls=nil,dual_number_drift='e',dual_number_phase_reduction='z',
        implicit_algebraic_coefficients=coefficients,passive_drift_pure_controls=passive,
        passive_marker_boundary=str(boundary),guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        universal_prime_congruence_proved_by_samples=False,arbitrary_formal_drift_decided=False,
        all_ring_drift_annihilation_claimed=False,arithmetic_model_recognizer_implemented=False,
        full_ternary_GVC_claimed=False,general_formal_gap_decided=False,
        full_exceptional_MZ_decided=False,unrestricted_frontier_endpoint_changed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (HERE/'prime_drift_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['completed','checker_sha256','passive_marker_boundary','guard_refusals']}))

if __name__=='__main__':main()
