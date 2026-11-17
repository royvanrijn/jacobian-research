#!/usr/bin/env python3
"""Small exact controls for universal written height/degree and torus-obstruction proofs."""
from pathlib import Path
import json,hashlib,math
import sympy as S

def lcm_to(n):
    result=1
    for k in range(1,n+1):result=math.lcm(result,k)
    return result

def occupations(m,n):
    if n==1:yield (m,);return
    for i in range(m+1):
        for tail in occupations(m-i,n-1):yield (i,)+tail

def gaussian_moment(poly,variables):
    out=0
    for powers,c in S.Poly(S.expand(poly),*variables).terms():
        if any(a%2 for a in powers):continue
        out+=c*S.prod(S.factorial(a)/(2**(a//2)*S.factorial(a//2)) for a in powers)
    return S.simplify(out)

def sphere_moment(poly,variables):
    N=len(variables);out=0
    for powers,c in S.Poly(S.expand(poly),*variables).terms():
        if any(a%2 for a in powers):continue
        half=[a//2 for a in powers];m=sum(half)
        out+=c*S.prod(S.rf(S.Rational(1,2),a) for a in half)/S.rf(S.Rational(N,2),m)
    return S.simplify(out)

def normalized_marked_moment(P,Q,variables,n,max_total_degree=24):
    P=S.sympify(P);Q=S.sympify(Q)
    if any(q.atoms(S.Float) for q in [P,Q]):raise ValueError('exact coefficients required')
    if not isinstance(n,int) or n<1:raise ValueError('positive integer moment index required')
    terms=[(a,c) for a,c in S.Poly(P,*variables).terms() if c!=0];degrees={sum(a) for a,c in terms}
    if len(degrees)!=1 or next(iter(degrees))<1:raise ValueError('nonzero positive-degree homogeneous phase required')
    D=next(iter(degrees));E=S.Poly(Q,*variables).total_degree()
    if D*n+E>max_total_degree:raise ValueError('configured exact Wick degree limit;no asymptotic inference')
    N=len(variables);k=(D*n+1)//2;epsilon=(D*n)%2
    original=gaussian_moment(Q*P**n,variables);radial=2**k*S.rf(S.Rational(N,2),k)
    pieces={}
    for powers,c in S.Poly(Q,*variables).terms():
        e=sum(powers);pieces[e]=pieces.get(e,0)+c*S.prod(x**a for x,a in zip(variables,powers))
    sphere=0
    for e,q in pieces.items():
        if e%2!=epsilon:continue
        j=(e-epsilon)//2
        sphere+=2**j*S.rf(S.Rational(N,2)+k,j)*sphere_moment(q*P**n,variables)
    assert S.simplify(original/radial-sphere)==0
    return dict(phase_degree=D,marker_degree=int(E),n=n,k=k,moment=str(original),radial=str(radial),normalized=str(S.simplify(sphere)))

def exact_factorial_degree_gate(degree,slope,shift=0):
    if not isinstance(degree,int) or degree<1 or not isinstance(slope,int) or slope<0 or not isinstance(shift,int) or shift<0:raise ValueError('integer degree/slope/shift required')
    return dict(necessary_degree=2*slope,compatible=degree==2*slope,scope='necessary only for an exact all-index marked factorial law;law itself not verified')

def balanced_pair(a,b):
    if not isinstance(a,int) or not isinstance(b,int) or a<=0 or b>=0:raise ValueError('opposite nonzero integer torus weights required')
    g=math.gcd(a,-b);powers=(-b//g,a//g);assert a*powers[0]+b*powers[1]==0
    return dict(weights=[a,b],balanced_product_powers=list(powers),consequence='nonzero torus-fixed product;all-power joint purity contradicts retained neutral radical zero')

def main():
    here=Path(__file__).parent;denom=[]
    for N in [1,2,3,4,5]:
        count=0
        for m in range(4):
            bound=4**m*lcm_to(2*m+N)**(N+2)
            for alpha in occupations(m,N):
                q=S.prod(S.rf(S.Rational(1,2),a) for a in alpha)/S.rf(S.Rational(N,2),m)
                assert S.denom(S.cancel(bound*q))==1;count+=1
        denom.append(dict(dimension=N,through_half_degree=3,even_monomials_checked=count))
    lcm=[]
    for k in range(1,9):
        left=lcm_to(2*k);right=lcm_to(k)*int(S.binomial(2*k,k));assert right%left==0 and lcm_to(k)<=16**k
        lcm.append(dict(k=k,doubling_divisibility=True))
    x,y=S.symbols('x y');radial=[]
    for P,Q,variables in [(x**2,1+x*x,(x,)),(x**3,x+x*x,(x,)),(x*x+S.sqrt(2)*y*y,1,(x,y))]:
        for n in [1,2,3]:radial.append(normalized_marked_moment(P,Q,variables,n))
    sharp=[]
    for ell in [1,2,3]:
        P=((x*x+y*y)/2)**ell
        for n in [1,2]:
            q=normalized_marked_moment(P,1,(x,y),n);assert S.sympify(q['moment'])==S.factorial(ell*n)
            sharp.append(dict(slope=ell,n=n,phase_degree=2*ell,marked=q['moment']))
    norm=[]
    for n in [1,2,3]:
        q=normalized_marked_moment(x*x+S.sqrt(2)*y*y,1,(x,y),n);value=S.sympify(q['normalized'])
        conjugate=value.xreplace({S.sqrt(2):-S.sqrt(2)});norm_=S.simplify(value*conjugate)
        k=n;T=4**k*lcm_to(2*k+2)**4;assert S.denom(S.cancel(T*T*norm_))==1 and norm_!=0
        norm.append(dict(n=n,normalized=str(value),field_norm=str(norm_),cleared_norm_is_nonzero_integer=True))
    balances=[balanced_pair(a,b) for a,b in [(2,-2),(8,-12),(4,-10)]]
    gates=[exact_factorial_degree_gate(D,ell) for D,ell in [(4,2),(16,2),(24,2),(6,3),(3,1)]]
    assert [q['compatible'] for q in gates]==[True,False,False,True,False]
    t=S.Symbol('t');hilbert=[]
    for name,numerator,degrees,minimum in [('2O',36,[8,12,18],8),('2I',60,[12,20,30],12)]:
        series=S.series((1-t**numerator)/S.prod(1-t**q for q in degrees),t,0,minimum+1).removeO().expand()
        assert all(series.coeff(t,j)==0 for j in range(1,minimum)) and series.coeff(t,minimum)==1
        hilbert.append(dict(group=name,minimum_nonconstant_invariant_degree=minimum,all_ambient_markers_below_threshold_are_radial_after_Reynolds='retained invariant-degree and homogeneous-lift theorem'))
    rejects=0
    for P in [1+x,x+S.Float('1.0')*x*x]:
        try:normalized_marked_moment(P,1,(x,),1)
        except ValueError:rejects+=1
    try:normalized_marked_moment(x**4,1,(x,),7)
    except ValueError:rejects+=1
    assert rejects==3
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        sphere_denominator_controls=denom,lcm_controls=lcm,odd_even_radial_controls=radial,sharp_homogeneous_factorial_examples=sharp,
        number_field_norm_controls=norm,balanced_weight_controls=balances,exact_factorial_degree_gates=gates,exceptional_low_marker_degree_controls=hilbert,guard_refusals=rejects,
        normalized_fixed_marker_calculator_implemented=True,exact_necessary_degree_gate_implemented=True,balanced_weight_certificate_implemented=True,supporting_checker_hashes={},
        universal_growth_theorem_proved_by_finite_controls=False,numerical_asymptotic_fit_used=False,complex_coefficient_pointwise_lower_bound_claimed=False,
        full_exceptional_MZ_decision_implemented=False,general_invariance_or_marker_solver_implemented=False,new_exceptional_counterexample_constructed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (here/'homogeneous_moment_marker_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','sphere_denominator_controls','balanced_weight_controls','exact_factorial_degree_gates']}))

if __name__=='__main__':main()
