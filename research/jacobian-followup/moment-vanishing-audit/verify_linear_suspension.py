#!/usr/bin/env python3
"""Small exact controls for the written suspension/nesting classification."""
from pathlib import Path
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent
x,y,u,v,X,Y,a,b,t,z,Q=S.symbols('x y u v X Y a b t z Q')

def action(lam,f):
    lam=S.sympify(lam);f=S.sympify(f)
    assert not lam.has(S.Float) and not f.has(S.Float)
    assert S.Poly(f,x,y,u).total_degree()<=24
    return S.expand(sum(c*S.diff(f,x,alpha[0],y,alpha[1]) for alpha,c in S.Poly(lam,X,Y).terms()))

def iterate(lam,f,k):
    assert isinstance(k,int) and 0<=k<=6
    for _ in range(k):f=action(lam,f)
    return f

def height(lam,f):
    assert S.Poly(lam,X,Y).coeff_monomial(1)==0
    k=0
    while True:
        f=action(lam,f)
        if f==0:return k
        k+=1
        assert k<=24

def coefficients(f,vars):return [c for alpha,c in S.Poly(S.expand(f),*vars).terms() if c!=0]

def main():
    counter=[]
    for d in [3,4]:
        P=x**d+y**d;C=S.binomial(3*d,d);H=C-3*S.binomial(2*d,d)+3
        R=3*(X**(2*d)+Y**(2*d))-(3+C)*X**d*Y**d
        assert H>0 and action(R,P**3)==0
        assert action(R,P**2)==-2*math.factorial(d)**2*H
        assert action(R,x**d*P)==-math.factorial(d)**2*H
        hs=[height(R,P**j) for j in range(1,4)]
        lifted=[height(R,(P+u)**j) for j in range(1,4)]
        assert hs==[0,1,0] and lifted==[0,1,1]
        assert action(R,(P+u)**3)==-6*math.factorial(d)**2*H*u
        counter.append(dict(degree=d,symbol=str(R),H=int(H),square_output=str(action(R,P**2)),cube_output='0',sharp_marker_output=str(action(R,x**d*P)),heights=hs,lifted_heights=lifted,sharp_cutoff=str(S.Rational(1,d))))
    # One all-degree counting identity is checked at a small exact instance.
    # Universal H_d>0 is the written three-block inclusion-exclusion proof.
    assert S.binomial(9,3)-3*S.binomial(6,3)+3==27
    binomials=0
    for lam,P in [(X**2+Y**3,x*y+x),(X**2+Y**4,x*y+y*y)]:
        for m,k in [(2,1),(3,2),(3,3)]:
            lhs=iterate(lam,(P+u)**m,k)
            rhs=sum(S.binomial(m,j)*u**(m-j)*iterate(lam,P**j,k) for j in range(m+1))
            assert S.expand(lhs-rhs)==0;binomials+=1
    # Ordinary prefix ideals agree and genuinely retain a square-zero direction.
    lam=a*X**2+b*Y;P=x*x+y
    base=[];lift=[]
    for m in [1,2]:
        base+=coefficients(iterate(lam,P**m,m),(x,y))
        lift+=coefficients(iterate(lam,(P+u)**m,m),(x,y,u))
    gb=S.groebner(base,b,a);gl=S.groebner(lift,b,a)
    assert list(gb)==list(gl)==[2*a+b,a*a]
    assert gb.reduce(a)[1]==a and gb.reduce(a*a)[1]==0
    # Reversed-source gauge from the actual binomial outputs, including height drop.
    P=x**3+y**3;R=3*(X**6+Y**6)-87*X**3*Y**3
    def reverse(phase):
        ans=1
        for m in range(1,4):
            for k in range(m):ans+=iterate(R,phase**m,k)*t**m*z**(m-k)/(S.factorial(m)*S.factorial(k))
        return S.expand(ans)
    F=reverse(P);Fp=reverse(P+u)
    gauge=sum((t*z*u)**j/S.factorial(j) for j in range(4))
    assert all(S.expand(Fp-gauge*F).coeff(t,j)==0 for j in range(4))
    # A supplied scalar time equation transports by D_t -> D_t-zu.
    exact=S.exp(z*x*t+z*t*t);exactp=S.exp(t*z*u)*exact
    assert S.simplify(S.diff(exactp,t)-z*(x+2*t+u)*exactp)==0
    # Coefficientwise Laurent descent still works with no uniform u-pole bound.
    scale=sum(z**j/u**(j+1) for j in range(5))
    recovered=S.expand(scale).coeff(u,-1)
    assert recovered==1
    assert S.expand(scale*(z*(x+2*t+u)-z*u)-scale*z*(x+2*t))==0
    # Quadratic radial finite-differential certificates for powers2 and3.
    radial=[]
    for m in [2,3]:
        seq=[];f=(Q+2)**m
        for k in range(1,m+1):
            f=S.expand(4*Q*S.diff(f,Q,2)+4*S.diff(f,Q));seq.append(f)
        matrix=S.Matrix([[f.coeff(Q,j) for f in seq] for j in range(m)])
        rhs=S.Matrix([S.expand((Q+2)**(m-1)).coeff(Q,j) for j in range(m)])
        weights=matrix.inv()*rhs
        assert S.expand(sum(c*f for c,f in zip(weights,seq))-(Q+2)**(m-1))==0
        direct=(x*x+y*y+2)**m
        for f in seq:
            direct=S.expand(S.diff(direct,x,2)+S.diff(direct,y,2))
            assert S.expand(direct-f.subs(Q,x*x+y*y))==0
        radial.append(dict(power=m,coefficients=list(map(str,weights))))
    # Polynomial shear and the precise commuting derivations.
    P=x**3+y**3;f=x*x*v*v+x*v+v**3
    dx=lambda g:S.expand(S.diff(g,x)+S.diff(P,x)*S.diff(g,v))
    dy=lambda g:S.expand(S.diff(g,y)+S.diff(P,y)*S.diff(g,v))
    assert S.expand(dx(f).subs(v,u+P)-S.diff(f.subs(v,u+P),x))==0
    assert S.expand(dy(f).subs(v,u+P)-S.diff(f.subs(v,u+P),y))==0
    assert S.expand(dx(dy(f))-dy(dx(f)))==0
    assert S.Poly(f.subs(v,u+P),x,y,u).total_degree()<=3*S.Poly(f,x,y,v).total_degree()
    # The theorem's zero-constant hypothesis is necessary.
    assert iterate(1,S.Integer(0),1)==0 and iterate(1,u,1)==u
    guards=0
    for lam,f in [(X,S.Float(1)),(X,x**25)]:
        try:action(lam,f)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failure')
    deps=['verify_formal_apolar_gaps.py','verify_scalar_time_gaps.py','verify_gvc_false_gaps.py','verify_phase_power_quotients.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in deps},
        nonnesting_countercontrols=counter,binomial_controls=binomials,
        ordinary_prefix_basis=list(map(str,gb)),ordinary_nilpotent_retained=True,
        reversed_source_gauge_coefficients=4,scalar_time_gauge_pass=True,
        laurent_descent_unbounded_pole_model_pass=True,quadratic_radial_certificates=radial,
        polynomial_shear_controls=3,zero_constant_boundary_pass=True,guard_refusals=guards,
        original_power_annihilators_always_nested=False,pure_safe_heights_always_monotone=False,
        nested_suspension_solves_general_gap=False,general_scalar_equation_produced=False,
        least_marked_cutoffs_identified=False,unrestricted_variable_coefficient_GVC_claimed=False,
        full_exceptional_MZ_decided=False,unrestricted_frontier_endpoint_changed=False,
        universal_proof_by_samples=False,independent_replay=False,external_review=False,
        formal_verification=False,novelty_claimed=False)
    assert guards==2
    (HERE/'linear_suspension_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,counterdegrees=[3,4],binomial_controls=binomials,ordinary_nilpotents=True,quadratic_certificates=2,shear_controls=3,guards=guards)))

if __name__=='__main__':main()
