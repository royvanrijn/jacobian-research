#!/usr/bin/env python3
"""Sufficient moment certificates without inverse; small exact cotangent controls."""
from pathlib import Path
import hashlib,json
import sympy as S
from verify_polynomial_gaussian_symmetries import Quotient,energy,compose

def sufficient_moment_certificate(phi,x,B,epsilon=None,length=1):
    x=tuple(x);ring=Quotient(epsilon,length)
    polyvars=x+((epsilon,) if epsilon is not None else ())
    if len(phi)!=len(x) or len(set(x))!=len(x) or not x:raise ValueError('matching polynomial coordinates required')
    phi=[S.Poly(f,*polyvars,domain=S.QQ_I).as_expr() for f in phi]
    if B.shape!=(len(x),len(x)) or B!=B.T or B.det()==0:raise ValueError('invertible symmetric covariance required')
    j=ring.red(S.Matrix(phi).jacobian(x).det());ell,jinv,u=ring.normalized_log(j,x)
    en=energy(x,B);res=ring.red(compose(en,phi,x,ring)-en-ell)
    constant=all(ring.red(S.diff(res,q))==0 for q in x)
    return dict(unit_Jacobian_verified=True,jacobian=str(j),normalized_log=str(ell),residual=str(res),
        sufficient_criterion_passed=constant,polynomial_inverse_verified=False,
        point_injectivity_inferred=False,necessity_for_endomorphisms_claimed=False)

def cotangent_data(f,x,p):
    f=S.Matrix([S.expand(v) for v in f]);A=f.jacobian(x);j=S.expand(A.det())
    if j==0 or any(S.diff(j,v)!=0 for v in x):raise ValueError('constant nonzero base Jacobian required')
    inv=A.adjugate().applyfunc(lambda z:S.expand(z/j));assert (A*inv).applyfunc(S.expand)==S.eye(len(x))
    y=(inv.T*S.Matrix(p)).applyfunc(S.expand)
    return f,A,inv,y

def scaling_cotangent_constructor(f,x,p,t,z):
    f,A,inv,y=cotangent_data(f,x,p);n=len(x)
    if f.subs(dict.fromkeys(x,0))!=S.zeros(n,1) or A.subs(dict.fromkeys(x,0))!=S.eye(n) or S.expand(A.det())!=1:raise ValueError('tangent-identity normalized Keller base required')
    ft=S.Matrix([S.expand(sum(c*t**(sum(a)-1)*S.prod(v**k for v,k in zip(x,a)) for a,c in S.Poly(v,*x).terms())) for v in f])
    ft,at,ati,yt=cotangent_data(ft,x,p)
    euler=(at*S.Matrix(x)-ft-t*ft.diff(t)).applyfunc(S.expand);assert euler==S.zeros(n,1)
    outz=S.expand(z+(ft.diff(t).T*yt)[0]);phi=list(ft)+list(yt)+[t,outz]
    defect=S.expand((ft.T*yt)[0]+t*outz-sum(a*b for a,b in zip(x,p))-t*z);assert defect==0
    return dict(map=phi,position_family=ft,inverse_Jacobian=ati,euler_residual=euler,energy_defect=defect,jacobian_by_verified_block_identity=1)

def wick_weighted(f,x,p,weights):
    out=0;n=len(x)
    for a,c in S.Poly(S.expand(f),*(tuple(x)+tuple(p))).terms():
        if tuple(a[:n])!=tuple(a[n:]):continue
        out+=c*S.prod(S.factorial(k)/S.Integer(w)**k for k,w in zip(a[:n],weights))
    return S.simplify(out)

def truncate_position_degree(f,x,p,d):
    variables=tuple(x)+tuple(p)
    return S.expand(sum(c*S.prod(v**k for v,k in zip(variables,a)) for a,c in S.Poly(f,*variables).terms() if sum(a[:len(x)])<=d))

def main():
    x,y,z,a,b,c,t,h,e=S.symbols('x y z a b c t h e');pos=(x,y,z);mom=(a,b,c);u=1+x*y
    f=[x-S.Rational(3,2)*x*x*y-x**3*z/2,
       y+3*x*u*u*z+3*x*y*y*(4+3*x*y),
       u**3*z+y*y*u*(4+3*x*y)]
    f,A,inv,Y=cotangent_data(f,pos,mom);assert S.expand(A.det())==1 and A.subs(dict.fromkeys(pos,0))==S.eye(3)
    weights=(1,-1,-2);W=S.diag(*weights)
    euler=(A*W*S.Matrix(pos)-W*f).applyfunc(S.expand);assert euler==S.zeros(3,1)
    en=(S.Matrix(pos).T*W*S.Matrix(mom))[0]
    defect=S.expand((f.T*W*Y)[0]-en);assert defect==0
    B=S.zeros(6)
    for i,w in enumerate(weights):B[i,3+i]=B[3+i,i]=1/S.Integer(w)
    assert S.expand(energy(pos+mom,B)-en)==0
    points=[(0,0,-S.Rational(1,4)),(1,-S.Rational(3,2),S.Rational(13,2)),(-1,S.Rational(3,2),S.Rational(13,2))]
    common=S.Matrix([0,0,-S.Rational(1,4)]);collision=[]
    for pt in points:
        assert f.subs(dict(zip(pos,pt)))==common
        yi=Y.subs(dict(zip(pos+mom,pt+(0,0,0))));assert yi==S.zeros(3,1)
        collision.append([str(v) for v in pt+(0,0,0)])
    # Exact low-degree moment controls, with safe pre-truncation by the number of momentum factors.
    # Polynomial multiplication never decreases position degree, and Wick requires equal totals.
    phi=list(f)+list(Y);monomials=[(0,)*6]
    for i in range(6):v=[0]*6;v[i]=1;monomials.append(tuple(v))
    for i in range(3):
        for j in range(3):v=[0]*6;v[i]=1;v[3+j]=1;monomials.append(tuple(v))
    for i in range(3):
        for k in [2,3]:v=[0]*6;v[3+i]=k;monomials.append(tuple(v))
    for beta in [(1,1,0),(2,0,1),(2,2,0)]:monomials.append((0,0,0)+beta)
    monomials.append((1,0,0,2,1,0))
    controls=[]
    for alpha in monomials:
        d=sum(alpha[3:]);small=[truncate_position_degree(v,pos,mom,d) for v in phi]
        source=S.expand(S.prod(v**k for v,k in zip(small,alpha)))
        actual=wick_weighted(source,pos,mom,weights)
        target=S.prod(v**k for v,k in zip(pos+mom,alpha));expected=wick_weighted(target,pos,mom,weights)
        assert S.simplify(actual-expected)==0
        controls.append(dict(target_exponents=list(alpha),value=str(actual),exact_position_truncation=d))
    # Fixed complex linear conjugation to ordinary independent real variance-one Gaussian inputs.
    real=S.symbols('U1 V1 U2 V2 U3 V3');split=[]
    for i in range(3):split.append(real[2*i]+S.I*real[2*i+1])
    for i,w in enumerate(weights):split.append((real[2*i]-S.I*real[2*i+1])/(2*w))
    C=S.Matrix(split).jacobian(real);assert (C*C.T-B).applyfunc(S.simplify)==S.zeros(6)
    assert S.expand(energy(split,B)-sum(v*v for v in real)/2)==0
    standard_collision=[[str(v) for v in C.inv()*S.Matrix([S.sympify(v) for v in pt])] for pt in collision]
    assert len({tuple(v) for v in standard_collision})==3
    # Universal constructor is checked on a small exact triangular base, without a degree campaign.
    xx=(x,y);pp=(a,b);toy=scaling_cotangent_constructor([x+y*y,y],xx,pp,t,h)
    v6=(x,y,a,b,t,h);B6=S.zeros(6)
    for i,j in [(0,2),(1,3),(4,5)]:B6[i,j]=B6[j,i]=1
    tc=sufficient_moment_certificate(toy['map'],v6,B6);assert tc['sufficient_criterion_passed']
    assert not tc['polynomial_inverse_verified'] and not tc['point_injectivity_inferred']
    # The nilpotent logarithmic correction remains sufficient without supplying an inverse.
    Z,Q=S.symbols('Z Q');B2=S.Matrix([[0,1],[1,0]])
    nil=sufficient_moment_certificate([Z-e*Z*Z,Q+e*(Z*Q-1)],(Z,Q),B2,e,2)
    assert nil['sufficient_criterion_passed']
    tr=sufficient_moment_certificate([Z+1,Q],(Z,Q),B2);assert not tr['sufficient_criterion_passed']
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        normalized_foundational_map=[str(v) for v in f],base_Jacobian_determinant='1',inverse_Jacobian=[[str(v) for v in inv.row(i)] for i in range(3)],
        weighted_Euler_residual=[str(v) for v in euler],weights=list(weights),lifted_energy_defect=str(defect),lifted_Jacobian_by_block_identity=1,
        six_dimensional_collision=collision,collision_image=[str(v) for v in list(common)+[0,0,0]],
        covariance=[[str(v) for v in B.row(i)] for i in range(6)],standard_gaussian_conjugation=[[str(v) for v in C.row(i)] for i in range(6)],
        standard_complex_collision=standard_collision,finite_moment_controls=controls,
        scaling_constructor_control=tc,nilpotent_sufficient_control=nil,translation_not_certified=tr,
        six_dimensional_map_max_degree=max(S.Poly(v,*(pos+mom)).total_degree() for v in phi),
        sufficient_logarithmic_endomorphism_certificate_implemented=True,weighted_cotangent_constructor_implemented=True,scaling_cotangent_constructor_implemented=True,
        general_inverse_solver_implemented=False,general_safe_pullback_theorem_claimed=False,necessary_endomorphism_log_criterion_claimed=False,
        all_source_marker_transfer_without_inverse_claimed=False,new_scalar_GMC_counterexample_constructed=False,universal_theorem_verified_by_controls=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    dest=Path(__file__).with_name('gaussian_keller_endomorphism_verification.json');dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=True,checker_sha256=receipt['checker_sha256'],finite_moment_controls=len(controls),six_dimensional_map_max_degree=receipt['six_dimensional_map_max_degree'],collision_size=len(collision),Gaussian_standard_conjugation_checked=True,universal_lift_control=tc['sufficient_criterion_passed'])))

if __name__=='__main__':main()
