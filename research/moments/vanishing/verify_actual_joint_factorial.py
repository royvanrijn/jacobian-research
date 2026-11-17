#!/usr/bin/env python3
"""Exact small controls for actual image/boundary and factorial recomposition proofs."""
from pathlib import Path
import json,hashlib,math
import sympy as S

def factorial_recomposition(h,z,m,kappa=1,c=1,max_index=24):
    h=S.sympify(h);c=S.sympify(c)
    if h.atoms(S.Float) or c.atoms(S.Float):raise ValueError('exact coefficients required')
    if h==0 or h.subs(z,0)!=0 or not isinstance(m,int) or m<1 or not isinstance(kappa,int) or kappa<1 or c==0:raise ValueError('nonzero zero-constant polynomial and positive integer indices required')
    p=S.Poly(h,z);d=p.degree()
    if d*m>max_index:raise ValueError('configured finite coefficient work limit;no classification')
    coeff=S.Poly(S.expand(h**m),z)
    value=S.expand(sum(a*c**k[0]*S.factorial(kappa*k[0]) for k,a in coeff.terms()))
    leading=p.LC()**m*c**(d*m)*S.factorial(kappa*d*m)
    next_=p.nth(d-1);limit=S.exp(next_/(c*d*p.LC())) if kappa==1 else S.S.One
    return dict(moment=value,normalized=S.simplify(value/leading),predicted_asymptotic_limit=limit,degree=d)

def circular_real_wick(f,z,w,t):
    result=0
    for (a,b,c),q in S.Poly(S.expand(f),z,w,t).terms():
        if a!=b or c%2:continue
        real=S.factorial(c)/(2**(c//2)*S.factorial(c//2))
        result+=q*S.factorial(a)*real
    return S.simplify(result)

def main():
    here=Path(__file__).parent;u,v,U,X=S.symbols('u v U X');delta=1-2*u*v;rel=delta*U*U-1
    def red(expr):
        num,den=S.together(expr).as_numer_denom()
        return S.rem(S.Poly(num,U,domain=S.QQ.frac_field(u,v)),S.Poly(rel,U,domain=S.QQ.frac_field(u,v))).as_expr()
    x=(U-1)/v
    du=lambda f:S.diff(f,u)+S.diff(f,U)*v*U**3
    dv=lambda f:S.diff(f,v)+S.diff(f,U)*u*U**3
    identities={
        'quadratic_inverse_relation':v*delta*x*x+2*delta*x-2*u,
        'du_X':du(x)-U**3,
        'dv_X':dv(x)-x*x*(U+2)/2,
        'constant_in_D_X':du(x)-2*v*v*dv(x)-3*v*x-1,
        'rational_generator_is_source_position_square':delta*(1+v*x)**2-1,
        'x_squared_generation_control':x*x-2*dv(x)+2*u*x+4*u*u/delta*(U-1)}
    remainders={name:str(red(expr)) for name,expr in identities.items()};assert set(remainders.values())=={'0'}
    odd=[];rational=[]
    for ell in [1,3,5]:
        assert red(U**(ell+2)-U**ell-S.Rational(2,ell)*u*du(U**ell))==0;odd.append(ell)
    for k in [1,2,3]:
        b=delta**(-k)
        assert S.simplify(delta**(-k-1)-b-u*S.diff(b,u)/k)==0;rational.append(k)
    xp=(1/S.sqrt(delta)-1)/v;xm=-2/v-xp;principal=[]
    for k in [1,2,3,4]:
        expansion=S.series(S.diff(xm,v,k-1),v,0,1).removeO().expand()
        negative=S.Add(*(q for q in S.Add.make_args(expansion) if q.as_powers_dict().get(v,0)<0))
        expected=-2*(-1)**(k-1)*S.factorial(k-1)/v**k;assert S.simplify(negative-expected)==0
        plus=S.series(xp**k,v,0,1).removeO();assert not any(q.as_powers_dict().get(v,0)<0 for q in S.Add.make_args(plus))
        principal.append(dict(derivative_order=k-1,minus_principal_part=str(negative)))
    # Discriminant and the first-order square-root connection are exact identities.
    quadratic=v*delta*X*X+2*delta*X-2*u
    assert S.expand(S.discriminant(quadratic,X)-4*delta)==0
    assert S.simplify(v/delta+S.diff(delta,u)/(2*delta))==0
    assert S.simplify(u/delta+S.diff(delta,v)/(2*delta))==0
    z=S.Symbol('z');examples=[]
    for h,kappa,c in [(z+z*z,1,1),(z+z*z,2,S.Rational(1,2)),(z**3-z*z+2*z,1,1),(z*z+S.I*z,1,1),(z**3,3,2)]:
        rows=[]
        for m in [1,2,3]:
            q=factorial_recomposition(h,z,m,kappa,c)
            # Direct coefficient convolution is also an exact independently ordered sum.
            p=S.Poly(h,z);coeff={0:S.S.One}
            for j in range(m):
                nxt={}
                for a,b in coeff.items():
                    for (power,),value in p.terms():nxt[a+power]=S.expand(nxt.get(a+power,0)+b*value)
                coeff=nxt
            other=S.expand(sum(value*S.sympify(c)**index*S.factorial(kappa*index) for index,value in coeff.items()))
            assert S.expand(other-q['moment'])==0
            rows.append(dict(m=m,moment=str(q['moment']),normalized=str(q['normalized'])))
        examples.append(dict(h=str(h),kappa=kappa,c=str(c),predicted_limit=str(q['predicted_asymptotic_limit']),finite_controls=rows))
    for m in [1,2,3,4]:
        q=factorial_recomposition(z+z*z,z,m,2,S.Rational(1,2))
        bound=(1+S.Rational(2,m*m))**m-1
        assert 0<=q['normalized']-1<=bound
    Z,W,T=S.symbols('Z W T');P=(1+Z)*(W-S.Rational(1,2)*(2+Z)*T*T);base=[]
    for k in [1,2,3,4]:
        pure=circular_real_wick(P**k,Z,W,T);marked=circular_real_wick(Z*P**k,Z,W,T)
        assert pure==0 and marked==S.factorial(k);base.append(dict(k=k,pure=str(pure),marked=str(marked)))
    direct=[]
    hP=P+P*P
    for m in [1,2]:
        value=circular_real_wick(Z*hP**m,Z,W,T);q=factorial_recomposition(z+z*z,z,m)
        assert value==q['moment'];direct.append(dict(m=m,marked=str(value)))
    refusals=0
    for h,m in [(1+z,1),(S.Float('1.0')*z,1),(z*z,13)]:
        try:factorial_recomposition(h,z,m)
        except ValueError:refusals+=1
    assert refusals==3
    # Monomial annihilator controls;the union/normalization formula is proved separately.
    quotient_relation=(1,1);normalization=[]
    for N in [1,2,3]:
        colon=tuple(max(0,a-b) for a,b in zip(quotient_relation,(N,0)));assert colon==(0,1)
        normalization.append(dict(coefficient_ring='K[a,b]/(ab)',I_power=N,annihilator_generator='b'))
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        quadratic_field_identity_remainders=remainders,odd_generator_recursion_controls=odd,rational_generator_recursion_controls=rational,
        minus_boundary_principal_parts=principal,joint_rank=2,actual_D_generating_source_markers=['x','(1+x*y)**2'],
        defect_filtration='C/O=(B/O) direct_sum B*U;R/C=H_v;written nonsplit extension proof',
        factorial_recomposition_controls=examples,Long_finite_base_controls=base,Long_direct_recomposition_controls=direct,guard_refusals=refusals,
        normalization_annihilator_controls=normalization,
        exact_recomposition_coefficients_implemented=True,quadratic_field_identity_checker_implemented=True,
        supporting_checker_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ['verify_joint_safe_algebras.py','verify_nilpotent_inverse_sections.py','verify_weighted_cotangent_inverse.py','verify_gaussian_keller_endomorphisms.py','verify_polynomial_gaussian_symmetries.py']},
        universal_theorem_verified_by_finite_controls=False,numerical_asymptotics_used_as_proof=False,general_Weyl_presentation_or_kernel_solver_implemented=False,
        complete_individual_nonlinear_safety_classifier_implemented=False,new_Gaussian_moment_preserving_map_constructed=False,minimum_GMC_dimension_or_degree_improved=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (here/'actual_joint_factorial_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','joint_rank','quadratic_field_identity_remainders','Long_direct_recomposition_controls','actual_D_generating_source_markers']}))

if __name__=='__main__':main()
