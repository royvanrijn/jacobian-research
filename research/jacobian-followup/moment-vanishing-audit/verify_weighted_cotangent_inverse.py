#!/usr/bin/env python3
"""Exact inverse-branch and sparse Wick certificates, not an all-index numerical proof."""
from pathlib import Path
import hashlib,json
import sympy as S
from verify_gaussian_keller_endomorphisms import cotangent_data,scaling_cotangent_constructor

def occupations(total,n):
    if n==1:
        yield (total,);return
    for i in range(total+1):
        for tail in occupations(total-i,n-1):yield (i,)+tail

def multiply_to_cap(left,right,cap):
    out={}
    for a,c in left.items():
        for b,d in right.items():
            e=tuple(i+j for i,j in zip(a,b))
            if all(i<=j for i,j in zip(e,cap)):out[e]=out.get(e,0)+c*d
    return {a:S.expand(c) for a,c in out.items() if c!=0}

def linear_momentum_moment(vector,m,marker_exponent,x,covariance_diagonal):
    # P=sum vector_i(x)*p_i, Q=x^gamma.
    # Wick forces the occupation alpha in p to match x; componentwise truncation is exact.
    n=len(x);result=0
    terms=[dict(S.Poly(v,*x).terms()) for v in vector]
    for alpha in occupations(m,n):
        if any(i<j for i,j in zip(alpha,marker_exponent)):continue
        cap=tuple(i-j for i,j in zip(alpha,marker_exponent));prod={(0,)*n:S.S.One}
        for source,k in zip(terms,alpha):
            small={a:c for a,c in source.items() if all(i<=j for i,j in zip(a,cap))}
            for j in range(k):prod=multiply_to_cap(prod,small,cap)
        result+=prod.get(cap,0)*S.prod(d**k for d,k in zip(covariance_diagonal,alpha))
    return S.simplify(S.factorial(m)*result)

def generic_inverse_phase(f,x,y,t,z,direction):
    lift=scaling_cotangent_constructor(f,x,y,t,z);n=len(x)
    phi=lift['map'];Y=S.Matrix(phi[n:2*n]);ft=lift['position_family']
    phase=S.expand(z+((S.Matrix(direction)+ft.diff(t)).T*Y)[0])
    assert S.expand(phase-(phi[-1]+(S.Matrix(direction).T*Y)[0]))==0
    return phase,lift

def split_wick(f,x,p):
    out=0;n=len(x)
    for alpha,c in S.Poly(S.expand(f),*(tuple(x)+tuple(p))).terms():
        if alpha[:n]!=alpha[n:]:continue
        out+=c*S.prod(S.factorial(k) for k in alpha[:n])
    return S.simplify(out)

def inverse_line_certificate(f,x,direction,max_branch_order=16):
    n=len(x);F=S.Matrix(f);A=F.jacobian(x)
    if any(S.sympify(v).atoms(S.Float) for v in list(F)+list(direction)):raise ValueError('exact coefficients required')
    if F.subs(dict.fromkeys(x,0))!=S.zeros(n,1) or A.subs(dict.fromkeys(x,0))!=S.eye(n) or S.expand(A.det())!=1:raise ValueError('normalized Keller identities required')
    degree=max(1,max(int(S.Poly(v,*x).total_degree()) for v in F));bound=degree**(n-1)
    if bound>max_branch_order:raise ValueError('Bezout branch order exceeds configured work limit;no classification performed')
    s=S.Symbol('_branch_time');b=S.zeros(n,1);rhs=s*S.Matrix(direction);nonlinear=F-S.Matrix(x)
    trunc=lambda v:S.expand(S.series(S.expand(v),s,0,bound+1).removeO())
    for k in range(bound):b=(rhs-nonlinear.subs(dict(zip(x,b)),simultaneous=True)).applyfunc(trunc)
    residual=(F.subs(dict(zip(x,b)),simultaneous=True)-rhs).applyfunc(S.expand)
    return dict(polynomial_inverse_section=residual==S.zeros(n,1),classification_scope='characteristic-zero field directions',
        branch_degree_bound=bound,base_degree=degree,branch=[str(v) for v in b],composition_residual=[str(v) for v in residual],full_all_ring_converse_claimed=False)

def safe_direction_certificate(f,x,weights,direction,max_branch_order=16):
    n=len(x);F=S.Matrix(f);W=S.diag(*weights)
    if (F.jacobian(x)*W*S.Matrix(x)-W*F).applyfunc(S.expand)!=S.zeros(n,1):raise ValueError('weighted Euler identity required')
    result=inverse_line_certificate(f,x,W.inv()*S.Matrix(direction),max_branch_order)
    degree=result['base_degree'];bound=result['branch_degree_bound'];hbound=(n-1)*(degree-1)
    cp=(2*n-1)*(1+bound*max(hbound-1,0))
    result.update(certified_safe=result['polynomial_inverse_section'],sufficient_common_linear_cutoff=max(bound,cp))
    return result

def main():
    here=Path(__file__).parent;prior=json.loads((here/'gaussian_keller_endomorphism_verification.json').read_text())
    x,y,z,p1,p2,p3,s,U=S.symbols('x y z p1 p2 p3 s U');pos=(x,y,z);mom=(p1,p2,p3)
    f=[S.sympify(v) for v in prior['normalized_foundational_map']]
    f,A,inv,Y=cotangent_data(f,pos,mom);W=S.diag(1,-1,-2);D=W.inv();a=S.Matrix([1,-1,0]);vector=(inv*a).applyfunc(S.expand)
    assert inv==S.Matrix([[S.sympify(v) for v in row] for row in prior['inverse_Jacobian']])
    euler=(A*W*S.Matrix(pos)-W*f).applyfunc(S.expand);assert euler==S.zeros(3,1)
    # Exact determinant-cancellation matrix identity before substitution into a formal branch.
    target=S.Matrix(S.symbols('v1 v2 v3'));h=D*inv*target
    left=(S.eye(3)-h.jacobian(pos)).subs(dict(zip(target,W*f)),simultaneous=True).applyfunc(S.expand)
    right=(D*inv*W*A).applyfunc(S.expand);assert left==right
    determinant_product=S.simplify(D.det()*inv.det()*W.det()*A.det());assert determinant_product==1
    # Exact inverse branch in Q(s)[U]/((1-2s^2)U^2-1), with U(0)=1.
    relation=(1-2*s*s)*U*U-1;branch=[(U-1)/s,s,-s*s*(1+3*U)/(U*U)]
    reductions=[]
    for value,expected in zip(f,[s,s,0]):
        value=S.together(value.subs(dict(zip(pos,branch)),simultaneous=True)-expected)
        num,den=value.as_numer_denom();rem=S.rem(S.Poly(num,U,domain=S.QQ.frac_field(s)),S.Poly(relation,U,domain=S.QQ.frac_field(s))).as_expr()
        assert rem==0 and den!=0;reductions.append(str(rem))
    use=1/S.sqrt(1-2*s*s);egf=(use-1)/s
    pure=[];marked=[];amplitude_controls=[]
    covariance=[1,-1,-S.Rational(1,2)]
    for m in range(1,6):
        value=linear_momentum_moment(vector,m,(0,0,0),pos,covariance);assert value==0;pure.append(str(value))
        value=linear_momentum_moment(vector,m,(1,0,0),pos,covariance)
        expected=S.factorial(m)*S.expand(S.series(egf,s,0,m+1).removeO()).coeff(s,m)
        assert value==expected;marked.append(str(value))
    assert marked==['1','0','9','0','300']
    for gamma,branchvalue in [((0,1,0),s),((0,0,1),-s*s*(1+3*use)/(use*use)),((1,1,0),s*egf)]:
        for m in [1,2,3,4]:
            value=linear_momentum_moment(vector,m,gamma,pos,covariance)
            expected=S.factorial(m)*S.expand(S.series(branchvalue,s,0,m+1).removeO()).coeff(s,m)
            assert S.simplify(value-expected)==0;amplitude_controls.append(str(value))
    wrong_transpose=linear_momentum_moment(inv.T*a,2,(0,0,0),pos,covariance);assert wrong_transpose==-3
    wrong_sign=linear_momentum_moment(inv*S.Matrix([1,1,0]),3,(1,0,0),pos,covariance);assert wrong_sign==-9
    # Standard Gaussian target linear form is isotropic and source marker stays linear.
    real=S.symbols('U1 V1 U2 V2 U3 V3');C=S.Matrix([[S.sympify(v) for v in row] for row in prior['standard_gaussian_conjugation']])
    target_linear=(real[0]-S.I*real[1]+real[2]-S.I*real[3])/2
    coeff=S.Matrix([S.diff(target_linear,v) for v in real]);assert S.simplify((coeff.T*coeff)[0])==0
    assert S.expand((S.Matrix([0,0,0,1,-1,0]).T*C*S.Matrix(real))[0]-target_linear)==0
    # A small general normalized-Keller phase construction checks t=s and the signs in (12).
    t,hv,b1,b2=S.symbols('t h b1 b2');aa=(S.Integer(2),S.Integer(3))
    phase,lift=generic_inverse_phase([x+y*y,y],(x,y),(b1,b2),t,hv,aa)
    assert S.expand(phase-(hv+2*b1+3*b2+y*y*b1-6*t*y*b1))==0
    generic_controls=[]
    for marker,series in [(x,2*s-9*s**3),(y,3*s),(t,s)]:
        for m in [1,2,3]:
            value=split_wick(marker*phase**m,(x,y,t),(b1,b2,hv));expected=S.factorial(m)*S.expand(series).coeff(s,m)
            assert value==expected;generic_controls.append(str(value))
    safe_gate=safe_direction_certificate([x+y*y,y],(x,y),(2,1),(2,3))
    assert safe_gate['certified_safe'] and safe_gate['branch_degree_bound']==2
    assert safe_gate['composition_residual']==['0','0']
    general_gate=inverse_line_certificate([x+y*y+y**3,y],(x,y),(2,3));assert general_gate['polynomial_inverse_section'] and general_gate['branch_degree_bound']==3
    zero_gate=safe_direction_certificate([x+y*y,y],(x,y),(2,1),(0,0));assert zero_gate['certified_safe']
    limit_rejected=False
    try:safe_direction_certificate(list(f),pos,(1,-1,-2),(1,-1,0),max_branch_order=16)
    except ValueError as error:limit_rejected='work limit' in str(error)
    assert limit_rejected
    phase_expr=S.expand(sum(v*p for v,p in zip(vector,mom)));degree=int(S.Poly(phase_expr,*(pos+mom)).total_degree());assert degree==12
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        base_map=[str(v) for v in f],direction=['1','-1','0'],inverse_covariance_weights=['1','-1','-2'],
        phase_coefficient_vector=[str(S.factor(v)) for v in vector],phase_total_degree=degree,source_marker='x',
        target_standard_linear_phase=str(target_linear),source_standard_linear_marker='U1 + I*V1',standard_conjugation=prior['standard_gaussian_conjugation'],
        determinant_cancellation_matrix_checked=True,determinant_product=str(determinant_product),formal_branch_relation=str(relation),formal_branch=[str(v) for v in branch],branch_remainders=reductions,
        marked_egf='((1-2*s**2)**(-1/2)-1)/s',all_odd_index_formula='(2*k-1)!*binomial(2*k,k)/2**k for k>=1',
        finite_pure_controls=pure,finite_marked_controls=marked,finite_other_amplitude_controls=amplitude_controls,
        wrong_transpose_second_pure=str(wrong_transpose),wrong_covariance_sign_third_marked=str(wrong_sign),
        universal_phase_constructor_control=str(phase),universal_phase_small_controls=generic_controls,
        supporting_checker_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ['verify_gaussian_keller_endomorphisms.py','verify_polynomial_gaussian_symmetries.py']},
        guarded_safe_direction_gate_implemented=True,safe_direction_gate_controls=[safe_gate,zero_gate],general_inverse_line_control=general_gate,large_branch_gate_work_limit_enforced=limit_rejected,
        sparse_exact_moment_extractor_implemented=True,generic_inverse_phase_constructor_implemented=True,
        explicit_linear_failure_witness_constructed=True,forward_safety_implication_disproved_internally=True,
        universal_theorem_verified_by_finite_controls=False,all_index_numerical_sweep_performed=False,
        general_inverse_solver_implemented=False,complete_surviving_phase_classifier_implemented=False,
        new_Gaussian_moment_preserving_map_constructed=False,minimum_GMC_dimension_or_degree_improved=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    out=here/'weighted_cotangent_inverse_verification.json';out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','phase_total_degree','finite_marked_controls','branch_remainders','wrong_transpose_second_pure','wrong_covariance_sign_third_marked']}))

if __name__=='__main__':main()
