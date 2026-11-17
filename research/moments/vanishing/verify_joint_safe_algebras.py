#!/usr/bin/env python3
"""Finite exact controls for joint inverse gates and written nonclosure proofs."""
from pathlib import Path
import json,hashlib
import sympy as S
from verify_weighted_cotangent_inverse import occupations,multiply_to_cap

def marked_linear_momentum(vector,m,gamma,beta,x,covariance):
    n=len(x);out=0;terms=[dict(S.Poly(v,*x).terms()) for v in vector]
    for alpha in occupations(m,n):
        cap=tuple(alpha[i]+beta[i]-gamma[i] for i in range(n))
        if min(cap)<0:continue
        prod={(0,)*n:S.S.One}
        for source,k in zip(terms,alpha):
            small={a:c for a,c in source.items() if all(i<=j for i,j in zip(a,cap))}
            for j in range(k):prod=multiply_to_cap(prod,small,cap)
        out+=prod.get(cap,0)*S.factorial(m)*S.prod(S.factorial(alpha[i]+beta[i])/S.factorial(alpha[i])*covariance[i]**(alpha[i]+beta[i]) for i in range(n))
    return S.expand(out)

def joint_moment(v,w,i,j,gamma,beta,x,covariance):
    a,b=S.symbols('_joint_a _joint_b')
    value=marked_linear_momentum(a*v+b*w,i+j,gamma,beta,x,covariance)
    return S.simplify(value.coeff(a,i).coeff(b,j)/S.binomial(i+j,i))

def joint_tail_generators(f,x,V,t,N,max_order=8,max_parameters=3):
    F=S.Matrix(f);V=S.Matrix(V);n=len(x)
    if not isinstance(N,int) or N<1 or N>max_order or len(t)>max_parameters:raise ValueError('joint gate exceeds configured order/parameter limit;no classification')
    if V.shape!=(n,len(t)):raise ValueError('direction matrix shape mismatch')
    if any(S.sympify(q).atoms(S.Float) for q in list(F)+list(V)):raise ValueError('exact coefficients required')
    A=F.jacobian(x)
    if F.subs(dict.fromkeys(x,0))!=S.zeros(n,1) or A.subs(dict.fromkeys(x,0))!=S.eye(n) or S.expand(A.det())!=1:raise ValueError('normalized Keller identities required')
    trunc=lambda q:S.Add(*(c*S.prod(u**k for u,k in zip(t,a)) for a,c in S.Poly(S.expand(q),*t).terms() if sum(a)<=N))
    b=S.zeros(n,1);rhs=V*S.Matrix(t);nonlinear=F-S.Matrix(x)
    for k in range(N):b=(rhs-nonlinear.subs(dict(zip(x,b)),simultaneous=True)).applyfunc(trunc)
    residual=(F.subs(dict(zip(x,b)),simultaneous=True)-rhs).applyfunc(S.expand)
    generators=[]
    for value in residual:
        for a,c in S.Poly(value,*t).terms():
            if c!=0 and c not in generators:generators.append(c)
    return dict(branch=list(b),ordinary_generators=generators,residual=list(residual),total_order=N)

def main():
    here=Path(__file__).parent;x,y,z,u,v,U=S.symbols('x y z u v U');pos=(x,y,z)
    old=json.loads((here/'gaussian_keller_endomorphism_verification.json').read_text())
    F=S.Matrix([S.sympify(q) for q in old['normalized_foundational_map']]);Ainv=S.Matrix([[S.sympify(q) for q in row] for row in old['inverse_Jacobian']])
    W=[1,-1,-2];cov=[1,-1,-S.Rational(1,2)];vs=Ainv[:,0];vt=-Ainv[:,1]
    assert (F.jacobian(pos)*Ainv).applyfunc(S.expand)==S.eye(3)
    for vec,charge in [(vs,-1),(vt,1)]:
        for j,c in enumerate(vec):assert S.expand(sum(w*q*S.diff(c,q) for q,w in zip(pos,W))-W[j]*c-charge*c)==0
    degrees=[1+max(int(S.Poly(c,*pos).total_degree()) for c in vec) for vec in [vs,vt]];assert degrees==[12,10]
    relation=(1-2*u*v)*U*U-1;branch=[(U-1)/v,v,-v*v*(1+3*U)/(U*U)];remainders=[]
    for q,target in zip(F,[u,v,0]):
        value=S.together(q.subs(dict(zip(pos,branch)),simultaneous=True)-target);num,den=value.as_numer_denom()
        rem=S.rem(S.Poly(num,U,domain=S.QQ.frac_field(u,v)),S.Poly(relation,U,domain=S.QQ.frac_field(u,v))).as_expr()
        assert rem==0 and den!=0;remainders.append(str(rem))
    zero=(0,0,0);sharp=[]
    for label,vec,gamma,beta,expected in [('S',vs,zero,(0,0,1),-15),('T',vt,(0,0,1),zero,-8)]:
        actual=marked_linear_momentum(vec,2,gamma,beta,pos,cov);assert actual==expected
        sharp.append(dict(phase=label,marker='p3' if label=='S' else 'z',index=2,value=str(actual)))
    pure=[]
    for i,j in [(1,1),(2,1),(1,2)]:
        actual=joint_moment(vs,vt,i,j,zero,zero,pos,cov);assert actual==0
        pure.append(dict(powers=[i,j],value=str(actual)))
    product=[]
    for m in [1,2]:
        actual=joint_moment(vs,vt,m,m,(1,1,0),zero,pos,cov);expected=S.factorial(2*m)/2**m;assert actual==expected
        product.append(dict(index=m,marker='xy',value=str(actual)))
    # Formal binomial extraction agrees without expanding the degree22 product.
    coeffs=[]
    for m in [1,2,3]:
        c=S.binomial(2*m,m)/2**m;assert S.factorial(m)**2*c==S.factorial(2*m)/2**m
        coeffs.append(str(c))
    a,b,c,d,t1,t2=S.symbols('a b c d t1 t2');V=S.Matrix([[a,b],[c,d]]);t=(t1,t2)
    controls=[]
    for map_,N,expected in [([x+y*y,y],1,[c*c,c*d,d*d]),([x+y*y,y],2,[]),([x+y*y+y**3,y],2,[c**3,c*c*d,c*d*d,d**3]),([x+y*y+y**3,y],3,[])]:
        q=joint_tail_generators(map_,(x,y),V,t,N)
        basis=[s.as_expr() for s in S.groebner(q['ordinary_generators'],a,b,c,d).polys] if q['ordinary_generators'] else []
        assert basis==expected
        controls.append(dict(map=[str(v) for v in map_],order=N,branch=[str(v) for v in q['branch']],ordinary_ideal_basis=[str(v) for v in basis]))
    change=S.Matrix([[1,1],[0,1]])
    q=joint_tail_generators([x+y*y,y],(x,y),V*change,t,1)
    changed_basis=[v.as_expr() for v in S.groebner(q['ordinary_generators'],a,b,c,d).polys];assert changed_basis==[c*c,c*d,d*d]
    refusals=0
    for N in [0,9]:
        try:joint_tail_generators([x+y*y,y],(x,y),V,t,N)
        except ValueError:refusals+=1
    try:joint_tail_generators([x+S.Float('1.0')*y*y,y],(x,y),V,t,1)
    except ValueError:refusals+=1
    assert refusals==3
    # Tensor grid interpolation is tested only on a two-parameter quadratic toy.
    poly=a*t1*t1+b*t1*t2+c*t2*t2+d
    grid=[0,1,2];basis1=[S.prod((t1-j)/(i-j) for j in grid if j!=i) for i in grid];basis2=[q.subs(t1,t2) for q in basis1]
    interp=sum(poly.subs({t1:i,t2:j})*basis1[i]*basis2[j] for i in grid for j in grid);assert S.expand(interp-poly)==0
    sig,tau=S.symbols('sigma tau');h=sig*(1+sig*tau);assert S.expand(-sig*S.diff(h,sig)+tau*S.diff(h,tau)+h)==0
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        phase_degrees=degrees,product_phase_degree=sum(degrees),Euler_charges=[-1,1],joint_branch=[str(q) for q in branch],joint_branch_relation=str(relation),branch_remainders=remainders,
        sharp_individual_slope_controls=sharp,finite_joint_pure_controls=pure,finite_product_witness_controls=product,binomial_coefficient_controls=coeffs,
        ordinary_joint_gate_controls=controls,frame_change_preserves_ordinary_ideal=True,guard_refusals=refusals,toy_tensor_interpolation_passed=True,
        safe_mixed_nonlinear_example='S*(1+S*T),Euler charge-1,cutoff2',exact_safe_plane='lambda*S+mu*T safe iff lambda*mu nilpotent;written proof',
        exact_product_scaling_ideals='I_N=(e**(2*N+1)) for every integerN>=0;written proof',sharp_product_nilpotent_cutoff='real(L-1)/2;integer floor(L/2);written proof',
        guarded_multivariate_tail_generator_implemented=True,sparse_momentum_marker_extractor_implemented=True,
        large_joint_grid_enumerated=False,foundational_order49_branch_expanded=False,
        supporting_checker_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ['verify_weighted_cotangent_inverse.py','verify_nilpotent_inverse_sections.py','verify_gaussian_keller_endomorphisms.py','verify_polynomial_gaussian_symmetries.py']},
        universal_theorem_verified_by_finite_controls=False,all_index_numerical_sweep_performed=False,general_inverse_solver_implemented=False,
        complete_individual_nonlinear_safety_classifier_implemented=False,new_Gaussian_moment_preserving_map_constructed=False,minimum_GMC_dimension_or_degree_improved=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (here/'joint_safe_algebras_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','phase_degrees','product_phase_degree','branch_remainders','sharp_individual_slope_controls','finite_product_witness_controls']}))

if __name__=='__main__':main()
