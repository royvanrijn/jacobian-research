#!/usr/bin/env python3
"""Finite exact controls for written inverse-section and all-index ideal proofs."""
from pathlib import Path
import hashlib,json
import sympy as S
from verify_weighted_cotangent_inverse import inverse_line_certificate,linear_momentum_moment

def cutoff_bound(n,degree,nilpotence_order):
    if any(not isinstance(q,int) or q<1 for q in [n,degree,nilpotence_order]):raise ValueError('positive integer bounds required')
    B=degree**(n-1);H=(n-1)*(degree-1);M=H+degree
    k=(nilpotence_order-1).bit_length();E=B*M**k
    weighted=max(E,(2*n-1)*(1+E*max(H-1,0)))
    suspended_E=max(1,2*E-1);suspended_H=n*(2*degree-2)
    suspended=max(suspended_E,(2*n+1)*(1+suspended_E*max(suspended_H-1,0)))
    return dict(field_branch_bound=B,inverse_Jacobian_degree_bound=H,Newton_steps=k,corrected_branch_degree_bound=E,weighted_sufficient_slope=weighted,unweighted_suspension_sufficient_slope=suspended)

def tail_generators(f,x,direction,N,max_branch_order=16):
    if not isinstance(N,int) or N<1 or N>max_branch_order:raise ValueError('configured positive branch order required;no classification on refusal')
    F=S.Matrix(f);A=F.jacobian(x);n=len(x)
    if any(S.sympify(v).atoms(S.Float) for v in list(F)+list(direction)):raise ValueError('exact inputs required')
    if F.subs(dict.fromkeys(x,0))!=S.zeros(n,1) or A.subs(dict.fromkeys(x,0))!=S.eye(n) or S.expand(A.det())!=1:raise ValueError('normalized Keller identities required')
    s=S.Symbol('_section_time');b=S.zeros(n,1);rhs=s*S.Matrix(direction);nonlinear=F-S.Matrix(x)
    trunc=lambda q:S.expand(S.series(q,s,0,N+1).removeO())
    for k in range(N):b=(rhs-nonlinear.subs(dict(zip(x,b)),simultaneous=True)).applyfunc(trunc)
    residual=(F.subs(dict(zip(x,b)),simultaneous=True)-rhs).applyfunc(S.expand)
    generators=[]
    for v in residual:
        for c in S.Poly(v,s).all_coeffs():
            if c!=0 and c not in generators:generators.append(S.expand(c))
    return dict(branch=list(b),residual=list(residual),ordinary_generators=generators,branch_order=N)

def reduce_univariate(f,u,order):
    return S.Add(*(c*u**a[0] for a,c in S.Poly(S.expand(f),u).terms() if a[0]<order))

def compose_quotient(f,x,b,u,order):
    # Truncate after every multiplication; never expand a large unreduced substitution.
    out=S.S.Zero
    for powers,coefficient in S.Poly(f,*x).terms():
        term=coefficient
        for value,power in zip(b,powers):
            for j in range(power):term=reduce_univariate(term*value,u,order)
        out=reduce_univariate(out+term,u,order)
    return out

def nilpotent_newton(f,x,direction,seed,order,max_nilpotence_order=8,wrong_transpose=False):
    if not isinstance(order,int) or order<1 or order>max_nilpotence_order:raise ValueError('monogenic quotient order exceeds configured work limit')
    n=len(x);F=S.Matrix(f);A=F.jacobian(x)
    if F.subs(dict.fromkeys(x,0))!=S.zeros(n,1) or A.subs(dict.fromkeys(x,0))!=S.eye(n) or S.expand(A.det())!=1:raise ValueError('normalized Keller identities required')
    u=S.Symbol('_nil_time');inv=A.adjugate();inv=inv.T if wrong_transpose else inv
    b=S.Matrix([S.sympify(q).subs(S.Symbol('u'),u) for q in seed]);rhs=u*S.Matrix(direction)
    red=lambda q:reduce_univariate(q,u,order)
    compose=lambda q:compose_quotient(q,x,b,u,order)
    rows=[]
    for j in range((order-1).bit_length()+1):
        residual=(F.applyfunc(compose)-rhs).applyfunc(red)
        valuation=min([a[0] for value in residual for a,c in S.Poly(value,u).terms() if c!=0],default=order)
        rows.append(dict(step=j,error_order_at_least=valuation,residual=[str(v) for v in residual]))
        if residual==S.zeros(n,1):return dict(exact_section=True,quotient_order=order,branch=[str(v) for v in b],errors=rows)
        correction=(inv.applyfunc(compose)*residual).applyfunc(red);b=(b-correction).applyfunc(red)
    return dict(exact_section=False,quotient_order=order,branch=[str(v) for v in b],errors=rows)

def main():
    here=Path(__file__).parent;x,y,z,a,b,u,s=S.symbols('x y z a b u s');pos=(x,y,z)
    controls=[]
    for F,N,expected in [([x+y*y,y],1,[b*b]),([x+y*y,y],2,[]),([x+y*y+y**3,y],1,[b*b]),([x+y*y+y**3,y],2,[b**3]),([x+y*y+y**3,y],3,[])]:
        q=tail_generators(F,(x,y),(a,b),N)
        basis=list(S.groebner(q['ordinary_generators'],a,b).polys) if q['ordinary_generators'] else []
        assert [v.as_expr() for v in basis]==expected
        controls.append(dict(map=[str(v) for v in F],order=N,branch=[str(v) for v in q['branch']],ordinary_ideal_basis=[str(v.as_expr()) for v in basis]))
    # The ordinary first gate is (b^2),not its radical (b).
    assert controls[0]['ordinary_ideal_basis']==['b**2']
    rejected=[]
    for N in [0,17]:
        try:tail_generators([x+y*y,y],(x,y),(a,b),N)
        except ValueError:rejected.append(N)
    assert rejected==[0,17]
    prior=json.loads((here/'weighted_cotangent_inverse_verification.json').read_text())
    old=json.loads((here/'gaussian_keller_endomorphism_verification.json').read_text())
    F=S.Matrix([S.sympify(v) for v in prior['base_map']]);inv=S.Matrix([[S.sympify(v) for v in row] for row in old['inverse_Jacobian']])
    assert (F.jacobian(pos)*inv).applyfunc(S.expand)==S.eye(3)
    lifts=[]
    for order in [3,5]:
        q=nilpotent_newton(F,pos,(1,1,0),(u,u,0),order);assert q['exact_section']
        U=1/S.sqrt(1-2*s*s)
        expected=[(U-1)/s,s,-s*s*(1+3*U)/(U*U)]
        uu=S.Symbol('_nil_time')
        for value,expr in zip(q['branch'],expected):
            target=S.expand(S.series(expr,s,0,order).removeO()).subs(s,uu)
            assert S.expand(S.sympify(value)-target)==0
        valuations=[row['error_order_at_least'] for row in q['errors']]
        assert valuations[0]==2 and valuations[-1]==order
        assert all(right>=min(2*left,order) for left,right in zip(valuations,valuations[1:]))
        lifts.append(q)
    wrong=nilpotent_newton(F,pos,(1,1,0),(u,u,0),5,wrong_transpose=True)
    assert wrong['errors'][1]['error_order_at_least']==3 < 4
    refusal=False
    try:nilpotent_newton(F,pos,(1,1,0),(u,u,0),9)
    except ValueError:refusal=True
    assert refusal
    # Independent sparse Wick values check both parity witnesses at the retained small orders.
    vector=inv*S.Matrix([1,-1,0]);cov=[1,-1,-S.Rational(1,2)];moments=[]
    for m,gamma,expected in [(1,(1,0,0),1),(2,(0,0,1),-8),(3,(1,0,0),9),(4,(0,0,1),120)]:
        value=linear_momentum_moment(vector,m,gamma,pos,cov);assert value==expected
        moments.append(dict(moment_index=m,marker='x' if m%2 else 'z',value=str(value),scaling_factor='e**'+str(m)))
    # Binomial coefficients are controls only;their sign formula is proved in the note.
    coeffs=[]
    for k in range(1,5):
        c=-4 if k==1 else 5 if k==2 else -3*S.binomial(S.Rational(1,2),k-1)*(-2)**(k-1)
        expr=-s*s*(1-2*s*s)-3*s*s*S.sqrt(1-2*s*s)
        assert S.expand(S.series(expr,s,0,2*k+1).removeO()).coeff(s,2*k)==c and c!=0
        coeffs.append(str(c))
    bounds=[cutoff_bound(3,7,L) for L in [1,2,3]]
    assert [q['Newton_steps'] for q in bounds]==[0,1,2]
    assert [q['corrected_branch_degree_bound'] for q in bounds]==[49,931,17689]
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        ordinary_tail_ideal_controls=controls,nonradical_first_gate_retained=True,nilpotent_newton_controls=lifts,
        wrong_transpose_violates_error_squaring=wrong['errors'][1]['error_order_at_least']<4,wrong_transpose_errors=wrong['errors'],
        guard_refusals=dict(tail_orders=rejected,nilpotence_order9=refusal,foundational_field_branch49_not_computed=True),
        finite_parity_witness_controls=moments,even_branch_coefficient_controls=coeffs,explicit_bound_controls=bounds,
        exact_scaling_ideal_formula='I_N=(e**(N+1)) for every integer N>=0;written proof',sharp_nilpotence_slope='L-1 for exact nilpotence order L>=1;written proof',
        finite_ordinary_tail_generator_implemented=True,guarded_nilpotent_Newton_implemented=True,explicit_thickness_bound_calculator_implemented=True,
        supporting_checker_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ['verify_weighted_cotangent_inverse.py','verify_gaussian_keller_endomorphisms.py','verify_polynomial_gaussian_symmetries.py']},
        universal_theorem_verified_by_finite_controls=False,all_index_numerical_sweep_performed=False,general_inverse_solver_implemented=False,
        arbitrary_ideal_nilpotence_solver_implemented=False,new_Gaussian_moment_preserving_map_constructed=False,minimum_GMC_dimension_or_degree_improved=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (here/'nilpotent_inverse_sections_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','nonradical_first_gate_retained','wrong_transpose_violates_error_squaring','finite_parity_witness_controls','exact_scaling_ideal_formula']}))

if __name__=='__main__':main()
