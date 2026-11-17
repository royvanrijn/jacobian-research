#!/usr/bin/env python3
"""Exact Q(i) real-image filter and small transcription controls; not universal proof."""
from pathlib import Path
import hashlib,json
import sympy as S

def real_image_filter(expr,variables):
    variables=tuple(variables)
    if not variables or len(set(variables))!=len(variables):raise ValueError('distinct nonempty coordinate list required')
    p=S.Poly(expr,*variables,domain=S.QQ_I)
    # Never conjugate the variables: they are real coordinate functions even if symbols have no assumptions.
    conjugate=sum(S.conjugate(c)*S.prod(x**e for x,e in zip(variables,a)) for a,c in p.terms())
    minors=[]
    for i in range(len(variables)):
        for j in range(i+1,len(variables)):
            m=S.Poly(S.diff(p.as_expr(),variables[i])*S.diff(conjugate,variables[j])-S.diff(p.as_expr(),variables[j])*S.diff(conjugate,variables[i]),*variables,domain=S.QQ_I)
            if not m.is_zero:minors.append({'coordinates':[i,j],'minor':str(m.as_expr())})
    degree=int(p.total_degree()) if not p.is_zero else 0
    thin=not minors
    return dict(degree=degree,zero=p.is_zero,constant=degree==0,real_image_has_interior=bool(minors),thin_image=thin,
        nonzero_minors=minors,conditional_bad_monomial_degree_bound=(degree-1 if thin and degree>=2 else None),
        conditional_bound_requires_all_order_purity=True,purity_tested=False,safety_decided=False,witness_selected=False)

def gaussian(expr,variables):
    out=0
    for a,c in S.Poly(S.expand(expr),*variables).terms():
        if any(e%2 for e in a):continue
        out+=c*S.prod(S.factorial2(e-1) if e else 1 for e in a)
    return S.expand(out)

def main():
    X,Y,T=S.symbols('X Y T');I=S.I;Z=X+I*Y;W=(X-I*Y)/2;v=(X,Y,T)
    cases=[]
    def check(name,p,interior):
        a=real_image_filter(p,v);assert a['real_image_has_interior']==interior,name
        assert a['purity_tested'] is a['safety_decided'] is a['witness_selected'] is False
        cases.append(dict(name=name,**a));return a
    check('zero',0,False);check('nonzero constant',1+I,False)
    check('real linear',X,False);check('pure circular linear',Z,True)
    for k in [2,3,5]:check('circular power '+str(k),Z**k,True)
    h=X**2+Y+1
    for d in [2,3,4]:
        a=check('real seed composition '+str(d),h+I*h**d,False)
        assert a['degree']==2*d and a['conditional_bad_monomial_degree_bound']==2*d-1
        assert gaussian((h+I*h**d)**2,v)!=0 # explicitly nonpure: no invented thin pure counterexample
    check('unary complex phase',X+I*X**3,False)
    # Coordinate symbols need not be declared real; coefficient conjugation still fixes them.
    assert len(real_image_filter(Z,(X,Y))['nonzero_minors'])==1
    for k in [1,2,3]:
        p=W-Z**k*T-Z**(2*k+1)/2
        check('safe shear '+str(k),p,True)
        for m in [1,2,3]:assert gaussian(p**m,v)==0
    long=(1+Z)*(W-(2+Z)*T**2/2)
    check('retained quartic failure with interior',long,True)
    long_controls=[]
    for m in [1,2,3]:
        assert gaussian(long**m,v)==0
        a=gaussian(Z*long**m,v);assert a==S.factorial(m);long_controls.append(str(a))
    # Affine involution repair, both coordinate-choice branches.
    repairs=[]
    for lam,mu in [(1+I,2-I),(I,1+2*I),(1,3*I)]:
        alpha=S.simplify(S.conjugate(lam)/lam);beta=S.simplify(S.conjugate(mu)-alpha*mu)
        c=1+alpha if alpha!=-1 else I;d=S.simplify(S.conjugate(c)*beta/2)
        assert S.simplify(S.conjugate(c)*alpha-c)==0
        assert S.simplify(d-S.conjugate(d)-S.conjugate(c)*beta)==0
        a=S.Poly(S.expand(c*(lam*h+mu)+d),*v)
        assert all(S.simplify(z-S.conjugate(z))==0 for z in a.coeffs());repairs.append(str(alpha))
    # Free-module division transcription: each remainder has basis degree below d.
    U,A=S.symbols('U A');poly=U**3+I*U+2;division_controls=0
    for k in range(10):
        q,rem=S.div(U**k,poly-A,U)
        assert S.expand(U**k-q*(poly-A)-rem)==0
        assert S.degree(rem,U)<3;division_controls+=1
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        cases=cases,affine_involution_controls=repairs,module_division_controls=division_controls,long_marked_controls=long_controls,
        exact_Q_i_rank_filter_implemented=True,conditional_degree_bound_implemented=True,
        normalization_algorithm_implemented=False,all_order_purity_solver_implemented=False,general_holonomic_endpoint_solver_implemented=False,
        bad_monomial_selector_implemented=False,universal_theorem_verified_by_controls=False,independent_replay=False,external_review=False,formal_verification=False,
        thin_pure_Gaussian_example_constructed=False,interior_implies_safety=False)
    path=Path(__file__).with_name('gaussian_real_curve_verification.json');path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','module_division_controls','long_marked_controls']}))

if __name__=='__main__':main()
