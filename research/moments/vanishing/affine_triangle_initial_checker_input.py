#!/usr/bin/env python3
"""Exact auxiliary checks for affine triangle-leading cubics."""
from pathlib import Path
import hashlib,itertools,json
import sympy as S


def partitions(total, minimum=1):
    if total==0:yield []
    for first in range(minimum,total+1):
        for rest in partitions(total-first,first):yield [first]+rest


def main():
    directory=Path(__file__).resolve().parent
    u,v,z,t,w=S.symbols('u v z t w')
    b,c,e,f,h,k,A,B,C,D,E,m,lu,lv,lz,q0=S.symbols('b c e f h k A B C D E m lu lv lz q0')
    P=z*(u*v+h)+b*u*u+c*v*v+e*u+f*v+k
    Q=(A*u*u+B*v*v+m*z*z+2*C*u*v+2*D*u*z+2*E*v*z)/2+lu*u+lv*v+lz*z+q0
    gradient=[S.diff(t*P-Q,x) for x in (u,v,z)]
    homogeneous=[S.expand(w*w*g.subs({u:u/w,v:v/w,z:z/w},simultaneous=True)) for g in gradient]
    infinity=[]
    for axis,others in [(u,(v,z,w)),(v,(u,z,w)),(z,(u,v,w))]:
        jac=S.Matrix([g.subs(axis,1) for g in homogeneous]).jacobian(others)
        determinant=S.factor(jac.det().subs({x:0 for x in others}))
        assert determinant!=0
        infinity.append({'axis':str(axis),'Jacobian':str(determinant),'length_when_nonzero':1})
    R=b*u*u+e*u+k-f*h/u+c*h*h/u**2
    assert S.expand(u**3*S.diff(R,u))==2*b*u**4+e*u**3+f*h*u-2*c*h*h
    cases=[];positive=[]
    for H,Bb,Cc,Ee,Ff in itertools.product([0,1],repeat=5):
        support=[power for flag,power in [(Bb,2),(Ee,1),(Ff,-1),(Cc,-2)] if flag]
        if H:
            if not support:kind='conic'
            else:
                kind='finite';finite=max(support)-min(support)
                nT=1
                nu=0 if Bb else (1 if Ee else (3 if Ff else 4))
                nv=0 if Cc else (1 if Ff else (3 if Ee else 4))
        else:
            if not Ee and not Ff:kind='axes'
            elif not(Bb or Ee) or not(Cc or Ff):kind='conic'
            else:
                kind='finite';finite=int(Bb and Ee)+int(Cc and Ff)
                nT=3 if Ee and Ff else 4
                nu=int(not Bb);nv=int(not Cc)
        row={'nonzero_h_b_c_e_f':[H,Bb,Cc,Ee,Ff],'kind':kind}
        if kind=='finite':
            assert finite+nT+nu+nv==5
            row.update(finite_rank=finite,escapes_T_u_v=[nT,nu,nv]);cases.append(row)
        else:positive.append(row)
    assert len(cases)==23 and len(positive)==9
    assert sum(x['kind']=='conic' for x in positive)==5
    assert sum(x['kind']=='axes' for x in positive)==4

    U,V,K,ee,ff,cc,mm=S.symbols('U V K ee ff cc mm')
    escape_checks=[]
    for n,F0,point,relation,action in [
        (3,K*U*V+ee*U+ff*V-mm*K*K/2,{U:-ff/K,V:-ee/K},K**3-ee*ff/mm,-3*mm*K*K/2),
        (4,K*U*V+ee*U+cc*V*V-mm*K*K/2,{U:2*cc*ee/K**2,V:-ee/K},K**4+2*cc*ee**2/mm,-mm*K*K)]:
        gb=S.groebner([relation],K,domain=S.QQ.frac_field(ee,ff,cc,mm))
        def zero(expr):assert gb.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
        for variable in [U,V,K]:zero(S.diff(F0,variable).subs(point))
        zero(F0.subs(point)-action)
        zero(S.hessian(F0,(U,V,K)).det().subs(point)-n*mm*K*K)
        escape_checks.append({'n':n,'action':str(action),'Hessian':str(n*mm*K*K)})

    # The unnormalized confluent numerator frame has degree<5, even
    # when escape nodes separate at different rates.
    X=S.Symbol('X');frame_cases=[]
    for finite in range(5):
        escapes=5-finite
        for multiplicities in partitions(finite):
            nodes=list(range(1,len(multiplicities)+1));escape_nodes=list(range(10,10+escapes))
            polynomial=S.prod((X-node)**mu for node,mu in zip(nodes,multiplicities))*S.prod(X-node for node in escape_nodes)
            basis=[]
            for node,mu in zip(nodes,multiplicities):
                for jet in range(mu):basis.append(S.cancel(polynomial/(X-node)**(mu-jet)))
            for node in escape_nodes:basis.append(S.cancel(polynomial/(X-node)))
            matrix=S.Matrix(5,5,lambda i,j:S.Poly(basis[j],X).nth(i))
            assert matrix.det()!=0
            assert all(S.degree(poly,X)<=4 for poly in basis)
            frame_cases.append({'finite_rank':finite,'escapes':escapes,'multiplicities':multiplicities,'basis_determinant':str(matrix.det())})
    assert len(frame_cases)==12
    # Exact off-diagonal factorization underlying the scale-independent bound.
    xi,xj,a0,a1,a2,delta=S.symbols('xi xj a0 a1 a2 delta')
    polynomial=S.prod(X-node for node in [xi,xj,a0,a1,a2])
    numerator=S.cancel(polynomial/(X-xi))
    rhs=delta/(xj-xi)*S.prod(1+delta/(xj-node) for node in [a0,a1,a2])
    assert S.factor(numerator.subs(X,xj+delta)/S.diff(polynomial,X).subs(X,xj)-rhs)==0

    # Fixed Laurent block on the critical hyperbola, rank4.
    x,Ss,W,r,a,ellx,elly=S.symbols('x Ss W r a ellx elly')
    conic_P=W*(x*v+c*W+a)+k
    substitution={v:(r*Ss-a)/x,W:r*W,t:r**-2}
    assert S.expand((t*(conic_P-k)).subs(substitution,simultaneous=True))==Ss*W+c*W*W
    assert S.det(S.hessian(Ss*W+c*W*W,(Ss,W)))==-1
    psi=-A*x*x/2-B*a*a/(2*x*x)-ellx*x+a*elly/x
    j=S.Symbol('j',integer=True)
    expected=j*x**j-A*x**(j+2)-ellx*x**(j+1)-a*elly*x**(j-1)+B*a*a*x**(j-2)
    assert S.simplify(x*S.diff(x**j,x)+x*S.diff(psi,x)*x**j-expected)==0
    assert S.factor(conic_P.subs({x:0,v:0,W:-a/(2*c)})-(k-a*a/(4*c)))==0
    assert S.factor(S.hessian(conic_P,(x,v,W)).det().subs({x:0,v:0,W:-a/(2*c)})+a*a/(2*c))==0
    prior=json.loads((directory/'cuspidal_affine_profiles_verification.json').read_text())
    assert prior['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_positive_seed_checker_sha256':prior['checker_sha256'],
            'normal_form':'T*(uv+h)+b*u^2+c*v^2+e*u+f*v+k',
            'rank':5,'infinity':infinity,'finite_patterns':cases,'positive_dimensional_patterns':positive,
            'escape_checks':escape_checks,'confluent_numerator_frames':frame_cases,
            'critical_conic':'w*(xy+c*w+a)+k, a!=0: fixed Laurent rank4 with t^-1 seed; one Morse point if c!=0, one quadratic-action escape if c=0',
            'critical_axes':'e=f=h=0: uniform complex gradient bound for tP-Q and bounded critical actions kt+O(1); rank5 current-limit decay reuses the retained general carrier proof',
            'resolved':'Every ternary cubic whose leading form is a product of three linearly independent linear forms has a nonzero positive pure Gaussian moment, for arbitrary lower terms.',
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Exact 32 support-pattern partition and rank counts, local infinity Jacobians, escape actions/Hessians, confluent numerator algebra, Laurent recurrence and conic identities. Complete formal frames and current-limit/seed arguments are internal written mathematics.',
            'remaining':'Nonhomogeneous line-smooth-conic leading cubics, three concurrent distinct leading lines, repeated-line leading cubics; numerical uniform gates; coefficient schemes and lower-power GVC; independent review.'}
    (directory/'affine_triangle_profiles_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','rank','resolved','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
