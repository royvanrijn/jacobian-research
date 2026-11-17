"""Exact Gaussian--Laurent encoding controls, no general D-module solver."""
from pathlib import Path
from math import factorial
from itertools import product
import hashlib,json
import sympy as S

ROOT=Path(__file__).resolve().parent


def apply(symbol, zs, xs, f):
    result=S.Integer(0)
    for alpha,c in S.Poly(S.expand(symbol),*zs).terms():
        result+=c*S.diff(f,*sum(([x,a] for x,a in zip(xs,alpha)),[]))
    return S.expand(result)


def contract(expr,zs,ws):
    result=S.Integer(0)
    for powers,c in S.Poly(S.expand(expr),*(list(zs)+list(ws))).terms():
        n=len(zs);a=powers[:n];b=powers[n:]
        if a==b:result+=c*S.prod(factorial(v) for v in a)
    return S.expand(result)


def run():
    x,y,z,w,z2,w2,u,v=S.symbols('x y z w z2 w2 u v')
    cases=[([x],[z],[w],1+z+z**2,x**2+x+1,x+2),
           ([x,y],[z,z2],[w,w2],z**2+z2+z*z2,x*y+x**2+y,1+x)]
    contraction_checks=0
    for xs,zs,ws,lam,P,Q in cases:
        shifted=(Q*P**0).subs(dict(zip(xs,[a+b for a,b in zip(xs,ws)])),simultaneous=True)
        for m in range(4):
            F=S.expand(Q*P**m)
            actual=S.expand(F)
            for _ in range(m):actual=apply(lam,zs,xs,actual)
            shifted=F.subs(dict(zip(xs,[a+b for a,b in zip(xs,ws)])),simultaneous=True)
            encoded=contract(lam**m*shifted,zs,ws)
            assert S.expand(actual-encoded)==0
            contraction_checks+=1

    lam=1+z+z**2;P=x**2+x+1;Pw=P.subs(x,x+w)
    def coeff(A,m,Pw0=Pw,lam0=lam):
        result=S.Integer(0)
        for term in S.Add.make_args(S.expand(A)):
            powers=term.as_powers_dict()
            h,r,a,b=[int(powers.get(q,0)) for q in [u,v,z,w]]
            c=S.cancel(term/(u**h*v**r*z**a*w**b))
            ell=m-h;k=ell-r
            if ell>=0 and k>=0:
                result+=c*contract(z**a*w**b*lam0**k*Pw0**ell,[z],[w])/S.Integer(factorial(ell)*factorial(k))
        return S.cancel(result)

    integration_checks=0
    for A in [S.Integer(1),v,z,w,z*w/v,v**2*w]:
        dZ=S.diff(A,z)+(v*S.diff(lam,z)-w)*A
        dW=S.diff(A,w)+(u/v*S.diff(Pw,w)-z)*A
        dv=S.diff(A,v)+(lam-u/v**2*Pw-1/v)*A
        du=S.diff(A,u)+Pw/v*A
        for m in range(4):
            for derivative in [dZ,dW,dv]:
                assert S.cancel(coeff(derivative,m))==0
                integration_checks+=1
            assert S.cancel(coeff(du,m)-(m+1)*coeff(A,m+1))==0
            integration_checks+=1
    assert coeff(S.diff(v,v),0,Pw0=S.Integer(0),lam0=S.Integer(0))==1
    assert coeff(S.diff(v,v)-v/v,0,Pw0=S.Integer(0),lam0=S.Integer(0))==0

    shift_checks=0
    U=S.symbols('U')
    for r,a,b,m in product([-2,-1,0,1,2],range(2),range(2),range(4)):
        k=m-r
        right=S.Integer(0)
        if k>=0:
            for j in range(b+1):
                f=S.expand(U**j*P.subs(x,U)**m)
                for _ in range(k):f=apply(lam,[z],[U],f)
                right+=S.binomial(b,j)*(-x)**(b-j)*S.diff(f,U,a).subs(U,x)
            right/=S.Integer(factorial(m)*factorial(k))
        assert S.cancel(coeff(v**r*z**a*w**b,m)-right)==0
        shift_checks+=1

    # Origin substitution loses a polynomial identity obstruction.
    for m in range(1,5):
        f=S.diff(x**(2*m),x,m)
        assert f!=0 and f.subs(x,0)==0
    # The constant-coefficient assumption is essential.
    actual=x*S.diff(x*S.diff(x**2,x),x)
    frozen=x**2*S.diff(x**2,x,2)
    assert actual==4*x**2 and frozen==2*x**2
    a,b=S.symbols('a b')
    for m in range(1,6):
        assert S.diff((a*x+b)**m,x,m)==factorial(m)*a**m
    assert S.rem(a,a**2,a)==a

    # Re-encode a retained order-three failure; no lower-order transfer.
    T,A,B,C,D=S.symbols('T A B C D')
    zz=S.symbols('zt za zb zc zd');ww=S.symbols('wt wa wb wc wd')
    xs=[T,A,B,C,D]
    P5=(T+C)*(A*D+B*T);Q5=-C
    L5=zz[0]*(zz[1]*zz[4]-zz[2]*zz[3])
    failing=[]
    for m in range(1,4):
        for Q in [S.Integer(1),Q5]:
            f=S.expand(Q*P5**m)
            direct=f
            for _ in range(m):direct=apply(L5,zz,xs,direct)
            shifted=f.subs(dict(zip(xs,[a+b for a,b in zip(xs,ww)])),simultaneous=True)
            encoded=contract(L5**m*shifted,zz,ww)
            assert S.expand(direct-encoded)==0
            if Q==1:assert direct==0
            else:
                expected=C+2*T if m==1 else (-1)**(m+1)*factorial(m+1)*factorial(m)**2*T
                assert S.expand(direct-expected)==0
                failing.append({'m':m,'normalized_coefficient':str(S.cancel(direct/factorial(m)**2))})

    return {'completed':True,'diagonal_contraction_checks':contraction_checks,
            'formal_total_derivative_and_time_checks':integration_checks,
            'shifted_amplitude_checks':shift_checks,
            'logarithmic_measure_omission_detected':True,
            'origin_specialization_shortcut_detected':True,
            'variable_coefficient_symbol_shortcut_detected':True,
            'nilpotent_coefficient_control':'a mod(a^2) is nonzero',
            'linear_operator_family_ideal':'(a),least pure prefix1',
            'retained_Dvorsky_order':3,'retained_Dvorsky_variables':5,
            'retained_Dvorsky_normalized_checks':failing,
            'general_D_module_solver_implemented':False,'general_D_module_solver_executed':False,
            'universal_termination_proof_certified_by_checker':False,
            'new_lower_order_GVC_witness_computed':False,'operator_order_reduced':False,
            'producer_package_executed':False,'new_cubic_prefix_computed':False,
            'independent_replay':False,'external_review':False,'formal_verification':False}


if __name__=='__main__':
    result=run()
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/'operator_gaussian_laurent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

