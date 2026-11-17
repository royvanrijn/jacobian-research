#!/usr/bin/env python3
"""Exact first-amplitude algebra, not a certificate of the analytic argument."""
from pathlib import Path
import hashlib
import json
import sympy as S


def main():
    directory = Path(__file__).resolve().parent
    U,V,K,r = S.symbols('U V K r')
    a,b,c,d,e,m,A,B,C,D,E,lu,lv,q0 = S.symbols('a b c d e m A B C D E lu lv q0')
    xs = (U,V,K)
    relation = K**4+2*c*e**2/m
    field = S.QQ.frac_field(a,b,c,d,e,m,A,B,C,D,E,lu,lv,q0)
    gb = S.groebner([relation],K,domain=field)

    def zero(expression):
        numerator = S.cancel(expression).as_numer_denom()[0]
        assert gb.reduce(S.expand(numerator))[1] == 0, expression

    # Full affine Gaussian, after imposing gamma=0, in the exact quartic chart.
    terms = [K*U*V+e*U+c*V**2-m*K**2/2,
             d*V**3-d*e*m*K/(2*c),
             b*U**2-E*V*K-q0,
             -D*K*U-lv*V,
             a*U**3-B*V**2/2-lu*U,
             -C*U*V,
             -A*U**2/2]
    scaled = sum(r**i*term for i,term in enumerate(terms))
    u,v,T,t,k = S.symbols('u v T t k')
    P=T*u*v+a*u**3+d*v**3+b*u**2+c*v**2+e*u+k
    Q=(A*u**2+B*v**2+m*T**2+2*C*u*v+2*D*u*T+2*E*v*T)/2
    Q+=lu*u+lv*v+d*e*m*T/(2*c)+q0
    actual=S.expand(((t*(P-k)-Q)*r**2).subs(
        {t:r**-4,u:r**2*U,v:r*V,T:K/r},simultaneous=True))
    assert S.expand(actual-scaled)==0
    assert S.expand(scaled.subs({V:-V,K:-K},simultaneous=True)-scaled.subs(r,-r))==0
    p0={U:2*c*e/K**2,V:-e/K}
    F0,F1,F2,F3=terms[:4]
    for x in xs:zero(S.diff(F0,x).subs(p0))
    zero(F0.subs(p0)+m*K**2)
    zero(F1.subs(p0))
    H=S.hessian(F0,xs).subs(p0)
    zero(H.det()-4*m*K**2)
    z1=S.Matrix([-d*e**2/K**3,-d*e**2/(2*c*K**2),-d*e/(2*c)])
    g1=S.Matrix([S.diff(F1,x).subs(p0) for x in xs])
    for x in H*z1+g1:zero(x)
    g2=S.Matrix([S.diff(F2,x).subs(p0) for x in xs])
    c1=F3.subs(p0)+(g2.T*z1)[0]+(z1.T*S.hessian(F1,xs).subs(p0)*z1)[0]/2
    c1+=sum(S.diff(F0,xs[i],xs[j],xs[k]).subs(p0)*z1[i]*z1[j]*z1[k]/6
            for i in range(3) for j in range(3) for k in range(3))
    Hdot=S.hessian(F1,xs).subs(p0)+S.Matrix(3,3,lambda i,j:
        sum(S.diff(F0,xs[i],xs[j],xs[k]).subs(p0)*z1[k] for k in range(3)))
    trace=(H.inv()*Hdot).trace()
    Xi=-4*D*c**4+4*b*c**3*d*m+2*c**3*lv+c**2*d+d**3*e**2*m
    chi=e*Xi/(2*K*c**3)
    zero(trace+d*e/(K*c))
    zero(c1-e*(Xi-c**2*d)/(2*K*c**3))
    zero(c1-trace/2-chi)
    assert S.simplify(chi.subs(K,-K)+chi)==0

    # Separate exact Taylor substitution checks on both roots of K^2.
    # These are regressions, not a replacement for the symbolic identities.
    point={a:1,b:0,c:1,d:1,e:1,m:-2,A:1,B:1,C:0,D:0,E:0,lu:0,lv:0,q0:0}
    examples=[]
    for root in [S.Integer(1),S.I]:
        pp={**point,K:root}
        Hn=H.subs(pp);zn=z1.subs(pp)
        g2n=g2.subs(pp)
        F1H=S.hessian(F1,xs).subs(p0).subs(pp)
        third=S.Matrix([sum(S.diff(F0,xs[i],xs[j],xs[k]).subs(pp)*zn[j]*zn[k]/2
                             for j in range(3) for k in range(3)) for i in range(3)])
        z2=-Hn.inv()*(g2n+F1H*zn+third)
        initial=S.Matrix([p0.get(x,x).subs(pp) for x in xs])
        shifted={x:initial[i]+r*zn[i]+r*r*z2[i] for i,x in enumerate(xs)}
        polynomial=S.Poly(S.expand(scaled.subs(point).subs(shifted,simultaneous=True)),r)
        assert S.simplify(polynomial.nth(1))==0
        assert S.simplify(polynomial.nth(3)-c1.subs(pp))==0
        determinant=S.Poly(S.expand(S.hessian(scaled,xs).subs(point).subs(shifted,simultaneous=True).det()),r)
        assert S.simplify(determinant.nth(1)/determinant.nth(0)-trace.subs(pp))==0
        examples.append({'K':str(root),'Xi':str(Xi.subs(point)),
                         'critical_action_linear_r':str(S.simplify(polynomial.nth(3))),
                         'log_determinant_linear_r':str(S.simplify(trace.subs(pp))),
                         'normalized_amplitude_linear_r':str(S.simplify(chi.subs(pp)))})
    # Both Xi opens and Xi=0 are compatible with a nondegenerate Gaussian.
    metric=S.Matrix([[A,C,D],[C,B,E],[D,E,m]])
    assert metric.subs(point).det()==-2 and Xi.subs(point)==-1
    boundary={**point,lv:S.Rational(1,2)}
    assert Xi.subs(boundary)==0 and metric.subs(boundary).det()==-2
    s=S.Symbol('s')
    values=S.Poly((27*s**2+4)*(27*s-4),s)
    assert values.degree()==3 and S.gcd(values,values.diff()).degree()==0
    prior=json.loads((directory/'nodal_morse_dominance_verification.json').read_text())
    assert prior['completed']
    output={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_morse_checker_sha256':prior['checker_sha256'],
            'gamma':'d*e*m/(2*c)-l_T = 0',
            'Xi':str(Xi),'chi':str(chi),
            'critical_action_r_coefficient':str(e*(Xi-c**2*d)/(2*K*c**3)),
            'log_Hessian_determinant_r_coefficient':str(-d*e/(K*c)),
            'saddle_shift':[str(x) for x in z1],
            'pair_parity':'F_r(U,-V,-K)=F_-r(U,V,K); opposite normalized periods have first coefficients +chi,-chi',
            'pair_orders_in_r':[5,6],
            'exact_root_regressions':examples,
            'nonempty_example':{'a':1,'b':0,'c':1,'d':1,'e':1,'m':-2,'l_T':-1,'l_v':0,'Xi':-1,'Gaussian_determinant':-2,
                                'finite_critical_value_polynomial':str(values.as_expr())},
            'unresolved_boundary_example':{'same_coefficients_except_l_v':'1/2','Xi':0,'Gaussian_determinant':-2,'pure_null_witness':False},
            'analytic_proof_certified_by_checker':False,
            'proof_boundary':'Only exact chart, parity, stationary action and determinant coefficients, and specified regressions are checked. Parameter-uniform Morse expansion, global cut matching and grouped dominance are written analysis.',
            'remaining':'Morse/distinct-value quartic simultaneous gamma=Xi=0; finite critical degeneracies or value collisions; other cubic strata, uniform gates, coefficient schemes and lower-power GVC.'}
    (directory/'nodal_quartic_resonance_verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({key:output[key] for key in ['completed','Xi','pair_orders_in_r','exact_root_regressions','analytic_proof_certified_by_checker']}))


if __name__=='__main__':
    main()
