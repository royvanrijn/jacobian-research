#!/usr/bin/env python3
"""Exact algebra for the nodal critical profiles, not a period certificate."""
from pathlib import Path
import hashlib
import json
import time
import sympy as S


def run():
    started = time.monotonic()
    here = Path(__file__).resolve().parent
    u, v, z, t, r, U, V, K, w = S.symbols('u v z t r U V K w')
    a, d, h, m, e, f = S.symbols('a d h m e f', nonzero=True)
    b, c, j, k = S.symbols('b c j k')
    aa, bb, ab, az, bz, lu, lv, lz, q0 = S.symbols('aa bb ab az bz lu lv lz q0')
    q = (aa*u*u+bb*v*v+m*z*z+2*ab*u*v+2*az*u*z+2*bz*v*z)/2+lu*u+lv*v+lz*z+q0
    G = a*u**3+d*v**3+b*u**2+c*v**2+j*u*v+e*u+f*v+k
    P = z*(u*v+h)+G
    # Mixed cubic terms are first absorbed in z; the quadratic uv term is too.
    assert S.expand(P.subs(z,z-j)-(z*(u*v+h)+G-j*u*v-j*h)) == 0
    P = P.subs(j,0)
    phase = t*P-q
    gradient = [S.diff(phase,x) for x in (u,v,z)]
    hessian = S.hessian(phase,(u,v,z))
    # Projective gradient intersection at the unique leading node.
    projective = []
    for eq in gradient:
        poly = S.Poly(eq,u,v,z)
        homog = sum(coeff*u**nu[0]*v**nu[1]*z**nu[2]*w**(2-sum(nu))
                    for nu,coeff in poly.terms())
        projective.append(S.expand(homog.subs(z,1)))
    jac = S.Matrix(projective).jacobian((u,v,w)).subs({u:0,v:0,w:0})
    assert S.expand(jac.det()-t*t*m) == 0
    # The centered conic has six finite critical points with multiplicity.
    R = S.expand(G.subs(j,0).subs(v,-h/u))
    sextic = S.expand(u**4*S.diff(R,u))
    assert S.Poly(sextic,u).degree()==6
    assert S.Poly(sextic,u).LC()==3*a
    assert S.Poly(sextic,u).coeff_monomial(1)==3*d*h**3
    assert S.expand(sextic-(3*a*u**6+2*b*u**5+e*u**4+f*h*u**2-2*c*h*h*u+3*d*h**3))==0
    # Exact Gaussian elimination, including metric cross terms and means.
    H = u*v+h
    binary = q.subs(z,0)
    linear = az*u+bz*v+lz
    effective = t*P.subs(z,0)-binary+(t*H-linear)**2/(2*m)
    assert S.expand(phase-(-m*(z-(t*H-linear)/m)**2/2+effective))==0
    s = S.Symbol('s')
    chart = effective.subs({t:1/r,v:(r*s-h)/u},simultaneous=True)
    conic_width = S.factor(S.diff(chart,s,2))
    assert S.simplify(conic_width.subs(r,0)-1/m)==0
    profiles = []
    # h!=0: unique Laurent branch at infinity, including metric cross terms.
    sub = {t:r**-1,u:r*U,v:r*V,z:K/r}
    eqs = [S.expand(eq.subs(sub,simultaneous=True)*r) for eq in gradient]
    leading = [S.expand(eq.subs(r,0)) for eq in eqs]
    assert leading == [K*V+e-az*K,K*U+f-bz*K,h-m*K]
    lead_jac = S.Matrix(leading).jacobian((U,V,K))
    assert S.expand(lead_jac.det()-m*K*K)==0
    action = S.expand(phase.subs(sub,simultaneous=True)*r*r).subs(r,0)
    assert S.expand(action-(h*K-m*K*K/2))==0
    assert S.simplify(action.subs(K,h/m)-h*h/(2*m))==0
    assert S.expand(hessian.det().subs(sub,simultaneous=True)*r**4).subs(r,0)==m*K*K
    profiles.append({'case':'h!=0','finite_Milnor_sum':6,'escape_count':1,
                     'T_scale':'t','root_equation':'K=h/m','action_power':'2',
                     'leading_action':'h**2/(2*m)','hessian_power':'4',
                     'leading_hessian':'m*K**2','scaled_gradient':[str(x) for x in leading]})
    # h=0, e!=0: lowest nonzero term of G(0,v) determines all escapes.
    # Symmetry u<->v supplies the e=0,f!=0 branches.
    for degree in (1,2,3):
        n=degree+2
        g = S.Symbol('g',nonzero=True)
        axis = g*v**degree+sum(S.Symbol('q'+str(i))*v**i for i in range(degree+1,4))
        p = z*u*v+a*u**3+b*u**2+e*u+axis+k
        s = t*p-q
        grad = [S.diff(s,x) for x in (u,v,z)]
        sub = {t:r**-n,u:r**degree*U,v:r*V,z:K/r}
        eqs = [S.expand(eq.subs(sub,simultaneous=True)*r**power)
               for eq,power in zip(grad,(n,3,1))]
        assert all(S.Poly(eq,r).is_univariate for eq in eqs)
        leading = [S.expand(eq.subs(r,0)) for eq in eqs]
        assert leading == [K*V+e,K*U+degree*g*V**(degree-1),U*V-m*K]
        relation = K**n-(-1)**(degree+1)*degree*g*e**degree/m
        solutions = {V:-e/K,U:-degree*g*(-e)**(degree-1)/K**degree}
        quotient = S.groebner([relation],K)
        def reduce_rational(expr):
            num = S.factor(S.together(expr).as_numer_denom()[0])
            return S.factor(quotient.reduce(num)[1])
        lead_jac = S.Matrix(leading).jacobian((U,V,K))
        assert reduce_rational((lead_jac.det()-n*m*K*K).subs(solutions,simultaneous=True))==0
        action = S.expand((s-t*k).subs(sub,simultaneous=True)*r*r).subs(r,0)
        target = -S.Rational(n,2*degree)*m*K*K
        assert reduce_rational((action-target).subs(solutions,simultaneous=True))==0
        det_lead = S.expand(S.hessian(s,(u,v,z)).det().subs(sub,simultaneous=True)*r**(2*n+2)).subs(r,0)
        assert reduce_rational((det_lead-n*m*K*K).subs(solutions,simultaneous=True))==0
        profiles.append({'case':'h=0, e!=0, lowest v degree='+str(degree),
                         'finite_Milnor_sum':7-n,'escape_count':n,
                         'T_scale':'t**(1/'+str(n)+')','root_equation':str(relation)+'=0',
                         'action_power':str(S.Rational(2,n)),
                         'action_convention':'S=t*k+C*t**(2/n)+lower fractional powers','leading_action':str(target),'hessian_power':str(S.Rational(2*n+2,n)),
                         'leading_hessian':str(n*m*K*K),
                         'scaled_gradient':[str(x) for x in leading]})
    # Positive-dimensional P-critical boundary: finite off-axis values and axis limits.
    center_P = P.subs({h:0,e:0,f:0})
    for point,value in [({u:-2*b/(3*a),v:0,z:0},4*b**3/(27*a*a)+k),
                        ({u:0,v:-2*c/(3*d),z:0},4*c**3/(27*d*d)+k)]:
        assert all(S.simplify(S.diff(center_P,x).subs(point))==0 for x in (u,v,z))
        assert S.simplify(center_P.subs(point)-value)==0
    transverse = S.hessian(center_P,(u,v)).subs({u:0,v:0})
    assert S.expand(transverse.det()-(4*b*c-z*z))==0
    # Explicit degree-seven finite-gradient quotient: one regression, not Bezout's proof.
    model = (z*(u*v+1)+u**3+v**3)
    model_phase = 2*model-(u*u+v*v+z*z)/2
    gb = S.groebner([S.diff(model_phase,x) for x in (u,v,z)],u,v,z)
    lm = [poly.LM(order=gb.order).exponents for poly in gb.polys]
    standard = [nu for nu in __import__('itertools').product(range(8),repeat=3)
                if not any(all(x>=y for x,y in zip(nu,lead)) for lead in lm)]
    assert len(standard)==7
    receipt = {'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'projective_infinity_Jacobian':str(jac.det()),'finite_phase_critical_length':7,
               'conic_critical_sextic':str(sextic),'escape_profiles':profiles,
               'compact_conic_chart_limit_Phi_ss':'1/m',
               'swapped_profiles':'Exchange u,v; a,d; b,c; e,f and metric entries.',
               'center_line_boundary':{'conditions':'h=e=f=0','escaping_critical_points':0,
                                      'possible_axis_limits':'(m*T+lz)*(T**2-4*b*c)=0',
                                      'off_axis_values':['k+4*b**3/(27*a**2)','k+4*c**3/(27*d**2)']},
               'example_standard_monomials':[list(nu) for nu in standard],
               'proof_boundary':'The universal Bezout count, conservation of finite Milnor multiplicities, exhaustiveness from degree seven, and the center-line boundedness/limit argument remain written proofs. Symbolic identities and a single exact gradient quotient are not their formal verification.',
               'analytic_boundary':'No new GMC class excluded. Counts, Puiseux actions and local Morse Hessians do not by themselves prove a uniform sectorial decomposition or eliminate a constant period. The center-line Gaussian-weighted analysis remains open.',
               'elapsed_seconds':round(time.monotonic()-started,6)}
    (here/'nodal_gaussian_critical_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({key:receipt[key] for key in ['completed','finite_phase_critical_length','elapsed_seconds']}))


if __name__=='__main__':
    run()
