#!/usr/bin/env python3
"""Auxiliary exact algebra for NODAL_MORSE_DOMINANCE.md, not analytic QA."""
from pathlib import Path
import hashlib
import json
import sympy as S


def main():
    p=Path(__file__).resolve().parent
    x1,y1,x2,y2,lam1,lam2,sigma,eps=S.symbols('x1 y1 x2 y2 lam1 lam2 sigma eps',real=True)
    alpha1,alpha2=x1+S.I*y1,x2+S.I*y2
    lam0=S.Symbol('lam0',positive=True)
    alpha0=-(lam1*alpha1+lam2*alpha2)/lam0
    delta=alpha2-alpha1
    cross=x1*y2-y1*x2
    normal=sigma*S.I*S.conjugate(delta)
    assert S.expand(S.re(normal*alpha1)-sigma*cross)==0
    assert S.expand(S.re(normal*alpha2)-sigma*cross)==0
    assert S.simplify(S.re(normal*alpha0)+(lam1+lam2)*sigma*cross/lam0)==0
    perturbed=normal+eps*S.conjugate(delta)
    assert S.expand(S.re(perturbed*alpha2)-S.re(perturbed*alpha1)-eps*((x2-x1)**2+(y2-y1)**2))==0
    triangle=[-1-S.I,2-S.I,-1+2*S.I]
    ray=(27-33*S.I)/10
    projections=[S.expand(S.re(ray*z)) for z in triangle]
    assert projections==[-6,S.Rational(21,10),S.Rational(39,10)]
    u,v,z,t,r,K,U,V=S.symbols('u v z t r K U V')
    a,d,m=S.symbols('a d m',nonzero=True)
    b,c,e,f,k=S.symbols('b c e f k')
    aa,bb,ab,az,bz,lu,lv,lz,q0=S.symbols('aa bb ab az bz lu lv lz q0')
    Q=(aa*u*u+bb*v*v+m*z*z+2*ab*u*v+2*az*u*z+2*bz*v*z)/2+lu*u+lv*v+lz*z+q0
    center=z*u*v+a*u**3+d*v**3+b*u*u+c*v*v+k
    off_axis=[]
    for point,expected in [({u:-2*b/(3*a),v:0,z:0},8*b**3/(9*a*a)),
                           ({u:0,v:-2*c/(3*d),z:0},8*c**3/(9*d*d))]:
        H=S.hessian(center,(u,v,z)).subs(point)
        assert S.factor(H.det()-expected)==0
        off_axis.append({'point':[str(point[x]) for x in (u,v,z)],'Hessian_determinant':str(expected)})
        # First critical displacement and action correction with arbitrary Q.
        g=S.Matrix([S.diff(Q,x).subs(point) for x in (u,v,z)])
        shift=H.inv()*g
        assert (H*shift-g).applyfunc(S.simplify)==S.zeros(3,1)
        variables=(u,v,z)
        shifted={x:point[x]+r*shift[i] for i,x in enumerate(variables)}
        action=S.expand((center/r-Q).subs(shifted,simultaneous=True))
        correction=S.expand(action).coeff(r,1)
        assert S.simplify(correction+(g.T*shift)[0]/2)==0
    # The only paired escape positive-action ambiguity in the four profiles.
    P=z*u*v+a*u**3+d*v**3+b*u*u+c*v*v+e*u+k
    phase=t*P-Q
    sub={t:r**-4,u:r*r*U,v:r*V,z:K/r}
    scaled=S.expand((phase-t*k).subs(sub,simultaneous=True)*r*r)
    first=S.expand(scaled).coeff(r,1)
    assert S.expand(first-(d*V**3-lz*K))==0
    numerator=S.together(first.subs(V,-e/K)-(d*e*m/(2*c)-lz)*K).as_numer_denom()[0]
    relation=K**4+2*c*e*e/m
    assert S.groebner([relation],K).reduce(S.expand(numerator))[1]==0
    # A nonempty open test point in the six-finite-critical-point profile.
    s=S.Symbol('s')
    sextic=3*u**6+2*u**5+u**4+3*u**2-4*u+3
    value_equation=u**6+u**5+u**4-s*u**3-3*u**2+2*u-1
    assert S.gcd(sextic,S.diff(sextic,u))==1
    value_poly=S.Poly(S.resultant(sextic,value_equation,u),s).primitive()[1]
    assert value_poly.degree()==6 and S.gcd(value_poly,value_poly.diff()).degree()==0
    finite_examples=[{'profile':'h!=0','finite_count':6,
                     'P':'T*(uv+1)+u^3+v^3+u^2+2*v^2+u+3*v',
                     'critical_polynomial':str(sextic),'critical_value_polynomial':str(value_poly.as_expr())}]
    for name,poly,count in [('h=0, e*f!=0',(27*s*s+4)*(27*s*s+32),4),
                            ('h=0, e!=0, f=0, c!=0',(27*s*s+4)*(27*s-4),3),
                            ('h=0, e!=0, f=c=0',27*s*s+4,2)]:
        poly=S.Poly(poly,s)
        assert poly.degree()==count and S.gcd(poly,poly.diff()).degree()==0
        finite_examples.append({'profile':name,'finite_count':count,
                                'critical_value_polynomial':str(poly.as_expr())})
    center_prior=json.loads((p/'nodal_center_line_decay_verification.json').read_text())
    critical_prior=json.loads((p/'nodal_gaussian_critical_verification.json').read_text())
    assert center_prior['completed'] and critical_prior['completed']
    assert [q['escape_count'] for q in critical_prior['escape_profiles']]==[1,3,4,5]
    assert all(row['finite_count']+profile['escape_count']==7
               for row,profile in zip(finite_examples,critical_prior['escape_profiles']))
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_center_checker_sha256':center_prior['checker_sha256'],
            'prior_critical_checker_sha256':critical_prior['checker_sha256'],
            'strict_triangle':{'values':[str(x) for x in triangle],'direction':str(ray),
                               'ordered_real_parts':[str(x) for x in projections],
                               'pure_null_witness':False},
            'universal_triangle_ray':'Re(n*alpha1)=Re(n*alpha2)=sigma*cross>0; Re(n*k)<0; perturb by epsilon*conj(alpha2-alpha1) to split growth rates',
            'off_axis_Morse_checks':off_axis,
            'quartic_escape_resonance':'gamma=d*e*m/(2*c)-l_T; opposite K share their entire positive action when gamma=0',
            'generic_finite_examples':finite_examples,
            'center_line_class':'All h=e=f=0 members, arbitrary k, transverse quadratics and affine Gaussian means',
            'finite_class':'P has finite affine critical set, all its points are Morse with pairwise distinct P-values; additionally gamma!=0 on the quartic escape profile, after the stated variable swap',
            'analytic_proof_certified_by_checker':False,
            'proof_boundary':'Triangle geometry, Hessians, finite critical/action formulas, quartic coefficient and exact squarefree example polynomials are checked. Global polynomial carriers, flat cut bases, local Morse matching, parameter-uniform tails and dominance remain internal written analysis.',
            'remaining':'Finite-critical nodal degeneracies, coincident affine critical values or the quartic gamma=0 resonance; other singular-leading cubics, numerical gates, coefficient schemes and lower-power GVC.'}
    (p/'nodal_morse_dominance_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','center_line_class','finite_class','analytic_proof_certified_by_checker']}))


if __name__=='__main__':
    main()
