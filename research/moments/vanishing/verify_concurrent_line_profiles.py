#!/usr/bin/env python3
"""Exact profile algebra for three concurrent leading lines; not a proof assistant."""
from pathlib import Path
import hashlib,json
import sympy as S


def main():
    directory=Path(__file__).resolve().parent
    u,v,z,t,r,U,V,W,Ss=S.symbols('u v z t r U V W Ss')
    a,b,c,d,q,f,k,m,A,B,C,D,E,lu,lv,lz,q0=S.symbols('a b c d q f k m A B C D E lu lv lz q0')
    F=a*u**3+b*u*u*v+c*u*v*v+d*v**3
    Q=(A*u*u+B*v*v+m*z*z+2*C*u*v+2*D*u*z+2*E*v*z)/2+lu*u+lv*v+lz*z+q0
    H=S.Symbol('H')
    J=D*u+E*v+lz
    binaryQ=(A*u*u+B*v*v+2*C*u*v)/2+lu*u+lv*v+q0
    full=t*(z*H+F+q*v*v+f*v+k)-Q
    completed=full.subs(z,z+(t*H-J)/m)
    effective=t*(F+q*v*v+f*v+k)-binaryQ+(t*H-J)**2/(2*m)
    assert S.expand(completed-(-m*z*z/2+effective))==0
    # Schur complement remains nondegenerate whenever Q does and m!=0.
    metric=S.Matrix([[A,C,D],[C,B,E],[D,E,m]])
    schur=S.Matrix([[A-D*D/m,C-D*E/m],[C-D*E/m,B-E*E/m]])
    assert S.factor(metric.det()-m*schur.det())==0
    # Binary squarefree cubic gradient has length4, illustrated in its
    # exhaustive projective normal form uv(u+v).
    model=u*v*(u+v)
    gb=S.groebner([S.diff(model,u),S.diff(model,v)],u,v)
    assert [str(poly.as_expr()) for poly in gb.polys]==['u**2 - v**2','2*u*v + v**2','v**3']
    assert len([1,u,v,v*v])==4

    # H=u not dividing F3: remove uv^2, so F3=a*u^3+b*u^2*v+d*v^3.
    zeta=S.Symbol('zeta')
    R=a+b*zeta+d*zeta**3
    U0=-1/(3*m*R);V0=zeta*U0
    phi=a*U**3+b*U*U*V+d*V**3+U*U/(2*m)
    relation=3*d*zeta*zeta+b
    gr=S.groebner([relation],zeta,domain=S.QQ.frac_field(a,b,d,m))
    def zero(expr):assert gr.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    point={U:U0,V:V0}
    zero(S.diff(phi,U).subs(point));zero(S.diff(phi,V).subs(point))
    zero(phi.subs(point)-U0**2/(6*m))
    zero(S.hessian(phi,(U,V)).det().subs(point)+6*d*zeta*U0/m)
    assert S.factor(S.resultant(R,relation,zeta))==d*d*(27*a*a*d+4*b**3)
    # b=0: the one nonzero escaping point has A2 type.
    degenerate=phi.subs(b,0)
    primary={U:-1/(3*a*m),V:0}
    assert S.factor(degenerate.subs(primary))==1/(54*a*a*m**3)
    assert S.factor(S.diff(degenerate,U,2).subs(primary))==-1/m
    assert S.diff(degenerate,V,3)==6*d
    unfold_orders=[]
    for order in [1,2,3,4,5]:
        ratio=S.Rational(2)-S.Rational(3*order,4)
        perturb=S.Rational(order)-S.Rational(8,3)
        if order<=2:assert ratio>0;kind='two shrinking Gaussian charts'
        else:assert perturb>0;kind='fixed cubic rank2 chart'
        unfold_orders.append({'ord_unfolding':order,'kind':kind,'Gaussian_width_over_separation_power':str(ratio),'fixed_cubic_perturbation_power':str(perturb)})

    # H=u divides F3: remove u^2v and retain nonzero a,c.
    Pf=z*u+a*u**3+c*u*v*v+q*v*v+f*v+k
    phase=t*Pf-Q
    grad=[S.diff(phase,x) for x in [u,v,z]]
    # q!=0: two secondary escapes, in addition to one primary and one finite.
    res=[S.expand(g.subs({t:r**-2,u:U,v:V/r,z:W/r**2},simultaneous=True)*r**power)
         for g,power in zip(grad,[4,3,2])]
    lead=[S.expand(expr).coeff(r,0) for expr in res]
    assert all(S.expand(x-y)==0 for x,y in zip(lead,[W+c*V*V,2*V*(U*c+q),U-W*m]))
    point={U:-q/c,W:-q/(c*m)}
    rel=V*V-q/(m*c*c)
    gq=S.groebner([rel],V,domain=S.QQ.frac_field(q,c,m))
    def zeroq(expr):assert gq.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for expression in lead:zeroq(expression.subs(point))
    det=S.Matrix(lead).jacobian([U,V,W]).det()
    zeroq(det.subs(point)-4*q)
    leading_action=W*U+c*U*V*V+q*V*V-m*W*W/2
    zeroq(leading_action.subs(point)-q*q/(2*m*c*c))

    # q=0,f!=0: three secondary escapes.
    res=[S.expand(g.subs(q,0).subs({t:r**-3,u:r*U,v:V/r,z:W/r**2},simultaneous=True)*r**power)
         for g,power in zip(grad,[5,3,2])]
    lead=[S.expand(expr).coeff(r,0) for expr in res]
    assert all(S.expand(x-y)==0 for x,y in zip(lead,[W+c*V*V,2*U*V*c+f,U-W*m]))
    point={U:-m*c*V*V,W:-c*V*V}
    gf=S.groebner([V**3-f/(2*m*c*c)],V,domain=S.QQ.frac_field(f,c,m))
    def zerof(expr):assert gf.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for expression in lead:zerof(expression.subs(point))
    zerof(S.Matrix(lead).jacobian([U,V,W]).det().subs(point)-6*m*c*c*V*V)
    action=W*U+c*U*V*V+f*V-m*W*W/2
    zerof(action.subs(point)-3*f*V/4)
    # Fixed parabolic critical curve at q=f=0 has quartic restricted Q.
    curveQ=S.expand(Q.subs({u:0,z:-c*v*v}))
    assert S.Poly(curveQ,v).degree()==4 and S.Poly(curveQ,v).LC()==m*c*c/2
    # Exact normal scaling around that curve.
    parabolic=S.expand((t*(Pf-k)).subs({q:0,f:0,t:r**-2,u:r*U,z:r*Ss-c*v*v},simultaneous=True))
    assert parabolic==Ss*U+a*r*U**3
    # Primary escape shared by every factor-H profile.
    phif=a*U**3+c*U*V*V+U*U/(2*m)
    primary={U:-1/(3*a*m),V:0}
    assert S.factor(phif.subs(primary))==1/(54*a*a*m**3)
    assert S.factor(S.hessian(phif,(U,V)).det().subs(primary))==2*c/(3*a*m*m)
    # Positive zero-exponential scalar seed powers at finite and parabolic modes.
    positive=[S.Rational(1)]+[1+S.Rational(j,degree) for degree in [2,3] for j in range(1,degree)]
    assert all(power>0 for power in positive)
    prior=json.loads((directory/'affine_triangle_profiles_verification.json').read_text())
    assert prior['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_numerator_frame_checker_sha256':prior['checker_sha256'],
            'rank':4,'infinity_length':4,'constant_H_schur_complement':'det Q=m*det Schur',
            'nonfactor_H':'finite A_mu rank2 plus escape rank2; b!=0 gives two Morse escapes, b=0 gives a complete A2 unfolding split',
            'A2_unfolding_threshold':unfold_orders,
            'factor_H_profiles':[
                {'condition':'q!=0','finite':1,'secondary_escapes':2,'primary_escape':1,'secondary_action':'q^2/(2*m*c^2)*t^2'},
                {'condition':'q=0,f!=0','finite':0,'secondary_escapes':3,'primary_escape':1,'secondary_action':'3*f*V_i/4*t^(4/3)'},
                {'condition':'q=f=0','fixed_parabolic_block':3,'primary_escape':1,'zero_exponential_seed_power':'t^-1'}],
            'primary_action':'t^4/(54*a^2*m^3)',
            'resolved':'Every ternary cubic with three concurrent distinct leading lines is GMC-safe. Nonisotropic-vertex cases and F2(vertex)!=0 cases are pure-nonvanishing; isotropic F2(vertex)=0 uses the retained W-linear theorem.',
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Exact Gaussian completion, binary Jacobian algebra, squarefree resultant, actions and Hessians, secondary rescalings, quartic restriction, positive powers and A2 threshold arithmetic. Formal normal forms, complete frames and the seed argument are written mathematics.',
            'remaining':'Nonhomogeneous line-smooth-conic and repeated-line leading cubic GMC; numerical global pure gates; coefficient schemes; lower-power GVC; independent review/formalization.'}
    (directory/'concurrent_line_profiles_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','rank','resolved','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
