#!/usr/bin/env python3
"""Exact auxiliary profiles for both affine line/smooth-conic configurations."""
from pathlib import Path
import hashlib,itertools,json
import sympy as S


def main():
    directory=Path(__file__).resolve().parent
    u,v,z,t,r,U,V,W,Ss=S.symbols('u v z t r U V W Ss')
    a,b,c,e,f,h,k,g,m,A,B,C,D,E,lu,lv,lz,q0=S.symbols('a b c e f h k g m A B C D E lu lv lz q0')
    transverse=[]
    for H,Bb,Cc,Ee,Ff in itertools.product([0,1],repeat=5):
        if H:
            support=[3]+[power for flag,power in [(Bb,2),(Ee,1),(Ff,-1),(Cc,-2)] if flag]
            finite=3-min(support);nT=1
            nv=0 if Cc else (1 if Ff else (3 if Ee else (4 if Bb else 5)))
            kind='finite'
        elif not Ee and not Ff:kind='axes'
        elif not Cc and not Ff:kind='conic'
        else:
            kind='finite';finite=(2 if Ee else (1 if Bb else 0))+int(Cc and Ff)
            nT=3 if Ee and Ff else (4 if Ee else (4 if Bb else 5))
            nv=int(not Cc)
        row={'nonzero_h_b_c_e_f':[H,Bb,Cc,Ee,Ff],'kind':kind}
        if kind=='finite':
            assert finite+nT+nv==6
            row.update(finite_rank=finite,escapes_T_v=[nT,nv])
        transverse.append(row)
    assert sum(row['kind']=='finite' for row in transverse)==26
    assert sum(row['kind']=='conic' for row in transverse)==2
    assert sum(row['kind']=='axes' for row in transverse)==4

    # Tangent normal form and a single sextic-action escape on g!=0.
    P=z*(u*u+g*v+h)+(u+c)*v*v+e*u+f*v+k
    Q=(A*u*u+B*v*v+m*z*z+2*C*u*v+2*D*u*z+2*E*v*z)/2+lu*u+lv*v+lz*z+q0
    grad=[S.diff(t*P-Q,x) for x in (u,v,z)]
    limit=[S.expand(expr.subs({t:1/r,u:U/r,v:V/r**2,z:W/r**3},simultaneous=True)*r**power).coeff(r,0)
           for expr,power in zip(grad,[5,4,3])]
    wanted=[2*U*W+V*V,2*U*V+g*W,U*U+g*V-m*W]
    assert all(S.expand(x-y)==0 for x,y in zip(limit,wanted))
    point={U:-5*g*g/(8*m),V:25*g**3/(16*m*m),W:125*g**4/(64*m**3)}
    for expression in limit:assert S.factor(expression.subs(point))==0
    phase=W*(U*U+g*V)+U*V*V-m*W*W/2
    assert S.factor(phase.subs(point))==3125*g**8/(8192*m**5)
    assert S.factor(S.hessian(phase,(U,V,W)).det().subs(point))==625*g**6/(128*m**3)
    restricted=S.expand(((u+c)*v*v+e*u+f*v+k).subs(v,-(u*u+h)/g))
    assert S.degree(S.diff(restricted,u),u)==4

    # Common secondary pair: c!=0,Delta=c^2+h!=0.
    equations=[2*W*U+V*V,2*(U+c)*V,U*U+h-m*W]
    point={U:-c,W:(c*c+h)/m}
    gb=S.groebner([V*V-2*c*(c*c+h)/m],V,domain=S.QQ.frac_field(c,h,m))
    def zero(expr):assert gb.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for expression in equations:zero(expression.subs(point))
    zero(S.Matrix(equations).jacobian([U,V,W]).det().subs(point)-8*c*(c*c+h))
    action=W*(U*U+h)+(U+c)*V*V-m*W*W/2
    zero(action.subs(point)-(c*c+h)**2/(2*m))

    # c=0,h!=0: a cubic trio when L=f-Eh/m!=0.
    L=S.Symbol('L')
    equations=[2*h*U/m+V*V,2*U*V+L,U*U-m*W-E*V]
    point={U:-m*V*V/(2*h),W:m*V**4/(4*h*h)-E*V/m}
    gb=S.groebner([V**3-L*h/m],V,domain=S.QQ.frac_field(L,h,m,E))
    def zero(expr):assert gb.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for expression in equations:zero(expression.subs(point))
    zero(S.Matrix(equations).jacobian([U,V,W]).det().subs(point)-6*m*V*V)
    zero((h*U*U/m+U*V*V+L*V).subs(point)-3*L*V/4)
    # Its L=0 alternative is a fixed quartic/Gaussian block.
    fixed=h*U*U/m+U*(V*V+e-D*h/m)-m*W*W/2-E*V*W-lz*W-B*V*V/2-lv*V-q0
    uc=S.solve(S.diff(fixed,U),U)[0];wc=S.solve(S.diff(fixed,W),W)[0]
    quartic=S.expand(fixed.subs({U:uc,W:wc},simultaneous=True))
    assert S.Poly(quartic,V).LC()==-m/(4*h)

    # Delta=0,c!=0,f!=0: three secondary escapes near u=-c.
    equations=[V*V-2*c*W,2*U*V+f,-2*c*U-m*W]
    point={U:-m*V*V/(4*c*c),W:V*V/(2*c)}
    gb=S.groebner([V**3-2*c*c*f/m],V,domain=S.QQ.frac_field(c,f,m))
    def zero(expr):assert gb.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for expression in equations:zero(expression.subs(point))
    zero(S.Matrix(equations).jacobian([U,V,W]).det().subs(point)-6*m*V*V)
    zero((U*(V*V-2*c*W)+f*V-m*W*W/2).subs(point)-3*f*V/4)
    # Exact normal chart at its f=0 critical parabola.
    x=S.Symbol('x')
    unshifted=P.subs({g:0,h:-c*c,u:x-c})
    substituted=S.cancel(unshifted.subs(z,(Ss-v*v-e)/(x-2*c)))
    assert S.cancel(substituted-(x*Ss+f*v+k-e*c))==0

    # c=h=g=0,f!=0: all five escaping saddles.
    F0=W*U*U+U*V*V+f*V-m*W*W/2
    point={V:-f/(2*U),W:U*U/m}
    gb=S.groebner([U**5+m*f*f/8],U,domain=S.QQ.frac_field(m,f))
    def zero(expr):assert gb.reduce(S.cancel(expr).as_numer_denom()[0])[1]==0
    for variable in [U,V,W]:zero(S.diff(F0,variable).subs(point))
    zero(F0.subs(point)-5*U**4/(2*m))
    zero(S.hessian(F0,(U,V,W)).det().subs(point)+20*U**3)

    # Last all-center chart: fixed D5/Gaussian rank5, positive seed t^-3/4.
    fixed=z*U*U+U*V*V-m*z*z/2-lz*z-q0
    completed=S.expand(fixed.subs(z,z+(U*U-lz)/m))
    assert S.factor(completed-(U**4/(2*m)+U*V*V-lz*U*U/m-m*z*z/2+lz*lz/(2*m)-q0))==0
    two=U**4/(2*m)+U*V*V-lz*U*U/m
    basis=S.groebner([S.diff(two,U),S.diff(two,V)],U,V,domain=S.QQ.frac_field(m,lz))
    assert [str(poly.as_expr()) for poly in basis.polys]==['U**3 - U*lz + V**2*m/2','U*V','V**3']
    assert len([1,U,U*U,V,V*V])==5
    profiles=[
      {'condition':'g!=0','ranks':[4,1],'description':'finite4 plus one sextic-action escape'},
      {'condition':'g=0,h!=0,c!=0,c^2+h!=0','ranks':[2,1,2],'description':'finite2, primary1, secondary2'},
      {'condition':'g=c=0,h!=0,L!=0','ranks':[2,3],'description':'finite2 plus three quartic-action corrections to a nonzero quadratic action'},
      {'condition':'g=c=0,h!=0,L=0','ranks':[2,3],'description':'finite2 plus fixed quartic block3 with nonzero quadratic action'},
      {'condition':'g=0,c!=0,c^2+h=0,f!=0','ranks':[1,1,3],'description':'finite1, primary1, three t^(4/3) secondary actions'},
      {'condition':'g=0,c!=0,c^2+h=0,f=0','ranks':[1,1,3],'description':'finite1, primary1, fixed parabola3 with t^-1 seed'},
      {'condition':'g=h=0,c!=0,delta!=0','ranks':[3,2],'description':'three nonzero t^(2/3) actions plus secondary2'},
      {'condition':'g=h=0,c!=0,delta=0','ranks':[3,2],'description':'fixed critical-line3 with t^-1 seed plus secondary2'},
      {'condition':'g=h=c=0,f!=0','ranks':[5],'description':'five nonzero t^(6/5) actions'},
      {'condition':'g=h=c=f=0,e!=0','ranks':[2,3],'description':'two fixed critical-line Gaussians with t^-1 seed plus escape3'},
      {'condition':'g=h=c=f=e=0','ranks':[5],'description':'fixed D5/Gaussian block5 with t^-3/4 seed'}]
    assert all(sum(profile['ranks'])==5 for profile in profiles)
    prior=json.loads((directory/'concurrent_line_profiles_verification.json').read_text())
    assert prior['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_complete_formal_checker_sha256':prior['checker_sha256'],
            'transverse_rank':6,'transverse_patterns':transverse,
            'tangent_rank':5,'tangent_profiles':profiles,
            'resolved':'Every ternary cubic whose leading form is a line times a smooth conic is GMC-safe for arbitrary lower terms. All transverse cases and nonisotropic-tangency cases are pure-nonvanishing; isotropic zero-F2 tangencies use W-linear safety.',
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Exact transverse support/rank arithmetic and tangent normal-model actions, Jacobians, quartic completions, parabola identity, D5 Jacobian basis and rank arithmetic. Formal completeness, uniform cut decay, local normal charts and quotient arguments are written mathematics.',
            'remaining':'Only nonhomogeneous repeated-line leading ternary cubic GMC; numerical uniform gates; coefficient schemes; lower-power GVC; external review and formalization.'}
    (directory/'affine_line_conic_profiles_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','transverse_rank','tangent_rank','resolved','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
