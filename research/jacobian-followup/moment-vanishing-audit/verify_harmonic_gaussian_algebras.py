#!/usr/bin/env python3
"""Exact certificates for harmonic polynomial algebras; bounded controls are not proofs."""
from pathlib import Path
import hashlib,json
import sympy as S

def delta(f,x):return S.expand(sum(S.diff(f,t,2) for t in x))
def gamma(f,g,x):return S.expand(sum(S.diff(f,t)*S.diff(g,t) for t in x))
def heat(f,x):
    out=S.expand(f);term=out;j=0
    while term!=0:
        j+=1;term=delta(term,x)/(2*j);out+=term
    return S.expand(out)
def expectation(f,x):return S.expand(heat(f,x).subs(dict.fromkeys(x,0)))
def coefficient_conjugate(f,x):
    p=S.Poly(f,*x,domain=S.QQ_I)
    return S.expand(sum(S.conjugate(c)*S.prod(t**a for t,a in zip(x,b)) for b,c in p.terms()))
def scalar_certificate(expr,x):
    x=tuple(x)
    if not x or len(x)!=len(set(x)):raise ValueError('distinct nonempty coordinates required')
    p=S.Poly(expr,*x,domain=S.QQ_I);f=p.as_expr();lap=delta(f,x);null=gamma(f,f,x)
    origin=f.subs(dict.fromkeys(x,0));r=min((sum(a) for a,c in p.terms()),default=0) if not p.is_zero else None
    certified=lap==0 and null==0 and origin==0
    return dict(harmonic=lap==0,isotropic_gradient=null==0,centered=origin==0,zero=p.is_zero,
        laplacian=str(lap),gradient_square=str(null),origin_order=r,
        certified_pure_and_safe=certified,cutoff_coefficient=(str(S.Rational(1,r)) if certified and r else None),
        zero_phase_all_positive_outputs_zero=certified and p.is_zero)

def tuple_certificate(exprs,x):
    x=tuple(x);ps=[S.Poly(f,*x,domain=S.QQ_I).as_expr() for f in exprs]
    laps=[delta(f,x) for f in ps];gram=[[gamma(f,g,x) for g in ps] for f in ps]
    return dict(polynomial_algebra_harmonic=all(t==0 for t in laps) and all(t==0 for row in gram for t in row),
        laplacians=[str(t) for t in laps],pointwise_gradient_gram=[[str(t) for t in row] for row in gram])

def fixed_isotropic_gradient_test(exprs,x):
    # Every coefficient vector of every gradient lies in the span of real-point gradient values,
    # by polynomial identity duality. Their finite Gram matrix tests uniform isotropy exactly.
    vecs=[]
    for f in exprs:
        bymon={}
        for i,t in enumerate(x):
            for a,c in S.Poly(S.diff(f,t),*x,domain=S.QQ_I).terms():
                if c!=0:bymon.setdefault(a,[S.S.Zero]*len(x))[i]=c
        vecs.extend(bymon.values())
    nonzero=[]
    for i,v in enumerate(vecs):
        for j,w in enumerate(vecs[i:],i):
            z=S.simplify(sum(a*b for a,b in zip(v,w)))
            if z!=0:nonzero.append(dict(i=i,j=j,value=str(z)))
    return dict(uniform_gradient_span_isotropic=not nonzero,coefficient_vectors=len(vecs),nonzero_gram_entries=nonzero)

def main():
    x=S.symbols('x1:9');I=S.I;z=x[0]+I*x[1];u=x[2]+I*x[3];w=x[4]+I*x[5];wb=x[4]-I*x[5];t=x[6]
    h=S.expand(z*z*w+2*z*u*t-u*u*wb);v7=x[:7]
    hcert=scalar_certificate(h,v7);assert hcert['certified_pure_and_safe'] and hcert['origin_order']==3 and hcert['cutoff_coefficient']=='1/3'
    cross=fixed_isotropic_gradient_test([h],v7);assert not cross['uniform_gradient_span_isotropic']
    points=[dict.fromkeys(v7,0),dict.fromkeys(v7,0)];points[0][x[0]]=1;points[1][x[2]]=1
    grads=[[S.diff(h,q).subs(pt) for q in v7] for pt in points]
    crossvalue=S.expand(sum(a*b for a,b in zip(*grads)));assert crossvalue==-2
    compositions=[]
    for label,f,order in [('H',h,3),('H+H^2',h+h*h,3),('H^2',h*h,6),('H^3',h**3,9)]:
        c=scalar_certificate(f,v7);assert c['certified_pure_and_safe'] and c['origin_order']==order
        compositions.append(dict(name=label,**c))
    v=x[6]+I*x[7];vb=x[6]-I*x[7]
    family=[S.expand(z*v+u*w),S.expand(z*wb-u*vb)]
    fc=tuple_certificate(family,x);assert fc['polynomial_algebra_harmonic']
    assert not fixed_isotropic_gradient_test(family,x)['uniform_gradient_span_isotropic']
    f=family[0]+family[1]**2+family[0]*family[1]
    assert scalar_certificate(f,x)['certified_pure_and_safe']
    mixed=0
    for m in [1,2,3]:
        hm=S.expand(h**m);assert delta(hm,v7)==0
        for q in [S.S.One,x[0],x[6],x[0]**2,x[0]*x[6]]:
            assert expectation(S.expand(q*hm),v7)==0;mixed+=1
    # Independently check sharp monomial extraction at the first level.
    hp=S.Poly(h,*v7);a,c=hp.terms()[0];q=S.prod(b**e for b,e in zip(v7,a))
    marked=expectation(q*h,v7);assert marked==c*S.prod(S.factorial(e) for e in a) and marked!=0
    sharp_sums=[]
    for m in [1,2,3]:
        terms=S.Poly(S.expand(h**m),*v7).terms()
        norm=S.simplify(sum(S.prod(S.factorial(e) for e in a)*c*S.conjugate(c) for a,c in terms))
        assert norm>0;sharp_sums.append(str(norm))
    # Heat/log triangular inverse on nonharmonic polynomials, with rational coefficients.
    A,B=x[:2];inverse_controls=0
    for f in [A**4+A*A*B+2*B**2,A**3+I*A*B**2,A*A+B*B,(A+I*B)+(A+I*B)**2]:
        g=S.expand(f);log=0;k=0
        while g!=0:
            g=S.expand(heat(g,(A,B))-g);k+=1
            log+=2*S.Rational((-1)**(k+1),k)*g
        assert S.expand(log-delta(f,(A,B)))==0;inverse_controls+=1
    Y1,Y2=S.symbols('y1 y2');p=(A+I*B)+(A+I*B)**2
    shifted=S.expand(p.subs({A:A+Y1,B:B+Y2},simultaneous=True)-p.subs({A:Y1,B:Y2},simultaneous=True))
    for m in [1,2,3]:assert expectation(shifted**m,(A,B))==0
    # Ordinary origin-only tests are weaker: T is harmonic but has second moment one.
    real=scalar_certificate(A,(A,));assert real['harmonic'] and not real['isotropic_gradient'] and not real['certified_pure_and_safe']
    assert expectation(A*A,(A,))==1
    # A retained safe nonlinear shear falls outside the harmonic-power certificate.
    Z=A+I*B;W=(A-I*B)/2;T=x[2];shear=S.expand(W-Z*T-Z**3/2)
    sc=scalar_certificate(shear,(A,B,T));assert sc['harmonic'] and not sc['isotropic_gradient']
    assert S.expand(gamma(shear,shear,(A,B,T))+2*(T+Z*Z))==0
    for m in [1,2,3]:assert expectation(shear**m,(A,B,T))==0
    assert scalar_certificate(0,(A,))['zero_phase_all_positive_outputs_zero']
    assert not scalar_certificate(1,(A,))['certified_pure_and_safe']
    # Quotient arithmetic retains the coefficient e, instead of setting it to zero.
    e=S.symbols('e');mod=lambda f:S.rem(S.expand(f),e*e,e)
    nf=e*A*A;nil_gamma=mod(gamma(nf,nf,(A,)));nil_lap=mod(delta(nf,(A,)));nil_moment=mod(expectation(nf,(A,)))
    assert nil_gamma==0 and nil_lap==2*e and nil_moment==e
    # Harmonic cores also give polynomial differential outputs; this is not a Gaussian-only check.
    gvc_zero=0;gvc_sharp=[]
    for m in [1,2,3]:
        for q in [S.S.One,W,W*W]:
            if m>S.Poly(q,A,B).total_degree():
                out=S.expand(q*Z**(3*m))
                for k in range(m):out=delta(out,(A,B))
                assert out==0;gvc_zero+=1
        out=S.expand(W**m*Z**(3*m))
        for k in range(m):out=delta(out,(A,B))
        expected=2**m*S.factorial(m)*S.prod(3*m-k for k in range(m))*Z**(2*m)
        assert S.expand(out-expected)==0 and out!=0;gvc_sharp.append(str(S.factorial(m)*2**m*S.prod(3*m-k for k in range(m))))
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        harmonic_GVC_zero_controls=gvc_zero,harmonic_GVC_boundary_coefficients=gvc_sharp,
        scalar_compositions=compositions,published_cubic_certificate=hcert,nonuniform_cross_point_pairing=str(crossvalue),
        coefficient_gradient_gram=cross,published_quadratic_tuple_certificate=fc,mixed_zero_controls=mixed,
        sharp_first_marker_value=str(marked),positive_Fischer_controls=sharp_sums,heat_log_inverse_controls=inverse_controls,
        centered_translation_controls=3,safe_shear_outside_harmonic_class=sc,
        square_zero_mutation=dict(gradient_square=str(nil_gamma),laplacian=str(nil_lap),first_pure=str(nil_moment)),
        exact_Q_i_scalar_certificate_implemented=True,exact_Q_i_tuple_certificate_implemented=True,
        origin_order_cutoff_implemented=True,coefficient_gradient_gram_implemented=True,
        universal_theorem_verified_by_controls=False,general_safe_phase_classifier_implemented=False,
        general_harmonic_morphism_classification_implemented=False,field_gradient_only_shortcut_used=False,
        all_ring_gradient_only_criterion_valid=False,ordinary_origin_pure_prefix_two_claimed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    dest=Path(__file__).with_name('harmonic_gaussian_algebra_verification.json');dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','nonuniform_cross_point_pairing','mixed_zero_controls','positive_Fischer_controls','square_zero_mutation']}))

if __name__=='__main__':main()
