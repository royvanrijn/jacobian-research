"""Exact replay excluding the endpoint v=0 of L3 in the chart w=1."""
from math import comb

from sage.all import GF, QQ, PolynomialRing, matrix, vector

from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows


def homogeneous(polynomial, degree):
    return polynomial.parent()({e: c for e, c in polynomial.dict().items() if sum(e) == degree})


def l3_jets(rows):
    parameter = PolynomialRing(QQ, "v")
    v_base = parameter.gen()
    field = parameter.fraction_field()
    ring = PolynomialRing(field, names=("L", "H", "tau"), order="degrevlex")
    l, h, tau = ring.gens()
    order = 8
    expansions = []
    for k, row in enumerate(rows, 7):
        blocks = [ring.zero() for _ in range(order + 1)]
        for (r, ell, h_power, v, w_power), coefficient in row.dict().items():
            tau_power, remainder = divmod(k - 1 - (2*r + ell + 2*h_power + v + 2*w_power), 2)
            assert not remainder and tau_power >= 0
            if r + tau_power > order:
                continue
            for a in range(min(ell, order - r - tau_power) + 1):
                for b in range(min(h_power, order - r - tau_power - a) + 1):
                    blocks[r] += (
                        coefficient * comb(ell, a) * comb(h_power, b) * 8**(ell - a)
                        * 17**(h_power - b) * v_base**(v + ell - a)
                        * l**a * h**b * tau**tau_power
                    )
        expansions.append(blocks)

    def truncate(f, degree):
        return ring({e: c for e, c in f.dict().items() if sum(e) <= degree})

    def evaluate(blocks, r, degree):
        result = ring.zero()
        power = ring.one()
        for block in blocks:
            result += truncate(block * power, degree)
            power = truncate(power * r, degree)
            if not power:
                break
        return result

    pivot = expansions[0][1].constant_coefficient()
    assert pivot == -QQ(12829696)/194828780125125
    solution = ring.zero()
    for degree in range(1, order + 1):
        correction = -evaluate(expansions[0], solution, degree)/pivot
        assert correction == homogeneous(correction, degree)
        solution += correction
    residuals = [evaluate(blocks, solution, order) for blocks in expansions]
    assert not residuals[0] and all(not truncate(f, 2) for f in residuals)
    return field, residuals[1:]


def verify_l3_vzero_tauleast(field,jets):
    v=field.gen()
    residue=GF(331)
    R=PolynomialRing(residue,"h")
    h=R.gen()
    pH=h**3-69*h**2-57*h+56
    cubic=homogeneous(jets[1],3)
    c_mod=R({e[1]:residue(field(a)(0)) for e,a in cubic.dict().items() if e[0]==0})
    coeff=R({e[1]:residue(field(a)(0)) for e,a in homogeneous(jets[0],3).dict().items() if e[0]==1})
    assert c_mod.monic()==pH and coeff==124*h**2+100*h-101
    assert pH.gcd(coeff)==1
    # When tau is not least, S9 forces H out of the least-valuation direction.
    assert [ (tuple(e),residue(field(a)(0))) for e,a in cubic.dict().items()
            if not e[2] and residue(field(a)(0)) ]==[((0,3,0),residue(-4))]

    assert all(not residue(field(a)(0)) for e,a in homogeneous(jets[0],3).dict().items() if e[0]!=1)

    main=PolynomialRing(field,names=("Y","K"),order="lex")
    y,k=main.gens()
    def normalized(i,degree):
        f=homogeneous(jets[i],degree)
        parity=int(i%2==0)
        ans=main.zero()
        for (a,b,c),coef in f.dict().items():
            value=coef*v**(a-parity)
            assert value.denominator().degree()==0
            ans+=value*k**a*y**b
        return ans

    # The parity needed to compare jets after dividing by v.
    for i,f in enumerate(jets):
        for (a,b,c),q in f.dict().items():
            coeff=field(q)
            assert coeff.denominator().degree()==0
            assert all((a+j)%2==int(i%2==0) for j in coeff.numerator().dict())
    cubics=[normalized(i,3) for i in range(len(jets))]
    Q3=PolynomialRing(field,"k")
    root=Q3.gen()
    exact_H,exact_K=main.ideal(cubics).groebner_basis()
    modulus = main.hom([0,root],Q3)(exact_K)
    assert modulus.is_monic() and modulus.degree()==3
    algebra=Q3.quotient(modulus,"k")
    theta=algebra.gen()
    exact_y=-main.hom([0,theta],algebra)(exact_H)
    assert all(not main.hom([exact_y,theta],algebra)(f) for f in cubics)

    reduced_ring=PolynomialRing(residue,names=("Y","K"),order="lex")
    Y,K=reduced_ring.gens()
    def modpoly(p,target):
        return target({e:residue(field(c)(0)) for e,c in p.dict().items()})
    reduced_cubics=[modpoly(f,reduced_ring) for f in cubics]
    reduced_ideal=reduced_ring.ideal(reduced_cubics)
    for f in [Y-98*K**2+101*K-62, K**3-136*K**2-27*K+75]:
        cert=f.lift(reduced_ideal)
        assert sum(a*b for a,b in zip(cert,reduced_cubics,strict=True))==f

    r3=PolynomialRing(residue,"k")
    rq=r3.gen()
    pK=rq**3-136*rq**2-27*rq+75
    assert pK.is_squarefree() and modpoly(modulus,r3)==pK
    leading=r3.quotient(pK,"k")
    K0=leading.gen()
    Y0=98*K0**2-101*K0+62
    assert leading(modpoly(exact_y.lift(),r3))==Y0
    for f in cubics:
        modpoly(f,reduced_ring)
    # Ensure integrality for arbitrary v in the maximal ideal, not just values at 0.
    for f in (modulus,exact_y.lift(),*cubics,*(normalized(i,4) for i in range(3))):
        for c in f.dict().values():
            c=field(c); d=c.denominator()
            assert d[0]
            for a in (c.numerator()/d[0]).list()+(d/d[0]).list():
                residue(a)
    ev=reduced_ring.hom([Y0,K0],leading)
    matrix_pivot=matrix(leading,[[ev(f.derivative(x)) for x in (Y,K)] for f in reduced_cubics[:2]])
    norm=lambda a:pK.resultant(a.lift())
    assert norm(matrix_pivot.det())==-191
    gradients=[ev(modpoly(normalized(i,4),reduced_ring)) for i in range(3)]
    correction=matrix_pivot.inverse()*vector(leading,[-c for c in gradients[:2]])
    obstruction=gradients[2]+sum(ev(reduced_cubics[2].derivative(x))*d for x,d in zip((Y,K),correction,strict=True))
    assert norm(obstruction)==68
    print("At 331: tau-least chart at the v=0 endpoint of L3 excluded.")


def integral_coefficients(value, value_at_zero=False):
    """Reduce regular parameter coefficients modulo 331."""
    field=value.parent()
    d=value.denominator()
    assert d[0]
    unit=d[0]
    n=value.numerator()/unit
    denominator=d/unit
    if value_at_zero:
        for a in n.list()+denominator.list():
            GF(331)(a)
        return GF(331)(n[0])
    assert denominator.degree()==0
    return {i:GF(331)(a) for i,a in n.dict().items()}


def require_membership(ideal, polynomial):
    multipliers=polynomial.lift(ideal)
    assert sum(a*b for a,b in zip(multipliers,ideal.gens(),strict=True))==polynomial


def blow_up_lleast(field,jets):
    """
    Put s=L, X=H/L, Z=tau/L, and divide by s^3.
    """
    exact=PolynomialRing(QQ,names=("X","Z","v","s"))
    x,z,v,s=exact.gens()
    output=[]
    for f in jets:
        g=exact.zero()
        for (a,b,c),p in f.dict().items():
            p=field(p)
            assert p.denominator().degree()==0
            for n,q in (p.numerator()/p.denominator()[0]).dict().items():
                GF(331)(q)
                g += q*v**n*x**b*z**c*s**(a+b+c-3)
        output.append(g)
    R=PolynomialRing(GF(331),names=exact.variable_names(),order="degrevlex")
    X,Z,v,s=R.gens()
    forms=[homogeneous(R(f),min(map(sum,f.dict()))) for f in output]
    ideal=R.ideal([g.subs({s:0}) for g in forms])
    for target in (X-8*v,Z):
        ex=next(i for i in range(1,8) if not (target**i).reduce(ideal.groebner_basis()))
        require_membership(ideal,target**ex)
    return output


def sextic_data():
    field=GF(331)
    ring=PolynomialRing(field,"z")
    z=ring.gen()
    sextic=z**6+33*z**5+47*z**4-113*z**3-154*z**2+124*z+101
    graph=139*z**5+74*z**4+159*z**3+148*z**2-14*z-165
    assert sextic.is_squarefree()
    return sextic,graph


def check_sextic(leading,eta_parameter):
    """
    Obtain exact simple centers, with all coefficients known integral.
    The eta chart has parameter eta integral; the other chart has v in
    the maximal ideal.
    """
    S=leading[0].parent()
    K=S.base_ring()
    modparam=PolynomialRing(GF(331),"u")
    u=modparam.gen()
    bar=PolynomialRing(GF(331),names=("x","z","u"))
    X,Z,U=bar.gens()
    sextic,graph=sextic_data()
    def reduce_multivariate(f):
        output={}
        for e,c in f.dict().items():
            if eta_parameter:
                for j,a in integral_coefficients(K(c)).items():
                    output[(*e,j)]=a
            else:
                output[(*e,0)]=integral_coefficients(K(c),True)
        return bar(output)
    Br=bar.ideal([reduce_multivariate(f) for f in leading])
    require_membership(Br,sextic(Z))
    require_membership(Br,X-graph(Z)-8*U*int(eta_parameter))
    exact_graph,exact_sextic=S.ideal(leading).groebner_basis()
    P=PolynomialRing(K,"z")
    z=P.gen()
    assert exact_graph.degree(S.gen(0))==1
    assert exact_graph.monomial_coefficient(S.gen(0))==1
    assert exact_sextic.degree(S.gen(0))==0
    E=S.hom([0,z],P)(exact_sextic)
    assert E.degree()==6 and E.is_monic()
    A=P.quotient(E,"theta")
    theta=A.gen()
    a=-S.hom([0,theta],A)(exact_graph)
    evaluate=S.hom([a,theta],A)
    assert all(not evaluate(f) for f in leading)
    # Reduce exact centers and jets whose coefficients are rational in
    # the parameter into the corresponding finite algebra.
    reduced_param=PolynomialRing(GF(331),"u").fraction_field()
    u=reduced_param.gen()
    V=PolynomialRing(reduced_param,"z")
    z1=V.gen()
    p1=V(sextic.list())
    W=V.quotient(p1,"theta")
    t1=W.gen()
    def reduce_poly(f):
        if eta_parameter:
            return V({e:sum(c*u**j for j,c in integral_coefficients(K(a)).items()) for e,a in f.dict().items()})
        return V({e:integral_coefficients(K(a),True) for e,a in f.dict().items()})
    assert reduce_poly(E)==p1
    assert W(reduce_poly(a.lift()))==W(V(graph.list()))+8*u*int(eta_parameter)
    def reduce_algebra(a):
        return W(reduce_poly(a.lift()))
    return E,a,evaluate,reduce_algebra,p1


def correction_jets(normal,eta_parameter):
    E,a,evaluate,reduce_algebra,p1=check_sextic([f[0] for f in normal],eta_parameter)
    A=a.parent()
    X=normal[0][0].parent().gens()
    J=matrix(A,[[evaluate(row[0].derivative(x)) for x in X] for row in normal[:2]])
    assert p1.resultant(reduce_algebra(J.det()).lift())==26
    inv=J.inverse()
    C=-inv*vector(A,[evaluate(row[1]) for row in normal[:2]])
    def first(row):
        return evaluate(row[1])+sum(c*evaluate(row[0].derivative(x)) for x,c in zip(X,C,strict=True))
    if eta_parameter:
        for c in C:
            reduce_algebra(c)
        assert not first(normal[2])
        def second(row):
            f0,f1,f2=row
            return (evaluate(f2)+sum(c*evaluate(f1.derivative(x)) for x,c in zip(X,C,strict=True))
                    + sum(C[i]*C[j]*evaluate(f0.derivative(X[i]).derivative(X[j]))/2
                          for i in range(2) for j in range(2)))
        C2=-inv*vector(A,[second(row) for row in normal[:2]])
        for c in C2:
            reduce_algebra(c)
        obstruction=second(normal[2])+sum(c*evaluate(normal[2][0].derivative(x))
                                        for x,c in zip(X,C2,strict=True))
    else:
        v=A.base_ring().gen()
        for c in C:
            reduce_algebra(c/v**2)
        obstruction=first(normal[2])/v**2
    assert p1.resultant(reduce_algebra(obstruction).lift())==182


def check_sleast(polys):
    # Put X=sx, Z=sz, and v=s eta. Keep degrees 0,1,2 of
    # the equations normalized in s.
    parameter=PolynomialRing(QQ,"eta").fraction_field()
    eta=parameter.gen()
    R=PolynomialRing(parameter,names=("x","z"),order="lex")
    x,z=R.gens()
    normal=[]
    for f in polys:
        d=min(map(sum,f.dict()))
        normal.append([R(sum(c*x**i*z**j*eta**v for (i,j,v,t),c in f.dict().items()
                          if i+j+v+t==d+q)) for q in range(3)])
    correction_jets(normal,True)


def check_vleast(polys):
    # X=v(8+A), Z=vB, s=vT, with all A,B,T,v in the maximal ideal.
    R=PolynomialRing(QQ,names=("A","B","v","T"))
    A,B,v,T=R.gens()
    transformed=[]
    quadratics=[]
    for f in polys:
        d=min(map(sum,f.dict()))
        g=sum(c*(8+A)**i*B**j*v**(i+j+vp+t-d)*T**t for (i,j,vp,t),c in f.dict().items())
        assert all(i+j+t>=2 for i,j,vp,t in g.dict())
        transformed.append(g)
        quadratics.append(R({e:c for e,c in g.dict().items() if e[0]+e[1]+e[3]==2}))
    rows=[transformed[0]]
    for g,q in zip(transformed[1:],quadratics[1:],strict=True):
        ratio=R(q/quadratics[0])
        assert q==ratio*quadratics[0]
        assert all(not i and not j and not t for i,j,vp,t in ratio.dict())
        for c in ratio.dict().values():
            GF(331)(c)
        rows.append(g-ratio*transformed[0])
    for i,f in enumerate(rows):
        d=2+(i>0)
        assert all(a+b+t>=d and (vp>=2 or a+b+t==d) for a,b,vp,t in f.dict())
    # Force T least in the normal variables.
    Mod=PolynomialRing(GF(331),names=("A","B","T"))
    a,b,t=Mod.gens()
    forms=[]
    for i,f in enumerate(rows):
        d=2+(i>0)
        forms.append(Mod({(a,b,t):GF(331)(c) for (a,b,vp,t),c in f.dict().items()
                          if a+b+t==d and not vp}))
    for leading in (a,b):
        require_membership(Mod.ideal([f.subs({leading:1,t:0}) for f in forms]),Mod.one())
    # Rescale A=Tx and B=Tz. The constant leading equations have a
    # simple exact lift and the first nonconstant terms are in Tv².
    parameter=PolynomialRing(QQ,"v").fraction_field()
    vv=parameter.gen()
    S=PolynomialRing(parameter,names=("x","z"),order="lex")
    x,z=S.gens()
    normal=[]
    for i,f in enumerate(rows):
        d=2+(i>0)
        normal.append([S(sum(c*x**a*z**b*vv**vp for (a,b,vp,t),c in f.dict().items()
                          if a+b+t==d+q)) for q in range(3)])
    correction_jets(normal,False)


def verify_vzero_endpoint(rows):
    field,jets=l3_jets(rows)
    verify_l3_vzero_tauleast(field,jets)
    polynomials=blow_up_lleast(field,jets)
    check_sleast(polynomials)
    check_vleast(polynomials)
    print("At 331: all normal directions at the v=0 endpoint of L3 are excluded.")


if __name__ == "__main__":
    verify_vzero_endpoint(reduced_residual_rows(5, 13, even_linear=True, primitive=False))
