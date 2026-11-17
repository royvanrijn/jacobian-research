"""Exact controls for multibase quadratic-fibre primitive vacua.

Unrestricted scope follows from the companion written proof. These are small
independent Wick/source controls and exact structural rejection identities.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import sympy as sp


def zero(expr):
    assert sp.cancel(sp.together(expr)) == 0, expr


def direct(phase, variables, pairs, markers, order):
    base={powers:Fraction(int(value.p),int(value.q))
          for powers,value in sp.Poly(sp.expand(phase),*variables).terms()}
    power={(0,)*len(variables):Fraction(1)}
    answer={marker:[] for marker in markers}
    for m in range(order+1):
        if m:
            updated={}
            for left,lvalue in power.items():
                for right,rvalue in base.items():
                    powers=tuple(a+b for a,b in zip(left,right))
                    updated[powers]=updated.get(powers,Fraction(0))+lvalue*rvalue
            power={powers:value for powers,value in updated.items() if value}
        for marker in markers:
            value=Fraction(0)
            for powers,coefficient in power.items():
                exponents=tuple(a+b for a,b in zip(powers,marker))
                if any(exponents[i]!=exponents[pairs+i] for i in range(pairs)):
                    continue
                if any(exponent%2 for exponent in exponents[2*pairs:]):
                    continue
                factor=math.prod(math.factorial(exponents[i]) for i in range(pairs))
                factor*=math.prod(math.prod(range(1,exponent,2))
                                  for exponent in exponents[2*pairs:])
                value+=coefficient*factor
            value/=math.factorial(m)
            answer[marker].append(sp.Rational(value.numerator,value.denominator))
    return answer


def source_period(root, alpha, source, v, beta, eta, gamma):
    expression=sp.prod(entry**power for entry,power in zip(root,alpha))*source
    for variable,power in zip(v,beta):
        expression=sp.diff(expression,variable,power)
    for variable,power in zip(eta,gamma):
        expression=sp.diff(expression,variable,power)
    return sp.factor(expression.subs({variable:0 for variable in (*v,*eta)}))


def compare(phase,variables,pairs,markers,root,source,v,eta,order,t):
    actual=direct(phase,variables,pairs,markers,order)
    for marker in markers:
        expected=source_period(root,marker[:pairs],source,v,marker[pairs:2*pairs],eta,marker[2*pairs:])
        expected=sp.series(expected,t,0,order+1).removeO().expand()
        for m in range(order+1):
            zero(expected.coeff(t,m)-actual[marker][m])


def matrix_angular_core_controls():
    t,z,w,y1,y2,x1,x2=sp.symbols("t z w y1 y2 x1 x2")
    variables=[z,w,y1,y2]
    def wick(expr):
        answer=0
        for (a,b,c,d),coefficient in sp.Poly(sp.expand(expr),*variables).terms():
            if a==b and c%2==0 and d%2==0:
                answer+=coefficient*sp.factorial(a)*sp.factorial2(c-1)*sp.factorial2(d-1)
        return sp.expand(answer)

    # Nonzero isotropic gamma is allowed in an ordinary vector fibre.
    R=z*w
    phase=R*(y1+sp.I*y2)
    source=sp.exp((x1*x1+x2*x2)/2)/(1-t*(x1+sp.I*x2))
    markers=[(sp.Integer(1),0,0),(y1,1,0),(y2,0,1),
             (y1*y1,2,0),(y1*y2,1,1),(y2*y2,0,2)]
    for marker,a,b in markers:
        expected=sp.diff(source,x1,a,x2,b).subs({x1:0,x2:0}).expand()
        for order in range(5):
            actual=wick(marker*phase**order)/math.factorial(order)
            zero(actual-expected.coeff(t,order))
    radial_source=sp.exp((x1*x1+x2*x2)/2)/(1-t*(x1+sp.I*x2))**2
    expected=sp.diff(radial_source,x1).subs({x1:0,x2:0})
    for order in range(5):
        zero(wick(R*y1*phase**order)/math.factorial(order)-expected.expand().coeff(t,order))

    # lambda R^2 has the exact M0(2lambda t) unit, including factor two.
    for order in range(5):
        coefficient=wick(R**(2*order))/math.factorial(order)
        zero(coefficient-2**order*sp.factorial2(2*order-1))

    # Pure unit does not imply a finite image in the vector-fibre case.
    # A=I, eta=-2, gamma=(1,i), u=(-2,2i) has f=1 but a pole marker.
    d=1-2*t
    gamma=sp.Matrix([1,sp.I])
    u=sp.Matrix([-2,2*sp.I])
    covariance=sp.eye(2)/d
    zero((gamma.T*covariance*gamma)[0])
    zero((u.T*covariance*u)[0])
    L=1+2*t-t*t*(gamma.T*covariance*u)[0]
    zero(L-1/d)
    zero(d*d*L*L-1)
    marked=(t*covariance*u)[0]+(t*covariance*gamma/L)[0]
    zero(marked+t*(1+2*t)/(1-2*t))
    zero(marked-(1+t-1/d))
    unsafe=(-2*R+R*(y1+sp.I*y2)+y1*y1+y2*y2-2*y1+2*sp.I*y2)
    marked_series=sp.series(marked,t,0,4).removeO().expand()
    for order in range(4):
        zero(wick(unsafe**order)/math.factorial(order)-(1 if order==0 else 0))
        zero(wick(y1*unsafe**order)/math.factorial(order)-marked_series.coeff(t,order))
        for a in range(4):
            radial=sp.factorial(a)*d**a
            zero(wick(R**a*unsafe**order)/math.factorial(order)-sp.expand(radial).coeff(t,order))
    return {
        "kind":"small exact quadratic-angular-core controls",
        "isotropic_gamma":{"A":0,"gamma":["1","i"],"phase":"ZW(Y1+iY2)","unit_and_Y_marker_orders":[0,4],"balanced_R_Y1_orders":[0,4]},
        "quadratic_R_factor":{"phase":"(ZW)^2","unit_formula":"M0(2t)","orders":[0,4]},
        "pure_unsafe_vector_boundary":{"A":"I2","eta":-2,"gamma":["1","i"],"u":["-2","2i"],"unit":"1","Y1":"1+t-(1-2t)^-1","S_powers":"a!*(1-2t)^a, a=0..3","direct_orders":[0,3],"full_image":"K[t,(1-2t)^-1] by the written all-marker formula"},
        "scope":"Finite Wick and exact source identities; no divergence or finite-image classification is inferred from these finite controls.",
    }


def gamma_shape_controls():
    """Bounded exact source controls for positive integer Gamma shapes."""
    t,x=sp.symbols("gamma_t gamma_source")
    d=1-2*t
    L=1+2*t+4*t*t/d
    zero(L-1/d)
    marked_sources={}
    for r in range(1,4):
        source=sp.exp(-2*t*x/d+x*x/(2*d))/(d*(L-t*x/d)**r)
        zero(source.subs(x,0)-d**(r-1))
        for a in range(4):
            balanced=sp.rf(r,a)/(d*L**(r+a))
            zero(balanced-sp.rf(r,a)*d**(r+a-1))
        for D in range(r+2):
            marked=sp.diff(source,x,D).subs(x,0)
            q=sp.cancel(marked*d**(D-r+1))
            explicit=0
            for b in range(D//2+1):
                for c in range(D-2*b+1):
                    a=D-2*b-c
                    explicit+=(sp.factorial(D)*(-2*t)**a*sp.rf(r,c)*t**c*d**(b+c)
                               /(sp.factorial(a)*2**b*sp.factorial(b)*sp.factorial(c)))
            zero(q-explicit)
            assert sp.Poly(q,t) is not None
            zero(q.subs(t,sp.Rational(1,2))-(-1)**D)
            denominator=sp.denom(sp.cancel(marked))
            assert sp.degree(denominator,t)==max(0,D-r+1)
            marked_sources[(r,D)]=sp.cancel(marked)

    # An independent six-variable Wick check at shape two, only small orders.
    z1,z2,w1,w2,y1,y2=sp.symbols("gz1 gz2 gw1 gw2 gy1 gy2")
    variables=[z1,z2,w1,w2,y1,y2]
    S=z1*w1+z2*w2
    phase=-2*S+S*(y1+sp.I*y2)+y1*y1+y2*y2-2*y1+2*sp.I*y2
    def wick(expr):
        answer=0
        for powers,coefficient in sp.Poly(sp.expand(expr),*variables).terms():
            a,b,c,e,f,g=powers
            if a==c and b==e and f%2==0 and g%2==0:
                answer+=(coefficient*sp.factorial(a)*sp.factorial(b)
                         *sp.factorial2(f-1)*sp.factorial2(g-1))
        return sp.expand(answer)
    expected_y=sp.series(marked_sources[(2,2)],t,0,4).removeO().expand()
    for order in range(4):
        power=sp.expand(phase**order)
        zero(wick(power)/sp.factorial(order)-sp.expand(d).coeff(t,order))
        zero(wick(y1*y1*power)/sp.factorial(order)-expected_y.coeff(t,order))
        zero(wick(S*power)/sp.factorial(order)-sp.expand(2*d*d).coeff(t,order))
    return {
        "kind":"bounded exact Gamma-shape source and Wick controls",
        "shapes":[1,2,3],"balanced_S_powers":[0,3],
        "unit":"(1-2t)^(r-1)","balanced_source":"(r)_a*(1-2t)^(r+a-1)",
        "Y1_powers":"D=0..r+1",
        "pole_factor":"Phi[Y1^D]=(1-2t)^(r-1-D)*q_D(t), q_D(1/2)=(-1)^D",
        "shape_two_direct_Wick":{"variables":6,"moment_orders":[0,3],"markers":["1","Y1^2","S"]},
        "scope":"Finite exact controls only. The written Gamma-shape and all-marker arguments prove the unrestricted image and least-bad-marker claims.",
    }


def checks():
    t,z1,z2,w1,w2,y=sp.symbols("t z1 z2 w1 w2 y")
    v1,v2,eta=sp.symbols("v1 v2 eta")
    H=sp.Matrix([1+z2*z2,1])
    root=sp.Matrix([v1+t*(1+(t+v2)**2),t+v2])
    branch=sp.Matrix([t+t**3,t])
    F=z1-z2-z2**3
    Froot=sp.expand(F.subs({z1:root[0],z2:root[1]}))
    zero(Froot-(v1-v2-v2*(t+v2)**2))
    J=sp.eye(2)-t*H.jacobian([z1,z2]).subs({z1:branch[0],z2:branch[1]})
    assert J.det()==1
    assert all(sp.Poly(entry,t) is not None for entry in J.inv())
    markers4=[(0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0),
              (0,0,0,1),(1,0,0,1),(0,0,0,2),(2,0,0,0)]
    phase=w1*(1+z2*z2)+w2+F
    compare(phase,[z1,z2,w1,w2],2,markers4,root,sp.exp(t*Froot),[v1,v2],[],3,t)

    # Polynomial, nonnilpotent A(Z) that vanishes only on the source branch.
    Dsource=1-2*t*Froot
    Gsource=Froot-root[1]/2
    source=sp.exp(t*Gsource+(t+eta)**2/(2*Dsource))/sp.sqrt(Dsource)
    phase=w1*(1+z2*z2)+w2+F*y*y+y+F-z2/2
    markers5=[(0,0,0,0,0),(1,0,0,0,0),(0,0,1,0,0),
              (0,0,0,1,0),(0,0,0,0,1),(0,0,0,0,2),
              (1,0,0,1,0),(0,0,0,2,0)]
    compare(phase,[z1,z2,w1,w2,y],2,markers5,root,source,[v1,v2],[eta],3,t)

    # A nontrivial two-base determinant and feedback root, checked at small order.
    den=1-t*t
    feedback=sp.Matrix([(v1+t+t*v2)/den,(v2+t*v1+t*t)/den])
    source=sp.exp(t*feedback[0]**2)/den
    phase=w1*(1+z2)+w2*z1+z1*z1
    compare(phase,[z1,z2,w1,w2],2,markers4[:5],feedback,source,[v1,v2],[],3,t)

    # Real polynomial matrix A that is not pointwise nilpotent.
    z,w,y1,y2,v,e1,e2=sp.symbols("z w y1 y2 v e1 e2")
    A=sp.Matrix([[0,-sp.Rational(1,2)],[-sp.Rational(1,2),-z/2]])
    Db=sp.eye(2)-2*t*A.subs(z,t)
    assert Db.det()==1
    assert Db.inv()==sp.Matrix([[1+t*t,-t],[-t,1]])
    Dv=sp.eye(2)-2*t*A.subs(z,t+v)
    es=sp.Matrix([e1,e2])
    source=sp.exp((es.T*Dv.inv()*es)[0]/2)/sp.sqrt(Dv.det())
    phase=w-y1*y2-z*y2*y2/2
    markers=[(0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,2,0),
             (0,0,1,1),(0,0,0,2),(1,1,0,0)]
    compare(phase,[z,w,y1,y2],1,markers,[t+v],source,[v],[e1,e2],4,t)
    zero(source_period([t+v],(0,),source,[v],(0,),[e1,e2],(2,0))-(1+t*t))

    # Exact structural rejection cases; no classification is inferred from jets.
    r=t/(1-t)
    zero(r-t*(1+r))
    assert sp.denom(sp.cancel(r))==1-t or sp.factor(sp.denom(sp.cancel(r)))==t-1
    # For H=Z, r=0, J=1-t, centred marker Z^mW^m has m!/(1-t)^(m+1).
    for m in range(5):
        actual=sp.diff((v/(1-t))**m/(1-t),v,m).subs(v,0)
        zero(actual-sp.factorial(m)/(1-t)**(m+1))
    # For H=1,G=Z^2, R=t^3 is polynomial, so the finite image is exp(t^3)O,
    # but the unit is nonpolynomial and safety fails.
    exponential_source=sp.exp(t*(t+v)**2)
    for m in range(4):
        normalized=sp.simplify(sp.diff(exponential_source,v,m).subs(v,0)/sp.exp(t**3))
        assert sp.Poly(normalized,t) is not None

    # Long's exact pure-but-unsafe source: one linear base marker detects it.
    Along=-(1+z)*(2+z)/2
    Dlong=sp.cancel(1-2*t*Along.subs(z,r))
    zero(Dlong-(1-t)**-2)
    zero((1-t)**2*Dlong-1)

    # Sharp finite base-marker bound for retained polynomial seeds (1-t)^M.
    b=sp.symbols("b")
    for M in [1,2]:
        N=2*M+4
        C=sum(sp.binomial(N,j+1)*b**j for j in range(N))
        rb=t/(1-t)
        fibre=1+t*(1+rb)*C.subs(b,rb)
        zero(fibre-(1-t)**-N)
        f=(1-t)**M
        zero(f**2*(1-t)**4*fibre-1)
        for j in range(M+1):
            assert sp.denom(sp.cancel(f*rb**j))==1
        bad=sp.cancel(f*rb**(M+1))
        assert sp.degree(sp.denom(bad),t)==1

    m=sp.symbols("m",integer=True,nonnegative=True)
    divergent=sp.factorial(2*m)/sp.factorial(m)
    zero(sp.combsimp(divergent.subs(m,m+1)/divergent)-2*(2*m+1))

    # Mixed zero-root W coefficient: independently test the exact shifted-Y
    # marker formula from the additional all-degree finite-image theorem.
    mz,mw,my=sp.symbols("mixed_z mixed_w mixed_y")
    mixed_markers=[(0,d,0) for d in range(4)]
    actual=direct(my+mz*mz*mw*my,[mz,mw,my],1,mixed_markers,6)
    for d in range(4):
        shifted_moment=sum(sp.binomial(d,j)*sp.factorial2(j-1)*t**(d-j)
                           for j in range(0,d+1,2))
        normalized=sp.factorial(2*d)/sp.factorial(d)*t**d*shifted_moment
        assert sp.degree(normalized,t)==2*d
        assert sp.Poly(normalized,t).LC()==sp.factorial(2*d)/sp.factorial(d)
        expected=sp.series(sp.exp(t*t/2)*normalized,t,0,7).removeO().expand()
        for order in range(7):
            zero(expected.coeff(t,order)-actual[(0,d,0)][order])
    gamma_seed=direct(mz*mw*my,[mz,mw,my],1,[(0,0,0)],6)[(0,0,0)]
    for order in range(7):
        expected=0 if order%2 else sp.factorial2(order-1)
        zero(gamma_seed[order]-expected)

    return {
        "kind":"exact structural identities and bounded independent Wick/source controls",
        "triangular_multibase":{"source_orders":[0,3],"markers":len(markers4),"sharp_all_order_cutoff":"3 deg Q; attained by Z1^d"},
        "variable_scalar_fibre":{"source_orders":[0,3],"markers":len(markers5),"A":"Z1-Z2-Z2^3; nonzero away from branch"},
        "coupled_feedback":{"source_orders":[0,3],"markers":5,"detJ":"1-t^2"},
        "nonnilpotent_matrix_fibre":{"source_orders":[0,4],"markers":len(markers),"D_branch":[["1","t"],["t","1+t^2"]]},
        "exact_boundaries":["nonpolynomial source root","inverse-Jacobian pole ladder","nonzero polynomial generator exponent","Long pure unsafe linear-base marker","sharp base-marker degree L+1 for L=1,2"],
        "encoding_obstruction":"(ZW)^2 unit coefficient ratio 2(2m+1), radius zero",
        "quadratic_angular_core":matrix_angular_core_controls(),
        "gamma_shapes":gamma_shape_controls(),
        "zero_root_mixed_W_affine":{
            "phase":"Y+Z^2 W Y","W_marker_powers":[0,3],"moment_orders":[0,6],
            "normalized_formula":"(2d)!/d!*t^d*E[(Y+t)^d]",
            "gamma_control":{"phase":"ZWY","unit_coefficient_t^(2m)":"(2m-1)!!","moment_orders":[0,6]},
            "scope":"Finite exact Wick controls. Scalar divergence and the all-degree finite-image classification use the written proof."},
        "scope":"The written residue, determinant and source-jet proofs establish all-degree classification. These controls do not infer safety or rank from finite moment prefixes.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=checks()
    here=Path(__file__).resolve()
    note=here.with_name("GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md")
    receipt={
        "status":"passed","checked_on":"2026-11-17",
        "scope":"Exact finite controls supporting the written all-degree finite-image and primitive-vacuum classifications for multibase polynomial quadratic fibres, mixed zero-root W-affine phases, and quadratic angular cores with positive-weight perturbations and positive integer Gamma shapes.",
        "command":".venv/bin/python research/moments/vanishing/verify_gaussian_multibase_quadratic_fibre_primitivity.py --output research/moments/vanishing/gaussian_multibase_quadratic_fibre_primitivity_verification.json",
        "environment":{"python":platform.python_version(),"sympy":sp.__version__},
        "checker_sha256":hashlib.sha256(here.read_bytes()).hexdigest(),
        "proof_sha256":hashlib.sha256(note.read_bytes()).hexdigest(),
        "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False,
    }
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
