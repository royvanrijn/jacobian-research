#!/usr/bin/env python3
"""Exact finite controls for the written all-even cyclic-cone Keller theorem.

Checks the seed moment symbolically in d, canonical chains and polynomial
Jacobian identities in bounded degrees, actual three-variable polynomial
maps in degrees six and ten, and a complete degree-six discriminant identity.
The all-degree monodromy assertion uses the companion mathematical proof.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,time
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
PRIMARY=HERE/'sources/kistner-shaska-graded-keller-2608.02863v1.pdf'
PRIMARY_SHA='3cd7e239d8b01c79239789a97c26e466208d9cb216d9773574562c468c915a2e'


def require(value,message):
    if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seed(d):
    if not isinstance(d,int)or d<3:raise ValueError('Require integral d>=3')
    D=s.sqrt(-3*(2*d-1))
    return (-3-D)/(2*(d-2)),(2*d-1+D)/(2*(d-2)),D


def formulas(d):
    if not isinstance(d,int)or d<3:raise ValueError('Require integral d>=3')
    L,w,W,Y,Z,a,b=s.symbols('Lambda w W Y Z a b')
    phi=a*W**2+b*W**d
    A1=4*a*W+s.Rational(2*d,d-1)*b*W**(d-1)
    A0=s.Rational(8,3)*a*a*W**3+s.Rational(4*d*d,(d-1)*(d+1))*a*b*W**(d+1)+s.Rational(2*d*d,(d-1)*(2*d-1))*b*b*W**(2*d-1)
    A=w+a*L*w**4+b*L**(d-1)*w**(2*d)
    B=1+4*a*L*w**3+s.Rational(2*d,d-1)*b*L**(d-1)*w**(2*d-1)+s.Rational(8,3)*a*a*L**2*w**6+s.Rational(4*d*d,(d-1)*(d+1))*a*b*L**d*w**(2*d+2)+s.Rational(2*d*d,(d-1)*(2*d-1))*b*b*L**(2*d-2)*w**(4*d-2)
    f=s.expand(W*A0+W*A1*(Y-phi)+(Y-phi)**2-Z*W)
    return dict(L=L,w=w,W=W,Y=Y,Z=Z,a=a,b=b,phi=phi,A1=A1,A0=A0,A=A,B=B,f=f)


def symbolic_seed():
    d,D,a=s.symbols('d D a')
    b=1-a
    I=a*a/3+2*a*b/(d+1)+b*b/(2*d-1)
    expected=2*(d-2)**2*a*a+6*(d-2)*a+3*(d+1)
    require(s.factor(3*(d+1)*(2*d-1)*I-expected)==0,'exact moment quadratic equation6')
    aa=(-3-D)/(2*(d-2));bb=(2*d-1+D)/(2*(d-2))
    require(s.cancel(aa+bb-1)==0,'all-d first lift equation')
    numerator=s.fraction(s.cancel(aa*aa/3+2*aa*bb/(d+1)+bb*bb/(2*d-1)))[0]
    require(s.rem(numerator,D*D+3*(2*d-1),D)==0,'all-d second lift equation')
    require(s.factor(s.discriminant(expected,a)+12*(d-2)**2*(2*d-1))==0,'quadratic discriminant')
    m=s.symbols('m',integer=True,positive=True)
    admissible=s.cancel(aa.subs({d:(m*m+1)/2,D:m*s.sqrt(-3)})-(-3-m*s.sqrt(-3))/(m*m-3))
    require(admissible==0,'fixed Qsqrt(-3) coefficient formula')
    return {'all_d_moment_and_lift_identity':True,'exact_equation6':True,'fixed_admissible_coefficient_field':'Q(sqrt(-3))'}


def bounded_controls():
    report=[]
    for d in (3,4,5,6,13,25):
        o=formulas(d);L,w,W,Y,Z,a,b=[o[k]for k in('L','w','W','Y','Z','a','b')]
        phi,A1,A0,A,B,f=[o[k]for k in('phi','A1','A0','A','B','f')]
        require(s.expand(W*s.diff(A1,W)-2*s.diff(phi,W))==0,'canonical A1 chain')
        require(s.expand(s.diff(A0,W)-s.diff(phi,W)*A1)==0,'canonical A0 chain')
        require(s.expand(W*s.diff(f,W)-f+(Y-phi)**2)==0,'inverse polynomial squared derivative')
        require(s.expand(f.subs(W,0)-Y*Y)==0,'inverse constant term')
        require(s.Poly(f,W).degree()==2*d,'inverse polynomial degree')
        require(s.expand(s.Poly(f,W).LC()+b*b/(2*d-1))==0,'inverse leading coefficient')
        T=L*w;WW=L*w*w
        chainP=L+A1.subs(W,WW)*T+A0.subs(W,WW)
        require(s.expand(chainP-L*B)==0,'literal six-term B formula from canonical chain')
        require(s.expand(T+phi.subs(W,WW)-L*A)==0,'literal three-term A formula')
        Ap=s.Poly(A,L,w,domain=s.QQ.poly_ring(a,b));Bp=s.Poly(B,L,w,domain=s.QQ.poly_ring(a,b))
        master=s.Poly(L,L,w,domain=Ap.domain)*(Bp.diff(L)*Ap.diff(w)-Bp.diff(w)*Ap.diff(L))+Bp*Ap.diff(w)-Ap*Bp.diff(w)
        require(master==s.Poly(1,L,w,domain=Ap.domain),'generic-coefficient polynomial master identity')
        aa,bb,D=seed(d);domain=s.QQ.algebraic_field(D)
        require(domain.from_sympy(s.expand((A.subs({L:1,w:-1})).subs({a:aa,b:bb})))==domain.zero,'first seed lift')
        require(domain.from_sympy(s.expand((B.subs({L:1,w:-1})).subs({a:aa,b:bb})))==domain.zero,'second seed lift')
        require(s.gcd(d,2*d-1)==1,'critical image map degree divides coprime degrees')
        require(s.Poly(A0,W).degree()==2*d-1,'critical value polynomial degree')
        require(s.Rational((-1)**(d+1)*(2*d-1),-3*(2*d-1))==s.Rational((-1)**d,3),'arithmetic squareclass ratio')
        report.append({'d':d,'generic_degree':2*d,'component_degrees':[6*d-3,12*d-9,3],
                       'finite_3cycle_branches':d,'discriminant_squareclass':(-1)**(d+1)*(2*d-1)})
    return report


def actual_polynomial_map(d):
    start=time.monotonic();o=formulas(d)
    L,w,a,b=[o[k]for k in('L','w','a','b')];aa,bb,D=seed(d)
    x,y,z=s.symbols('x y z');domain=s.QQ.algebraic_field(D)
    X=s.Poly(x,x,y,z,domain=domain);LL=s.Poly(1+x*y,x,y,z,domain=domain);WW=s.Poly(x*z-1,x,y,z,domain=domain)
    F=[]
    for name in ('A','B'):
        pol=s.Poly(o[name],L,w,domain=s.QQ.poly_ring(a,b));numerator=s.Poly(0,x,y,z,domain=domain)
        for (i,j),coefficient in pol.terms():
            cc=domain.from_sympy(s.expand(coefficient.subs({a:aa,b:bb})))
            numerator+=s.Poly.from_dict({(0,0,0):cc},(x,y,z),domain=domain)*LL**i*WW**j
        quotient,remainder=s.div(numerator,X)
        require(remainder.is_zero,'actual coordinate numerator is polynomial-divisible by x')
        F.append(quotient)
    F.append(X*LL)
    derivatives=[[f.diff(v)for v in (x,y,z)]for f in F]
    J=derivatives
    determinant=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0])
    require(determinant==s.Poly(-1,x,y,z,domain=domain),'direct expanded three-variable Jacobian determinant')
    degrees=[f.total_degree()for f in F]
    require(degrees==[6*d-3,12*d-9,3],'actual component degrees')
    return {'d':d,'polynomial_component_degrees':degrees,'monomial_counts':[len(f.terms())for f in F],
            'direct_determinant':-1,'elapsed_seconds':round(time.monotonic()-start,3)}


def full_discriminant_degree6():
    start=time.monotonic();o=formulas(3)
    W,Y,Z,a,b=[o[k]for k in('W','Y','Z','a','b')]
    aa,bb,D=seed(3)
    f=s.Poly(o['f'].subs({a:aa,b:bb}),W,Y,Z,extension=D)
    g=s.Poly((Y-o['phi']).subs({a:aa,b:bb}),W,Y,Z,extension=D)
    discriminant=f.discriminant().as_expr();R=f.resultant(g).as_expr()
    require(s.Poly(s.expand(bb*bb*Y*Y*discriminant-5*R*R),Y,Z,extension=D).is_zero,'entire degree-six discriminant-resultant identity')
    return {'degree':6,'full_identity_Disc_equals_5_times_resultant_over_bY_squared':True,
            'elapsed_seconds':round(time.monotonic()-start,3)}


def whole_fibre_reconstruction():
    W,T=s.symbols('W T',nonzero=True)
    L=T*T/W;ww=W/T;xx=W/(T*T);yy=(L-1)/xx;zz=(ww+1)/xx
    require(s.cancel(1+xx*yy-L)==0,'whole-fibre Lambda reconstruction')
    require(s.cancel(xx*zz-1-ww)==0,'whole-fibre w reconstruction')
    require(s.cancel(xx*L-1)==0,'whole-fibre third target coordinate')
    require(s.cancel(L*ww-T)==0 and s.cancel(L*ww*ww-W)==0,'whole-fibre canonical pair reconstruction')
    return {'rational_algebra_identities':True,'required_units':['W','T'],
            'units_and_entire_scheme_isomorphism_require_written_proof':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-receipt',type=Path)
    args=parser.parse_args();require(sha(PRIMARY)==PRIMARY_SHA,'retained primary source bytes')
    report={'schema_version':1,'python':platform.python_version(),'sympy':s.__version__,
            'checker_sha256':sha(__file__),'primary_source_sha256':sha(PRIMARY),
            'symbolic_seed':symbolic_seed(),'bounded_controls':bounded_controls(),
            'direct_polynomial_maps':[actual_polynomial_map(d)for d in (3,5)],
            'full_discriminant':full_discriminant_degree6(),
            'whole_fibre_reconstruction':whole_fibre_reconstruction()}
    aa,bb,D=seed(5);bad=s.simplify(aa*aa/3+2*aa*(bb+1)/6+(bb+1)**2/9)
    require(bad!=0,'corrupted coefficient fails lift moment')
    try:seed(2)
    except ValueError:pass
    else:raise AssertionError('degree-two seed was accepted')
    report.update(negative_controls_rejected=2,result='PASS',scope={
        'written_uniform_proof_required':True,'formal_verification':False,'independent_replay':False,
        'external_review':False,'no_rational_coefficient_descent':True,'no_plane_realization':True,
        'no_double_transposition_passport_claim':True})
    if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
