#!/usr/bin/env sage-python
"""Construct the old-degree3 triangle pencil and interpolate its mod-p j-map.

One divisor O+P6+P8 on11952, zero P7. At most64 finite samples at1009,
49 good genus-one samples determine j of degree<=24; 300-second cap.
Any pole collision, failed genus/point computation or degree loss remains
unresolved. This prototype does not assert a302 parent or a Q(t) MW basis.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time
from sage.all import EllipticCurve,FunctionField,GF,PolynomialRing,QQ,ZZ,gcd,lcm,matrix,vector

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
OUT=ROOT/'artifacts/generated-results/elkies-k3-curve302-triangle-6-8-mw14-mod1009-v1.json'
LOCAL=ROOT/'artifacts/local/elkies-k3/curve302-triangle-6-8-mw14'
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS


def make_pencil(field):
    d=json.loads(SOURCE.read_text());R=PolynomialRing(field,'u');u=R.gen();K=R.fraction_field()
    def rat(v):return K(R([field(QQ(c)) for c in v['numerator_coefficients_low_to_high']])/R([field(QQ(c)) for c in v['denominator_coefficients_low_to_high']]))
    A=R([field(QQ(c)) for c in d['weierstrass_model']['A_coefficients_low_to_high']]);B=R([field(QQ(c)) for c in d['weierstrass_model']['B_coefficients_low_to_high']])
    delta=4*A**3+27*B**2;assert delta.degree()==24 and delta.gcd(delta.derivative())==1
    E=EllipticCurve(K,[A,B]);P,Q,Z=[E(rat(d['sections']['records'][i]['X']),rat(d['sections']['records'][i]['Y'])) for i in [5,7,6]]
    def meet_o(P):
        h=P[0].denominator();h=h//h.gcd(h.derivative());assert h.degree()==1;return -h[0]/h[1]
    rp,rq,lam=meet_o(P),meet_o(Q),meet_o(P-Q);assert len({rp,rq,lam})==3
    def principal(P,r):
        a=(u-r)*(-P[1]/P[0]);return a(r)/(u-r)+a.derivative()(r)
    kp,kq=principal(P,rp),principal(Q,rq)
    c1=(lam-rp)/((u-rp)*(u-lam));c2=-(lam-rq)/((u-rq)*(u-lam));c0=(kp(lam)-kq(lam))/(u-lam)-c1*kp-c2*kq
    zero_parameter=c1*(Z[1]+P[1])/(Z[0]-P[0])+c2*(Z[1]+Q[1])/(Z[0]-Q[0])+c0
    assert max(zero_parameter.numerator().degree(),zero_parameter.denominator().degree())==1
    T=PolynomialRing(K,'s');s=T.gen();X=PolynomialRing(T,'x');x=X.gen()
    den=(x-P[0])*(x-Q[0]);ny=c1*(x-Q[0])+c2*(x-P[0]);n0=c1*P[1]*(x-Q[0])+c2*Q[1]*(x-P[0])+c0*den
    f,rem=((s*den-n0)**2-(x**3+A*x+B)*ny**2).quo_rem(den);assert not rem and f.degree()==3
    coefficients=[f[i] for i in range(4)]
    common=lcm([c.denominator() for poly in coefficients for c in poly.list()])
    coefficients=[poly*common for poly in coefficients]
    content=gcd([c.numerator() for poly in coefficients for c in poly.list()])
    coefficients=[poly/content for poly in coefficients]
    assert all(c.denominator()==1 for poly in coefficients for c in poly.list())
    return K,f,coefficients,zero_parameter,(c1,c2,c0),(rp,rq,lam)


def fibre_j(field,K,f,value):
    K0=FunctionField(field,'u');XX=PolynomialRing(K0,'x');x=XX.gen()
    poly=sum(K0(f[i](field(value)))*x**i for i in range(4))
    if poly.degree()!=3 or not poly.is_irreducible():return None
    L=K0.extension(poly/poly.leading_coefficient(),'x')
    if L.genus()!=1:return None
    place=L.get_place(1)
    if place is None:return None
    D=place.divisor();xx=next(g for g in (2*D).basis_function_space() if g.valuation(place)==-2)
    yy=next(g for g in (3*D).basis_function_space() if g.valuation(place)==-3)
    V,frm,to=(6*D).function_space();monomials=[L(1),xx,yy,xx**2,xx*yy,xx**3,yy**2]
    kernel=matrix(field,[to(g) for g in monomials]).left_kernel();assert kernel.dimension()==1
    rel=kernel.basis()[0];c0,cx,cy,cx2,cxy,cx3,cy2=rel
    assert cy2 and cx3 and sum(c*g for c,g in zip(rel,monomials))==0
    scale=-cx3/cy2;e=EllipticCurve(field,[cxy/cy2,-cx2/cy2,cy*scale/cy2,-cx*scale/cy2,-c0*scale**2/cy2])
    return int(e.j_invariant())


def build(prime=1009):
    started=time.monotonic();field=GF(prime);K,f,coefficients,zparam,c,intersections=make_pencil(field)
    samples=[];skipped=[];LOCAL.mkdir(parents=True,exist_ok=True)
    for value in range(64):
        j=fibre_j(field,K,f,value)
        if j is None:skipped.append(value)
        else:samples.append([value,j])
        (LOCAL/f'j_samples{prime}.json').write_text(json.dumps({'samples':samples,'skipped':skipped},sort_keys=True)+'\n')
        print('SAMPLE',value,j,'good',len(samples),'seconds',round(time.monotonic()-started,1),flush=True)
        if len(samples)==52:break
    assert len(samples)>=49
    R=PolynomialRing(field,'s');s=R.gen()
    rows=[[field(a)**i for i in range(25)]+[-field(j)*field(a)**i for i in range(25)] for a,j in samples[:49]]
    kernel=matrix(field,rows).right_kernel();assert kernel.dimension()>=1
    v=next(v for v in kernel.basis() if any(v[25:]));n=R(list(v[:25]));d=R(list(v[25:]));common=n.gcd(d);n=n//common;d=d//common
    assert all(n(a)==field(j)*d(a) and d(a) for a,j in samples)
    degree=max(n.degree(),d.degree());assert degree<=24
    target=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    assert field(target.discriminant())!=0
    comparison=n-field(target.j_invariant())*d
    roots=[int(a) for a in field if comparison(a)==0]
    infinity=comparison[24]==0
    return {'schema':'curve302.triangle-6-8-mw14-modular-probe.v1',
        'status':f'EXCLUDED_MOD{prime}' if degree==24 and not roots and not infinity else 'MODULAR_SURVIVOR_OR_DEGREE_LOSS',
        'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'prime':prime,'divisor':'O+P6+P8','zero':'P7','old_degree':3,
        'intersections':list(map(int,intersections)),
        'pencil_coefficients_c1_c2_c0':list(map(str,c)),'zero_parameter':str(zparam),
        'trigonal_coefficients_x_s_u_low_to_high':[[list(map(int,a.numerator().list())) for a in poly.list()] for poly in coefficients],
        'samples':samples,'skipped_samples':skipped,'interpolation_samples':49,
        'j_numerator_coefficients_low_to_high':list(map(int,n.list())),
        'j_denominator_coefficients_low_to_high':list(map(int,d.list())),
        'j_degree':int(degree),'target_j_mod_p':int(field(target.j_invariant())),
        'comparison_coefficients_low_to_high':list(map(int,comparison.list())),
        'finite_target_roots':roots,'infinity_possible':bool(infinity),
        'boundary':'Finite-field construction and exact genus-one fibre conversions. Rational interpolation uses the degree<=24 bound from the K3 elliptic pencil. A degree24 no-root result is an inverse witness only after the characteristic-zero triangle pencil and its good reduction are certified. No302 parent or explicit characteristic-zero Weierstrass MW basis is asserted by this prototype.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');parser.add_argument('--prime',type=int,choices=[1009,1013,1021],default=1009);args=parser.parse_args();signal.alarm(300)
    OUT=OUT.with_name(f'elkies-k3-curve302-triangle-6-8-mw14-mod{args.prime}-v1.json')
    result=build(args.prime)
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)
