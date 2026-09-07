#!/usr/bin/env sage-python
"""Recover the Kihara control and test all rational parameters against302.

One published arithmetic rank>=14 family, exact symbolic construction,
180-second cap. Rank14--16 candidate-subspace intake is an exact matrix
audit, not a claim that their ranks are generic ranks.
"""
import argparse
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, prod

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
from kihara_rank14 import kihara_specialization
from ecsearch.kihara import kihara_rank14_replay,verify_kihara_rank14_manifest
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS

ART=ROOT/'artifacts/generated-results/elliptic-curves'
CONTROL=ART/'kihara_rank14_t2_v1.json'
SELECTION=ART/'low_height_mw_sublattices_v1_302_selection.json'
OUT=ART/'curve302_inverse_kihara_and_rank14_16_intake_v1.json'


def digest(p):return sha256(p.read_bytes()).hexdigest()


def model():
    R=PolynomialRing(QQ,'t');t=R.gen();S=PolynomialRing(R,'x');x=S.gen()
    p=t**2*(8+3*t**2);q=-6*(2+t**2)*(4+t**2)
    U=4*(2+t**2)*(2304+2400*t**2+928*t**4+150*t**6+9*t**8)*(1152+1632*t**2+860*t**4+201*t**6+18*t**8)
    a=[R(0),(2*p**2+p*q+2*q**2)**2,2*(p+q)**2*(2*p**2+p*q+q**2),q**2*(4*p**2-p*q+4*q**2),p*(2*p-q)*(2*p**2+4*p*q+5*q**2),4*p**4+8*p**3*q+9*p**2*q**2-2*p*q**3+2*q**4]
    # x'=t*x, y'=t^6*y removes the printed u=U/t denominator.
    roots=[t*ai+sign*U for ai in a for sign in [-1,1]]
    f=prod(x-bi for bi in roots);g=x**6
    for k in range(5,-1,-1):g+=(f[6+k]-(g*g)[6+k])/2*x**k
    quartic=g*g-f;assert quartic.degree()==4
    assert all(quartic(bi)==g(bi)**2 for bi in roots)
    coefficients=[quartic[i] for i in range(5)]
    # Independent rational-only implementation of the source formula.
    for value in [-3,-1,1,2,3]:
        source=kihara_specialization(Fraction(value))
        assert all(coefficients[i](value)==QQ(source.quartic_coefficients[i])*QQ(value)**(12-i) for i in range(5))
    e,d,c,b,a=coefficients
    I=12*a*e-3*b*d+c*c
    J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    j=R.fraction_field()(6912*I**3/(4*I**3-J**2))
    numerator=j.numerator();denominator=j.denominator()
    assert numerator.gcd(denominator)==1 and numerator.degree()==432 and denominator.degree()==426
    return R,coefficients,I,J,numerator,denominator


def intake():
    d=json.loads(SELECTION.read_text())
    core=min((e for e in d['finalists'] if e['rank']==17),key=lambda e:Decimal(e['determinant']))
    B=matrix(ZZ,core['primitive_basis_rows']);assert B.rank()==17
    result=[]
    for e in d['finalists']:
        if e['rank'] not in [14,15,16]:continue
        C=matrix(ZZ,e['primitive_basis_rows']);rank=e['rank'];assert C.rank()==rank
        assert all(abs(n)==1 for n in C.smith_form()[0].diagonal())
        union=int(B.stack(C).rank())
        result.append({'rank':rank,'candidate_index':e['candidate_index'],
            'primitive_basis_rows':[list(map(int,v)) for v in C.rows()],
            'intersection_rank_with_tested_rank17_core':17+rank-union,
            'contained_in_tested_rank17_core':union==17,'selection_channel':e['channel']})
    assert len(result)==11 and sum(not e['contained_in_tested_rank17_core'] for e in result)==5
    return result


def build():
    R,quartic,I,J,numerator,denominator=model();t=R.gen()
    manifest=json.loads(CONTROL.read_text());verify_kihara_rank14_manifest(manifest)
    replay=kihara_rank14_replay(Fraction(manifest['specialization']['parameter_t']))
    control=EllipticCurve(QQ,list(map(QQ,replay.weierstrass_coefficients)))
    comparison=numerator-control.j_invariant()*denominator
    assert all(comparison[i]==0 for i in range(1,comparison.degree()+1,2))
    RZ=PolynomialRing(QQ,'z');h=RZ([comparison[2*i] for i in range(comparison.degree()//2+1)])
    roots=h.roots(QQ,multiplicities=False);matches=[]
    for z in roots:
        if z<0 or not z.is_square():continue
        for value in sorted({z.sqrt(),-z.sqrt()}):
            a=-27*I(value);b=-27*J(value)
            if 4*a**3+27*b**2==0:continue
            fibre=EllipticCurve(QQ,[a,b]);maps=fibre.isomorphisms(control)
            if maps:matches.append({'t':str(value),'isomorphisms_u_r_s_t':[list(map(str,iso.tuple())) for iso in maps]})
    assert roots==[QQ(4)] and {e['t'] for e in matches}=={'-2','2'}
    target=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    target_j=target.j_invariant();f=numerator*target_j.denominator()-denominator*target_j.numerator()
    f=PolynomialRing(ZZ,'t')(f*f.denominator());f=f//f.content()
    degree=max(numerator.degree(),denominator.degree());assert f.degree()==degree
    prime=101;RP=PolynomialRing(GF(prime),'t');fp=RP(f)
    values=[int(fp(value)) for value in GF(prime)]
    infinity=int(fp[degree])
    assert infinity and all(values) and fp.gcd(RP.gen()**prime-RP.gen())==1
    paths=[Path(__file__),CAS/'kihara_rank14.py',ROOT/'elliptic-curves/ecsearch/kihara.py',
           CAS/'icarm_curve302.py',CONTROL,SELECTION]
    return {'schema':'curve302.inverse-kihara-rank-intake.v1','status':'PASS_CONTROL_AND_FULL_KIHARA_LINE_EXCLUSION',
        'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in paths},
        'limits':{'families':1,'seconds':180,'target_prime':prime,'subspaces':11},
        'source':'https://doi.org/10.3792/pjaa.77.50',
        'quartic_scaled_coefficients_low_to_high':[list(map(str,c.list())) for c in quartic],
        'jacobian_model':'y^2=x^3-27 I(t)x-27 J(t), I,J the binary quartic invariants',
        'j_numerator_coefficients_low_to_high':list(map(str,numerator.list())),
        'j_denominator_coefficients_low_to_high':list(map(str,denominator.list())),
        'j_degree':int(degree),'control_rational_roots_in_t_squared':list(map(str,roots)),
        'control_isomorphism_matches':matches,'control_rank14_manifest_replayed':True,
        'target_j':str(target_j),'comparison_primitive_integer_coefficients_low_to_high':list(map(str,f.list())),
        'projective_no_root_witness':{'prime':prime,'degree':int(degree),'finite_values_in_residue_order':values,'infinity_value':infinity},
        'rank14_16_intake':intake(),
        'boundary':'No rational parameter of this one published Kihara rank-at-least14 family has j equal to302, including the projective limits. This is not an exclusion of every MW14 parent. Candidate subgroup ranks14,15,16 are not certified generic ranks; five candidate spaces extend outside the previously tested rank17 core. No parent of302 is recovered.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();signal.alarm(180)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])
