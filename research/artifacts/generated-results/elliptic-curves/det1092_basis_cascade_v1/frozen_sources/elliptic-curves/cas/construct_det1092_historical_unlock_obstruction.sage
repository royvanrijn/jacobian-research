#!/usr/bin/env sage-python
"""Verified application: fixed historical first-centre RR obstruction.

Only generic parent/sections and the saved pre-search centre word are inputs.
The exceptional point, its chart coordinate and splitting outcomes are not
construction inputs. Selection of this one historical centre is retrospective.
Limit: one 19x20 RR system, one polynomial factorization, zero point searches.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,prime_range
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_historical_unlock_obstruction_v1.json'
WORD=[1,-1,-2,1,0,-1,1,-1,1,-1,1,0,1,0,0,-1,0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    path=ART/'curve302_recovered_mw17_parent_v1.json'
    loader=CAS/'load_curve302_recovered_parent.sage'
    source=CAS/'construct_curve302_parent_cheapest_lattice_bisection.sage'
    parent=json.loads(path.read_text());rr=runpy.run_path(str(source))
    E,basis,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](path)
    R=E.base_ring().ring();w=vector(ZZ,WORD);G=matrix(QQ,parent['generic_height_gram'])
    assert w*G*w==10
    trace=-sum((n*p for n,p in zip(w,basis)),E(0))
    relation=rr['primitive_kernel_relation'](trace,R)
    residual=rr['residual_quadratic'](E,trace,relation)
    c,b,a=residual.list();disc=R(b*b-4*a*c)
    factors=disc.factor();q=R(factors.unit());h=R(1)
    for f,e in factors:q*=f**(e%2);h*=f**(e//2)
    assert disc==h*h*q and q.degree()==2 and q.gcd(q.derivative()).degree()==0
    zero=QQ(q(0));assert zero and not zero.is_square()
    # An exact small-prime obstruction needs no factorization of q(0).
    witness=None
    for p in prime_range(3,200):
        if zero.denominator()%p==0:continue
        residue=int(zero.numerator()%p)*pow(int(zero.denominator()%p),-1,int(p))%int(p)
        if pow(residue,(int(p)-1)//2,int(p))==int(p)-1:
            witness={'prime':int(p),'residue':residue};break
    assert witness is not None
    assert a(0) and relation['f2'](0) and h(0) and E.discriminant()(0)
    return {'classification':'verified application','status':'PASS_HISTORICAL_FIRST_CENTRE_NONSPLIT',
        'selection':'Fixed generic-only word from historical first gaining chart; retrospective centre selection, no exceptional point input.',
        'limits':{'RR_systems':1,'matrix_shape':[19,20],'polynomial_factorizations':1,'nonsquare_prime_bound':199,'point_searches':0,'parameter_sweeps':0},
        'centre_word':WORD,'trace_word':list(map(int,-w)),'generic_norm':10,
        'line_coefficients':{k:rr['polynomial_record'](relation[k]) for k in ['f0','f1','f2']},
        'residual_coefficients':[rr['function_record'](v) for v in residual.list()],
        'q_coefficients':rr['polynomial_record'](q),'h_coefficients':rr['polynomial_record'](h),
        'zero_q':str(zero),'nonsquare_witness':witness,'geometric_genus':0,
        'splitting_condition':'Off zeros/poles of a,f2,h and the parent discriminant, q(t) is a square in Q.',
        'lift':{'x':'(-b+h*s)/(2*a)','y':'-(f0+f1*x)/f2','base':'s^2=q(t)'},
        'conclusion':'No rational point of this residual bisection lies above t=0. In particular it cannot supply the historical first rational exceptional point. This excludes one specified multisection, not all multisections or translated points.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [path,loader,source,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],d['nonsquare_witness'],flush=True)
