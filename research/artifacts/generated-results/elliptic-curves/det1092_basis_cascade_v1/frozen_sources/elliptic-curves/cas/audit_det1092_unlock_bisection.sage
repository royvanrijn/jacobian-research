#!/usr/bin/env sage-python
"""Verified application: one additional RR bisection, fixed by the first V3 gain.

Retrospective selection, bounded to orbit127449 and its 19x20 RR system.
This is not part of any search policy or the eight-fibre pilot.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_initial_unlock_bisection_obstruction_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    parent_path=ART/'curve302_recovered_mw17_parent_v1.json';audit_path=ART/'det1092_initial_unlock_construction_audit_v1.json'
    source=CAS/'construct_curve302_parent_cheapest_lattice_bisection.sage';loader=CAS/'load_curve302_recovered_parent.sage'
    parent=json.loads(parent_path.read_text());audit=json.loads(audit_path.read_text());rr=runpy.run_path(str(source))
    E,basis,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](parent_path)
    # Literal specialization, not a numerical recognition of generic origin.
    generic=[(36*p[0](0)+15,108*(2*p[1](0)+p[0](0)+1)) for p in basis]
    oldseed=ROOT/'artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2/curve302-generic17/seed.json'
    newseed=ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3/replay-M17/epoch-00/selection.json'
    assert generic==[tuple(map(QQ,p)) for p in json.loads(oldseed.read_text())['points']]
    assert generic==[tuple(map(QQ,p)) for p in json.loads(newseed.read_text())['basis']]
    w=vector(ZZ,audit['autonomous_V3']['centre_word_in_generic17']);G=matrix(QQ,parent['generic_height_gram'])
    assert w*G*w==10
    trace=-sum((n*p for n,p in zip(w,basis)),E(0));R=E.base_ring().ring()
    relation=rr['primitive_kernel_relation'](trace,R);residual=rr['residual_quadratic'](E,trace,relation)
    Z=PolynomialRing(QQ,'x');zero=Z([c(0) for c in residual.list()]);point=list(map(QQ,audit['autonomous_V3']['point_literal302']))
    incidence=sum(relation[k](0)*v for k,v in zip(['f0','f1','f2'],[1,*point]))
    split=zero.discriminant().is_square()
    assert not split and incidence!=0
    return {'classification':'verified application','status':'PASS_EXACT_FIRST_CENTRE_BISECTION_NONSPLIT_OBSTRUCTION',
       'orbit':127449,'selection':'Retrospective: the norm10 generic centre of the completed V3 initial gain; exactly one additional RR equation.',
       'limits':{'additional_RR_systems':1,'matrix_shape':[19,20],'point_searches':0,'parameter_sweeps':0},
       'generic_origin':{'historical_17_section_specializations_exact':True,'V3_17_section_specializations_exact':True},
       'trace_word':list(map(int,-w)),'generic_norm':10,
       'line_coefficients':{k:rr['polynomial_record'](relation[k]) for k in ['f0','f1','f2']},
       'residual_coefficients':[rr['function_record'](c) for c in residual.list()],
       'zero_residual_coefficients':rr['polynomial_record'](zero),'zero_discriminant':str(zero.discriminant()),
       'zero_split_over_Q':bool(split),'winning_point_line_incidence':str(incidence),
       'conclusion':'The winning pointed chart and the unique genus-zero surface bisection share a trace centre but are different constructions. The first rational302 gain does not lie on this bisection. Its residual fibre at0 is nonsplit; it cannot supply a rational302 point via this cover.',
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [parent_path,audit_path,source,loader,oldseed,newseed,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(d['status'],flush=True)
