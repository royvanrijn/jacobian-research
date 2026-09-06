#!/usr/bin/env python3
"""A bounded local test for a constant fourth-root field and its sign kernel."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prime_range
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'FOURTH_LIFT_SHARED_FIELD_PROTOCOL.json'
NORMS=r.OUT/'rank_jump_fourth_lift_on_relation_norms_v1.json'
FIELDS=r.OUT/'rank_jump_fourth_lift_on_relation_verification_v1.json'
BASE=r.OUT/'rank_jump_triple_translate_controls_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_fourth_lift_shared_field_v1.json'


def compute():
    src=r.read(NORMS);fields=r.read(FIELDS);base=r.read(BASE)
    for d in (src,fields,base):
        for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    c=next(c for row in src['rows'] if row['relation_index']==1 for c in row['components'] if c['parameter_degree']==11)
    R=PolynomialRing(QQ,'t');f=R(c['parameter_factor']);q=R(src['fourth_cover']['form']);N=QQ(c['fourth_value_norm']);D=f.discriminant()
    assert f.degree()==11 and f.is_irreducible() and f.resultant(q)==N
    assert N and D and not N.is_square() and not D.is_square()
    br=next(row for row in base['rows'] if row['index']==1);assert br['nonlinear_Galois_group']=='S11'
    patterns=[]
    for p in (73,79):
        k=GF(p);rp=PolynomialRing(k,'t');ff=rp(f)
        assert ff.degree()==11 and ff.is_squarefree()
        ds=sorted(int(g.degree()) for g,e in ff.factor() for _ in range(int(e)))
        patterns.append({'prime':p,'degrees':ds})
    assert [x['degrees'] for x in patterns]==[[11],[1,2,3,5]]
    witness=None;tested=[]
    for p in prime_range(3,200):
        if any(a.denominator()%p==0 for a in f.list()+[N,D]) or N.numerator()%p==0 or D.numerator()%p==0:continue
        k=GF(p);rp=PolynomialRing(k,'t');ff=rp(f);qq=rp(q)
        assert ff.degree()==11 and ff.is_squarefree()
        roots=sorted(ff.roots(multiplicities=False),key=int);tests=[]
        for t in roots:
            value=qq(t)/k(N)
            if value:
                tests.append({'root':int(t),'value':int(value),'is_square':bool(value.is_square())})
                if not value.is_square() and witness is None:
                    witness={'prime':int(p),'root':int(t),'F_derivative_at_root':int(ff.derivative()(t)),
                        'q_at_root':int(qq(t)),'norm_residue':int(k(N)),
                        'q_over_norm_residue':int(value),'Legendre_value':int(value**((p-1)//2))}
        tested.append({'prime':int(p),'root_tests':tests})
        if witness:break
    same=bool((N/D).is_square())
    ratio_witness=None
    if same:ratio_witness={'rational_square_root':str((N/D).sqrt())}
    else:
        ratio=N/D
        if ratio<0:ratio_witness={'place':'real','sign':-1}
        else:
            for p in prime_range(3,200):
                if ratio.numerator()%p and ratio.denominator()%p and not GF(p)(ratio).is_square():
                    ratio_witness={'prime':int(p),'nonsquare_unit_residue':int(GF(p)(ratio))};break
    status='PASS' if witness and ratio_witness else 'UNKNOWN'
    sign_rank=(10 if same else 11) if status=='PASS' else 'UNKNOWN'
    return {'schema':'rank-jump.fourth-lift-shared-field.v1','status':status,'layer':'solubility',
        'parameter_polynomial':list(map(str,f.list())),'fourth_polynomial':list(map(str,q.list())),
        'norm':str(N),'parameter_discriminant':str(D),'norm_and_discriminant_same_squareclass':same,
        'norm_discriminant_ratio_witness':ratio_witness,'base_S11_patterns':patterns,
        'q_over_norm_local_obstruction':witness,'prime_tests':tested,
        'degree22_field_has_quadratic_subfield':False if witness else 'UNKNOWN',
        'normal_closure_sign_kernel_rank':sign_rank,
        'normal_closure_Galois_group':('C2 wr S11' if not same else 'index-two signed permutation group: product(signs)=sign(permutation)') if status=='PASS' else 'UNKNOWN',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (NORMS,FIELDS,BASE,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print(d['status'],d['normal_closure_Galois_group'],d['q_over_norm_local_obstruction'])
