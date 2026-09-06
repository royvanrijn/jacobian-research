#!/usr/bin/env python3
"""Norms for the fourth native lift on three fixed finite relation schemes."""
import argparse
from pathlib import Path
from sage.all import QQ,PolynomialRing
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'FOURTH_LIFT_ON_RELATION_PROTOCOL.json'
FORMS=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
SCHEMES=[r.OUT/f'rank_jump_triple_translate_control_{i}_v1.json' for i in range(3)]
VERIFIED=r.OUT/'rank_jump_triple_translate_controls_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_fourth_lift_on_relation_norms_v1.json'


def compute():
    inp=r.read(FORMS);case=next(c for c in inp['cases'] if c['id']=='08234-003')
    fourth=next(c for c in case['covers'] if c['label']=='orbit-13109')
    R=PolynomialRing(QQ,'t');q=R(fourth['form']);rows=[]
    for i,path in enumerate(SCHEMES):
        scheme=r.read(path)
        for p,sha in scheme['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
        factors=[]
        for fdata in scheme['factorization']:
            f=R(fdata['coefficients']);assert fdata['multiplicity']==1 and f.is_monic() and f.is_irreducible()
            assert f.gcd(q)==1
            norm=f.resultant(q);square=bool(norm.is_square())
            row={'parameter_factor':list(map(str,f.list())),'parameter_degree':int(f.degree()),
                 'fourth_value_norm':str(norm),'norm_is_square':square,
                 'fourth_lift_over_residue_field':'NO' if not square else 'YES' if f.degree()==1 else 'UNKNOWN'}
            if square:row['norm_square_root']=str(norm.sqrt())
            if f.degree()==1:
                t=-f[0];assert norm==q(t);row['parameter']=str(t)
                if square:row['fourth_root']=str(q(t).sqrt())
            factors.append(row)
        rows.append({'relation_index':i,'components':factors})
    paths=[FORMS,VERIFIED,PROTOCOL,Path(__file__),HERE/'retrospective.py']+SCHEMES
    return {'schema':'rank-jump.fourth-lift-on-relation-norms.v1','status':'PASS','layer':'solubility',
        'fourth_cover':fourth,'rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    for row in d['rows']:print(row['relation_index'],[(c['parameter_degree'],c['norm_is_square'],c['fourth_lift_over_residue_field']) for c in row['components']])
