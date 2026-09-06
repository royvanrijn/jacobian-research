#!/usr/bin/env python3
"""Exact local epsilon signs refine one fixed Frobenius moment bound."""
import argparse
from math import prod
from pathlib import Path
from sage.all import GF, PolynomialRing
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_TWIST_MOMENT_PARITY_PROTOCOL.json'
SOURCE=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
OUTPUT=r.OUT/'rank_jump_native_twist_moment_parity_v1.json'


def compute():
    data=r.read(SOURCE)
    for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(GF(131),'t');A,B,q=[R(x) for x in data['geometry']['modular_coefficients']]
    delta=-16*(4*A**3+27*B**2);factors=list(delta.factor())
    assert prod(f**e for f,e in factors)*delta.leading_coefficient()==delta
    rows=[]
    for f,e in factors:
        assert e==1 and f.is_irreducible()
        norms=[f.resultant(864*B),f.resultant(864*B*q**3)]
        assert all(norms)
        chars=[int(a**65) for a in norms];chars=[-1 if a==130 else a for a in chars]
        assert all(a in(-1,1) for a in chars)
        rows.append({'factor':list(map(int,f.list())),'degree':int(f.degree()),'norms_minus_c6':list(map(int,norms)),
                     'split_characters':chars,'local_root_numbers':[-a for a in chars]})
    assert sum(x['degree'] for x in rows)==24
    result=[]
    for j,old in enumerate(data['rows']):
        W=prod(row['local_root_numbers'][j] for row in rows)
        upper=old['arithmetic_generic_rank_upper_bound'];refined=max(m for m in range(upper+1) if (-1)**m==W)
        assert old['arithmetic_generic_rank_lower_bound']<=refined
        result.append({'id':old['id'],'global_root_number':W,'central_multiplicity_parity':'odd' if W==-1 else 'even',
                       'moment_bound':upper,'parity_refined_upper_bound':refined})
    assert result[0]['global_root_number']==-1 and result[0]['parity_refined_upper_bound']==17
    return {'schema':'rank-jump.native-twist-moment-parity.v1','status':'PASS','multiplicative_places':rows,
            'twisted_additive_root_number_product':1,'additive_argument':'Two total geometric I0* places: product chi_residue(-1)=chi_131(-1)^2=1.',
            'infinity_root_number':1,'rows':result,
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,SOURCE,Path(__file__),HERE/'retrospective.py')},
            'boundary':'Parity of the functional-equation central multiplicity only. The rank bound uses rank<=analytic multiplicity over a finite function field and good surface reduction, not a number-field parity conjecture.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=compute()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    for row in result['rows']:print(row['id'],'sign',row['global_root_number'],'upper',row['parity_refined_upper_bound'])
