#!/usr/bin/env python3
"""Independent positive control exactly in characteristic equal to degree."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from verify_projective13_covers import check_homology,check_augmentation

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'characteristic_degree_control.json').read_text())
    for p,digest in data['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    free,torsion=check_homology(data);assert free==0 and torsion==[13]*26
    _,dimension=check_augmentation(data['augmentation_F13'],data['generator_images'],data['presentation_words'])
    assert dimension==36 and data['base_F13_betti']==3 and len(torsion)-3==data['gain_F13']==23
    assert data['split_formula_prediction']==12+23==35<dimension
    rejects=[]
    bad=deepcopy(data);bad['H1_torsion_factors']=bad['H1_torsion_factors'][:-1]
    try:check_homology(bad)
    except AssertionError:rejects.append('omitted_C13_factor')
    else:raise AssertionError('accepted wrong positive control homology')
    bad=deepcopy(data['augmentation_F13']);bad['kernel_basis']=bad['kernel_basis'][:-1]
    try:check_augmentation(bad,data['generator_images'],data['presentation_words'])
    except AssertionError:rejects.append('incomplete_nonsplit_kernel')
    else:raise AssertionError('accepted incomplete nonsplit kernel')
    files=['characteristic_degree_control.json','verify_characteristic_degree_control.py',
        'verify_projective13_covers.py','verify_core_subgroup_homology.py','verify_integer_lifts.py']
    receipt={'status':'PASS independent positive degree-characteristic control','degree':13,'characteristic':13,
        'H1_torsion_factors':torsion,'base_b1':3,'cover_b1':26,'gain':23,'augmentation_Z1':36,
        'split_prediction':35,'rank_lower_bound':3,'rejected_corruptions':rejects,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Finite positive control; the all-field rank inequality and every-odd-prime family are written uniform arguments, not inferred from this single sample.'}
    (HERE/'characteristic_degree_control_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','gain','augmentation_Z1','split_prediction','rejected_corruptions']},indent=2))


if __name__=='__main__':main()
