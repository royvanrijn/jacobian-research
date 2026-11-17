#!/usr/bin/env python3
"""Positive nonsplit-characteristic control for the all-field Betti gain bound."""
import hashlib
import json
from pathlib import Path
from audit_projective13_covers import homology,augmentation_cocycles

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'characteristic_degree_control_plan.json').read_text())
    images=[tuple(list(range(1,13))+[0]),tuple(range(13)),tuple(range(13))]
    data=homology(images,plan['presentation'])
    assert data['H1_free_rank']==0 and data['H1_torsion_factors']==[13]*26
    module=augmentation_cocycles(images,plan['presentation'],13)
    assert module['Z1_dimension']==36 and module['H1_dimension']==25 and module['fixed_dimension']==1
    files=['characteristic_degree_control_plan.json','audit_characteristic_degree_control.py',
        'audit_projective13_covers.py','audit_core_subgroup_homology.py','integer_lifts.py','verify_integer_lifts.py']
    report={'status':'CONSTRUCTED positive degree-characteristic control',**data,
        'augmentation_F13':module,'base_F13_betti':3,'cover_F13_betti':26,'gain_F13':23,
        'split_formula_prediction':35,'actual_augmentation_Z1':36,'rank_lower_bound':3,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'characteristic_degree_control.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'gain':23,'actual_Z1':36,'split_formula_prediction':35},indent=2))


if __name__=='__main__':main()
