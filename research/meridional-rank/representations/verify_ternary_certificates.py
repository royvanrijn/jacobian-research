#!/usr/bin/env python3
"""Run and seal the four independent ternary certificate replays.

--receipts-only checks already-completed exact replays and refreshes document
seals, without repeating finite calculations after a layout-only change.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
CHECKERS=['verify_ternary_transvection_control.py','verify_ternary_reflection_control.py',
          'verify_ternary_bases.py','verify_ternary_cover_modules.py']
RECEIPTS=['ternary_transvection_control_verification.json','ternary_reflection_control_verification.json',
          'ternary_base_verification.json','ternary_cover_module_verification.json']


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--receipts-only',action='store_true')
    args=parser.parse_args()
    if not args.receipts_only:
        for script in CHECKERS:
            completed=subprocess.run([sys.executable,str(HERE/script)],cwd=HERE,capture_output=True,text=True)
            if completed.returncode:sys.stdout.write(completed.stdout);sys.stderr.write(completed.stderr);raise SystemExit(completed.returncode)
            print(json.dumps({'checker':script,'status':'PASS'}),flush=True)
    data=[]
    for file in RECEIPTS:
        receipt=json.loads((HERE/file).read_text());assert receipt['status'].startswith('PASS')
        for p,digest in receipt['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest,(file,p)
        data.append(receipt)
    first,second,bases,covers=data
    assert first['checks']['root_dimension']==second['checks']['root_dimension']==4
    assert [c['tuple_orbits'] for c in bases['classes']]==[21632,34084]
    assert all(b['all_admissible_affine_lifts_principal'] for c in bases['classes'] for row in c['targets'] for b in row['base_orbits'])
    assert covers['representation_orbits']==16 and covers['module_tests']==32
    assert all(m['Z1_dimension']==4 and m['H1_dimension']==1 for r in covers['results'] for m in r['modules'])
    projective=covers['projective_index13_control'];assert projective['homology_gain']==15 and projective['sharp_gain_threshold']==13
    reject_count=len(first['rejected_corruptions'])+int(second['wrong_root_coefficient_rejected'])+len(bases['rejected_corruptions'])+int(covers['false_coefficient_rejected'])+int(covers['incomplete_kernel_rejected'])
    assert reject_count==10
    files=CHECKERS+RECEIPTS+['verify_ternary_certificates.py','TERNARY_AFFINE_AND_COVER_CERTIFICATES.tex',
        'output/pdf/TERNARY_AFFINE_AND_COVER_CERTIFICATES.pdf']
    report={'status':'PASS scoped complete ternary affine and cover package',
        'targets':16,'rank_one_classes':2,'tuple_orbits':[21632,34084],'raw_tuples':[1124864,1601613],
        'all_rank_one_admissible_affine_lifts_principal':True,
        'full_GL3_F3_representation_orbits':16,'knots_with_full_GL3_F3_image':8,
        'induced_cover_quotient':'SL3(F3),order5616','ordinary_cover_module_tests':32,
        'ordinary_Z1_dimension':4,'required_Z1_dimension':7,'control_root_dimensions':[4,4],
        'projective_index13_control':projective,'rejected_corruptions':reject_count,
        'new_target_closures':0,'full_MRC_proved':False,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'All rank-one GL3(F3) bases on the unchanged original16, and natural/dualF3 ordinary cover modules through the16 full-image representations. Finite inputs independently replayed; uniform cocycle/fixed-space, core-cover and Kurosh arguments require the written proof. Unformalized and externally unreviewed. No distinct-kernel, all-field/all-module or worldwide novelty claim; target13point homology untested.'}
    (HERE/'ternary_certificate_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','targets','full_GL3_F3_representation_orbits','ordinary_cover_module_tests','rejected_corruptions','new_target_closures']},indent=2))


if __name__=='__main__':main()
