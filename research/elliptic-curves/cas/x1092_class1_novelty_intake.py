#!/usr/bin/env python3
"""Closed intake for the frozen parent while a novelty verifier is unavailable.

No generic replacement is inferred from Agent 1's fixed-number-field recipe
checker. In particular, certificate status strings are never mathematics.
There is deliberately no launch function or automatic positive adapter.
"""
import argparse
import json
from pathlib import Path
from freeze_x1092_class1_arithmetic import check, DEST, sha, write_once

def decision(candidate=None):
    # No candidate can cross this gate until an applicable exact verifier has
    # actually been reviewed, bound to this parent, integrated and tested.
    return {'schema':'x1092.class1.novelty-intake.v1',
      'status':'UNAVAILABLE_APPLICABLE_PROSPECTIVE_NOVELTY_CHECKER',
      'candidate_supplied':candidate is not None,
      'gates':{'prospective_class_novelty':'UNAVAILABLE',
               'cover_construction':'BLOCKED', 'cover_solubility':'BLOCKED',
               'new_rational_point':'BLOCKED', 'independent_quotient_direction':'BLOCKED'},
      'next_action':'WAIT_FOR_APPLICABLE_EXACT_NOVELTY_VERIFIER',
      'allowed_class_indices':[1],'parameter_panel_enabled':False,'V3_release':False,
      'reason':'Agent 1 confirms no prospective novelty verifier exists yet. The fixed-number-field principal-recipe adapter is neither a Q(s) interface nor a novelty certificate.',
      'future_contract':[
        'Pin the implemented checker and its declared domain; do not infer a Q(s) theorem from number-field calibration.',
        'Bind the exact cubic, marked basis, inherited image, support, candidate and proof dependencies to the frozen parent.',
        'Replay principal-square identities and an applicable novelty certificate modulo the complete inherited image.',
        'Only a prospective novelty PASS may dispatch bounded cover construction and solubility work.',
        'Replay a rational point on the exact fibre and its quotient-independence certificate before allowing V3.'
      ]}

def seal():
    check()
    inputs=DEST/'generic-arithmetic.json'
    data=json.loads(inputs.read_text())
    if data['parent_sha256'] != sha(DEST/'parent.json') or data['freeze_sha256'] != sha(DEST/'freeze.json'):
        raise ValueError('generic input binding mismatch')
    producer=Path(__file__).with_name('prepare_x1092_class1_generic_arithmetic.sage')
    if data['producer_sha256'] != sha(producer):
        raise ValueError('generic arithmetic producer changed')
    result={'schema':'x1092.class1.generic-input-seal.v1',
      'generic_arithmetic_sha256':sha(inputs),'producer_sha256':sha(producer),
      'freeze_sha256':sha(DEST/'freeze.json')}
    write_once(DEST/'generic-input-seal.json',(json.dumps(result,indent=2,sort_keys=True)+'\n').encode())
    return result

def status(candidate=None):
    check()
    old=json.loads((DEST/'generic-input-seal.json').read_text())
    if old != seal():
        raise ValueError('generic input seal mismatch')
    result=decision(candidate)
    result['bindings']={'parent_sha256':sha(DEST/'parent.json'),
                        'generic_arithmetic_sha256':sha(DEST/'generic-arithmetic.json'),
                        'intake_sha256':sha(Path(__file__))}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seal',action='store_true');p.add_argument('--candidate',type=Path)
    a=p.parse_args()
    if a.seal:seal()
    candidate=json.loads(a.candidate.read_text()) if a.candidate else None
    result=status(candidate)
    write_once(DEST/'intake-status.json',(json.dumps(status(),indent=2,sort_keys=True)+'\n').encode())
    print(json.dumps(result,indent=2))
