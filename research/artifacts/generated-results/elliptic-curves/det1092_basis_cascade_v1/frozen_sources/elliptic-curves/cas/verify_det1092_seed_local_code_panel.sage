#!/usr/bin/env sage-python
"""Join independently replayed generic-only codes and first-seed evaluation."""
import hashlib,json,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_seed_local_code_v5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def run():
    protocol=read(DIR/'protocol.json');provenance(protocol)
    roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
    rows=[];paths=[DIR/'protocol.json',Path(__file__)]
    for i in range(9):
        path=DIR/f'case-{i:02d}-replay.json';r=read(path);provenance(r)
        assert r['status']=='PASS_INDEPENDENT_GENERIC_LOCAL_CODE'
        assert r['case']==protocol['cases'][i]
        assert r['case']['parameter']==roster['cases'][i]['parameter']
        assert r['all_generic_local_images_full']
        rows.append({'index':i,'label':roster['cases'][i]['label'],
            'parameter':r['case']['parameter'],
            'local_product_dimension':r['product_dimension'],
            'generic_joint_image_rank':r['generic_joint_rank'],
            'compatibility_check_dimension':r['compatibility_check_dimension'],
            'all_individual_local_images_filled_by_generic_sections':True,
            'extra_point_existence':'NOT_TESTED_BY_THIS_GENERIC_ONLY_COMPUTATION'})
        paths.append(path)
    seed_path=DIR/'first-seed-replay.json';seed=read(seed_path);provenance(seed)
    assert seed['status']=='PASS_INDEPENDENT_FIRST_SEED_LOCAL_COMPATIBILITY_DEFECT'
    assert seed['smallest_separator']['places']==[7,19,23,29,167]
    assert seed['smallest_separating_place_count_in_fixed_set']==5
    paths.append(seed_path)
    result={'status':'PASS_NINE_GENERIC_CODES_AND_HISTORICAL_FIRST_SEED',
        'classification':'verified application and new diagnostic deduction',
        'cases':rows,'first_seed':{'syndrome':seed['syndrome'],
            'generic_joint_rank':17,'with_seed_rank':18,
            'minimum_separating_place_count_in_fixed_set':5,
            'separating_places':[7,19,23,29,167],
            'local_seed_characters':[0,0,1,1,1]},
        'exact_obstructions':[
            'Every individual local point image is already generic on all nine addresses.',
            'Scale0290097 has zero joint compatibility quotient at the frozen places, so this footprint cannot detect any new point there.',
            'Scale0748009 has six checks versus five on302; check counts are not a rank or rational-seed criterion.',
            'Nonzero syndrome certifies a supplied rational point only after the generic saturation gate; it does not produce a rational point or establish a Selmer class.'
        ],
        'input_boundary':'Generic equations and17 sections only for code construction. Exactly the historical first302 point enters evaluation; no later point or V3 artifact.',
        'scope_limits':protocol['limits'],
        'failure_history':['v1 matrix coercion','v2 unavailable Sturm method',
            'v3 seed-only matrix coercion and control denominator clearing',
            'v4 case03 and06 time caps; v5 rational scaling-only retry'],
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    path=DIR/'panel-replay.json'
    if path.exists():assert path.read_text()==payload
    else:path.write_text(payload)
    print('PASS_NINE_GENERIC_CODES_AND_HISTORICAL_FIRST_SEED',flush=True)
    for r in rows:print(r['label'],r['local_product_dimension'],r['generic_joint_image_rank'],r['compatibility_check_dimension'])
if __name__=='__main__':
    signal.alarm(25)
    run()
