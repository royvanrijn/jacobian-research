#!/usr/bin/env sage-python
"""Replay all ten independent finite-ring image certificates under one cap."""
import hashlib,json,runpy,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_images_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    script=Path(__file__).with_name('verify_det1092_rr_dyadic_images.sage')
    check=runpy.run_path(str(script),run_name='independent_finite_ring')['verify']
    results=[];paths=[script,OUT/'protocol.json',OUT/'replay-failure-01.json']
    for i in range(10):
        result=check(i);p=OUT/('case-%02d-finite-replay.json'%i);paths.append(p)
        assert result['status']=='PASS_INDEPENDENT_COMPLETE_DYADIC_KUMMER_IMAGE'
        results.append({'case_index':i,'fake_image_rank':result['generic_fake_image_rank'],
                        'true_image_rank':result['generic_true_image_rank_including_D0'],
                        'generic_divisor_indices':result['independent_generic_divisor_indices'],
                        'include_D0':result['include_inherited_D0'],
                        'square_residue_count':len(result['unit_square_residues_mod8']),
                        'nonsquare_tests':len(result['nonsquare_tests'])})
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_TEN_COMPLETE_DYADIC_KUMMER_IMAGES',
            'cases':results,'all_images_spanned_by_inherited_subgroups':True,
            'nonsquare_tests':sum(r['nonsquare_tests'] for r in results),
            'limits':{'wall_seconds':25,'cases':10,'mod4_residue_vectors_per_case':4096,
                      'multiplier_candidates_per_case':64,'global_class_groups':0,
                      'full_discriminant_factorizations':0,'point_searches':0,'pilot_changes':0},
            'scope':'All true/fake2-adic images are complete. The proxy preserves these local squareclasses with proved error margins, not the global sextic field or global Selmer group. Full global descent remains uncomputed.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();out=OUT/'panel-finite-replay.json'
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if out.exists():assert out.read_text()==payload
    else:out.write_text(payload)
    print(result['status'],flush=True)
