#!/usr/bin/env sage-python
"""Replay both independent full-Picard comparisons under one25-second cap."""
import hashlib,json,runpy,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_genus1_picard_panel_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    scripts=[CAS/'verify_det1092_genus1_picard_image.sage',CAS/'verify_det1092_genus1_picard_controls.sage']
    for script in scripts:runpy.run_path(str(script))['verify']()
    paths=[ART/'det1092_genus1_picard_image_v1/replay.json',ART/'det1092_genus1_picard_controls_v1/replay.json']
    first,controls=[json.loads(p.read_text()) for p in paths]
    assert first['first_augmented_rank']==13 and len(controls['cases'])==9
    assert all(c['generic_Picard_image_rank']==12 and c['rank_with_marked_class']==13 for c in controls['cases'])
    result={'status':'PASS_FULL_PICARD_FIRST_SEED_OBSTRUCTION_AND_NINE_CONTROLS',
      'classification':'new exact arithmetic and specificity obstructions; independent replay',
      'inherited_rank':12,'first_augmented_rank':13,'generic_point_control_augmented_ranks':[13]*9,
      'conclusion':'The first carrier point is outside the full inherited rational span, but so are all nine generic-old-point controls. Carrier-class non-genericity is not an original-fibre seed discriminator. Inherited carrier group operations cannot generate the first point.',
      'full_Selmer_groups':'NOT_COMPUTED','nonzero_Sha_classes':'NOT_CLAIMED',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in scripts+paths+[Path(__file__)]}}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(result['status'],flush=True)
if __name__=='__main__':
    signal.alarm(25);main()
