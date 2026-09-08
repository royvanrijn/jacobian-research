#!/usr/bin/env sage-python
"""Export equations without exceptional coordinates or control outcomes.

Selection of cover01 remains retrospective; sanitizing an execution payload
does not make its selection prospective.25-second cap, no point evaluation.
"""
import hashlib,json,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves';DIR=ART/'det1092_norm8_seed_cover_v2'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def main():
    parent_path=ART/'curve302_recovered_mw17_parent_v1.json'
    parent=read(parent_path);proof=read(DIR/'replay.json')
    assert proof['status']=='PASS_INDEPENDENT_POSITIVE_RANK_GENUS_ONE_FIRST_SEED_COVER'
    for path,digest in proof['inputs'].items():assert sha(ROOT/path)==digest
    for i in range(2):
        cover_path=DIR/f'cover-{i:02d}.json';cover=read(cover_path)
        payload={'schema':'det1092-genus-one-carrier-equations-v1',
          'classification':cover['classification'],
          'selection_boundary':'Cover00 is generic-selected. Cover01 is first-witness-calibrated. Neither label is changed by removing point coordinates from execution inputs.',
          'cover_index':i,'cover_parameter_z':cover['z'],
          'quartic_coefficients':cover['quartic_coefficients'],'maps':cover['maps'],
          'parent_a_invariants':parent['a_invariants'],
          'generic_point_functions':parent['basis_weierstrass_coordinates'],
          'infinite_base_point_source':proof['covers'][i]['explicit_base_generator'],
          'function_field_rank_lower_bound':18,
          'finite_certificate_primes':[17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197],
          'no_two_torsion_prime':31,
          'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
                    [cover_path,parent_path,DIR/'replay.json',Path(__file__)]}}
        retain(DIR/f'equation-only-cover-{i:02d}.json',payload)
    print('EXPORTED_TWO_EQUATION_ONLY_CARRIER_PAYLOADS',flush=True)
if __name__=='__main__':
    signal.alarm(25);main()
