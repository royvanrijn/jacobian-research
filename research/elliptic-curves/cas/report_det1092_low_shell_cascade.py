#!/usr/bin/env python3
"""Terminal report for the target-blind determinant-1092 cascade deployment."""
import argparse, hashlib, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
CAS=ROOT/'elliptic-curves/cas';BATCH=LOCAL/'det1092-low-shell-cascade-v3';OUT=ART/'det1092_low_shell_cascade_v1.json'
OLD=LOCAL/'det1092-record-scale-points-v1/ledger.json';A1=ART/'a1_mw16_target_free_parameter_search_h300_v1.json';BASIS=ART/'det1092_low_shell_basis_randomization_v1.json'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def build():
    protocol,ledger,basis=read(BATCH/'protocol.json'),read(BATCH/'ledger.json'),read(BASIS)
    if ledger['status']!='PASS' or ledger['completed_boxes']!=147 or len(ledger['rows'])!=3:raise ArithmeticError('low-shell deployment is not terminal')
    rows=[];inputs=[BATCH/'protocol.json',BATCH/'ledger.json',BASIS,OLD,A1,Path(__file__)]
    for row in ledger['rows']:
        if row['status']!='PASS' or row['rank_lower_bound']!=17 or row['stop_reason']!='NO_CERTIFIED_GAIN' or len(row['waves'])!=1:raise ArithmeticError('low-shell terminal row differs')
        wave=row['waves'][0];directory=BATCH/row['id']/wave['id'];certificate=read(directory/'certification-ledger.json')
        if certificate['status']!='PASS' or certificate['rank_lower_bound']!=17 or certificate['odd_modulus_ranks']!={'3':17,'5':17}:raise ArithmeticError('low-shell finite certificate differs')
        inputs += [directory/'maps.json',directory/'result.json',directory/'certification-ledger.json']
        rows.append({'id':row['id'],'score_stratum':row['stratum'],'initial_rank':17,'terminal_rank_lower_bound':17,'completed_boxes':49,'stop_reason':row['stop_reason'],'mod2_mod3_mod5':[17,17,17],'mod2_cloud':wave['cloud_path']})
    old=read(OLD)
    if old['status']!='PASS' or old['completed_boxes']!=2352 or len(old['rows'])!=48 or any(row['rank_lower_bound']!=17 for row in old['rows']):raise ArithmeticError('earlier determinant-1092 negative exposure differs')
    a1=read(A1)
    if a1['status']!='PASS_COMPLETE_TARGET_FREE_A1_MW16_PARAMETER_SEARCH' or a1['completed_candidate_count']!=104 or a1['positive_candidate_count']!=0:raise ArithmeticError('existing A1/MW16 terminal negative exposure differs')
    if basis['status']!='PASS_COORDINATE_SENSITIVITY_AUDIT' or basis['policy_is_basis_invariant'] is not False:raise ArithmeticError('basis sensitivity conclusion differs')
    return {'schema':'elliptic-curves.det1092-low-shell-cascade-report.v1','status':'PASS_COMPLETE_TARGET_BLIND_LOW_SHELL_DEPLOYMENT','inputs':{str(path.relative_to(ROOT)):sha(path) for path in inputs},'new_low_shell_deployment':{'frozen_policy':protocol['adaptive_policy'],'selection_boundary':protocol['selection'],'rows':rows,'complete_negative_exposure':True},'matched_existing_parent_exposures':{'det1092_deep_centre_48':{'completed_boxes':old['completed_boxes'],'initial_rank':17,'terminal_rank_lower_bound':17,'adaptive_waves_entered':0},'a1_mw16_target_free_104':{'completed_candidates':a1['completed_candidate_count'],'generic_rank':16,'positive_candidates':a1['positive_candidate_count'],'adaptive_waves_entered':0}},'basis_randomization':{'canonical_norm_transport_preserved':True,'policy_is_basis_invariant':False,'selected_coset_intersection_size':basis['selected_intersection_size'],'consequence':'Do not claim intrinsic basis invariance for coordinate-hashed sampling. Replace the schedule by canonical parity addressing before a next deployment.'},'rank32_status':'No certified gain occurred, so no subgroup enlargement or later wave was authorized. Rank at least 32 remains open.','claim_boundary':'The new three-fibre low-shell deployment is target-blind and complete. The two matched existing exposures are terminal negative controls, not inputs to its selection or execution. No bounded miss is an upper rank bound.','reproducing_command':'python3 elliptic-curves/cas/report_det1092_low_shell_cascade.py --check'}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();payload=build();rendered=json.dumps(payload,indent=2,sort_keys=True)+'\n'
    if args.check:
        if OUT.read_text()!=rendered:raise ArithmeticError('low-shell cascade report did not replay')
    else:
        if OUT.exists():raise FileExistsError('preserve terminal low-shell cascade report')
        OUT.write_text(rendered)
    print('DET1092 LOW-SHELL CASCADE|rows=3|boxes=147|gain=0|basis_invariant=false',flush=True)
if __name__=='__main__':main()
