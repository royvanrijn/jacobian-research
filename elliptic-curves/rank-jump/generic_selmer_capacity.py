#!/usr/bin/env python3
"""Point-masked geometric capacity and historical parent geometry, not a selector."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fixed_field_transfer_geometry as geometry
import fresh_governing_panel as panel

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'GENERIC_SELMER_CAPACITY_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_generic_selmer_capacity_inputs_v1.json'
PROVENANCE=r.OUT/'rank_jump_generic_selmer_capacity_provenance_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_selmer_capacity_geometry_v1.json'
COMPARISON=r.OUT/'rank_jump_generic_selmer_capacity_comparison_v1.json'
VERIFICATION=r.OUT/'rank_jump_generic_selmer_capacity_verification_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-generic-selmer-capacity-v1'


def bind(paths):
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    sources={}
    def read(p):
        sources.update(bind([p]));return r.read(p)
    families=read(geometry.INPUT)['families']
    lineage=read(r.ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json')
    rep=lineage['representative']
    historic={str(x['curve_id']):x for x in lineage['target_isomorphisms'] if x['chart']=='norm12-orbit-074d9'}
    families.append({'family':'historic-lineage-074d9','A':rep['A_coefficients_low_to_high'],
        'B':rep['B_coefficients_low_to_high'],'irreducibility_witness_parameter':historic['356']['parameter']})
    parent=next(x for x in read(r.ROOT/'elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json')['parents'] if x['parent_id']=='curve398-p16875')
    families.append({'family':parent['parent_id'],'A':parent['pencil']['A_coefficients_low_to_high'],
        'B':parent['pencil']['B_coefficients_low_to_high'],'irreducibility_witness_parameter':parent['target_parameter']})
    labels={x['token']:x for x in read(panel.MANIFEST)['rows']}
    old={x['token']:x['local'] for x in read(panel.OUTPUT)['rows']}
    for x in read(r.OUT/'rank_jump_fresh_governing_completion_v1.json')['rows']:
        if x.get('local'):old[x['token']]=x
    cases=[]
    for x in read(panel.INPUT)['cases']:
        token=x['token'];label=labels[token];family=label['family'];t=label['parameter']
        if family=='published-R17':family='historic-lineage-074d9';t=historic[label['id'].split('-')[-1]]['parameter']
        if family=='curve398-p16875':t=parent['target_parameter']
        cases.append({**x,'family':family,'parameter':t,
            'independence_blocks':old[token]['independence_blocks'],
            'independence_signatures':old[token]['independence_signatures']})
    r.write_new(INPUT,{'schema':'rank-jump.generic-selmer-capacity-inputs.v1','families':families,'cases':cases})
    sources.update(bind([INPUT,PROTOCOL,Path(__file__)]))
    r.write_new(PROVENANCE,{'schema':'rank-jump.generic-selmer-capacity-provenance.v1','bindings':sources,
        'boundary':'Export projects equations and explicitly generic coordinates. It does not project any exceptional point or retained rank. Historical lineage supplies family identification only.'})


def worker(family):
    # Reuse an immutable worker without changing its source or original outputs.
    geometry.INPUT=INPUT;geometry.PROTOCOL=PROTOCOL
    result=geometry.worker(family)
    result['bindings']=bind([INPUT,PROTOCOL,Path(__file__),Path(geometry.__file__),Path(r.__file__)])
    r.write_new(WORK/(family+'.json'),result)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    prior=r.read(geometry.OUTPUT);rows=list(prior['rows']);have={x['family'] for x in rows}
    for family in [x['family'] for x in r.read(INPUT)['families'] if x['family'] not in have]:
        path=WORK/(family+'.json')
        if not path.exists():
            error=None
            with (WORK/(family+'.log')).open('x') as log:
                try:
                    p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--family',family],stdout=log,stderr=log,timeout=30)
                    if p.returncode:error='worker failure'
                except subprocess.TimeoutExpired:error='bounded worker timeout'
            if error:r.write_new(path,{'family':family,'status':'UNKNOWN','reason':error})
        row=r.read(path);rows.append(row)
        print(family,row['status'],row.get('geometric_transposition_branch_points'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.generic-selmer-capacity-geometry.v1','rows':rows,
        'bindings':bind([INPUT,PROTOCOL,Path(__file__),geometry.OUTPUT])})


def compare():
    measured=r.read(VERIFICATION);assert measured['status']=='PASS'
    cases={x['token']:x for x in measured['cases']};rows=[]
    for label in r.read(panel.MANIFEST)['rows']:
        x=cases[label['token']];R=label['retained_rank_lower_bound'];d=x['geometric_selmer_dimension']
        deficit=max(0,R-d)
        rows.append({'token':x['token'],'id':label['id'],'family':x['family'],
            'retained_rank_lower_bound':R,'rank_status':'lower bound; exact rank UNKNOWN',
            'generic_mod_two_dimension':x['generic_mod_two_dimension'],
            'geometric_selmer_dimension':d,'arithmetic_global_pool_dimension_upper_bound':d,
            'extra_pool_capacity_upper_bound':d-x['generic_mod_two_dimension'],
            'rational_fibre_dimensions_outside_pool_lower_bound':deficit,
            'global_block_local_obstruction_rank_lower_bound':deficit,
            'good_geometric_base_support_lower_bound_if_no_bad_obstructions':(deficit+1)//2,
            'additional_arithmetic_quotient_CT':'UNKNOWN'})
    r.write_new(COMPARISON,{'schema':'rank-jump.generic-selmer-capacity-comparison.v1','rows':rows,
        'bindings':bind([Path(__file__),VERIFICATION,panel.MANIFEST]),
        'boundary':'Rank labels joined after equation/generic-only verification. Outside-pool dimensions are necessities deduced from labels, not point-free measurements of fibre ranks or additional classes. Generic capacities do not discriminate within a family. Support bounds apply to proposed Q(t) blocks, not to all possible auxiliary constructions.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','worker','capture','compare']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':worker(args.family)
    else:globals()[args.mode]()
