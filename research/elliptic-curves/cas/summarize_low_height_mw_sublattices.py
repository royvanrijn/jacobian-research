#!/usr/bin/env python3
"""Bind both negative calibrations and the successful arithmetic replays."""
import argparse
from pathlib import Path
import search_low_height_mw_sublattices as b


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    paths=[b.OUT/'low_height_mw_sublattices_replay_v1.json',b.OUT/'common_cover_generated_lattices_v1.json',
           b.OUT/'record_rank17_core_candidates_v1.json']
    for prefix in ('low_height_mw_sublattices_v1','common_cover_mw_sublattices_v1'):
        paths.extend(b.OUT/f'{prefix}_{suffix}.json' for suffix in ('protocol','calibration','245_selection','302_selection'))
    replay=b.read(paths[0]);generated=b.read(paths[1]);old=b.read(paths[2])['curves'][-1]
    assert replay['status']=='PASS_ARITHMETIC_REPLAY_CALIBRATION_FAILED'
    assert generated['status']=='PASS_GENERATED_SUBGROUP_INDEX_AND_REDUCTION_REPLAY'
    rows=[]
    target=b.read(b.OUT/'low_height_mw_sublattices_v1_302_selection.json')
    for r in range(8,21):
        p=min((p for p in target['finalists'] if p['rank']==r),key=lambda p:float(p['determinant']))
        rows.append({k:p[k] for k in ('rank','candidate_index','determinant','minimum','ball_support','lll_heights')})
    new=next(p for p in target['finalists'] if p['candidate_index']==next(x['candidate_index'] for x in rows if x['rank']==17))
    oldrows=list(map(list,zip(*old['saturated_basis_columns_in_public_point_coordinates'])))
    comparison=b.run_gp(f'A={b.gp_matrix(new["primitive_basis_rows"])};B={b.gp_matrix(oldrows)};print(mathnf(A~)==mathnf(B~));')
    assert comparison==['1']
    report={'status':'COMPLETED_BOUNDED_SEARCHES_FAILED_RECOVERY_CALIBRATION',
        'rank302_best_primitive_determinants_by_rank':rows,
        'old_rank17_candidate_exactly_rediscovered':True,
        'calibrations':{p:b.read(b.OUT/f'{p}_calibration.json')['status'] for p in ('low_height_mw_sublattices_v1','common_cover_mw_sublattices_v1')},
        'existing_control_replay':replay['legacy_positive_control']['status'],
        'finalists_replayed':sum(c['finalists_verified'] for c in replay['curves']),
        'generic_family_identification':'UNKNOWN','root_lattice_complement':'UNKNOWN',
        'boundary':'The two new algorithms do not meet the requested Fermigier recovery calibration. Candidate spaces and known-subgroup combinations are retained without provenance or rank promotion.',
        'inputs':{str(p.relative_to(b.ROOT)):b.digest(p) for p in (*paths,Path(__file__))}}
    path=b.OUT/'low_height_mw_sublattices_summary_v1.json'
    if args.check:
        assert b.read(path)==report
    else:
        b.save(path,report)
    print(report['status'])


if __name__=='__main__':
    main()
