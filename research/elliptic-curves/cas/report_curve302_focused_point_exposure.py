#!/usr/bin/env python3
"""Fail-closed final focused302 results, combined cloud and height-domain audit."""
import argparse
from pathlib import Path
from fractions import Fraction as F
import curve302_focused_point_exposure_v2 as batch
from research_runtime.store import checkpoint

ROOT,ART,D=batch.ROOT,batch.ART,batch.BATCH;cert=batch.cert
OUT=ART/'curve302_focused_point_exposure_v2.json'

def expected():
    p=batch.protocol();driver=cert.read(D/'run-ledger.json');search=cert.read(D/'ledger.json')
    union=cert.read(D/'union-replay/ledger.json')
    assert driver['status']==search['status']==union['status']=='PASS'
    assert len(search['maps'])==len(search['rows'])==len(driver['proofs'])==4
    paths=[D/n for n in ['protocol.json','intake.json','run-protocol.json','run-ledger.json','ledger.json',
        'seed-replay/result.json','seed-replay/protocol.json','union-replay/ledger.json','union-replay/protocol.json']]
    inputs={};rows=[]
    for row in p['rows']:
        folder=D/row['id'];seed=cert.read(folder/'seed.json');result=cert.read(folder/'result.json')
        cp=cert.read(folder/'certification-protocol.json');proof=cert.read(folder/'certification-ledger.json')
        assert proof['status']=='PASS' and proof['completed_boxes']==49 and len(result['charts'])==49
        assert all(c['search']['status']=='bounded_search_complete' for c in result['charts'])
        assert proof['initial_rank']==row['initial_rank'] and result['rank_lower_bound']<=proof['rank_lower_bound']
        assert all(v==proof['rank_lower_bound'] for v in proof['odd_modulus_ranks'].values())
        for name,h in {**cp['sources'],**cp['inputs'],**proof['certificates']}.items():assert cert.hashed(ROOT/name)==h
        inputs.update(proof['certificates'])
        entry=next(r for r in search['rows'] if r['id']==row['id']);geometry=next(r for r in search['maps'] if r['id']==row['id'])
        stages=[geometry['supervision']]+[s['supervision'] for s in entry['stages']]+[s['supervision'] for s in proof['stages']]
        assert all(s['outcome']=='completed' and s['returncode']==0 for s in stages)
        a,b=map(F,seed['curve'][3:]);j=1728*4*a**3/(4*a**3+27*b*b)
        rows.append(dict(id=row['id'],role=row['role'],initial_rank=row['initial_rank'],
            rank_lower_bound=proof['rank_lower_bound'],relative_seed_gain=proof['rank_lower_bound']-row['initial_rank'],
            completed_boxes=49,point_count=proof['point_count'],j_bits=[abs(j.numerator).bit_length(),j.denominator.bit_length()],
            completed_square_quartic_height_floor=max(1,(abs(j.numerator).bit_length()-20)//6,(j.denominator.bit_length()-11)//6),
            searched_quartic_coefficient_bit_range=[min(c['search']['maximum_coefficient_bits'] for c in result['charts']),max(c['search']['maximum_coefficient_bits'] for c in result['charts'])],
            point_search_cpu_ms=sum(c['search']['search_cpu_ms'] or 0 for c in result['charts']),
            point_worker_seconds=next(s['supervision']['wall_seconds'] for s in entry['stages'] if s['name']=='worker'),
            total_stage_seconds=sum(s['wall_seconds'] for s in stages),finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']}))
        paths += [folder/n for n in ['seed.json','maps.json','result.json','certification-protocol.json','certification-ledger.json']]
    for name,h in union['files'].items():assert cert.hashed(ROOT/name)==h
    up=cert.read(D/'union-replay/protocol.json');assert all(cert.hashed(ROOT/n)==h for n,h in up['sources'].items())
    assert up['input_run_ledger_sha256']==cert.hashed(D/'run-ledger.json')
    inputs.update(union['files']);visibility=cert.read(ART/'curve302_focused_visibility_v1.json')
    for r in visibility['rows']:
        q=D/r['id']/'public-visibility-after-search.json';assert cert.hashed(q)==r['details_sha256'];paths.append(q)
    failed=batch.LOCAL/'curve302-focused-point-exposure-v1/run-ledger.json';failure=cert.read(failed)
    assert failure['status']=='FAILED_OR_CENSORED' and len(failure['stages'])==1
    wasted=sum(s['supervision']['wall_seconds'] for s in failure['stages']);paths.append(failed)
    total=sum(r['total_stage_seconds'] for r in rows)+wasted
    total+=sum(s['supervision']['wall_seconds'] for s in driver['stages'] if not s['name'].startswith('proof-'))
    total+=sum(s['supervision']['wall_seconds'] for s in union['stages'])
    chartpath=ART/'det1092_reduced_parameter_chart_v1/chart-search.json';chart=cert.read(chartpath)
    m=chart['selected']['state']['parameter_matrix'];anchor=-F(m[1])/F(m[0]);paths.append(chartpath)
    inputs.update({str(q.relative_to(ROOT)):cert.hashed(q) for q in paths})
    return dict(schema='elliptic-curves.curve302-focused-point-exposure.v2',status='PASS',rows=rows,
        completed_boxes=196,combined_public31_cloud_lower_bound=union['rank_lower_bound'],
        combined_public31_cloud_points=union['point_count'],certified_gain_above_public31=max(0,union['rank_lower_bound']-31),
        sibling_rank_gains=sum(r['relative_seed_gain'] for r in rows if r['id'].startswith('d1092')),
        total_stage_seconds=total,failed_intake_seconds=wasted,visibility=visibility['rows'],
        parent_height_diagnostic=dict(curve302_reduced_parameter=str(anchor),parameter_height=max(abs(anchor.numerator),anchor.denominator),
            curve302_j_bits=rows[0]['j_bits'],previous_sibling_selection_j_numerator_band=[256,319],
            scope='Post-freeze diagnostic only. The earlier height8 panel and its score band do not sample302-sized arithmetic territory. No sibling parameter or centre is selected from the anchor address.'),
        inputs=inputs,source_sha256=cert.hashed(Path(__file__)),
        boundary='Masked generic17 control and user-directed public31 search are distinct endpoints. Two sibling fibres were selected by certified19 bounds after the unchanged six-fibre comparison, then received new classes outside their old17 span. All maps precede point searches; no validation-prime or public-exceptional-point input enters sibling selection/geometry. Full31 union checks all302 search outputs. Matching finite bounds do not give an exact rank or prove rational-span membership of every point. V1 stopped at a string-versus-Fraction seed comparison before maps or points; V2 preserves it and charges its cost. No automatic larger scan or next wave.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();r=expected()
    if a.check:assert cert.read(OUT)==r
    else:assert not OUT.exists();checkpoint(OUT,r)
    print('PASS196 focused boxes',[(x['id'],x['rank_lower_bound']) for x in r['rows']],
          'combined302',r['combined_public31_cloud_lower_bound'],'seconds',r['total_stage_seconds'],flush=True)
